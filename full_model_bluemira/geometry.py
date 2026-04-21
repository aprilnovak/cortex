#!/usr/bin/env python3
"""
geometry.py
===========
Builds the OpenMC geometry, mixes materials, computes chunk mappings
and structural origin maps.

Imported by neutronics_model.py, depletion_model.py, and post-processing scripts.

Material definitions live in inputs.py (materials) and eudemo_materials.py (recipes).

Exported names
--------------
model, pydagmc_model, dagmc_universe, dagmc_universe_cells, all_cells, sector_cell
neutron_source_rate, s_in_y, ev_to_joule
build_breeder_chunks, ob_by_key, ib_by_key, n_breeder, OB_CHUNK_SIZE, IB_CHUNK_SIZE
cell_ids_all
info, all_surface_ids, external_surface_ids, internal_surface_ids
materials  (the low-level project materials module, for Ed() etc.)
MIX_RECIPES_NORM, structure_material_list
"""

from __future__ import annotations

import json
import math
import re
import sys
from enum import Enum
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import openmc
import pydagmc

# ──────────────────────────────────────────────────────────────────────────────
# User configuration
# ──────────────────────────────────────────────────────────────────────────────
import inputs as cfg
import materials

# ──────────────────────────────────────────────────────────────────────────────
# Create output directories
# ──────────────────────────────────────────────────────────────────────────────
for _d in [cfg.NEUTRONICS_RUN_DIR, cfg.NEUTRONICS_RESULTS_DIR, cfg.DEPLETION_RUN_DIR]:
    _d.mkdir(parents=True, exist_ok=True)

# ──────────────────────────────────────────────────────────────────────────────
# Source / scaling constants
# ──────────────────────────────────────────────────────────────────────────────
ev_to_joule: float         = 1.60218e-19
_convert_ev_fusion: float  = ev_to_joule * cfg.EV_PER_FUSION
_section_power: float      = cfg.TOTAL_FUSION_POWER_W / cfg.NUMBER_OF_SECTORS
neutron_source_rate: float = _section_power / _convert_ev_fusion
s_in_y: float              = 365.0 * 24.0 * 60.0 * 60.0

# ──────────────────────────────────────────────────────────────────────────────
# DAGMC load
# ──────────────────────────────────────────────────────────────────────────────
if not cfg.DAGMC_MODEL_FILE.is_file():
    raise FileNotFoundError(f"DAGMC file not found: {cfg.DAGMC_MODEL_FILE}")

dagmc_universe = openmc.DAGMCUniverse(filename=str(cfg.DAGMC_MODEL_FILE))
pydagmc_model  = pydagmc.Model(str(cfg.DAGMC_MODEL_FILE))

openmc.reserve_ids([v.id for v in pydagmc_model.volumes],  cls=openmc.Cell)
openmc.reserve_ids([s.id for s in pydagmc_model.surfaces], cls=openmc.Surface)

model = openmc.Model()

# ──────────────────────────────────────────────────────────────────────────────
# Surface-info helpers
# ──────────────────────────────────────────────────────────────────────────────
class _Orientation(Enum):
    FORWARD = 1
    REVERSE = -1

def dagmc_bounding_box(pydagmc_model, vol_id: int):
    volume = pydagmc_model.volumes_by_id[vol_id]
    triangle_coords = volume.triangle_coords
    min_coords = np.min(triangle_coords, axis=0)
    max_coords = np.max(triangle_coords, axis=0)
    return openmc.BoundingBox(min_coords, max_coords)

def dagmc_volume_surface_info(
    pydagmc_model,
    volume_ids: List[int],
    ) -> Tuple[dict, List[int], List[int], List[int]]:
    volume_ids_set = set(volume_ids)
    result = {}

    surface_ids_all = set()
    surface_ids_external = set()
    surface_ids_internal = set()

    for vol_id in volume_ids:
        volume = pydagmc_model.volumes_by_id[vol_id]

        all_list = []
        external_list = []
        internal_list = []

        for surface in volume.surfaces:
            parent_volumes = surface.senses
            if parent_volumes and parent_volumes[0].id == vol_id:
                orientation = _Orientation.FORWARD
            else:
                orientation = _Orientation.REVERSE

            adj_all = [v.id for v in surface.volumes if v.id != vol_id]
            adj_ext = [aid for aid in adj_all if aid not in volume_ids_set]

            entry = {
                "surface_id": surface.id,
                "orientation": orientation,
                "adjacent_volumes_all": adj_all,
                "adjacent_volumes_external": adj_ext,
            }

            all_list.append(entry)
            surface_ids_all.add(surface.id)

            if not adj_all or adj_ext:
                external_list.append(entry)
                surface_ids_external.add(surface.id)
            else:
                internal_list.append(entry)
                surface_ids_internal.add(surface.id)

        result[vol_id] = {
            "all_surfaces": all_list,
            "external_surfaces": external_list,
            "internal_surfaces": internal_list,
        }

    return (
        result,
        sorted(surface_ids_all),
        sorted(surface_ids_external),
        sorted(surface_ids_internal),
    )

# ──────────────────────────────────────────────────────────────────────────────
# GEOMETRY  –  tokamak vs slab
# ──────────────────────────────────────────────────────────────────────────────
def _azimuthal_plane(theta_deg, boundary_type=None, name=None, surface_id=None):
    """
    Generate a CSG plane to cut the tokamak into a wedge. The unit normals of
    these planes are generated such that the correct region definitions later
    in this file will include the small-angle space between the two planes.
    """
    theta = math.radians(theta_deg)
    if abs(theta_deg - 90.0) < 1e-6:
        return openmc.XPlane(boundary_type=boundary_type, name=name, surface_id=surface_id)
    return openmc.Plane(a=-math.tan(theta), b=1.0, c=0.0, d=0.0,
                        boundary_type=boundary_type, name=name, surface_id=surface_id)

def _build_tokamak_geometry() -> openmc.Cell:
    theta0 = cfg.THETA0_DEG
    theta1 = theta0 + 360.0 / cfg.NUMBER_OF_SECTORS
    cut_lo = _azimuthal_plane(theta0, boundary_type="reflective", surface_id=100_001)
    cut_hi = _azimuthal_plane(theta1, boundary_type="reflective", surface_id=100_002)
    z_min  = openmc.ZPlane(z0=-3000, boundary_type="vacuum",      surface_id=100_003)
    z_max  = openmc.ZPlane(z0=+2000, boundary_type="vacuum",      surface_id=100_004)
    r_out  = openmc.ZCylinder(r=3000, boundary_type="vacuum",     surface_id=100_005)
    region = (+cut_lo & -cut_hi) & (+z_min & -z_max) & (-r_out)
    return openmc.Cell(
        cell_id=10_000, region=region, fill=dagmc_universe, name="sector_container"
    )

# Slab cut functions:
# This works for now, but it can break quite easily if we decide to
# use a different breeder chunk (such as IB_1_b4).
# TBD: Improve to avoid these crashes. It may not be needed if dagmc updated.
def _get_breeder_reflective_cuts(
    volume_id: int, base_surface_id: int = 50_000
) -> list:
    """
    Find 2 Z-dominant + 2 Y-dominant reflective cut planes for a DAGMC volume.
    Uses pymoab.
    """
    from pymoab import core as mb_core, types as mb_types

    mb = mb_core.Core()
    mb.load_file(str(cfg.DAGMC_MODEL_FILE))

    surf_map, _, _, _ = dagmc_volume_surface_info([volume_id])
    vol_data          = surf_map[volume_id]

    id_tag  = mb.tag_get_handle("GLOBAL_ID")
    dim_tag = mb.tag_get_handle("GEOM_DIMENSION")
    s_sets  = mb.get_entities_by_type_and_tag(0, mb_types.MBENTITYSET, dim_tag, [2])

    z_cands: list = []
    y_cands: list = []

    for entry in vol_data["all_surfaces"]:
        sid = int(entry["surface_id"])
        ori = entry["orientation"]

        target = next(
            (s for s in s_sets if int(mb.tag_get_data(id_tag, s)[0][0]) == sid), None
        )
        if target is None:
            continue

        ents  = mb.get_entities_by_handle(target, True)
        tris  = [e for e in ents if mb.type_from_handle(e) == mb_types.MBTRI]
        if not tris:
            continue

        conn   = mb.get_connectivity(tris).reshape(-1, 3)
        coords = mb.get_coords(conn.flatten()).reshape(-1, 3, 3)
        norms  = np.cross(coords[:, 1] - coords[:, 0], coords[:, 2] - coords[:, 0])
        nmag   = np.linalg.norm(norms, axis=1)
        valid  = nmag > 1e-14
        if not np.any(valid):
            continue

        avg_n = np.mean(norms[valid] / nmag[valid, None], axis=0)
        avg_n /= np.linalg.norm(avg_n)
        d_val  = float(np.dot(avg_n, coords[valid].mean(axis=(0, 1))))

        plane  = openmc.Plane(
            a=float(avg_n[0]), b=float(avg_n[1]), c=float(avg_n[2]),
            d=d_val, boundary_type="reflective",
        )
        region = -plane if ori == _Orientation.FORWARD else +plane
        area   = float(0.5 * np.sum(nmag))
        rec    = {"id": sid, "area": area, "region": region,
                  "normal": avg_n, "d": d_val, "plane": plane}

        if abs(avg_n[2]) > 0.7:
            z_cands.append(rec)
        elif abs(avg_n[1]) > 0.7:
            y_cands.append(rec)

    poloidal = sorted(z_cands, key=lambda x: x["area"], reverse=True)[:2]
    toroidal = sorted(y_cands, key=lambda x: x["area"], reverse=True)[:2]
    final    = poloidal + toroidal

    for i, cut in enumerate(final):
        cut["plane"].id = base_surface_id + i

    return final


def _build_slab_geometry() -> openmc.Cell:
    r_in  = openmc.ZCylinder(r=cfg.SLAB_R_IN_CM,  boundary_type="vacuum")
    r_out = openmc.ZCylinder(r=cfg.SLAB_R_OUT_CM, boundary_type="vacuum")
    cuts  = _get_breeder_reflective_cuts(cfg.SLAB_TARGET_VOL_ID)
    region = +r_in & -r_out
    for cut in cuts:
        region &= cut["region"]
    return openmc.Cell(
        cell_id=10_000, region=region, fill=dagmc_universe, name="breeder_cell"
    )


if cfg.SIM_TYPE == "tokamak":
    sector_cell = _build_tokamak_geometry()
elif cfg.SIM_TYPE == "slab":
    sector_cell = _build_slab_geometry()
else:
    raise ValueError(
        f"Unknown SIM_TYPE: {cfg.SIM_TYPE!r}. Choose 'tokamak' or 'slab'."
    )

model.geometry = openmc.Geometry(root=[sector_cell])

# ──────────────────────────────────────────────────────────────────────────────
# MATERIALS  –  pulled directly from cfg + eudemo_materials.py
# ──────────────────────────────────────────────────────────────────────────────

# Per-breeder material attributes from cfg
armor_material:      openmc.Material           = cfg.armor_material
structural_material: openmc.Material           = cfg.structural_material
coolant_material:    openmc.Material           = cfg.coolant_material
breeder_material:    openmc.Material           = cfg.breeder_material
multiplier_material: Optional[openmc.Material] = cfg.multiplier_material
vv_material:         openmc.Material           = cfg.vv_material

structure_material_list: List[openmc.Material] = [
    structural_material, armor_material, vv_material
]

# Volume-fraction recipes from eudemo_materials.py
from eudemo_materials import build_mix_recipes, BreederType   

MIX_RECIPES_OBJ: dict = build_mix_recipes(cfg)

# Normalise volume fractions 
def _normalize_recipes(recipes: dict, tol: float = 1e-14) -> dict:
    out: dict = {}
    for name, comp in recipes.items():
        items = [(m, float(vf)) for m, vf in comp.items() if float(vf) > 0.0]
        s     = sum(vf for _, vf in items)
        if s <= tol:
            raise ValueError(f"Recipe '{name}' volume fractions sum to zero.")
        out[name] = {m: vf / s for m, vf in items}
    return out

MIX_RECIPES_NORM = _normalize_recipes(MIX_RECIPES_OBJ)

# Build mixed materials 
def _build_mixed_materials(recipes: dict) -> dict:
    mixed: dict = {}
    for name, comp in recipes.items():
        m = openmc.Material.mix_materials(
            list(comp.keys()), list(comp.values()), percent_type="vo", name=name
        )
        mixed[name] = m
    return mixed

mixed_materials    = _build_mixed_materials(MIX_RECIPES_NORM)
model.materials    = openmc.Materials(list(mixed_materials.values()))

# ──────────────────────────────────────────────────────────────────────────────
# CHUNK / RADIAL BIN HELPERS
# ──────────────────────────────────────────────────────────────────────────────
_OB_KEY_RE = re.compile(r"^OB_(\d+)_b(\d+)$")
_IB_KEY_RE = re.compile(r"^IB_(\d+)_b(\d+)$")

def parse_key(key: str) -> Tuple[str, int, int]:
    key = key.strip()
    m   = _OB_KEY_RE.match(key)
    if m:
        return "OB", int(m.group(1)), int(m.group(2))
    m = _IB_KEY_RE.match(key)
    if m:
        return "IB", int(m.group(1)), int(m.group(2))
    raise ValueError(f"Not a valid chunk key: {key!r} (expected 'OB_<r>_b<b>' or 'IB_<r>_b<b>')")

def _interp_three_point(t: float, start: float, mid: float, end: float) -> float:
    """Piecewise-linear through (start -> mid -> end) with mid at t=0.5, t in [0,1]."""
    t = float(t)
    if t <= 0.5:
        a = t / 0.5
        return float(start + a * (mid - start))
    a = (t - 0.5) / 0.5
    return float(mid + a * (end - mid))

def _bins_from_thicknesses(
    thicknesses: np.ndarray, *, start_cm: float = 0.0
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    w = np.asarray(thicknesses, float).ravel()
    if np.any(w <= 0.0):
        raise ValueError(f"Non-positive layer thickness: {w}")
    edges     = float(start_cm) + np.concatenate(([0.0], np.cumsum(w)))
    centroids = 0.5 * (edges[:-1] + edges[1:])
    return centroids, w, edges

def _known_thicknesses(geom: dict, *, side: str) -> np.ndarray:
    """
    Returns thicknesses for the "known" part BEFORE the last breeder thickness:
      [Armor, FW, breeder_1..breeder_(n_layers-1)]

    offsets_cm are assumed measured from start of breeder stack (after Armor+FW),
    so we do NOT include Armor/FW inside offsets; we prepend them explicitly.
    """
    s        = side.lower()
    armor_cm = float(geom["armor_cm"] - geom["fw_cm"])
    fw_cm    = float(geom["fw_cm"])
    n_layers = int(geom[f"{s}_n_layers"])
    offsets  = np.asarray(geom[f"{s}_offsets_cm"], float).ravel()

    need = n_layers - 1
    if offsets.size != need:
        raise ValueError(f"{side}: expected {need} offsets, got {offsets.size}")

    breeder    = np.empty(need, float)
    breeder[0] = abs(offsets[0])
    if need > 1:
        breeder[1:] = np.abs(np.diff(offsets))

    return np.concatenate(([armor_cm, fw_cm], breeder))

def _breeder_t(geom: dict, *, b: int) -> float:
    """
    Map breeder index b=1..n_breeder to a poloidal parameter t in [0,1].
    """
    n_breeder = int(geom["n_breeder"])
    if not (1 <= b <= n_breeder):
        raise ValueError(f"b={b} out of range 1..{n_breeder}")

    fr = geom.get("fractions", None)
    if fr is not None:
        fr = np.asarray(fr, dtype=float).ravel()
        if fr.size != max(n_breeder - 1, 0):
            raise ValueError(f"fractions length {fr.size} != n_breeder-1 ({n_breeder-1})")
        if np.any(fr <= 0.0) or np.any(fr >= 1.0) or np.any(np.diff(fr) <= 0.0):
            raise ValueError(f"fractions must be strictly increasing in (0,1): {fr}")

        bp = np.concatenate(([0.0], fr, [1.0]))
        return float(0.5 * (bp[b - 1] + bp[b]))

    if n_breeder < 2:
        return 0.5
    return float((b - 1) / (n_breeder - 1))

def _make_radial_bins(
    *, key: str, geom: dict, gap_cm: float, start_cm: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    OB_LAST = (87.62, 96.825, 89.67)
    OB_VV   = (76.5,  110.0,  84.63)
    IB_LAST = (84.11, 75.5,   79.65)
    IB_VV   = (68.33, 60.0,   60.0)

    side, _, b = parse_key(key)
    s_l        = side.lower()
    n_layers   = int(geom[f"{s_l}_n_layers"])

    known = _known_thicknesses(geom, side=s_l)
    t     = _breeder_t(geom, b=b)

    tL, mL, bL = OB_LAST if side == "OB" else IB_LAST
    tV, mV, bV = OB_VV   if side == "OB" else IB_VV

    interp_last = _interp_three_point(t, tL, mL, bL)
    vv          = _interp_three_point(t, tV, mV, bV)

    offsets      = np.asarray(geom[f"{s_l}_offsets_cm"], float).ravel()
    last_offset  = abs(offsets[-1]) + float(geom["armor_cm"])
    last_breeder = interp_last - last_offset

    if last_breeder <= 0.0:
        raise ValueError(
            f"{key}: last_breeder={last_breeder:.4g} cm ≤ 0 "
            f"(interp_last={interp_last:.4g}, last_offset={last_offset:.4g})"
        )

    thicknesses = np.concatenate([known, [last_breeder, gap_cm + vv]])
    expected    = n_layers + 3
    if thicknesses.size != expected:
        raise RuntimeError(
            f"{key}: {thicknesses.size} thickness bins built, expected {expected}"
        )

    return _bins_from_thicknesses(thicknesses, start_cm=start_cm)

def build_breeder_chunks(
    INPUT_JSON: Path,
    *,
    default_chunk_key: str   = "OB_1_b6",
    gap_cm:            float = cfg.CHUNK_GAP_CM,
    start_cm:          float = cfg.CHUNK_START_CM,
) -> dict:
    """
    Parse the geometry JSON and return a dict with all chunk/bin helpers.

    """
    if not INPUT_JSON.is_file():
        raise FileNotFoundError(f"Geometry JSON not found: {INPUT_JSON}")

    with open(INPUT_JSON, encoding="utf-8") as f:
        data = json.load(f)

    inv  = data["inventory"]
    geom = data["geometry"]

    BLANKET_N   = int(inv["BLANKET_FINAL"])
    ob_n_layers = int(geom["ob_n_layers"])
    ib_n_layers = int(geom["ib_n_layers"])
    n_breeder   = int(geom["n_breeder"])
    ob_regions  = int(geom["ob_regions"])
    ib_regions  = int(geom["ib_regions"])

    OB_CHUNK_SIZE = ob_n_layers + 3
    IB_CHUNK_SIZE = ib_n_layers + 3

    cell_ids_all = list(range(1, BLANKET_N + 1))
    OB_TOTAL     = OB_CHUNK_SIZE * n_breeder * ob_regions

    if OB_TOTAL > BLANKET_N:
        raise ValueError(f"OB_TOTAL={OB_TOTAL} > BLANKET_N={BLANKET_N}")

    ob_ids = cell_ids_all[:OB_TOTAL]
    ib_ids = cell_ids_all[OB_TOTAL:]

    ob_chunks = [ob_ids[i:i + OB_CHUNK_SIZE] for i in range(0, len(ob_ids), OB_CHUNK_SIZE)]
    ib_chunks = [ib_ids[i:i + IB_CHUNK_SIZE] for i in range(0, len(ib_ids), IB_CHUNK_SIZE)]

    ob_by_key: Dict[str, List[int]] = {}
    k = 0
    for r in range(1, ob_regions + 1):
        for b in range(1, n_breeder + 1):
            ob_by_key[f"OB_{r}_b{b}"] = ob_chunks[k]; k += 1

    ib_by_key: Dict[str, List[int]] = {}
    k = 0
    for r in range(1, ib_regions + 1):
        for b in range(1, n_breeder + 1):
            ib_by_key[f"IB_{r}_b{b}"] = ib_chunks[k]; k += 1

    ALL_KEYS = list(ob_by_key) + list(ib_by_key)

    armor_cell_ids: List[int] = []
    for chunk in ob_by_key.values():
        if len(chunk) > 0:
            armor_cell_ids.append(int(chunk[0]))  # armor is always the first cell
    for chunk in ib_by_key.values():
        if len(chunk) > 0:
            armor_cell_ids.append(int(chunk[0]))

    def cell_ids_for_key(key: str) -> List[int]:
        side, _, _ = parse_key(key)
        return ob_by_key[key] if side == "OB" else ib_by_key[key]

    def radial_bins_for_key(key: str):
        return _make_radial_bins(key=key, geom=geom, gap_cm=gap_cm, start_cm=start_cm)

    if default_chunk_key in ob_by_key:
        selected = ob_by_key[default_chunk_key]
    elif default_chunk_key in ib_by_key:
        selected = ib_by_key[default_chunk_key]
    else:
        raise KeyError(
            f"default_chunk_key={default_chunk_key!r} not found. "
            f"Valid keys (sample): {ALL_KEYS[:6]} ..."
        )

    return {
        "data":                    data,
        "geom":                    geom,
        "cell_ids_all":            cell_ids_all,
        "ob_by_key":               ob_by_key,
        "ib_by_key":               ib_by_key,
        "armor_cell_ids":          armor_cell_ids,
        "ALL_KEYS":                ALL_KEYS,
        "OB_CHUNK_SIZE":           OB_CHUNK_SIZE,
        "IB_CHUNK_SIZE":           IB_CHUNK_SIZE,
        "n_breeder":               n_breeder,
        "selected_chunk_key":      default_chunk_key,
        "selected_chunk_cell_ids": selected,
        "cell_ids_for_key":        cell_ids_for_key,
        "radial_bins_for_key":     radial_bins_for_key,
    }

# Build chunks at import time using the INPUT_JSON from inputs.py
_chunk_result = build_breeder_chunks(
    cfg.INPUT_JSON, default_chunk_key=cfg.ALBEDO_CHUNK_KEY
)
ob_by_key     = _chunk_result["ob_by_key"]
ib_by_key     = _chunk_result["ib_by_key"]
n_breeder     = _chunk_result["n_breeder"]
OB_CHUNK_SIZE = _chunk_result["OB_CHUNK_SIZE"]
IB_CHUNK_SIZE = _chunk_result["IB_CHUNK_SIZE"]
cell_ids_all  = _chunk_result["cell_ids_all"]

# Surface info for the albedo / current-tally chunk
_albedo_cells = _chunk_result["selected_chunk_cell_ids"]

info, all_surface_ids, external_surface_ids, internal_surface_ids = (
    dagmc_volume_surface_info(pydagmc_model, _albedo_cells)
)

