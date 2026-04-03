#!/usr/bin/env python3

"""
Slab model based on Bluemira/Tokamak wedge dagmc model
"""

from __future__ import annotations

# built-in modules
from enum import Enum
from typing import Dict, Tuple, List, Iterable, Set, Optional
import math
import re
import json
from pathlib import Path

import openmc
import numpy as np
import h5py
from openmc_plasma_source import tokamak_source
from pymoab import core, types
import pydagmc

import os
import sys

# =============================================================================
# PATH SETUP
# =============================================================================

SCRIPT_DIR = Path.cwd()        # cortex/slab
PROJECT_ROOT = SCRIPT_DIR.parent                    # cortex
BLUEMIRA_DIR = PROJECT_ROOT / "detailed_bluemira_hcll"  # cortex/detailed_bluemira

SLAB_RUN_DIR = (SCRIPT_DIR / "neutronics_run")
SLAB_RUN_DIR.mkdir(parents=True, exist_ok=True)

SLAB_MODEL_XML = SLAB_RUN_DIR / "model.xml"

# Add materials folder (if it lives in cortex/materials)
module_path = (PROJECT_ROOT / "materials")
sys.path.append(str(module_path))
import materials

# =============================================================================
# USER INPUTS
# =============================================================================

_DAGMC_MODEL_FILE = (BLUEMIRA_DIR / "eudemo_hcll.h5m")
INPUT_JSON = (BLUEMIRA_DIR / "EUDEMO_HCLL_inputs.json")

# Fixed-source run configuration
# Read the surface source written by the detailed_bluemira neutronics run
SURF_SOURCE_FILE = (BLUEMIRA_DIR / "neutronics_run"  / "surface_source.h5")

TARGET_VOL_ID = 78

# Radial vacuum cylinders (cm)
R_IN_CM = 1140.0
R_OUT_CM = 1375.0

# -----------------------------------------------------------------------------
# SURFACE INFO HELPERS
# -----------------------------------------------------------------------------
class Orientation(Enum):
    FORWARD = 1
    REVERSE = -1

def dagmc_bounding_box(pydagmc_model, volume_id):
    volume = pydagmc_model.volumes_by_id[volume_id]
    triangle_coords = volume.triangle_coords
    min_coords = np.min(triangle_coords, axis=0)
    max_coords = np.max(triangle_coords, axis=0)
    return openmc.BoundingBox(min_coords, max_coords)

def dagmc_volume_surface_info(pydagmc_model, volume_ids):
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
                orientation = Orientation.FORWARD
            else:
                orientation = Orientation.REVERSE

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

def get_breeder_reflective_cuts(
    mb: core.Core,
    pydagmc_model: pydagmc.Model,
    volume_id: int,
    plane_surface_ids=None,   # <-- NEW: list/tuple of 4 ints, or None
    base_surface_id=50_000    # <-- NEW: used if plane_surface_ids is None
    ):
    """
    Finds 2 Z-dominant and 2 Y-dominant reflective planes for a given volume.
    Returns list of dicts: {id, area, region, normal, d, plane}.
    Also assigns OpenMC surface IDs to the 4 cut planes (reflective).
    """
    surface_info_map, _, _, _ = dagmc_volume_surface_info(pydagmc_model, [volume_id])
    volume_data = surface_info_map[int(volume_id)]

    id_tag = mb.tag_get_handle("GLOBAL_ID")
    dim_tag = mb.tag_get_handle("GEOM_DIMENSION")
    s_sets = mb.get_entities_by_type_and_tag(0, types.MBENTITYSET, dim_tag, [2])

    z_candidates = []
    y_candidates = []

    for surf_entry in volume_data["all_surfaces"]:
        sid = int(surf_entry["surface_id"])
        orientation = surf_entry["orientation"]

        target_s_list = [s for s in s_sets if int(mb.tag_get_data(id_tag, s)[0][0]) == sid]
        if not target_s_list:
            continue
        target_s = target_s_list[0]

        s_entities = mb.get_entities_by_handle(target_s, True)
        tris = [e for e in s_entities if mb.type_from_handle(e) == types.MBTRI]
        if not tris:
            continue

        conn = mb.get_connectivity(tris).reshape(-1, 3)
        coords = mb.get_coords(conn.flatten()).reshape(-1, 3, 3)
        v1 = coords[:, 1] - coords[:, 0]
        v2 = coords[:, 2] - coords[:, 0]
        normals = np.cross(v1, v2)
        norm_vals = np.linalg.norm(normals, axis=1)

        valid = norm_vals > 1e-14
        if not np.any(valid):
            continue

        avg_n = np.mean(normals[valid] / norm_vals[valid][:, np.newaxis], axis=0)
        avg_n /= np.linalg.norm(avg_n)

        surface_centroid = coords[valid].mean(axis=(0, 1))
        d_openmc = float(np.dot(avg_n, surface_centroid))

        # create the plane (ID assigned later, after we know the final 4)
        p = openmc.Plane(
            a=float(avg_n[0]), b=float(avg_n[1]), c=float(avg_n[2]),
            d=d_openmc, boundary_type="reflective"
        )

        region = -p if orientation == Orientation.FORWARD else +p

        area = float(0.5 * np.sum(norm_vals))
        entry = {
            "id": sid,
            "area": area,
            "region": region,
            "normal": avg_n,
            "d": d_openmc,
            "plane": p,  # <-- NEW: keep a direct handle to the plane
        }

        if abs(avg_n[2]) > 0.7:
            z_candidates.append(entry)
        elif abs(avg_n[1]) > 0.7:
            y_candidates.append(entry)

    poloidal_cuts = sorted(z_candidates, key=lambda x: x["area"], reverse=True)[:2]
    toroidal_cuts = sorted(y_candidates, key=lambda x: x["area"], reverse=True)[:2]
    final_cuts = poloidal_cuts + toroidal_cuts

    # ---- NEW: assign OpenMC surface IDs to the *final 4* planes ----
    if plane_surface_ids is None:
        plane_surface_ids = [base_surface_id + i for i in range(len(final_cuts))]

    if len(plane_surface_ids) != len(final_cuts):
        raise ValueError(
            f"plane_surface_ids must have length {len(final_cuts)}, got {len(plane_surface_ids)}"
        )

    for cut, new_sid in zip(final_cuts, plane_surface_ids):
        cut["plane"].id = int(new_sid)

    return final_cuts

# =============================================================================
# INITIALIZE MODEL
# =============================================================================
model = openmc.Model()

# =============================================================================
# GEOMETRY
# =============================================================================
dagmc_universe = openmc.DAGMCUniverse(filename=str(_DAGMC_MODEL_FILE))
pydagmc_model = pydagmc.Model(str(_DAGMC_MODEL_FILE))

mb = core.Core()
mb.load_file(str(_DAGMC_MODEL_FILE))

# reserve IDs
openmc.reserve_ids([v.id for v in pydagmc_model.volumes], cls=openmc.Cell)
openmc.reserve_ids([s.id for s in pydagmc_model.surfaces], cls=openmc.Surface)

r_out = openmc.ZCylinder(r=float(R_OUT_CM), boundary_type="vacuum")
r_in = openmc.ZCylinder(r=float(R_IN_CM), boundary_type="vacuum")

final_cuts = get_breeder_reflective_cuts(mb, pydagmc_model, TARGET_VOL_ID,
                                        base_surface_id=50_000) # or  plane_surface_ids=[9001, 9002, 9003, 9004]
if len(final_cuts) != 4:
    print(f"[warn] expected 4 cuts, got {len(final_cuts)} (continuing)")

sector_region = (+r_in & -r_out)
for cut in final_cuts:
    sector_region &= cut["region"]

sector_cell = openmc.Cell(
    cell_id=10_000,
    fill=dagmc_universe,
    region=sector_region,
    name="breeder_cell",
)
model.geometry = openmc.Geometry([sector_cell])

# -----------------------------------------------------------------------------
# MATERIALS
# -----------------------------------------------------------------------------
# TODO: add citations for where these densities come from
# TODO: need to review materials.py for correctness for all materials
ccz     = materials.CuCrZr(8.9)
h       = materials.Helium(0.0001785)
nb3sn   = materials.Nb3Sn(5.7)
epoxy   = materials.Epoxy(1.207)
bronze  = materials.Bronze(8.8775)
nbti    = materials.NbTi(6.538)
c       = materials.Cu(8.96)

# tungsten, density based on PNNL material compendium value
w       = materials.W(19.3)

ss304_b4 = materials.ss304_b4(7.8)
ss316   = materials.ss316Ln_ig(7.93)
# Plasma Region
plasma = openmc.Material()
plasma.set_density("g/cm3", 1e-6)
plasma.add_element("H", 1.0)

# ----------------------------------------------
# CASE MATERIALS
# USER STRUCTURAL OR ARMOR MATERIAL + COOLANT
# ----------------------------------------------
# Armor (default tungsten if not provided by user)
armor_material = materials.W(19.3)

# Structural material (Eurofer if not provided by user)
structural_material = materials.eurofer97(7.87)

# Coolant material (dependent on the breeder type) (water for WCLL)
coolant_material = materials.Helium(0.00498)  # ~8 MPa and 500 C

# Breeder material (dependent on the breeder type) (PbLi for WCLL)
breeder_material = materials.PbLi(0.90, 9.8)

# -----------------------------------------------------------------------------
# Structural / breeder / coolant definitions
# -----------------------------------------------------------------------------
# TODO: very unclear what these lists are. Why are the other solid materials not included in structure_material_list?
structure_material_list: list[openmc.Material] = [structural_material, armor_material, ss316]
breeder_material_list: list[openmc.Material] = [breeder_material]
coolant_material_list: list[openmc.Material] = [coolant_material]

# -----------------------------------------------------------------------------
# MIX RECIPES
# (easier to find nuclides ratios with this)
# (ADD any openmc.Materials.mix_materials() to this dict)
# TODO: add citation for where these numbers come from
# -----------------------------------------------------------------------------
MIX_RECIPES_OBJ: dict[str, dict[openmc.Material, float]] = {
    # Armor
    "Armor": {armor_material: 1.00},
    # FW
    "First_Wall":  {coolant_material: 0.049, structural_material: 0.951},
    # IB
    "ib_layer_1":  {breeder_material: 0.883, coolant_material: 0.002, structural_material: 0.115},
    "ib_layer_2":  {breeder_material: 0.878, coolant_material: 0.003, structural_material: 0.119},
    "ib_layer_3":  {breeder_material: 0.878, coolant_material: 0.003, structural_material: 0.119},
    "ib_layer_4":  {breeder_material: 0.313, coolant_material: 0.118, structural_material: 0.569},
    "ib_layer_5":  {breeder_material: 0.901, coolant_material: 0.009, structural_material: 0.090},
    "ib_layer_6":  {breeder_material: 0.104, coolant_material: 0.058, structural_material: 0.838},
    "ib_layer_7":  {breeder_material: 0.010, coolant_material: 0.861, structural_material: 0.129},
    "ib_layer_8":  {breeder_material: 0.010, coolant_material: 0.876, structural_material: 0.114},
    "ib_layer_9":  {breeder_material: 0.010, coolant_material: 0.011, structural_material: 0.979},
    "ib_layer_10": {breeder_material: 0.010, coolant_material: 0.011, structural_material: 0.979},    
    # OB
    "ob_layer_1":  {breeder_material: 0.883, coolant_material: 0.002, structural_material: 0.115},
    "ob_layer_2":  {breeder_material: 0.878, coolant_material: 0.003, structural_material: 0.119},
    "ob_layer_3":  {breeder_material: 0.878, coolant_material: 0.003, structural_material: 0.119},
    "ob_layer_4":  {breeder_material: 0.313, coolant_material: 0.118, structural_material: 0.569},
    "ob_layer_5":  {breeder_material: 0.901, coolant_material: 0.009, structural_material: 0.090},
    "ob_layer_6":  {breeder_material: 0.104, coolant_material: 0.058, structural_material: 0.838},
    "ob_layer_7":  {breeder_material: 0.010, coolant_material: 0.861, structural_material: 0.129},
    "ob_layer_8":  {breeder_material: 0.010, coolant_material: 0.876, structural_material: 0.114},
    "ob_layer_9":  {breeder_material: 0.010, coolant_material: 0.011, structural_material: 0.979},
    "ob_layer_10": {breeder_material: 0.010, coolant_material: 0.011, structural_material: 0.979},    
    # Divertor
    "Divertor":    {ccz: 0.00552, c: 0.00438, structural_material: 0.5238, armor_material: 0.01026, coolant_material: 0.45604},
    # VV
    "VV_IB":       {ss316: 0.6, coolant_material: 0.4},
    "VV_OB":       {ss316: 0.6, coolant_material: 0.4},
    "VV_ports_all":{ss316: 0.6, coolant_material: 0.4},
    "VV_port_Fill":{ss316: 0.6, coolant_material: 0.4},
    # Coils
    "PC_PFC":      {nbti: 0.02895, c: 0.1169, epoxy: 0.18, bronze: 0.0735, h: 0.1682, ss316: 0.43245},
    "TFcoil":      {nb3sn: 0.02895, c: 0.1169, epoxy: 0.18, bronze: 0.0735, h: 0.1682, ss316: 0.43245},
    # Plasma Region
    "Plasma_Region": {plasma: 1.00},
    # Cryostat
    "Cryostat": {ss316: 1.00},
    # Shielding
    "RadiationShield_all": {ss304_b4: 1.00},
}

# -----------------------------------------------------------------------------
# Build mixed materials and assign model.materials (VO ONLY)
# -----------------------------------------------------------------------------

def normalize_mix_recipes_obj(
    recipes_obj: dict[str, dict[openmc.Material, float]],
    *,
    drop_nonpositive: bool = True,
    tol: float = 1e-14,
    ) -> dict[str, dict[openmc.Material, float]]:
    """
    Normalize volume % of each materials components in the mixture of [mix_recipes_obj]

    - If drop_nonpositive: remove components with vf <= 0 before normalizing.
    - Raises if the sum is <= 0 after filtering.
    """
    out: dict[str, dict[openmc.Material, float]] = {}

    for name, comp in recipes_obj.items():
        if comp is None or len(comp) == 0:
            raise ValueError(f"[MIX_RECIPES_OBJ] Recipe '{name}' is empty.")

        items = []
        for m, vf in comp.items():
            vf = float(vf)
            if drop_nonpositive and vf <= 0.0:
                continue
            items.append((m, vf))

        s = sum(vf for _, vf in items)
        if s <= tol:
            raise ValueError(f"[MIX_RECIPES_OBJ] Recipe '{name}' fractions sum <= 0 (sum={s}).")

        out[name] = {m: vf / s for m, vf in items}

    return out


def build_and_set_model_materials_from_obj_recipes_vo(
    model: openmc.Model,
    *,
    recipes_obj: dict[str, dict[openmc.Material, float]],
    ) -> dict[str, openmc.Material]:
    """
    VO-only builder: interprets recipe fractions as volume fractions.

    This routine returns a dictionary of the material names, followed by the
    openmc.Material definition for each material in the geometry.
    """
    mixed: dict[str, openmc.Material] = {}

    for name, comp in recipes_obj.items():
        mats = list(comp.keys())
        fracs = [v for v in comp.values()]

        m = openmc.Material.mix_materials(mats, fracs, percent_type="vo", name=name)
        mixed[name] = m

    model.materials = openmc.Materials(list(mixed.values()))
    return mixed

MIX_RECIPES_OBJ_NORM = normalize_mix_recipes_obj(MIX_RECIPES_OBJ)

mixed = build_and_set_model_materials_from_obj_recipes_vo(
    model,
    recipes_obj=MIX_RECIPES_OBJ_NORM,
)

# -----------------------------------------------------------------------------
# SETTINGS
# -----------------------------------------------------------------------------
model.settings = openmc.Settings()
model.settings.photon_transport = True
model.settings.run_mode = "fixed source"
model.settings.surf_source_read = {"path": str(SURF_SOURCE_FILE)}

# Set TALLY_CONVERGENCE_THRESHOLD to 0.01 (1%) or 0.001 (0.1%)
TALLY_CONVERGENCE_THRESHOLD = 0.01

model.settings.batches = 10           # minimum batches before triggers are checked
model.settings.trigger_active = True
model.settings.trigger_batch_interval = 10   # check triggers every N batches
model.settings.particles = 1_000_000
model.settings.trigger_max_batches = 2000     # hard ceiling

# -----------------------------------------------------------------------------
# DAGMC volume sync so cells have volumes
# -----------------------------------------------------------------------------
openmc.Cell.reset_ids()
openmc.Surface.reset_ids()

model.init_lib(output=False)
model.sync_dagmc_universes()
model.finalize_lib()

openmc.reserve_ids([c_id for c_id in model.geometry.get_all_cells()], cls=openmc.Cell)
openmc.reserve_ids([s_id for s_id in model.geometry.get_all_surfaces()], cls=openmc.Surface)

# Apply volumes/bounding boxes to dagmc universe cells
dagmc_universe_cells = dagmc_universe.get_all_cells()
for volume in pydagmc_model.volumes:
    dagmc_universe_cells[volume.id].volume = volume.volume
    dagmc_universe_cells[volume.id].bounding_box = dagmc_bounding_box(pydagmc_model, volume.id)

all_cells = model.geometry.get_all_cells()  

# -----------------------------------------------------------------------------
# IMPORT JSON GEOMETRY INFO (cell IDs)
# Also provide functions to create centroids and bin wedges
# Added in the neutronics_model.py to simplify other scripts
# -----------------------------------------------------------------------------
_OB_KEY_RE = re.compile(r"^OB_(\d+)_b(\d+)$")
_IB_KEY_RE = re.compile(r"^IB_(\d+)_b(\d+)$")

def parse_key(key: str) -> Tuple[str, int, int]:
    key = key.strip()
    m = _OB_KEY_RE.match(key)
    if m:
        return "OB", int(m.group(1)), int(m.group(2))
    m = _IB_KEY_RE.match(key)
    if m:
        return "IB", int(m.group(1)), int(m.group(2))
    raise ValueError(f"Not a chunk key: {key!r} (expected 'OB_<r>_b<b>' or 'IB_<r>_b<b>')")

def interp_three_point(t: float, start: float, mid: float, end: float) -> float:
    """Piecewise-linear through (start -> mid -> end) with mid at t=0.5, t in [0,1]."""
    t = float(t)
    if t <= 0.5:
        a = t / 0.5
        return float(start + a * (mid - start))
    a = (t - 0.5) / 0.5
    return float(mid + a * (end - mid))

def build_bins_from_thicknesses(thicknesses_cm: np.ndarray, *, start_cm: float = 0.0):
    w = np.asarray(thicknesses_cm, dtype=float).ravel()
    if np.any(w <= 0.0):
        raise ValueError(f"Non-positive thickness encountered: {w}")
    edges = float(start_cm) + np.concatenate(([0.0], np.cumsum(w)))
    centroids = 0.5 * (edges[:-1] + edges[1:])
    return centroids, w, edges

def build_known_thicknesses_from_offsets(geom: dict, *, side: str) -> np.ndarray:
    """
    Returns thicknesses for the "known" part BEFORE the last breeder thickness:
      [Armor, FW, breeder_1..breeder_(n_layers-1)]

    offsets_cm are assumed measured from start of breeder stack (after Armor+FW),
    so we do NOT include Armor/FW inside offsets; we prepend them explicitly.
    """
    side = side.lower()
    armor_cm = float(geom["armor_cm"] - geom["fw_cm"])
    fw_cm = float(geom["fw_cm"])
    n_layers = int(geom[f"{side}_n_layers"])
    offsets = np.asarray(geom[f"{side}_offsets_cm"], dtype=float).ravel()

    if armor_cm <= 0.0:
        raise ValueError(f"{side}: armor_cm must be positive (got {armor_cm})")
    if fw_cm <= 0.0:
        raise ValueError(f"{side}: fw_cm must be positive (got {fw_cm})")

    need = n_layers - 1
    if offsets.size != need:
        raise ValueError(f"{side}: expected {need} offsets for n_layers={n_layers}, got {offsets.size}")

    # breeder thicknesses from offsets: [|o0|, |o1-o0|, ...]
    breeder = np.empty(need, dtype=float)
    breeder[0] = abs(offsets[0])
    if need > 1:
        breeder[1:] = np.abs(np.diff(offsets))
    if np.any(breeder <= 0.0):
        raise ValueError(f"{side}: non-positive breeder thicknesses from offsets: {breeder}")

    return np.concatenate(([armor_cm, fw_cm], breeder))

def breeder_t_from_fractions(geom: dict, *, b: int) -> float:
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

def make_radial_bins_for_key(
    *,
    key: str,
    geom: dict,
    gap_cm: float = 2.0,
    start_cm: float = 0.0,
    ):
    """
    Build centroids/widths/edges for either OB_* or IB_* key.

    thicknesses are:
      [Armor, FW, breeder_1..breeder_(n_layers-1), breeder_last, (gap + vv)]

    Correction:
      - offsets are from breeder start -> last_offset must include (armor_cm + fw_cm)
      - t is based on b using fractions (poloidal discretization)
    """
    OB_LAST_LAYER_TOP_MID_BOT = (87.62, 96.825, 89.67)
    OB_VV_TOP_MID_BOT         = (76.5, 110.0, 84.63)

    IB_LAST_LAYER_TOP_MID_BOT = (84.11, 75.5, 79.65)
    IB_VV_TOP_MID_BOT         = (68.33, 60.0, 60.0)
    side, _, b = parse_key(key)
    side_l = side.lower()

    n_layers = int(geom[f"{side_l}_n_layers"])
    offsets = np.asarray(geom[f"{side_l}_offsets_cm"], dtype=float).ravel()
    if offsets.size == 0:
        raise ValueError(f"{key}: {side_l}_offsets_cm is empty.")

    known = build_known_thicknesses_from_offsets(geom, side=side_l)

    t = breeder_t_from_fractions(geom, b=b)

    if side == "OB":
        topL, midL, botL = OB_LAST_LAYER_TOP_MID_BOT
        topV, midV, botV = OB_VV_TOP_MID_BOT
        interp_last = interp_three_point(t, topL, midL, botL)
        vv = interp_three_point(t, topV, midV, botV)
    else:
        topL, midL, botL = IB_LAST_LAYER_TOP_MID_BOT
        topV, midV, botV = IB_VV_TOP_MID_BOT
        interp_last = interp_three_point(t, topL, midL, botL)
        vv = interp_three_point(t, topV, midV, botV)

    pre_breeder_shift = float(geom["armor_cm"]) # geom["armor_cm"] => fw + armor offset
    last_offset = float(abs(offsets[-1]) + pre_breeder_shift)

    last_breeder = float(interp_last - last_offset)
    if last_breeder <= 0.0:
        raise ValueError(
            f"{key}: last breeder thickness <= 0: "
            f"interp_last={interp_last:.6g} cm, last_offset={last_offset:.6g} cm "
            f"(t={t:.6f}, b={b}/{int(geom['n_breeder'])})"
        )

    thicknesses = np.concatenate([known, [last_breeder, float(gap_cm) + float(vv)]])

    expected = n_layers + 3
    if thicknesses.size != expected:
        raise RuntimeError(f"{key}: thickness bins={thicknesses.size}, expected {expected}")

    return build_bins_from_thicknesses(thicknesses, start_cm=start_cm)

def build_breeder_chunks(
    INPUT_JSON: Path,
    *,
    default_equatorial_ob_key: str = "OB_1_b6",
    gap_cm: float = 2.0,
    start_cm: float = 0.0,
    ) -> dict:
    with INPUT_JSON.open("r", encoding="utf-8") as f:
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

    OB_TOTAL = OB_CHUNK_SIZE * n_breeder * ob_regions
    if OB_TOTAL > BLANKET_N:
        raise ValueError(f"OB_TOTAL={OB_TOTAL} exceeds BLANKET_N={BLANKET_N}")

    ob_ids = cell_ids_all[:OB_TOTAL]
    ib_ids = cell_ids_all[OB_TOTAL:]

    ob_chunks = [ob_ids[i:i + OB_CHUNK_SIZE] for i in range(0, len(ob_ids), OB_CHUNK_SIZE)]
    ib_chunks = [ib_ids[i:i + IB_CHUNK_SIZE] for i in range(0, len(ib_ids), IB_CHUNK_SIZE)]

    ob_by_key: Dict[str, List[int]] = {}
    k = 0
    for r in range(1, ob_regions + 1):
        for b in range(1, n_breeder + 1):
            ob_by_key[f"OB_{r}_b{b}"] = ob_chunks[k]
            k += 1

    ib_by_key: Dict[str, List[int]] = {}
    k = 0
    for r in range(1, ib_regions + 1):
        for b in range(1, n_breeder + 1):
            ib_by_key[f"IB_{r}_b{b}"] = ib_chunks[k]
            k += 1

    ALL_KEYS = list(ob_by_key.keys()) + list(ib_by_key.keys())

    equatorial_ob_cell_ids = ob_by_key.get(default_equatorial_ob_key, [])

    # convenience callables
    def cell_ids_for_key(key: str) -> List[int]:
        side, _, _ = parse_key(key)
        return ob_by_key[key] if side == "OB" else ib_by_key[key]

    def radial_bins_for_key(key: str):
        return make_radial_bins_for_key(key=key, geom=geom, gap_cm=gap_cm, start_cm=start_cm)

    return {
        "data": data,
        "geom": geom,
        "cell_ids_all": cell_ids_all,
        "ob_by_key": ob_by_key,
        "ib_by_key": ib_by_key,
        "ALL_KEYS": ALL_KEYS,
        "OB_CHUNK_SIZE": OB_CHUNK_SIZE,
        "IB_CHUNK_SIZE": IB_CHUNK_SIZE,
        "equatorial_ob_key": default_equatorial_ob_key,
        "equatorial_ob_cell_ids": equatorial_ob_cell_ids,
        "cell_ids_for_key": cell_ids_for_key,
        "radial_bins_for_key": radial_bins_for_key,
    }

_cells_chunk = build_breeder_chunks(INPUT_JSON, default_equatorial_ob_key="OB_1_b6")

cell_ids = _cells_chunk["cell_ids_for_key"]("OB_1_b6")
cell_ids_equatorial_ob = _cells_chunk["equatorial_ob_cell_ids"]

info, all_surface_ids, external_surface_ids, internal_surface_ids = dagmc_volume_surface_info(
    pydagmc_model,
    cell_ids
)
OB_KEY = "OB_1_b6"

#present = sorted(int(cid) for cid in dagmc_universe.get_all_cells().keys())
#print(f"[present] dagmc_universe cells present: {len(present)}")
#print(f"[present] min/max: {present[0]}..{present[-1]}")
#print(f"[present] sample: {present[:20]}")

ob_cells = _cells_chunk["cell_ids_for_key"](OB_KEY)
#ob_present = [cid for cid in ob_cells if cid in present]
#print(f"[chunk] {OB_KEY}: {len(ob_cells)} from JSON, {len(ob_present)} present in wrapper")
#print(f"[chunk] present ids: {ob_present}")
# -----------------------------------------------------------------------------
# TALLIES
# -----------------------------------------------------------------------------

# Source intensity -> scaling
total_power = 2e9  # 2000 MW
number_sectors = 16
section_power = total_power / number_sectors
ev_to_joule = 1.60218e-19
ev_fusion = 17.6e6
convert_e = ev_to_joule * ev_fusion
neutron_source_rate = section_power / convert_e
s_in_y = (365 * 24 * 60 * 60)

# Check surface_source.h5 file:
if not SURF_SOURCE_FILE.is_file():
    raise FileNotFoundError(f"Surface source file not found: {SURF_SOURCE_FILE}")

with h5py.File(str(SURF_SOURCE_FILE), "r") as f:
    bank = f["source_bank"]

    particles = bank["particle"][:]   # particle type array

total = len(particles)

n_neutrons = np.sum(particles == 0)
n_photons  = np.sum(particles == 1)

ratio_neutrons =  n_neutrons / total if total > 0 else 0.0
ratio_photons  =  n_photons  / total if total > 0 else 0.0

neutron_ratio_source = (1.0 / ratio_neutrons)

pct_neutrons = 100 * ratio_neutrons
pct_photons = 100 * ratio_photons

print("\n--- Surface Source Particle Composition ---\n")

print(f"Source file : {SURF_SOURCE_FILE}")
print(f"Total particles : {total:,}\n")

print(f"Neutrons : {n_neutrons:,}  ({pct_neutrons:.2f} %)")
print(f"Photons  : {n_photons:,}  ({pct_photons:.2f} %)\n")
print(f"Neutron ratio source scaling: {neutron_ratio_source:,}")

print("Sanity check:")
print(f"Sum = {n_neutrons + n_photons:,}")
# ------------------------------------------

# Filters
cell_filter = openmc.CellFilter(cell_ids)
particle_filter = openmc.ParticleFilter(bins=["neutron", "photon"])
t_surf_filter = openmc.SurfaceFilter(external_surface_ids)
n_particle_filter = openmc.ParticleFilter(bins=["neutron"])
energies = openmc.mgxs.GROUP_STRUCTURES["CCFE-709"]
energy_filter = openmc.EnergyFilter(energies)

model.tallies = openmc.Tallies()

# Neutron and photon flux
flux_tally = openmc.Tally()
flux_tally.filters = [cell_filter, particle_filter, energy_filter]
flux_tally.scores = ["flux"]
model.tallies.append(flux_tally)

# Adding total flux_total_tally (in OB_1_b6) for trigger only
flux_tally_total = openmc.Tally()
flux_tally_total.filters = [cell_filter, n_particle_filter]
flux_tally_total.scores = ["flux"]
flux_tally_total.triggers = [
    openmc.Trigger(trigger_type="rel_err", threshold=TALLY_CONVERGENCE_THRESHOLD)
]
model.tallies.append(flux_tally_total)


# Turn it off current tallies for albedo in the slab model
DO_ALBEDO = False
if DO_ALBEDO:
    t_current_tally = openmc.Tally()
    t_current_tally.filters = [t_surf_filter, particle_filter]
    t_current_tally.scores = ["current"]
    model.tallies.append(t_current_tally)

    p_current_tallies: dict[int, openmc.Tally] = {}
    for cid in cell_ids_equatorial_ob:
        ocell = dagmc_universe_cells[cid]
        surf_ids_for_cell = [int(s["surface_id"]) for s in info.get(cid, {}).get("all_surfaces", [])]
        if not surf_ids_for_cell:
            continue

        cell_from_filter = openmc.CellFromFilter([ocell])
        surf_filter = openmc.SurfaceFilter(surf_ids_for_cell)

        p_current_tally = openmc.Tally()
        p_current_tally.filters = [cell_from_filter, surf_filter, particle_filter]
        p_current_tally.scores = ["current"]

        model.tallies.append(p_current_tally)
        p_current_tallies[cid] = p_current_tally

heating_tally = openmc.Tally()
heating_tally.filters = [cell_filter]
heating_tally.scores = ["heating"]
model.tallies.append(heating_tally)

# -----------------------------------------------------------------------------
# Define Structural_materials nuclides fractions and totals
# -----------------------------------------------------------------------------
_ATOMS_PER_BARNCM_TO_ATOMS = 1.0e24  # (barn*cm)^-1 * cm^3 * 1e24 = atoms

def build_structural_nuclides(structural_materials: List[openmc.Material]) -> List[str]:
    """
    Sorted list of nuclides present in the provided *structural* materials.

    This is a whitelist used to:
      - restrict per-nuclide tallies (damage-energy, H/He production) to only
        nuclides that exist in structural materials set.
    """
    s: Set[str] = set()
    for m in structural_materials:
        # Use the same API as the later accumulation (atom densities).
        for nuc in m.get_nuclide_atom_densities().keys():
            s.add(str(nuc))
    return sorted(s)

def build_structural_maps_vo(
    model: openmc.Model,
    *,
    cell_ids: List[int],
    mix_recipes_obj: Dict[str, Dict[openmc.Material, float]],
    structural_materials: List[openmc.Material],
    ) -> Tuple[
    List[str],                    # structural_nuclides (whitelist)
    Dict[int, Dict[str, float]],  # cell_nuclide_atoms (total atoms per nuclide)
    Dict[int, Dict[str, float]],  # cell_struct_nuclide_atoms (struct atoms per nuclide)
    Dict[int, float],             # cell_total_atoms_struct (sum struct atoms)
    Dict[int, Dict[str, float]],  # cell_struct_origin_frac (struct/total per nuclide)
    ]:
    """
    For each cell:
        total_by_nuc(n)  = sum_i [ dens_i(n) * vol * vf_i * 1e24 ] over ALL components i
        struct_by_nuc(n) = sum_{i in structural materials} [ dens_i(n) * vol * vf_i * 1e24 ]
        f_struct_origin(n) = struct_by_nuc(n) / total_by_nuc(n)

    Important:
      - We only track nuclides in `structural_nuclides` (whitelist), defined as
        nuclides present in the *structural materials list*.
      - Structural-origin is determined by which *component material* contributed,
        via membership in `structural_materials_set`.
    """
    structural_materials_set: Set[openmc.Material] = set(structural_materials)

    structural_nuclides: List[str] = build_structural_nuclides(structural_materials)
    structural_nuclide_set: Set[str] = set(structural_nuclides)

    all_cells = model.geometry.get_all_cells()

    cell_nuclide_atoms: Dict[int, Dict[str, float]] = {}
    cell_struct_nuclide_atoms: Dict[int, Dict[str, float]] = {}
    cell_total_atoms_struct: Dict[int, float] = {}
    cell_struct_origin_frac: Dict[int, Dict[str, float]] = {}

    for cid in map(int, cell_ids):
        cell = all_cells[cid]
        vol = float(cell.volume or 0.0)

        cell_nuclide_atoms[cid] = {}
        cell_struct_nuclide_atoms[cid] = {}
        cell_total_atoms_struct[cid] = 0.0
        cell_struct_origin_frac[cid] = {}

        if vol <= 0.0 or not isinstance(cell.fill, openmc.Material):
            continue

        fill_mat: openmc.Material = cell.fill
        fill_name = str(fill_mat.name or "")

        recipe = mix_recipes_obj.get(fill_name)
        if recipe is None:
            raise KeyError(
                f"[mix_recipes_obj] Missing recipe for fill material name '{fill_name}' "
                f"(cell {cid})."
            )

        mats = list(recipe.keys())
        vfs = [float(recipe[m]) for m in mats]

        total_by_nuc: Dict[str, float] = {}
        struct_by_nuc: Dict[str, float] = {}

        for m_i, vf_i in zip(mats, vfs):
            if vf_i <= 0.0:
                continue

            nd_i = m_i.get_nuclide_atom_densities()
            scale = vol * vf_i * _ATOMS_PER_BARNCM_TO_ATOMS

            # Structural-origin classification is by 'component material'
            is_struct = (m_i in structural_materials_set)

            for nuc, dens in nd_i.items():
                nuc = str(nuc)

                # only nuclides present in structural materials list.
                if nuc not in structural_nuclide_set:
                    continue

                a = float(dens) * scale
                total_by_nuc[nuc] = total_by_nuc.get(nuc, 0.0) + a
                if is_struct:
                    struct_by_nuc[nuc] = struct_by_nuc.get(nuc, 0.0) + a

        # Structural-origin fraction per nuclide
        frac_by_nuc: Dict[str, float] = {}
        for nuc, n_tot in total_by_nuc.items():
            n_str = float(struct_by_nuc.get(nuc, 0.0))
            frac_by_nuc[nuc] = min((n_str / n_tot), 1.0) if n_tot > 0.0 else 0.0

        # store
        cell_nuclide_atoms[cid] = total_by_nuc
        cell_struct_nuclide_atoms[cid] = struct_by_nuc
        cell_total_atoms_struct[cid] = float(sum(struct_by_nuc.values()))
        cell_struct_origin_frac[cid] = frac_by_nuc

    return (
        structural_nuclides,
        cell_nuclide_atoms,
        cell_struct_nuclide_atoms,
        cell_total_atoms_struct,
        cell_struct_origin_frac,
    )


structural_nuclides, cell_nuclide_atoms, cell_struct_nuclide_atoms, cell_total_atoms_struct, cell_struct_origin_frac = (
    build_structural_maps_vo(
        model,
        cell_ids=cell_ids,
        mix_recipes_obj=MIX_RECIPES_OBJ_NORM,
        structural_materials=structure_material_list,
    )
)

# -----------------------------------------------------------------------------
# DPA + gas production tallies per cell
# -----------------------------------------------------------------------------
cell_filters_by_cid = {int(cid): openmc.CellFilter([int(cid)]) for cid in cell_ids}
structural_nuclide_set = set(structural_nuclides)

dpa_gas_tallies: Dict[int, openmc.Tally] = {}

for cid in cell_ids:
    cid = int(cid)

    # NOTE: cell_nuclide_atoms was built using the whitelist already, so
    # this membership filter is redundant but harmless. Keep it for clarity.
    nuclides_in_cell = sorted(
        nuc for nuc in cell_nuclide_atoms.get(cid, {}).keys()
        if nuc in structural_nuclide_set
    )
    if not nuclides_in_cell:
        continue

    tg = openmc.Tally()
    tg.filters = [cell_filters_by_cid[cid], n_particle_filter]
    tg.scores = [
        "damage-energy",
        "H1-production", "H2-production", "H3-production",
        "He3-production", "He4-production",
    ]
    tg.nuclides = nuclides_in_cell

    model.tallies.append(tg)
    dpa_gas_tallies[cid] = tg

# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------
model.export_to_model_xml(path=SLAB_MODEL_XML)

# TODO: is this necessary? How would these come to exist? Suggest to remove if not needed
# remove redundant defaults
# (REPLY): I think the model.init_lib(output=False) is creating these outputs (not the tallies)
redundant_files = ["geometry.xml", "materials.xml", "settings.xml", "tallies.xml"]
for f in redundant_files:
    if os.path.exists(f):
        os.remove(f)
