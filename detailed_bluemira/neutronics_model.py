#!/usr/bin/env python3

"""
Generates a wedge of the EU DEMO tokamak, based on a CAD file generated from the
BlueMira/Process codes.
"""

from __future__ import annotations

# built-in modules
from enum import Enum
from typing import Dict, Tuple, List, Iterable, Set, Optional
import math
import re
import json
from pathlib import Path
import os
import sys

import openmc
import numpy as np
from openmc_plasma_source import tokamak_source
import pydagmc

# -----------------------------------------------------------------------------
# Inputs / paths
# -----------------------------------------------------------------------------
BASE_DIR = Path.cwd()

_DAGMC_MODEL_FILE = (BASE_DIR / "eudemo_f_1_27a.h5m")
INPUT_JSON = (BASE_DIR / "Tokamak_inputs.json")

RUN_DIR = (BASE_DIR / "neutronics_run")
RUN_DIR.mkdir(parents=True, exist_ok=True)

NEUTRONICS_MODEL_XML = RUN_DIR / "model.xml"
SURFACE_SOURCE_FILE = RUN_DIR / "surface_source.h5"

model = openmc.Model()

# Load materials compositions
module_path = (BASE_DIR.parent / "materials")
sys.path.append(str(module_path))
import materials

# -----------------------------------------------------------------------------
# GEOMETRY
# -----------------------------------------------------------------------------
dagmc_universe = openmc.DAGMCUniverse(filename=str(_DAGMC_MODEL_FILE))
pydagmc_model = pydagmc.Model(str(_DAGMC_MODEL_FILE))

# reserve IDs
openmc.reserve_ids([v.id for v in pydagmc_model.volumes], cls=openmc.Cell)
openmc.reserve_ids([s.id for s in pydagmc_model.surfaces], cls=openmc.Surface)

def azimuthal_plane(theta_deg, boundary_type=None, name=None, surface_id=None):
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

number_sectors = 16
theta0 = 0.0
theta1 = theta0 + 360 / number_sectors

# generate the CSG wedge to place the DAGMC model within
cut_lo = azimuthal_plane(theta0, boundary_type="reflective", surface_id=100_001)
cut_hi = azimuthal_plane(theta1, boundary_type="reflective", surface_id=100_002)

z_min = openmc.ZPlane(z0=-3000, boundary_type="vacuum", surface_id=100_003)
z_max = openmc.ZPlane(z0=+2000, boundary_type="vacuum", surface_id=100_004)
r_out = openmc.ZCylinder(r=3000, boundary_type="vacuum", surface_id=100_005)

sector_region = (+cut_lo & -cut_hi) & (+z_min & -z_max) & (-r_out)
sector_cell = openmc.Cell(
    cell_id=10_000,
    region=sector_region,
    fill=dagmc_universe,
    name="sector_container",
)
model.geometry = openmc.Geometry(root=[sector_cell])

# -----------------------------------------------------------------------------
# MATERIALS
# -----------------------------------------------------------------------------
# TODO: add citations for where these densities come from
# TODO: need to review materials.py for correctness for all materials
ss316   = materials.ss316Ln_ig(7.93)
ccz     = materials.CuCrZr(8.9)
h       = materials.Helium(0.0001785)
nb3sn   = materials.Nb3Sn(5.7)
epoxy   = materials.Epoxy(1.207)
bronze  = materials.Bronze(8.8775)
nbti    = materials.NbTi(6.538)
c       = materials.Cu(8.96)
ss304_b4 = materials.ss304_b4(7.8)

# tungsten, density based on PNNL material compendium value
w       = materials.W(19.3)

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
coolant_material = materials.Water(0.866)

# Breeder material (dependent on the breeder type) (PbLi for WCLL)
breeder_material = materials.PbLi(0.90, 9.8)

# -----------------------------------------------------------------------------
# Structural / breeder / coolant definitions
# -----------------------------------------------------------------------------
# TODO: very unclear what these lists are. Why are the other solid materials not included in structure_material_list?
structure_material_list: list[openmc.Material] = [structural_material, armor_material, ss316]
breeder_material_list: list[openmc.Material] = [breeder_material]
coolant_material_list: list[openmc.Material] = [coolant_material]

# ARMOR DEFINITION (not a mixture)
#armor = armor_material
#armor.name = "Armor"
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
    "First_Wall":  {armor_material: 0.0027, coolant_material: 0.14268, structural_material: 0.85462},
    # IB
    "ib_layer_1":  {breeder_material: 0.833, coolant_material: 0.025, structural_material: 0.139},
    "ib_layer_2":  {breeder_material: 0.858, coolant_material: 0.018, structural_material: 0.124},
    "ib_layer_3":  {breeder_material: 0.8132, coolant_material: 0.0158, structural_material: 0.171},
    "ib_layer_4":  {breeder_material: 0.8132, coolant_material: 0.0158, structural_material: 0.171},
    "ib_layer_5":  {breeder_material: 0.8132, coolant_material: 0.0158, structural_material: 0.171},
    "ib_layer_6":  {breeder_material: 0.8132, coolant_material: 0.0158, structural_material: 0.171},
    "ib_layer_7":  {breeder_material: 0.427, coolant_material: 0.016, structural_material: 0.558},
    "ib_layer_8":  {coolant_material: 0.486, structural_material: 0.514},
    # OB
    "ob_layer_1":  {breeder_material: 0.833, coolant_material: 0.025, structural_material: 0.139},
    "ob_layer_2":  {breeder_material: 0.858, coolant_material: 0.018, structural_material: 0.124},
    "ob_layer_3":  {breeder_material: 0.8132, coolant_material: 0.0158, structural_material: 0.171},
    "ob_layer_4":  {breeder_material: 0.8132, coolant_material: 0.0158, structural_material: 0.171},
    "ob_layer_5":  {breeder_material: 0.8132, coolant_material: 0.0158, structural_material: 0.171},
    "ob_layer_6":  {breeder_material: 0.8132, coolant_material: 0.0158, structural_material: 0.171},
    "ob_layer_7":  {breeder_material: 0.427, coolant_material: 0.016, structural_material: 0.558},
    "ob_layer_8":  {coolant_material: 0.486, structural_material: 0.514},
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
# SOURCE (openmc-plasma-source)
# -----------------------------------------------------------------------------
my_source = tokamak_source(
    angles=(0.0, math.pi/8),
    elongation=1.739,
    ion_density_centre=6.3e19, #6.8e19,
    ion_density_pedestal=5.355e19, #5.78e19,
    ion_density_peaking_factor=1.0,
    ion_density_separatrix=0.945e19, #1.02e19,
    ion_temperature_centre=23.7e3,
    ion_temperature_pedestal=5.5e3,
    ion_temperature_separatrix=0.1e3,
    ion_temperature_peaking_factor=1.45,
    ion_temperature_beta=2.0,
    major_radius=840.67,
    minor_radius=300.2,
    pedestal_radius=0.94 * 300.2,
    mode="H",
    shafranov_factor=0.44789,
    triangularity=0.333,
    fuel={"D": 0.5, "T": 0.5},
)

# -----------------------------------------------------------------------------
# SETTINGS
# -----------------------------------------------------------------------------

model.settings = openmc.Settings()
model.settings.dagmc = True
model.settings.photon_transport = True
model.settings.run_mode = "fixed source"
model.settings.source = my_source

# Set TALLY_CONVERGENCE_THRESHOLD to 0.01 (1%) or 0.001 (0.1%)
TALLY_CONVERGENCE_THRESHOLD = 0.1

# Choose initial batches * particles per batch > 15-20 * max_particles
model.settings.batches = 10           
model.settings.trigger_active = True
model.settings.trigger_batch_interval = 5   # check triggers every N batches
model.settings.particles = 100_000
model.settings.trigger_max_batches = 1000     # hard ceiling

# output particle track, selected at random
import random
model.settings.track = [(1, 1, random.randint(1, model.settings.particles))]

# TODO: change 245 and 56 to not be hard-coded
_WRITE_SOURCE = False
if _WRITE_SOURCE:
    model.settings.surf_source_write = {
        "surface_ids": [245],
        "max_particles": 1_000_000,
        "cellto": 56,
    }
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
    default_chunk_key: str = "OB_1_b6",
    gap_cm: float = 2.0,
    start_cm: float = 0.0,
) -> dict:
    with INPUT_JSON.open("r", encoding="utf-8") as f:
        data = json.load(f)

    inv = data["inventory"]
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

    if default_chunk_key in ob_by_key:
        selected_chunk_cell_ids = ob_by_key[default_chunk_key]
    elif default_chunk_key in ib_by_key:
        selected_chunk_cell_ids = ib_by_key[default_chunk_key]
    else:
        raise KeyError(
            f"default_chunk_key={default_chunk_key!r} not found. "
            f"Valid keys include: {ALL_KEYS[:10]}{' ...' if len(ALL_KEYS) > 10 else ''}"
        )

    # -------------------------
    # Armor cell ids (ALL)
    # -------------------------
    ob_armor_idx = ob_n_layers
    ib_armor_idx = ib_n_layers

    armor_cell_ids = []

    for chunk in ob_by_key.values():
        if len(chunk) > ob_armor_idx:
            armor_cell_ids.append(int(chunk[ob_armor_idx]))

    for chunk in ib_by_key.values():
        if len(chunk) > ib_armor_idx:
            armor_cell_ids.append(int(chunk[ib_armor_idx]))

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
        "armor_cell_ids": armor_cell_ids,
        "ALL_KEYS": ALL_KEYS,
        "OB_CHUNK_SIZE": OB_CHUNK_SIZE,
        "IB_CHUNK_SIZE": IB_CHUNK_SIZE,
        "selected_chunk_key": default_chunk_key,
        "selected_chunk_cell_ids": selected_chunk_cell_ids,
        "cell_ids_for_key": cell_ids_for_key,
        "radial_bins_for_key": radial_bins_for_key,
    }

chunk_cells = build_breeder_chunks(INPUT_JSON, default_chunk_key="OB_1_b6")

cell_ids = chunk_cells["cell_ids_all"]
cell_ids_selected_chunk = chunk_cells["selected_chunk_cell_ids"]
armor_cell_ids = chunk_cells["armor_cell_ids"]

info, all_surface_ids, external_surface_ids, internal_surface_ids = dagmc_volume_surface_info(
    pydagmc_model, cell_ids_selected_chunk
)

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

# Filters
cell_filter = openmc.CellFilter(cell_ids)
chunk_cell_filter = openmc.CellFilter(cell_ids_selected_chunk) # If only using OB_1_b6 to check trigger
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
flux_tally_total.filters = [chunk_cell_filter, n_particle_filter]
flux_tally_total.scores = ["flux"]
flux_tally_total.triggers = [
    openmc.Trigger(trigger_type="rel_err", threshold=TALLY_CONVERGENCE_THRESHOLD)
]
model.tallies.append(flux_tally_total)

# TODO: why is this only looking at the neutrons? I guess we are only computing the albedos for the neutrons?
# (REPLY) I have been checked neutrons only. But I agree this should have been more in depth explored with photons.
# I will introduce Photons analysis after 2/26/2026.
t_current_tally = openmc.Tally()
t_current_tally.filters = [t_surf_filter, particle_filter]
t_current_tally.scores = ["current"]
model.tallies.append(t_current_tally)

p_current_tallies: dict[int, openmc.Tally] = {}
for cid in cell_ids_selected_chunk:
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
model.export_to_model_xml(path=NEUTRONICS_MODEL_XML)

# TODO: is this necessary? How would these come to exist? Suggest to remove if not needed
# remove redundant defaults
# (REPLY): I think the model.init_lib(output=False) is creating these outputs (not the tallies)
redundant_files = ["geometry.xml", "materials.xml", "settings.xml", "tallies.xml"]
for f in redundant_files:
    if os.path.exists(f):
        os.remove(f)

# Check if any of the source sites overlap with the material regions; this can be commented
# out to make the model run faster but is helpful to make sure the plasma source is
# behaving as we expect
_CHECK_SOURCE = False
if _CHECK_SOURCE:
    openmc.lib.init(output=False, args=[str(NEUTRONICS_MODEL_XML)])
    n_samples = 100000
    particles = openmc.lib.sample_external_source(n_samples=n_samples)

    in_cells = {}
    for p in particles:
        c = openmc.lib.find_cell([p.r[0], p.r[1], p.r[2]])
        i = c[0].id
        if (i not in in_cells):
            in_cells[i] = 1
        else:
            in_cells[i] += 1

    print('\nPercent of source sites in each cell: ')
    for k, v in in_cells.items():
        print("Cell : ", k, " % Sites: ", v/n_samples * 100)
    openmc.lib.finalize()

    # end check on source site overlaps
