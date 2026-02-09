#!/usr/bin/env python3
"""
neutronics_model.py 
"""

from __future__ import annotations

# built-in modules
from enum import Enum
from typing import Dict, Tuple, List, Iterable
import math
import re
import json
from pathlib import Path

import openmc
import numpy as np
from openmc_plasma_source import tokamak_source
import pydagmc

import os
import sys
module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "materials"))
sys.path.append(module_path)
import materials  

# -----------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------
_DAGMC_MODEL_FILE = "eudemo_f_1_27a.h5m"
INPUT_JSON = Path("Tokamak_inputs.json")

model = openmc.Model()

# -----------------------------------------------------------------------------
# GEOMETRY 
# -----------------------------------------------------------------------------
dagmc_universe = openmc.DAGMCUniverse(filename=_DAGMC_MODEL_FILE)
pydagmc_model = pydagmc.Model(str(dagmc_universe.filename))

# reserve IDs
openmc.reserve_ids([v.id for v in pydagmc_model.volumes], cls=openmc.Cell)
openmc.reserve_ids([s.id for s in pydagmc_model.surfaces], cls=openmc.Surface)

def azimuthal_plane(theta_deg, boundary_type=None, name=None, surface_id=None):
    theta = math.radians(theta_deg)
    if abs(theta_deg - 90.0) < 1e-6:
        return openmc.Plane(a=1.0, b=0.0, c=0.0, d=0.0,
                            boundary_type=boundary_type, name=name, surface_id=surface_id)
    a = -math.tan(theta)
    b = 1.0
    c = 0.0
    d = 0.0
    return openmc.Plane(a=a, b=b, c=c, d=d,
                        boundary_type=boundary_type, name=name, surface_id=surface_id)

number_sectors = 16
theta0 = 0.0
theta1 = theta0 + 360 / number_sectors

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
ss316   = materials.ss316Ln_ig(7.93)
ccz     = materials.CuCrZr(8.9)
eurofer = materials.eurofer97(7.87)
w       = materials.W(19.3)
water   = materials.Water(0.866)
h       = materials.Helium(0.0001785)
nb3sn   = materials.Nb3Sn(5.7)
epoxy   = materials.Epoxy(1.207)
bronze  = materials.Bronze(8.8775)
nbti    = materials.NbTi(6.538)
c       = materials.Cu(8.96)
PbLi    = materials.PbLi(0.90, 9.8)

# Plasma Region
plasma = openmc.Material(name="Plasma_Region")
plasma.set_density("g/cm3", 1e-6)
plasma.add_element("He", 1.0)
# Armor
armor = w
armor.name = "Armor" 
# Cryostat
CS = materials.ss316Ln_ig(7.93)
CS.name = "Cryostat"
# Radiation Shielding
RS = materials.ss304_b4(7.8)
RS.name = "RadiationShield_all"

# -----------------------------------------------------------------------------
# MIX RECIPES
# (easier to find nuclides ratios with this)
# -----------------------------------------------------------------------------
MIX_RECIPES_OBJ: dict[str, dict[openmc.Material, float]] = {
    # FW
    "First_Wall":  {w: 0.0027, water: 0.14268, eurofer: 0.85462},
    # IB
    "ib_layer_1":  {PbLi: 0.85, water: 0.05, eurofer: 0.10},
    "ib_layer_2":  {PbLi: 0.85, water: 0.05, eurofer: 0.10},
    "ib_layer_3":  {PbLi: 0.85, water: 0.05, eurofer: 0.10},
    "ib_layer_4":  {PbLi: 0.85, water: 0.05, eurofer: 0.10},
    "ib_layer_5":  {PbLi: 0.85, water: 0.05, eurofer: 0.10},
    "ib_layer_6":  {PbLi: 0.85, water: 0.05, eurofer: 0.10},
    "ib_layer_7":  {PbLi: 0.85, water: 0.05, eurofer: 0.10},
    "ib_layer_8":  {water: 0.486, eurofer: 0.514},
    # OB
    "ob_layer_1":  {PbLi: 0.85, water: 0.05, eurofer: 0.10},
    "ob_layer_2":  {PbLi: 0.85, water: 0.05, eurofer: 0.10},
    "ob_layer_3":  {PbLi: 0.85, water: 0.05, eurofer: 0.10},
    "ob_layer_4":  {PbLi: 0.85, water: 0.05, eurofer: 0.10},
    "ob_layer_5":  {PbLi: 0.85, water: 0.05, eurofer: 0.10},
    "ob_layer_6":  {PbLi: 0.85, water: 0.05, eurofer: 0.10},
    "ob_layer_7":  {PbLi: 0.85, water: 0.05, eurofer: 0.10},
    "ob_layer_8":  {water: 0.486, eurofer: 0.514},
    # Divertor
    "Divertor":    {ccz: 0.00552, c: 0.00438, eurofer: 0.5238, w: 0.01026, water: 0.45604},
    # VV
    "VV_IB":       {ss316: 0.6, water: 0.4},
    "VV_OB":       {ss316: 0.6, water: 0.4},
    "VV_ports_all":{ss316: 0.6, water: 0.4},
    "VV_port_Fill":{ss316: 0.6, water: 0.4},
    # Coils
    "PC_PFC":      {nbti: 0.02895, c: 0.1169, epoxy: 0.18, bronze: 0.0735, h: 0.1682, ss316: 0.43245},
    "TFcoil":      {nb3sn: 0.02895, c: 0.1169, epoxy: 0.18, bronze: 0.0735, h: 0.1682, ss316: 0.43245},
}

# -----------------------------------------------------------------------------
# Build mixed materials and assign model.materials
# -----------------------------------------------------------------------------
def build_and_set_model_materials_from_obj_recipes(
    model: openmc.Model,
    *,
    recipes_obj: dict[str, dict[openmc.Material, float]],
    static_materials: list[openmc.Material] | None = None,
    basis: str = "vo",          # "vo" or "wo"
    normalize: bool = True,
) -> dict[str, openmc.Material]:
    """
    Build mixed materials from object-keyed recipes and assign them to model.materials.

    """
    import inspect

    def _normalize(vals: list[float]) -> list[float]:
        s = float(sum(vals))
        if s <= 0.0:
            raise ValueError("Fractions sum to <= 0")
        return [float(v) / s for v in vals]

    def _weight_to_volume_fracs(mats: list[openmc.Material], w: list[float]) -> list[float]:
        # v_i ∝ w_i / rho_i
        v_unnorm: list[float] = []
        for m, wi in zip(mats, w):
            rho = getattr(m, "density", None)
            if rho is None:
                raise ValueError(
                    f"[MIX_RECIPES] basis='wo' requested but material '{m.name}' has no density set."
                )
            rho = float(rho)
            if rho <= 0.0:
                raise ValueError(f"[MIX_RECIPES] material '{m.name}' has non-positive density: {rho}")
            v_unnorm.append(float(wi) / rho)
        return _normalize(v_unnorm)

    static_materials = list(static_materials or [])
    mixed: dict[str, openmc.Material] = {}

    # Detect whether this OpenMC supports 'basis' kwarg on mix_materials
    mix_sig = inspect.signature(openmc.Material.mix_materials)
    supports_basis = ("basis" in mix_sig.parameters)

    basis = str(basis).strip().lower()
    if basis not in ("vo", "wo"):
        raise ValueError(f"basis must be 'vo' or 'wo' (got {basis!r})")

    for name, comp in recipes_obj.items():
        mats = list(comp.keys())
        vals = [float(v) for v in comp.values()]

        s = float(sum(vals))
        if s <= 0.0:
            raise ValueError(f"[MIX_RECIPES] {name}: fractions sum to {s}")

        if normalize and abs(s - 1.0) > 1e-12:
            vals = [v / s for v in vals]

        # If OpenMC supports basis, let OpenMC do it.
        # Otherwise, force everything to volume fractions before mixing.
        if supports_basis:
            mixed[name] = openmc.Material.mix_materials(
                mats, vals, basis=basis, name=name
            )
        else:
            if basis == "wo":
                vals_vo = _weight_to_volume_fracs(mats, vals)
            else:
                vals_vo = _normalize(vals) if normalize else vals  # keep as provided

            mixed[name] = openmc.Material.mix_materials(
                mats, vals_vo, name=name
            )

        # Ensure the name is present (necessary for DAGMC matching)
        mixed[name].name = name

    model.materials = openmc.Materials(static_materials + list(mixed.values()))
    return mixed

# no mixture materials
NO_MIXTURE = [armor, CS, RS, plasma]
# Build all materials and assign to model.materials
mixed = build_and_set_model_materials_from_obj_recipes(
    model,
    recipes_obj=MIX_RECIPES_OBJ,
    static_materials=NO_MIXTURE,
    basis="vo",
)

# -----------------------------------------------------------------------------
# Structural / breeder / coolant definitions 
# -----------------------------------------------------------------------------
structure_material_list: list[openmc.Material] = [eurofer, w, ss316]
breeder_material_list: list[openmc.Material] = [PbLi]
coolant_material_list: list[openmc.Material] = [water]

# -----------------------------------------------------------------------------
# SOURCE (openmc-plasma-source)
# -----------------------------------------------------------------------------
my_source = tokamak_source(
    angles=(0.0, math.pi/8),
    elongation=1.739,
    ion_density_centre=6.8e19,
    ion_density_pedestal=5.78e19,
    ion_density_peaking_factor=1,
    ion_density_separatrix=1.02e19,
    ion_temperature_centre=23.7e3,
    ion_temperature_pedestal=5.5e3,
    ion_temperature_separatrix=0.1e3,
    ion_temperature_peaking_factor=8.06,
    ion_temperature_beta=6,
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
model.settings.batches = 10
model.settings.particles = 100_000
model.settings.run_mode = "fixed source"
model.settings.source = my_source

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
    armor_cm = float(geom["armor_cm"])
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

    pre_breeder_shift = float(geom["armor_cm"] + geom["fw_cm"])
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

chunk_cells = build_breeder_chunks(INPUT_JSON, default_equatorial_ob_key="OB_1_b6")

cell_ids = chunk_cells["cell_ids_all"]
cell_ids_equatorial_ob = chunk_cells["equatorial_ob_cell_ids"]

info, all_surface_ids, external_surface_ids, internal_surface_ids = dagmc_volume_surface_info(
    pydagmc_model, cell_ids_equatorial_ob
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
particle_filter = openmc.ParticleFilter(bins=["neutron", "photon"])
t_surf_filter = openmc.SurfaceFilter(external_surface_ids)
n_particle_filter = openmc.ParticleFilter(bins=["neutron"])
energies = openmc.mgxs.GROUP_STRUCTURES["CCFE-709"]
energy_filter = openmc.EnergyFilter(energies)
unit_lethargy = [np.log(energies[i + 1] / energies[i]) for i in range(len(energies) - 1)]

model.tallies = openmc.Tallies()

flux_tally = openmc.Tally()
flux_tally.filters = [cell_filter, particle_filter, energy_filter]
flux_tally.scores = ["flux"]
model.tallies.append(flux_tally)

flux_tally_total = openmc.Tally()
flux_tally_total.filters = [cell_filter, particle_filter]
flux_tally_total.scores = ["flux"]
model.tallies.append(flux_tally_total)

t_current_tally = openmc.Tally()
t_current_tally.filters = [t_surf_filter, n_particle_filter]
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
    p_current_tally.filters = [cell_from_filter, surf_filter, n_particle_filter]
    p_current_tally.scores = ["current"]

    model.tallies.append(p_current_tally)
    p_current_tallies[cid] = p_current_tally

heating_tally = openmc.Tally()
heating_tally.filters = [cell_filter]
heating_tally.scores = ["heating"]
model.tallies.append(heating_tally)

# -----------------------------------------------------------------------------
# STRUCTURAL NUCLIDE FRACTIONS 
# -----------------------------------------------------------------------------
def compute_structural_nuclide_fractions(
    model: openmc.Model,
    *,
    cell_ids: list[int],
    mix_recipes: dict[str, dict[openmc.Material, float]],
    structure_materials: set[openmc.Material],
) -> Tuple[Dict[int, Dict[str, float]], Dict[int, float], Dict[int, Dict[str, float]]]:
    """
    Build per-cell structural-nuclide atoms / totals / fractions.

    - Uses cell.fill.name ONLY to select a mix recipe (because fill is the MIXED material).
    - Recipe components are openmc.Material objects; we DO NOT rely on component names.
    - structure_materials is a set of openmc.Material objects (identity match).
    """
    all_cells = model.geometry.get_all_cells()

    cell_struct_nuclide_atoms: Dict[int, Dict[str, float]] = {}
    cell_total_atoms_struct: Dict[int, float] = {}
    cell_struct_nuclide_frac: Dict[int, Dict[str, float]] = {}

    for cid in cell_ids:
        cid = int(cid)
        cell = all_cells[cid]
        vol = float(cell.volume or 0.0)

        struct_atoms: Dict[str, float] = {}

        if vol <= 0.0 or not isinstance(cell.fill, openmc.Material):
            cell_struct_nuclide_atoms[cid] = {}
            cell_total_atoms_struct[cid] = 0.0
            cell_struct_nuclide_frac[cid] = {}
            continue

        fill_mat: openmc.Material = cell.fill
        fill_name = str(fill_mat.name or "")

        # mixture case
        if fill_name in mix_recipes:
            for comp_mat, vf in mix_recipes[fill_name].items():
                if comp_mat not in structure_materials:
                    continue
                nd = comp_mat.get_nuclide_atom_densities()
                for nuc, dens in nd.items():
                    struct_atoms[str(nuc)] = struct_atoms.get(str(nuc), 0.0) + float(dens) * vol * float(vf) * 1e24

        # single-material case
        else:
            if fill_mat in structure_materials:
                nd = fill_mat.get_nuclide_atom_densities()
                for nuc, dens in nd.items():
                    struct_atoms[str(nuc)] = float(dens) * vol * 1e24

        total_atoms = float(sum(struct_atoms.values()))
        frac = {n: a / total_atoms for n, a in struct_atoms.items()} if total_atoms > 0.0 else {}

        cell_struct_nuclide_atoms[cid] = struct_atoms
        cell_total_atoms_struct[cid] = total_atoms
        cell_struct_nuclide_frac[cid] = frac

    return cell_struct_nuclide_atoms, cell_total_atoms_struct, cell_struct_nuclide_frac

# Build globals used by post-processing
cell_struct_nuclide_atoms, cell_total_atoms_struct, cell_struct_nuclide_frac = (
    compute_structural_nuclide_fractions(
        model,
        cell_ids=cell_ids,
        mix_recipes=MIX_RECIPES_OBJ,
        structure_materials=set(structure_material_list),
    )
)

# -----------------------------------------------------------------------------
# Structural nuclides list 
# -----------------------------------------------------------------------------
solid_nuclides: set[str] = set()
for mat in structure_material_list:
    # get_nuclide_densities returns dict nuclide -> density
    for name in mat.get_nuclide_densities().keys():
        solid_nuclides.add(str(name))
solid_nuclides = sorted(solid_nuclides)

# -----------------------------------------------------------------------------
# Atom totals for whole cell + solid-only
# -----------------------------------------------------------------------------
all_cells = model.geometry.get_all_cells()
all_nuclides: set[str] = set()

cell_nuclide_atoms: Dict[int, Dict[str, float]] = {}
cell_total_atoms: Dict[int, float] = {}
cell_total_atoms_solid: Dict[int, float] = {}
cell_atomic_ratio: Dict[int, Dict[str, float]] = {}

for cid in cell_ids:
    cell = all_cells[int(cid)]
    if not isinstance(cell.fill, openmc.Material):
        cell_nuclide_atoms[int(cid)] = {}
        cell_total_atoms[int(cid)] = 0.0
        cell_total_atoms_solid[int(cid)] = 0.0
        cell_atomic_ratio[int(cid)] = {}
        continue

    all_nuclides.update(cell.fill.get_nuclides())

    # OpenMC helper gives atoms directly if volume is set
    nuc_atoms = cell.fill.get_nuclide_atoms(volume=cell.volume)
    tot_atoms = float(sum(nuc_atoms.values()))

    cell_nuclide_atoms[int(cid)] = {str(k): float(v) for k, v in nuc_atoms.items()}
    cell_total_atoms[int(cid)] = tot_atoms

    cell_atomic_ratio[int(cid)] = (
        {nuc: atoms / tot_atoms for nuc, atoms in cell_nuclide_atoms[int(cid)].items()}
        if tot_atoms > 0.0 else {}
    )

    solid_total_atoms = sum(
        atoms for nuc, atoms in cell_nuclide_atoms[int(cid)].items()
        if nuc in solid_nuclides
    )
    cell_total_atoms_solid[int(cid)] = float(solid_total_atoms)

# -----------------------------------------------------------------------------
# STRUCTURAL ORIGIN FRACTION PER NUCLIDE
# -----------------------------------------------------------------------------
# For each nuclide in solid_nuclides:
#   f_struct_of_nuclide = (atoms of that nuclide coming from structural materials)
#                         / (total atoms of that nuclide in the FULL mixture)
# Also store atom percentage of each nuclide in the FULL mixture.
# -----------------------------------------------------------------------------

cell_struct_origin_frac: Dict[int, Dict[str, float]] = {}     # {cid: {nuc: N_struct(n)/N_total(n)}}
cell_nuclide_atom_percent: Dict[int, Dict[str, float]] = {}  # {cid: {nuc: 100*N_total(n)/N_total_all}}

for cid in cell_ids:
    cid = int(cid)

    total_by_nuc = cell_nuclide_atoms.get(cid, {})         # full mixture atoms per nuclide
    struct_by_nuc = cell_struct_nuclide_atoms.get(cid, {}) # struct-only atoms per nuclide

    # total atoms of *all* nuclides in the mixture
    total_atoms_all = float(cell_total_atoms.get(cid, 0.0))

    # fraction of each "structural nuclide" that comes from structure
    frac_map: Dict[str, float] = {}
    atompct_map: Dict[str, float] = {}

    for nuc in solid_nuclides:
        n_tot = float(total_by_nuc.get(nuc, 0.0))
        n_str = float(struct_by_nuc.get(nuc, 0.0))

        # solid_nuclide_fraction from solid_materials
        frac_map[nuc] = (n_str / n_tot) if n_tot > 0.0 else 0.0

        # atom% of that nuclide in the full mixture
        atompct_map[nuc] = (100.0 * n_tot / total_atoms_all) if total_atoms_all > 0.0 else 0.0

    cell_struct_origin_frac[cid] = frac_map
    cell_nuclide_atom_percent[cid] = atompct_map

# -----------------------------------------------------------------------------
# DPA + gas production tallies per cell
# -----------------------------------------------------------------------------
cell_filters_by_cid = {int(cid): openmc.CellFilter([int(cid)]) for cid in cell_ids}
solid_set = set(solid_nuclides)

dpa_gas_tallies: Dict[int, openmc.Tally] = {}

for cid in cell_ids:
    cid = int(cid)

    nuclides_in_cell = sorted(
        nuc for nuc in cell_nuclide_atoms[cid].keys()
        if nuc in solid_set
    )
    if not nuclides_in_cell:
        continue

    tg = openmc.Tally()
    tg.filters = [cell_filters_by_cid[cid], n_particle_filter]
    tg.scores = [
        "damage-energy",
        "H1-production", "He3-production", "He4-production",
    ]
    tg.nuclides = nuclides_in_cell

    model.tallies.append(tg)
    dpa_gas_tallies[cid] = tg

# -----------------------------------------------------------------------------
# Export
# -----------------------------------------------------------------------------
model.export_to_xml()