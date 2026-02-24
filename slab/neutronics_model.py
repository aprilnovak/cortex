#!/usr/bin/env python3
"""
plot_dagmc_sector_and_surface_source.py

DAGMC sector (with 4 reflective cut planes) + overlay surface_source.h5 points
WITHOUT shifting the point cloud.

Outputs:
  plot_xy.png
  plot_xz.png
  plot_yz.png
"""

from __future__ import annotations

import os
import sys
from enum import Enum
from pathlib import Path
from typing import Dict, Tuple, List, Iterable, Set, Optional

import numpy as np
import h5py
import openmc
from pymoab import core, types
import pydagmc

module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "materials"))
if module_path not in sys.path:
    sys.path.append(module_path)
import materials  

# =============================================================================
# USER INPUTS
# =============================================================================
_DAGMC_MODEL_FILE = "eudemo_f_1_27a.h5m"
# Fixed-source run configuration
SURF_SOURCE_FILE = "surface_source.h5" 

TARGET_VOL_ID = 66

# Radial vacuum cylinders (cm)
R_IN_CM = 1140.0
R_OUT_CM = 1375.0

# source plotting controls
N_SOURCE_PLOT = 20_000
DPI = 300

# SlicePlot controls 
PLOT_ORIGIN = (1250.0, 5.0, 0.0)
PLOT_WIDTH = (1000.0, 1000.0)     # (u_width, v_width) for each basis
PLOT_PIXELS = (3000, 3000)

# =============================================================================
# Read surface-source positions (plotting)
# =============================================================================
def load_positions_cm(path: Path) -> np.ndarray:
    with h5py.File(path, "r") as f:
        if "source_bank" not in f:
            raise KeyError(f"{path} missing 'source_bank'")
        sb = f["source_bank"][()]
        if "r" not in sb.dtype.fields:
            raise KeyError(f"'source_bank' fields={list(sb.dtype.fields)} missing 'r'")
        r = sb["r"]

        if r.dtype.kind in ("f", "d"):
            arr = np.asarray(r, dtype=float)
            if arr.ndim != 2 or arr.shape[1] != 3:
                raise ValueError(f"Expected (N,3), got {arr.shape}")
            return arr

        if r.dtype.fields and {"x", "y", "z"}.issubset(r.dtype.fields):
            return np.column_stack((r["x"], r["y"], r["z"])).astype(float)

        raise TypeError(f"Unrecognized dtype for r: {r.dtype}")

def summarize_positions(r_xyz: np.ndarray) -> dict:
    mins = r_xyz.min(axis=0)
    maxs = r_xyz.max(axis=0)
    means = r_xyz.mean(axis=0)
    return {
        "N": int(r_xyz.shape[0]),
        "xmin": float(mins[0]), "xmax": float(maxs[0]), "xmean": float(means[0]),
        "ymin": float(mins[1]), "ymax": float(maxs[1]), "ymean": float(means[1]),
        "zmin": float(mins[2]), "zmax": float(maxs[2]), "zmean": float(means[2]),
    }

def overlay_points_on_png(png_path: Path, basis: str, origin, width, points_xyz: np.ndarray, out_path: Path):
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg

    img = mpimg.imread(str(png_path))
    ny, nx = img.shape[0], img.shape[1]

    pts = np.asarray(points_xyz, dtype=float)

    if basis == "xy":
        u, v = pts[:, 0], pts[:, 1]
        u0, v0 = origin[0], origin[1]
        wu, wv = width[0], width[1]
    elif basis == "xz":
        u, v = pts[:, 0], pts[:, 2]
        u0, v0 = origin[0], origin[2]
        wu, wv = width[0], width[1]
    elif basis == "yz":
        u, v = pts[:, 1], pts[:, 2]
        u0, v0 = origin[1], origin[2]
        wu, wv = width[0], width[1]
    else:
        raise ValueError(f"Unknown basis: {basis}")

    umin, umax = u0 - wu / 2.0, u0 + wu / 2.0
    vmin, vmax = v0 - wv / 2.0, v0 + wv / 2.0

    xpix = (u - umin) / (umax - umin) * (nx - 1)
    ypix = (1.0 - (v - vmin) / (vmax - vmin)) * (ny - 1)

    fig, ax = plt.subplots(figsize=(nx / DPI, ny / DPI), dpi=DPI)
    ax.imshow(img)
    ax.scatter(xpix, ypix, s=1, marker=".", alpha=0.7)
    ax.set_axis_off()
    fig.tight_layout(pad=0)
    fig.savefig(str(out_path), dpi=DPI, bbox_inches="tight", pad_inches=0)
    plt.close(fig)

def make_slice_plot(plot_id: int, basis: str, origin, width, pixels) -> openmc.SlicePlot:
    p = openmc.SlicePlot(plot_id=plot_id)
    p.basis = basis
    p.origin = origin
    p.width = width
    p.pixels = pixels
    p.color_by = "cell"
    return p

# =============================================================================
# DAGMC helpers
# =============================================================================
class Orientation(Enum):
    FORWARD = 1
    REVERSE = -1

def dagmc_volume_surface_info(pydagmc_model: pydagmc.Model, volume_ids: List[int]):
    volume_ids_set = set(int(v) for v in volume_ids)
    result: Dict[int, dict] = {}
    surface_ids_all = set()
    surface_ids_external = set()
    surface_ids_internal = set()

    for vol_id in volume_ids:
        vol_id = int(vol_id)
        volume = pydagmc_model.volumes_by_id[vol_id]
        all_list, external_list, internal_list = [], [], []

        for surface in volume.surfaces:
            parent_volumes = surface.senses
            orientation = Orientation.FORWARD if (parent_volumes and parent_volumes[0].id == vol_id) else Orientation.REVERSE

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
            if (not adj_all) or adj_ext:
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

    return result, sorted(surface_ids_all), sorted(surface_ids_external), sorted(surface_ids_internal)

def dagmc_bounding_box(pydagmc_model, volume_id):
    volume = pydagmc_model.volumes_by_id[volume_id]
    triangle_coords = volume.triangle_coords
    min_coords = np.min(triangle_coords, axis=0)
    max_coords = np.max(triangle_coords, axis=0)
    return openmc.BoundingBox(min_coords, max_coords)

def get_breeder_reflective_cuts(mb: core.Core, pydagmc_model: pydagmc.Model, volume_id: int):
    """
    Finds 2 Z-dominant and 2 Y-dominant reflective planes for a given volume.
    Returns list of dicts: {id, area, region, normal, d}.
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

        p = openmc.Plane(a=float(avg_n[0]), b=float(avg_n[1]), c=float(avg_n[2]),
                         d=d_openmc, boundary_type="reflective")

        region = -p if orientation == Orientation.FORWARD else +p

        area = float(0.5 * np.sum(norm_vals))
        entry = {"id": sid, "area": area, "region": region, "normal": avg_n, "d": d_openmc}

        if abs(avg_n[2]) > 0.7:
            z_candidates.append(entry)
        elif abs(avg_n[1]) > 0.7:
            y_candidates.append(entry)

    poloidal_cuts = sorted(z_candidates, key=lambda x: x["area"], reverse=True)[:2]
    toroidal_cuts = sorted(y_candidates, key=lambda x: x["area"], reverse=True)[:2]
    return poloidal_cuts + toroidal_cuts

# -----------------------------------------------------------------------------
# Normalize and create openmc.Material.Mixture
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

# -----------------------------------------------------------------------------
# Define Structural_materials nuclides fractions and totals
# -----------------------------------------------------------------------------
_ATOMS_PER_BARNCM_TO_ATOMS = 1.0e24  # (barn*cm)^-1 * cm^3 * 1e24 = atoms

def build_structural_nuclides(structural_materials: List[openmc.Material]) -> List[str]:
    """
    Sorted list of nuclides present in the provided *structural* materials.

    This is a whitelist used to:
      - restrict per-nuclide tallies (damage-energy, H/He production) to only
        nuclides that exist in your structural materials set.
      - restrict book-keeping for f_struct_origin(n) to the same whitelist.
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

# =============================================================================
# Build OpenMC model: DAGMC in a cut sector cell
# =============================================================================

model = openmc.Model()

mb = core.Core()
mb.load_file(_DAGMC_MODEL_FILE)

dagmc_universe = openmc.DAGMCUniverse(filename=_DAGMC_MODEL_FILE)
pydagmc_model = pydagmc.Model(str(dagmc_universe.filename))

openmc.reserve_ids([v.id for v in pydagmc_model.volumes], cls=openmc.Cell)
openmc.reserve_ids([s.id for s in pydagmc_model.surfaces], cls=openmc.Surface)

r_out = openmc.ZCylinder(r=float(R_OUT_CM), boundary_type="vacuum")
r_in = openmc.ZCylinder(r=float(R_IN_CM), boundary_type="vacuum")

final_cuts = get_breeder_reflective_cuts(mb, pydagmc_model, TARGET_VOL_ID)
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

# base materials
# TODO: add citations for where these densities come from
# TODO: need to review materials.py for correctness for all materials
ss316   = materials.ss316Ln_ig(7.93)
ccz     = materials.CuCrZr(8.9)

# tungsten, density based on PNNL material compendium value
w       = materials.W(19.3)

h       = materials.Helium(0.0001785)
nb3sn   = materials.Nb3Sn(5.7)
epoxy   = materials.Epoxy(1.207)
bronze  = materials.Bronze(8.8775)
nbti    = materials.NbTi(6.538)
c       = materials.Cu(8.96)
ss304_b4 = materials.ss304_b4(7.8)

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

# --------------------------------------------
# Structural / breeder / coolant definitions 
# --------------------------------------------
structure_material_list: list[openmc.Material] = [structural_material, armor_material, ss316]
breeder_material_list: list[openmc.Material] = [breeder_material]
coolant_material_list: list[openmc.Material] = [coolant_material]

# --------------------------------------------
# MIX RECIPES
# --------------------------------------------
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

MIX_RECIPES_OBJ_NORM = normalize_mix_recipes_obj(MIX_RECIPES_OBJ)

mixed = build_and_set_model_materials_from_obj_recipes_vo(
    model,
    recipes_obj=MIX_RECIPES_OBJ_NORM,
)

# -----------------------------------------------------------------------------
# SETTINGS
# -----------------------------------------------------------------------------
model.settings = openmc.Settings()
model.settings.dagmc = True
model.settings.photon_transport = True
model.settings.batches = 10
model.settings.particles = 2_000_000
model.settings.run_mode = "fixed source"
model.settings.surf_source_read = {'path': SURF_SOURCE_FILE}

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

#info, all_surface_ids, external_surface_ids, internal_surface_ids = dagmc_volume_surface_info(
#pydagmc_model, cell_ids_equatorial_ob)

# =============================================================================
# TALLIES
# =============================================================================
# Source intensity -> scaling (keep for later postprocessing)
total_power = 2e9
number_sectors = 16
section_power = total_power / number_sectors
ev_to_joule = 1.60218e-19
ev_fusion = 17.6e6
convert_e = ev_to_joule * ev_fusion
neutron_source_rate = section_power / convert_e
s_in_y = (365 * 24 * 60 * 60)

# Find correct dagmc cell_ids
cell_ids = [56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66]
xcentroids_ob = (0.1, 1.1, 6.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 84.4, 156.0)
widths = (0.2, 1.8, 8.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 28.0, 110)

# Filters
cell_filter = openmc.CellFilter([int(c) for c in cell_ids])
particle_filter = openmc.ParticleFilter(bins=["neutron", "photon"])
n_particle_filter = openmc.ParticleFilter(bins=["neutron"])

energies = openmc.mgxs.GROUP_STRUCTURES["CCFE-709"]
energy_filter = openmc.EnergyFilter(energies)
unit_lethargy = [float(np.log(energies[i + 1] / energies[i])) for i in range(len(energies) - 1)]

model.tallies = openmc.Tallies()

# Flux (MG)
flux_tally = openmc.Tally()
flux_tally.filters = [cell_filter, particle_filter, energy_filter]
flux_tally.scores = ["flux"]
model.tallies.append(flux_tally)

# Flux (total)
flux_tally_total = openmc.Tally()
flux_tally_total.filters = [cell_filter, particle_filter]
flux_tally_total.scores = ["flux"]
model.tallies.append(flux_tally_total)

# Heating
heating_tally = openmc.Tally()
heating_tally.filters = [cell_filter]
heating_tally.scores = ["heating"]
model.tallies.append(heating_tally)

# -----------------------------------------------------------------------------
# Define Structural_materials nuclides fractions and totals
# -----------------------------------------------------------------------------
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
model.export_to_model_xml(path="neutronics_model.xml")

# remove redundant defaults
redundant_files = ["geometry.xml", "materials.xml", "settings.xml", "tallies.xml"]
for f in redundant_files:
    if os.path.exists(f):
        os.remove(f)


# PLOT
# load points 
r = load_positions_cm(SURF_SOURCE_FILE)
s = summarize_positions(r)

print(f"[surf_source] N={s['N']}")
print(f"[surf_source] x min/max/mean = {s['xmin']:.3f}, {s['xmax']:.3f}, {s['xmean']:.3f} cm")
print(f"[surf_source] y min/max/mean = {s['ymin']:.3f}, {s['ymax']:.3f}, {s['ymean']:.3f} cm")
print(f"[surf_source] z min/max/mean = {s['zmin']:.3f}, {s['zmax']:.3f}, {s['zmean']:.3f} cm")

r_plot = r  

if N_SOURCE_PLOT is not None and r_plot.shape[0] > N_SOURCE_PLOT:
    idx = np.random.default_rng(123).choice(r_plot.shape[0], size=N_SOURCE_PLOT, replace=False)
    r_plot = r_plot[idx]
print("[plot] using points:", r_plot.shape)


# --- plots ---
plots = openmc.Plots()

p_xy = make_slice_plot(
    plot_id=1, basis="xy",
    origin=PLOT_ORIGIN,
    width=PLOT_WIDTH,
    pixels=PLOT_PIXELS,
)
plots.append(p_xy)

p_xz = make_slice_plot(
    plot_id=2, basis="xz",
    origin=PLOT_ORIGIN,
    width=PLOT_WIDTH,
    pixels=PLOT_PIXELS,
)
plots.append(p_xz)

model.plots = plots

model.export_to_xml()
openmc.plot_geometry(output=False)

in_xy = Path("plot_1.png")
in_xz = Path("plot_2.png")

overlay_points_on_png(in_xy, "xy", p_xy.origin, p_xy.width, r_plot, Path("plot_xy.png"))
overlay_points_on_png(in_xz, "xz", p_xz.origin, p_xz.width, r_plot, Path("plot_xz.png"))

print("Wrote: plot_xy.png, plot_xz.png")

