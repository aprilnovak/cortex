#!/usr/bin/env python3
"""
neutronics_model.py
===================
Assembles the OpenMC Model (settings + tallies) and writes model.xml.

Imports from:
  inputs.py   : user configuration flags
  geometry.py : model, dagmc_universe, pydagmc_model, chunk helpers

Exports (for post-processing scripts)
--------------------------------------
model, all_cells
neutron_source_rate, s_in_y, ev_to_joule
flux_tally, heating_tally, t_current_tally, p_current_tallies, dpa_gas_tallies
structural_nuclides, cell_nuclide_atoms, cell_struct_nuclide_atoms
cell_total_atoms_struct, cell_struct_origin_frac, cell_volumes
"""

from __future__ import annotations

import json
import random
from pathlib import Path
from typing import Dict, List, Set, Tuple
import numpy as np

import openmc

import inputs as cfg
import geometry as geo

# ──────────────────────────────────────────────────────────────────────────────
# Aliases from geometry
# ──────────────────────────────────────────────────────────────────────────────
model         = geo.model
ob_by_key     = geo.ob_by_key
ib_by_key     = geo.ib_by_key
cell_ids_all  = geo.cell_ids_all
info          = geo.info
external_surface_ids = geo.external_surface_ids

# ──────────────────────────────────────────────────────────────────────────────
# SOURCE
# ──────────────────────────────────────────────────────────────────────────────
model.settings = openmc.Settings()
#model.settings.dagmc            = True # Not necessary anymore?
model.settings.photon_transport = cfg.DO_PHOTON_TRANSPORT
model.settings.run_mode         = "fixed source"

if cfg.SIM_TYPE == "tokamak":
    # Plasma source (tokamak_neutron_source)
    from tokamak_neutron_source import (
        FluxMap, FractionalFuelComposition,
        TokamakNeutronSource, TransportInformation,
    )
    from tokamak_neutron_source.profile import ParabolicPedestalProfile
    from tokamak_neutron_source.reactions import Reactions

    temperature_profile = ParabolicPedestalProfile(
        2.37754249767383570e+01, 5.5, 0.1, 2.0, 1.45, 0.94
    )
    density_profile = ParabolicPedestalProfile(
        9.83393828196113777e+19, 5.86476334188244419e+19,
        3.44986078934261350e+19, 1.0, 2.0, 0.94,
    )
    density_profile.set_scale(
        6.30197378312059699e+19 / 7.47634856031986811e+19
    )

    my_source = TokamakNeutronSource(
        transport=TransportInformation.from_parameterisations(
            ion_temperature_profile=temperature_profile,
            fuel_density_profile=density_profile,
            rho_profile=__import__("numpy").linspace(0, 1, 30),
            fuel_composition=FractionalFuelComposition(D=0.5, T=0.5),
        ),
        source_type=[Reactions.D_T, Reactions.D_D],
        flux_map=FluxMap.from_eqdsk(str(cfg.EQUILIBRIUM_FILE)),
        cell_side_length=0.05,
    )
    #dt_fusion_power = 1.91480972334663875e+09
    #my_source.normalise_fusion_power(dt_fusion_power)
    model.settings.source = my_source.to_openmc_source(
            start_angle=0.0,
            end_angle=np.radians(22.5),
            )

    # Surface source write (optional) 
    if cfg.DO_WRITE_SURFACE_SOURCE:
        model.settings.surf_source_write = {
            "surface_ids":  [cfg.SURFACE_SOURCE_SURFACE_ID],
            "max_particles": cfg.SURFACE_SOURCE_MAX_PARTICLES,
            "cellto":        cfg.SURFACE_SOURCE_CELL_TO,
        }

    neutron_ratio_source: float = 1.0

else:
    # Surface source read (slab / fast_slab)
    import h5py
    import numpy as np

    # load surface_source.h5 file
    if not cfg.SURFACE_SOURCE_FILE.is_file():
        raise FileNotFoundError(
            f"Surface source file not found: {cfg.SURFACE_SOURCE_FILE}"
        )
    model.settings.surf_source_read = {"path": str(cfg.SURFACE_SOURCE_FILE)}

    # Check ratio of neutrons on source file
    with h5py.File(str(cfg.SURFACE_SOURCE_FILE), "r") as f:
        particles = f["source_bank"]["particle"][:]
    total      = len(particles)
    n_neutrons = int((particles == 0).sum())
    neutron_ratio_source: float = total / n_neutrons if n_neutrons > 0 else 1.0
    print(
        f"Surface source: {total} particles, {n_neutrons} neutrons "
        f"→ neutron_ratio_source = {neutron_ratio_source:.4f}"
    )

# ──────────────────────────────────────────────────────────────────────────────
# Run settings
# ──────────────────────────────────────────────────────────────────────────────
model.settings.particles = cfg.PARTICLES_PER_BATCH

if cfg.USE_TRIGGER:
    model.settings.batches                = cfg.BATCHES
    model.settings.trigger_active         = True
    model.settings.trigger_batch_interval = cfg.TRIGGER_BATCH_INTERVAL
    model.settings.trigger_max_batches    = cfg.TRIGGER_MAX_BATCHES
else:
    model.settings.batches        = cfg.FIXED_BATCHES
    model.settings.trigger_active = False

#if cfg.TRACKS:
#    if cfg.SIM_TYPE == "tokamak":
#        model.settings.track = [
#            (1, 1, random.randint(1, cfg.PARTICLES_PER_BATCH))
#        ]

# ──────────────────────────────────────────────────────────────────────────────
# DAGMC volume sync  (requires settings to be defined first)
# ──────────────────────────────────────────────────────────────────────────────
openmc.Cell.reset_ids()
openmc.Surface.reset_ids()

model.init_lib(output=False)
model.sync_dagmc_universes()
model.finalize_lib()

openmc.reserve_ids(list(model.geometry.get_all_cells()),    cls=openmc.Cell)
openmc.reserve_ids(list(model.geometry.get_all_surfaces()), cls=openmc.Surface)

dagmc_universe_cells: Dict[int, openmc.Cell] = geo.dagmc_universe.get_all_cells()
for volume in geo.pydagmc_model.volumes:
    dagmc_universe_cells[volume.id].volume       = volume.volume
    dagmc_universe_cells[volume.id].bounding_box = geo.dagmc_bounding_box(geo.pydagmc_model, volume.id)


all_cells: Dict[int, openmc.Cell] = model.geometry.get_all_cells()

# ──────────────────────────────────────────────────────────────────────────────
# ============================ DEBUGGING =======================================
# Export cell/material/volume mapping for depletion (may not be necessary anymore)
# ──────────────────────────────────────────────────────────────────────────────
_cell_map = []
for cid, cell in all_cells.items():
    mat = cell.fill
    if not isinstance(mat, openmc.Material):
        continue
    _cell_map.append({
        "cell_id":       int(cid),
        "mat_id":        int(mat.id),
        "mat_name":      str(mat.name or ""),
        "volume_cm3":    float(cell.volume or 0.0),
    })

_cell_map_path = cfg.NEUTRONICS_RESULTS_DIR / "cell_material_volume_map.json"
_cell_map_path.parent.mkdir(parents=True, exist_ok=True)
with open(_cell_map_path, "w", encoding="utf-8") as _f:
    json.dump(_cell_map, _f, indent=2)
print(f"[neutronics_model] Cell/material/volume map written to: {_cell_map_path}")
# ============================ DEBUGGING =======================================

# ──────────────────────────────────────────────────────────────────────────────
# STRUCTURAL MAPS  →  exported to JSON for neutronics_post.py
# ──────────────────────────────────────────────────────────────────────────────
_ATOMS_PER_BARNCM = 1.0e24

def _structural_nuclides(structural_materials: List[openmc.Material]) -> List[str]:
    """
    Sorted list of nuclides present in the provided *structural* materials.

    This is a whitelist used to:
      - restrict per-nuclide tallies (damage-energy, H/He production) to only
        nuclides that exist in structural materials set.
    """
    s: Set[str] = set()
    for m in structural_materials:
        for nuc in m.get_nuclide_atom_densities().keys():
            s.add(str(nuc))
    return sorted(s)

def build_structural_maps_vo(
    model: openmc.Model,
    *,
    cell_ids:             List[int],
    mix_recipes_norm:     Dict[str, Dict[openmc.Material, float]],
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

    structural_nuclides: List[str] = _structural_nuclides(structural_materials)
    structural_nuclide_set: Set[str] = set(structural_nuclides)

    all_cells = model.geometry.get_all_cells()

    cell_nuclide_atoms:        Dict[int, Dict[str, float]] = {}
    cell_struct_nuclide_atoms: Dict[int, Dict[str, float]] = {}
    cell_total_atoms_struct:   Dict[int, float]            = {}
    cell_struct_origin_frac:   Dict[int, Dict[str, float]] = {}

    for cid in map(int, cell_ids):
        cell = all_cells[cid]
        vol  = float(cell.volume or 0.0)

        cell_nuclide_atoms[cid]        = {}
        cell_struct_nuclide_atoms[cid] = {}
        cell_total_atoms_struct[cid]   = 0.0
        cell_struct_origin_frac[cid]   = {}

        if vol <= 0.0 or not isinstance(cell.fill, openmc.Material):
            continue

        fill_mat:  openmc.Material = cell.fill
        fill_name: str             = str(fill_mat.name or "")

        recipe = mix_recipes_norm.get(fill_name)
        if recipe is None:
            raise KeyError(
                f"No recipe for fill material name '{fill_name}' (cell {cid}). "
                f"Ensure the name matches a key in {cfg.BREEDER_TYPE} mix_recipes_norm."
            )

        total_by_nuc:  Dict[str, float] = {}
        struct_by_nuc: Dict[str, float] = {}

        for m_i, vf_i in recipe.items():
            if vf_i <= 0.0:
                continue

            scale     = vol * vf_i * _ATOMS_PER_BARNCM
            is_struct = (m_i in structural_materials_set)

            for nuc, dens in m_i.get_nuclide_atom_densities().items():
                nuc = str(nuc)
                if nuc not in structural_nuclide_set:
                    continue

                a = float(dens) * scale
                total_by_nuc[nuc] = total_by_nuc.get(nuc, 0.0) + a
                if is_struct:
                    struct_by_nuc[nuc] = struct_by_nuc.get(nuc, 0.0) + a

        # Structural-origin fraction per nuclide
        frac_by_nuc: Dict[str, float] = {
            nuc: min(float(struct_by_nuc.get(nuc, 0.0)) / n_tot, 1.0) if n_tot > 0.0 else 0.0
            for nuc, n_tot in total_by_nuc.items()
        }

        cell_nuclide_atoms[cid]        = total_by_nuc
        cell_struct_nuclide_atoms[cid] = struct_by_nuc
        cell_total_atoms_struct[cid]   = float(sum(struct_by_nuc.values()))
        cell_struct_origin_frac[cid]   = frac_by_nuc

    return (
        structural_nuclides,
        cell_nuclide_atoms,
        cell_struct_nuclide_atoms,
        cell_total_atoms_struct,
        cell_struct_origin_frac,
    )

(
    structural_nuclides,
    cell_nuclide_atoms,
    cell_struct_nuclide_atoms,
    cell_total_atoms_struct,
    cell_struct_origin_frac,
) = build_structural_maps_vo(
    model,
    cell_ids             = geo.cell_ids_all,
    mix_recipes_norm     = geo.MIX_RECIPES_NORM,
    structural_materials = geo.structure_material_list,
)

# Save information for neutronics_post.py
def export_structural_maps_json(output_path: Path) -> None:
    """
    export/save structural maps + cell volumes to JSON so that
    neutronics_post.py can load them without importing neutronics_model.py.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cell_volumes = {
        str(cid): float(cell.volume or 0.0)
        for cid, cell in all_cells.items()
    }

    payload = {
        "breeder_type":              cfg.BREEDER_TYPE,
        "structural_nuclides":       structural_nuclides,
        "cell_volumes":              cell_volumes,
        "cell_nuclide_atoms":        {str(k): v for k, v in cell_nuclide_atoms.items()},
        "cell_struct_nuclide_atoms": {str(k): v for k, v in cell_struct_nuclide_atoms.items()},
        "cell_total_atoms_struct":   {str(k): v for k, v in cell_total_atoms_struct.items()},
        "cell_struct_origin_frac":   {str(k): v for k, v in cell_struct_origin_frac.items()},
    }

    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    print(f"[neutronics_model] Structural maps written to {output_path}")


# Write structural maps JSON 
STRUCTURAL_MAPS_JSON = cfg.NEUTRONICS_RESULTS_DIR / "structural_maps.json"
export_structural_maps_json(STRUCTURAL_MAPS_JSON)

# ──────────────────────────────────────────────────────────────────────────────
# TALLIES
# ──────────────────────────────────────────────────────────────────────────────
model.tallies = openmc.Tallies()

# Albedo chunk cells 
_albedo_cells: List[int] = (
    ob_by_key[cfg.ALBEDO_CHUNK_KEY]
    if cfg.ALBEDO_CHUNK_KEY in ob_by_key
    else ib_by_key[cfg.ALBEDO_CHUNK_KEY]
)

# Energy structure
energies = openmc.mgxs.GROUP_STRUCTURES["CCFE-709"]
e_filter = openmc.EnergyFilter(energies)

# Particle filters
particle_filter = openmc.ParticleFilter(["neutron", "photon"])
n_filter        = openmc.ParticleFilter(["neutron"])

# Cell filters 
_slab_cells = (
    ob_by_key[cfg.SLAB_CHUNK_KEY]
    if cfg.SLAB_CHUNK_KEY in ob_by_key
    else ib_by_key[cfg.SLAB_CHUNK_KEY]
)

_tally_cells = _slab_cells if cfg.SIM_TYPE in ("slab" ,"fast_slab") else cell_ids_all
cell_filter_all = openmc.CellFilter(_tally_cells)

cell_filter_trigger = openmc.CellFilter(
    ob_by_key.get(cfg.TRIGGER_CHUNK_KEY, []) or
    ib_by_key.get(cfg.TRIGGER_CHUNK_KEY, [])
)

# 1. Flux spectrum tally (neutron and photon) 
flux_tally = openmc.Tally()
flux_tally.filters = [cell_filter_all, particle_filter, e_filter]
flux_tally.scores  = ["flux"]
model.tallies.append(flux_tally)

# 2. Trigger tally (convergence monitoring)
if cfg.USE_TRIGGER: 
    flux_tally_trigger = openmc.Tally()
    flux_tally_trigger.filters = [cell_filter_trigger, n_filter]
    flux_tally_trigger.scores  = ["flux"]
    flux_tally_trigger.triggers = [
        openmc.Trigger(
            trigger_type="rel_err",
            threshold=cfg.TALLY_CONVERGENCE_THRESHOLD,
        )
    ]
    model.tallies.append(flux_tally_trigger)

# 3. Heating tally (all cells) 
heating_tally = openmc.Tally()
heating_tally.filters = [cell_filter_all]
heating_tally.scores  = ["heating"]
model.tallies.append(heating_tally)

# 4. Total current tally 
# (tokamak only, for albedo and calculating
#  _J_in in neutron_source_rate for slab models) 
t_current_tally:   openmc.Tally | None        = None
p_current_tallies: Dict[int, openmc.Tally]    = {}

if cfg.SIM_TYPE == "tokamak":
    t_surf_filter   = openmc.SurfaceFilter(external_surface_ids)
    t_current_tally = openmc.Tally()
    t_current_tally.filters = [t_surf_filter, particle_filter]
    t_current_tally.scores  = ["current"]
    model.tallies.append(t_current_tally)

    # 5. Per-cell partial current tallies (albedo chunk)
    for cid in _albedo_cells:
        ocell          = dagmc_universe_cells[cid]
        surf_ids_local = [
            int(s["surface_id"])
            for s in info.get(cid, {}).get("all_surfaces", [])
        ]
        if not surf_ids_local:
            continue
        cell_from_f = openmc.CellFromFilter([ocell])
        surf_f      = openmc.SurfaceFilter(surf_ids_local)
        pt = openmc.Tally()
        pt.filters = [cell_from_f, surf_f, particle_filter]
        pt.scores  = ["current"]
        model.tallies.append(pt)
        p_current_tallies[cid] = pt


# 6. DPA + gas production tallies (per cell, structural nuclides only)
structural_nuclide_set = set(structural_nuclides)
cell_filters_by_cid = {
    int(cid): openmc.CellFilter([int(cid)]) for cid in _tally_cells
}

dpa_gas_tallies: Dict[int, openmc.Tally] = {}

for cid in _tally_cells:
    cid = int(cid)
    nuclides_in_cell = sorted(
        nuc for nuc in cell_nuclide_atoms.get(cid, {})
        if nuc in structural_nuclide_set
    )
    if not nuclides_in_cell:
        continue

    tg = openmc.Tally()
    tg.filters  = [cell_filters_by_cid[cid], n_filter]
    tg.scores   = [
        "damage-energy",
        "H1-production", "H2-production", "H3-production",
        "He3-production", "He4-production",
    ]
    tg.nuclides = nuclides_in_cell
    model.tallies.append(tg)
    dpa_gas_tallies[cid] = tg

# 7 TBD (n,gamma)
# 8 flux in the VV_port_filling? (tokamak only)

# ──────────────────────────────────────────────────────────────────────────────
# EXPORT model.xml
# ──────────────────────────────────────────────────────────────────────────────
NEUTRONICS_MODEL_XML = cfg.NEUTRONICS_RUN_DIR / "model.xml"
model.export_to_model_xml(path=NEUTRONICS_MODEL_XML)

# Write neutron_ratio_source for neutronics_post.py
_run_meta = {
    "neutron_ratio_source":  neutron_ratio_source,
    "sim_type":              cfg.SIM_TYPE,
    "breeder_type":          cfg.BREEDER_TYPE,
    "tally_ids": {
        "flux_spectrum":   flux_tally.id,
        "heating":         heating_tally.id,
        "total_current":   t_current_tally.id if t_current_tally else None,
        "partial_current": {str(cid): t.id for cid, t in p_current_tallies.items()},
        "dpa_gas":         {str(cid): t.id for cid, t in dpa_gas_tallies.items()},
    },
}
with open(cfg.NEUTRONICS_RESULTS_DIR / "run_meta.json", "w", encoding="utf-8") as _f:
    json.dump(_run_meta, _f, indent=2)
print(f"[neutronics_model] run_meta.json written to: {cfg.NEUTRONICS_RESULTS_DIR}")

# Remove stray XML files that model.init_lib() created
for redundant_files in ["geometry.xml", "materials.xml", "settings.xml", "tallies.xml"]:
    p = Path(".") / redundant_files
    if p.exists():
        p.unlink()

# ──────────────────────────────────────────────────────────────────────────────
# OPTIONAL: source site cell-distribution check
# ──────────────────────────────────────────────────────────────────────────────
# Check if any of the source sites overlap with the material regions; this can be commented
# out to make the model run faster but is helpful to make sure the plasma source is
# behaving as we expect

if cfg.DO_CHECK_SOURCE:
    openmc.lib.init(output=False, args=[str(NEUTRONICS_MODEL_XML)])
    n_samples  = 100_000
    particles  = openmc.lib.sample_external_source(n_samples=n_samples)
    in_cells:  Dict[int, int] = {}
    for p in particles:
        c = openmc.lib.find_cell([p.r[0], p.r[1], p.r[2]])
        i = c[0].id
        in_cells[i] = in_cells.get(i, 0) + 1
    print("\nPercent of source sites in each cell:")
    for k, v in sorted(in_cells.items()):
        print(f"  Cell {k:6d}: {v / n_samples * 100:.2f} %")
    openmc.lib.finalize()

print(f"\n[neutronics_model] model.xml written to: {NEUTRONICS_MODEL_XML}")