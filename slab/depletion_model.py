#!/usr/bin/env python3
"""
depletion_model.py

Runs:
  (1) D1S shutdown dose tally for a single breeder chunk (OB_KEY)
  (2) Depletion microXS + flux (fixed-source, KEEPING surface_source)
      using Strategy B: set VOID material in every non-target DAGMC cell

Key idea:
  DAGMCUniverse still contains all DAGMC volume-cells. Even if your wrapper
  "sector_cell" clips the geometry, OpenMC can still encounter/initialize
  materials for many DAGMC cells. If any DAGMC cell references a material ID
  that is not present, you'll get:
     "Material with ID 'X' not found for DAGMC cell Y"
  or crashes.

So:
  - Keep surf_source_read (you said you don't want to remove it)
  - Force all non-target material-filled DAGMC cells to a known VOID material
  - Only keep real materials in the 11 OB_1_b6 cells
"""

from __future__ import annotations

from pathlib import Path
import os
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import openmc
import openmc.data
import openmc.deplete
from openmc.deplete import d1s

# -----------------------------------------------------------------------------
# materials module path (your setup)
# -----------------------------------------------------------------------------
module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "materials"))
if module_path not in sys.path:
    sys.path.append(module_path)
import materials  

# -----------------------------------------------------------------------------
# Import from neutronics_model.py
# -----------------------------------------------------------------------------
from neutronics_model import (
    neutron_source_rate,
    sector_cell,
    model,
    pydagmc_model, 
    all_cells,
    build_breeder_chunks,
)

# -----------------------------------------------------------------------------
# USER INPUTS
# -----------------------------------------------------------------------------
path_file = Path(__file__).resolve().parent
INPUT_JSON = (path_file / "Tokamak_inputs.json").resolve()
OB_KEY = "OB_1_b6"

# Chain file (base)
chain_path = (Path(__file__).resolve().parent.parent / "depletion_chain" / "chain_endfb80_sfr.xml").resolve()
if not chain_path.exists():
    raise FileNotFoundError(f"Chain file not found: {chain_path}")
chain = openmc.deplete.Chain.from_xml(str(chain_path))

# Register for both Python-side and OpenMC executable
openmc.config["chain_file"] = str(chain_path)
model.settings.depletion = {"chain_file": str(chain_path)}

# -----------------------------------------------------------------------------
# Build chunk mapping (from JSON), then restrict to cells PRESENT in model
# -----------------------------------------------------------------------------
_cells_chunk = build_breeder_chunks(INPUT_JSON, default_equatorial_ob_key=OB_KEY)
ob_by_key = _cells_chunk["ob_by_key"]
OB_CHUNK_SIZE = _cells_chunk["OB_CHUNK_SIZE"]

if OB_KEY not in ob_by_key:
    raise KeyError(f"{OB_KEY} not found in ob_by_key. Available (sample): {list(ob_by_key)[:10]}")

# Cells from JSON mapping
ob_cells_json = [int(x) for x in ob_by_key[OB_KEY]]

# What cells actually exist in this *reduced* model geometry
present_cell_ids = set(int(c.id) for c in model.geometry.get_all_cells().values())
# Restrict to those present
ob_cells = [cid for cid in ob_cells_json if cid in present_cell_ids]

print(f"[chunk] {OB_KEY}: {len(ob_cells_json)} from JSON, {len(ob_cells)} present in wrapper")
print(f"[chunk] present ids: {ob_cells}")
if len(ob_cells) == 0:
    raise RuntimeError(f"{OB_KEY}: none of the JSON chunk cells are present in model.geometry")

# Radial binning (centroids/widths/edges) for OB_KEY
centroids_cm, widths_cm, edges_cm = _cells_chunk["radial_bins_for_key"](OB_KEY)

# Sanity checks: binning returns per-layer for one OB chunk
if len(ob_cells_json) != OB_CHUNK_SIZE:
    raise ValueError(f"{OB_KEY}: JSON chunk size {len(ob_cells_json)} != OB_CHUNK_SIZE {OB_CHUNK_SIZE}")
if len(centroids_cm) != OB_CHUNK_SIZE:
    raise ValueError(f"{OB_KEY}: len(centroids_cm)={len(centroids_cm)} != OB_CHUNK_SIZE={OB_CHUNK_SIZE}")

# Map cell_id -> radial center 
center_by_cell_full = dict(zip(ob_cells_json, map(float, centroids_cm)))

# -----------------------------------------------------------------------------
# Time grids / source rates
# -----------------------------------------------------------------------------
s_to_h = 3600.0
y_to_s = 24.0 * 365.0 * s_to_h
to_μSv = 1e-6
to_mSv = 1e-9

# irradiation (5 years)
irradiation_time_y = np.array([5.0])
irradiation_time = (irradiation_time_y * y_to_s).tolist()

# cooling (1e-8 y to 1000 y)
timesteps_years = np.concatenate([
    np.logspace(-8, -4, 6),
    np.logspace(-4,  0, 8)[1:],   # drop 1e-4 duplicate
    np.logspace( 0,  3, 14)[1:],  # drop 1e0 duplicate
])
cooling_times = (timesteps_years * y_to_s).tolist()

timesteps = irradiation_time + cooling_times

SURFACE_SOURCE_POWER_RATIO = 7.171062e-02
constant_power_ratio = 0.3 * SURFACE_SOURCE_POWER_RATIO
source_rates = [constant_power_ratio * neutron_source_rate] * len(irradiation_time) + [0.0] * len(cooling_times)

# -----------------------------------------------------------------------------
# Volumes from DAGMC cells (already synced in neutronics_model.py)
# -----------------------------------------------------------------------------
vol_by_cell = {int(cid): float(cell.volume) for cid, cell in all_cells.items() if getattr(cell, "volume", None) is not None}

# -----------------------------------------------------------------------------
# D1S: Build dose tally for ONLY the present OB cells
# -----------------------------------------------------------------------------
dose_cells_all = list(ob_cells)
cell_filter2 = openmc.CellFilter(dose_cells_all)

# Per-cell volumes for normalization
vol_by_cell_all = {}
for cid in dose_cells_all:
    v = vol_by_cell.get(int(cid), None)
    if v is None:
        raise KeyError(f"Missing volume for cell {cid} in vol_by_cell (all_cells).")
    vol_by_cell_all[int(cid)] = float(v)

# Save the current tallies to restore later
orig_tallies = [t for t in model.tallies]

# Dose coefficients + tally
energies, pSv_cm2 = openmc.data.dose_coefficients(particle="photon", geometry="AP")
dose_filter = openmc.EnergyFunctionFilter(energies, pSv_cm2, interpolation="cubic")
photon_filter = openmc.ParticleFilter("photon")

dose_tally = openmc.Tally(name="dose tally")
dose_tally.filters = [dose_filter, photon_filter, cell_filter2]
dose_tally.scores = ["flux"]

# Use only dose tally for D1S run
model.tallies = [dose_tally]

# D1S settings
model.settings.photon_transport = True
model.settings.use_decay_photons = True

nuclides = d1s.prepare_tallies(model)
factors = d1s.time_correction_factors(nuclides, timesteps, source_rates)

print("---------------------------")
print("Performing D1S run")
print("---------------------------")

statepoint = model.run(output=False)

with openmc.StatePoint(statepoint) as sp:
    tally = sp.get_tally(name="dose tally")

corrected_tallies = [d1s.apply_time_correction(tally, factors, i + 1) for i in range(len(timesteps))]

print("Displaying cell dose rates")

profiles = []
for t_cool, ctally in zip(timesteps[1:], corrected_tallies[1:]):
    d1s_df = ctally.get_pandas_dataframe()

    if "cell" in d1s_df.columns:
        cell_col = "cell"
    elif "cell_id" in d1s_df.columns:
        cell_col = "cell_id"
    else:
        raise KeyError(f"Could not find a cell column in tally dataframe. Columns: {list(d1s_df.columns)}")

    d1s_df["cell_volume"] = d1s_df[cell_col].map(vol_by_cell_all)
    if d1s_df["cell_volume"].isna().any():
        missing = d1s_df.loc[d1s_df["cell_volume"].isna(), cell_col].unique().tolist()
        raise KeyError(f"Missing volumes for cells: {missing}")

    d1s_df["μSv/h"] = d1s_df["mean"] * (s_to_h * to_μSv) / d1s_df["cell_volume"]
    d1s_df["mSv/h"] = d1s_df["mean"] * (s_to_h * to_mSv) / d1s_df["cell_volume"]

    df_prof = d1s_df[d1s_df[cell_col].isin(dose_cells_all)].copy()
    df_prof["centers"] = df_prof[cell_col].map(center_by_cell_full)

    if df_prof["centers"].isna().any():
        missing = df_prof.loc[df_prof["centers"].isna(), cell_col].unique().tolist()
        raise KeyError(f"Missing centroids for OB cells: {missing}")

    df_prof = df_prof.sort_values("centers")
    profiles.append({"t_s": float(t_cool), "df": df_prof.copy()})

    print(f"Cooling {t_cool:.3e} s | {OB_KEY} rows={len(df_prof)}")

# Plot a few profiles
if len(profiles) == 0:
    raise RuntimeError("profiles is empty: OB rows were never captured.")

ncurves = min(10, len(profiles))
idxs = np.linspace(0, len(profiles) - 1, ncurves, dtype=int)

plt.figure()
for i in idxs:
    t_s_i = profiles[i]["t_s"]
    dfp = profiles[i]["df"]
    plt.plot(dfp["centers"], dfp["μSv/h"], label=f"{t_s_i:.1e} s")
plt.grid()
plt.yscale("log")
plt.ylim(1e-6, 1e12)
plt.ylabel("Shutdown Dose (μSv/h)")
plt.xlabel("Radial Position [cm]")
plt.title(f"D1S spatial profile: {OB_KEY}")
plt.legend()
plt.savefig(f"sdr_profile_{OB_KEY}.png", dpi=300)
plt.show()



# -----------------------------------------------------------------------------
# Depletion
# -----------------------------------------------------------------------------
print("--------------------------------")
print("Performing depletion")
print("--------------------------------")
# -----------------------------------------------------------------------------
# 1. Prepare Geometry & Materials
# -----------------------------------------------------------------------------
# Differentiate only what is necessary
model.differentiate_mats("match cell", depletable_only=True)

# Synchronize the model materials after differentiation
model.materials = openmc.Materials(model.geometry.get_all_materials().values())

# -----------------------------------------------------------------------------
# 2. Build Target Lists (Cells 55-66)
# -----------------------------------------------------------------------------
target_ids = sorted([int(cid) for cid in ob_cells]) # Ensure sorted order
deplete_mats = []
dagmc_cell_ids = [] # Store this alongside materials

all_cells = model.geometry.get_all_cells()

for cid in target_ids:
    cell = all_cells.get(cid)
    if cell and isinstance(cell.fill, openmc.Material):
        mat = cell.fill
        
        # Ensure volume is set
        if mat.volume is None:
            mat.volume = getattr(cell, 'volume', None)
            
        if mat.volume is None:
            print(f"[Warning] Cell {cid} has no volume assigned. Results will be 0.")

        mat.depletable = True
        
        deplete_mats.append(mat)
        dagmc_cell_ids.append(cid)

if not deplete_mats:
    raise RuntimeError(f"No materials found in cells {target_ids}. Check cell IDs.")

print(f"Targeting {len(deplete_mats)} cells/materials for depletion.")

# -----------------------------------------------------------------------------
# 3. Chain Reduction & MicroXS
# -----------------------------------------------------------------------------
initial_nuclides = model.geometry.get_all_nuclides()
reduced_chain = chain.reduce(initial_nuclides, level=2)
reduced_chain.export_to_xml("bluemira_chain.xml")

# Run MicroXS calculation
fluxes, micros = openmc.deplete.get_microxs_and_flux(
    model,
    deplete_mats, 
    chain_file="bluemira_chain.xml",
    run_kwargs={"output": True} # Threads=1 helps debug SEGV
)

# -----------------------------------------------------------------------------
# Operator + Integrator
# -----------------------------------------------------------------------------
operator = openmc.deplete.IndependentOperator(
    deplete_mats,
    fluxes,
    micros,
    chain_file="bluemira_chain.xml",
    normalization_mode="source-rate",
)
operator.output_dir = "r2s/activation"

integrator = openmc.deplete.PredictorIntegrator(
    operator,
    timesteps,
    source_rates=source_rates,
)
integrator.integrate()

# -----------------------------------------------------------------------------
# Post-processing: results + cell-material map
# -----------------------------------------------------------------------------
results = openmc.deplete.Results("r2s/activation/depletion_results.h5")

rows = []
for cid, mat in zip(dagmc_cell_ids, deplete_mats):
    rows.append({"cell_id": int(cid), "mat_id": int(mat.id), "material_name": mat.name})

Path("r2s/activation").mkdir(parents=True, exist_ok=True)
df_map = pd.DataFrame(rows)
df_map.to_csv("r2s/activation/cell_material_map.csv", index=False)
print("Written r2s/activation/cell_material_map.csv")

print("[done] depletion_model.py finished successfully")