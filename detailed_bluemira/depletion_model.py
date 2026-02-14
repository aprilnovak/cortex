from pathlib import Path

from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import openmc
import openmc.data
import openmc.deplete
from openmc.deplete import d1s

import sys
import os
module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "materials"))
sys.path.append(module_path)
import materials 

# ------------------------------------------------------------------
# Import from neutronics_model.py
# ------------------------------------------------------------------
from neutronics_model import (
    neutron_source_rate,
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

chunk_cells = build_breeder_chunks(INPUT_JSON, default_equatorial_ob_key=OB_KEY)

ob_by_key = chunk_cells["ob_by_key"]
ib_by_key = chunk_cells["ib_by_key"]
OB_CHUNK_SIZE = chunk_cells["OB_CHUNK_SIZE"]
IB_CHUNK_SIZE = chunk_cells["IB_CHUNK_SIZE"]

cell_ids = chunk_cells["cell_ids_all"]
cell_ids_equatorial_ob = chunk_cells["equatorial_ob_cell_ids"]

# ------------------------------------------------------------------
# Helper: bounding box
# ------------------------------------------------------------------
def dagmc_bounding_box(pydagmc_model, volume_id):
    """Returns the bounding box of a given volume in a PyDAGMC model."""
    volume = pydagmc_model.volumes_by_id[volume_id]
    triangle_coords = volume.triangle_coords
    min_coords = np.min(triangle_coords, axis=0)
    max_coords = np.max(triangle_coords, axis=0)
    return openmc.BoundingBox(min_coords, max_coords)


# ------------------------------------------------------------------
# Resolve chain file (ENDF/B-VIII.0)
# ------------------------------------------------------------------
chain_path = (Path(__file__).resolve().parent.parent / "depletion_chain" / "chain_endfb80_sfr.xml").resolve()
if not chain_path.exists():
    raise FileNotFoundError(f"Chain file not found: {chain_path}")

# Parse once (sanity check)
chain = openmc.deplete.Chain.from_xml(str(chain_path))

# Register for both Python-side and the OpenMC executable
openmc.config["chain_file"] = str(chain_path)
model.settings.depletion = {"chain_file": str(chain_path)}

# ------------------------------------------------------------------
# Time grids / source rates
# ------------------------------------------------------------------
s_to_h = 3600
y_to_s = 24 * 365 * s_to_h
to_μSv = 1e-6
to_mSv = 1e-9

# irradiation (5 years)
irradiation_time_y = np.array([5.0])
irradiation_time = (irradiation_time_y * y_to_s).tolist()

# cooling (1e-8 y to 1000 y)
timesteps_years = np.concatenate([
    np.logspace(-8, -4, 6),
    np.logspace(-4,  0, 8)[1:],  # drop 1e-4
    np.logspace( 0,  3, 14)[1:],  # drop 1e1
])
cooling_times = (timesteps_years * y_to_s).tolist()

timesteps = irradiation_time + cooling_times
constant_power_ratio = 0.3
source_rates = [constant_power_ratio * neutron_source_rate] * len(irradiation_time) + [0.0] * len(cooling_times)

# ------------------------------------------------------------------
# Volumes from DAGMC cells
# ------------------------------------------------------------------
vol_by_cell = {cid: cell.volume for cid, cell in all_cells.items()}

# TO BE UPDATED:
xcentroids_ob = (0.1, 1.1, 6.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 84.4, 156.0)
xcentroids_ib = (0.1, 1.1, 6.0, 15.0, 25.0, 35.0, 45.0, 55.0, 65.0, 73.8, 108.0)
# ------------------------------------------------------------------
# Identify plasma and vv port-fill cells (last two cells)
# ------------------------------------------------------------------
vols = list(pydagmc_model.volumes)

vols = list(pydagmc_model.volumes)
if len(vols) < 2:
    raise RuntimeError(f"PyDAGMC model has only {len(vols)} volumes; expected >= 2.")

plasma_vol_id = vols[-2].id
vvportfill_vol_id = vols[-1].id
print("vol_ids for sdr= ", plasma_vol_id, "and", vvportfill_vol_id)

# ------------------------------------------------------------------
# Add OB_1_b6 to D1S tally set
# ------------------------------------------------------------------
OB_KEY = "OB_1_b6"
if OB_KEY not in ob_by_key:
    raise KeyError(f"{OB_KEY} not found in ob_by_key. Available (sample): {list(ob_by_key)[:10]}")

ob_1_b6_cells = ob_by_key[OB_KEY]

# xcentroids_ob is per-layer for one OB chunk: [Armor, FW, OB1..OBn, VV]
if len(xcentroids_ob) != OB_CHUNK_SIZE:
    raise ValueError(
        f"len(xcentroids_ob)={len(xcentroids_ob)} but OB_CHUNK_SIZE={OB_CHUNK_SIZE}. "
        "xcentroids_ob must have one entry per item in the OB chunk "
        "(Armor, FW, OB layers..., VV)."
    )
if len(ob_1_b6_cells) != OB_CHUNK_SIZE:
    raise ValueError(
        f"len(ob_1_b6_cells)={len(ob_1_b6_cells)} but OB_CHUNK_SIZE={OB_CHUNK_SIZE}. "
        "OB_1_b6 chunk size mismatch."
    )

# Map OB_1_b6 cell_id -> radial center [cm]
center_by_cell_ob6 = dict(zip(ob_1_b6_cells, xcentroids_ob))

# Cells included in the dose tally
dose_cells_time = [plasma_vol_id, vvportfill_vol_id]      # time series
dose_cells_prof = list(ob_1_b6_cells)                     # spatial profile
dose_cells_all = dose_cells_time + dose_cells_prof

cell_filter2 = openmc.CellFilter(dose_cells_all)

# Per-cell volumes for normalization for all dose cells
vol_by_cell_all = {}
for cid in dose_cells_all:
    v = vol_by_cell.get(cid, None)
    if v is None:
        raise KeyError(f"Missing volume for cell {cid} in vol_by_cell (all_cells).")
    vol_by_cell_all[cid] = float(v)

# Save the current set of tallies to re-apply during activation/cooling (R2S later)
orig_tallies = [t for t in model.tallies]

# ------------------------------------------------------------------
# Dose coefficients and D1S tally
# ------------------------------------------------------------------
energies, pSv_cm2 = openmc.data.dose_coefficients(particle="photon", geometry="AP")
dose_filter = openmc.EnergyFunctionFilter(energies, pSv_cm2, interpolation="cubic")
photon_filter = openmc.ParticleFilter("photon")

dose_tally = openmc.Tally(name="dose tally")
dose_tally.filters = [dose_filter, photon_filter, cell_filter2]
dose_tally.scores = ["flux"]

model.tallies = [dose_tally]

# ------------------------------------------------------------------
# D1S run
# ------------------------------------------------------------------
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

# Apply time correction for each timestep index
corrected_tallies = []
for i in range(len(timesteps)):
    corrected_tallies.append(d1s.apply_time_correction(tally, factors, i + 1))

print("Displaying cell dose rates")

# Storage:
time_s = []  # cooling times (seconds), corresponds to timesteps[1:]
dose_time_by_cell = {plasma_vol_id: [], vvportfill_vol_id: []}
profiles = []  # list of {"t_s":..., "df":...} for OB profile

for t_cool, ctally in zip(timesteps[1:], corrected_tallies[1:]):
    d1s_df = ctally.get_pandas_dataframe()

    # identify cell column name
    if "cell" in d1s_df.columns:
        cell_col = "cell"
    elif "cell_id" in d1s_df.columns:
        cell_col = "cell_id"
    else:
        raise KeyError(f"Could not find a cell column in tally dataframe. Columns: {list(d1s_df.columns)}")

    # normalize per cell
    d1s_df["cell_volume"] = d1s_df[cell_col].map(vol_by_cell_all)
    if d1s_df["cell_volume"].isna().any():
        missing = d1s_df.loc[d1s_df["cell_volume"].isna(), cell_col].unique().tolist()
        raise KeyError(f"Missing volumes for cells: {missing}")

    d1s_df["μSv/h"] = d1s_df["mean"] * (s_to_h * to_μSv) / d1s_df["cell_volume"]
    d1s_df["mSv/h"] = d1s_df["mean"] * (s_to_h * to_mSv) / d1s_df["cell_volume"]

    # -------------------------
    # (A) time series: plasma & vvportfill
    # -------------------------
    df_time = d1s_df[d1s_df[cell_col].isin(dose_cells_time)].copy()
    got = set(df_time[cell_col].tolist())
    want = set(dose_cells_time)
    if got != want:
        raise RuntimeError(f"Time-series rows mismatch. Got {sorted(got)}, want {sorted(want)}")

    time_s.append(float(t_cool))
    for cid in dose_cells_time:
        dose_time_by_cell[cid].append(float(df_time.loc[df_time[cell_col] == cid, "μSv/h"].iloc[0]))

    # -------------------------
    # (B) spatial profile: OB_1_b6
    # -------------------------
    df_prof = d1s_df[d1s_df[cell_col].isin(ob_1_b6_cells)].copy()
    df_prof["centers"] = df_prof[cell_col].map(center_by_cell_ob6)

    if df_prof["centers"].isna().any():
        missing = df_prof.loc[df_prof["centers"].isna(), cell_col].unique().tolist()
        raise KeyError(f"Missing centroids for OB_1_b6 cells: {missing}")

    df_prof = df_prof.sort_values("centers")
    profiles.append({"t_s": float(t_cool), "df": df_prof.copy()})

    print(
        f"Cooling {t_cool:.3e} s | plasma={dose_time_by_cell[plasma_vol_id][-1]:.3e} μSv/h | "
        f"vvpf={dose_time_by_cell[vvportfill_vol_id][-1]:.3e} μSv/h | "
        f"{OB_KEY} rows={len(df_prof)}"
    )

# ------------------------------------------------------------------
# Plot 1A: time series for plasma
# ------------------------------------------------------------------
plt.figure()
plt.plot(time_s, dose_time_by_cell[plasma_vol_id], label="Plasma (D1S)")
plt.grid()
plt.xscale("log")
plt.yscale("log")
plt.ylabel("Shutdown Dose (μSv/h)")
plt.xlabel("Cooling Time [s]")
plt.legend()
plt.savefig("sdr_time_plasma.png", dpi=300, bbox_inches="tight")
plt.show()

# ------------------------------------------------------------------
# Plot 1B: time series for VV port-fill
# ------------------------------------------------------------------
plt.figure()
plt.plot(time_s, dose_time_by_cell[vvportfill_vol_id], label="VV port-fill (D1S)")
plt.grid()
plt.xscale("log")
plt.yscale("log")
plt.ylabel("Shutdown Dose (μSv/h)")
plt.xlabel("Cooling Time [s]")
plt.legend()
plt.savefig("sdr_time_vvportfill.png", dpi=300, bbox_inches="tight")
plt.show()

# ------------------------------------------------------------------
# Plot 2: spatial profiles for OB_1_b6 at a few representative cooling times
# ------------------------------------------------------------------
if len(profiles) == 0:
    raise RuntimeError("profiles is empty: OB_1_b6 rows were never captured.")

ncurves = min(10, len(profiles))
idxs = np.linspace(0, len(profiles) - 1, ncurves, dtype=int)

plt.figure()
for i in idxs:
    t_s_i = profiles[i]["t_s"]
    dfp = profiles[i]["df"]
    plt.plot(dfp["centers"], dfp["μSv/h"], label=f"{t_s_i:.1e} s")

plt.grid()
plt.yscale("log")
plt.ylim(1e-6,1e12)
plt.ylabel("Shutdown Dose (μSv/h)")
plt.xlabel("Radial Position [cm]")
plt.title(f"D1S spatial profile: {OB_KEY}")
plt.legend()
plt.savefig(f"sdr_profile_{OB_KEY}.png", dpi=300)
plt.show()

# ------------------------------------------------------------------
# Depletion 
# ------------------------------------------------------------------
print("--------------------------------")
print("Performing depletion")
print("--------------------------------")

# Make sure the model's materials list matches the geometry
model.materials = openmc.Materials(list(model.geometry.get_all_materials().values()))
model.geometry.determine_paths()

model.settings.photon_transport = False
model.settings.use_decay_photons = False

# place unique materials in each cell, to be activated individually
model.differentiate_mats("match cell", depletable_only=False)

# include original tallies during activation and cooling
model.tallies = orig_tallies

# Prepare list cells_ids/mat for post processing of depletion
cells = list(model.geometry.get_all_cells().values())
dagmc_cell_ids = []
deplete_cells = []
deplete_mats = []

for cid in sorted(all_cells.keys()):
    cell = all_cells[cid]
    mat = cell.fill

    # Skip void or universes that is not a Material (complementary volume)
    if not isinstance(mat, openmc.Material):
        continue

    dagmc_cell_ids.append(cid)
    deplete_cells.append(cell)

    # Ensure material is depletable and has a volume
    if mat.volume is None:
        mat.volume = cell.volume
    mat.depletable = True
    # MAYBE REMOVE PLASMA_REGION AND VV_port_Fill as depletable

    deplete_mats.append(mat)

# Mapping for post-processing 
cell_to_mat = {cid: str(mat.id) for cid, mat in zip(dagmc_cell_ids, deplete_mats)}
mat_to_cell = {str(mat.id): cid for cid, mat in zip(dagmc_cell_ids, deplete_mats)}
mat_id_to_name = {str(mat.id): (mat.name or f"material_{mat.id}") for mat in deplete_mats}

# ---- Build + use reduced chain everywhere below ----
initial_nuclides = model.geometry.get_all_nuclides()
reduced_chain = chain.reduce(initial_nuclides, level=2)

bluemira_chain = Path("bluemira_chain.xml").resolve()
reduced_chain.export_to_xml(str(bluemira_chain))
print(f"[info] Wrote reduced chain: {bluemira_chain}")

openmc.config["chain_file"] = str(bluemira_chain)
model.settings.depletion = {"chain_file": str(bluemira_chain)}

fluxes, micros = openmc.deplete.get_microxs_and_flux(
    model,
    deplete_mats,
    chain_file=str(bluemira_chain),
    run_kwargs={"output": False},
)

operator = openmc.deplete.IndependentOperator(
    deplete_mats,
    fluxes,
    micros,
    chain_file=str(bluemira_chain),
    normalization_mode="source-rate",
)
operator.output_dir = "r2s/activation"

integrator = openmc.deplete.PredictorIntegrator(
    operator,
    timesteps,
    source_rates=source_rates,
)
integrator.integrate()

# ------------------------------------------------------------------
# Post-processing: Results
# ------------------------------------------------------------------
results = openmc.deplete.Results("r2s/activation/depletion_results.h5")

# -----------------------------------------------------------------
# EXPORT CELL-MAT MAPPING FOR POST PROCESSING
# -----------------------------------------------------------------
rows = []
for cid, mat in zip(dagmc_cell_ids, deplete_mats):
    rows.append({"cell_id": cid, "mat_id": int(mat.id), "material_name": mat.name})

df_map = pd.DataFrame(rows)
df_map.to_csv("r2s/activation/cell_material_map.csv", index=False)
print("Written r2s/activation/cell_material_map.csv")

# -----------------------------------------------------------------
# Introduce R2S decay-gamma run (placeholder)
# -----------------------------------------------------------------
# 
