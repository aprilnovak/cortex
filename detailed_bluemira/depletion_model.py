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
module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'materials'))
sys.path.append(module_path)
import materials

from neutronics_model import neutron_source_rate
from neutronics_model import model, pydagmc_model, all_cells, xcentroids_ob, cell_ids, cell_filter

# Resolve chain file (Using ENDF/B-VIII.0)
chain_path = (Path(__file__).resolve().parent.parent / "depletion_chain" / "chain_endfb80_sfr.xml").resolve()

if not chain_path.exists():
    raise FileNotFoundError(f"Chain file not found: {chain_path}")

# Parse once (sanity check)
chain = openmc.deplete.Chain.from_xml(str(chain_path))

# Register for both Python-side and the OpenMC executable
openmc.config['chain_file'] = str(chain_path)
model.settings.depletion = {'chain_file': str(chain_path)}

# Load cell volumes (already set from pydagmc)
vol_by_cell = {cid: cell.volume for cid, cell in all_cells.items()}

#print(f'Neutron source rate: {neutron_source_rate}')
# conversion factors
s_to_h = 3600
y_to_s = 24 * 365 * s_to_h
to_μSv = 1e-6
to_mSv = 1e-9

#irradiation_time (5 years)
irradiation_time_y = np.array([1.0, 1.0, 1.0, 1.0, 1.0])
irradiation_time = (irradiation_time_y * y_to_s).tolist()
# cooling  (1e-9 y to 1000 y)
timesteps_years = np.logspace(-9, 3, num=13)
cooling_times = (timesteps_years * y_to_s).tolist()

timesteps = irradiation_time + cooling_times
source_rates = [neutron_source_rate] * len(irradiation_time) + [0.0] * len(cooling_times)

# save the current set of tallies to be re-applied
# during R2S calculations
orig_tallies = [t for t in model.tallies]

energies, pSv_cm2 = openmc.data.dose_coefficients(
    particle='photon', geometry='AP'
)

# cubic interpolation recommended by ICRP
dose_filter = openmc.EnergyFunctionFilter(
    energies, pSv_cm2, interpolation='cubic'
)

# ------------------------------------------------------------------
# D1S
# ------------------------------------------------------------------
# Adding D1S as soon as DAGMC model is available


# ------------------------------------------------------------------
# Depletion -> *change that for R2S later
# ------------------------------------------------------------------
# THIS MIGHT CHANGE WHEN WE INTRODUCE R2S
# (Can we optimize to run less transport?)
print('--------------------------------')
print('Performing depletion')
print('--------------------------------')

model.settings.photon_transport = False
model.settings.use_decay_photons = False

# place unique materials in each cell, to be activated individually
model.differentiate_mats('match cell', depletable_only=False)

# include original tallies during activation and cooling
model.tallies = orig_tallies

# Prepare list cells_ids/mat for post processing of depletion
# (materials ids will change with depletion)
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

    deplete_mats.append(mat)

# Mapping for post-processing
# Since post-processing is in another file, I am saving the mapping information
# damgc_cell_ids -> deplete_mats_id in a CSV file.
# I would appreciate smarter/more efficient suggestions 
cell_to_mat = {
    cid: str(mat.id)
    for cid, mat in zip(dagmc_cell_ids, deplete_mats)
}

mat_to_cell = {
    str(mat.id): cid
    for cid, mat in zip(dagmc_cell_ids, deplete_mats)
}

mat_id_to_name = {
    str(mat.id): (mat.name or f"material_{mat.id}")
    for mat in deplete_mats
}

# REDUCED_CHAIN
initial_nuclides = model.geometry.get_all_nuclides()
reduced_chain = chain.reduce(initial_nuclides, level=5) # previous simulation was on 4
reduced_chain.export_to_xml("bluemira_chain.xml")
bluemira_chain = Path("bluemira_chain.xml").resolve()

# compute nuclide fluxes and microscopic cross-sections
fluxes, micros = openmc.deplete.get_microxs_and_flux(
    model,
    deplete_mats,                 
    chain_file=bluemira_chain,
    run_kwargs={'output': False},
)

# perform neutron activation and produces a depletion_results.h5 file
# As discussed, using IndepentOperator
operator = openmc.deplete.IndependentOperator(
    deplete_mats,
    fluxes,
    micros,
    chain_file=bluemira_chain,
    normalization_mode='source-rate',
)
operator.output_dir = 'r2s/activation'

integrator = openmc.deplete.PredictorIntegrator(
    operator,
    timesteps,
    source_rates=source_rates
)
integrator.integrate()

###########################################################
# Alternative (commented out for now)
###########################################################
# model.deplete(
#     timesteps,
#     source_rates=source_rates,
#     output=False,
#     method="cf4",
#     directory="r2s/activation",
#     final_step=False,
#     operator_kwargs={
#         "normalization_mode": "source-rate",
#         "chain_file": bluemira_chain,
#         "reduce_chain_level": 5,
#         "reduce_chain": True,
#     },
# )

# ------------------------------------------------------------------
# Post-processing: Results
# ------------------------------------------------------------------
results = openmc.deplete.Results("r2s/activation/depletion_results.h5")

# -----------------------------------------------------------------
# Introduce  R2S Decay Gamma run
# -----------------------------------------------------------------
# MISSING TO DEFINE A SPATIAL REGION FOR GAMMA SOURCE


# -----------------------------------------------------------------
# EXPORT CELL-MAT MAPPING FOR POST PROCESSING
# -----------------------------------------------------------------
# Save mapping cell_id -> mat_id during depletion for post-processing
rows = []
for cid, mat in zip(dagmc_cell_ids, deplete_mats):
    rows.append({
        "cell_id": cid,
        "mat_id": int(mat.id),
        "material_name": mat.name,
    })

df_map = pd.DataFrame(rows)
df_map.to_csv("r2s/activation/cell_material_map.csv", index=False)
print("Written r2s/activation/cell_material_map.csv")

