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

from bluemira_full_setup import neutron_source_rate
from bluemira_full_setup import model, all_cells, xcentroids_ob, cell_ids, cell_filter

# Resolve chain path relative to this file (adjust folders as needed)
chain_path = (Path(__file__).resolve().parent.parent / "depletion_chain" / "chain_endfb71_sfr.xml").resolve()
# ^ change filename to one that actually exists on your machine

if not chain_path.exists():
    raise FileNotFoundError(f"Chain file not found: {chain_path}")

# Parse once (sanity check)
chain = openmc.deplete.Chain.from_xml(str(chain_path))

# Register for both Python-side and the OpenMC executable
openmc.config['chain_file'] = str(chain_path)
model.settings.depletion = {'chain_file': str(chain_path)}

# compute cell volumes (already set from DAGMC volume calc)
vol_by_cell = {cid: cell.volume for cid, cell in all_cells.items()}

print(f'Neutron source rate: {neutron_source_rate}')

# conversion factors
to_hours = 3600
to_μSv = 1e-6
to_mSv = 1e-9

irradiation_time = 1E7  # s
# cooling at ~12 d, 7.6 mo, 3.2 y
cooling_times = [1e6, 2e7, 1e8]  # s

timesteps = [irradiation_time] + cooling_times
source_rates = [neutron_source_rate] + [0.0] * len(cooling_times)

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

photon_filter = openmc.ParticleFilter('photon')
#dose_cell_filter = openmc.CellFilter(cell_ids)
#print("cell_filter bins:", dose_cell_filter.bins)

# Is the cell_filter(cell_ids) correct for dose_tally?
dose_tally = openmc.Tally(name='dose tally')
dose_tally.filters = [dose_filter, photon_filter, cell_filter]
dose_tally.scores = ['flux']

model.tallies = [dose_tally]

model.settings.photon_transport = True
model.settings.use_decay_photons = True


nuclides = d1s.prepare_tallies(model)
factors = d1s.time_correction_factors(nuclides,
                                      timesteps,
                                      source_rates)

print('---------------------------')
print(f'Performing D1S run')
print('---------------------------')

statepoint = model.run(output=False)

with openmc.StatePoint(statepoint) as sp:
    tally = sp.get_tally(name='dose tally')

corrected_tallies = []
for i, t in enumerate(timesteps):
    corrected_tally = d1s.apply_time_correction(tally, factors, i+1)
    corrected_tallies.append(corrected_tally)

x_by_cell = dict(zip(cell_ids, xcentroids_ob))
for t_cool, t_corr in zip(timesteps[1:], corrected_tallies):
    print('---------------------------')
    print(f'Cooling Time: {t_cool:.3e} s')
    print('---------------------------')

    d1s_df = t_corr.get_pandas_dataframe()

    # Attach volume and dose rate only to the actually-present cells
    d1s_df['volume']  = d1s_df['cell'].map(vol_by_cell)
    d1s_df['μSv/h']   = d1s_df['mean'] * (to_hours * to_μSv) / d1s_df['volume']
    d1s_df['mSv/h']   = d1s_df['mean'] * (to_hours * to_mSv) / d1s_df['volume']
    d1s_df['centers'] = d1s_df['cell'].map(x_by_cell)
    print(d1s_df)
    


# plot the D1S results
plt.plot(d1s_df['centers'], d1s_df['μSv/h'], label='Direct One-Step')
plt.grid()
plt.ylabel('Shutdown Dose - d1s (μSv/h)')
plt.xlabel('Radial Position [cm]')
plt.legend()
plt.savefig('sdr_d1s.png')
plt.close()

# ------------------------------------------------------------------
# R2S / depletion setup
# ------------------------------------------------------------------
print('--------------------------------')
#print('Performing R2S Activation run')
print('Performing depletion')
print('--------------------------------')

model.settings.photon_transport = False
model.settings.use_decay_photons = False

# place unique materials in each cell, to be activated individually
model.differentiate_mats('match cell', depletable_only=False)

# Prepare list cells_ids/mat for post processing of depletion
dagmc_cell_ids = []
deplete_cells = []
deplete_mats = []

for cid in sorted(all_cells.keys()):
    cell = all_cells[cid]
    mat = cell.fill

    # Skip void / universes / anything that is not a Material
    if not isinstance(mat, openmc.Material):
        continue

    dagmc_cell_ids.append(cid)
    deplete_cells.append(cell)

    # Ensure material is depletable and has a volume
    if mat.volume is None:
        mat.volume = cell.volume
    mat.depletable = True

    deplete_mats.append(mat)

# include original tallies during activation and cooling
model.tallies = orig_tallies

# Load the depletion chain
# Load all nuclides in the geometry
initial_nuclides = list(model.geometry.get_all_nuclides())

# Using this to avoid crash with C0 not in the chain_file
# Should be easier than I think to fix this
chain_nuclides = set(chain.nuclide_dict)
valid = sorted(set(initial_nuclides) & chain_nuclides)
missing = sorted(set(initial_nuclides) - chain_nuclides)

if missing:
    print("Ignoring nuclides not in chain:", ", ".join(missing))

# REDUCED_CHAIN
reduced_chain = chain.reduce(valid, level=5)
reduced_chain.export_to_xml("bluemira_chain.xml")
bluemira_chain = Path("bluemira_chain.xml").resolve()

fluxes, micros = openmc.deplete.get_microxs_and_flux(
    model,
    deplete_mats,                 # materials as domains
    chain_file=bluemira_chain,
    run_kwargs={'output': False},
)

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
integrator.integrate(final_step=False, output=False)

###########################################################
# Alternative: model.deplete(...) (commented out for now)
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
# DECAY HEAT, ACTIVATION, MATERIAL TRANSMUTATION
# -----------------------------------------------------------------
# Map cell id -> depletion index
index_by_cid = {cid: i for i, cid in enumerate(dagmc_cell_ids)}

# ------------------------------------------------------------------
# Decay Heat 
# ------------------------------------------------------------------

def get_cell_decay_heat_profile(results,
                                cell_ids,
                                all_cells,
                                xcentroids,
                                time_index=-1,
                                units='W/cm3'):
    """
    Radial profile of decay heat per cell.

    Returns
    -------
    times : np.ndarray
    x : np.ndarray
    heat : np.ndarray
    """
    x = []
    heat = []

    for i, cid in enumerate(cell_ids):
        cell = all_cells[cid]
        mat = cell.fill         
        vol = cell.volume      

        times, heat_t = results.get_decay_heat(
            mat=mat,           
            units=units,
            by_nuclide=False,
            volume=vol         
        )

        x.append(xcentroids[i])
        heat.append(heat_t[time_index])

    return np.array(times), np.array(x), np.array(heat)

# Last timestep, volumetric heat
times, x, heat = get_cell_decay_heat_profile(
    results,
    cell_ids=cell_ids,
    all_cells=all_cells,
    xcentroids=xcentroids_ob,
    time_index=-1,
    units='W/cm3'
)

plt.figure()
plt.semilogy(x, heat, marker='o')
plt.xlabel('Radial position [cm]')
plt.ylabel('Decay heat [W/cm$^3$]')
plt.grid(True, which='both')
plt.tight_layout()
plt.savefig('radial_decay_heat.png')
plt.close()

# ------------------------------------------------------------------
# Activity
# ------------------------------------------------------------------

def get_cell_activity_profile(results,
                              cell_ids,
                              all_cells,
                              xcentroids,
                              time_index=-1,
                              units='Bq/cm3'):
    """
    Radial profile of activity per cell.

    Returns
    -------
    times : np.ndarray
    x : np.ndarray
    activity : np.ndarray
    """
    x = []
    act = []

    for i, cid in enumerate(cell_ids):
        cell = all_cells[cid]
        mat = cell.fill
        vol = cell.volume

        times, act_t = results.get_activity(
            mat=mat,
            units=units,       
            by_nuclide=False,
            volume=vol
        )

        x.append(xcentroids[i])
        act.append(act_t[time_index])

    return np.array(times), np.array(x), np.array(act)

times, x, activity = get_cell_activity_profile(
    results,
    cell_ids=cell_ids,
    all_cells=all_cells,
    xcentroids=xcentroids_ob,
    time_index=-1,
    units='Bq/cm3'
)

plt.figure()
plt.semilogy(x, activity, marker='s')
plt.xlabel('Radial position [cm]')
plt.ylabel('Activity [Bq/cm$^3$]')
plt.grid(True, which='both')
plt.tight_layout()
plt.savefig('radial_activity.png')
plt.close()

# ------------------------------------------------------------------
# Transmutation list -> CSV
# ------------------------------------------------------------------

def export_cell_compositions_to_csv(results,
                                    burn_step_index,
                                    cell_ids,
                                    dagmc_cell_ids,
                                    deplete_mats,
                                    xcentroids,
                                    fname_prefix="cell_composition"):
    """
    Export depleted nuclide compositions per cell into a CSV using
    the depletion results directly (results.get_atoms).

    Returns
    -------
    df : pandas.DataFrame
    """
    # Resolve negative index (e.g. -1 -> last step)
    n_steps = len(results)
    if burn_step_index < 0:
        step_idx = n_steps + burn_step_index
    else:
        step_idx = burn_step_index

    if step_idx < 0 or step_idx >= n_steps:
        raise IndexError(f"burn_step_index {burn_step_index} is out of range for {n_steps} depletion steps")

    step0 = results[0]
    nuc_list = list(step0.index_nuc.keys())

    # Map cell id -> index into deplete_mats
    index_by_cid = {cid: i for i, cid in enumerate(dagmc_cell_ids)}

    rows = []

    for i, cid in enumerate(cell_ids):
        if cid not in index_by_cid:
            print(f"Skipping cell {cid}: not in depletable cell list")
            continue

        mat = deplete_mats[index_by_cid[cid]]

        # Loop over all nuclides and grab atoms at this step
        for nuc in nuc_list:
            times, atoms = results.get_atoms(mat, nuc,
                                             nuc_units='atoms',
                                             time_units='s')
            n_at_step = atoms[step_idx]
            
            # ADD FILTER TO VERY LOW MATERIAL

            if n_at_step == 0.0:
                continue

            rows.append({
                "cell_id": cid,
                "xcentroid_cm": xcentroids[i],
                "material_id": mat.id,
                "nuclide": nuc,
                "atoms": n_at_step
            })

    df = pd.DataFrame(rows)
    fname = f"{fname_prefix}_step{burn_step_index}.csv"
    df.to_csv(fname, index=False)
    print(f"Written {fname} with {len(df)} rows")
    return df


comp_df = export_cell_compositions_to_csv(
    results,
    burn_step_index=-1,
    cell_ids=cell_ids,
    dagmc_cell_ids=dagmc_cell_ids,
    deplete_mats=deplete_mats,
    xcentroids=xcentroids_ob,
    fname_prefix="ob_cells"
)

