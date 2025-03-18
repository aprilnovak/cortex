from pathlib import Path

from matplotlib import pyplot as plt

import openmc
from openmc.deplete import d1s

from tile import neutron_source_rate
from tile import model, cell_filter, frontal_side, cell_volume, xcentroids

# compute cell volumes
bbox = model.geometry.bounding_box
all_cells = list(model.geometry.get_all_cells().values())
cell_vc = openmc.VolumeCalculation(all_cells, 10_000_000, *bbox)
model.settings.volume_calculations = [cell_vc]
model.calculate_volumes(apply_volumes=True, output=False)

# all cell volumes are equal in this case
# cell_volume = all_cells[0].volume


print(f'Neutron source rate: {neutron_source_rate}')
# taken from the ITER-SA2 irradiation scenario
#

# conversion factors
to_hours = 3600
to_μSv = 1e-6
to_mSv = 1e-9

irradiation_time = 1E7 # s

# cooling at ~12 d, 7.6 mo, 3.2 y
cooling_times = [1e6, 2e7, 1e8] # s

timesteps = [irradiation_time] + cooling_times
source_rates = [neutron_source_rate] + [0.0] * len(cooling_times)

# Tallies

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

dose_tally = openmc.Tally(name='dose tally')
dose_tally.filters = [dose_filter, photon_filter, cell_filter]
dose_tally.scores = ['flux']

# adding a mesh tally for visualization
mesh = openmc.RegularMesh.from_domain(model.geometry.root_universe,
                                      dimension=(5, 5, 5))

mesh_filter = openmc.MeshFilter(mesh)
neutron_filter = openmc.ParticleFilter('neutron')

mesh_tally = openmc.Tally()
mesh_tally.filters = [mesh_filter, neutron_filter]
mesh_tally.scores = ['flux']

model.settings.photon_transport = True
model.settings.use_decay_photons = True

model.tallies = [dose_tally, mesh_tally]

# this should only modify the dose tally b/c
# it has a particle filter with one 'photon' bin
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

print('Displaying cell dose rates')
for t_cool, t in zip(timesteps[1:], corrected_tallies):
    print('---------------------------')
    print(f'Cooling Time: {t_cool:.3e} s')
    print('---------------------------')
    d1s_df = t.get_pandas_dataframe()
    d1s_df['μSv/h'] = d1s_df['mean'] * (to_hours * to_μSv) / cell_volume
    d1s_df['mSv/h'] = d1s_df['mean'] * (to_hours * to_mSv) / cell_volume
    d1s_df['centers'] = xcentroids
    print(d1s_df)

# R2S setup
print('--------------------------------')
print(f'Performing R2S Activation run')
print('--------------------------------')

model.settings.photon_transport = False
model.settings.use_decay_photons = False

# perform a volume calculation for
bbox = model.geometry.bounding_box

# place unique materials in each cell, to be activated individually
model.differentiate_mats('match cell', depletable_only=False)

# mark all materials as depletable for now
for m in model.geometry.get_all_materials().values():
    m.depletable = True

model.materials = []

# include original tallies during
# activation and cooling
model.tallies = orig_tallies

# load the depletion chain
chain = openmc.deplete.Chain.from_xml(openmc.config['chain_file'])
initial_nuclides = model.geometry.get_all_nuclides()
reduced_chain = chain.reduce(initial_nuclides, level=5)
reduced_chain.export_to_xml('tungsten_chain.xml')

tungsten_chain = Path('tungsten_chain.xml').resolve()


cells = list(model.geometry.get_all_cells().values())
cell_volume = cells[0].volume

# compute nuclide fluxes and microscopic cross-sections
fluxes, micros = openmc.deplete.get_microxs_and_flux(model,
                                                     cells,
                                                     chain_file=tungsten_chain,
                                                     run_kwargs={'output': False})

# perform neutron activation and produces a depletion_results.h5 file
activation_materials = list(model.geometry.get_all_materials().values())
operator = openmc.deplete.IndependentOperator(activation_materials,
                                              fluxes,
                                              micros,
                                              chain_file=tungsten_chain,
                                              normalization_mode='source-rate')
operator.output_dir = 'r2s/activation'

# integrator = openmc.deplete.PredictorIntegrator(operator,
#                                                 timesteps,
#                                                 source_rates=source_rates)
# integrator.integrate(final_step=False, output=False)


###########################################################
# could also do this. It's simpler, but less efficient
# since it involves additional transport steps
###########################################################
model.deplete(
    timesteps,
    source_rates=source_rates,
    output=False,
    # method="predictor",  # predictor is a simple but quick method
    method="cf4",  # CF4Integrator is an accurate but slower method
    directory="r2s/activation",
    final_step=False,
    operator_kwargs={
        "normalization_mode": "source-rate",  # needed as this is a fixed source simulation
        "chain_file": tungsten_chain,
        "reduce_chain_level": 5,
        "reduce_chain": True,
    },
)

# load depletion/activation results
results = openmc.deplete.Results("r2s/activation/depletion_results.h5")

print('--------------------------------')
print(f'Performing R2S Decay Gamma run')
print('--------------------------------')

# update the dose filter
dose_tally.filters = [dose_filter, photon_filter, cell_filter]
model.tallies = [dose_tally]
model.settings.photon_transport = True

# generate decay photon sources for each cooling time
for i_cool, (t_cool, src_rate) in enumerate(zip(timesteps, source_rates)):

    if src_rate != 0:
        continue

    step_photon_sources = []

    for cell in model.geometry.get_all_material_cells().values():
        material = cell.fill
        if not material.depletable:
            continue

        # retrieve activated material composition from depletion results
        actiavted_material = results[i_cool].get_material(str(material.id))

        energy = actiavted_material.get_decay_photon_energy(units='Bq')

        strength = energy.integral()

        if strength == 0:
            continue

        space = openmc.stats.Box(*cell.bounding_box)
        source = openmc.IndependentSource(
            space=space,
            energy=energy,
            particle="photon",
            constraints= {'domains': [cell]},
            strength=strength
        )

        step_photon_sources.append(source)

    # TODO: should technically also use activated material compositions in transport here
    model.settings.source = step_photon_sources

    statepoint = model.run(cwd=f'r2s/cooling_{t_cool:.3e}_s', output=False)

    with openmc.StatePoint(statepoint) as sp:
        sp_dose = sp.get_tally(name=dose_tally.name)
        r2s_df = sp_dose.get_pandas_dataframe()

    # adjust units of dose
    print('---------------------------')
    print(f'Cooling Time: {t_cool:.3e} s')
    print('---------------------------')
    r2s_df['μSv/h'] = r2s_df['mean'] * (to_hours * to_μSv) / cell_volume
    r2s_df['mSv/h'] = r2s_df['mean'] * (to_hours * to_mSv) / cell_volume
    r2s_df['centers'] = xcentroids
    print(r2s_df)

# plot the D1S and R2S results on the same plot

plt.plot(d1s_df['centers'], d1s_df['μSv/h'], label='Direct One-Step')
plt.plot(r2s_df['centers'], r2s_df['μSv/h'], label='Rigorous Two-Step')
plt.grid()
plt.ylabel('Shutdown Dose (μSv/h)')
plt.xlabel('Radial Position [cm]')
plt.legend()
plt.savefig('sdr.png')
plt.show()