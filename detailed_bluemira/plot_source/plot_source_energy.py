import numpy as np
import openmc
import openmc.lib
from tempfile import TemporaryDirectory
import matplotlib
import matplotlib.pyplot as plt
import math

matplotlib.rcParams.update({'font.size': 12})

# Plots the NEUTRON plasma source in tandem with the surface source for comparison

# Number of samples to use when sampling the source
n_samples = 10000000

from openmc_plasma_source import tokamak_source

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

def sample_initial_particles(this, prn_seed: int = None):
    """smaples particles from the source.

    Args:
        this: The openmc source, settings or model containing the source to plot
        n_samples: The number of source samples to obtain.
        prn_seed: The pseudorandom number seed.
    """
    with TemporaryDirectory() as tmpdir:
        model = openmc.Model()

        materials = openmc.Materials()
        model.materials = materials

        sph = openmc.Sphere(r=99999999999, boundary_type="vacuum")
        cell = openmc.Cell(region=-sph)
        geometry = openmc.Geometry([cell])
        model.geometry = geometry

        settings = openmc.Settings()
        settings.particles = 1
        settings.batches = 1
        settings.source = this
        model.settings = settings

        model.export_to_model_xml()

        openmc.lib.init(output=False)
        particles = openmc.lib.sample_external_source(
            n_samples=n_samples, prn_seed=prn_seed
        )
        openmc.lib.finalize()

        p = []
        for i in particles:
          if str(i.particle) == 'neutron':
            p.append(i)

    return p

# For each source, calculate the area under the main peak to get a rough sense of
# how different the spectrum is. This is screened by looking at bins with more
# than 1e-3 probability AND above 10e6 MeV
integral_tokamak = 0
integral_slab = 0

# initialise a new source object - openmc plasma source (for the tokamak)
my_source = my_source
data = sample_initial_particles(my_source, prn_seed=1)
e_values = [particle.E for particle in data]

energy_filter = openmc.EnergyFilter.from_group_structure('CCFE-709')
bins = np.unique(energy_filter.bins.flatten())
probability, bin_edges = np.histogram(e_values, bins=bins, density=True)

energy = bin_edges[:-1]
unit_lethargy = energy_filter.lethargy_bin_width
e = np.unique(bin_edges)

p_times_e = [probability[i] * (e[i+1] - e[i]) for i in range(len(probability))]
plt.step(e[:-1], p_times_e, label='From Plasma', color='k')

for i in range(len(p_times_e)):
  if e[i] > 10e6 and p_times_e[i] > 1e-3:
    integral_tokamak += p_times_e[i]

import csv
with open('tokamak_source.csv', 'w') as csvfile:
  writer = csv.writer(csvfile)
  for i in range(len(p_times_e)):
    writer.writerow([e[i], p_times_e[i]])

# initialise a new source object - surface source file
my_source = openmc.FileSource('../surface_source.h5')

data = sample_initial_particles(my_source, prn_seed=1)
e_values = [particle.E for particle in data]

energy_filter = openmc.EnergyFilter.from_group_structure('CCFE-709')
bins = np.unique(energy_filter.bins.flatten())
probability, bin_edges = np.histogram(e_values, bins=bins, density=True)

energy = bin_edges[:-1]
unit_lethargy = energy_filter.lethargy_bin_width
e = np.unique(bin_edges)

p_times_e = [probability[i] * (e[i+1] - e[i]) for i in range(len(probability))]
plt.step(e[:-1], p_times_e, label='Entering Armor (WCLL)')

with open('slab_source.csv', 'w') as csvfile:
  writer = csv.writer(csvfile)
  for i in range(len(p_times_e)):
    writer.writerow([e[i], p_times_e[i]])

for i in range(len(p_times_e)):
  if e[i] > 10e6 and p_times_e[i] > 1e-3:
    integral_slab += p_times_e[i]

print('Tokamak, area under peak: ', integral_tokamak)
print('Slab, area under peak: ', integral_slab)

plt.xlim([1e-3, 1e8])
plt.ylim([1e-7, 1e1])
plt.xscale('log')
plt.xlabel('Neutron Energy (eV)')
plt.ylabel('Fraction of Neutrons Born in Bin')
plt.yscale('log')
plt.legend(loc='upper left', ncol=2)
plt.grid()
plt.savefig('source_energy.png', bbox_inches="tight")

plt.xlim([1e7, 2e7])
plt.savefig('source_energy_zoom.png', bbox_inches="tight")
plt.close()
