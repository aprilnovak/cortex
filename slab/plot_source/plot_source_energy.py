import numpy as np
import openmc
import openmc.lib
from tempfile import TemporaryDirectory
import matplotlib.pyplot as plt
import math

from openmc_plasma_source import tokamak_source

def sample_initial_particles(this, n_samples: int = 1000, prn_seed: int = None):
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

    return particles

# initialise a new source object
my_source = openmc.FileSource('../surface_source.h5')

data = sample_initial_particles(my_source, n_samples=1000000, prn_seed=1)

e_values = [particle.E for particle in data]

# Calculate pdf for source energies
energy_filter = openmc.EnergyFilter.from_group_structure('CCFE-709')
bins = np.unique(energy_filter.bins.flatten())
probability, bin_edges = np.histogram(e_values, bins=bins, density=True)

# scaling by strength
energy = bin_edges[:-1]

unit_lethargy = energy_filter.lethargy_bin_width
e = np.unique(bin_edges)

# plot the tokamak source on top
import csv
import os
if (not os.path.exists('../../detailed_bluemira/plot_source/tokamak_source.csv')):
  raise ValueError("Go run the plot_source_energy.py script in ../../detailed_bluemira/plot_soure first!")

# load in the tokamak source
et = []
pt = []
with open('../../detailed_bluemira/plot_source/tokamak_source.csv') as csvfile:
  reader = csv.reader(csvfile)
  for row in reader:
    et.append(float(row[0]))
    pt.append(float(row[1]))

plt.step(et, pt, label='Neutrons From Plasma')

p_times_e = [probability[i] * (e[i+1] - e[i]) for i in range(len(probability))]

plt.step(e[:-1], p_times_e, label='Neutrons Entering Armor Front')

plt.legend()
plt.xlim([1e-3, 1e8])
plt.xscale('log')
plt.xlabel('Neutron Energy (eV')
plt.ylabel('Fraction of Neutrons Born in Bin')
plt.yscale('log')
plt.grid()
plt.savefig('source_energy.png', bbox_inches="tight")
plt.close()

# count fraction of neutrons under the 14.1 MeV peak (above 10 MeV)
frac_tokamak = 0
frac_slab = 0
for i in range(1, len(et)):
  if et[i] > 10e6:
    frac_tokamak += pt[i - 1]
for i in range(1, len(e)):
  if e[i] > 10e6:
    frac_slab += p_times_e[i - 1]
print('Slab fraction: ', frac_slab)
print('Tokamak fraction: ', frac_tokamak)
