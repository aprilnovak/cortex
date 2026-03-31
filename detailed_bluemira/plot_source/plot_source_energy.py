import numpy as np
import openmc
import openmc.lib
from tempfile import TemporaryDirectory
import matplotlib.pyplot as plt
import math

from tokamak_neutron_source import (
FluxMap,
FractionalFuelComposition,
TokamakNeutronSource,
TransportInformation,
)
from tokamak_neutron_source.profile import ParabolicPedestalProfile
from tokamak_neutron_source.reactions import Reactions

# I took all these hard-coded numbers from your OUT.DAT
# I haven't written this as a script, because this is an old PROCESS version, and
# it's a huge pain with their variable name changes...
temperature_profile = ParabolicPedestalProfile(2.37754249767383570e+01, 5.5, 0.1, 2.0, 1.45, 0.94) # [keV]
density_profile = ParabolicPedestalProfile(9.83393828196113777e+19, 5.86476334188244419e+19, 3.44986078934261350e+19, 1.0, 2.0, 0.94)
density_profile.set_scale(6.30197378312059699e+19/7.47634856031986811e+19)
rho_profile = np.linspace(0, 1, 30)

my_source = TokamakNeutronSource(
transport=TransportInformation.from_parameterisations(
ion_temperature_profile=temperature_profile,
fuel_density_profile=density_profile,
rho_profile=rho_profile,
fuel_composition=FractionalFuelComposition(D=0.5, T=0.5),
),
source_type=[Reactions.D_T, Reactions.D_D],
flux_map=FluxMap.from_eqdsk("../../equilibrium_eqdsk.json"),
cell_side_length=0.05,
)
my_source = my_source.to_openmc_source()

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
my_source = my_source

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

p_times_e = [probability[i] * (e[i+1] - e[i]) for i in range(len(probability))]

plt.step(e[:-1], p_times_e)
plt.xlim([1e-3, 1e8])
plt.xscale('log')
plt.xlabel('Neutron Energy (eV')
plt.ylabel('Fraction of Neutrons Born in Bin')
plt.yscale('log')
plt.grid()
plt.savefig('source_energy.png', bbox_inches="tight")
plt.close()

# write data to CSV so that we can plot it alongside the slab source
import csv
with open('tokamak_source.csv', 'w') as csvfile:
  writer = csv.writer(csvfile)
  for i in range(len(p_times_e)):
    writer.writerow([e[i], p_times_e[i]])
