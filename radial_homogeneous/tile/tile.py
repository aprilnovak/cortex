import openmc
import openmc.stats
import numpy as np
import matplotlib.pyplot as plt

import sys
import os
module_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..', 'materials'))
sys.path.append(module_path)
import materials

model = openmc.Model()

# OpenMC model of a bare tile with dimensions as given in the ARPA-E FOA.
# This is the simplest model possible, and will be both (i) a lower bound
# on our runtime and (ii) the baseline case we compare all materials against.

# TODO: numbers are about 5x lower than expected. Is my source rate wrong?

thickness = 5                      # [cm] thickness of the region
frontal_side = 50                  # [cm] side length of the tile facing the plasma
ncells = 10                        # number of cells in the radial direction
nwl = 1e6                          # [W/m2] neutron wall loading
e_per_neutron = 14.07e6            # [eV] energy carried by each neutron

dx = thickness / ncells
cell_volume = frontal_side**2 * dx
ev_to_joule = 1.60218e-19

neutron_source_rate = nwl / (100**2)   # W/cm2
neutron_source_rate /= ev_to_joule     # eV/s/cm2
neutron_source_rate /= e_per_neutron   # neutrons/s/cm2
neutron_source_rate *= frontal_side**2 # neutrons/s

# get materials; also fetch the number of each element's atom in the cell
# for later DPA normalization
t = materials.W(19.254)
t.volume = cell_volume
atoms_of_each_element = materials.atoms_of_each_element(t)
nuclides_of_each_element = materials.nuclides_for_each_element(t)

model.materials = openmc.Materials([t])

# define geometry; space is divided into cells for space-dependent tallies
# using tracklength estimators for better statistics
xplanes = [openmc.XPlane(x0=x) for x in np.linspace(0, thickness, ncells + 1)]
xcentroids = np.linspace(dx / 2, thickness - dx/2, ncells)

# TODO: figure out if I can have an actual surface source which lies on the
# boundary of a model; right now, I have to shift the model so that the plane
# source at x=0 is actually inside a cell
shift = 1e-6
xplanes[0].x0 = -shift

# TODO: need to think about the boundary condition on the incident face, it
# is not actually going to be vacuum
xplanes[0].boundary_type = 'vacuum'

xplanes[-1].boundary_type = 'vacuum'

ybot = openmc.YPlane(y0=0, boundary_type = 'reflective')
ytop = openmc.YPlane(y0=frontal_side, boundary_type = 'reflective')
zbot = openmc.ZPlane(z0=0, boundary_type = 'reflective')
ztop = openmc.ZPlane(z0=frontal_side, boundary_type = 'reflective')
chop = +ybot & -ytop & +zbot & -ztop

tile_cells = []
for i in range(ncells):
  tile_cells.append(openmc.Cell(region=+xplanes[i] & -xplanes[i+1] & chop, fill=t))

root_universe = openmc.Universe(cells=tile_cells)
model.geometry = openmc.Geometry(root=root_universe)

# define the neutron source
# TODO: use a Gaussian energy distribution instead of monoenergetic
energy_distribution = openmc.stats.Discrete([e_per_neutron], [1])

# TODO: use an isotropic angular distribution instead of all moving to the right
angle_distribution = openmc.stats.Monodirectional()

# TODO: determine if necessary to have a volume source or if this plane source
# is sufficient; if plane source is sufficient, is the correct angle distribution
# isotropic, or is it more like a white boundary?
x_distribution = openmc.stats.Discrete([0], [1])
y_distribution = openmc.stats.Uniform(a=0, b=frontal_side)
z_distribution = openmc.stats.Uniform(a=0, b=frontal_side)
space_distribution = openmc.stats.CartesianIndependent(x=x_distribution, y=y_distribution, z=z_distribution)

model.settings = openmc.Settings()
model.settings.source = openmc.IndependentSource(space=space_distribution, energy=energy_distribution, angle=angle_distribution)

# set other model settings
model.settings.particles = 500
model.settings.photon_transport = True
model.settings.batches = 100
model.settings.run_mode = 'fixed source'

# create general filters to be re-used across the tallies
cell_filter = openmc.CellFilter(tile_cells)
n_particle_filter = openmc.ParticleFilter(bins=['neutron'])
p_particle_filter = openmc.ParticleFilter(bins=['photon'])
energies = openmc.mgxs.GROUP_STRUCTURES['CCFE-709']
energy_filter = openmc.EnergyFilter(energies)
unit_lethargy = [np.log(energies[i+1]/energies[i]) for i in range(len(energies)-1)]

# add tallies for solution quantities
model.tallies = openmc.Tallies()

# flux score, with an energy filter
n_flux_tally = openmc.Tally()
n_flux_tally.filters = [cell_filter, n_particle_filter, energy_filter]
n_flux_tally.scores = ['flux']
model.tallies.append(n_flux_tally)

p_flux_tally = openmc.Tally()
p_flux_tally.filters = [cell_filter, p_particle_filter, energy_filter]
p_flux_tally.scores = ['flux']
model.tallies.append(p_flux_tally)

# for the dpa tally, we tally on a per-element basis, because each may have a different
# value of Ed. So, we need to create one dpa tally for each element in the material
dpa_tallies = []
for key in nuclides_of_each_element:
  dpa_tally = openmc.Tally()
  dpa_tally.scores = ['damage-energy']
  dpa_tally.filters = [cell_filter]
  dpa_tally.nuclides = nuclides_of_each_element[key]
  dpa_tallies.append(dpa_tally)
  model.tallies.append(dpa_tally)

statepoint = model.run()
with openmc.StatePoint(statepoint) as sp:
  n_tally = sp.get_tally(id=n_flux_tally.id)
  neutron_flux = n_tally.get_reshaped_data()

  for c in range(ncells):
    scaling = 1 / cell_volume * neutron_source_rate
    plt.loglog(energies[:-1], neutron_flux[c].flatten() * scaling / unit_lethargy, label='Cell {}'.format(c))

  plt.legend()
  plt.grid()
  plt.ylabel('Neutron Flux Per Unit Lethargy [1/cm$^2$/s]')
  plt.xlabel('Energy [eV]')
  plt.xlim([1e-2, 20e6])
  plt.savefig('n_flux_spectrum.png')
  plt.close()

  p_tally = sp.get_tally(id=p_flux_tally.id)
  photon_flux = p_tally.get_reshaped_data()

  for c in range(ncells):
    scaling = 1 / cell_volume * neutron_source_rate
    plt.loglog(energies[:-1], photon_flux[c].flatten() * scaling / unit_lethargy, label='Cell {}'.format(c))

  plt.legend()
  plt.grid()
  plt.ylabel('Photon Flux [1/cm$^2$/s/eV]')
  plt.xlabel('Energy [eV]')
  plt.xlim([1e-2, 20e6])
  plt.savefig('p_flux_spectrum.png')
  plt.close()

  # now, just plot the total fluxes by integrating over energy
  total_neutron_flux = np.zeros(ncells)
  for c in range(ncells):
    total_neutron_flux[c] += np.sum(neutron_flux[c].flatten()) * scaling

  plt.semilogy(xcentroids, total_neutron_flux, label='Neutron flux')

  total_photon_flux = np.zeros(ncells)
  for c in range(ncells):
    total_photon_flux[c] += np.sum(photon_flux[c].flatten()) * scaling

  plt.semilogy(xcentroids, total_photon_flux, label='Photon flux')
  plt.legend()
  plt.grid()
  plt.ylabel('Flux [1/cm$^2$/s]')
  plt.xlabel('Radial Position [cm]')
  plt.savefig('flux.png')
  plt.close()

  # create radial plots of the dpa; each of the tallies is for a particular element
  dpa_per_y = np.zeros(ncells)
  for tally in dpa_tallies:
    dpa_tally = sp.get_tally(id=tally.id)

    # TODO: damage-energy is only for neutrons, right? Believe so, the photon
    # bins return nothing
    dpa = dpa_tally.summation(nuclides=dpa_tally.nuclides).mean.flatten()

    # get the Ed for this element
    Ed = materials.Ed(materials.element(dpa_tally.nuclides))

    displacements_per_source = 0.8 * dpa / (2 * Ed)
    displacements_per_s = displacements_per_source * neutron_source_rate
    displacements_per_y = displacements_per_s * (365 * 24 * 60 * 60)
    displacements_per_y_per_all_atoms = displacements_per_y / materials.atoms(t)

    for i in range(ncells):
      dpa_per_y[i] += displacements_per_y_per_all_atoms[i]

  plt.semilogy(xcentroids, dpa_per_y)
  plt.grid()
  plt.ylabel('DPA/y')
  plt.xlabel('Radial Position [cm]')
  plt.savefig('dpa.png')
  plt.close()
