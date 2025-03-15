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
model.settings.particles = 1000
model.settings.photon_transport = True
model.settings.batches = 100
model.settings.run_mode = 'fixed source'

# create general filters to be re-used across the tallies
cell_filter = openmc.CellFilter(tile_cells)
particle_filter = openmc.ParticleFilter(bins=['neutron', 'photon'])

# add tallies for solution quantities
model.tallies = openmc.Tallies()

cell_tally = openmc.Tally()
cell_tally.filters = [cell_filter, particle_filter]
cell_tally.scores = ['flux', 'heating']
model.tallies.append(cell_tally)

# for the dpa tally, we tally on a per-element basis, because each may have a different
# value of Ed. So, we need to create one dpa tally for each element in the material
dpa_tallies = []
for key in nuclides_of_each_element:
  dpa_tally = openmc.Tally()
  dpa_tally.scores = ['damage-energy']
  dpa_tally.filter = cell_filter
  dpa_tally.nuclides = nuclides_of_each_element[key]
  dpa_tallies.append(dpa_tally)
  model.tallies.append(dpa_tally)

statepoint = model.run()
with openmc.StatePoint(statepoint) as sp:
  tally = sp.get_tally(id=cell_tally.id)
  neutron_flux = tally.get_slice(scores=['flux'], filters=[type(particle_filter)], filter_bins=[('neutron',)])
  photon_flux = tally.get_slice(scores=['flux'], filters=[type(particle_filter)], filter_bins=[('photon',)])
  neutron_flux = neutron_flux.mean.flatten() / cell_volume * neutron_source_rate
  photon_flux = photon_flux.mean.flatten() / cell_volume * neutron_source_rate

  # create radial plots of the flux
  plt.semilogy(xcentroids, neutron_flux, label='Neutron flux')
  plt.semilogy(xcentroids, photon_flux, label='Photon flux')
  plt.legend()
  plt.grid()
  plt.ylabel('Flux [1/cm$^2$/s]')
  plt.xlabel('Radial Position [cm]')
  plt.savefig('scores.png')
  plt.close()

  # create radial plots of the dpa
  # TODO: damage-energy is only for neutrons, right?
  dpa = tally.get_slice(scores=['damage-energy'], filters=[type(particle_filter)], filter_bins=[('neutron',)])

  # TODO: get actual Ed for each nuclide, will need to generalize this for alloys.
  Ed = 90
  displacements_per_source = 0.8 * dpa.mean.flatten() / (2 * Ed)
  displacements_per_s = displacements_per_source * neutron_source_rate
  displacements_per_y = displacements_per_s * (365 * 24 * 60 * 60)
  displacements_per_atom_per_y = displacements_per_y / atoms_of_element['W']
  print(displacements_per_atom_per_y)

  plt.semilogy(xcentroids, displacements_per_atom_per_y)
  plt.grid()
  plt.ylabel('DPA/y')
  plt.xlabel('Radial Position [cm]')
  plt.savefig('dpa.png')
  plt.close()
