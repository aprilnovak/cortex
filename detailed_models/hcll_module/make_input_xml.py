#!/usr/bin/env python3

import openmc
import materials as fusion_mats
import matplotlib.pyplot as plt
import numpy as np

steel = fusion_mats.ss316(7.8)
steel.id=1
steel.name="SS"

#Lithium lead
lithium_lead = fusion_mats.PbLi(0.6, 1.593)
lithium_lead.id=2
lithium_lead.name="PbLi"

tungsten = fusion_mats.W(19.30)
tungsten.id=3
tungsten.name="W"

pyne_mats = [ steel, lithium_lead, tungsten ]
mats = openmc.Materials(pyne_mats)
mats.export_to_xml()
print("Created materials.xml")

colors = ["green", "blue", "orange"]
mats_dict = dict(zip(pyne_mats, colors))

p_xz = openmc.Plot()
p_xz.width = (1.5,0.85)
p_xz.pixels = (1000, 600)
p_xz.origin = (0.,0.412,0.)
p_xz.color_by = 'material'
p_xz.colors = mats_dict
p_xz.basis = 'xz'
plots = openmc.Plots([p_xz])
plots.export_to_xml()
print("Created plots.xml")

dagmc_univ = openmc.DAGMCUniverse(filename="dagmc.h5m")
geometry = openmc.Geometry(root=dagmc_univ)
geometry.export_to_xml()

# Specify spatial distribution for source location
# Specify each axis direction independently
x_disbn = openmc.stats.Uniform(-23.1,23.1)
y_disbn = openmc.stats.Uniform(-83.5,83.5)
z_disbn = openmc.stats.Discrete(x=[69.5],p=[1.0])
spatial_disbn = openmc.stats.CartesianIndependent(x_disbn,y_disbn,z_disbn)

# Specify an angular distribution for flux directionality
angle_disbn = openmc.stats.Monodirectional([0.,0.,-1.])

# Specify an energy distribution: 14 MeV neutrons
energy_disbn = openmc.stats.Discrete(x=[14.0e6],p=[1.0])

source = openmc.Source(space=spatial_disbn,
                       angle=angle_disbn,
                       energy=energy_disbn,
                       particle='neutron')

# Create settings object
settings = openmc.Settings()
settings.source = source
settings.run_mode = 'fixed source'
settings.photon_transport = False

settings.batches = 10
settings.particles = 10000

# add a volume calculation
nl = 70
dx = 69.5 / nl
vol_calcs = []
x = []
for l in range(nl):
  zmax = 69.5 - l*dx
  zmin = 69.5 - (l+1)*dx
  x.append(0.5*(zmin+zmax))
  print('Range in z: ', zmin, ' to ', zmax)
  vc = openmc.VolumeCalculation(pyne_mats, samples=1_000_000, lower_left=(-23.1, -83.5, zmin), upper_right=(23.1, 83.5, zmax))
  vol_calcs.append(vc)
settings.volume_calculations = vol_calcs
layer_volume = 23.1*2*83.5*2*69.5/nl

settings.export_to_xml()
print("Created settings.xml")

# Unstructured mesh to calculate tallies upon
meshname = "HCLL-cm.e"
umesh = openmc.UnstructuredMesh(meshname, library='libmesh')
mesh_filter = openmc.MeshFilter(umesh)

# Tallies
tally = openmc.Tally()
tally.filters = [mesh_filter]
tally.scores = ['heating-local', 'flux']
tally.estimator = 'collision'
tallies = openmc.Tallies([tally])
tallies.export_to_xml()
print("Created tallies.xml")

openmc.calculate_volumes()
steel = np.zeros(len(vol_calcs))
breeder = np.zeros(len(vol_calcs))
tungsten = np.zeros(len(vol_calcs))
helium = np.zeros(len(vol_calcs))

idx = 1
for ii, v in enumerate(vol_calcs):
  i = len(vol_calcs) - 1 - ii
  v.load_results('volume_' + str(idx) + '.h5')
  steel[i] = (v.volumes[1].n+v.volumes[3].n)/layer_volume
  breeder[i] = v.volumes[2].n/layer_volume
  #tungsten[i] = v.volumes[3].n/layer_volume
  leftover = (layer_volume - v.volumes[1].n - v.volumes[2].n - v.volumes[3].n)/layer_volume
  helium[i] = max(leftover, 0.0)
  idx += 1

layers = [3, 9, 10, 10, 3, 9, 3, 10, 10, 3]
fig, ax = plt.subplots()

left = 0
for l in layers:
  ax.vlines(left+l, color='k', ymin=1.0, ymax=1.1)
  left += l

print('Total volume: ', 23.1*2*83.5*2*69.5)
print('Volume of each layer: ', layer_volume)
print('dx: ', dx)
plt.bar(x, steel, label='Structural', color='magenta', width=1)
plt.bar(x, breeder, bottom=steel, label='PbLi', color='cyan', width=1)
plt.bar(x, helium, bottom=breeder+steel, label='Helium', color='blue', width=1)
#plt.bar(x, tungsten, bottom=helium+breeder+steel, label='First Wall', color='magenta')
plt.legend(loc='center left')
plt.xlabel('Distance from Plasma')
plt.xlim([0, 70])
plt.ylim([0.0, 1.1])
plt.ylabel('Fraction of Volume')
plt.savefig('hcll_layers.png', bbox_inches="tight")
plt.close()


tungsten = list(reversed(tungsten))
steel = list(reversed(steel))
helium = list(reversed(helium))
breeder = list(reversed(breeder))

start = 0
for j in layers:
  #print(tungsten[start:start+j])
  print('\nLayer of thickness ', j)
  print('Structural: ', np.sum(steel[start:start+j])/j)
  print('Coolant: ', np.sum(helium[start:start+j])/j)
  print('Breeder: ', np.sum(breeder[start:start+j])/j)
  print('Sum: ', (np.sum(steel[start:start+j]))/j+np.sum(helium[start:start+j])/j+np.sum(breeder[start:start+j])/j)
  start = start+j


