import openmc

model = openmc.Model()
#import geometry
#dagmc_univ = openmc.DAGMCUniverse(filename='sphere_fine.h5m')
dagmc_univ = openmc.DAGMCUniverse(filename='sphere_coarse.h5m')
model.geometry = openmc.Geometry(root=dagmc_univ)

# Attribute materials to regions
# outside material (ss316)
mat1 = openmc.Material(name="mat1")
mat1.add_element('C', 0.0003, 'wo')
mat1.add_element('Mn', 0.02, 'wo')
mat1.add_element('Si', 0.01, 'wo')
mat1.add_element('P', 0.0005, 'wo')
mat1.add_element('S', 0.0002, 'wo')
mat1.add_element('Cr', 0.185, 'wo')
mat1.add_element('Mo', 0.025, 'wo')
mat1.add_element('Ni', 0.13, 'wo')
mat1.add_element('Fe', 0.629, 'wo')
mat1.set_density("g/cm3", 8.0)

model.materials = openmc.Materials([mat1])

# define simulation settings
model.settings.dagmc = True
model.settings.batches = 10
model.settings.particles = 10000
model.settings.run_mode = "fixed source"

# define the neutron source
source = openmc.Source()
source.space = openmc.stats.Point((0.0, 0.0, 0.0))  # (x, y, z) in cm
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([14.08e6], [1.0])
model.settings.source = source

# provide tallies
sphere_ss316_filter = openmc.CellFilter(1)
tally_ss316 = openmc.Tally(name="ss316_tally")
tally_ss316.filters = [sphere_ss316_filter]
tally_ss316.scores = ['flux', 'total']  # or 'fission', 'absorption', etc.

tallies = openmc.Tallies([tally_ss316])
model.tallies = tallies

model.export_to_xml()








