import openmc

model = openmc.Model()
# import geometry
dagmc_univ = openmc.DAGMCUniverse(filename='brick_1.h5m')
model.geometry = openmc.Geometry(root=dagmc_univ)

# Attribute materials to regions
ss316 = openmc.Material(name="mat1")
ss316.add_element('C', 0.0003, 'wo')
ss316.add_element('Mn', 0.02, 'wo')
ss316.add_element('Si', 0.01, 'wo')
ss316.add_element('P', 0.0005, 'wo')
ss316.add_element('S', 0.0002, 'wo')
ss316.add_element('Cr', 0.185, 'wo')
ss316.add_element('Mo', 0.025, 'wo')
ss316.add_element('Ni', 0.13, 'wo')
ss316.add_element('Fe', 0.629, 'wo')
ss316.set_density("g/cm3", 8.0)

model.materials = openmc.Materials([ss316])

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
brick_cell_filter = openmc.CellFilter(1)
tally = openmc.Tally(name="brick_tally")
tally.filters = [brick_cell_filter]
tally.scores = ['flux', 'total']  # or 'fission', 'absorption', etc.
tallies = openmc.Tallies([tally])
model.tallies = tallies

model.export_to_xml()