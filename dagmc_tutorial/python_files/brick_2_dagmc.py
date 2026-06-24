import openmc

model = openmc.Model()
# import geometry
dagmc_univ_bounded = openmc.DAGMCUniverse(filename='brick_2.h5m').bounded_universe()
model.geometry = openmc.Geometry(root=dagmc_univ_bounded)

# Attribute nuclides/elements to dagmc materials
# brick 1 (inside material) - water
water = openmc.Material(name="mat1")
water.add_nuclide('H1', 2.0)
water.add_nuclide('O16', 1.0)
water.set_density('g/cm3', 1.0)

# brick 2 (outside material) - ss316
ss316 = openmc.Material(name="mat2")
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

model.materials = openmc.Materials([water, ss316])

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

#Tallies
brick_water_filter = openmc.CellFilter(1)
tally_water = openmc.Tally(name="water_tally")
tally_water.filters = [brick_water_filter]
tally_water.scores = ['flux', 'total']  # or 'fission', 'absorption', etc.

brick_ss316_filter = openmc.CellFilter(2)
tally_ss316 = openmc.Tally(name="ss316_tally")
tally_ss316.filters = [brick_ss316_filter]
tally_ss316.scores = ['flux', 'total']  # or 'fission', 'absorption', etc.

tallies = openmc.Tallies([tally_water, tally_ss316])
model.tallies = tallies

model.export_to_xml()


