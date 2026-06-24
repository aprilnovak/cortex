import openmc

model = openmc.Model()
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

# Surfaces
# Inner cube (1x1x1)
x1_min = openmc.XPlane(x0=-0.5)
x1_max = openmc.XPlane(x0=0.5)
y1_min = openmc.YPlane(y0=-0.5)
y1_max = openmc.YPlane(y0=0.5)
z1_min = openmc.ZPlane(z0=-0.5)
z1_max = openmc.ZPlane(z0=0.5)

# Outer cube (2x2x2)
x2_min = openmc.XPlane(x0=-1.0, boundary_type='vacuum')
x2_max = openmc.XPlane(x0=1.0, boundary_type='vacuum')
y2_min = openmc.YPlane(y0=-1.0, boundary_type='vacuum')
y2_max = openmc.YPlane(y0=1.0, boundary_type='vacuum')
z2_min = openmc.ZPlane(z0=-1.0, boundary_type='vacuum')
z2_max = openmc.ZPlane(z0=1.0, boundary_type='vacuum')

# Regions
# Inner cube: material 1
inner_region = +x1_min & -x1_max & +y1_min & -y1_max & +z1_min & -z1_max
# Outer hollow shell: material 2
outer_shell_region = (+x2_min & -x2_max & +y2_min & -y2_max & +z2_min & -z2_max) & ~inner_region
# Cells
inner_cell = openmc.Cell(name='Inner Cube', fill=water, region=inner_region)
outer_cell = openmc.Cell(name='Outer Shell', fill=ss316, region=outer_shell_region)
# Define the universe
root_universe = openmc.Universe(cells=[inner_cell, outer_cell])
geometry = openmc.Geometry(root_universe)
model.geometry = geometry

# define simulation settings
model.settings.batches = 10
model.settings.particles = 10000
model.settings.run_mode = "fixed source"

# define the neutron source
source = openmc.Source()
source.space = openmc.stats.Point((0.0, 0.0, 0.0))  # (x, y, z) in cm
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.Discrete([14.08e6], [1.0])
model.settings.source = source

# Tallies
brick_water_filter = openmc.CellFilter(inner_cell)
tally_water = openmc.Tally(name="water_tally")
tally_water.filters = [brick_water_filter]
tally_water.scores = ['flux', 'total']  # or 'fission', 'absorption', etc.

brick_ss316_filter = openmc.CellFilter(outer_cell)
tally_ss316 = openmc.Tally(name="ss316_tally")
tally_ss316.filters = [brick_ss316_filter]
tally_ss316.scores = ['flux', 'total']  # or 'fission', 'absorption', etc.

tallies = openmc.Tallies([tally_water, tally_ss316])
model.tallies = tallies

model.export_to_xml()
