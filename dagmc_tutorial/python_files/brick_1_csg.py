import openmc

model = openmc.Model()

# Attribute materials to regions
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

# Define surfaces for the cube: x=0 to 1, y=0 to 1, z=0 to 1
x0 = openmc.XPlane(x0=-0.5, boundary_type='vacuum')
x1 = openmc.XPlane(x0=0.5, boundary_type='vacuum')
y0 = openmc.YPlane(y0=-0.5, boundary_type='vacuum')
y1 = openmc.YPlane(y0=0.5, boundary_type='vacuum')
z0 = openmc.ZPlane(z0=-0.5, boundary_type='vacuum')
z1 = openmc.ZPlane(z0=0.5, boundary_type='vacuum')

# Create a region inside all six surfaces
region = +x0 & -x1 & +y0 & -y1 & +z0 & -z1

# define boundary conditions
cell = openmc.Cell(region =region, fill=mat1)
universe = openmc.Universe(cells=[cell])
geometry = openmc.Geometry(universe)
model.geometry= geometry

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

# provide tallies
brick_cell_filter = openmc.CellFilter(1)
tally = openmc.Tally(name="brick_general_tally")
tally.filters = [brick_cell_filter]
tally.scores = ['flux', 'total']  # or 'fission', 'absorption', etc.
tallies = openmc.Tallies([tally])
model.tallies = tallies

model.export_to_xml()
