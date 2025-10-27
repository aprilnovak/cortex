import openmc

model = openmc.Model()

ss316 = openmc.Material(name="ss316")
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

# implement neutron source
# define boundary conditions
sphere = openmc.Sphere(r=0.5, boundary_type='vacuum')
cell = openmc.Cell(region =-sphere, fill=ss316)
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

sphere_ss316_filter = openmc.CellFilter(cell)
tally_ss316 = openmc.Tally(name="ss316_tally")
tally_ss316.filters = [sphere_ss316_filter]
tally_ss316.scores = ['flux', 'total']  # or 'fission', 'absorption', etc.

tallies = openmc.Tallies([tally_ss316])
model.tallies = tallies

model.export_to_xml()







