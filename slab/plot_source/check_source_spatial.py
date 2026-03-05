import numpy as np
import openmc
import openmc.lib
from tempfile import TemporaryDirectory
import math
import numpy as np
import matplotlib.pyplot as plt

# check that the source does not have any points being cut off outside
# the limits of the CSG planes

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
my_source = openmc.FileSource('../surface_source.h5')

ns = 1000000
data = sample_initial_particles(my_source, n_samples=ns, prn_seed=1)

# the four planes defining the chunk; all regions are in the negative halfspace
def plane1(x,y,z):
  return -0.110682 *x + -0.014572 *y + +0.993749 *z + 68.566475
def plane2(x,y,z):
  return -0.159268 *x + -0.020968 *y + -0.987013 *z + 123.101395
def plane3(x,y,z):
  return -0.130526 *x + +0.991445 *y + +0.000000 *z - 0.999806
def plane4(x,y,z):
  return -0.000000 *x + -1.000000 *y + -0.000000 *z + 0.000000
#Plane 1
#  Surface ID: 50000
#  Normal vector (a,b,c): (-0.110682, -0.014572, +0.993749)
#  Plane equation: -0.110682 x + -0.014572 y + +0.993749 z = -68.566475
#  Region operator: -p (inside)
#
#Plane 2
#  Surface ID: 50001
#  Normal vector (a,b,c): (-0.159268, -0.020968, -0.987013)
#  Plane equation: -0.159268 x + -0.020968 y + -0.987013 z = -123.101395
#  Region operator: -p (inside)
#
#Plane 3
#  Surface ID: 50002
#  Normal vector (a,b,c): (-0.130526, +0.991445, +0.000000)
#  Plane equation: -0.130526 x + +0.991445 y + +0.000000 z = +0.999806
#  Region operator: -p (inside)
#
#Plane 4
#  Surface ID: 50003
#  Normal vector (a,b,c): (-0.000000, -1.000000, -0.000000)
#  Plane equation: -0.000000 x + -1.000000 y + -0.000000 z = -0.000000
#  Region operator: -p (inside)

xo = []
yo = []
zo = []

for p in data:
  x = p.r[0]
  y = p.r[1]
  z = p.r[2]

  in_plane1 = plane1(x,y,z) < 0
  in_plane2 = plane2(x,y,z) < 0
  in_plane3 = plane3(x,y,z) < 0
  in_plane4 = plane4(x,y,z) < 0

  if (not in_plane1 or not in_plane2 or not in_plane3 or not in_plane4):
    xo.append(x)
    yo.append(y)
    zo.append(z)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.scatter(xo, yo, zo)
plt.savefig('outside.png')
plt.close()

print('Percent of points outside CSG box: ', len(xo) / len(data) * 100)
