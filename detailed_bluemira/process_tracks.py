import openmc
import os

n_tracks = 10
for i in range(n_tracks):
  os.system("python neutronics_model.py")
  os.system("openmc neutronics_model.xml")

  tracks = openmc.Tracks('tracks.h5')
  tracks.write_to_vtk()

  os.system("mv tracks_0.vtp tokamak_tracks/tracks_" + str(i) + ".vtp")
