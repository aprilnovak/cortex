#!/usr/bin/env python3

from openmc import Material, Materials
import openmc
import os
import subprocess

from make_materials import make_materials_xml
from make_source import make_plane_source

def createMaterials():
    return make_materials_xml()

def createPlots(mats_dict):
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

def createGeometry():
    dagmc_univ = openmc.DAGMCUniverse(filename="dagmc.h5m")
    geometry = openmc.Geometry(root=dagmc_univ)
    geometry.export_to_xml()

def createSettings( suppressOutput=True ):

    src = make_plane_source()

    # Create settings object
    settings = openmc.Settings()
    settings.source = src
    settings.run_mode = 'fixed source'
    settings.photon_transport = False

    # Turn on dagmc
    settings.dagmc = True
    settings.batches = 10
    settings.particles = 1000000

    ## Interpolate temperatures
    #settings.temperature['method'] = 'interpolation'

    settings.export_to_xml()
    print("Created settings.xml")

def createTallies():

    # Unstructured mesh to calculate tallies upon
    #meshname = "hcpb-mod.h5m"
    #umesh = openmc.UnstructuredMesh(meshname, library='moab')
    meshname = "HCLL-cm.e"
    umesh = openmc.UnstructuredMesh(meshname, library='libmesh')
    mesh_filter = openmc.MeshFilter(umesh)

    # Tallies
    tally = openmc.Tally()
    tally.filters = [mesh_filter]
    tally.scores = ['heating-local', 'flux']
    #tally.estimator = 'tracklength'
    tally.estimator = 'collision'
    tallies = openmc.Tallies([tally])
    tallies.export_to_xml()
    print("Created tallies.xml")
