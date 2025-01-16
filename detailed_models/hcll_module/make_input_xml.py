#!/usr/bin/env python3

from make_settings import *

mats_dict = createMaterials()
createGeometry()
createSettings()
createTallies()

# Make a plot...
createPlots(mats_dict)
