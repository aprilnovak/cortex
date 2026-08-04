# ----------------------------------------------------------------------
# Name the component

component_name = 'W Sheet'

# ----------------------------------------------------------------------
# Specify the material
#
# Materials can be selected from a directory (specified below)
# containing material input files, or a new material can be defined. To
# select a material from the directory, make sure to use the same name
# as specified in the material input file. To define a new material not
# found in the directory, add the name below and define the material
# parameters when prompted. You will be prompted to enter the material
# composition and the material-to-process compatibility relative cost
# coefficients for each process used in the component. Alternatively,
# materials can also be added directly to the directory. New materials 
# are not saved between runs. To access a material during a later run,
# you must manually upload the material to the directory.

# This component is a functionally graded W/Eurofer97 armor part, modeled
# as one effective blended material rather than a single pure material.
material_mode = 'blend'

# Two materials to blend (must be resolvable from the materials database
# and/or material_dir below)
material_a = 'W'
material_b = 'EROFER97'

# Fraction of material_a in the blend, and the basis it's given in
# ('volume' or 'mass')
blend_fraction_a = 0.75   # 75% tungsten, 25% Eurofer97, by volume
blend_basis = 'volume'

# A blended material has no material-to-process compatibility data of its
# own (it isn't safely inferable from the two source materials), so it
# must be supplied here. Reusing W's prior Spray Deposition Cmp (1.1) as
# an engineering placeholder for this tungsten-dominant blend.
blend_cmp = {'Spray Deposition': 1.1}

# Specify the directory containing material input files
material_dir = 'inputs_materials'

# ----------------------------------------------------------------------
# List the fabrication processes
#
# Processes can be selected from a list of pre-defined processes or a
# custom process can be specified. A custom process may be selected
# from a directory (specified below) containing process input files,
# or a new custom process can be defined. To select a custom process
# from the directory, make sure to use the same name as specified in
# the process input file. To define a new custom process not found in
# the directory, add the name below and define the process parameters
# when prompted. You will be prompted to enter the cost correlation
# coefficients for each process and the material-to-process
# compatibility relative cost coefficients for each material used in
# the component.

processes = ['Spray Deposition']

# Specify the directory containing process input files
process_dir = 'inputs_processes'

# Relative cost coefficients (Cc, Cs, Ct, Cf) and component volume are
# looked up from a named geometry in the geometries database instead of
# being hardcoded here. Cc/Cs/Ct/Cf default to 1.0 for any listed process
# not covered by the geometry's coefficient maps.

geometry = 'Plasma Facing Surface'

# Specify the directory containing custom geometry input files (or None
# to use only the geometries database)
geometry_dir = None

# List the scrap coefficients for each of the listed processes.
#
# Scrap coefficients determine how much more material, relative to part
# volume, will be rerquired to produce the component when acounting for
# the scrap fraction of each process.
Wc  =       [1.0]

# ----------------------------------------------------------------------
# Specify production quantity

production_qty = 1000

