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

material = 'W'

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

# List the relative cost coefficients for each of the listed processes.
#
# Relative cost coefficients determine how much more expensive it will 
# be to produce a component with more demanding features than the 'ideal 
# design'.

# Relative cost associated with producing components of different 
# geometrical complexity
Cc =        [1.0]

# Relative cost associated with size considerations and achieving 
# component section reductions/thickness
Cs =        [1.0]

# Relative cost associated with obtaining a specified tolerance
Ct =        [1.0]

# Relative cost associated with obtaining a specified a specified 
# surface finish
Cf =        [1.0]

# List the scrap coefficients for each of the listed processes.
#
# Scrap coefficients determine how much more material, relative to part 
# volume, will be rerquired to produce the component when acounting for 
# the scrap fraction of each process.
Wc  =       [1.0]

# ----------------------------------------------------------------------
# Specify component volume (mm3)

volume_mm3 = 3000000*.75  # 75% tungsten, 25% eurofer97     

# ----------------------------------------------------------------------
# Specify production quantity

production_qty = 1000 

