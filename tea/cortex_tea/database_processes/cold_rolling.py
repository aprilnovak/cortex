"""
Cold Rolling


########################################################################
DESCRIPTION
########################################################################

########################################################################
PARAMETERS
########################################################################
------------------------------------------------------------------------
name : str
    Name of the manufacturing process (e.g., 'CNC', 'HIP')

------------------------------------------------------------------------
tooling_level, equipment_level, time_level : str
    One of process_def.COST_LEVELS each (Low, Low-Medium, Medium,
    Medium-High, High, High-Very High, Very High). alphaT and beta are
    computed from these by process_def.build_process().

------------------------------------------------------------------------
description : str, optional
    A brief description of the process. This field is optional and
    provides additional context about the manufacturing method.

    Example:
        "Computer numerical control machining"
"""

name = 'Cold Rolling'

tooling_level = 'Medium'
equipment_level = 'High-Very High'
time_level = 'Low'

description = ""