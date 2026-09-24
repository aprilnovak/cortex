"""
Example Material


########################################################################
DESCRIPTION
########################################################################
Example input file.

########################################################################
PARAMETERS
########################################################################
------------------------------------------------------------------------
name : str
    Name of the material. Must be unique

------------------------------------------------------------------------
density : float 
    Density of the material in g/cm³

------------------------------------------------------------------------
composition : dict
    Dictionary defining the composition of the material. The keys are 
    element symbols (e.g., 'Fe', 'C', 'Cr'), and the values are 
    dictionaries with the following structure:

    composition[symbol] = {'wt': wt_percent, 'type': type}
    
    Where:
        symbol : str
            String representing an element (e.g., 'Fe', 'C', 'Cr').

        wt_percent : float
            Numeric value representing the proportion of the element in 
            weight percent (wt%).

        type : str
            String indicating the source of the element in the material:
                * 'alloy'    : Intentionally added via purchased 
                               material during alloying/processing

                * 'residual' : Present without requiring direct material 
                               purchase. This includes elements both 
                               intentionally and incidentally introduced 
                               during processing.

------------------------------------------------------------------------
remainder_element : str
    Element symbol representing the remainder element to balance 
    composition to 100 wt%. This element is not specified in the
    composition dictionary and is assumed to be the balance of the 
    material after accounting for the explicitly listed elements.

------------------------------------------------------------------------
Cmp_map : dict
    Dictionary mapping process names to material to process 
    compatability relative cost coefficients. The keys are process names
    (e.g., 'CNC', 'Hot Rolling'), and the values are the relative cost 
    of processing the material using the specified process compared to 
    an ideal material. 
    
    Cmp_map[process] = Cmp

    Where:
        Process: str
            String representing a manufacturing process (e.g., 'CNC', 
            'Hot Rolling'). All potetinal process options are listed 
            below, but only a subset will be used for each material
            based on user input.
                * CNC
                * Hot Rolling
                * Cold Rolling
                * HIP
                * Spray Deposition
                * Electron Beam

        Cmp: float
            Numeric (float) value representing the relative cost of 
            processing the material using the specified process compared
            to an ideal material. A value of 1 indicates the same cost
            as the ideal, while values greater than 1 indicate higher 
            costs.
"""

name = 'Example Material'

density = 6

composition = {
    # alloying elements
    'Cr': {'wt': 4.1, 'type': 'alloy'},
    'Ti': {'wt': 4.3, 'type': 'alloy'},

    # residuals / impurities
    'O':  {'wt': 0.035,  'type': 'residual'},
    'N':  {'wt': 0.022,  'type': 'residual'},
    'C':  {'wt': 0.020,  'type': 'residual'},
    'Si': {'wt': 0.087,  'type': 'residual'},
    'S':  {'wt': 0.002,  'type': 'residual'},
    'P':  {'wt': 0.004,  'type': 'residual'},
    'Nb': {'wt': 0.010,  'type': 'residual'},
    'Mo': {'wt': 0.010,  'type': 'residual'},
}

remainder_element = 'V'

import cost_variables as defaults

Cmp_map ={                    
    # Replace 'defaults.Cmp_map_1[process]' with the desired 'Cmp'
    'CNC': defaults.Cmp_map_1['CNC'],
    'Hot Rolling': defaults.Cmp_map_1['Hot Rolling'],
    'Cold Rolling': defaults.Cmp_map_1['Cold Rolling'],
    'HIP': defaults.Cmp_map_1['HIP'],
    'Spray Deposition': defaults.Cmp_map_1['Spray Deposition'],
    'Electron Beam': defaults.Cmp_map_1['Electron Beam']
}