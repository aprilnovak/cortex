"""
Oxide Dispersion Strengthened Steel


########################################################################
DESCRIPTION
########################################################################
Example. Cmp values have not be validated.

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

name = 'ODS'

density = 8

composition = {
    # TODO validate composition
    # alloying elements
    'Cr': {'wt': 13.5,   'type': 'alloy'}, 
    'Mn': {'wt': 0.0917, 'type': 'alloy'},
    'W':  {'wt': 0.99,   'type': 'alloy'},
    'Y':  {'wt': 0.146,  'type': 'alloy'}, 
    'Ti': {'wt': 0.16,   'type': 'alloy'},

    # residuals / impurities
    'Ni': {'wt': 0.0599, 'type': 'residual'},
    'H':  {'wt': 0.0057, 'type': 'residual'},
    'N':  {'wt': 0.0078, 'type': 'residual'},
    'O':  {'wt': 0.117,  'type': 'residual'},
    'C':  {'wt': 0.0164, 'type': 'residual'}
}

remainder_element = 'Fe'

import cost_variables as defaults

Cmp_map ={                    
    'CNC': 4,
    'Hot Rolling': 2,
    'Cold Rolling': 2,
    'HIP': 1.1,
    'Spray Deposition': 1.1,
    'Electron Beam': 1
}