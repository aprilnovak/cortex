"""
Ditungsten Carbide (W2C, Ni binder)


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

name = 'W2C'

density = 16.767

composition = {
    # alloying elements
    'C':  {'wt': 4.0,  'type': 'alloy'},     # 3.8 - 4.2
    'Ni': {'wt': 4.0,  'type': 'alloy'},     # 3.0 - 5.0

    # residuals / impurities (max limits)
    'Cr': {'wt': 0.1,  'type': 'residual'},
    'V':  {'wt': 0.1,  'type': 'residual'},
    'Fe': {'wt': 0.05, 'type': 'residual'},
    'Ti': {'wt': 0.1,  'type': 'residual'},
    'Co': {'wt': 0.1,  'type': 'residual'},
    # W balance: 91.55
}

remainder_element = 'W'

from cortex_tea import cost_variables as defaults

Cmp_map ={                    
    # Replace 'defaults.Cmp_map_1[process]' with the desired 'Cmp'
    # TODO: Validate values
    'HIP': 2,
    'Spray Deposition': 1.5
}