"""
EROFER97 test


########################################################################
DESCRIPTION
########################################################################
Test variant of EROFER97 with every composition element tagged 'alloy'
(no 'residual' entries) - used to verify that cost calculation includes
every element when none of them are residual.

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

name = 'EROFER97 test'

density = 8

composition = {
    # alloying elements
    'Cr': {'wt': 9.0,   'type': 'alloy'},
    'W':  {'wt': 1.1,   'type': 'alloy'},
    'Mn': {'wt': 0.4,   'type': 'alloy'},
    'V':  {'wt': (0.25 + 0.15) / 2, 'type': 'alloy'},  # average = 0.20
    'Ta': {'wt': 0.12,  'type': 'alloy'},
    'B':  {'wt': 0.002, 'type': 'alloy'},


    # residuals / impurities - tagged 'alloy' here (test material)
    'C':  {'wt': 0.11,  'type': 'alloy'},
    'O':  {'wt': 0.01,  'type': 'alloy'},
    'N':  {'wt': 0.03,  'type': 'alloy'},
    'P':  {'wt': 0.005, 'type': 'alloy'},
    'S':  {'wt': 0.005, 'type': 'alloy'},
    'Nb': {'wt': 0.005, 'type': 'alloy'},
    'Mo': {'wt': 0.005, 'type': 'alloy'},
    'Ni': {'wt': 0.01,  'type': 'alloy'},
    'Cu': {'wt': 0.01,  'type': 'alloy'},
    'Al': {'wt': 0.01,  'type': 'alloy'},
    'Ti': {'wt': 0.02,  'type': 'alloy'},
    'Si': {'wt': 0.05,  'type': 'alloy'},
    'Co': {'wt': 0.01,  'type': 'alloy'},
    'As': {'wt': 0.05 / 4, 'type': 'alloy'},
    'Sn': {'wt': 0.05 / 4, 'type': 'alloy'},
    'Sb': {'wt': 0.05 / 4, 'type': 'alloy'},
    'Zr': {'wt': 0.05 / 4, 'type': 'alloy'}
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
