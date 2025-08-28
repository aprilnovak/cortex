""" 
This file contains the TEA code.

"""
from cost_variables import material_cost_database

def TEA_Lv1(material_masses):
    """
    Calculate part cost using the material cost method.

    Parameters
    ----------
    material_masses : dict
        Keys are element symbols (e.g., "Fe"), values are mass fractions.

    Returns
    -------
    total_cost : float
        Part cost in $/kg.
    """


    total_cost = 0.0
    for symbol, mass in material_masses.items():
        if symbol in material_cost_database:
            unit_cost = material_cost_database[symbol].cost * mass
            total_cost += unit_cost
        else:
            print(f"Warning: {symbol} not found in database.")
    return total_cost


