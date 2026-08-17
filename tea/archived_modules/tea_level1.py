""" 
This file contains the TEA code.

"""
#import openmc
from cost_variables import material_cost_database

def TEA_Lv1(material):
    """
    Calculate part cost using the material cost method.

    Parameters
    ----------
    material : openmc.Material class
        Nuclide types can be accessed via material.nuclides[i].name
        Nuclide fraction in teh material can be accessed via material.nuclides[i].percent
        Nuclide percent type can be accessed via material.nuclides[i].percent_type

    Returns
    -------
    total_cost : float
        Part cost in $/kg.
    """
 # Build a dictionary of material masses in weight percent
    material_masses = {}
    for nuclide in material.nuclides:
        # Extract atomic symbol (letters only) from nuclide name
        symbol = ''.join(filter(str.isalpha, nuclide.name))
        if nuclide.percent_type == 'wo':
            mass = nuclide.percent/100
        else:
            raise ValueError(f"Unsupported percent_type '{nuclide.percent_type}' for nuclide {nuclide.name}. Only 'wo' is supported.")
        material_masses[symbol] = material_masses.get(symbol, 0.0) + mass

 # Calculate total cost
    total_cost = 0.0
    for symbol, mass in material_masses.items():
        if symbol in material_cost_database:
            unit_cost = material_cost_database[symbol].cost * mass
            total_cost += unit_cost
        else:
            print(f"Warning: {symbol} not found in database.")
    return total_cost


