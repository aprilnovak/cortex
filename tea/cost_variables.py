""" 
This database contains the raw material cost data for materials common in fusion reactors.
"""

class MaterialCost:
    def __init__(self, name, element, cost):
        """
        Initialize a material with basic techno-economic attributes.

        Parameters:
        - name (str): Common name of the material
        - element (str): Primary element or alloy base
        - cost (float): Cost per unit (e.g., $/kg)
        """
        self.name = name
        self.element = element
        self.cost = cost

    def __str__(self):
        return f"{self.name} ({self.element}): ${self.cost}/kg"

material_cost_database = {
    'B' : MaterialCost(name='Boron', element='B', cost=550),
    'Cr' : MaterialCost(name="Chromium", element="Cr", cost=12.0),
    "Fe" : MaterialCost(name="Iron", element="Fe", cost=1.2),
    "Mn" : MaterialCost(name="Manganese", element="Mn", cost=2.2),
    "Ni" : MaterialCost(name="Nickel", element="Ni", cost=22.0),
    "Ta" : MaterialCost(name="Tantalum", element="Ta", cost=350.0),
    "Ti" : MaterialCost(name="Titanium", element="Ti", cost=15.0),
    "V" : MaterialCost(name="Vanadium", element="V", cost=20.0),
    "W" : MaterialCost(name="Tungsten", element="W", cost=40.0),
    "Y" : MaterialCost(name="Yttrium", element="Y", cost=33.0)    
}

# The map below was based on generic stainless steel
Cmp_map_1 ={                    
    'CNC': 4,
    'Hot Rolling': 2,
    'Cold Rolling': 2,
    'HIP': 1.1,
    'Spray Deposition': 1.1,
    'Electron Beam': 1
}