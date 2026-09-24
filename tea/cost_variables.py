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
    "Al" : MaterialCost(name="Aluminum", element="Al", cost=3.0),
    "As" : MaterialCost(name="Arsenic", element="As", cost=0.8),
    'B' : MaterialCost(name='Boron', element='B', cost=550), # check
    "C" : MaterialCost(name="Carbon", element="C", cost=10), # check
    "Co" : MaterialCost(name="Cobalt", element="Co", cost=56),
    'Cr' : MaterialCost(name="Chromium", element="Cr", cost=11.0),
    "Cu" : MaterialCost(name="Copper", element="Cu", cost=14),
    "Fe" : MaterialCost(name="Iron", element="Fe", cost=1.2),
    "Mg" : MaterialCost(name="Magnesium", element="Mg", cost=2.4),
    "Mn" : MaterialCost(name="Manganese", element="Mn", cost=1.9),
    "Mo" : MaterialCost(name="Molybdenum", element="Mo", cost=60),
    "N" : MaterialCost(name="SNitrogen", element="N", cost=0),
    "Nb" : MaterialCost(name="Niobium", element="Nb", cost=90),
    "Ni" : MaterialCost(name="Nickel", element="Ni", cost=21.0),
    "O" : MaterialCost(name="Oxygen", element="O", cost=0),
    "P": MaterialCost(name="Posphorous", element="P", cost=0),
    "Pb" : MaterialCost(name="Lead", element="Pb", cost=2.0),
    "S" : MaterialCost(name="Sulfur", element="S", cost=0),
    "Sb" : MaterialCost(name="Antimony", element="Sb", cost=18),
    "Si" : MaterialCost(name="Silicon", element="Si", cost=1.2),
    "Sn" : MaterialCost(name="Tin", element="Sn", cost=54),
    "Ta" : MaterialCost(name="Tantalum", element="Ta", cost=400.0),
    "Ti" : MaterialCost(name="Titanium", element="Ti", cost=15.0), # check
    "V" : MaterialCost(name="Vanadium", element="V", cost=21.0),
    "W" : MaterialCost(name="Tungsten", element="W", cost=50.0),
    "Y" : MaterialCost(name="Yttrium", element="Y", cost=33.0),
    "Zn" : MaterialCost(name="Zinc", element="Zn", cost=3.4),
    "Zr" : MaterialCost(name="Zirconium", element="Zr", cost=85)
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