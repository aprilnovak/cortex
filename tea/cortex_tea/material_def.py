"""
Material Definition


########################################################################
DESCRIPTION
########################################################################
Defines the 'Material' class and functions to build materials from input
data.

########################################################################
CLASSES
########################################################################
Material(name: str)
    Represents a material.

########################################################################
FUNCTIONS
########################################################################
------------------------------------------------------------------------
build_material(name, density, composition, remainder_element, Cmp_map)
------------------------------------------------------------------------
load_materials(material_dir)
"""

# ######################################################################
# CLASSES
# ######################################################################
# ----------------------------------------------------------------------
# Material class

class Material:
    """
    Represents a material.

    
    ########################################################################
    ATTRIBUTES
    ########################################################################
    ------------------------------------------------------------------------
    name : str
        Name of the material

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

                    * 'rem'      : Remainder element to balance composition
                                to 100 wt%
    ------------------------------------------------------------------------
    density : float
        Density of the material (g/cm³)
    
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
    
    ------------------------------------------------------------------------
    cost : float
        Material cost ($/kg) calculated from composition

    ########################################################################
    METHODS
    ########################################################################
    ------------------------------------------------------------------------
    add_element(symbol: str, wt_percent: float, type: str)
        Add or set an element wt fraction.  

    ------------------------------------------------------------------------ 
    add_Cmp(process: str, Cmp: float)
        Add or set a relative cost coefficient for a given process.

    ------------------------------------------------------------------------
    set_density(density: float)
        Set the density of the material.

    ------------------------------------------------------------------------
    get_total_wt()
        Calculate total mass currently specified.

    ------------------------------------------------------------------------
    calculate_cost_breakdown(cost_db)
        Per-element $/kg contribution to material cost.

    ------------------------------------------------------------------------
    calculate_cost(cost_db)
        Calculate material cost ($/kg) from composition.

    ------------------------------------------------------------------------
    get_Cmp(process: str)
        Look up relative cost coefficient for a given material and process.
    """

    def __init__(self, name: str):
        self.name = name
        self.composition = {}   # {wt: float, type: str}
        self.density = None     # g/cm^3
        self.Cmp_map = {}
        self.cost = None        # $/kg

    def add_element(self, symbol: str, wt_percent: float, type: str):
        """
        Add or set an element wt fraction.


        ########################################################################
        DESCRIPTION
        ########################################################################
        Adds and element to the 'composition' attribute of a 'Material' object

        ########################################################################
        PARAMETERS
        ########################################################################
        ------------------------------------------------------------------------
        symbol : str
            String representing an element (e.g., 'Fe', 'C', 'Cr').

        ------------------------------------------------------------------------
        wt_percent : float
            Numeric value representing the proportion of the element in 
            weight percent (wt%).

        ------------------------------------------------------------------------
        type : str
            String indicating the source of the element in the material:
                * 'alloy'    : Intentionally added via purchased 
                               material during alloying/processing

                * 'residual' : Present without requiring direct material 
                               purchase. This includes elements both 
                               intentionally and incidentally introduced 
                               during processing.

                * 'rem'      : Remainder element to balance composition
                               to 100 wt%
        """

        # If element already exists, warn user of override
        # TODO: Ask user if they want to proceed
        if symbol in self.composition:
            print(f"Warning: element '{symbol}' already exists. Overriding " f"{self.composition[symbol]} wt with new value of {wt_percent}")

        # Override / add element in composition
        self.composition[symbol] = {'wt': wt_percent, 'type': type}

    def add_Cmp(self, process: str, Cmp: float):
        """
        Add or set a relative cost coefficient for a given process.

        
        ########################################################################
        PARAMETERS
        ########################################################################
        ------------------------------------------------------------------------
        process : str
            Name of a process (e.g., CNC)

        ------------------------------------------------------------------------
        Cmp : float
            Relative cost associated with material to process compatability
        """

        # If process already exists, warn user of override
        # TODO: Ask user if they want to proceed
        if process in self.Cmp_map:
            print(f"Warning: process '{process}' already exists. Overriding " f"{self.Cmp_map[process]} Cmp with new value of {Cmp}")

        # Override process in the map
        self.Cmp_map[process] = {'Cmp': Cmp }

    def set_density(self, density: float):
        """
        Set the density of the material.
        """
        self.density = density

    def get_total_wt(self):
        """
        Calculate total mass currently specified.
        """
        return sum(el_data['wt'] for el_data in self.composition.values())
  
    def calculate_cost_breakdown(self, cost_db):
        """
        Per-element $/kg contribution to material cost. Only 'alloy' and
        'rem' elements contribute; 'residual' elements are free (present
        without requiring direct material purchase) and always report 0.0.


        ########################################################################
        PARAMETERS
        ########################################################################
        ------------------------------------------------------------------------
        cost_db : dict
            Dictionary of elements and associated costs

        ########################################################################
        RETURNS
        ########################################################################
        ------------------------------------------------------------------------
        rows : list of dict
            One entry per element: {'element', 'wt', 'type', 'dollar_per_kg'}
        missing_elements : list of str
            Cost-contributing elements absent from cost_db
        """

        rows = []
        missing_elements = []

        for element, data in self.composition.items():
            wt_percent = data['wt']
            el_type = data['type']
            contributes = el_type in ('alloy', 'rem')

            dollar_per_kg = 0.0
            if contributes:
                if element not in cost_db:
                    missing_elements.append(element)
                else:
                    dollar_per_kg = (wt_percent / 100.0) * cost_db[element].cost

            rows.append({'element': element, 'wt': wt_percent, 'type': el_type, 'dollar_per_kg': dollar_per_kg})

        return rows, missing_elements

    def calculate_cost(self, cost_db):
        """
        Calculate material cost ($/kg) from composition.


        ########################################################################
        PARAMETERS
        ########################################################################
        ------------------------------------------------------------------------
        cost_db : dict
            Dictionary of elements and associated costs
        """

        rows, missing_elements = self.calculate_cost_breakdown(cost_db)
        total_cost = sum(row['dollar_per_kg'] for row in rows)

        if missing_elements:
            warn_block(f"Warning: Missing cost data for elements {missing_elements} in {self.name}")

        return total_cost

    def get_Cmp(self, process: str) -> float:
        """
        Look up relative cost coefficient for a given process.

        
        ########################################################################
        PARAMETERS
        ########################################################################
        ------------------------------------------------------------------------
        process : str
            Name of a process (e.g., 'CNC')

        ########################################################################
        RETURNS
        ########################################################################
        ------------------------------------------------------------------------
        Cmp : float
            Relative cost coefficient for the material and process
        """

        if process not in self.Cmp_map:
            raise KeyError(
                f"Process '{process}' not found in Cmp_map for material '{self.name}'."
            )

        Cmp = self.Cmp_map[process]['Cmp']

        return Cmp

# ######################################################################
# FUNCTIONS
# ######################################################################
# ----------------------------------------------------------------------
# Function to build a 'Material' object from input data.

from cost_variables import material_cost_database as cost_db

def build_material(name, density, composition, remainder_element, Cmp_map):
    """
    Build a 'Material' object from input data.


    ########################################################################
    PARAMETERS
    ########################################################################
    ------------------------------------------------------------------------
    name : str
        Name of the material (e.g., 'V4Cr4Ti')

    ------------------------------------------------------------------------
    density : float
        Material density (g/cm³)

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
    # ------------------------
    # Create 'Material' object
    mat = Material(name)

        # ----------------------
    # Verify that input data
    if density <= 0:
        print_error(f"Error loading material, '{name}'. Material skipped. Density must be greater than zero.")
        return None

       
    for element in composition:

        if element.strip() == "":
            print_error(f"Error loading material, '{name}'. Material skipped. Element symbol cannot be empty.")
            return None

    total_wt = 0

    for element, data in composition.items():

        # Validate weight %
        try:
            data["wt"] = positive_float(data["wt"])
        except Exception:
            print_error(f"Error loading material, '{name}'. Material skipped. Weight % for element, '{element}', must be a positive number.")
            return None

        # Validate type
        if data.get("type") not in ["alloy", "residual"]:
            print_error(f"Error loading material, '{name}'. Material skipped. Type for element, '{element}', must be 'alloy' or 'residual'.")
            return None

        total_wt += data["wt"]

    if total_wt > 100:
        print_error(f"Error loading material, '{name}'. Material skipped. Sum of alloying and residual weight cannot exceed 100%.")
        return None

    if remainder_element in composition:
        print_error(f"Error loading material, '{name}'. Material skipped. Remainder element, '{remainder_element}', is already defined in composition.")
        return None

    if remainder_element.strip() == "":
        print_error(f"Error loading material, '{name}'. Material skipped. Element symbol cannot be empty.")
        return None

    for pname in Cmp_map:

        try:
            Cmp_map[pname] = positive_float(Cmp_map[pname])
        except Exception:
            print_error(f"Error loading material, '{name}'. Material skipped. Cmp for process, '{pname}', must be a positive number.")
            return None

    # ------------
    # Add elements
    for el, el_data in composition.items():
        mat.add_element(el, el_data['wt'], el_data['type'])

    # Add remainder element
    remainder = 100 - total_wt
    mat.add_element(remainder_element, remainder, type="rem")

    # -----------
    # Set density
    mat.set_density(density)

    # -----------------------
    # Calculate material cost
    mat.cost = mat.calculate_cost(cost_db)

    # ----------------------------------------
    # Define material to process compatability
    for process, val in Cmp_map.items():
        mat.add_Cmp(process, val)

    # ----------------------------
    # Display material information
    print("\n➤ ", name)
    print(f"       Alloying elements (wt%): {total_wt:.2f}%")
    print(f"       Remainder ({remainder_element}) (wt%):   {' ' if len(remainder_element) == 1 else ''} {remainder}%")
    print(f"       Material cost:           ${mat.cost:.2f}/kg")

    return mat

# ----------------------------------------------------------------------
# Function to blend two materials into one effective material.

def blend_materials(name, material_a, material_b, fraction_a, basis="volume"):
    """
    Blend two materials into a single effective 'Material' object.

    Models a part made from two materials (e.g. a functionally graded W/
    Eurofer97 armor part) as one homogeneous material with rule-of-mixtures
    density and mass-weighted composition/cost. This is NOT a spatially
    resolved gradient - 'Component' only supports a single bulk material.


    ########################################################################
    PARAMETERS
    ########################################################################
    ------------------------------------------------------------------------
    name : str
        Name of the blended material

    ------------------------------------------------------------------------
    material_a, material_b : Material
        The two materials to blend

    ------------------------------------------------------------------------
    fraction_a : float
        Fraction of material_a in the blend (0 < fraction_a < 1), in the
        units given by 'basis'. material_b makes up the remainder.

    ------------------------------------------------------------------------
    basis : str
        'volume' (default): fraction_a is a volume fraction; density is
            blended by the volume-weighted rule of mixtures.
        'mass': fraction_a is a mass fraction; density is blended by the
            mass-weighted (harmonic) rule of mixtures.

    ########################################################################
    RETURNS
    ########################################################################
    ------------------------------------------------------------------------
    Material or None
        The blended material, with an empty Cmp_map (process compatibility
        for a novel blend isn't safely inferable from the source materials
        and must be supplied by the caller, same as a custom material).
        Returns None (and prints an error) if inputs are invalid.
    """

    if basis not in ("volume", "mass"):
        print_error(f"Error blending material, '{name}'. Blend skipped. Basis must be 'volume' or 'mass'.")
        return None

    try:
        fraction_a = positive_float(fraction_a)
    except Exception:
        print_error(f"Error blending material, '{name}'. Blend skipped. Fraction of material_a must be a positive number.")
        return None

    if fraction_a >= 1:
        print_error(f"Error blending material, '{name}'. Blend skipped. Fraction of material_a must be less than 1.")
        return None

    # ------------------------------------
    # Blended density and mass fractions
    if basis == "volume":
        mass_a = fraction_a * material_a.density
        mass_b = (1 - fraction_a) * material_b.density
        density = mass_a + mass_b
        w_a = mass_a / density
        w_b = mass_b / density
    else:
        w_a = fraction_a
        w_b = 1 - fraction_a
        density = 1 / (w_a / material_a.density + w_b / material_b.density)

    # ------------------------------------
    # Blended composition (mass-weighted average of each element's wt%;
    # a source material's remainder element becomes a concrete alloy
    # component in the blend rather than a remainder placeholder)
    elements = set(material_a.composition) | set(material_b.composition)
    composition = {}

    for element in elements:
        data_a = material_a.composition.get(element)
        data_b = material_b.composition.get(element)
        wt_a = data_a["wt"] if data_a else 0.0
        wt_b = data_b["wt"] if data_b else 0.0

        is_alloy = (data_a and data_a["type"] in ("alloy", "rem")) or (data_b and data_b["type"] in ("alloy", "rem"))
        el_type = "alloy" if is_alloy else "residual"

        composition[element] = {"wt": w_a * wt_a + w_b * wt_b, "type": el_type}

    # ------------------------------------
    # Build blended material
    mat = Material(name)
    mat.set_density(density)

    for element, data in composition.items():
        mat.add_element(element, data["wt"], data["type"])

    mat.cost = mat.calculate_cost(cost_db)

    print("\n➤ ", name)
    print(f"       Blend:                    {fraction_a*100:.1f}% {material_a.name} / {(1-fraction_a)*100:.1f}% {material_b.name} (by {basis})")
    print(f"       Density:                  {density:.3f} g/cm^3")
    print(f"       Material cost:            ${mat.cost:.2f}/kg")

    return mat

# ----------------------------------------------------------------------
# Function to load custom materials from materials directory.

import os
import importlib.util

def load_materials(material_dir):
    """
    Load custom materials from materials directory.


    ########################################################################
    PARAMETERS
    ########################################################################
    ------------------------------------------------------------------------
    material_dir : str
        Name of folder containing material input files

    ########################################################################
    RETURNS
    ########################################################################
    ------------------------------------------------------------------------
    materials : 
        List of 'Material' objcets from materials directory
    """
    materials = []
    
    if not os.path.isdir(material_dir):
        print(f"\nDirectory '{material_dir}' not found. Skipping material loading.")
        return materials

    for filename in os.listdir(material_dir):

        if not filename.endswith(".py"):
            continue

        filepath = os.path.join(material_dir, filename)
        module_name = os.path.splitext(filename)[0]

        # Load module
        spec = importlib.util.spec_from_file_location(module_name, filepath)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        # Extract data
        name = module.name
        density = module.density
        composition = module.composition
        remainder_element = module.remainder_element
        Cmp_map = module.Cmp_map

        # Build material
        mat = build_material(name, density, composition, remainder_element, Cmp_map)

        materials.append(mat)

    return materials

# ----------------------------------------------------------------------
# Functions to format error checking
import textwrap
def print_error(msg, width=72):

    print("\n        ---------------------- ✖ ERROR ✖ -----------------------        ")

    wrapped = textwrap.wrap(msg, width - 20)

    for line in wrapped:
        print(f"        |{line.center(width - 18)}|")

    # blank spacer line
    print(f"        |{' ' * (width -18)}|")
    print("        " + "-" * (width - 16) + "")

def positive_float(value):
    v = float(value)
    if v <= 0:
        raise ValueError("Please enter a positive numeric value.")
    return v

def warn_block(messages, width=72):
    print("\n        --------------------- ⚠ WARNING ⚠ ----------------------        ")

    # Allow single string or list
    if isinstance(messages, str):
        messages = [messages]

    for msg in messages:
        wrapped = textwrap.wrap(msg, width - 20)

        for line in wrapped:
            print(f"        |{line.center(width - 18)}|       ")

        print("        |                                                      |        ")

    print("        --------------------------------------------------------        ")