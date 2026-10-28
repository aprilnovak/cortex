"""
Component Geometry Definition


########################################################################
DESCRIPTION
########################################################################
Defines the 'Geometry' class and functions to build a component's
geometry and look up the design-dependent relative cost coefficients
(Cc, Cs, Ct, Cf) for each fabrication process.

A geometry captures:
    * volume_mm3            - component volume (mm^3)
    * section_thickness_mm  - representative section/wall thickness (mm)
    * shape_class           - shape classification band (e.g. 'C1'); see
                               Figs. 3.8/3.9 of the DFM cost reference,
                               where 'C' denotes a flat/thin wall section
                               component and the number is its complexity
                               band (1 = simplest)
    * tolerance_mm          - tightest critical total tolerance (mm)
    * surface_finish_um_ra  - tightest critical surface finish (um Ra)

and, for each fabrication process, the resulting relative cost
coefficients:
    * Cc - relative cost associated with geometrical (shape) complexity
    * Cs - relative cost associated with size/section thickness
    * Ct - relative cost associated with tolerance
    * Cf - relative cost associated with surface finish

Coefficients are looked up per process rather than derived from a
continuous regression: digitizing the full family of Cc/Cs/Ct/Cf curves
against these geometry parameters isn't supported by any data we have
for our specific processes (Hot Rolling, Cold Rolling, HIP, CNC,
Electron Beam, Diffusion Bonding, Spray Deposition) - most of them
aren't covered by the generic DFM cost reference at all. Only
pre-defined ('database') geometries are supported for now; there is no
free-form 'custom geometry' input yet.

########################################################################
CLASSES
########################################################################
Geometry(name: str)
    Represents a component's geometry.

########################################################################
FUNCTIONS
########################################################################
------------------------------------------------------------------------
build_geometry(name, volume_mm3, section_thickness_mm, shape_class,
               tolerance_mm, surface_finish_um_ra, Cc_map, Cs_map,
               Ct_map, Cf_map)
------------------------------------------------------------------------
load_geometries(geometry_dir)
"""

# ######################################################################
# CLASSES
# ######################################################################
# ----------------------------------------------------------------------
# Geometry class

class Geometry:
    """
    Represents a component's geometry.


    ########################################################################
    ATTRIBUTES
    ########################################################################
    ------------------------------------------------------------------------
    name : str
        Name of the geometry (e.g., 'Plasma Facing Surface')

    ------------------------------------------------------------------------
    volume_mm3 : float
        Component volume (mm^3)

    ------------------------------------------------------------------------
    section_thickness_mm : float
        Representative section/wall thickness (mm)

    ------------------------------------------------------------------------
    shape_class : str
        Shape classification band (e.g. 'C1')

    ------------------------------------------------------------------------
    tolerance_mm : float
        Tightest critical total tolerance (mm)

    ------------------------------------------------------------------------
    surface_finish_um_ra : float
        Tightest critical surface finish (um Ra)

    ------------------------------------------------------------------------
    Cc_map, Cs_map, Ct_map, Cf_map : dict
        Dictionaries mapping process name to the relative cost
        coefficient for that process, given this geometry.

        Cc_map[process] = Cc, etc.

        Where:
            process : str
                String representing a manufacturing process (e.g.,
                'CNC', 'Hot Rolling').

            Cc / Cs / Ct / Cf : float
                Numeric (float) value representing the relative cost of
                producing this geometry with the specified process
                compared to an ideal design. A value of 1 indicates the
                same cost as the ideal, while values greater than 1
                indicate higher costs.

    ########################################################################
    METHODS
    ########################################################################
    ------------------------------------------------------------------------
    add_Cc(process, Cc) / add_Cs(process, Cs) / add_Ct(process, Ct) / add_Cf(process, Cf)
        Add or set a relative cost coefficient for a given process.

    ------------------------------------------------------------------------
    get_Cc(process) / get_Cs(process) / get_Ct(process) / get_Cf(process)
        Look up a relative cost coefficient for a given process.

    ------------------------------------------------------------------------
    get_coefficients(process)
        Look up Cc, Cs, Ct, and Cf for a given process.
    """

    def __init__(self, name: str):
        self.name = name
        self.volume_mm3 = None
        self.section_thickness_mm = None
        self.shape_class = None
        self.tolerance_mm = None
        self.surface_finish_um_ra = None
        self.Cc_map = {}
        self.Cs_map = {}
        self.Ct_map = {}
        self.Cf_map = {}

    def _add_coefficient(self, coeff_map: dict, process: str, value: float, label: str):
        # If process already exists, warn user of override
        # TODO: Ask user if they want to proceed
        if process in coeff_map:
            print(f"Warning: process '{process}' already exists. Overriding " f"{coeff_map[process]} {label} with new value of {value}")

        coeff_map[process] = value

    def add_Cc(self, process: str, Cc: float):
        """
        Add or set the shape complexity coefficient (Cc) for a process.
        """
        self._add_coefficient(self.Cc_map, process, Cc, "Cc")

    def add_Cs(self, process: str, Cs: float):
        """
        Add or set the size/section thickness coefficient (Cs) for a process.
        """
        self._add_coefficient(self.Cs_map, process, Cs, "Cs")

    def add_Ct(self, process: str, Ct: float):
        """
        Add or set the tolerance coefficient (Ct) for a process.
        """
        self._add_coefficient(self.Ct_map, process, Ct, "Ct")

    def add_Cf(self, process: str, Cf: float):
        """
        Add or set the surface finish coefficient (Cf) for a process.
        """
        self._add_coefficient(self.Cf_map, process, Cf, "Cf")

    def _get_coefficient(self, coeff_map: dict, process: str, label: str) -> float:
        if process not in coeff_map:
            raise KeyError(
                f"Process '{process}' not found in {label}_map for geometry '{self.name}'."
            )
        return coeff_map[process]

    def get_Cc(self, process: str) -> float:
        """
        Look up the shape complexity coefficient (Cc) for a given process.
        """
        return self._get_coefficient(self.Cc_map, process, "Cc")

    def get_Cs(self, process: str) -> float:
        """
        Look up the size/section thickness coefficient (Cs) for a given process.
        """
        return self._get_coefficient(self.Cs_map, process, "Cs")

    def get_Ct(self, process: str) -> float:
        """
        Look up the tolerance coefficient (Ct) for a given process.
        """
        return self._get_coefficient(self.Ct_map, process, "Ct")

    def get_Cf(self, process: str) -> float:
        """
        Look up the surface finish coefficient (Cf) for a given process.
        """
        return self._get_coefficient(self.Cf_map, process, "Cf")

    def get_coefficients(self, process: str) -> dict:
        """
        Look up Cc, Cs, Ct, and Cf for a given process.


        ########################################################################
        RETURNS
        ########################################################################
        ------------------------------------------------------------------------
        dict
            {"Cc": Cc, "Cs": Cs, "Ct": Ct, "Cf": Cf} for the given process.
        """
        return {
            "Cc": self.get_Cc(process),
            "Cs": self.get_Cs(process),
            "Ct": self.get_Ct(process),
            "Cf": self.get_Cf(process),
        }

# ######################################################################
# FUNCTIONS
# ######################################################################
# ----------------------------------------------------------------------
# Function to build a 'Geometry' object from input data.

def build_geometry(name, volume_mm3, section_thickness_mm, shape_class,
                    tolerance_mm, surface_finish_um_ra,
                    Cc_map, Cs_map, Ct_map, Cf_map, verbose=True):
    """
    Build a 'Geometry' object from input data.


    ########################################################################
    PARAMETERS
    ########################################################################
    ------------------------------------------------------------------------
    name : str
        Name of the geometry (e.g., 'Plasma Facing Surface')

    ------------------------------------------------------------------------
    volume_mm3 : float
        Component volume (mm^3)

    ------------------------------------------------------------------------
    section_thickness_mm : float
        Representative section/wall thickness (mm)

    ------------------------------------------------------------------------
    shape_class : str
        Shape classification band (e.g. 'C1')

    ------------------------------------------------------------------------
    tolerance_mm : float
        Tightest critical total tolerance (mm)

    ------------------------------------------------------------------------
    surface_finish_um_ra : float
        Tightest critical surface finish (um Ra)

    ------------------------------------------------------------------------
    Cc_map, Cs_map, Ct_map, Cf_map : dict
        Dictionaries mapping process name to the relative cost
        coefficient for that process, given this geometry.
    """

    # ------------------------
    # Validate name
    if not isinstance(name, str) or name.strip() == "":
        print_error(f"Error loading geometry, '{name}'. Geometry skipped. Geometry name must be a non-empty string.")
        return None

    # ------------------------
    # Validate volume
    try:
        volume_mm3 = positive_float(volume_mm3)
    except Exception:
        print_error(f"Error loading geometry, '{name}'. Geometry skipped. Volume must be a positive number.")
        return None

    # ------------------------
    # Validate section thickness
    try:
        section_thickness_mm = positive_float(section_thickness_mm)
    except Exception:
        print_error(f"Error loading geometry, '{name}'. Geometry skipped. Section thickness must be a positive number.")
        return None

    # ------------------------
    # Validate shape class
    if not isinstance(shape_class, str) or shape_class.strip() == "":
        print_error(f"Error loading geometry, '{name}'. Geometry skipped. Shape class must be a non-empty string.")
        return None

    # ------------------------
    # Validate tolerance
    try:
        tolerance_mm = positive_float(tolerance_mm)
    except Exception:
        print_error(f"Error loading geometry, '{name}'. Geometry skipped. Tolerance must be a positive number.")
        return None

    # ------------------------
    # Validate surface finish
    try:
        surface_finish_um_ra = positive_float(surface_finish_um_ra)
    except Exception:
        print_error(f"Error loading geometry, '{name}'. Geometry skipped. Surface finish must be a positive number.")
        return None

    # ------------------------------------------
    # Validate coefficient maps
    for map_name, coeff_map in [("Cc_map", Cc_map), ("Cs_map", Cs_map), ("Ct_map", Ct_map), ("Cf_map", Cf_map)]:
        for process, value in coeff_map.items():
            try:
                coeff_map[process] = positive_float(value)
            except Exception:
                print_error(f"Error loading geometry, '{name}'. Geometry skipped. {map_name} value for process, '{process}', must be a positive number.")
                return None

    # -----------------------
    # Create 'Geometry' object
    geom = Geometry(name)
    geom.volume_mm3 = volume_mm3
    geom.section_thickness_mm = section_thickness_mm
    geom.shape_class = shape_class
    geom.tolerance_mm = tolerance_mm
    geom.surface_finish_um_ra = surface_finish_um_ra
    geom.Cc_map = dict(Cc_map)
    geom.Cs_map = dict(Cs_map)
    geom.Ct_map = dict(Ct_map)
    geom.Cf_map = dict(Cf_map)

    # ------------------------
    # Display info
    if verbose:
        print("\n➤  ", name)
        print(f"       Volume:              {volume_mm3} mm^3")
        print(f"       Section thickness:   {section_thickness_mm} mm")
        print(f"       Shape class:         {shape_class}")
        print(f"       Tolerance:           {tolerance_mm} mm")
        print(f"       Surface finish:      {surface_finish_um_ra} um Ra")

    return geom

# ----------------------------------------------------------------------
# Function to load geometries from a geometries directory.

import os
import importlib.util

def load_geometries(geometry_dir):
    """
    Load geometries from a directory of geometry input files.


    ########################################################################
    PARAMETERS
    ########################################################################
    ------------------------------------------------------------------------
    geometry_dir : str
        Name of folder containing geometry input files

    ########################################################################
    RETURNS
    ########################################################################
    ------------------------------------------------------------------------
    geometries :
        List of 'Geometry' objects from the geometries directory
    """
    geometries = []

    if not os.path.isdir(geometry_dir):
        print(f"\nDirectory '{geometry_dir}' not found. Skipping geometry loading.")
        return geometries

    for filename in os.listdir(geometry_dir):

        if not filename.endswith(".py"):
            continue

        filepath = os.path.join(geometry_dir, filename)
        module_name = os.path.splitext(filename)[0]

        # Load module
        spec = importlib.util.spec_from_file_location(module_name, filepath)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        # Extract data
        name = module.name
        volume_mm3 = module.volume_mm3
        section_thickness_mm = module.section_thickness_mm
        shape_class = module.shape_class
        tolerance_mm = module.tolerance_mm
        surface_finish_um_ra = module.surface_finish_um_ra
        Cc_map = module.Cc_map
        Cs_map = module.Cs_map
        Ct_map = module.Ct_map
        Cf_map = module.Cf_map

        # Build geometry
        geom = build_geometry(
            name, volume_mm3, section_thickness_mm, shape_class,
            tolerance_mm, surface_finish_um_ra,
            Cc_map, Cs_map, Ct_map, Cf_map,
        )

        geometries.append(geom)

    return geometries

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
