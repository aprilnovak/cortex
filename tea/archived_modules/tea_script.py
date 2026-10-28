"""
Main program


########################################################################
DESCRIPTION
########################################################################
The main execution script. It handles command-line input, orchestrates 
the workflow, performs input validation, and coordinates interaction 
between all modules.
"""

import sys
import textwrap
from material_def import load_materials
from material_def import build_material
from material_def import blend_materials
from component_def import load_geometries
from tea1 import Component
import importlib
from types import SimpleNamespace
import argparse
import process_def as process
import difflib

GEOMETRIES_DB_DIR = "geometries_database"

# ----------------------------------------------------------------------
# Formating

def print_error(msg, width=72):

    print("\n        ---------------------- ✖ ERROR ✖ -----------------------        ")

    wrapped = textwrap.wrap(msg, width - 20)

    for line in wrapped:
        print(f"        |{line.center(width - 18)}|")

    # blank spacer line
    print(f"        |{' ' * (width -18)}|")
    print("        " + "-" * (width - 16) + "")

def print_invalid(msg, width=72):
    print("    |")
    print("    |     ------------------ ✖ INVALID INPUT ✖ -------------------      ")

    wrapped = textwrap.wrap(msg, width - 20)

    for line in wrapped:
        print(f"    |     |{line.center(width - 18)}|")

    # blank spacer line
    print(f"    |     |{' ' * (width -18)}|")
    print("    |     " + "-" * (width - 16) + "")
    print("    |")


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


print("\n═══✦════════════════════════════════════════════════════════════════✦═══")
print(  "                  CORTEX Design for Manufacture Toolkit                  ")
print(  "═══✦════════════════════════════════════════════════════════════════✦═══\n")

# ######################################################################
# LOAD INPUTS
# ######################################################################

# ----------------------------------------------------------------------
# Custom ValueError functions
def positive_float(value):
    v = float(value)
    if v <= 0:
        raise ValueError("Please enter a positive numeric value.")
    return v


def positive_int(value):
    v = int(value)
    if v <= 0:
        raise ValueError("Value must be a positive integer.")
    return v

# ----------------------------------------------------------------------
# Load file
print("──────────────────────── Loading toolkit inputs ────────────────────────\n")
# ----------------------
# Look for an input file

parser = argparse.ArgumentParser()

parser.add_argument(
    "-i",
    "--input",
    dest="input_file_name",
    default=None,
    help="Name of input file (without .py)"
)

args = parser.parse_args()

# --------------------------------------
# Try to import input module if provided

inputs = None

if args.input_file_name:
    name = args.input_file_name.replace(".py", "")
    module_name = f"inputs_components.{name}"

    try:
        inputs = importlib.import_module(module_name)
        print(f"Input file,'{module_name}', loaded.")

    except ModuleNotFoundError:
        print(f"Input file, '{module_name}', not found.\n")

else:
    print("No input file specified\n")

# ----------------------------------------------------------------------
# If import failed OR no file given, prompt user

if inputs is None:
    print("Switching to manual input mode.")
    print("\n    ----------------- Please enter required inputs ---------------------")

    # Build a simple object to store user inputs
    inputs = SimpleNamespace()

    # Custom materials directory input
    print("    |")
    inputs.material_dir = input("    | Material directory (or leave blank): ")
    if inputs.material_dir.strip() == "":
        inputs.material_dir = None

    # Material name input
    # Cannot be empty
    while True:
        print("    |")
        inputs.material = input("    | Material name: ").strip()
        if inputs.material == "":
            print_invalid("Material name cannot be empty. Please enter a name.") # doesn't break loop and material name is overriden on next pass
        else:
            break

    # Processes list input
    # Cannot be empty
    inputs.processes = []
    print("    |")
    print("    | Enter processes (type 'done' to finish):")

    while True:
        pname = input("    |     Process name: ").strip()

        if pname.lower() == "done":
            break

        if pname == "":
            print_invalid("Process name cannot be empty. Please enter a name.")
            continue # goes back to beginng of loop skipping adding process to list

        inputs.processes.append(pname)

    # Custom process directory input
    print("    |")
    inputs.process_dir = input("    | Process directory (or leave blank): ")
    if inputs.process_dir.strip() == "":
        inputs.process_dir = None

    # Component name input
    # Cannot be empty
    while True:
        print("    |")
        inputs.component_name = input("    | Component name: ").strip()
        if inputs.component_name == "":
            print_invalid("Component name cannot be empty. Enter a name.") # doesn't break loop and material name is overriden on next pass
        else:
            break

    # Volume input
    # Must be positive numeric value
    while True:
        try:
            print("    |")
            inputs.volume_mm3 = positive_float(input("    | Volume (mm^3): "))
            break
        except ValueError as e:
            print_invalid(f"Please enter a positive number. Try again.")

    # Production quantity input
    while True:
        try:
            print("    |")
            inputs.production_qty = positive_int(input("    | Production quantity: "))
            break
        except ValueError as e:
            print_invalid(f"Please enter a positive integer. Try again.")

    # Cost factors
    n_proc = len(inputs.processes)
    print("    |")
    print(f"    | Enter cost factors for each of the processes:")

    def get_list_input(name):
        values = []
        for pname in inputs.processes:
            while True:
                try:
                    val = positive_float(input(f"    |     {name} for process, '{pname}': "))
                    values.append(val)
                    break
                except ValueError:
                    print_invalid(f"Please enter a positive number. Try again.")
        return values

    inputs.Cc = get_list_input("Cc")
    inputs.Cs = get_list_input("Cs")
    inputs.Ct = get_list_input("Ct")
    inputs.Cf = get_list_input("Cf")
    inputs.Wc = get_list_input("Wc")
    print("    |")
    print("    --------------------------------------------------------------------")

# ----------------------------------------------------------------------
# Verify file inputs validity

allowed_attrs = {
    "material",
    "material_dir",
    "material_mode",
    "material_a",
    "material_b",
    "blend_fraction_a",
    "blend_basis",
    "blend_cmp",
    "processes",
    "process_dir",
    "component_name",
    "geometry",
    "geometry_dir",
    "volume_mm3",
    "production_qty",
    "Cc",
    "Cs",
    "Ct",
    "Cf",
    "Wc"
}

# Get attributes from module (exclude built-ins
input_attrs = {a for a in vars(inputs).keys() if not a.startswith("__")}

# Find unknown attributes
unknown_attrs = input_attrs - allowed_attrs

if unknown_attrs:
    msgs = []

    for attr in unknown_attrs:
        suggestion = difflib.get_close_matches(attr, allowed_attrs, n=1)

        if suggestion:
            msgs.append(
                f"Warning: input file attribute, '{attr}', not recognized. "
                f"Did you mean '{suggestion[0]}'?"
            )
        else:
            msgs.append(
                f"Warning: input file attribute, '{attr}', not recognized "
                f"and will not be used."
            )

    warn_block(msgs)

    

# If blending two materials, validate the blend attributes now and give
# 'inputs.material' a placeholder name so the existing single-material
# checks/lookup below keep working unchanged; the real blended material
# is built and registered under this name in the MATERIALS section.
material_mode = getattr(inputs, "material_mode", "single")

if material_mode == "blend":
    missing = [a for a in ("material_a", "material_b", "blend_fraction_a") if not hasattr(inputs, a)]
    if missing:
        print_error(textwrap.fill(
            f"Blended material requires {', '.join(missing)}. Edit input file, '{module_name}', and try again."
        ))
        print("\n════════════════════════ CALCULATION TERMINATED ════════════════════════\n")
        sys.exit(1)
    inputs.material = f"{inputs.material_a}/{inputs.material_b} Blend"

# Material name must not be empty
if not hasattr(inputs, "material") or not str(inputs.material).strip():
    print_error(textwrap.fill(
        f"Material name cannot be empty. Edit attribute, 'material_name', in input file, '{module_name}', and try again."
    ))
    print("\n════════════════════════ CALCULATION TERMINATED ════════════════════════\n")
    sys.exit(1)

# Component name must not be empty
if not hasattr(inputs, "component_name") or not str(inputs.component_name).strip():
    print_error(textwrap.fill(
        f"Component name cannot be empty. Edit attribute, 'component_name', in input file, '{module_name}', and try again."
    ))
    print("\n════════════════════════ CALCULATION TERMINATED ════════════════════════\n")
    sys.exit(1)

# Processes must exist and not be empty
if not hasattr(inputs, "processes") or len(inputs.processes) == 0:
    print_error(textwrap.fill(
        f"At least one process must be specified. Edit attribute 'processes', in input file, '{module_name}', and try again."
    ))
    print("\n════════════════════════ CALCULATION TERMINATED ════════════════════════\n")
    sys.exit(1)

# Ensure no empty process names
for pname in inputs.processes:
    if not str(pname).strip():
        print_error(textwrap.fill(
            f"Process names cannot be empty. Edit attribute 'processes', in input file, '{module_name}', and try again."
        ))
        print("\n════════════════════════ CALCULATION TERMINATED ════════════════════════\n")
        sys.exit(1)

# Resolve geometry (optional) and, from it, default volume and any
# per-process cost coefficients (Cc/Cs/Ct/Cf) not explicitly given in the
# input file. Explicit values in the input file always take precedence.
geometry = None

if getattr(inputs, "geometry", None):
    geometry_lookup = {}

    db_geometries = load_geometries(GEOMETRIES_DB_DIR)
    geometry_lookup.update({g.name: g for g in db_geometries if g is not None})

    geometry_dir = getattr(inputs, "geometry_dir", None)
    if geometry_dir:
        custom_geometries = load_geometries(geometry_dir)
        geometry_lookup.update({g.name: g for g in custom_geometries if g is not None})

    if inputs.geometry not in geometry_lookup:
        print_error(textwrap.fill(
            f"Geometry, '{inputs.geometry}', not found. Edit attribute, 'geometry', in input file, '{module_name}', and try again."
        ))
        print("\n════════════════════════ CALCULATION TERMINATED ════════════════════════\n")
        sys.exit(1)

    geometry = geometry_lookup[inputs.geometry]

if not hasattr(inputs, "volume_mm3"):
    if geometry is None:
        print_error(textwrap.fill(
            f"Volume must be specified (directly via 'volume_mm3', or indirectly via 'geometry'). Edit input file, '{module_name}', and try again."
        ))
        print("\n════════════════════════ CALCULATION TERMINATED ════════════════════════\n")
        sys.exit(1)
    inputs.volume_mm3 = geometry.volume_mm3

for attr_name, geometry_map_attr in [("Cc", "Cc_map"), ("Cs", "Cs_map"), ("Ct", "Ct_map"), ("Cf", "Cf_map")]:
    manual = getattr(inputs, attr_name, None)
    geom_map = getattr(geometry, geometry_map_attr) if geometry is not None else {}
    setattr(inputs, attr_name, [
        geom_map[pname] if pname in geom_map else (manual[i] if manual is not None else 1.0)
        for i, pname in enumerate(inputs.processes)
    ])

# Volume must be positive
try:
    inputs.volume_mm3 = positive_float(inputs.volume_mm3)
except Exception:
    print_error(textwrap.fill(
        f"Volume must be a positive numeric value. Edit attribute, 'volume_mm3', in input file, '{module_name}', and try again."
    ))
    print("\n════════════════════════ CALCULATION TERMINATED ════════════════════════\n")
    sys.exit(1)

# Production quantity must be positive integer
try:
    inputs.production_qty = positive_int(inputs.production_qty)
except Exception:
    print_error(textwrap.fill(
        f"Production quantity must be a positive integer. Edit attribute, 'production_qty', in input file, '{module_name}', and try again."
    ))
    print("\n════════════════════════ CALCULATION TERMINATED ════════════════════════\n")
    sys.exit(1)

# Cost factors must match number of processes
n = len(inputs.processes)
for arr in [inputs.Cc, inputs.Cs, inputs.Ct, inputs.Cf, inputs.Wc]:
    if len(arr) != n:
        print_error(textwrap.fill(f"Cost factors must match number of processes. Edit Input file, '{module_name}', and try again."))
        print("\n════════════════════════ CALCULATION TERMINATED ════════════════════════\n")
        sys.exit(1)

# ######################################################################
# MATERIALS
# ######################################################################

print("\n────────────────────────── Loading materials ───────────────────────────")

# Build process lookup
material_lookup = {}

# ----------------------------------------------------------------------
# Load defualt materials

print(f"\nPre-defined materials from materials database:")

database_dir = "materials_database"

db_materials = load_materials(database_dir)

# Add to lookup
material_lookup.update({m.name: m for m in db_materials})

# ----------------------------------------------------------------------
# Load materials from directory

if inputs.material_dir:
    print(f"\nCustom materials from '{inputs.material_dir}' directory:")
    custom_materials = load_materials(inputs.material_dir) 
    valid_materials = [m for m in custom_materials if m is not None]
    material_lookup.update({m.name: m for m in valid_materials})

else:
    print("\nNo custom materials directory specified.")

# ----------------------------------------------------------------------
# If blending two materials, build the blend now and register it under
# 'inputs.material' (set above) so the lookup below resolves to it
# unchanged, exactly like any other pre-loaded material.

if material_mode == "blend":

    if inputs.material_a not in material_lookup:
        print_error(textwrap.fill(f"Material, '{inputs.material_a}', not found in loaded materials for blending."))
        print("\n════════════════════════ CALCULATION TERMINATED ════════════════════════\n")
        sys.exit(1)

    if inputs.material_b not in material_lookup:
        print_error(textwrap.fill(f"Material, '{inputs.material_b}', not found in loaded materials for blending."))
        print("\n════════════════════════ CALCULATION TERMINATED ════════════════════════\n")
        sys.exit(1)

    blended = blend_materials(
        name=inputs.material,
        material_a=material_lookup[inputs.material_a],
        material_b=material_lookup[inputs.material_b],
        fraction_a=inputs.blend_fraction_a,
        basis=getattr(inputs, "blend_basis", "volume"),
    )

    if blended is None:
        print_error(textwrap.fill(f"Failed to blend materials '{inputs.material_a}' and '{inputs.material_b}'."))
        print("\n════════════════════════ CALCULATION TERMINATED ════════════════════════\n")
        sys.exit(1)

    # A blend has no Cmp_map of its own; apply any compatibility
    # overrides explicitly supplied by the input file
    for pname, Cmp in getattr(inputs, "blend_cmp", {}).items():
        blended.add_Cmp(pname, Cmp)

    material_lookup[blended.name] = blended

# ----------------------------------------------------------------------
# Lookup or define material

# ------------------------
# Case 1: Found in lookup
if inputs.material in material_lookup:
    material = material_lookup[inputs.material]

# ----------------------------------------------------------
# Case 2: material cannot be found, prompt user to define it
else:
    print(f"\nMaterial, '{inputs.material}', not found in loaded materials.")
    print("\n    ------------------ Please define the material ----------------------")

    # --------------------------
    # Prompt user for properties
    while True:
        try:
            print("    |")
            density = positive_float(input("    | Enter density (g/cm^3): "))
            break 
        except ValueError:
            print_invalid("Please enter a positive number. Try again.")

    # ---------------------------
    # Prompt user for composition

    # Alloy and residual elements
    composition = {}

    print("    |")
    print("    | Enter Alloying Elements (type 'done' to finish):")

    while True:
        # Get symobol
        # Cannot be empty
        while True:
            element = input("    |     Element symbol (e.g., Fe, Cr): ").strip()
            if element == "":
                print_invalid("Element symbol cannot be empty. Please enter a symbol.") # doesn't break loop and element is overriden on next pass
            else:
                break
        

        if element.lower() == "done":
            break

        # Get valid weight %
        while True:
            try:
                wt = positive_float(input(f"    |         Weight % of {element}: "))
                break
            except ValueError:
                print_invalid("Please enter a postitive number. Try again")

        # Get valid type
        while True:
            el_type = input("    |         Type (alloy/residual): ").strip().lower()

            if el_type in ["alloy", "residual"]:
                break
            else:
                print_invalid("Invalid type. Must be 'alloy' or 'residual'")

        # Store result
        composition[element] = {"wt": wt, "type": el_type}

    # Remainder element
    print("    |")

    while True:
        remainder_element = input("    | Enter remainder element (e.g., Fe): ").strip()
        if remainder_element == "":
            print_invalid("Element symbol cannot be empty. Please enter a symbol.") # doesn't break loop and element is overriden on next pass
        else:
            break

    # ------------------------------------
    # Prompt user for Cmp for each process
    Cmp_map = {}
    print("    |")
    print("    | Enter Cmp values for each process:")

    for pname in inputs.processes:
        while True:
            try:
                Cmp = positive_float(input(f"    |    Cmp for process, '{pname}': "))
                break
            except ValueError:
                print_invalid(f"Please enter a postitive number. Try again")

        Cmp_map[pname] = Cmp

    print("    |")
    print("    --------------------------------------------------------------------")

    # --------------
    # Build material
    material = build_material(
        name=inputs.material,
        density=density,
        composition=composition,
        remainder_element=remainder_element,
        Cmp_map=Cmp_map
    )









# ######################################################################
# PROCESSES
# ######################################################################

print("\n──────────────────── Loading fabrication processes ─────────────────────")

processes = []

# Build process lookup
process_lookup = {}

# ----------------------------------------------------------------------
# Load defualt processes
print(f"\nPre-defined processes from process database:")

database_dir = "processes_database"

db_processes = process.load_processes(database_dir)

# Add to lookup
process_lookup.update({p.name: p for p in db_processes})

# ----------------------------------------------------------------------
# Load from custom directory (if provided)

if inputs.process_dir:
    print(f"\nCustom processes from '{inputs.process_dir}' directory:")

    custom_processes = process.load_processes(inputs.process_dir)
    valid_processes = [p for p in custom_processes if p is not None]

    # Add to process lookup
    process_lookup.update({p.name: p for p in valid_processes})

    
    
else:
    print("\nNo custom processes directory specified.")

# ----------------------------------------------------------------------
# Lookup or define processes from inputs file

for pname in inputs.processes:

    # -----------------------
    # Case 1: Found in lookup
    if pname in process_lookup:
        
        proc = process_lookup[pname]
        #Tagging these processes as loaded (used later in verification)
        proc._source = "loaded"
        processes.append(proc)

        continue

    # ---------------------------------------------------------
    # Case 2: Process cannot be found, prompt user to define it
    print(f"\nProcess, '{pname}', not found in loaded processes.")
    print("\n    ------------------- Please define the process ----------------------")
    print("    |")

    # ---------------------------------
    # Prompt user for cost correlations

    # Get valid alphaT
    while True:
        try:
            alphaT = positive_float(input(f"    | Enter αT for process, '{pname}': "))
            break
        except ValueError:
            print_invalid("Please enter a positive number.")

    # Get valid beta
    while True:
        try:
            print("    |")
            beta = positive_float(input(f"    | Enter β for process, '{pname}': "))
            break
        except ValueError:
            print_invalid("Please enter a positive number.")

    print("    |")
    print("    -------------------------------------------------------------------")
    proc = process.build_process(
        name=pname,
        alphaT=alphaT,
        beta=beta
    )

    # Tagging these processes as user defined (used later in verification)
    proc._source = "user"

    processes.append(proc)

# ######################################################################
# VERIFICATION
# ######################################################################
# ----------------------------------------------------------------------
# Ensure Cmp coverage for all material–process combinations

print("\n─────────────────────── Verifying input validity ───────────────────────\n")

for proc in processes:

    if proc.name in material.Cmp_map:
        continue  # already defined

    # If process came from loaded processes
    is_loaded_process = proc._source == "loaded"

    if is_loaded_process:
        print_error(textwrap.fill(f"According to the loaded input files, material, '{material.name}', is not compatible with process, '{proc.name}'."))
        print("\n══════════════════════════ CALCULATION TERMINATED ══════════════════════")
        print()
        sys.exit(0)

    print(textwrap.fill(f"Material, '{material.name}', is missing Cmp value for newly added process, '{proc.name}'.\n"))

    print("\n        -------------------- Add Cmp value ---------------------        \n")
    # Otherwise prompt user
    while True:
        val = input(f"        Enter Cmp for material, '{material.name}': ").strip()

        try:
            Cmp = float(val)
            material.add_Cmp(proc.name, Cmp)
            break
        except ValueError:
            print("Invalid input. Enter a number.")

    print("\n        --------------------------------------------------------        \n")

print("Inputs are valid.")
# ######################################################################
# BUILD COMPONENT
# ######################################################################

component = Component(
    name=inputs.component_name,
    volume_mm3=inputs.volume_mm3,
    material=material,
    processes=processes,
    production_qty=inputs.production_qty,
    Cc=inputs.Cc,
    Cs=inputs.Cs,
    Ct=inputs.Ct,
    Cf=inputs.Cf,
    Wc=inputs.Wc
)

# ######################################################################
# GENERATE SUMMARY
# ######################################################################

summary = component.manufacturing_cost()

# ---------------
# Display summary
print("\n══════════════════════════ COMPONENT SUMMARY ═══════════════════════════")
print(f"\nComponent: {summary['Component']}")
print(f"\nMaterials: ")
print(f"    {summary['Material']}: Cost = ${material.cost:.2f}/kg, Mass = {summary['Mass']}, Wc = {summary['Wc']}")

print("\nProcesses:")
for p in summary["Processes"]:
    print(
        f"    {p['Process']}: "
        f"Pc = {p['Pc']:.2f}, Rc = {p['Rc']:.3f}, Cost = ${p['Cost']:.2f}"
    )

print("\nCost Breakdown:")
for k, v in summary["Cost Breakdown"].items():
    print(f"    {k}: ${v:.2f}")

print("\nTotal Manufacturing Cost:")
print(f"    ${summary['Total cost']:.2f} per component")
print(f"    ${summary['Unit cost']:.2f} per kg\n")
