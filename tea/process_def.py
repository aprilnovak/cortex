"""
Process Definition


########################################################################
DESCRIPTION
########################################################################
Defines the 'Process' class and functions to build processes from input
data.

########################################################################
CLASSES
########################################################################
Process(name: str, alphaT: float, beta: float)
    Represents a manufacturing process.

########################################################################
FUNCTIONS
########################################################################
"""

# ######################################################################
# CLASSES
# ######################################################################
# ----------------------------------------------------------------------
# Process class

class Process:
    """
    Represents a manufacturing process.

    
    ########################################################################
    DESCRIPTION
    ########################################################################
    A process is defined by the following equation:

    Pc = α*T + β/N

    Where ideal process cost (Pc) is determined from operating and overhead
    cost per sec (α), time (T), and tooling cost (β).

    ########################################################################
    ATTRIBUTES
    ########################################################################
    ------------------------------------------------------------------------
    name : str
        Name of the process (e.g., 'CNC', 'Hot Rolling')
    
    ------------------------------------------------------------------------
    alphaT : float
        Product of operating and overhead cost per sec (α) and time (T)

    ------------------------------------------------------------------------
    beta : float
        Tooling cost (β)

    ------------------------------------------------------------------------
    description : str
        Optional description of the process 

    ########################################################################
    METHODS
    ########################################################################
    ------------------------------------------------------------------------
    basic_processing_cost(production_qty: int) -> float
        Calculate basic processing cost (Pc)
    """

    def __init__(self, name: str, alphaT: float, beta: float, description: str = ""):
        
        self.name = name
        self.alphaT = alphaT 
        self.beta = beta                 
        self.description = description

    def basic_processing_cost(self, production_qty: int) -> float:
        """
        Calculate basic processing cost (Pc).

        
        ########################################################################
        DESCRIPTION
        ########################################################################
        Cost is calculated with the following equation:

        Pc = α*T + β/N

        Where basic process cost (Pc) is determined from operating and overhead
        cost per sec (α), time (T), and tooling cost (β).

        ########################################################################
        PARAMETERS
        ########################################################################
        ------------------------------------------------------------------------
        production_qty : int
            Production quantity (N)
         
        ########################################################################
        RETURNS
        ########################################################################
        ------------------------------------------------------------------------
        Pc : float
            Basic processing cost (USD)
        """

        Pc = self.alphaT + self.beta / production_qty

        return Pc

    def __repr__(self):
        return f"<Process: {self.name}>"
    
# ######################################################################
# FUNCTIONS
# ######################################################################
# ----------------------------------------------------------------------
# Function to build a 'Process' object from input data.

def build_process(name, alphaT, beta, description="",verbose=True):
    """
    Build a 'Process' object from input data.

    
    ########################################################################
    PARAMETERS
    ########################################################################
    ------------------------------------------------------------------------
    name : str
        Process name (e.g., 'CNC')

    ------------------------------------------------------------------------
    alphaT : float
        Time-related cost coefficient

    ------------------------------------------------------------------------
    beta : float
        Fixed or scaling cost parameter

    ------------------------------------------------------------------------
    description : str, optional
        Description of the process
    ------------------------------------------------------------------------
    """
    
    # ------------------------
    # Validate name
    if not isinstance(name, str) or name.strip() == "":
        print_error(f"Error loading proces, '{name}'. Process skipped. Process name must be a non-empty string.")
        return None

    # ------------------------
    # Validate alphaT
    try:
        alphaT = positive_float(alphaT)
    except Exception:
        print_error(f"Error loading proces, '{name}'. Process Skipped. αT must be a positive number.")
        return None

    # ------------------------
    # Validate beta
    try:
        beta = positive_float(beta)
    except Exception:
        print_error(f"Error loading proces, '{name}'. Process Skipped. β must be a positive number.")
        return None

    # ---------------------
    # Create process object
    proc = Process(
        name=name,
        alphaT=alphaT,
        beta=beta,
        description=description
    )

    # ------------------------
    # Display info
    if verbose:
        print("\n➤  ", name)
        print(f"       alphaT:      {alphaT}")
        print(f"       beta:        {beta}")
        if description:
            print(f"       Description: {description}")

    return proc

# ----------------------------------------------------------------------
# Function to load custom processes from processes directory.

import os
import importlib.util

def load_processes(process_dir):
    """
    Load custom processes from processes directory.


    ########################################################################
    PARAMETERS
    ########################################################################
    ------------------------------------------------------------------------
    process_dir : str
        Name of folder containing process input files

    ########################################################################
    RETURNS
    ########################################################################
    ------------------------------------------------------------------------
    processes : 
        List of 'Process' objcets from processes directory
    """

    processes = []

    for filename in os.listdir(process_dir):

        if not filename.endswith(".py"):
            continue

        filepath = os.path.join(process_dir, filename)
        module_name = os.path.splitext(filename)[0]

        # Load module
        spec = importlib.util.spec_from_file_location(module_name, filepath)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        # Extract data (required fields)
        name = module.name
        alphaT = module.alphaT
        beta = module.beta

        # Optional description
        description = getattr(module, "description", "")

        # Build process object
        proc = build_process(name, alphaT, beta, description)

        processes.append(proc)


    return processes

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