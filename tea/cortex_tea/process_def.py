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

    ------------------------------------------------------------------------
    tooling_level, equipment_level, time_level : str or None
        The COST_LEVELS entries alphaT/beta were computed from, if the
        process was built from categorical levels rather than direct
        alphaT/beta values. None otherwise.

    ########################################################################
    METHODS
    ########################################################################
    ------------------------------------------------------------------------
    basic_processing_cost(production_qty: int) -> float
        Calculate basic processing cost (Pc)
    """

    def __init__(self, name: str, alphaT: float, beta: float, description: str = "",
                 tooling_level=None, equipment_level=None, time_level=None):

        self.name = name
        self.alphaT = alphaT
        self.beta = beta
        self.description = description
        self.tooling_level = tooling_level
        self.equipment_level = equipment_level
        self.time_level = time_level

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
# CATEGORICAL COST LEVELS

# Instead of specifying alphaT/beta directly, a process can be described by
# a Low..Very High level for tooling cost, equipment cost, and processing
# time. Each level maps to a value below (from the GUI prototype's cost
# table), and the three combine into alphaT/beta as:
#
#   alphaT = (equipment_cost / (5*24*365*0.5)) * (1.05**5) * processing_time
#   beta   = tooling_cost
#
# where the denominator amortizes equipment cost over 5 years at 50%
# utilization, and (1.05**5) escalates it over that period.

COST_LEVELS = ["Low", "Low-Medium", "Medium", "Medium-High", "High", "High-Very High", "Very High"]

TOOLING_COST_BY_LEVEL = dict(zip(COST_LEVELS, [1e2, 1e3, 3e3, 1e4, 3e4, 1e5, 1e6]))
EQUIPMENT_COST_BY_LEVEL = dict(zip(COST_LEVELS, [1e3, 1e4, 3e4, 1e5, 3e5, 1e6, 1e7]))
PROCESSING_TIME_BY_LEVEL = dict(zip(COST_LEVELS, [1e-2, 1e-1, 2e0, 3e1, 1e2, 3e2, 1e3]))

EQUIPMENT_UTILIZATION_SEC = 5 * 24 * 365 * 0.5  # 5 yr amortization at 50% utilization
EQUIPMENT_ESCALATION = 1.05 ** 5


def compute_alphaT_beta(tooling_level, equipment_level, time_level):
    """
    Convert Low..Very High tooling cost, equipment cost, and processing
    time levels into (alphaT, beta).


    ########################################################################
    PARAMETERS
    ########################################################################
    ------------------------------------------------------------------------
    tooling_level : str
        One of COST_LEVELS (e.g. 'Medium')

    ------------------------------------------------------------------------
    equipment_level : str
        One of COST_LEVELS

    ------------------------------------------------------------------------
    time_level : str
        One of COST_LEVELS

    ########################################################################
    RETURNS
    ########################################################################
    ------------------------------------------------------------------------
    (alphaT, beta) : tuple of float
    """
    equipment_cost = EQUIPMENT_COST_BY_LEVEL[equipment_level]
    tooling_cost = TOOLING_COST_BY_LEVEL[tooling_level]
    processing_time = PROCESSING_TIME_BY_LEVEL[time_level]

    alphaT = (equipment_cost / EQUIPMENT_UTILIZATION_SEC) * EQUIPMENT_ESCALATION * processing_time
    beta = tooling_cost

    return alphaT, beta


# ----------------------------------------------------------------------
# Function to build a 'Process' object from input data.

def build_process(name, alphaT=None, beta=None, tooling_level=None, equipment_level=None,
                   time_level=None, description="", verbose=True):
    """
    Build a 'Process' object from input data.

    Either provide alphaT and beta directly, or provide tooling_level,
    equipment_level, and time_level (each one of COST_LEVELS) and alphaT/beta
    will be computed via compute_alphaT_beta().


    ########################################################################
    PARAMETERS
    ########################################################################
    ------------------------------------------------------------------------
    name : str
        Process name (e.g., 'CNC')

    ------------------------------------------------------------------------
    alphaT : float, optional
        Time-related cost coefficient. Required unless tooling_level,
        equipment_level, and time_level are provided instead.

    ------------------------------------------------------------------------
    beta : float, optional
        Fixed or scaling cost parameter. Required unless tooling_level,
        equipment_level, and time_level are provided instead.

    ------------------------------------------------------------------------
    tooling_level : str, optional
        One of COST_LEVELS. Used with equipment_level and time_level to
        compute alphaT/beta when they aren't given directly.

    ------------------------------------------------------------------------
    equipment_level : str, optional
        One of COST_LEVELS.

    ------------------------------------------------------------------------
    time_level : str, optional
        One of COST_LEVELS.

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

    # ------------------------------------------------
    # Compute alphaT/beta from categorical levels if
    # they weren't provided directly.
    if alphaT is None or beta is None:
        if tooling_level is None or equipment_level is None or time_level is None:
            print_error(
                f"Error loading proces, '{name}'. Process skipped. Provide either alphaT and beta, "
                f"or tooling_level, equipment_level, and time_level."
            )
            return None

        for level_name, level in [
            ("tooling_level", tooling_level),
            ("equipment_level", equipment_level),
            ("time_level", time_level),
        ]:
            if level not in COST_LEVELS:
                print_error(
                    f"Error loading proces, '{name}'. Process skipped. {level_name} '{level}' "
                    f"must be one of {COST_LEVELS}."
                )
                return None

        alphaT, beta = compute_alphaT_beta(tooling_level, equipment_level, time_level)

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
        description=description,
        tooling_level=tooling_level,
        equipment_level=equipment_level,
        time_level=time_level,
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

        # Extract data. A process file provides either alphaT/beta directly,
        # or tooling_level/equipment_level/time_level for build_process to
        # compute them from.
        name = module.name
        alphaT = getattr(module, "alphaT", None)
        beta = getattr(module, "beta", None)
        tooling_level = getattr(module, "tooling_level", None)
        equipment_level = getattr(module, "equipment_level", None)
        time_level = getattr(module, "time_level", None)

        # Optional description
        description = getattr(module, "description", "")

        # Build process object
        proc = build_process(
            name, alphaT, beta,
            tooling_level=tooling_level,
            equipment_level=equipment_level,
            time_level=time_level,
            description=description,
        )

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