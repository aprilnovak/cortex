#!/usr/bin/env python3
"""
inputs.py
=========
Central configuration file for all tokamak / slab neutronics simulations.

Edit this file to switch between:
  - Simulation type  : "tokamak", "slab" , or "fast_slab"
  - Breeder type     : "WCLL", "HCLL", or "HCPB"
  - Environment      : "CORTEX" or "TESTING"
  - Run options      : triggers, surface-source write/read, albedo, etc.
"""

from pathlib import Path
import os
import sys
import shutil

# =============================================================================
# SIMULATION TYPE
# =============================================================================
SIM_TYPE     = "tokamak"      # "tokamak"  |  "slab" | "fast_slab"

# =============================================================================
# BREEDER TYPE
# =============================================================================
BREEDER_TYPE = "HCPB"         # "WCLL"  |  "HCLL"  |  "HCPB"

# =============================================================================
# ENVIRONMENT FLAG
# =============================================================================
# FOR SLAB MODEL we need two files generated from tokamak (surface_source and current_armor)
# TESTING -> will load from previous tokamak run/folder
# CORTEX  -> will load from file saved in the same slab folder
#            (if not available, it will copy from tokamak run/folder)
ENVIRONMENT  = "TESTING"      # "CORTEX"  |  "TESTING"

# =============================================================================
# VALIDATION
# =============================================================================
if SIM_TYPE not in ("tokamak", "slab", "fast_slab"):
    sys.exit(f"[inputs] Unknown SIM_TYPE '{SIM_TYPE}'. "
             f"Choose from: tokamak, slab")

if BREEDER_TYPE not in ("WCLL", "HCLL", "HCPB"):
    sys.exit(f"[inputs] Unknown BREEDER_TYPE '{BREEDER_TYPE}'. "
             f"Choose from: WCLL, HCLL, HCPB")

if ENVIRONMENT not in ("CORTEX", "TESTING"):
    sys.exit(f"[inputs] Unknown ENVIRONMENT '{ENVIRONMENT}'. "
             f"Choose from: CORTEX, TESTING")

# =============================================================================
# DIRECTORY LAYOUT
# =============================================================================
BASE_DIR = Path.cwd()
PROJECT_ROOT = BASE_DIR.parent             

# Simulation-type subdirectory (all run/result files go here)
SIM_DIR = BASE_DIR / SIM_TYPE / BREEDER_TYPE
# e.g. BASE_DIR/tokamak/WCLL/  or  BASE_DIR/slab/HCPB/

# Run / result directories — all nested under SIM_DIR 
NEUTRONICS_RUN_DIR     = SIM_DIR / "neutronics_run"
NEUTRONICS_RESULTS_DIR = SIM_DIR / "neutronics_results"
DEPLETION_RUN_DIR      = SIM_DIR / "depletion_run"
DEPLETION_RESULTS_DIR  = SIM_DIR / "depletion_results"
SUMMARY_RESULTS_CSV = NEUTRONICS_RESULTS_DIR / "OB_1_b6" / "summary_results.csv"

# Important input files — all live in PROJECT_ROOT/dagmc_files/ 
DAGMC_DIR = PROJECT_ROOT / "dagmc_files"

_JSON_FILENAMES = {
    "WCLL": "EUDEMO_WCLL_inputs.json",
    "HCLL": "EUDEMO_HCLL_inputs.json",
    "HCPB": "EUDEMO_HCPB_inputs.json",
}

_DAGMC_FILENAMES = {
    "WCLL": "eudemo_wcll.h5m",
    "HCLL": "eudemo_hcll.h5m",
    "HCPB": "eudemo_hcpb.h5m",
}

_FAST_SLAB_DAGMC_FILENAMES = {
    "WCLL": "eudemo_wcll_fast_slab.h5m",
    "HCLL": "eudemo_hcll_fast_slab.h5m",
    "HCPB": "eudemo_hcpb_fast_slab.h5m",
}

if SIM_TYPE == "fast_slab":
    DAGMC_MODEL_FILE = DAGMC_DIR / _FAST_SLAB_DAGMC_FILENAMES[BREEDER_TYPE]
else:
    DAGMC_MODEL_FILE = DAGMC_DIR / _DAGMC_FILENAMES[BREEDER_TYPE]

INPUT_JSON       = DAGMC_DIR / _JSON_FILENAMES[BREEDER_TYPE]
EQUILIBRIUM_FILE = DAGMC_DIR  / "equilibrium_eqdsk.json"
MATERIALS_MODULE_DIR = PROJECT_ROOT / "materials"  

# read from bashrc as an environment variable
#CHAIN_FILE           = BASE_DIR / "depletion_chain" / "chain_endfb80_sfr.xml"

# Depletion chain — path set in ~/.bashrc as $CHAIN_FILE
_chain_env = os.environ.get("OPENMC_CHAIN_FILE")
if not _chain_env:
    sys.exit("[inputs] Environment variable OPENMC_CHAIN_FILE is not set. "
             "Add it to your ~/.bashrc: export OPENMC_CHAIN_FILE=/path/to/chain.xml")
OPENMC_CHAIN_FILE = Path(_chain_env)
if not OPENMC_CHAIN_FILE.is_file():
    sys.exit(f"[inputs] CHAIN_FILE does not exist: {OPENMC_CHAIN_FILE}")

# =============================================================================
# TALLY SETTINGS
# =============================================================================
TRIGGER_CHUNK_KEY           = "OB_1_b6"
ALBEDO_CHUNK_KEY            = "OB_1_b6"
SLAB_CHUNK_KEY  = "OB_1_b6"   # DO NOT CHANGE — fixed by DAGMC slab geometry

# =============================================================================
# SLAB STATIC FILES
# =============================================================================
if ENVIRONMENT == "CORTEX":
    # Static files, one per breeder type, fixed location
    # Pre-copied once into slab/BREEDER_TYPE/ and never change between runs.
    SURFACE_SOURCE_FILE = SIM_DIR / "surface_source.h5"
    ARMOR_CURRENT_JSON  = SIM_DIR / "armor_current_neutron.json"

else:  # TESTING
    # Load directly from the tokamak simulation outputs
    SURFACE_SOURCE_FILE = (
        BASE_DIR / "tokamak" / BREEDER_TYPE / "neutronics_run" / "surface_source.h5"
    )
    ARMOR_CURRENT_JSON  = (
        BASE_DIR / "tokamak" / BREEDER_TYPE / "neutronics_results"
        / ALBEDO_CHUNK_KEY / "albedo" /"armor_current_neutron.json"
    )

# CORTEX + slab only: one-time copy from tokamak outputs if files missing
if ENVIRONMENT == "CORTEX" and SIM_TYPE in ("slab", "fast_slab"):
    _tokamak_surface_source = (
        BASE_DIR / "tokamak" / BREEDER_TYPE / "neutronics_run" / "surface_source.h5"
    )
    _tokamak_armor_current = (
        BASE_DIR / "tokamak" / BREEDER_TYPE / "neutronics_results"
        / ALBEDO_CHUNK_KEY / "albedo" / "armor_current_neutron.json"
    )
    _file_map = {
        SURFACE_SOURCE_FILE: _tokamak_surface_source,
        ARMOR_CURRENT_JSON:  _tokamak_armor_current,
    }
    missing = [dst for dst in _file_map if not dst.exists()]
    if missing:
        print(f"[inputs] CORTEX: slab static files missing for {BREEDER_TYPE}, "
              f"attempting one-time copy from tokamak results ...")
        SIM_DIR.mkdir(parents=True, exist_ok=True)
        for dst, src in _file_map.items():
            if dst.exists():
                continue
            if not src.exists():
                print(f"[inputs] WARNING: tokamak source not found, cannot copy:\n"
                      f"         src : {src}\n"
                      f"         dst : {dst}")
                continue
            shutil.copy2(src, dst)
            print(f"[inputs] Copied {src.name}:\n"
                  f"         src : {src}\n"
                  f"         dst : {dst}")

# =============================================================================
# MATERIALS 
# (FOR CORTEX, UPDATE HERE WITH USER'S MATERIALS)
# =============================================================================
if str(MATERIALS_MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MATERIALS_MODULE_DIR))
import materials

print("materials file:", materials.__file__)
print("SIM_TYPE:", SIM_TYPE)
print("BREEDER_TYPE:", BREEDER_TYPE)

# Materials common to all breeder types 
armor_material      = materials.W(19.3)
structural_material = materials.eurofer97(7.87)
vv_material = materials.ss316Ln_ig(7.93)

# Breeder-specific materials 
if BREEDER_TYPE == "WCLL":
    coolant_material    = materials.Water(0.866)       # pressurised water
    breeder_material    = materials.PbLi(0.90, 9.8)    # 90% Li6 enriched
    multiplier_material = None

elif BREEDER_TYPE == "HCLL":
    coolant_material    = materials.Helium(0.00498)    # ~8 MPa, 500 °C
    breeder_material    = materials.PbLi(0.90, 9.8)    # 90% Li6 enriched
    multiplier_material = None

elif BREEDER_TYPE == "HCPB":
    coolant_material    = materials.Helium(0.00498)     # ~8 MPa, 500 °C
    breeder_material    = materials.kalos_cb(2.52)      # Li4SiO4 pebbles
    multiplier_material = materials.be12ti(2.25)        # Be12Ti neutron multiplier

#print("vv_material:", vv_material)
#print("structural_material:", structural_material)
#print("coolant_material:", coolant_material)
#print("breeder_material:", breeder_material)

# =============================================================================
# TOKAMAK GEOMETRY
# =============================================================================
NUMBER_OF_SECTORS = 16
THETA0_DEG        = 0.0

# =============================================================================
# SLAB GEOMETRY
# =============================================================================
# SLAB GEOMETRY
if SIM_TYPE == "fast_slab":      # (DO NOT CHANGE)
    SLAB_TARGET_VOL_ID = 13      # VV in equatorial chunk in stripped dagmc
else:
    SLAB_TARGET_VOL_ID = 78      # VV in full dagmc Model 

SLAB_R_IN_CM       = 1130.0  # DO NOT CHANGE 
SLAB_R_OUT_CM      = 1380.0  # DO NOT CHANGE 

# =============================================================================
# REACTOR POWER & SOURCE SCALING
# =============================================================================
TOTAL_FUSION_POWER_W = 2.0e9
EV_PER_FUSION        = 17.6e6

# =============================================================================
# CHUNK / RADIAL BIN SETTINGS
# =============================================================================
CHUNK_GAP_CM   = 2.0
CHUNK_START_CM = 0.0

# =============================================================================
# OPENMC RUN SETTINGS
# =============================================================================
if SIM_TYPE == "tokamak":
    TALLY_CONVERGENCE_THRESHOLD = 0.1 # 10%
    BATCHES                     = 15
    PARTICLES_PER_BATCH    = 1_000_000   
else:
    TALLY_CONVERGENCE_THRESHOLD = 0.01 # 1%
    BATCHES                     = 10
    PARTICLES_PER_BATCH    = 100_000 


USE_TRIGGER            = True   # False = fixed batches, no convergence check
TRIGGER_BATCH_INTERVAL = 5
TRIGGER_MAX_BATCHES    = 2000

# if USE_TRIGGER = False
FIXED_BATCHES          = 10    # only used when USE_TRIGGER = False

# Depletion
OPENMP_THREADS         = None    # None for max available
DEPLETION_PROCESSES    = None    # None for max available

# =============================================================================
# FEATURE TOGGLES
# =============================================================================
DO_PHOTON_TRANSPORT          = True
DO_WRITE_SURFACE_SOURCE      = True # TBD: Add if statement to "False" if CORTEX?
SURFACE_SOURCE_SURFACE_ID    = 287 # DO NOT CHANGE 
SURFACE_SOURCE_CELL_TO       = 66  # DO NOT CHANGE 
SURFACE_SOURCE_MAX_PARTICLES = 1_000_000


if SIM_TYPE == "tokamak":
    DO_ALBEDO       = True
    DO_CHECK_SOURCE = True
    #TRACKS = True
else:
    DO_ALBEDO       = False
    DO_CHECK_SOURCE = False
    #TRACKS = False

# Depletion
IRRADIATION_YEARS    = 5.0
CONSTANT_POWER_RATIO = 0.3
REDUCED_CHAIN_LEVEL = 5

# =============================================================================
# POST-PROCESSING CHUNK SELECTION
# =============================================================================
NEUTRONICS_KEYS_TO_PROCESS = [
    "OB_1_b6", # DO NOT REMOVE THIS ONE ! TBD: add safe guards for slab tokamak/slab interaction
]

if "OB_1_b6" not in NEUTRONICS_KEYS_TO_PROCESS:
    sys.exit(
        "[inputs] NEUTRONICS_KEYS_TO_PROCESS must include 'OB_1_b6' — "
    )

# Initial desired layers for CORTEX SUMMARY
# Can be added others, such as "Breeder layer 1"
SUMMARY_TARGET_LAYERS = [
    "Armor",
    "First_Wall",
    "Vacuum Vessel",
]

POLOIDAL_LAYERS = [
    ("Armor",      0),
    ("First_Wall", 1),
    ("VV",        -1),
]

DEPLETION_CHUNKS = {
    "OB_1_b6": None,
}

DEPLETION_IDX_TO_PLOT = (0, 3, 6, 9, 12, 15, 18, 21, 24, 27)

# =============================================================================
# REFERENCE (GOLD) RUN COMPARISON
# =============================================================================
REFERENCE_LABEL       = "W_armor_eurofer_structural"   # None = disabled
REFERENCE_PLOT_SUFFIX = "_vs_ref"                      # appended to overlay figs


# Improvements ideas
# A) Introduce Random RAY
# B) Fine tune number of particles
# C) Clear code to Cortex to improve performance:
#   1 - Remove albedo tallies and calculations.
#   2 - Merge neutronics/post/depletion/post in same script.
#       2.1 - Avoid re-loading neutronics_model.py
#       2.2 - Remove storing/loading info mapping across files (tallies, material, cell, volume)
#       2.3 - Add safe guard for ONLY OB_1_b6 for slab model. (add warning that is ignoring others)
#
# LIMITATIONS:
#   Current tokamak model allows any breeder chunk tally. However Slab model only works
#   for OB_1_b6. To improve for any given breeder chunk we must first improve:
#      1 - _get_breeder_reflective_cuts
#      2 - TBD