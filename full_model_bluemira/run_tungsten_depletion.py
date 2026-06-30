import random
import os

# Runs the depletion_model.py script 'ns' times, each with a different random number
# seed. That random number seed is then applied in the neutronics_model.py to change
# tokamak armor material composition stochastically. This script was used to
# do the stochastic sampling of tungsten composition for depletion in the WCLL blanket
# for the FS&T paper.

# number of samples to run
ns = 4

for i in range(ns):
  seed = random.randint(1, 100000)
  os.environ["OPENMC_MATERIAL_SEED"] = str(seed)
  print("Running with seed ", os.getenv('OPENMC_MATERIAL_SEED'))
  os.system("python depletion_model.py")
  os.system("python depletion_post.py")

  # the output files are copied into a new folder - "WCLL" will need to be replaced with whatever blanket
  # you are actually running
  os.system('mkdir -p tokamak/WCLL/tungsten_depletion_run/seed_' + str(seed))
  os.system('cp -r tokamak/WCLL/depletion_results/* tokamak/WCLL/tungsten_depletion_run/seed_' + str(seed))
