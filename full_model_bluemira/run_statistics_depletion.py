import random
import os

# Runs the depletion_model.py 'ns' times, each with a different random number seed. This
# was used for the FS&T paper to justify that the propagated uncertainty through depletion
# is sufficiently small.

# number of samples to run
ns = 10

for i in range(ns):
  seed = random.randint(1, 100000)
  os.environ["OPENMC_RNG_SEED"] = str(seed)
  print("Running with seed ", os.getenv('OPENMC_RNG_SEED'))
  os.system("python depletion_model.py")
  os.system("python depletion_post.py")

  # the output files are copied into a new folder - "WCLL" will need to be replaced with whatever blanket
  # you are actually running
  os.system('mkdir -p tokamak/WCLL/statistics_run/seed_' + str(seed))
  os.system('cp -r tokamak/WCLL/depletion_results/* tokamak/WCLL/statistics_run/seed_' + str(seed))
