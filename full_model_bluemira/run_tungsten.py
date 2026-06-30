import random
import os

# Runs the neutronics_model.py script 'ns' times, each with a different random number
# seed. That random number seed is then applied in the inputs.py to change
# eurofer structural material composition stochastically. This script was used to
# do the stochastic sampling of eurofer composition for transport (not depletion) in the WCLL blanket
# for the FS&T paper. This script will launch jobs on Improv, wait until results are obtained,
# then keep launching jobs sequentially.

# number of samples to run
ns = 4

for i in range(ns):
  seed = random.randint(1, 100000)
  os.environ["OPENMC_MATERIAL_SEED"] = str(seed)
  print("Running with seed ", os.getenv('OPENMC_MATERIAL_SEED'))
  os.system("python neutronics_model.py")

  # the output files are copied into a new folder - "WCLL" will need to be replaced with whatever blanket
  # you are actually running
  os.chdir('tokamak/WCLL/neutronics_run')
  os.system('qsub job_improv')

  # check if statepoint has been created with n > 10 every 2 minutes
  import glob
  import time
  wait_time = 120
  while True:
    time.sleep(wait_time)
    files = glob.glob("statepoint*")
    print(files)
    if (len(files)) > 1:
      print('Finished with statepoint ', files[0])
      break

  os.chdir('../../../')
  os.system('python neutronics_post.py')

  # the output files are copied into a new folder - "WCLL" will need to be replaced with whatever blanket
  # you are actually running
  os.system('mkdir -p tokamak/WCLL/tungsten_run/seed_' + str(seed))
  os.system('cp -r tokamak/WCLL/neutronics_results/* tokamak/WCLL/tungsten_run/seed_' + str(seed))
