import random
import os

for i in range(15):
  seed = random.randint(1, 100000)
  os.environ["OPENMC_MATERIAL_SEED"] = str(seed)
  print("Running with seed ", os.getenv('OPENMC_MATERIAL_SEED'))
  os.system("python neutronics_model.py")

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
  os.system('mkdir -p tokamak/WCLL/tungsten_run/seed_' + str(seed))
  os.system('cp -r tokamak/WCLL/neutronics_results/* tokamak/WCLL/tungsten_run/seed_' + str(seed))
