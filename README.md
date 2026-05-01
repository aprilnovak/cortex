This repository is to hold all OpenMC models and files generated during
the CORTEX ARPA-E project...

For any large files (like meshes), [put them on Box](https://uofi.app.box.com/folder/303093159824).

## Instructions for building BlueMira with EUDEMO model enabled

```
# delete old environments and do a fresh clone if you get messed up
conda remove -n bluemira --all
rm -rf bluemira
rm -rf process
git clone git@github.com:Fusion-Power-Plant-Framework/bluemira.git
git clone https://github.com/ukaea/PROCESS.git
cd bluemira

# create bluemira conda environment assuming you already have conda installed
sed s/".*python.*"/"  - python="3.11/g ./conda/environment.yml > ./conda/tmp_env.yml
mamba env create -f ./conda/tmp_env.yml -n bluemira

# enter bluemira environment and add some missing dependencies
conda activate bluemira
conda install pythonocc-core
pip install fast_ctd@git+https://github.com/Fusion-Power-Plant-Framework/fast_ctd@main

# install process
cd $HOME
cd PROCESS
pip install setuptools\<65
git checkout 48e3910bed862af1e159e7a46a31b0acacd1676e
cmake -G Ninja -S . -B build -DRELEASE=TRUE
cmake --build build

# check out the bluemira commit with the working model
cd ../bluemira
pip install -e .\[dagmc\]
git checkout 75307027c22364a3ee7710351f050c9d88b475ff
cd eudemo
ln -s $(pwd)/../data/materials config/materials_data
```
## Instructions for Installing MOAB

Provide CMake arguments during MOAB installation to enable HDF5 format and be able to load exodus files.

```
CMAKE_ARGS="-DENABLE_HDF5=ON -DENABLE_NETCDF=ON" python -m pip install .
```

## Instructions for Improv

Put these in your bashrc before you begin anything (change for XS location as appropriate).

```
export OPENMC_CROSS_SECTIONS=$HOME_DIRECTORY_SYM_LINK/cross_sections/endfb-vii.1-hdf5/cross_sections.xml
export ENABLE_DAGMC=true
export PATH=$PATH:$DIRECTORY_WHERE_YOU_HAVE_CARDINAL/cardinal/build/openmc/bin
```

Example job script:

```
#!/bin/bash -l

# Usage:
# 1. Copy to the directory where you have your files
# 2. Update any needed environment variables and input file names in this script
# 3. qsub job_improv

#PBS -A Radiant
#PBS -l select=20:ncpus=128:mpiprocs=128
#PBS -l walltime=01:00:00
#PBS -q compute
#PBS -j oe
#PBS -N cardinal

#PBS -m bea

module purge
module load gcc/11.4.0
module load openmpi/4.1.6-gcc-11.4.0-pbs
module load cmake/3.27.4
module load perl/5.38.0-gcc-11.4.0
module load miniforge3/25.3.0

# Revise for your Cardinal repository location
DIRECTORY_WHERE_YOU_HAVE_CARDINAL=$HOME

# This is needed because your home directory on Improv is actually a symlink
HOME_DIRECTORY_SYM_LINK=$(realpath -P $DIRECTORY_WHERE_YOU_HAVE_CARDINAL)
export NEKRS_HOME=$HOME_DIRECTORY_SYM_LINK/cardinal/install

# Revise for your cross sections location
export OPENMC_CROSS_SECTIONS=$HOME_DIRECTORY_SYM_LINK/cross_sections/endfb-vii.1-hdf5/cross_sections.xml

export CARDINAL_DIR=$HOME_DIRECTORY_SYM_LINK/cardinal

# The name of the input file you want to run
input_file=nek.i

# Moving into the working directory (where the job script was launched).
echo "Working directory: $PBS_O_WORKDIR"
cd $PBS_O_WORKDIR

# Run a Cardinal case
mpirun -np 20 openmc -s 128 neutronics_model.xml > logfile
```

To build Cardinal with the modules shown in the example job script above, add these three lines
in `cardinal/config/moab.mk`:

```
       -DENABLE_FORTRAN=OFF \
       -DENABLE_EIGEN3=ON \
       -DENABLE_PYMOAB=ON \
```

Compile Cardinal following the without-conda instructions (here)[https://cardinal.cels.anl.gov/without_conda.html]. After you have obtained the cardinal executable, the last step is to install pymoab.

```
cd cardinal/build/moab
pip install .
```

## TOKAMAK - SLAB INSTRUCTIONS

In order to run the slab model, you must first run the Tokamak model to generate two important files: ```surface_source.h5``` and ```armor_current_neutron.json```. 

To generate these files, run the Tokamak model: 
1. Run ```python neutronics_model.py``` in the respective directory ```./detailed_bluemira(_hcll, _wcll, _hcpb)```.
2. Execute ```openmc``` command inside ```./neutronics_run/``` folder.
3. Run the ```python neutronics_post.py```.

Ideally we will have one example ```surface_source.h5``` and the ```armor_current_neutron.json``` provided to be used in all slab cases for a breeder type. This would remove the necessity of running the detailed_model ahead of the slab. However, we need to test that this source don't change drastically with different materials in the same breeder type.
