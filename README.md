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
