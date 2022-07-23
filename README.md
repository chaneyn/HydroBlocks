HydroBlocks
==========

HydroBlocks relies on a number python libraries. To make this straightforward use conda (http://conda.pydata.org/miniconda.html).

## 1. Clone HydroBlocks:
```
git clone -b dev_noemi https://github.com/chaneyn/HydroBlocks.git
cd HydroBlocks
```

## 2. Install the libray dependencies and compile the model:
Option 1: By creating a conda environment from the yml file:
```
conda env create -f yml/HBenv.yml
source activate HBenv
python setup.py
cd ..
```

Option 2: By installing the dependencies yourself:
```
conda create -n HBenv -y
source activate HBenv
conda install -c conda-forge netcdf4 gdal geos jpeg scikit-learn numpy=1.23 scipy h5py matplotlib cartopy mpi4py zarr opencv gfortran pandas numba
python -m pip install git+https://github.com/chaneyn/geospatialtools.git
python setup.py
cd ..
```

## 3. Run the model on a test dataset:
```
wget https://www.dropbox.com/s/w10u8ocghk3oe21/HB_sample_nv.tar.gz?dl=0
tar -xvzf HB_sample_nv.tar.gz
cd HB_sample
python ../HydroBlocks/Preprocessing/Driver.py metadata.json
python ../HydroBlocks/HydroBlocks/Driver.py metadata.json 
```

## 4. Plot results 
```
python plot.py
```

