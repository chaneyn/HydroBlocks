HydroBlocks
==========

HydroBlocks relies on a number python libraries. To make this straightforward use conda (http://conda.pydata.org/miniconda.html) and the intel repository to install all the packages. Here are the steps to install the model.


# 1. Install the libray dependencies:
Option 1: By creating a conda environment from yml file:
```
conda env create -f yml/HBenv.yml
source activate HBenv
```

Option 2: By install the dependencies yourself:
```
conda create -n HBenv -y
source activate HBenv
conda install -c conda-forge netcdf4 gdal geos jpeg scikit-learn numpy scipy h5py matplotlib cartopy mpi4py zarr opencv gfortran pandas rasterio xarray
python -m pip install git+https://github.com/chaneyn/geospatialtools.git
```

# 2. Install HydroBlocks Model:
```
git clone -b dev_noemi https://github.com/chaneyn/HydroBlocks.git
cd HydroBlocks
python setup.py 
cd ..
```

# 3.Run the model on a test dataset:
```
wget https://www.dropbox.com/s/k7su7af5dk1l2vf/HB_sample.tar.gz?dl=0
tar -xvzf HB_sample.tar.gz
cd HB_sample
python ../HydroBlocks/Preprocessing/Driver.py metadata.json
python ../HydroBlocks/HydroBlocks/Driver.py metadata.json 
```

# 4. Plot results 
```
python plot.py
```

