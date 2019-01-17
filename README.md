HydroBlocks
==========

HydroBlocks relies on a number python libraries. To make this straightforward use conda (http://conda.pydata.org/miniconda.html) and the intel repository to install all the packages. Here are the steps to install the model.


Create a conda environment for HydroBlocks:
```
conda update conda
conda create -n HydroBlocks -c intel python=3.6 anaconda
source activate HydroBlocks
```

Install HydroBlocks dependencies from intel channel:
```
conda install -c intel gdal netcdf4 geos xerces-c jpeg scikit-image scikit-learn numpy pandas h5py kealib gcc libgcc python=3.6 mpi4py 
conda install conda-forge::gdal=2.2.4 --force --update-dependencies
```

Install HydroBlocks:
```
git clone https://github.com/chaneyn/HydroBlocks.git
cd HydroBlocks
python setup.py 
cd ..
```

Install geospatial tools:
```
git clone https://github.com/chaneyn/geospatialtools.git
cd geospatialtools
python setup.py install
cd ..
```

To run the model on a test dataset:
```
wget https://www.dropbox.com/s/k7su7af5dk1l2vf/HB_sample.tar.gz?dl=0
tar -xvzf HB_sample.tar.gz
cd HB_sample
python ../HydroBlocks/Preprocessing/Driver.py metadata.json
python ../HydroBlocks/HydroBlocks/Driver.py metadata.json 
```

```
source deactivate 
```

