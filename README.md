HydroBlocks
==========

HydroBlocks relies on a number python libraries. To make this straightforward use conda (http://conda.pydata.org/miniconda.html) and the intel repository to install all the packages. Here are the steps to install the model.


Create a conda environment for HydroBlocks:
```
conda update conda
conda create -n HydroBlocks -c conga-forge python=3.6 anaconda
source activate HydroBlocks
```

Install HydroBlocks dependencies from intel channel:
```
conda install -c conda-forge python=3.6 netcdf4 h5py gdal geos gcc libgcc-ng numpy-base numpy scipy pandas scikit-image scikit-learn mpi4py jpeg kealib xerces-c
```

Install HydroBlocks:
```
git clone https://github.com/chaneyn/HydroBlocks.git
cd HydroBlocks
git checkout dev_noemi
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

