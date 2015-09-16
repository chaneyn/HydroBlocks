HydroBloks
==========

HydroBloks relies on a number python libraries. To make this straightforward use conda (http://conda.pydata.org/miniconda.html) to install all the packages. Here are the steps to install the model.

1. conda install gdal
2. conda install scikit-learn
3. conda install netCDF4
4. python setup.py

Note: If you have the intel mkl library on your machine, we recommend that you set the mkl_flag variable in setup.py to True.

To run the model on a test dataset, do the following:

1. wget http://hydrology.princeton.edu/~nchaney/HydroBloks/LittleWashita.tar.gz
2. untar -xvzf LittleWashita.tar.gz
3. cd LittleWashita
4. python ../HydroBloks/Preprocessing/Driver.py metadata_littlewashita.json
5. python ../HydroBloks/HydroBloks/Driver.py metadata_littlewashita.json


