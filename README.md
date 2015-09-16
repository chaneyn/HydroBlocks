HydroBloks
==========

HydroBloks relies on a number python libraries. To make this straightforward use conda (http://conda.pydata.org/miniconda.html) to install all the packages. Here are the steps to install the model.

1. conda install gdal
2. conda install scikit-learn
3. conda install netCDF4
4. Compile HydroBloks by entering "python setup.py" at the command line.

Note: If you have the intel mkl library on your machine, we recommend that you set the mkl_flag variable in setup.py to True.

To run the model on a test dataset, follow the following steps

wget http://hydrology.princeton.edu/~nchaney/HydroBloks/LittleWashita.tar.gz
untar -xvzf LittleWashita.tar.gz
cd LittleWashita
python ../HydroBloks/Preprocessing/Driver.py metadata_littlewashita.json
python ../HydroBloks/HydroBloks/Driver.py metadata_littlewashita.json


