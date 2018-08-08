HydroBlocks
==========

HydroBlocks relies on a number python libraries. To make this straightforward use conda (http://conda.pydata.org/miniconda.html) and the intel repository (e.g., conda -c intel numpy) to install all the packages. Here are the steps to install the model.

1. conda install gdal
2. conda install scikit-learn
3. conda install netCDF4
4. git clone https://github.com/chaneyn/HydroBlocks.git
5. cd HydroBlocks
6. python setup.py
7. cd ..
8. git clone https://github.com/chaneyn/geospatialtools.git
9. cd geospatialtools
10. python setup.py install
11. cd ..

Note: If you have the intel mkl library on your machine, we recommend that you set the mkl_flag variable in setup.py to True.

To run the model on a test dataset, do the following:

1. wget https://www.dropbox.com/s/tw4z4rf9ol6p24x/HB_sample.tar.gz?dl=0
2. tar -xvzf HB_sample.tar.gz
3. cd HB_sample
4. python ../HydroBlocks/Preprocessing/Driver.py metadata.json
5. python ../HydroBlocks/HydroBlocks/Driver.py metadata.json


