HydroBloks
==========

HydroBloks relies on a number of packages. Please make sure you have the following packages installed when running the model. The steps outlined below to install these packages are for a Centos 6 machine. 

1. Numpy (Python)
   pip install numpy
2. Scipy (Python) 
   pip install scipy
3. Scikit-Learn (Python)
   pip install scikit-learn
4. netCDF4 (Python) 
   pip install netCDF4
5. gdal (Python) 
   pip install gdal
6. pygrads (Python) 
   wget http://sourceforge.net/projects/opengrads/files/python-grads/1.1.9/pygrads-1.1.9.tar.gz/download
   tar -xvzf pygrads-1.1.9.tar.gz
   cd pygrads-1.1.9
   python setup.pt install
4. Openblas 
   yum install openblas-devel
5. Netcdf 
   yum install netcdf-devel
6. grads
   yum install grads

Enter into the directory Tools/Tools and insert "python compile.py" to prepare the model tools library

Now we can compile the different model libraries by following these steps:

1. cd Model/pyNoahMP
   python make.py
   cd ..
2. Model/pyDTopmodel
   python Model/pyDTopmodel/src
   python compile.py

*Note - Before compiling dynamic topmodel, you will need to ensure that the intel MKL library is installed.

To run the model on a test dataset, please download and untar the following file within the HydroBloks directory:

http://hydrology.princeton.edu/~nchaney/HydroBloks/Test.tar.gz

To ensure that the data can be prepared for HydroBloks, please add the following line at the end of the esri_extra.wkt file which you can find generally at /usr/local/share/gdal/esri_extra.wkt

102039,PROJCS["USA_Contiguous_Albers_Equal_Area_Conic_USGS_version",GEOGCS["GCS_North_American_1983",DATUM["D_North_American_1983",SPHEROID["GRS_1980",6378137.0,298.257222101]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Albers_Conic_Equal_Area"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",-96.0],PARAMETER["Standard_Parallel_1",29.5],PARAMETER["Standard_Parallel_2",45.5],PARAMETER["Latitude_Of_Origin",23.0],UNIT["Meter",1.0]]


