HydroBlocks
==========

HydroBlocks relies on a number python libraries. To make this straightforward use conda (http://conda.pydata.org/miniconda.html). Here are the steps to install the model.


1. Clone and install the dev_nate branch of HydroBlocks
```
git clone --single-branch --branch dev_nate https://github.com/chaneyn/HydroBlocks.git
cd HydroBlocks
python setup.py 
cd ..
```

2. Create a conda environment named HB from the yml file. Note that the only current yml file in the repository is for a linux64 machine. 
```
conda update conda
conda env create -f yml/HB_linux64.yml
source activate HB
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

