HydroBlocks
==========

The following steps walk through how to install HydroBlocks and run it over the SGP site in Oklahoma

**1. Clone the the HBrouting branch of HydroBlocks.**

```
git clone --single-branch --branch HBrouting https://github.com/chaneyn/HydroBlocks.git
cd HydroBlocks
```

**2. Create a conda environment named HB from the yml file. Note that the only current yml file in the repository is for a linux64 machine.** 

```
conda install conda==4.8.5
conda env create -f yml/HB_linux64.yml
source activate HB
```

**3. Install the HBrouting branch of HydroBlocks.**

```
python setup.py 
cd ..
```

**4. Download the SGP site data and run the model.**

```
wget http://hydrology.cee.duke.edu/HydroBlocks/SGP_OK_1deg.tar.gz
tar -xvzf SGP_OK_1deg.tar.gz
vi SGP_OK_1deg/experiments/json/baseline.json
Set the variable rdir to the absolute path of SGP_OK_1deg
mpirun -n 16 ./HydroBlocks/HB -m SGP_OK_1deg/experiments/json/baseline.json -t preprocess
mpirun -n 16 ./HydroBlocks/HB -m SGP_OK_1deg/experiments/json/baseline.json -t model
```

