HydroBlocks
==========

The following steps walk through how to install HydroBlocks

**1. Clone the HydroBlocks repository.**

```
git clone https://github.com/chaneyn/HydroBlocks.git
cd HydroBlocks
```

**2. Create a conda environment named HB from the spec-file. Note that the only current spec-file in the repository is for a linux64 machine.** 

```
conda create --name HB --file spec-file.txt
source activate HB
pip install git+https://github.com/chaneyn/geospatialtools@dev_nate
pip install psutil==5.9.4
```

**3. Install HydroBlocks.**

```
python setup.py 
```

