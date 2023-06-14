HydroBlocks
==========

The following steps walk through how to install HydroBlocks

**1. Clone the the HBrouting branch of HydroBlocks.**

```
git clone --single-branch --branch dev_laura https://github.com/chaneyn/HydroBlocks.git
cd HydroBlocks
```

**2. Create a conda environment named HB from the spec-file. Note that the only current spec-file in the repository is for a linux64 machine.** 

```
conda create --name HB --file spec-file.txt
source activate HB
pip install psutil==5.9.4

```

**3. Install the dev_laura branch of HydroBlocks.**

```
python setup.py 
```

