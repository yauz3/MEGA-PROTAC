<h1> MEGA PROTAC: PROTAC-Mediated Ternary Complex Formation Pipeline based on MEGADOCK with Sequential Filtering integrated with Rank Aggregation. </h1>


**Pre-installation:**


*1*- MEGA DOCK Installation:
Before proceeding with the installation, ensure that MEGA DOCK is installed on your system. It is highly recommended to follow the original documentation provided by the developers for detailed instructions:

Please refer to their documentation to install MEGA DOCK correctly before continuing with this project.

[MEGA DOCK GitHub Repository
](https://github.com/akiyamalab/MEGADOCK)

Note: Before running the project, ensure that the library PATH located on "../bin/MEGADOCK/Makefile" and compiler settings are properly configured on your system.


*2*- DockQ Installation:

To install DockQ, please follow the protocol outlined in their official documentation:

[DockQ GitHub Repository
](https://github.com/bjornwallner/DockQ)

*3*- PIZSA Installation:
To intall PIZSA, please follow the protocol outlined in their official documentation:

[PIZSA Website
](http://cospi.iiserpune.ac.in/pizsa/Download/Download.html)

The three program files are located in the *BIN* directory.

<h2>How to set up the environment</h2>

*4* Make sure FCC and Vorquma paths are correct on your system
FCC_path="../fcc/scripts"
voronota_path="../voronota_1.22.3149"


**1- Conda Installation**:
-We have provided an Anaconda environment file for easy set-up. If you do not have Anaconda installed, you can get Miniconda from [HERE](https://docs.anaconda.com/free/miniconda/).

Then, install ```bash mamba```:
```bash
conda install -c conda-forge mamba
```

If an individual encounters difficulties in implementing mamba in their primary setting, proceed to establish a fresh environment for mamba by following these steps:
```bash
conda create -n mamba-env -c conda-forge mamba
conda activate mamba-env
```

```bash
conda create -n mamba-env -c conda-forge mamba
conda activate mamba-env
```




bunlar yüklenip değiştirilmesi lazım.

