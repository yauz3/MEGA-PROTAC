**Pre-installation:**

First, create a **./bin** dictionary and put [MEGADOCK](https://github.com/akiyamalab/MEGADOCK), [DockQ](https://github.com/bjornwallner/DockQ) and [PIZSA](http://cospi.iiserpune.ac.in/pizsa/Download/Download.html) into **./bin** dictionary.

*1*- MEGA DOCK Installation:
Before proceeding with the installation, ensure that MEGA DOCK is installed on your system. It is highly recommended to follow the original documentation provided by the developers for detailed instructions:
@@ -13,42 +14,37 @@ Please refer to their documentation to install MEGA DOCK correctly before contin
](https://github.com/akiyamalab/MEGADOCK)

Note: Before running the project, ensure that the library PATH is located on "./bin/MEGADOCK/Makefile" and that compiler settings are properly configured on your system.
Put the MEGA DOCK file into the *bin* dictionary.
Put the MEGA DOCK file into the **./bin** dictionary.

HINT: megadock and megadock-gpu should be in the folder. Otherwise, the MEGA PROTAC won't work because of failure installation for MEGADOCK.

*2*- DockQ [1] Installation:

*2*- DockQ [1] Installation: 

To install DockQ, please follow the protocol outlined in their official documentation:

[DockQ GitHub Repository
](https://github.com/bjornwallner/DockQ)
Put the DOCKQ file into the *bin* dictionary.
Put the DOCKQ file into the **./bin** dictionary.

HINT: DockQ has been used in the validation step. MEGA PROTAC can be used without DockQ, but it is necessary for validation on test tests.

*3*- PIZSA [2] Installation:
To install PIZSA, please follow the protocol outlined in their official documentation:

[PIZSA Website
](http://cospi.iiserpune.ac.in/pizsa/Download/Download.html)

Put the PIZSA file into the *bin* dictionary.


*4*- Install [Voromqa](https://github.com/kliment-olechnovic/voronota), [FCC](https://github.com/haddocking/fcc).
Make sure FCC [3] and VoroMQA [4] paths are correct on the files:

"../uti/clustering_fast.py"

"../uti/formation_searh_space_small.py"

"../uti/ranking_of_cluster_faster.py"
Put the PIZSA file into the **./bin** dictionary.

The paths should be arranged before a run:
HINT: **./uti/ranking_PIZSA_cluster_rank.py** should be worked in the PIZSA environment.
Be careful about the **err_dir** path in **./bin/run_PIZSA.py**

FCC_path="../fcc/scripts"
*4*- Install [Voromqa](https://github.com/kliment-olechnovic/voronota); and put the Voromqa [4] file into **./bin**

voronota_path="../voronota_1.22.3149"
HINT: The "voromqa" binary should be after the successful installation of Voromqa.

fcc="../fcc"
*5*- Install [FCC](https://github.com/haddocking/fcc), and put FCC [3] into **./bin**


<h2>How to set up the environment</h2>
@@ -71,32 +67,28 @@ conda activate mamba-env
**2- Conda Env Formation**:
- To establish a conda environment, please execute the following code:
```bash
<<<<<<< HEAD
mamba env create -f environment.yml
=======
mamba env create -n mega_protact -f environment.yml
>>>>>>> c5778318 (read_me)
```

- PIZSA requires Python2, which may lead to conflicts. The second environment has been specifically developed to mitigate that issue.

```bash
mamba env create -n pizsa -f environment_pizsa.yml
mamba env create -f environment_pizsa.yml
```

**3- Activate Conda Env**:
- Before running MEF-AlloSite, enable the conda environment.
```bash
mamba activate mega_protact
mamba activate mega_protac
```
It should be noted that the terminal may need to be closed and re-opened to activate the mega_protact environment.
It should be noted that the terminal may need to be closed and re-opened to activate the mega_protac environment.

Next, proceed to install supplementary packages using the pip command as follows:

```bash
pip install -r requirements.txt
```
NOTE: Please make sure installation into the *mega_protact* environment. 
NOTE: Please make sure installation into the *mega_protac* environment. 

https://github.com/thelahunginjeet/pyrankagg/tree/master

@@ -210,28 +202,6 @@ python3 python3 STEP_14_PROTAC_docking.py


**References**

[1] Basu, Sankar, and Björn Wallner. "DockQ: a quality measure for protein-protein docking models." PloS one 11.8 (2016): e0161879.

[2] Roy, Ankit A., et al. "Protein Interaction Z Score Assessment (PIZSA): an empirical scoring scheme for evaluation of protein–protein interactions." Nucleic acids research 47.W1 (2019): W331-W337.

[3] Rodrigues, João PGLM, et al. "Clustering biomolecular complexes by residue contacts similarity." Proteins: Structure, Function, and Bioinformatics 80.7 (2012): 1810-1817.

[4] Olechnovič, Kliment, and Česlovas Venclovas. "VoroMQA: Assessment of protein structure quality using interatomic contact areas." Proteins: Structure, Function, and Bioinformatics 85.6 (2017): 1131-1145.

[5] Michaud‐Agrawal, Naveen, et al. "MDAnalysis: a toolkit for the analysis of molecular dynamics simulations." Journal of computational chemistry 32.10 (2011): 2319-2327.

[6] Mitternacht, Simon. "FreeSASA: An open source C library for solvent accessible surface area calculations." F1000Research 5 (2016).

[7] O'Boyle, Noel M., et al. "Open Babel: An open chemical toolbox." Journal of cheminformatics 3 (2011): 1-14.

[8] Weng, Gaoqi, et al. "Integrative modeling of PROTAC-mediated ternary complexes." Journal of Medicinal Chemistry 64.21 (2021): 16271-16281.

[9] Rao, Arjun, et al. "Bayesian optimization for ternary complex prediction (BOTCP)." Artificial Intelligence in the Life Sciences 3 (2023): 100072.
