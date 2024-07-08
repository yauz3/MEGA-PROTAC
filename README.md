<h1> MEGA PROTAC: PROTAC-Mediated Ternary Complex Formation Pipeline based on MEGADOCK with Sequential Filtering integrated with Rank Aggregation. </h1>


**Pre-installation:**

First, create a **./bin** dictionary and put [MEGADOCK](https://github.com/akiyamalab/MEGADOCK), [DockQ](https://github.com/bjornwallner/DockQ) and [PIZSA](http://cospi.iiserpune.ac.in/pizsa/Download/Download.html) into **./bin** dictionary.

*1*- MEGA DOCK Installation:
Before proceeding with the installation, ensure that MEGA DOCK is installed on your system. It is highly recommended to follow the original documentation provided by the developers for detailed instructions:

Please refer to their documentation to install MEGA DOCK correctly before continuing with this project.

[MEGA DOCK GitHub Repository
](https://github.com/akiyamalab/MEGADOCK)

Note: Before running the project, ensure that the library PATH is located on "./bin/MEGADOCK/Makefile" and that compiler settings are properly configured on your system.
Put the MEGA DOCK file into the **./bin** dictionary.


*2*- DockQ [1] Installation:

To install DockQ, please follow the protocol outlined in their official documentation:

[DockQ GitHub Repository
](https://github.com/bjornwallner/DockQ)
Put the DOCKQ file into the **./bin** dictionary.

*3*- PIZSA [2] Installation:
To install PIZSA, please follow the protocol outlined in their official documentation:

[PIZSA Website
](http://cospi.iiserpune.ac.in/pizsa/Download/Download.html)

Put the PIZSA file into the **./bin** dictionary.


*4*- Install [Voromqa](https://github.com/kliment-olechnovic/voronota), [FCC](https://github.com/haddocking/fcc).
Make sure FCC [3] and VoroMQA [4] paths are correct on the files:

"./uti/clustering_fast.py"

"./uti/formation_searh_space_small.py"

"./uti/ranking_of_cluster_faster.py"

The paths should be arranged before a run:

FCC_path="../fcc/scripts"

voronota_path="../voronota_1.22.3149"

fcc="../fcc"


<h2>How to set up the environment</h2>


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

**2- Conda Env Formation**:
- To establish a conda environment, please execute the following code:
```bash
mamba env create -f environment.yml
```

- PIZSA requires Python2, which may lead to conflicts. The second environment has been specifically developed to mitigate that issue.

```bash
mamba env create -n pizsa -f environment_pizsa.yml
```

**3- Activate Conda Env**:
- Before running MEF-AlloSite, enable the conda environment.
```bash
mamba activate mega_protact
```
It should be noted that the terminal may need to be closed and re-opened to activate the mega_protact environment.

Next, proceed to install supplementary packages using the pip command as follows:

```bash
pip install -r requirements.txt
```
NOTE: Please make sure installation into the *mega_protact* environment. 

https://github.com/thelahunginjeet/pyrankagg/tree/master

```bash
python2 setup.py build
python2 setup.py install
```

**4 - Make predictions**

MEGA PROTAC consists of 14 Python scripts designed to investigate ternary structures for PROTACs. While some scripts may have similar functionalities, we have streamlined the process to facilitate the straightforward execution of the MEGA PROTAC analysis.


- The MEGA PROTAC initiates the MEGA DOCK process to identify the initial search space, which can also be referred to as seeds.
```bash
python3 STEP_1_MEGA_DOCK_Docking.py
```

- MEGA PROTAC filtrates protein-protein complexes based on their ligand distance.
```bash
python3 STEP_2_Ligand_Filtration_MEGA_DOCK_Seeds.py
```

- The second filtration has been implemented by protein-based docking, utilising MDAnalysis [5], SASA [6], Energy [7], and PIZSA [2].
```bash
python3 STEP_3_Filtrate_and_Rank_MEGA_DOCK_Seeds.py
```

NOTE: if you faced "ZeroDivisionError: integer division or modulo by zero" error. The PIZSA score probably has not worked properly. To solve this problem, go for:
```bash
mamba activate pızsa
```

Then, try to execute 
```bash
python2 ranking_PIZSA_cluster_rank.py -p -i
```
you can find help --help option

- Following filtration, MEGA PROTAC employs rank aggregation to prioritise protein complexes and identify the most favourable candidate for subsequent investigation. Here, the script selects only the **top 10** proteins for the sake of representation. These selected proteins have been translated to enlarge the search space.

```bash
python3 STEP_4_Take_Top_N_with_ligand_from_MEGA_DOCK_Seeds.py
```

- The same ligand filtration has been used to filtrate translated proteins.

```bash
python3 STEP_5_Ligand_Filtration_Translated.py
```

- The same protein-based filtration has been used to filtrate translated proteins.

```bash
python3 STEP_6_Filtrate_and_Rank_Translation.py
```

- In the seventh step, rank aggregation has been used to select the **top 10** translated protein to search on rotation. 

```bash
python3 STEP_7_Take_Top_N_Translated_Complex_with_ligand.py
```

- The same ligand filtration has been used to filtrate rotated proteins.

```bash
python3 STEP_8_Ligand_Filtration_For_Translated.py
```

- The same protein-based filtration has been used to filtrate rotated proteins.

```bash
python3 STEP_9_Filtrate_and_Rank_Rotation.py
```

- The selected proteins have been clustered using FCC.

```bash
python3 STEP_10_Clustering.py
```

- The selected clusters have been filtered based on energy.

```bash
python3 STEP_11_Cluster_Filtering.py
```

- The remaining protein has been reclustered.

```bash
python3 STEP_12_reClustering.py
```

- The reclustered proteins have been reranked using our rank aggregation.

```bash
python3 STEP_13_Cluster_ranking_clusters.py
```
- The final step is local docking, which can be done using any other molecular docking program, such as ZDock.
- MEGA PROTAC defaults to perform docking PROTAC into the first-ranked protein in the first cluster. Kindly note that the parameters have been reorganised to demonstrate the functioning of MEGA PROTAC. Hence, the chosen value of 10 must be increased to a minimum of 200, as specified in the paper.

```bash
python3 python3 STEP_14_PROTAC_docking.py
```

**Acknowledgements**

- We express our heartfelt gratitude to the creators of all the software components that constitute the PROTAC-Model [8] pipeline.
- The authors of BOTCP [9] are acknowledged for their generous contribution of data.



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
