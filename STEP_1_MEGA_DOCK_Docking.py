# -*- coding: utf-8 -*-
# 28/06/2024
# Author: Sadettin Y. Ugurlu & David McDonald

import pdb
import time
from pymol import cmd
from uti import mega_dock_optimisation
import sys
import glob
import os
import shutil
from statistics import mean
import pandas as pd
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBParser, PDBIO, Select
import math
import subprocess
import re
import numpy as np
from Bio import PDB
import random
from tqdm import tqdm

# Setting the seed for the random number generator to ensure reproducibility
# The seed value (42) is arbitrary but consistent, so the sequence of random numbers
# generated will be the same each time this code is run
random.seed(42)

# Get the current working directory
current_path = os.getcwd()

# MEGA DOCK path to execute the program
MEGA_DOCK_PATH = f"{current_path}/bin/MEGADOCK"


# Check_out_output_folder
folder_path = "outputs/Megadock_OUT"

# If the folder has not been there, create it.
if not os.path.exists(folder_path):
    os.makedirs(folder_path)
    print(f"{folder_path} formed!")
else:
    print(f"{folder_path} has been already there.")



def mega_dock_protein_protein(protein, number):
    # the paths for input files, which obtained by BOTCP study
    # https://www.sciencedirect.com/science/article/pii/S2667318523000168
    test_path=f"{current_path}/test_inputs/{protein}"
    os.chdir(test_path)
    ####################################################################################################################
    # PREPARE INPUTS: Input proteins should have ligand. Pymol has been used to prepare inputs
    cmd.reinitialize()
    cmd.load(f"{protein}_receptor.pdb")
    cmd.load("warhead.pdb")
    cmd.save(f"{protein}_receptor_input.pdb")
    cmd.remove("all")
    ##########################################
    cmd.reinitialize()
    cmd.load(f"{protein}_target.pdb")
    cmd.load("binder.pdb")
    cmd.save(f"{protein}_target_input.pdb")
    cmd.remove("all")
    # the paths for prepared inputs
    receptor=f"{current_path}/test_inputs/{protein}/{protein}_target_input.pdb"
    target=f"{current_path}/test_inputs/{protein}/{protein}_receptor_input.pdb"
    ####################################################################################################################
    # Execute MEGA DOCK to have intial search space
    mega_dock_optimisation.megadock(
        target=receptor,
        ligand=target,
        MEGA_DOCK_PATH=MEGA_DOCK_PATH,
        output_number=number,
        ouput_file_path=f"{current_path}/outputs/Megadock_OUT",
        # FFT_point=1,  # None
        voxel=3,  # while voxels becomes larger, it reduces times and lose performance for protein-protein docking
        degree=12,
        rotation=54000,
        # electrostatic=1,
        # hydrophobic=1,
        receptor_pel=-1000,
        ligand_pel=10,
        score=1
    )
    return f"{current_path}/outputs/Megadock_OUT/{protein}_target_input",

# MEGA DOCK provides separeted ligand-protein and target-protein files. Therefore, the complexes have been
# formed by using Pymol.
def complex_preparetion_with_ligand(protein, output_number):
    """directory_path = (f"{current_path}/outputs/Megadock_OUT/{protein}_target_input_0")
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

    input_path = f"/media/yavuz/New Volume/megadock_out/{protein}"
    source_path = f"/home/yavuz/yavuz_proje/protac/test_inputs/{protein}/{protein}_target_input.pdb"
    target_path = f"/home/yavuz/yavuz_proje/protac/outputs/Megadock_OUT/{protein}_target_input_0/{protein}_target_input_0.pdb"
    shutil.copy(source_path, target_path)"""

    source_path = f"{current_path}/test_inputs/{protein}/{protein}_receptor_input.pdb"
    target_path = f"{current_path}/outputs/Megadock_OUT/{protein}_target_input/{protein}_receptor_input.pdb"
    shutil.copy(source_path, target_path)

    target_path = f"{current_path}/outputs/Megadock_OUT/{protein}_target_input"
    os.chdir(target_path)
    print("Complexes have been merging using Pymol.")
    time.sleep(3)
    for i in tqdm(range(1, output_number), desc="Complexes construction", unit="  Output PDB file"):
        cmd.reinitialize()
        cmd.load(f"{protein}_target_input.pdb")
        cmd.load(f"ligand_{i}.pdb")
        cmd.load(f"{protein}_receptor_input.pdb")
        cmd.align(f"{protein}_receptor_input", f"ligand_{i}")
        cmd.remove(f"ligand_{i}")
        cmd.save(f"{protein}_target_input_{i}_complex.pdb")
        cmd.remove(f"{protein}_target_input")
        cmd.save(f"{protein}_target_input_{i}_ligand.pdb")
        cmd.remove("all")
    print("Done!")
def run_mega_dock(protein, output_number):
    # run MEGA DOCK.
    mega_dock_protein_protein(
        protein=protein, # the name of input protein
        number=output_number)  # the number of outputs
    # MEGA DOCK provides ligand protein pose without target protein.
    # Therefore, the ligand pose should be located into target protein.
    """complex_preparetion_with_ligand(
        output_path=f"{current_path}/outputs/Megadock_OUT/{protein}_target_input",
        number=output_number
    )"""
    complex_preparetion_with_ligand(
        protein=protein,
        output_number=output_number
    )

pdb_listem= [
    "5t35-DA", "5t35-HE", "6bn7-BC", "6boy-BC", "6hax-BA",
    "6hax-FE", "6hay-BA", "6hay-FE", "6hr2-BA", "6hr2-FE",
    "6sis-DA", "6sis-HE", "6w7o-CA", "6w7o-DB", "6w8i-DA",
    "6w8i-EB", "6w8i-FC", "6zhc-AD", "7jto-LB", "7jtp-LA",
    "7khh-CD", "7q2j-CD"
]
pdb_listem= ["5t35-DA"]
for protein in pdb_listem:
    run_mega_dock(protein,5000) # The defualt output number is 5000.
    # To perform a more exhaustive search:
    # 1) Increase the number of outputs.
    # 2) Execute MEGADOCK with different parameters.
