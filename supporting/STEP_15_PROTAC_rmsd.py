# -*- coding: utf-8 -*-
# 28/06/2024
# Author: Sadettin Y. Ugurlu & David McDonald

import pdb
import time
from pymol import cmd
from script import protein_uti
from script import pocket_identification
from script import psvina
from script import mega_dock_optimisation
from script import lighdock
from script import zdock
from script import vina_uti
from script import protein_protein_selection
from script import clustering_fast
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
from script import uti_rotate_and_translate
from script import rotate_and_translate
from Bio import PDB
import random

# Setting the seed for the random number generator to ensure reproducibility
# The seed value (42) is arbitrary but consistent, so the sequence of random numbers
# generated will be the same each time this code is run
random.seed(42)

# Get the current working directory
current_path = os.getcwd()


pdb_listem=["7q2j-CD", "5t35-DA", "5t35-HE", "6bn7-BC", "6boy-BC", "6hax-BA",
    "6hax-FE", "6hay-BA", "6hay-FE", "6hr2-BA", "6hr2-FE",
    "6sis-DA", "6sis-HE", "6w7o-CA", "6w7o-DB", "6w8i-DA",
    "6w8i-EB", "6w8i-FC", "6zhc-AD", "7jto-LB", "7jtp-LA",
    "7khh-CD",]
pdb_listem=["6w8i-DA"]

protein_protein_list=["5t35-DA_target_input_0_122_center_x_-3-0_y_3-0_z_-4-5_center_x_-15-0_y_-10-0_z_-5-0_center",
                      "5t35-HE_target_input_0_144_center_x_-3-0_y_3-0_z_-4-5_center_x_-15-0_y_-10-0_z_-5-0_center",
                      "6sis-DA_target_input_0_1806_center_x_-1-5_y_-1-5_z_1-5_center_x_-5-0_y_-5-0_z_0-0_center",
                      "6sis-HE_target_input_0_153_center_x_-3-0_y_3-0_z_-4-5_center_x_-15-0_y_-10-0_z_-5-0_center",
                      "7jtp-LA_target_input_0_1128_center_x_-4-5_y_3-0_z_4-5_center_x_-5-0_y_-5-0_z_-15-0_center",
                      "7khh-CD_target_input_0_2441_center_x_1-5_y_-3-0_z_0-0_center_x_10-0_y_-5-0_z_15-0_center"]
for protein in protein_protein_list:
    directory_path = f"{current_path}/outputs/Megadock_OUT/{protein}"
    os.chdir(directory_path)
    ligand_list=glob.glob("*ligand*.pdb")
    rmsd_list=[]
    for ligand in ligand_list:
        ligand_name=protein.split("_")[0]
        input_path = f"{current_path}/Publication_Data_Final/Publication_Data/Input_data/native_TC/{ligand_name}_protac.pdb"
        output_path = f"{directory_path}/{ligand}"
        #output_path="/home/yavuz/deneme.pdb"
        completed_process = subprocess.run(["obrms", input_path, output_path], capture_output=True, text=True,
                                           check=True)
        completed_process = subprocess.run(["obrms", "-m", input_path, output_path], capture_output=True, text=True,
                                           check=True)
        # Print the captured output
        output_lines = completed_process.stdout.strip().split("\n")
        rmsd_line = output_lines[-1]  # The last line contains the RMSD value
        rmsd_value = float(rmsd_line.split()[-1])  # Extract the RMSD value
        # Print the RMSD value
        print("ligand:",ligand)
        print("RMSD value:", rmsd_value)
        rmsd_list.append(rmsd_value)
        #exit()
    print(rmsd_list)




