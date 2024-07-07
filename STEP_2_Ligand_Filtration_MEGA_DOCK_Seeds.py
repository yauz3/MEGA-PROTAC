# -*- coding: utf-8 -*-
# 28/06/2024
# Author: Sadettin Y. Ugurlu & David McDonald


import pdb
import time
from pymol import cmd
from uti import protein_uti
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
import time
import json
from tqdm import tqdm
from uti import uti_rotate_and_translate

# Start the timer by recording the current time
start_time = time.time()

# Setting the seed for the random number generator to ensure reproducibility
# The seed value (42) is arbitrary but consistent, so the sequence of random numbers
# generated will be the same each time this code is run
random.seed(42)

# Get the current working directory
current_path = os.getcwd()



# all Ternary structures
pdb_listem=[ "5t35-DA", "5t35-HE", "6bn7-BC", "6boy-BC", "6hax-BA",
    "6hax-FE", "6hay-BA", "6hay-FE", "6hr2-BA", "6hr2-FE",
    "6sis-DA", "6sis-HE", "6w7o-CA", "6w7o-DB", "6w8i-DA",
    "6w8i-EB", "6w8i-FC", "6zhc-AD", "7jto-LB", "7jtp-LA",
    "7khh-CD", "7q2j-CD"]

# selected the first one to demonstrates how MEGA PORTAC is working
pdb_listem=[ "5t35-DA"]
# with the help of for loop, MEGA PROTAC can be used for entire input files
for protein in pdb_listem:
    # find the input file
    input_path=f"{current_path}/outputs/Megadock_OUT/{protein}_target_input"
    os.chdir(input_path)
    # make a list using all COMPLEX.PDB files.
    pdb_list=glob.glob("*complex.pdb")
    print("The number of found prepared complexes:")
    print(len(pdb_list))
    output_dictionary={}
    filtered_list=[]
    # get the ligand protein chain ID to find PREDICTED WARHEAD location from PDB file.
    ligand_chain = protein_uti.get_chain_list(f"{current_path}/outputs/Megadock_OUT/{protein}_target_input/{protein}_target_input_1_ligand")[0]
    print("Ligand_based_filtration has been started")
    time.sleep(5)
    for complex in tqdm(pdb_list, desc="Ligand-Based Filtration", unit="pdb"):
        distance = uti_rotate_and_translate.distance_filtration(protein=protein,
                                       input_path=f"{current_path}/outputs/Megadock_OUT/{protein}_target_input",
                                       complex=complex,
                                       ligand_chain=ligand_chain)
        output_dictionary[complex]=distance
        # the limits to filter MEGA DOCK seeds
        if 3 < distance < 20:
            filtered_list.append(complex)
    print("\n\n\nReport:")
    print("Input number:", len(pdb_list))
    print("Remaning protein number:",len(filtered_list))
    with open(f'ligand_filtration.json', 'w') as f:
        json.dump(filtered_list, f)

end_time = time.time()
elapsed_time = end_time - start_time
print("Run time:", elapsed_time, "seconds")