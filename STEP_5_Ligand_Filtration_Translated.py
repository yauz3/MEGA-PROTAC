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

# Start the timer by recording the current time
start_time = time.time()

# Setting the seed for the random number generator to ensure reproducibility
# The seed value (42) is arbitrary but consistent, so the sequence of random numbers
# generated will be the same each time this code is run
random.seed(42)

# Get the current working directory
current_path = os.getcwd()


def distance_filtration(
        protein,
        input_path,
        complex,
        ligand_chain):
    os.chdir(input_path)

    # these are fixed points
    warhead_coordinates = protein_uti.identify_centre_of_mass(f"{current_path}/test_inputs/{protein}/binder.pdb")

    """print("pdb_file")
    print(complex)"""
    if "_complex.pdb" in complex:
        ligand_name = complex.replace("_complex.pdb", "")
    cmd.reinitialize()
    cmd.load(f"{ligand_name}_ligand.pdb")
    #start_point_1 = cmd.centerofmass(selection=f"chain {ligand_chain}")
    anchor_coordinate = cmd.centerofmass(selection=f"het")
    cmd.remove("all")
    dist = math.sqrt(sum([(a - b) ** 2 for a, b in zip(warhead_coordinates, anchor_coordinate)]))
    #print(dist)
    return dist




pdb_listem=[ "5t35-DA", "5t35-HE", "6bn7-BC", "6boy-BC", "6hax-BA",
    "6hax-FE", "6hay-BA", "6hay-FE", "6hr2-BA", "6hr2-FE",
    "6sis-DA", "6sis-HE", "6w7o-CA", "6w7o-DB", "6w8i-DA",
    "6w8i-EB", "6w8i-FC", "6zhc-AD", "7jto-LB", "7jtp-LA",
    "7khh-CD", "7q2j-CD"]

pdb_listem=[ "5t35-DA"]
for protein in pdb_listem:
    input_path=f"{current_path}/outputs/Megadock_OUT/{protein}_translate"
    os.chdir(input_path)
    pdb_list=glob.glob("*x*y*z*_complex.pdb")
    print("The number of found prepared complexes:")
    print(len(pdb_list))
    output_dictionary={}
    filtered_list=[]
    ligand_chain = protein_uti.get_chain_list(f"{current_path}/outputs/Megadock_OUT/{protein}_target_input/{protein}_target_input_1_ligand")[0]
    time.sleep(4)
    for complex in tqdm(pdb_list, desc="Ligand-Based Filtration", unit="pdb"):
        distance = distance_filtration(protein=protein,
                                       input_path=f"{current_path}/outputs/Megadock_OUT/{protein}_translate",
                                       complex=complex,
                                       ligand_chain=ligand_chain)
        output_dictionary[complex]=distance
        if 3 < distance < 20:
            filtered_list.append(complex)
    print("\n\n\nReport:")
    print("Ligand_based_filtration has been finished!")
    print("Input number:", len(pdb_list))
    print("Remaning protein number:", len(filtered_list))
    with open(f'ligand_filtration_translated.json', 'w') as f:
        json.dump(filtered_list, f)

end_time = time.time()
elapsed_time = end_time - start_time
print("Run time:", elapsed_time, "seconds")