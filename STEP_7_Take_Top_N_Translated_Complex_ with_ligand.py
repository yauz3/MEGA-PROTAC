# -*- coding: utf-8 -*-
# 28/06/2024
# Author: Sadettin Y. Ugurlu & David McDonald


import pdb
import time
from pymol import cmd
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
from uti import formation_searh_space_small

start_time = time.time()

# Setting the seed for the random number generator to ensure reproducibility
# The seed value (42) is arbitrary but consistent, so the sequence of random numbers
# generated will be the same each time this code is run
random.seed(42)

# Get the current working directory
current_path = os.getcwd()




def take_top_n_translated(protein, output_number):
    input_path = f"{current_path}/outputs/Megadock_OUT/{protein}_translate"
    os.chdir(input_path)
    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_translate/rank_aggregation.json", 'r') as f:
        ranking_output = json.load(f)
    index = 1
    new_sorted_dict={}
    selected_proteins=[]
    for key,value in ranking_output.items():
        new_sorted_dict[key] = index
        if index <= output_number:
            selected_proteins.append(key)
        index += 1
    print("selected")
    print(selected_proteins)
    #selected_proteins=["6hay-FE_target_input_0_4841_center_x_0-0_y_-4-5_z_3-0_center","6hay-FE_target_input_0_419_center_x_0-0_y_-4-5_z_3-0_center"]
    optmisation_path = f"{current_path}/outputs/Megadock_OUT/{protein}_rotate"

    # Dosya yolu var mı kontrol et
    if not os.path.exists(optmisation_path):
        # Klasörü oluştur
        os.makedirs(optmisation_path)
        print(f"{optmisation_path} oluşturuldu.")
    else:
        print(f"{optmisation_path} zaten var.")

    optimisation_input=[]
    for sel_pro in selected_proteins:
        source= f"{sel_pro}.pdb"
        optimisation_input.append(source)
        target=f"{optmisation_path}/{sel_pro}.pdb"
        shutil.copy(source,target)
    os.chdir(optmisation_path)
    formation_searh_space_small.rotational_grid_search(protein,optimisation_input,optmisation_path)


pdb_listem=[ "5t35-DA", "5t35-HE", "6bn7-BC", "6boy-BC", "6hax-BA",
    "6hax-FE", "6hay-BA", "6hay-FE", "6hr2-BA", "6hr2-FE",
    "6sis-DA", "6sis-HE", "6w7o-CA", "6w7o-DB", "6w8i-DA",
    "6w8i-EB", "6w8i-FC", "6zhc-AD", "7jto-LB", "7jtp-LA",
    "7khh-CD", "7q2j-CD"]


pdb_listem=["5t35-DA"]
for protein in pdb_listem:
    take_top_n_translated(protein,20)

end_time = time.time()
elapsed_time = end_time - start_time
print("İşlem süresi:", elapsed_time, "saniye")