# -*- coding: utf-8 -*-
# 28/06/2024
# Author: Sadettin Y. Ugurlu & David McDonald


import pdb
import time
from pymol import cmd
from uti import clustering_fast
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

start_time = time.time()

# Setting the seed for the random number generator to ensure reproducibility
# The seed value (42) is arbitrary but consistent, so the sequence of random numbers
# generated will be the same each time this code is run
random.seed(42)

# Get the current working directory
current_path = os.getcwd()



pdb_listem=["5t35-DA", "5t35-HE", "6bn7-BC", "6boy-BC", "6hax-BA",
    "6hax-FE", "6hay-BA", "6hay-FE", "6hr2-BA", "6hr2-FE",
    "6sis-DA", "6sis-HE", "6w7o-CA", "6w7o-DB", "6w8i-DA",
    "6w8i-EB", "6w8i-FC", "6zhc-AD", "7jto-LB", "7jtp-LA",
    "7khh-CD","7q2j-CD", ]
pdb_listem=["5t35-DA"]
for protein in pdb_listem:
    input_path = f"{current_path}/outputs/Megadock_OUT/{protein}_rotate"
    os.chdir(input_path)
    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_rotate/vorquma_dic_final.json", 'r') as f:
        vorquma_dic_cleaned = json.load(f)
    pdb_without_pdb=list(vorquma_dic_cleaned.keys())


    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_rotate/group_filtered.json",
              'r') as f:
        cluster_filtered_proteins = json.load(f)
    pdb_without_pdb=cluster_filtered_proteins

    input_file_name=clustering_fast.prepare_input_file(input_path=input_path, file_list=pdb_without_pdb)
    cluster_output_name = "cluster_5_3_model_final"  # Adding chunk number to the output name
    cluster_5_3_model = clustering_fast.cluster(input_path=input_path, output_name=cluster_output_name,
                                                min_number=2, threshold=0.5)
    #########################


end_time = time.time()
elapsed_time = end_time - start_time
print("İşlem süresi:", elapsed_time, "saniye")