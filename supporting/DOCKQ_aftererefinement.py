# -*- coding: utf-8 -*-
# 28/06/2024
# Author: Sadettin Y. Ugurlu & David McDonald

import pdb
import time
from pymol import cmd
from uti import protein_uti
from uti import pocket_identification
from uti import psvina
from uti import mega_dock_optimisation
from uti import lighdock
from uti import zdock
from uti import vina_uti
from uti import protein_protein_selection
from uti import dockq_rmsd
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

protein_protein_list=codes = [
    "5t35-DA", "5t35-HE", "6bn7-BC", "6boy-BC", "6hax-BA",
    "6hax-FE", "6hay-BA", "6hay-FE", "6hr2-BA", "6hr2-FE",
    "6sis-DA", "6sis-HE", "6w7o-CA", "6w7o-DB", "6w8i-DA",
    "6w8i-EB", "6w8i-FC", "6zhc-AD", "7jto-LB", "7jtp-LA",
    "7khh-CD", "7q2j-CD"
]

# Get the current working directory
current_path = os.getcwd()

protein_protein_list=["5t35-DA"]

max_dockq_score=0
output_dictionary={}
for protein in protein_protein_list:
    true_ternary=f"{current_path}/Publication_Data_Final/Publication_Data/Input_data/native_TC"
    predicted=f"{current_path}/Publication_Data_Final/Publication_Data/Refinement_Publication_data"
    os.chdir(
        f"{current_path}/Publication_Data_Final/Publication_Data/Refinement_Publication_data/{protein}")
    proteins = glob.glob("complex*protein_protein*")
    rmsd_list=[]
    for pro in proteins:
        # print(protein)
        try:
            fnat, fnonnat, irms, lrms, dockq = dockq_rmsd.rmsd_calculation(f"{true_ternary}/{protein}_Native.pdb",
                                                                           f"{current_path}/Publication_Data_Final/Publication_Data/Refinement_Publication_data/{protein}/{pro}"
                                                                           )
            if float(dockq) > 0.22:
                print(pro)
                print("fnat", fnat)
                print("fnonnat", fnonnat)
                print("irms", irms)
                print("lrms", lrms)
                print("dockq", dockq)
                rmsd_list.append(dockq)
                if dockq > max_dockq_score:
                    max_dockq_score = dockq

        except:
            continue
    max_rmsd = max(rmsd_list)
    output_dictionary[protein] = max_rmsd
print(max_dockq_score)