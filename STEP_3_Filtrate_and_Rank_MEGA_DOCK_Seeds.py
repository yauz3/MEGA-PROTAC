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
from uti import ranking_of_cluster_faster
import time
import json
import rank_aggregation
from tqdm import tqdm


start_time = time.time()

# Setting the seed for the random number generator to ensure reproducibility
# The seed value (42) is arbitrary but consistent, so the sequence of random numbers
# generated will be the same each time this code is run
random.seed(42)

# Get the current working directory
current_path = os.getcwd()

#  In order to keep their outputs from main outputs, PIZSA output should be stored at different folder.
pizsa_path = f"{current_path}/outputs/Megadock_OUT/pizsa"
# Check if the directory exists
if not os.path.exists(pizsa_path):
    # If it doesn't exist, create the directory
    os.makedirs(pizsa_path)

# all Ternary structures
pdb_listem=[ "5t35-DA", "5t35-HE", "6bn7-BC", "6boy-BC", "6hax-BA",
    "6hax-FE", "6hay-BA", "6hay-FE", "6hr2-BA", "6hr2-FE",
    "6sis-DA", "6sis-HE", "6w7o-CA", "6w7o-DB", "6w8i-DA",
    "6w8i-EB", "6w8i-FC", "6zhc-AD", "7jto-LB", "7jtp-LA",
    "7khh-CD", "7q2j-CD"]

# selected the first one to demonstrates how MEGA PORTAC is working
pdb_listem=[  "5t35-DA"]

for protein in pdb_listem:
    # find the input file
    input_path = f"{current_path}/outputs/Megadock_OUT/{protein}_target_input"
    os.chdir(input_path)
    # get the filtered MEGA DOCK outputs from output json file.
    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_target_input/ligand_filtration.json",
              'r') as f:
        pdb_list = json.load(f)

    # get rid of extension of files
    pdb_without_pdb=[]
    for i in pdb_list:
        pdb_without_pdb.append(i.replace(".pdb",""))
    #the number of input file
    print("Filtered proteins have been cleaning from het atoms")
    time.sleep(3)
    for i in tqdm(pdb_without_pdb, desc="Removing HET atoms", unit="pdb"):
        cmd.reinitialize()
        cmd.load(f"{i}.pdb")
        cmd.remove("het")
        cmd.save(f"{i}.pdb")
    #############################################################################################################################
    print("MDA based filtration has been started!")
    mda_filtered_list, mda_dic = ranking_of_cluster_faster.MDA_filtering(pdb_without_pdb, input_path, 0.2)
    with open(f'mda_dic.json', 'w') as f:
        json.dump(mda_dic, f)
    #############################################################################################################################
    print("SASA based filtration has been started!")
    sasa_filtered_list,sasa_dic = ranking_of_cluster_faster.sasa_filtering(mda_filtered_list, input_path,0.5) # 10-10
    with open(f'sasa_dic.json', 'w') as f:
        json.dump(sasa_dic, f)
    #############################################################################################################################
    print("Energy based filtration has been started!")
    energy_filtered_list,energy_dic = ranking_of_cluster_faster.energy_filtering_new(protein,sasa_filtered_list, input_path)

    with open(f'energy_dic.json', 'w') as f:
        json.dump(energy_dic, f)
    #############################################################################################################################
    print("PIZSA based filtration has been started!")
    # For MEGA DOCK seed filtration, PIZSA threshold was 0.5. Since higher threshold may be reason to lose promissing
    # candidates.
    # If the output number is much higher than 5000 MEGA DOCK outputs, or use 10 times different MEGA DOCK parameters, the theshold should be 1.
    pizsa_filtered_list,pizsa_dic = ranking_of_cluster_faster.pizsa_filtering_with_threshold(energy_filtered_list, input_path, 0.5)
    with open(f'pizsa_dic.json', 'w') as f:
        json.dump(pizsa_dic, f)
    # because of the same reasons dicussed above, the stability check has not been applied for 5000 MEGA DOCK outputs.
    """#########################################################
    os.chdir(pizsa_path)
    pizsa_output_files=[]
    for pdb_file_name in pizsa_filtered_list:
        pizsa_output_files.append(f"{pdb_file_name}._scores.txt{pdb_file_name}._scores.txt")
    stability_filtered=ranking_of_cluster_faster.pizsa_stability_filtration(pizsa_output_files)
    print("stability_filrted")
    print(len(stability_filtered))
    print(stability_filtered)"""
    #############################################################################################################################
    print("Vorquma assessment has been started!")
    # Vorquma is used in rank aggregation to order proteins. In other words, complex having the low vorquma score
    # will be located at the bottom of ranks, so they will be filtered because of rank aggregation.
    # Therefore, any protein has not been filtered based on vorquma.
    # On the other hand, 0.0 input can be use for any user to filtrate some proteins.
    # For example, 0.2 will filtrates 20% of proteins
    os.chdir(input_path)
    vorquma_filtered_list, vorquma_dic = ranking_of_cluster_faster.vormqa_filtering(pizsa_filtered_list,
                                                                                    input_path,
                                                                                    0.0)  # 0
    # write output into json file
    with open(f'vorquma_dic_final.json', 'w') as f:
        json.dump(vorquma_dic, f)

    # get the final filtered protein list
    filtered_without_extension = list(vorquma_dic.keys())

    # check protein and filter out proteins have been filtered in other stages.
    # after this filtration, any one can use different rank aggregation options, such as use all ranking powers.
    vorquma_dic_cleaned = {key: vorquma_dic[key] for key in filtered_without_extension}
    mda_cleaned = {key: mda_dic[key] for key in filtered_without_extension}
    sasa_cleaned = {key: sasa_dic[key] for key in filtered_without_extension}
    energy_cleaned = {key: energy_dic[key] for key in filtered_without_extension}
    # make sure energy values are reverse. Since rank aggregation select the highest value, but lowest energy is favorubale
    for key in energy_cleaned:
        energy_cleaned[key] = -energy_cleaned[key]
    pizsa_cleaned= {key: pizsa_dic[key] for key in filtered_without_extension}

    # MEGA PROTAC uses only two of them in rank aggregation to increase all over performace
    scorelist = [
        vorquma_dic_cleaned,
        #mda_cleaned,
        sasa_cleaned,
        #energy_cleaned,
        #pizsa_cleaned
    ]
    # write out scorelist into a json file
    with open(f'scorelist.json', 'w') as f:
        json.dump(scorelist, f)

    # lowest score first here, lowest, 1st should be perfect result
    vorqma_pizsa = rank_aggregation.perform_rank_aggregation(scorelist)
    print("Filtration has been finished!")
    print("######################")
    print("Ranked list")
    print(vorqma_pizsa)
    with open(f'rank_aggregation.json', 'w') as f:
        json.dump(vorqma_pizsa, f)
    print("\n\n\n")
    print("######################")
    print("Summary:")
    print("The number of protein after ligand-based filtration:",len(pdb_without_pdb))
    print("MDA filtration number:", len(mda_filtered_list))
    print("SASA filtration number:", len(sasa_filtered_list))
    print("Energy filtration number:", len(energy_filtered_list))
    print("PIZSA filtration number:",len(pizsa_filtered_list))
    print("Vorquma assessment number:",len(vorquma_filtered_list))

end_time = time.time()
elapsed_time = end_time - start_time
print("Run Time:", elapsed_time, "Seconds")
