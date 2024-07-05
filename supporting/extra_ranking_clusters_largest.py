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
from uti import uti_rotate_and_translate
from uti import rotate_and_translate
from Bio import PDB
import random
from script import ranking_of_cluster_faster
import time
import json
import rank_aggregation

start_time = time.time()

# Setting the seed for the random number generator to ensure reproducibility
# The seed value (42) is arbitrary but consistent, so the sequence of random numbers
# generated will be the same each time this code is run
random.seed(42)

# Get the current working directory
current_path = os.getcwd()

def get_filtered_protein(protein):
    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_target_input/{protein}_scorelist.json", 'r') as f:
        score_list = json.load(f)
    vorquma_dic=score_list[-1]
    keys_list = list(vorquma_dic.keys())
    return keys_list





def group_values(protein):
    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_rotate/vorquma_dic_final.json", 'r') as f:
        vorquma_dic_cleaned = json.load(f)
    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_rotate/sasa_dic.json", 'r') as f:
        sasa_cleaned = json.load(f)
    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_rotate/mda_dic.json", 'r') as f:
        mda_cleaned = json.load(f)
    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_rotate/pizsa_dic.json", 'r') as f:
        pizsa_cleaned = json.load(f)
    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_rotate/energy_dic.json", 'r') as f:
        energy_cleaned = json.load(f)
    re_cluster = re.compile('Cluster (.*?) -> (.*?)\n', re.S)
    re_number = re.compile(f'(.*?)_{protein}_target_input_pro_comp_(\d+)')
    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_rotate/cluster_5_3_model_filtered", 'r') as cluster_file:
        cluster_lines = cluster_file.readlines()

    vorqmaaa_dictionary = {}
    sasa_dictionary = {}
    pizsa_dictionary = {}
    energy_dictionary= {}
    mda_dictionary = {}

    vorqmaaa_dictionary_0 = {}
    sasa_dictionary_0 = {}
    pizsa_dictionary_0 = {}
    energy_dictionary_0 = {}
    mda_dictionary_0 = {}

    vorqmaaa_dictionary_1 = {}
    sasa_dictionary_1 = {}
    pizsa_dictionary_1 = {}
    energy_dictionary_1 = {}
    mda_dictionary_1 = {}
    filtered_proteins={}
    for cluster_line in cluster_lines:
        voqruma_score_list = []
        sasa_score_list = []
        pizsa_list=[]
        mda_list=[]
        energy_list=[]
        re_line = re.search(re_cluster, cluster_line)
        if re_line:
            cluster_id, model = re_line.group(1), re_line.group(2)
            model_list = model.split()
            for pro in model_list:
                vorquma_score = vorquma_dic_cleaned.get(pro)
                sasa_score = sasa_cleaned.get(pro)
                pizsa_score=pizsa_cleaned.get(pro)
                mda_score=mda_cleaned.get(pro)
                energy_score=energy_cleaned.get(pro)
                mda_list.append(mda_score)
                energy_list.append(energy_score)
                voqruma_score_list.append(vorquma_score)
                sasa_score_list.append(sasa_score)
                pizsa_list.append(pizsa_score)
            vorqmaaa_dictionary[cluster_id] = max(voqruma_score_list)
            sasa_dictionary[cluster_id]=max(sasa_score_list)
            pizsa_dictionary[cluster_id]=-int(cluster_id)
            energy_dictionary[cluster_id]=max(energy_list)*-1
            mda_dictionary[cluster_id]=max(mda_list)

            vorqmaaa_dictionary_1[cluster_id] = min(voqruma_score_list)
            sasa_dictionary_1[cluster_id] = min(sasa_score_list)
            pizsa_dictionary_1[cluster_id] = min(pizsa_list)
            energy_dictionary_1[cluster_id] = min(energy_list)*-1
            mda_dictionary_1[cluster_id] = min(mda_list)

            vorqmaaa_dictionary_0[cluster_id] = mean(voqruma_score_list)
            sasa_dictionary_0[cluster_id] = mean(sasa_score_list)
            pizsa_dictionary_0[cluster_id] = mean(pizsa_list)
            energy_dictionary_0[cluster_id] = mean(energy_list)*-1
            mda_dictionary_0[cluster_id] = mean(mda_list)
            filtered_proteins[cluster_id]=model_list
    return vorqmaaa_dictionary, sasa_dictionary, pizsa_dictionary, energy_dictionary,mda_dictionary,vorqmaaa_dictionary_0,sasa_dictionary_0,pizsa_dictionary_0, energy_dictionary_0, mda_dictionary_0, vorqmaaa_dictionary_1,sasa_dictionary_1,pizsa_dictionary_1, energy_dictionary_1, mda_dictionary_1, filtered_proteins

pdb_listem=[ "5t35-DA", "5t35-HE", "6bn7-BC", "6boy-BC", "6hax-BA",
            "6hax-FE", "6hay-BA", "6hay-FE", "6hr2-BA", "6hr2-FE",
            "6sis-DA", "6sis-HE", "6w7o-CA", "6w7o-DB", "6w8i-DA",
            "6w8i-EB", "6w8i-FC", "6zhc-AD", "7jto-LB", "7jtp-LA",
            "7khh-CD", "7q2j-CD"]


for protein in pdb_listem:
    input_path = f"{current_path}/outputs/Megadock_OUT/{protein}_rotate"
    score_list=[group_values(protein)[0],
                group_values(protein)[1],
                group_values(protein)[2],
               ]
    #print(score_list)
    os.chdir(input_path)
    with open(f'group_scorelist_filtered_vor+sasa+largest.json', 'w') as f:
        json.dump(score_list, f)
    # lowest score first here, lowest, 1st should be perfect result
    vorqma_pizsa = rank_aggregation.perform_rank_aggregation(score_list)
    print(protein)
    print(vorqma_pizsa)
    with open(f'group_scorelist_filtered_vor+sasa+largest.json', 'w') as f:
        json.dump(vorqma_pizsa, f)
end_time = time.time()
elapsed_time = end_time - start_time
print("İşlem süresi:", elapsed_time, "saniye")