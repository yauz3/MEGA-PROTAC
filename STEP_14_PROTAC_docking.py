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
import json
import rank_aggregation


random.seed(42)

# Setting the seed for the random number generator to ensure reproducibility
# The seed value (42) is arbitrary but consistent, so the sequence of random numbers
# generated will be the same each time this code is run
random.seed(42)

# Get the current working directory
current_path = os.getcwd()

# MEGA DOCK path to execute the program
MEGA_DOCK_PATH = f"{current_path}/bin/MEGADOCK-master"


def protac_into_protein_protein(receptor, target):
    test_path=f"{current_path}/test_inputs"
    os.chdir(test_path)
    mega_dock_optimisation.megadock(
        target=target,
        ligand=receptor,
        MEGA_DOCK_PATH=MEGA_DOCK_PATH,
        output_number=1000,
        ouput_file_path=f"{current_path}/outputs/Megadock_OUT",
        # FFT_point=1,  # None
        voxel=1.2,  # 1 553  0.8 653
        degree=None,  # 6 fark etmiyor
        rotation=24,  # 54000, 24 fark etmiyor
        # electrostatic=1,  # 1
        # hydrophobic=1,  # 1
        receptor_pel=-45,  # -30
        ligand_pel=1,  # 1
        score=3
    )


def get_cluster_proteins(protein, file_input):
    output_dictionary={}
    input_path=f"{current_path}/outputs/Megadock_OUT/{protein}_rotate"
    os.chdir(input_path)
    re_cluster = re.compile('Cluster (.*?) -> (.*?)\n', re.S)

    with open(file_input, 'r') as cluster_file:
        cluster_lines = cluster_file.readlines()
    cluster_number=1
    for cluster_line in cluster_lines:
        re_line = re.search(re_cluster, cluster_line)
        if re_line:
            model = re_line.group(2)
            model_list = model.split()
            #print("model_list")
            #print(model_list)
            output_dictionary[cluster_number]=model_list
            cluster_number += 1
    return output_dictionary
#out=get_cluster_proteins("5t35-DA","cluster_5_3_model_final")


def ranks_protein_in_cluster(protein,input_protein_list):

    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_rotate/vorquma_dic_final.json",
                  'r') as f:
        vor_dic = json.load(f)
    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_rotate/sasa_dic.json",
              'r') as f:
        sasa_dic = json.load(f)
    vor_cleaned = {key: vor_dic[key] for key in input_protein_list}
    sasa_cleaned = {key: sasa_dic[key] for key in input_protein_list}
    score_list = [sasa_cleaned,
                  vor_cleaned]
    vorqma_pizsa = rank_aggregation.perform_rank_aggregation(score_list)
    #print(vorqma_pizsa)
    return vorqma_pizsa
#print(ranks_protein_in_cluster("5t35-DA",output_dictionary))

def find_protein(input_dictionary, max_value):
    result_list = []
    for key, val in input_dictionary.items():
        if val <= max_value:
            result_list.append(key)
    return result_list


def get_first_cluster(protein,cluster_order,protein_number):
    input_path = f"{current_path}/outputs/Megadock_OUT/{protein}_rotate"
    os.chdir(input_path)
    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_rotate/group_rank_aggregation_vor+sasa_selected.json", 'r') as f:
        ordered_clusters = json.load(f)
    #print(ordered_clusters)
    for key, value in ordered_clusters.items():
        if value == cluster_order:
            wanted_cluster_number=key

    clusters_proteins = get_cluster_proteins(protein, "cluster_5_3_model_final")
    proteins_in_selected_cluster = clusters_proteins.get(int(wanted_cluster_number))
    protein_orders=ranks_protein_in_cluster(protein,proteins_in_selected_cluster)
    ordered_proteins_for_PROTAC_docking=find_protein(protein_orders,protein_number)
    print(ordered_proteins_for_PROTAC_docking)
    for protein in ordered_proteins_for_PROTAC_docking:
        input_protein = protein.split("_")[0]
        # PROTAC docking
        protac_into_protein_protein(
            target=f"{current_path}/outputs/Megadock_OUT/{protein}_rotate",
            receptor=f"{current_path}/test_inputs/{input_protein}/{input_protein}_protac.pdb")

get_first_cluster(protein="5t35-DA",
                  cluster_order=1,
                  protein_number=1)








