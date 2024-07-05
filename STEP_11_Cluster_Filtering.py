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
import rank_aggregation

start_time = time.time()

# Setting the seed for the random number generator to ensure reproducibility
# The seed value (42) is arbitrary but consistent, so the sequence of random numbers
# generated will be the same each time this code is run
random.seed(42)

# Get the current working directory
current_path = os.getcwd()


def get_cluster_number(input_path):
    cluster_numbers = set()  # Benzersiz küme numaralarını saklamak için bir küme oluştur
    with open(f"{input_path}/cluster_5_3_model_before_filtration", 'r') as cluster_file:
        for line in cluster_file:
            # Her satırın başındaki "Cluster" ifadesini arayarak küme numarasını bul
            if line.startswith('Cluster'):
                cluster_number = line.split()[1]  # Küme numarasını al
                cluster_numbers.add(cluster_number)  # Küme numarasını küme içine ekle
    # Benzersiz küme numaralarını sıralayarak listeye dönüştür
    sorted_cluster_numbers = sorted(cluster_numbers, key=int)
    return sorted_cluster_numbers


def filtared_by_group(input_path, cluster_input,selected_proteins):
    re_cluster = re.compile('Cluster (.*?) -> (.*?)\n', re.S)
    with open(f"{input_path}/cluster_5_3_model_before_filtration", 'r') as cluster_file:
        cluster_lines = cluster_file.readlines()
    vorqmaaa_dictionary = {}
    proteins_dictionary = {}
    for cluster_line in cluster_lines:
        protein_list = []
        voqruma_score_list = []
        re_line = re.search(re_cluster, cluster_line)
        if re_line:
            cluster_id, model = re_line.group(1), re_line.group(2)
            if cluster_id in cluster_input:
                model_list = model.split()
                if any(element in selected_proteins for element in model_list):
                        proteins_dictionary[cluster_id] = model_list
    return proteins_dictionary


pdb_listem=[ "5t35-DA", "5t35-HE", "6bn7-BC",
             "6boy-BC", "6hax-BA", "6hax-FE",
             "6hay-BA", "6hay-FE", "6hr2-BA",
             "6hr2-FE", "6sis-DA", "6sis-HE",
             "6w7o-CA", "6w7o-DB", "6w8i-DA",
             "6w8i-EB", "6w8i-FC", "6zhc-AD",
             "7jto-LB", "7jtp-LA", "7khh-CD",
             "7q2j-CD"]

pdb_listem=["5t35-DA"]

for protein in pdb_listem:
    input_path=f"{current_path}/outputs/Megadock_OUT/{protein}_rotate"

    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_rotate/mda_dic.json",
              'r') as f:
        mda_dic = json.load(f)
    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_rotate/sasa_dic.json",
              'r') as f:
        sasa_dic = json.load(f)
    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_rotate/energy_dic.json",
              'r') as f:
        energy_dic = json.load(f)
    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_rotate/pizsa_dic.json",
              'r') as f:
        pizsa_dic = json.load(f)
    with open(f"{current_path}/outputs/Megadock_OUT/{protein}_rotate/vorquma_dic_final.json",
              'r') as f:
        vor_dic = json.load(f)

    filtered_keys = set(vor_dic.keys())
    pizsa_cleaned = {key: pizsa_dic[key] for key in filtered_keys}
    energy_cleaned = {key: energy_dic[key] for key in filtered_keys}
    for key in energy_cleaned:
        energy_cleaned[key] = -energy_cleaned[key]
    sasa_cleaned = {key: sasa_dic[key] for key in filtered_keys}
    mda_cleaned = {key: mda_dic[key] for key in filtered_keys}
    score_list=[energy_cleaned]
    vorqma_pizsa = rank_aggregation.perform_rank_aggregation(score_list)
    sorted_energy = dict(sorted(energy_cleaned.items(), key=lambda item: item[1], reverse=True))
    # sorted_energy is equal to vorqma_pizsa

    # Dictionary'nin ilk %10'unu alın
    percent_to_keep = 0.25
    num_to_keep = int(len(vorqma_pizsa) * percent_to_keep)
    # Slice kullanarak liste elde edin ve ardından dictionary'ye dönüştürün
    top_10_percent = dict(list(vorqma_pizsa.items())[:num_to_keep])
    print("total")
    print(len(vorqma_pizsa))
    print("top_10_percent")
    """print(top_10_percent)
    print(list(top_10_percent.keys()))"""
    print(len(top_10_percent))
    total_cluster_indexs = get_cluster_number(input_path)
    fitrated_groups=filtared_by_group(input_path, total_cluster_indexs,top_10_percent.keys())
    selected_proteins=[]
    for key,value in fitrated_groups.items():
        for prot in value:
            selected_proteins.append(prot)
    print(len(selected_proteins))
    os.chdir(input_path)
    with open(f'group_filtered.json', 'w') as f:
        json.dump(selected_proteins, f)
    #time.sleep(5)
end_time = time.time()
elapsed_time = end_time - start_time
print("İşlem süresi:", elapsed_time, "saniye")