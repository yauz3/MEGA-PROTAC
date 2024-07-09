# -*- coding: utf-8 -*-
# 28/06/2024
# Author: Sadettin Y. Ugurlu & David McDonald

import os
from uti import MDAnalysis
import glob
import subprocess
import numpy as np
import freesasa
import re
import pandas as pd
import os
import glob
from pymol import cmd
import rank_aggregation
import concurrent.futures
import json

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir)) # ..
FCC_path = os.path.join(PROJECT_ROOT, "bin/fcc/scripts")
voronota_path = os.path.join(PROJECT_ROOT, "bin/voronota")
fcc= os.path.join(PROJECT_ROOT, "bin/fcc")

########################################################################################################################
def MDAnalysis_filter(proteins,input_path):
    os.chdir(input_path)
    output_dictionary_MDA={}
    for pro in proteins:
        output_dictionary_MDA[pro]=float(MDAnalysis.quality_score_mda(f"{pro}.pdb"))
    sorted_dict_mda = dict(sorted(output_dictionary_MDA.items(), key=lambda item: item[1], reverse=True))
    return output_dictionary_MDA

def MDA_filtering(protein_list,input_path,filtering_ratio):
    if len(protein_list) > 500:
        num_chunks = 500  # Toplam işlem sayısı
    else:
        num_chunks = int(len(protein_list)/5)
    chunk_size = len(protein_list) // num_chunks  # Her bir işlem parçasının boyutu
    chunks = [protein_list[i:i + chunk_size] for i in range(0, len(protein_list), chunk_size)]

    # Paralel işlemler için ProcessPoolExecutor kullanımı
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = {}
        futures["mda"] = {executor.submit(MDAnalysis_filter, chunk, input_path) for chunk in chunks}

        scores = {}
        for func_name, future_set in futures.items():
            scores[func_name] = [future.result() for future in future_set]

        # Birleştirilmiş sonuçları kullanarak işlevlerin sonuçlarını oluşturun
    mda_filtered_1 = {}
    for score_dict in scores["mda"]:
        mda_filtered_1.update(score_dict)
    sorted_dict_mda = dict(sorted(mda_filtered_1.items(), key=lambda item: item[1], reverse=True))
    end_index = int(len(sorted_dict_mda) * (1 - filtering_ratio))
    selected = {k: v for i, (k, v) in enumerate(sorted_dict_mda.items()) if i < end_index}
    return list(selected.keys()), selected
########################################################################################################################

def sasa_feature(protein):
    structure = freesasa.Structure(protein)
    result = freesasa.calc(structure,
                           freesasa.Parameters({'algorithm': freesasa.LeeRichards,
                                                'n-slices': 100}))
    area_classes = freesasa.classifyResults(result, structure)

    areas = []
    for key in area_classes:
        areas.append(area_classes[key])
        # print(key, ": %.2f A2" % area_classes[key])

    # atom based-areas
    class DerivedClassifier(freesasa.Classifier):
        # this must be set explicitly in all derived classifiers
        purePython = True

        def classify(self, residueName, atomName):
            if re.match('\s*N', atomName):
                return 'Nitrogen'
            if re.match('\s*O', atomName):
                return 'Oxygen'
            if re.match('\s*S', atomName):
                return 'Sulfur'
            if re.match('\s*C', atomName):
                return 'Carbon'
            return 'Not-selected'

        def radius(self, residueName, atomName):
            if re.match('\s*N', atomName):  # Nitrogen
                return 1.6
            if re.match('\s*C', atomName):  # Carbon
                return 1.7
            if re.match('\s*O', atomName):  # Oxygen
                return 1.4
            if re.match('\s*S', atomName):  # Sulfur
                return 1.8
            return 0;  # everything else (Hydrogen, etc)

    classifier = DerivedClassifier()

    # use the DerivedClassifier to calculate atom radii
    structure = freesasa.Structure(protein, classifier)
    result = freesasa.calc(structure,
                           freesasa.Parameters({'algorithm': freesasa.LeeRichards,
                                                'n-slices': 100}))
    # use the DerivedClassifier to classify atoms
    area_classes = freesasa.classifyResults(result, structure, classifier)
    #print(result.totalArea(), areas[0],areas[1],area_classes.get("Nitrogen"),area_classes.get('Carbon'),area_classes.get('Oxygen'),area_classes.get('Sulfur'))
    return result.totalArea(), areas[0],areas[1]

def sasa_to_list(proteins, input_path):
    os.chdir(input_path)
    os.chdir(input_path)
    output_dictionary_sasa = {}
    for prot in proteins:
        total_area = sasa_feature(f"{prot}.pdb")[0]
        output_dictionary_sasa[prot] = total_area

    sorted_dict_sasa = dict(sorted(output_dictionary_sasa.items(), key=lambda item: item[1], reverse=True))
    return sorted_dict_sasa

def sasa_filtering(protein_list,input_path,filtering_ratio):
    if len(protein_list) > 500:
        num_chunks = 500  # Toplam işlem sayısı
    else:
        num_chunks = int(len(protein_list))
    chunk_size = len(protein_list) // num_chunks  # Her bir işlem parçasının boyutu
    chunks = [protein_list[i:i + chunk_size] for i in range(0, len(protein_list), chunk_size)]

    # Paralel işlemler için ProcessPoolExecutor kullanımı
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = {}
        futures["sasa"] = {executor.submit(sasa_to_list, chunk, input_path) for chunk in chunks}

        scores = {}
        for func_name, future_set in futures.items():
            scores[func_name] = [future.result() for future in future_set]

        # Birleştirilmiş sonuçları kullanarak işlevlerin sonuçlarını oluşturun
    sasa_filtered_1 = {}
    for score_dict in scores["sasa"]:
        sasa_filtered_1.update(score_dict)

    # en büyük değer filtreden sonra duruyor. bir sıkıntı yok.
    sorted_dict_sasa = dict(sorted(sasa_filtered_1.items(), key=lambda item: item[1], reverse=True))
    """print("sorted_dictionary_sasa")
    print(sorted_dict_sasa)"""
    end_index = int(len(sorted_dict_sasa) * (1 - filtering_ratio))
    selected = {k: v for i, (k, v) in enumerate(sorted_dict_sasa.items()) if i < end_index}
    """print(selected)"""
    return list(selected.keys()), selected
########################################################################################################################

def obenergy(protein,scoring_function="GAFF"):
    """
    avaliable scoring functions:
    GAFF :General Amber Force Field (GAFF)
    Ghemical :Ghemical force field
    MMFF94 :MMFF94 force field
    MMFF94s :MMFF94s force field
    UFF :Universal Force Field.

    """

    try:
        if scoring_function == "GAFF" or scoring_function == "Ghemical" or scoring_function == "UFF":
            command = f"obenergy -ff {scoring_function} {protein}.pdb"
            result = str(subprocess.check_output(command,
                                                 shell=True))
            energies = result.split(r"E N E R G Y")[1].split(" kJ")

        else:
            if scoring_function == "MMFF94":
                command = f"obenergy -ff  MMFF94 {protein}.pdb"
            else:
                command = f"obenergy {protein}.pdb"
            result = str(subprocess.check_output(command,
                                                 shell=True))
            energies = result.split(r"E N E R G Y")[1].split(" kcal")
        # print(energies)
        energy_list = []
        for energy in energies:
            if "=" in energy:
                #print(energy.split(" = ")[0].split("TOTAL")[1])
                energy_list.append(energy.split(" = ")[1])
    except:
        if scoring_function == "GAFF":
            energy_list = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
        elif scoring_function == "UFF":
            energy_list = [np.nan, np.nan, np.nan, np.nan, np.nan]
        elif scoring_function == "Ghemical":
            energy_list = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]
        elif scoring_function == "MMFF94s" or scoring_function == "MMFF94":
            energy_list = [np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan]

    return energy_list

def find_test_inputs_directory(start_path):
    path = start_path
    while path != os.path.dirname(path):  # Kök dizine ulaşana kadar devam et
        uti_path = os.path.join(path, "test_inputs")
        if os.path.exists(uti_path) and os.path.isdir(uti_path):
            return uti_path
        path = os.path.dirname(path)
    return None

def get_top_and_lowest_energy(protein):
    # Get the current working directory
    current_path = os.getcwd()
    input_files=find_test_inputs_directory(current_path)
    target_energy=obenergy(f"{input_files}/{protein}/{protein}_target",scoring_function="UFF")[-1]
    receptor_energy=obenergy(f"{input_files}/{protein}/{protein}_receptor",scoring_function="UFF")[-1]
    lowest_energy=(float(target_energy)+float(receptor_energy))/2
    higest_energy=(float(target_energy)+float(receptor_energy))*2
    print("Lowest energy:", lowest_energy)
    print("Highest energy:",higest_energy)
    return lowest_energy,higest_energy

def energy_to_list(proteins,scoring_function,input_path):
    # third filtered proteins
    os.chdir(input_path)
    output_dictionary_energy = {}
    for pro in proteins:
        energy = (obenergy(pro,
                           scoring_function=scoring_function)[-1])
        output_dictionary_energy[pro] = float(energy)
    return output_dictionary_energy

def energy_filtering_new(protein,protein_list,input_path):
    lowest_energy,top_energy=get_top_and_lowest_energy(protein)
    if len(protein_list) > 350:
        num_chunks = 350  # Toplam işlem sayısı
    else:
        num_chunks = int(len(protein_list))
    chunk_size = len(protein_list) // num_chunks  # Her bir işlem parçasının boyutu
    chunks = [protein_list[i:i + chunk_size] for i in range(0, len(protein_list), chunk_size)]

    # Paralel işlemler için ProcessPoolExecutor kullanımı
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = {}
        futures["obenergy"] = {executor.submit(energy_to_list, chunk, "UFF", input_path) for chunk in chunks}

        scores = {}
        for func_name, future_set in futures.items():
            scores[func_name] = [future.result() for future in future_set]

        # Birleştirilmiş sonuçları kullanarak işlevlerin sonuçlarını oluşturun
    energy_filtered_1 = {}
    for score_dict in scores["obenergy"]:
        energy_filtered_1.update(score_dict)

    #print(sorted_energy)
    output_list=[]
    energy_out={}
    for key, value in energy_filtered_1.items():
        if top_energy > value > lowest_energy:
            output_list.append(key)
            energy_out[key]=value
    return output_list,energy_out

########################################################################################################################

def find_uti_directory(start_path):
    path = start_path
    while path != os.path.dirname(path):  # Kök dizine ulaşana kadar devam et
        uti_path = os.path.join(path, "uti")
        if os.path.exists(uti_path) and os.path.isdir(uti_path):
            return uti_path
        path = os.path.dirname(path)
    return None

def pizsa_stability_filtration(pizsa_input):
    # Dosyayı aç ve satırları oku
    filtered_list=[]
    for pizsa_file in pizsa_input:
        protein_name=pizsa_file.split("._scores")[0]
        with open(pizsa_file, 'r') as file:
            lines = file.readlines()

            # Her bir satırı kontrol et
            for line in lines:
                # 'Binding Prediction' satırını bul
                if line.strip().startswith('Binding Prediction'):
                    # 'Binding Prediction' satırının değerini al
                    binding_prediction = line.split(':')[1].strip()
                    if binding_prediction == "Unstable association (NON-BINDER)":
                        filtered_list.append(protein_name)
    return filtered_list

def pizsa_score(first_filtered,input_path):
    current_path = os.getcwd()
    output_dictionary = {}
    for protein in first_filtered:
        uti_directory = find_uti_directory(current_path)
        os.chdir(uti_directory)
        try:
            process = subprocess.run(
                ["conda", "run", "-n", "pizsa", "python2", "ranking_PIZSA_cluster_rank.py", "--protein",
                 f"{protein}.pdb", "--input_path",
                 input_path],
                capture_output=True, text=True)
            output = process.stdout
            """print("outputoutputoutput")
            print(process)"""
            z_score = float(output.split()[-1])  # Split on spaces and take the last item
            pizsa_txt=output.split()[-4].replace(")","")
            #print("z_score", z_score)
            output_dictionary[protein] = z_score
        except:
            z_score=float(-55555)
            output_dictionary[protein] =z_score

    return output_dictionary


def pizsa_filtering_with_threshold(protein_list,input_path,lowest):
    if len(protein_list) > 350:
        num_chunks = 350  # Toplam işlem sayısı
    else:
        num_chunks = len(protein_list)
    chunk_size = len(protein_list) // num_chunks  # Her bir işlem parçasının boyutu
    chunks = [protein_list[i:i + chunk_size] for i in range(0, len(protein_list), chunk_size)]

    # Paralel işlemler için ProcessPoolExecutor kullanımı
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = {}
        futures["pizsa"] = {executor.submit(pizsa_score, chunk, input_path) for chunk in chunks}
        scores = {}
        for func_name, future_set in futures.items():
            scores[func_name] = [future.result() for future in future_set]

        # Birleştirilmiş sonuçları kullanarak işlevlerin sonuçlarını oluşturun
    pizsa_score_filtered_1 = {}
    for score_dict in scores["pizsa"]:
        pizsa_score_filtered_1.update(score_dict)
    sorted_pizsa = dict(sorted(pizsa_score_filtered_1.items(), key=lambda item: item[1], reverse=True))
    pizsa_filtered_list = []
    pizsa_dic={}
    for key, value in pizsa_score_filtered_1.items():
        if float(value) > lowest:
            pizsa_filtered_list.append(key)
            pizsa_dic[key]=value
    return pizsa_filtered_list,pizsa_dic

########################################################################################################################
def vormqa(protein_list,input_path):
    os.chdir(input_path)
    with open("pdb.list", "w") as pdb_list_file:
        for protein in protein_list:
            pdb_list_file.write(f"{protein}.pdb" + "\n")

    process = subprocess.run(
        ["conda", "run", "-n", "pizsa", "python2", f"{FCC_path}/make_contacts.py", "-f", "pdb.list", "-n", "8"],
        capture_output=True,
        text=True
    )
    #print("stdout:", process.stdout)
    #print("stderr:", process.stderr)
    os.system(r"sed -e 's/pdb/contacts/' pdb.list | sed -e '/^$/d' > pdb.contacts")

    output_dictionary = {}
    for protein in protein_list:
        command = f'{voronota_path}/voronota-voromqa --input {protein}.pdb --contacts-query \'--no-same-chain --no-solvent\' --print-energy-of-contacts-selection'

        output = os.popen(command).read()

        # Use regular expression to extract the first number
        match = re.search(r'(\d+\.\d+)', output)

        voromqa_score = float(match.group(1))
        output_dictionary[protein] = voromqa_score

    return output_dictionary

def vormqa_filtering(protein_list, input_path,filtering_ratio):
    if len(protein_list) > 350:
        num_chunks = 350  # Toplam işlem sayısı
    else:
        num_chunks = len(protein_list)
    chunk_size = len(protein_list) // num_chunks  # Her bir işlem parçasının boyutu
    chunks = [protein_list[i:i + chunk_size] for i in range(0, len(protein_list), chunk_size)]

    # Paralel işlemler için ProcessPoolExecutor kullanımı
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = {}
        futures["vormqa"] = {executor.submit(vormqa, chunk, input_path) for chunk in chunks}

        scores = {}
        for func_name, future_set in futures.items():
            scores[func_name] = [future.result() for future in future_set]


    vormqa_list_1 = {}
    for score_dict in scores["vormqa"]:
        vormqa_list_1.update(score_dict)

    sorted_vorqma = dict(sorted(vormqa_list_1.items(), key=lambda item: item[1], reverse=True))

    end_index = int(len(sorted_vorqma) * (1 - filtering_ratio))
    selected = {k: v for i, (k, v) in enumerate(sorted_vorqma.items()) if i < end_index}
    return list(selected.keys()), selected
########################################################################################################################



