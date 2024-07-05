import numpy as np
from scipy.spatial.transform import Rotation as R
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
import random
import concurrent.futures
from tqdm import tqdm
from uti import protein_uti


# Set the random seed.
random.seed(42)


def find_test_inputs_directory(start_path):
    path = start_path
    while path != os.path.dirname(path):  # Kök dizine ulaşana kadar devam et
        uti_path = os.path.join(path, "test_inputs")
        if os.path.exists(uti_path) and os.path.isdir(uti_path):
            return uti_path
        path = os.path.dirname(path)
    return None



def translate_three_axes_with_ligands(pdb,
                                       target,
                                       ligand,
                                       chain_id,
                                       input_path, x_distances, y_distances, z_distances):

    result_files = []  # Yeni bir liste oluştur
    time.sleep(1)
    for x_distance in x_distances:
        for y_distance in y_distances:
            for z_distance in z_distances:
                    os.chdir(input_path)
                    cmd.reinitialize()
                    cmd.load(f"{pdb}.pdb")
                    cmd.load(f"/home/yavuz/yavuz_proje/protac/test_new/{target.split('_')[0]}/{target}.pdb")
                    cmd.load(f"/home/yavuz/yavuz_proje/protac/test_new/{ligand.split('_')[0]}/{ligand}.pdb")
                    cmd.align(ligand, pdb)
                    cmd.align(target, pdb)
                    cmd.remove(pdb)
                    cmd.translate([x_distance,y_distance,z_distance], ligand)
                    x_distance_name = str(x_distance).replace(".", "-")
                    y_distance_name = str(y_distance).replace(".", "-")
                    z_distance_name = str(z_distance).replace(".", "-")
                    previous_name=str(pdb).replace("_complex","")
                    result_file = f"{previous_name}_x_{x_distance_name}_y_{y_distance_name}_z_{z_distance_name}_complex.pdb"
                    cmd.save(result_file)
                    cmd.remove(target)
                    cmd.save(f"{previous_name}_x_{x_distance_name}_y_{y_distance_name}_z_{z_distance_name}_ligand.pdb")
                    cmd.remove("all")
                    result_files.append(result_file)
    return result_files


def rotate_poi_three_axes_with_ligands(pdb,
                                       target,
                                       ligand,
                                       chain_id,
                                       input_path, x_angles, y_angles, z_angles, center_rotation=True):


    result_files = []  # Yeni bir liste oluştur

    for x_angle in x_angles:
        for y_angle in y_angles:
            for z_angle in z_angles:
                cmd.remove("het")
                if center_rotation:
                    os.chdir(input_path)
                    cmd.reinitialize()
                    cmd.load(f"{pdb}.pdb")
                    cmd.load(f"/home/yavuz/yavuz_proje/protac/test_new/{target.split('_')[0]}/{target}.pdb")
                    cmd.load(f"/home/yavuz/yavuz_proje/protac/test_new/{ligand.split('_')[0]}/{ligand}.pdb")
                    cmd.align(ligand, pdb)
                    cmd.align(target, pdb)
                    cmd.remove(pdb)
                    center = cmd.centerofmass(selection=ligand)

                    cmd.rotate("x", angle=x_angle, origin=center, selection=ligand)
                    cmd.rotate("y", angle=y_angle, origin=center, selection=ligand)
                    cmd.rotate("z", angle=z_angle, origin=center, selection=ligand)


                    x_angle_name = str(x_angle).replace(".", "-")
                    y_angle_name = str(y_angle).replace(".", "-")
                    z_angle_name = str(z_angle).replace(".", "-")
                    file_name=str(pdb).replace("_complex","")
                    result_file = f"{file_name}_x_{x_angle_name}_y_{y_angle_name}_z_{z_angle_name}_complex.pdb"
                    print("result_file")
                    print(result_file)
                    cmd.save(result_file)
                    cmd.remove(target)
                    cmd.save(f"{file_name}_x_{x_angle_name}_y_{y_angle_name}_z_{z_angle_name}_ligand.pdb")
                    cmd.remove("all")
                    result_files.append(result_file)
                else:
                    os.chdir(input_path)
                    cmd.reinitialize()
                    cmd.load(f"{pdb}.pdb")
                    cmd.load(f"/home/yavuz/yavuz_proje/protac/test_new/{pdb}/{target}.pdb")
                    cmd.load(f"/home/yavuz/yavuz_proje/protac/test_new/{pdb}/{ligand}.pdb")
                    cmd.align(ligand, pdb)
                    cmd.align(target, pdb)
                    cmd.remove(pdb)
                    cmd.rotate("x", angle=x_angle, selection=ligand)
                    cmd.rotate("y", angle=y_angle, selection=ligand)
                    cmd.rotate("z", angle=z_angle, selection=ligand)

                    result_file = f"{pdb}_x_{x_angle}_y_{y_angle}_z_{z_angle}_no_center.pdb"
                    cmd.save(result_file)
                    cmd.remove("all")
                    result_files.append(result_file)

    return result_files




# the function for distance based filtration. The function calculates
# distance between ligand (small molecules, warhead and anchor) by using their mass center.
def distance_filtration(
        protein,
        input_path,
        complex,
        ligand_chain):
    # go for input files
    os.chdir(input_path)
    # Get the current working directory
    current_path = os.getcwd()
    # these are fixed points
    test_input_path=find_test_inputs_directory(current_path)
    warhead_coordinates = protein_uti.identify_centre_of_mass(f"{test_input_path}/{protein}/binder.pdb")

    # _complex files contains TWO PROTEINS and THEIR SMALL MOLECULES, including warhead and anchor.
    if "_complex.pdb" in complex:
        ligand_name = complex.replace("_complex.pdb", "")
    # referesh Pymol
    cmd.reinitialize()
    # _ligand files contains ONE PROTEIN (LIGAND PROTEIN) and ITS ligand.
    cmd.load(f"{ligand_name}_ligand.pdb")
    # Pymol gives the centerofmass
    anchor_coordinate = cmd.centerofmass(selection=f"het")
    # remove everything in Pymol in order to inhibit open several PDB files.
    cmd.remove("all")
    # calculate the distance between warhead (Fix point) and anchor (Predicted point)
    dist = math.sqrt(sum([(a - b) ** 2 for a, b in zip(warhead_coordinates, anchor_coordinate)]))
    return dist
