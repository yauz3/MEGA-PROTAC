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
from uti import uti_rotate_and_translate
from Bio import PDB
import random
import re
import os
import os
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
from tqdm import tqdm

""""""
FCC_path="/home/yavuz/fcc/scripts"
voronota_path="/home/yavuz/voronota_1.22.3149"
""""""



def search_space_translate_with_ligand(protein,protein_list,input_path):
    # Get the current working directory
    current_path = os.getcwd()
    parent_path = os.path.dirname(current_path)
    ligand_chain = protein_uti.get_chain_list(f"{parent_path}/{protein}_target_input/{protein}_target_input_1_ligand")[0]
    """x_distances = list(np.arange(-6, 6, 1.5))
    y_distances = list(np.arange(-6, 6, 1.5))
    z_distances = list(np.arange(-6, 6, 1.5))"""
    # Define the precision
    precision = 10

    # Generate the distances
    x_distances = [round(x, precision) for x in np.arange(-4.5, 4.5 + 1e-10, 1.5)]
    y_distances = [round(y, precision) for y in np.arange(-4.5, 4.5 + 1e-10, 1.5)]
    z_distances = [round(z, precision) for z in np.arange(-4.5, 4.5 + 1e-10, 1.5)]
    target=f"{protein}_target_input"
    ligand=f"{protein}_receptor_input"
    os.chdir(input_path)
    #protein_list = [file.replace('.pdb', '') for file in protein_list]
    time.sleep(2)
    for pro in tqdm(protein_list, desc="Search space is enlarging", unit="pdb"):
        """uti_rotate_and_translate.prepare_with_ligands(pro.replace(".pdb",""),
                                                      target,
                                                      ligand,
                                                      input_path,
                                                      )"""
        uti_rotate_and_translate.translate_three_axes_with_ligands(pro.replace(".pdb",""),
                                                                      target,
                                                                      ligand,
                                                                      ligand_chain,
                                                                      input_path,
                                                                      x_distances,
                                                                      y_distances,
                                                                      z_distances,
                                                                      )


def rotational_grid_search(protein,protein_list,input_path):
    # take ligand chain ID
    current_path = os.getcwd()
    parent_path = os.path.dirname(current_path)
    chain_id = protein_uti.get_chain_list(f"{parent_path}/{protein}_target_input/{protein}_target_input_1_ligand")[
        0]
    precision=10
    target=f"{protein}_target_input"
    ligand=f"{protein}_receptor_input"
    x_angles = [round(x, precision) for x in np.arange(-15, 15 + 1e-10, 5)]
    y_angles = [round(y, precision) for y in np.arange(-15, 15 + 1e-10, 5)]
    z_angles = [round(z, precision) for z in np.arange(-15, 15 + 1e-10, 5)]
    print("angles")
    print(x_angles)
    print(y_angles)
    print(z_angles)
    """
    x_angles = list(range(-25, 25, 5))
    y_angles = list(range(-25, 25, 5))
    z_angles = list(range(-25, 25, 5))
    """
    """
    x_angles = list(range(-25, 25, 3))
    y_angles = list(range(-25, 25, 3))
    z_angles = list(range(-25, 25, 3))
    """
    time.sleep(5)
    for pro in tqdm(protein_list, desc="Search space is enlarging", unit="pdb"):
        uti_rotate_and_translate.rotate_poi_three_axes_with_ligands(pro.replace(".pdb",""),
                                                                                   target,
                                                                                   ligand,
                                                                                   chain_id,
                                                                                   input_path, x_angles, y_angles,
                                                                                   z_angles, center_rotation=True)