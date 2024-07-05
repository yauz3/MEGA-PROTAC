# -*- coding: utf-8 -*-
# 28/06/2024
# Author: Sadettin Y. Ugurlu & David McDonald

"""
https://github.com/akiyamalab/MEGADOCK


 Available options:
  -o filename    : set the output filename (default to "$R-$L.out")
  -O             : output docking detail files
  -N integer     : set the number of output predictions (default to 2000)
  -t integer     : set the number of predictions per each rotation (default to 1)
  -F integer     : set the number of FFT point (default to none)
  -v float       : set the voxel size (default to 1.2)
  -D             : set the 6 deg. (54000 angles) of rotational sampling
                   (default to none, 15 deg. (3600 angles) of rotational sampling)
  -r integer     : set the number of rotational sampling angles
                   (54000: 54000 angles, 1: 1 angles, 24: 24 angles, default to 3600 angles)
  -e float       : set the electrostatics term ratio (default to 1.0)
  -d float       : set the hydrophobic term ratio (default to 1.0)
  -a float       : set the rPSC receptor core penalty (default to -45.0)
  -b float       : set the rPSC ligand core penalty (default to 1.0)
  -f 1/2/3       : set function
                   (1: rPSC, 2: rPSC+Elec, 3: rPSC+Elec+RDE, default to 3)
  -h             : show this message


"""

import os
import subprocess
import sys
import shutil
import time
from pymol import cmd
import random
from tqdm import tqdm

# Setting the seed for the random number generator to ensure reproducibility
# The seed value (42) is arbitrary but consistent, so the sequence of random numbers
# generated will be the same each time this code is run
random.seed(42)

def megadock(
        target,
        ligand,
        MEGA_DOCK_PATH,
        output_number,
        ouput_file_path,
        #FFT_point,
        voxel,
        degree,
        rotation,
        #electrostatic,
        #hydrophobic,
        receptor_pel,
        ligand_pel,
        score

        ):
    only_target_name=target.split("/")[-1].split(".")[0]
    only_ligand_name = ligand.split("/")[-1].split(".")[0]
    print("Docking has been started!")
    try:
        os.mkdir(f"{ouput_file_path}/{only_target_name}")
    except:
        print("There is already output file in the dictionary.")
    print("########################################################")
    megadock_out=ouput_file_path
    ouput_file_path=f"{ouput_file_path}/{only_target_name}"
    os.chdir(MEGA_DOCK_PATH)
    os.system(f"./megadock-gpu -R {target} -L {ligand} -N {output_number} -o {ouput_file_path}/{only_target_name}-{only_ligand_name}.out "
              f" -r {rotation} -D {degree} -a {receptor_pel} -b {ligand_pel} -v {voxel} -f {score}")
    print("Docking is finished!")

    shutil.copy(f"{ligand}", f"{ouput_file_path}/{only_ligand_name}.pdb")
    try:
        shutil.copy(f"{target}", f"{ouput_file_path}/{only_target_name}.pdb")
    except:
        # it will work for ligand + protein-protein complex docking
        #print("buraaa")
        #print(f"{megadock_out}/{only_target_name.split('_')[0]}_rotate/{only_target_name}.pdb", f"{ouput_file_path}/{only_target_name}.pdb")
        shutil.copy(f"{megadock_out}/{only_target_name.split('_')[0]}_rotate/{str(only_target_name).replace('_rotate','')}.pdb", f"{ouput_file_path}/{only_target_name}.pdb")
    print("Ligand-protein structure have been forming using MEGA DOCK")
    time.sleep(3)
    for i in tqdm(range(1, output_number), desc="Complexes construction", unit="  Output PDB file"):
        os.chdir(MEGA_DOCK_PATH)
        os.system(f"./decoygen {ouput_file_path}/ligand_{i}.pdb {ouput_file_path}/{only_ligand_name}.pdb {ouput_file_path}/{only_target_name}-{only_ligand_name}.out {i}")


if __name__ == "__main__":
    megadock(
        target=target,
        ligand=ligand,
        MEGA_DOCK_PATH=MEGA_DOCK_PATH,
        output_number=output_number,
        ouput_file_path=ouput_file_path,
        #FFT_point=FFT_point,
        voxel=voxel,
        degree=degree,
        rotation=rotation,
        #electrostatic=electrostatic,
        #hydrophobic=hydrophobic,
        receptor_pel=receptor_pel,
        ligand_pel=ligand_pel,
        score=score
    )
