# -*- coding: utf-8 -*-
# 28/06/2024
# Author: Sadettin Y. Ugurlu & David McDonald


import os
import datetime
import math
import statistics
import subprocess
import sys
import textwrap
import timeit
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from pymol.cgo import *
from scipy.optimize import minimize
import os
from collections import OrderedDict
import math
import numpy as np
from pymol.cgo import *
from scipy.optimize import minimize
from pymol import cmd
import pandas as pd
import sys
import argparse
import timeit
import statistics
import textwrap
import datetime
import subprocess
import csv
import time
import scoria
import math
import numpy as np
from pymol.cgo import *
from pathlib import Path
from Bio.PDB import PDBParser
from scipy.optimize import minimize
from pymol.cgo import *
import pandas as pd
import sys
import argparse
import timeit
import statistics
import textwrap
import datetime
import os.path
import time
import csv
import os
import shutil
import numpy
import scoria
import numpy as np
import mdtraj as md
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBParser, PDBIO, Select
import glob
from Bio import PDB



def load_scoria_molecule(
        mol_path: str,
        verbose: bool = False,
):
    """Use the `scoria` library to build a molecule model from the 3D structure file.
    Parameters
    ----------
    mol_path : str
        The path to the 3D file. Will be converted to PDB format, if it is not already.
    verbose : bool, optional
        Flag to print updates to a console, by default False
    Returns
    -------
    scoria.Molecule
        A `scoria` molecule object constructed from the given file
    """

    if verbose:
        print("Building scoria molecule from", mol_path)

    # convert file if necessary
    if not mol_path.endswith(".pdb"):
        stem, ext = os.path.splitext(mol_path)
        ext = ext.replace(".", "")  # remove .
        mol_path = obabel_convert(
            input_format=ext,
            input_filename=mol_path,
            output_format="pdb",
            output_filename=stem + ".pdb")

    return scoria.Molecule(mol_path)


def identify_centre_of_mass(
        mol_path: str = None,
        mol: scoria.Molecule = None,
        precision: int = 3,
        geometric: bool = True,
        verbose: bool = False,
):
    """Use scoria to compute the center of mass / geometric center for a PDB file or pre-constructed scoria Molecule.
    Parameters
    ----------
    mol_path : str, optional
        Path of 3D molecule file, by default None
    mol : scoria.Molecule, optional
        Pre-constructed scoria Molecule, by default None
    precision : int, optional
        Rounding precision, by default 3
    geometric : bool, optional
        Flag to compute geometric mean, rather than center of mass, by default True
    verbose : bool, optional
        Flag to print updates to the console, by default False
    Returns
    -------
    tuple
        Computed x, y, z co-ordinates of the center
    """

    if verbose:
        print("Using scoria to compute center of molecule")

    if mol is None:  # construct molecule
        assert mol_path is not None and isinstance(mol_path, str)
        mol = load_scoria_molecule(mol_path)

    if not geometric:
        if verbose:
            print("Calculating mass center")
        try:
            center_x, center_y, center_z = mol.get_center_of_mass().data
        except Exception as e:
            print("Mass center failed")
            print("Exception was", e)
            geometric = True

    if geometric:
        if verbose:
            print("Calculating geometric center", )
        try:
            center_x, center_y, center_z = mol.get_geometric_center()
        except Exception as e:
            print("geometric center calculation error", e)
            return None

    # round using precision
    if isinstance(precision, int):
        center_x = round(center_x, precision)
        center_y = round(center_y, precision)
        center_z = round(center_z, precision)
    #print(center_x, center_y, center_z)
    return center_x, center_y, center_z

def get_chain_list(target_pdb):
    parser = PDBParser()
    structure_1 = parser.get_structure(f"{target_pdb}", f"{target_pdb}.pdb")
    # For each chain
    chain_list=[]
    for chain_1 in structure_1.get_chains():
        if " " != str(chain_1.id):
            chain_list.append (chain_1.id)
    print("Chain_list",chain_list)
    return chain_list







