# -*- coding: utf-8 -*-
# 28/06/2024
# Author: Sadettin Y. Ugurlu & David McDonald

import MDAnalysis as mda
import numpy as np

def quality_score_mda(pdb_file):
    # MDAnalysis ile PDB dosyasını yükle
    u = mda.Universe(pdb_file)

    # Analiz için seçim yapılabilir, örneğin, proteinin Cα atomları
    selection = u.select_atoms('protein and name CA')

    # Cα atomlarının konum vektörlerini bir listeye kaydedin
    positions = []
    for atom in selection:
        positions.append(atom.position)

    # Vektör normlarını hesaplayın
    norms = [np.linalg.norm(position) for position in positions]

    # Kalite skorunu hesaplayın
    quality_score = np.mean(norms)

    return quality_score

