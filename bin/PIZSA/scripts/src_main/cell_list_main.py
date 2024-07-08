"""
Copyright (C) 2013 Neelesh Soni <neelesh.soni@alumni.iiserpune.ac.in>, <neelrocks4@gmail.com>
This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

''' 
Main routine for Cell_list implementation for calculating pairwise distances
'''

# @Abhilesh modified - Flags for specifying proteins
def get_dist(input_pdb, cutoff, p1 = None, p2 = None):
    
    from readpdb import parse_pdb
    from mesh import create_mesh
    from assign_coord import assign_atomlist_to_mesh
    from assign_neighbor import assign_ngh
    from calculate_distance_matrix import calculate_distance_dm, calculate_distance_mm

    #Open the input file for parsing.#
    inf2=input_pdb;
    #parse_pdb function returns the details of all atoms in the pdb file with maximum and minimum value of coordinates in each direction. All these values gets stored in the 'param' list#
    param=parse_pdb(inf2);
    #close the file#
    #inf2.close();

    #Set the distance cutoff#
    #cutoff=8.0;

    #Set the coarse mesh cell length in angstrom#
    cmesh_len=cutoff;

    #function 'create_mesh' creates a 3D mesh such that all atoms should be contained within this mesh. Extra 'cutoff' length is added in each direction for boundary cells#
    cmesh=create_mesh(param,cmesh_len,cutoff);

    #getting the all atom coordinates and their identities like atom type, residue name, residue number and chain identifier. All this are stored int he 6th element of param list#
    atom_props=param[6];

    #Getting the number of total number cell in the mesh#
    totcell=cmesh[3]*cmesh[4]*cmesh[5];

    #initializing the cell_list variable. This stores the list of all atoms in each cell of the mesh#
    cell_list=[ [] for n in range(0,totcell) ]

    #Assign the atoms in each cell of the mesh and store them in cell_list#
    cell_list=assign_atomlist_to_mesh(cmesh,atom_props,cell_list);

    #get the neighbors of each cell in the mesh and store them in the neighbor dictionary "nghdict" #
    ngh_dict=assign_ngh(cmesh);

    #Calculate the dictance of all the atoms that are withion cutoff distance. This function stores the distances with atom identities in the outline string.#
    # @Abhilesh modified - If clause added
    if p1 == None and p2 == None:
        outline=calculate_distance_dm(cell_list,cmesh,ngh_dict,cutoff)
    else:
        outline=calculate_distance_mm(cell_list,cmesh,ngh_dict,cutoff,p1,p2)

    #Open the output file for writing the distances#
    dist_list = []
    #Write the outline variable to the file#
    for element in outline.split('\n')[:-1]:
        dist_list.append(element + '\n')

    return dist_list
