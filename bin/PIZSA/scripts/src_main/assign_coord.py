"""
Copyright (C) 2013 Neelesh Soni <neelesh.soni@alumni.iiserpune.ac.in>, <neelrocks4@gmail.com>
This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine assigns each atom in the 3D mesh. This is done through binary search on the x,y and z arrays of mesh points
'''

from binary_search import bsearch
def assign_atomlist_to_mesh(mesh,atom_props,cell_lst):
    
    #get the total number of mesh points in each direction#
    xn=mesh[3];
    yn=mesh[4];
    zn=mesh[5];
    
    #storing coordinates in x-y-z major. ie. z coordinates first then y then x
    for c in atom_props:
        #Do a binary search of current x,y,z coordinates in the corresponding x,y,z array of points#
        ix=bsearch(c[0],mesh[0])
        iy=bsearch(c[1],mesh[1])
        iz=bsearch(c[2],mesh[2])
        
        #z-y-x major
        index=zn*yn*ix + iy*zn +iz;
        
        #append the mesh cell with coordinate & its properties
        cell_lst[index].append(c);
        #Return cell list
    return cell_lst
