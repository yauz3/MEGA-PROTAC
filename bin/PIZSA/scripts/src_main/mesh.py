"""
Copyright (C) 2013 Neelesh Soni <neelesh.soni@alumni.iiserpune.ac.in>, <neelrocks4@gmail.com>
This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine create the 3D mesh of points. Mesh is created so that all the atoms in the pdb file will be inside the mesh
'''
import numpy as np
def create_mesh(param,loc,cutoff):
    
    #Geting the maximum and minimum values of each coordinates from the param argument#
    minx=param[0]
    maxx=param[1]
    miny=param[2]
    maxy=param[3]
    minz=param[4]
    maxz=param[5]
    
    #Add the cutoff to the maximum value and subtract from the minimum value. This is done to extend the mesh dimension by cutoff distance in each direction.#
    minx-=cutoff
    maxx+=cutoff
    miny-=cutoff
    maxy+=cutoff
    minz-=cutoff
    maxz+=cutoff
    
    #loc is the length of cell. Here loc is equal to the cutoff value. This can ve a variable#
    #Getting the total number of cells in each dimension#
    xn=int((maxx-minx)/loc);
    yn=int((maxy-miny)/loc);
    zn=int((maxz-minz)/loc);
    
    #Setting the origin as minimum value in each dimension#
    origin=[minx,miny,minz]
    
    #Generating the equidistant points in x, y and z direction#
    xarray = np.linspace(minx,maxx,xn)
    yarray = np.linspace(miny,maxy,yn)
    zarray = np.linspace(minz,maxz,zn)
    
    #Setting up the mesh variable#
    mesh=[xarray,yarray,zarray,xn,yn,zn,origin,loc]
    
    #Print the MEsh Statistics#
    #print "\nMesh with following no of mesh points in (x,y,z) created: ","(",xn,yn,zn,")"
    #print "Mesh Origin: ", "(",origin[0],",",origin[1],",",origin[2],")"
    #print "Length of square cell/cutoff: ", loc
    #print "Total no of cells: ",xn*yn*zn,"\n"
    
    #Return the Mesh#
    return mesh
