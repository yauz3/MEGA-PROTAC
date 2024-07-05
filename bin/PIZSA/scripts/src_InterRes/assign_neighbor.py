"""
Copyright (C) 2013 Neelesh Soni <neelesh.soni@alumni.iiserpune.ac.in>, <neelrocks4@gmail.com>
This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA.  If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine finds the neighbor cells for each cell in the mesh. Neighbor cells are stored in the nghdict variable
'''
def assign_ngh(mesh):
    #get the total number of points in each direction
    xn=mesh[3];
    yn=mesh[4];
    zn=mesh[5];
    
    #initializing the neighbor dictionary
    nghdict={};
    
    #Calculating the total number of cells
    totcell=xn*yn*zn;
    
    #Iterating through total number of cells
    for lstidx in range(0,totcell):
        
        #calculate the x,y,z index using the mesh indexes. This is possible because the mesh indexes are z then y then x major. Here we are unwinding the mesh in x,y and z direction
        r=lstidx%(yn*zn);
        #get mesh x index
        mix=lstidx/(yn*zn);
        #get mesh z index
        miz=r%zn;
        #get mesh y index
        miy=r/zn;
        
        #initializing the temp neighbor list variable
        temp_nghlst=[];
        
        #All permuations/neighbors are taken by adding -1,0,-1 values in each direction. here for each cell in the mesh, its neighbors are recorded in nghdict variable
        for idsx in -1,0,1:
            for idsy in -1,0,1:
                for idsz in -1,0,1:
                    #Get index in each dimension
                    xid,yid,zid=(mix+idsx),(miy+idsy),(miz+idsz);
                    #Check if each index is less than the maximum number of cells in that direction and if each index is greater than -1
                    if((xid < xn)&(xid > -1)):
                        if((yid < yn)&(yid > -1)):
                            if((zid < zn)&(zid > -1)):
                                #Recalculate the mesh index
                                index=zn*yn*xid + yid*zn +zid;
                                #Append the mesh index in the nighbor list 
                                temp_nghlst.append(index)
        
        #Add the temp neighbor list in the dictionary
        nghdict[lstidx]=temp_nghlst;
        
    return nghdict
