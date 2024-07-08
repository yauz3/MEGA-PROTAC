"""
Copyright (C) 2013 Neelesh Soni <neelesh.soni@alumni.iiserpune.ac.in>, <neelrocks4@gmail.com>
This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine reads the pdb file given as input and then return the coordinates and corresponding atom identities.
'''
def parse_pdb(inf):
    
    #read the pdb file#
    lines=inf;
    
    #initializing the maximum variable with negative infinity and  minimum value as positive infinity
    minx=10000000.0
    maxx=-10000000.0
    miny=10000000.0
    maxy=-10000000.0
    minz=10000000.0
    maxz=-10000000.0
    
    #Initiallizing the atom_props variable
    atom_props=[]
    
    #Iterating through each line
    for l in lines:
        #if line starts with 'ATOM' identifier, then continue
        if (l[0:4]=='ATOM'):
            #Getting the coordinates from pdb files
            xc=float(l[30:38]);
            yc=float(l[38:46]);
            zc=float(l[46:54]);
            
            #Getting the corresponding atom properties. Atom name, residue name,number, chain id
            # @Abhilesh modified - Hydrogen filter 
            if l.strip()[-1] != 'H':
                atom_type=l[12:16].strip();
                res_name=l[16:20].strip();
                res_num=l[22:26].strip();
                chain_num=l[21:22].strip();
            
            #Getting the max, min of x direction
                if minx>xc:
                    minx=xc;
                elif maxx<xc:
                    maxx=xc;
            
            #Getting the max, min of y direction
                if miny>yc:
                    miny=yc;
                elif maxy<yc:
                    maxy=yc;
            #Getting the max, min of z direction
                if minz>zc:
                    minz=zc;
                elif maxz<zc:
                    maxz=zc;
            
            #Append the atom_coord with each atom properties
                atom_props.append([xc,yc,zc,atom_type,res_name,res_num,chain_num]);

    #for element in atom_props:
    #   print element
    
    #Return            
    return minx,maxx,miny,maxy,minz,maxz,atom_props


    
