"""
Copyright (C) 2013 Neelesh Soni <neelesh.soni@alumni.iiserpune.ac.in>, <neelrocks4@gmail.com>
This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

'''
This routine do the binary search for x in the list nums
'''
def recBinSearch(x, nums, low, high):
    #Check if low value is greater than high value
    if low > high:
        return low
    #fing the mid point
    mid = (low + high) / 2
    
    item1 = nums[mid]
    item2 = nums[mid+1]
    #if range is found
    if ((item1 <= x) & (item2 > x)):
        return mid
    #Else do the binary search on left and right parts
    elif x < item1:
        return recBinSearch(x, nums, low, mid-1)
    else: 
        return recBinSearch(x, nums, mid+1, high)

'''
Call the binary search routine
'''
def bsearch(x, nums):
    return recBinSearch(x, nums, 0, len(nums)-1)


#below shows the test example
'''
from sortlist import qsort
list2=[[4,1],[56,6],[5,3],[76,7],[83,9],[2,12],[45,23],[67,4],[11,67],[10,7],[20,5],[98,9],[9,22],[6,7]]
list2=qsort(list2)
list2.append([10000000])
print list2
cord=25;rcord=30;lcord=20
print bsearch(7,list2)
print bsearch(99,list2)
'''
