1. USAGE
---------

1. To run the software package from the command line using default parameters, run the run_PIZSA.py script as:
	python run_PIZSA.py XXXX.pdb

2. The default distance threshold for defining an interaction is 4.0 Å. This can be changed to either 4.0, 6.0 or 8.0 Å by using the flag -d followed by 4.0, 6.0 or 8.0
	eg. python run_PIZSA.py XXXX.pdb -d 6.0


2. DEPENDENCIES
----------------

* PIZSA requires Python 2, Python2.7.15 and above but not Python3.x.
* PIZSA requires NumPy and SciPy.


3. SOURCE CODE
---------------

All the scripts for the package can be found in the 'scripts' directory.


4. LICENSE
-----------

Please refer to COPYING.txt for the full license or COPYING_LESSER.txt for a condensed version.
