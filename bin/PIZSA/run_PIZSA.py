#!/usr/bin/env python

"""
Copyright (C) 2015 Abhilesh Dhawanjewar <abhilesh7@gmail.com> and Ankit Roy <intellect.ankit@gmail.com>

This file is part of PIZSA.

PIZSA is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PIZSA is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with PIZSA. If not, see <http://www.gnu.org/licenses/>.
"""

__author__ = "Abhilesh Dhawanjewar, Ankit Roy and Neelesh Soni"
__copyright__ = "Copyright C 2015"
__license__ = "LGPL"
__version__ = "2.0"
__maintainer__ = "Abhilesh Dhawanjewar and Ankit Roy"
__email__ = "abhilesh7@gmail.com, intellect.ankit@gmail.com"


import scripts.src_main.main
import sys
import logging
import traceback
import time

start_time = time.time()

err_dir = "/home/yavuz/yavuz_proje/protac/deneme/PIZSA"

logging.basicConfig()

logger = logging.getLogger(__name__)

ch = logging.FileHandler(err_dir + "/pizsa_log.err", mode = "w", delay =True)

logger.addHandler(ch)

try:
	scripts.src_main.main.main()
except Exception:
	logger.exception("Stopping Execution.... Check traceback below")

end_time = time.time()

#print "Run Completed! Check the 'output' folder for results!"
#print "Time elapsed: ", round(end_time - start_time, 3), 's'
