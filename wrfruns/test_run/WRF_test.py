#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 11:48:27 2026

@author: u1301408
"""

import glob
from netCDF4 import Dataset
from pathlib import Path
import sys
from wrf import (getvar, ALL_TIMES)


# ==================================
# - Establish Relative File Path - 
# ==================================

current_file_directory = Path(__file__).resolve().parent
# parent_directory = current_file_directory.parent
sys.path.append(str(current_file_directory))


# =========================
# - Check Files for Data - 
# =========================

# create a list with all control experiment files and open them
collect_files = sorted(glob.glob(current_file_directory + 'wrfout*'))
open_files = [Dataset(file) for file in collect_files]

# find hourly precipitation vales (data essentially creating a staircase so you have to difference values from the previous timestamp to current)
accumulated_pr = getvar(open_files, 'RAINNC', timeidx = ALL_TIMES)
hourly_pr = accumulated_pr.diff(dim = 'Time') / 6

# call temperature data
temp =  getvar(open_files, 'T2', timeidx = ALL_TIMES)









