# %%

import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
from pathlib import Path
import sys
from wrf import (CoordPair, getvar, vertcross)


# ==================================
# - Establish Relative File Path - 
# ==================================

current_file_directory = Path(__file__).resolve().parent
# parent_directory = current_file_directory.parent
sys.path.append(str(current_file_directory))

from anomalies import data_access

# access all variables necessary for cross section
qvapor = data_access('QVAPOR', 7, 'd01', 'ctl')
ua = data_access('ua', 7, 'd01', 'ctl')
va = data_access('va', 7, 'd01', 'ctl')
ht = data_access('z', 7, 'd01', 'ctl')

# end points of the vertical cross section
start = CoordPair(lat = 35, lon = 72.5)
stop = CoordPair(lat = 33, lon = 76)

# rise and run of crooss section
dlat = stop.lat - start.lat
dlon = stop.lon - start.lon

# find angle of cross section
angle = np.arctan2(dlat, dlon)

# calculate what winds flow perpendicular to cross section
perpendicular_wind = -ua[0] * np.sin(angle) + va[0] * np.cos(angle)

# dot product to get water vapor transport 
WVT = qvapor[0] * perpendicular_wind

cross_section = vertcross(WVT, ht[0], wrfin = qvapor, start_point = start, end_point = stop, latlon = True, autolevels = 100)


# %%
