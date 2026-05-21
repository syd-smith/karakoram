# %%
import glob
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from netCDF4 import Dataset
import numpy as np
from pathlib import Path
import sys
from wrf import (CoordPair, getvar, to_np, vertcross)


# ==================================
# - Establish Relative File Path - 
# ==================================

current_file_directory = Path(__file__).resolve().parent
# parent_directory = current_file_directory.parent
sys.path.append(str(current_file_directory))

from anomalies import data_access
import from_savanna.nclcmaps as cmap

def load_variables(experiment):

    # access all variables necessary for cross section
    qvapor = data_access('QVAPOR', 7, 'd01', experiment)
    ua = data_access('ua', 7, 'd01', experiment)
    va = data_access('va', 7, 'd01', experiment)
    ht = data_access('z', 7, 'd01', experiment)

    return qvapor, ua, va, ht

noise_qvapor, noise_ua, noise_va, noise_ht = load_variables('noise')

# %%
def WVT(year, qvapor, ua, va, ht):
    q = qvapor[year] /(1 + qvapor[year]) # convert to specific humidity

    # end points of the vertical cross section
    start = CoordPair(lat = 35, lon = 72.5)
    stop = CoordPair(lat = 33, lon = 76)
        
    # rise and run of cross section
    dlat = stop.lat - start.lat
    dlon = (stop.lon - start.lon) * np.cos(dlat)
        
    # find angle of cross section
    angle = np.arctan2(dlat, dlon)

    # get total magnitude of wind
    magnitude = np.sqrt(ua**2 + va**2)
                        
    # calculate what winds flow perpendicular to cross section
    perpendicular_wind = ua[year] * np.sin(angle) + va[year] * np.cos(angle)
        
    # dot product to get water vapor transport 
    # WVT = np.sqrt((qvapor[0]*perp_ua)**2 + (qvapor[0]*perp_va)**2)
    wvt_dot = qvapor[year] * (perpendicular_wind)
        
    path = current_file_directory / f'wrfout_{experiment}/2016/'
    collect_files = sorted(glob.glob(str(path / f'wrfout_d01_2016-07*')))
    open_files = [Dataset(file) for file in collect_files]

    cross_section = vertcross(wvt_dot, ht[year], wrfin = open_files, start_point = start, end_point = stop, latlon = True, autolevels = 100)

    return cross_section

wvt = WVT(2, noise_qvapor, noise_ua, noise_va, noise_ht)


fig = plt.figure(figsize=(12, 6))
ax = plt.axes()

norm = mcolors.TwoSlopeNorm(vmin = -0.01, vcenter = 0, vmax = 0.05)
zvals = wvt['vertical'][0:39]
x_vals = np.arange(wvt.shape[1])

# Plot the cross-section data
# to_np converts the xarray back to a numpy array for plotting
wvt_contours = ax.contourf(x_vals, zvals, to_np(wvt.sel(vertical = slice(0, 12000))), cmap=cmap.cmap("MPL_BrBG"), norm = norm)

# xvals = np.arange(72, 81.5, 0.5)
# ax.set_xticks(xvals)
# ax.set_xticklabels(xvals, rotation = 45, ha = 'right')

yvals = np.arange(0, 12000, 1000)
ax.set_yticks(yvals)
ax.set_yticklabels(yvals)
ax.set_ylabel('Height (m)')

ax.set_xticks([])
ax.set_xticklabels([])

# Add a colorbar
plt.colorbar(wvt_contours, ax=ax, label='WVT (kg/m/s)')
plt.show()

# %%
