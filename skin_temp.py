# %%
"""
Author: Sydney Smith
Created: June 3, 2026
"""

import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from netCDF4 import Dataset
import numpy as np
from pathlib import Path
import sys
from wrf import (getvar, latlon_coords)


# ==================================
# - Establish Relative File Path - 
# ==================================

try:
    current_file_directory = Path(__file__).resolve().parent
except NameError:
    current_file_directory = Path().resolve().parent
# parent_directory = current_file_directory.parent
sys.path.append(str(current_file_directory))

import from_savanna.nclcmaps as cmap
from preprocessor.data_access import get_data


noise = get_data('TSK', 5, 'd01', 'noise')
ctl = get_data('TSK', 5, 'd01', 'ctl')

# calculate anomaly skin temperature over d01 between noise experiment and ctl
# note we want this to be close to 0
domain_anom = noise[0] - ctl[0]
print(domain_anom)

# select only the first time stamp of the month
noise_data = noise[0].isel(Time = 0)
ctl_data = ctl[0].isel(Time = 0)
anom_data = domain_anom.isel(Time = 0)


# make three panel figure of skin tmeperature in noise experiment and control
fig, axs = plt.subplots(nrows = 1, ncols = 3, subplot_kw = {'projection': ccrs.PlateCarree()}, figsize = (25, 20))

# flatten axs to a 1D array of 2 elements
flat_axs = axs.flatten()

# create labels for each plot
labels = ['a.', 'b.', 'c.']

# open WRF to pull elevation data
open_elevation = Dataset(current_file_directory / 'wrfout' / f'wrfout_ctl' / '2016' / f'wrfout_d01_2016-05-01_00:00:00')
# use getvar to access data
elevation_data = getvar(open_elevation, 'ter')
elevation_lat, elevation_lon = latlon_coords(elevation_data)

# set elevation bands to appear every 500 meters
elevation_levels = np.arange(0, 8500, 500)

for ax, exp, label in zip(flat_axs, [noise_data, ctl_data, anom_data], labels):
    if label == 'c.':
        chosen_cmap = plt.get_cmap(cmap.cmap('MPL_bwr'), 10) # Red for warm anomalies, Blue for cold
        # Center the anomaly levels cleanly around 0
        # max_abs = max(float(abs(exp.min())), float(abs(exp.max())))
        # levels = np.linspace(-max_abs, max_abs, 11) 
        levels = np.linspace(-0.01, 0.01, 11)
        cbar_format = '%.2f'
        norm = mcolors.TwoSlopeNorm(vmin = -0.01, vcenter = 0, vmax = 0.01)
        
    else:
        chosen_cmap = plt.get_cmap(cmap.cmap('sunshine_9lev'), 10)  # Standard clean yellow-orange-red
        # levels = np.linspace(250, 290, 11)
        levels = np.linspace(222, 282, 11)
        cbar_format = '%.0f'
        
        # define norm values to standardize colorbar to map colors
        norm = mcolors.BoundaryNorm(levels, ncolors = 10) # 15 bins means 14 edges
   
    ax.contourf(exp['XLONG'], exp['XLAT'], exp.values, levels = levels, norm = norm, cmap = chosen_cmap, transform = ccrs.PlateCarree(), extend = 'both')

    # add label to panel
    ax.text(0.02, 0.92, label, fontsize = 15, fontweight = 'bold', transform = ax.transAxes)

    # plot elevation bands on each panel
    elevations = ax.contour(elevation_lon, elevation_lat, elevation_data, levels = elevation_levels, colors = 'black', linewidths = 0.25, transform = ccrs.PlateCarree())
    
    # create a scalar mappable as a standin for contour so the colorbar remains standardized across different maps
    sm = mpl.cm.ScalarMappable(norm = norm, cmap = chosen_cmap)
    sm.set_array([])
    
    # colorbar function passed using the scalar mappable 
    cbar = fig.colorbar(sm, ax = ax, orientation = 'horizontal', extend = 'both', ticks = levels, boundaries = levels, aspect = 30, location = 'bottom', pad = 0.015, format = cbar_format)
    cbar.ax.tick_params(labelsize = 10)
    

# check out land classifications of red dots
domain = 'd01'
fpath = current_file_directory.parent.parent / 'husile' / 'karakoram' / 'model_result' / 'wps_k2' / f'geo_em.{domain}.nc'

wrf_file = Dataset(fpath)
land_use = getvar(wrf_file, 'LU_INDEX')

fig = plt.figure(figsize = (10, 8))
ax = plt.axes(projection = ccrs.PlateCarree())

# see land use categories here: https://www2.mmm.ucar.edu/wrf/users/docs/user_guide_V3/user_guide_V3.3/users_guide_chap3.htm#_Land_Use_and
# create custom color map for land use categories
colors = plt.cm.jet(np.linspace(0, 1, 25))
# call out number 16 as black for water bodies
colors[20] = [0, 0, 0, 1] 
cmap = mcolors.ListedColormap(colors)
norm = mcolors.BoundaryNorm(np.arange(0, 26, 1), cmap.N)

# plot land use data
contour = ax.contourf(land_use['XLONG_M'], land_use['XLAT_M'], land_use.values, levels = np.arange(0, 25, 1), cmap = cmap, norm = norm, transform = ccrs.PlateCarree())
plt.colorbar(contour, ax = ax,label = 'Land Use Classification', ticks = np.arange(0, 25, 1), orientation = 'horizontal', pad = 0.05)
plt.show()





# %%
