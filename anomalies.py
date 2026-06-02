"""
Author: Sydney Smith
Created: April 2, 2026
"""

# %%
import cartopy.crs as ccrs
import copy 
import glob
import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from netCDF4 import Dataset
import numpy as np
from pathlib import Path
import sys
from wrf import (ALL_TIMES, getvar, latlon_coords)
import xarray as xr


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
from data_access import five_yr_anom


# ==============
# - Constants -
# ==============


# ==============
# - Functions - 
# ==============

def color_bar(var_min, var_center, var_max, color, location, axs):
    """
    Function to simply create a colorbar for maps. Note that this colorbar
    won't be standardized using a scalar mappable to apply to multiple maps. 
    """
    
    ticks = np.linspace(var_min, var_max, num = 10)
    norm = mcolors.TwoSlopeNorm(vmin = var_min, vcenter = var_center, vmax = var_max)

    # create a scalar mappable as a standin for contour so the colorbar remains standardized across different maps
    sm = mpl.cm.ScalarMappable(norm = norm, cmap = color)
    sm.set_array([]) # makes sure no data is attached to the colorbar
    
    # specify the layout of the colorbar
    cbar = plt.colorbar(sm, ax = axs, orientation = 'vertical', pad = 0.015, aspect = 50, extend = 'both', ticks = ticks, location = location)
    cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))
    cbar.ax.tick_params(labelsize = 30)
    
    return cbar


def six_panel_fig(variable, month, domain, experiment, units, elevation = True, save_name = 'place holder', save = False):
    """
    Function used to create 
    """

    # define path for data call
    path = current_file_directory / 'preprocessor' / f'{variable}_analysis.nc'
    ds = xr.open_dataset(path)
    # calculate yearly anomalies
    anom_data = ds['exp'].sel(experiment = experiment, month = month, domain = domain) - ds['ctl'].sel(experiment = experiment, month = month, domain = domain)
    # pull five year average anomaly data
    anom = ds['anom'].sel(experiment = experiment, month = month, domain = domain)

    # set up figure
    fig, axs = plt.subplots(nrows = 2, ncols = 3, subplot_kw = {'projection': ccrs.PlateCarree()}, figsize = (15, 10))
    
    # create label names for each panel
    labels = ['2016', '2017', '2018', '2019', '2020', 'Mean']
    
    if elevation: 
        # open WRF to pull elevation data
        open_elevation = Dataset(current_file_directory / 'wrfout' / f'wrfout_{experiment}' / '2016' / f'wrfout_{domain}_2016-06-08_00:00:00')
        # use getvar to access data
        elevation_data = getvar(open_elevation, 'ter')
        elevation_lat, elevation_lon = latlon_coords(elevation_data)

        # set elevation bands to appear every 500 meters
        elevation_levels = np.arange(0, 8500, 500)

    if domain == 'd01':
        xvals = np.arange(72, 81.5, 0.5)
        yvals = np.arange(32.5, 39.5, 0.5)
    
    elif domain == 'd02':
        xvals = np.arange(75.2, 77.8, 0.25)
        yvals = np.arange(35, 37, 0.25)

    else:
        xvals = anom.XLONG[0, ::3].values
        yvals = anom.XLAT[::3, 0].values

    # define levels
    levels = np.linspace(anom.min(), anom.max(), 15) 
        
    # define norm values to standardize colorbar to map colors
    norm = mcolors.BoundaryNorm(levels, ncolors = plt.get_cmap(cmap.cmap('MPL_BrBG'), 14).N) # 15 bins means 14 edges
    # TODO: define how many bins there should be for the colorbar? should it be 14?

    # flatten axs to a 1D array of 6 elements
    flat_axs = axs.flatten()

    # add data to each panel
    for i, ax in enumerate(flat_axs[:-1]):
        # ax.set_extent([-143, -67.5, 20, 44.5])
        # TODO: define the right boundaries for d01 and d02

        # define calendar year
        year = 2016 + i
        # call data for that year
        data = anom_data.sel(year = year) 
        maps_row1 = ax.contourf(data['XLONG'], data['XLAT'], data.values, levels = levels, norm = norm, cmap = cmap.cmap('MPL_BrBG'), extend = 'both', transform = ccrs.PlateCarree())
    
        # add label to panel
        ax.text(0.02, 0.89, labels[i], fontsize = 30, fontweight = 'bold', transform = ax.transAxes)
        
        if elevation:
            # plot elevation bands on each panel
            elevations = ax.contour(elevation_lon, elevation_lat, elevation_data, levels = elevation_levels, colors = 'black', linewidths = 0.25, transform = ccrs.PlateCarree())
        
        # add x and y axis values
        ax.set_xticks(xvals, crs = ccrs.PlateCarree())
        ax.set_xticklabels(xvals, rotation = 45, ha = 'right')
        ax.set_yticks(yvals, crs = ccrs.PlateCarree())
        ax.set_yticklabels(yvals)

    # add anomaly data to last panel
    anom_map = flat_axs[-1].contourf(anom['XLONG'], anom['XLAT'], anom.values, levels = levels, norm = norm, cmap = cmap.cmap('MPL_BrBG'), extend = 'both', transform = ccrs.PlateCarree())
    
    # add label to five year average panel
    flat_axs[-1].text(0.02, 0.89, labels[-1], fontsize = 30, fontweight = 'bold', transform = flat_axs[-1].transAxes)
    
    if elevation:
        # add elevation bands to five year average panel
        elevations = flat_axs[-1].contour(elevation_lon, elevation_lat, elevation_data, levels = elevation_levels, colors = 'black', linewidths = 0.25, transform = ccrs.PlateCarree())
        
    # add x and y axis values
    flat_axs[-1].set_xticks(xvals, crs = ccrs.PlateCarree())
    flat_axs[-1].set_xticklabels(xvals, rotation = 45, ha = 'right')
    flat_axs[-1].set_yticks(yvals, crs = ccrs.PlateCarree())
    flat_axs[-1].set_yticklabels(yvals)
    
    # fine tuning control of colorbar size and placement
    axins = inset_axes(ax, width = '10%', height = '20%', loc = 'center right', bbox_to_anchor = (0.935, 0.08, 0.08, 1.5), bbox_transform = fig.transFigure, borderpad = 0)

    # create a scalar mappable as a standin for contour so the colorbar remains standardized across different maps
    sm = mpl.cm.ScalarMappable(norm = norm, cmap = cmap.cmap('MPL_BrBG'))
    sm.set_array([])
    
    #  colorbar function passed using the scalar mappable 
    cbar = fig.colorbar(sm, cax = axins, orientation = 'vertical', extend = 'both', ticks = levels, boundaries = levels, aspect = 50)
    cbar.ax.tick_params(labelsize = 20)
    cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter(units))

    if save:
        # all PNGs stored to anomaly_maps directory but ignored in Git
        save_path = current_file_directory.joinpath(save_name)
        plt.savefig(save_path, dpi = 400)
        
    plt.show()












# ctl = test_mfc('d02', 7, 'ctl')
# exp = test_mfc('d02', 7, 'MODISImproved')
# noise = test_mfc('d02', 7, 'noise')


# # termini experiment precip data
# anom_term, ctl_term_mean, exp_term_mean, ctl_term, exp_term = five_yr_anom('RAINNC', 7, 'd02', 'MODISImproved')
# anom_term_1, ctl_term_1_mean, exp_term_1_mean, ctl_term_1, exp_term_1 = five_yr_anom('RAINNC', 7, 'd01', 'MODISImproved')

# # noise experiment precip data
# anom_noise, ctl_noise_mean, exp_noise_mean, ctl_noise, exp_noise = five_yr_anom('RAINNC', 7, 'd02', 'noise')
# anom_noise_1, ctl_noise_1_mean, exp_noise_1_mean, ctl_noise_1, exp_noise_1 = five_yr_anom('RAINNC', 7, 'd01', 'noise')

# # %%
# # domain 2 four year (2016-2019) precip anomaly 
# d02_term_anomaly = plot_anom(anom_term, 'Anomaly of Four Year Precip Average', '%', 'MPL_BrBG', 'd02', elevation = True)
# d02_noise_anomaly = plot_anom(anom_noise, 'Anomaly of Four Year Precip Average', '%', 'MPL_BrBG', 'd02', elevation = True)
# # %%
# d01_term_anomaly = plot_anom(exp_term_1_mean - ctl_term_1_mean, 'Anomaly of Four Year Precip Average', 'mm / 6 hours', 'MPL_BrBG', 'd01', elevation = True)
# d01_noise_anomaly = plot_anom(exp_noise_1_mean - ctl_noise_1_mean, 'Anomaly of Four Year Precip Average', 'mm / 6 hours', 'MPL_BrBG', 'd01', elevation = True)


# # %%
# # domain 1 2016 precip anomaly
# precip_2016_anom = plot_anom((exp_noise_1[0]-ctl_noise_1[0]), '2016 Precipitation Anomaly', 'mm/ 6 hours', 'MPL_BrBG', 'd01', elevation = True)

# # domain 2 raw data
# precip_2016_ctl = plot_anom(ctl_noise[0], '2016 Precipitation - Control', 'mm/ 6 hours', 'MPL_PuBuGn', 'd02', elevation = True)
# precip_2016_exp = plot_anom(exp_noise[0], '2016 Precipitation - Noise Experiment', 'mm/ 6 hours', 'MPL_PuBuGn', 'd02', elevation = True)
# precip_2016_term = plot_anom(exp_term[0], '2016 Precipitation - Termini Experiment', 'mm/ 6 hours', 'MPL_PuBuGn', 'd02', elevation = True)

# # %%

# anom = (exp_d01[0]/ctl_d01[0])*100
# precip_2016_anom = plot_anom(anom, '2016 Precipitation Anomaly', '%', 'MPL_BrBG', 'd01', elevation = True)
# # precip_2016_ctl = plot_anom(ctl_d01[0], '2016 Precipitation - Control', 'mm/ 6 hours', 'MPL_PuBuGn', 'd01', elevation = True)
# # precip_2016_ctl = plot_anom(exp_d01[0], '2016 Precipitation - Experiment', 'mm/ 6 hours', 'MPL_PuBuGn', 'd01', elevation = True)



# %%


# def main(precip = True, mfc = True):
    
#     if precip:
#         # five years of raw data for control and experiment respectively
#         ctl_02 = data_access('RAINNC', 7, 'd01', 'ctl')
#         exp_02 = data_access('RAINNC', 7, 'd01', 'MODISImproved')
        
#         # five year anom
#         # pr_anom, ctl_mean, exp_mean = five_yr_anom('RAINNC', 7, 'd01')
#         # pr_anomaly = plot_anom(pr_anom, 'Anomaly of Five Year Precip Average', '%', 'MPL_BrBG')
#         # five_year = plot_anom(ctl_mean, 'Five Year Precipitation  Average', 'mm/ 6 hr', 'MPL_PuBuGn')
        
        
#         for index, year in enumerate(range(2016, 2021)):
#             # yearly_map = plot_anom(ctl_02[index], f'{year} Precipitation', 'mm/ 6 hrs', 'MPL_PuBuGn')
            
#             # find the anomaly from experiment to control for that given year
#             pr_anom = exp_02[index] - ctl_02[index]
#             pr_map = plot_anom(pr_anom, f'{year} Precip Anomalies', 'mm / 6hr', 'MPL_BrBG')
            
#             # find anomaly from given year to five year control average
#             # anom_from_ctl_mean = ctl_02[index] / ctl_mean * 100
#             # from_ctl_mean_map = plot_anom(anom_from_ctl_mean, f'{year} from Five Year Average - Control Precip', '%', 'MPL_BrBG')
            
#             # find anomaly from given year to five year experiment average
#             # anom_from_exp_mean = exp_02[index] / exp_mean * 100
#             # from_exp_mean_map = plot_anom(anom_from_exp_mean, f'{year} from Five Year Average - Experiemtn Precip', '%', 'MPL_BrBG')
    
#     if mfc:
#         ctl = test_mfc('d01', 7, 'ctl')
#         exp = test_mfc('d01', 7, 'MODISImproved')
        
#         combine_ctl = xr.concat(ctl, dim = 'Years')
#         ctl_coords = {'XLAT': combine_ctl.XLAT.isel(Years = 0), 'XLONG': combine_ctl.XLONG.isel(Years = 0)}
#         ctl_mean =  combine_ctl.mean(dim = 'Years').assign_coords(ctl_coords)
#         # five_year_mfc = plot_anom(ctl_mean, 'Five Year MFC Average', 'g kg-1 s-1', 'posneg_2')
        
        
#         # combine_exp = xr.concat(exp, dim = 'Years')
#         # exp_coords = {'XLAT': combine_exp.XLAT.isel(Years = 0), 'XLONG': combine_exp.XLONG.isel(Years = 0)}
#         # exp_mean =  combine_exp.mean(dim = 'Years').assign_coords(exp_coords)
        
#         for index, year in enumerate(range(2016, 2021)):
            
#             mfc_yearly = plot_anom(exp[index], f'{year} Moisture Flux Convergence at 350 hPa', 'g kg-1 s-1', 'posneg_2')
            
#             # anom_data = (exp[index] - ctl[index])
#             # mfc_anom = plot_anom(anom_data, f'{year} MFC Anomaly at 350 hPa', 'as difference', 'posneg_2')
    
# # if __name__ == '__main__':
# #     main(precip = False)

# # %%

# %%
