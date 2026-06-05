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
from typing import Literal
from wrf import (ALL_TIMES, getvar, latlon_coords, to_np)
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
from preprocessor.data_access import five_yr_anom


# ==============
# - Functions - 
# ==============

def six_panel_fig(
        variable: str, 
        month: Literal[5, 6, 7], 
        domain: Literal['d01', 'd02', 'd03'], 
        experiment: str, 
        cbar_label: str, 
        elevation: bool = True, 
        save_name: str = 'place holder.png', 
        save: bool = False
    ) -> plt.Figure:
        
    """
    Function used to create a six panel figure for to show monthly averages 
    of a variable over the five simulation plus a five year average. 
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
    # labels = ['2016', '2017', '2018', '2019', '2020', 'Mean']
    labels = ['a.', 'b.', 'c.', 'd.', 'e.', 'f.']
    
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
        maps = ax.contourf(data['XLONG'], data['XLAT'], data.values, levels = levels, norm = norm, cmap = cmap.cmap('MPL_BrBG'), extend = 'both', transform = ccrs.PlateCarree())
    
        # add label to panel
        ax.text(0.02, 0.89, labels[i], fontsize = 15, fontweight = 'bold', transform = ax.transAxes)
        
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
    cbar.ax.tick_params(labelsize = 20) # size of numbers on colorbar
    cbar.set_label(cbar_label, size = 25, weight = 'bold', labelpad = 15)
    # cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter())

    if save:
        # all PNGs stored to anomaly_maps directory but ignored in Git
        save_path = current_file_directory / 'figures' / f'{save_name}.png'
        plt.savefig(save_path, dpi = 400)
        
    plt.show()


def six_panel_transect(
        month: Literal[5, 6, 7], 
        experiment: str, 
        cbar_label: str = 'kg kg⁻¹ m s⁻¹', 
        save_name: str = 'place holder.png', 
        save: bool = False
    ) -> plt.Figure:
   
    """
    Function used to create a six panel figure of a vertical transect for 
    monthly average water vappor transport over the five simulation plus a 
    five year average. 
    """
   
    # define path for data call
    path = current_file_directory / 'preprocessor' / 'wvt_analysis.nc'
    ds = xr.open_dataset(path)
    # calculate yearly anomalies
    anom_data = ds['exp'].sel(experiment = experiment, month = month) - ds['ctl'].sel(experiment = experiment, month = month)
    # pull five year average anomaly data
    anom = ds['anom'].sel(experiment = experiment, month = month)

    # set up figure
    fig, axs = plt.subplots(nrows = 2, ncols = 3, figsize = (15, 10))

    # create label names for each panel
    # labels = ['2016', '2017', '2018', '2019', '2020', 'Mean']
    labels = ['a.', 'b.', 'c.', 'd.', 'e.', 'f.']

    # define levels
    levels = np.linspace(float(anom.min()), float(anom.max()), 15) 
        
    # define norm values to standardize colorbar to map colors
    norm = mcolors.BoundaryNorm(levels, ncolors = plt.get_cmap(cmap.cmap('MPL_BrBG'), 14).N) # 15 bins means 14 edges
    # TODO: define how many bins there should be for the colorbar? should it be 14?

    # flatten axs to a 1D array of 6 elements
    flat_axs = axs.flatten()

    # add data to each panel
    for i, ax in enumerate(flat_axs[:-1]):

        # ax.set_extent([-143, -67.5, 20, 44.5])
        # TODO: define the right boundaries for d01 and d02
        
        # set nan values to show as black
        ax.set_facecolor('black')

        # define calendar year
        year = 2016 + i
        # call data for that year
        data = anom_data.sel(year = year) 

        # restrict cross section to show up to 12000 m
        data_slice = data.sel(vertical = slice(0, 12000))
        zvals = data_slice['vertical']
        x_vals = np.arange(data_slice.shape[1])
        # map wvt data
        maps = ax.contourf(x_vals, zvals, to_np(data_slice), levels = levels, norm = norm, cmap = cmap.cmap('MPL_BrBG'), extend = 'both')
    
        # add label to panel
        ax.text(0.02, 0.89, labels[i], fontsize = 15, fontweight = 'bold', transform = ax.transAxes)
  
        # add y values that reflect the height in meters
        yvals = np.arange(0, 12000, 1000)
        ax.set_yticks(yvals)
        ax.set_yticklabels(yvals)
        ax.set_ylabel('Height (m)')

        # remove x axis labels
        ax.set_xticks([])
        ax.set_xticklabels([])

    ### repeat same process with five year average ###
    # restrict cross section to show up to 12000 m
    anom_slice = anom.sel(vertical = slice(0, 12000))
    zvals = anom_slice['vertical']
    x_vals = np.arange(anom_slice.shape[1])

    # add anomaly data to last panel
    anom_map = flat_axs[-1].contourf(x_vals, zvals, to_np(anom_slice), levels = levels, norm = norm, cmap = cmap.cmap('MPL_BrBG'), extend = 'both')
    
    # set nan values to black
    flat_axs[-1].set_facecolor('black')

    # add label to five year average panel
    flat_axs[-1].text(0.02, 0.89, labels[-1], fontsize = 15, fontweight = 'bold', transform = flat_axs[-1].transAxes)
    
    # add y values that reflect the height in meters
    yvals = np.arange(0, 12000, 1000)
    flat_axs[-1].set_yticks(yvals)
    flat_axs[-1].set_yticklabels(yvals)
    flat_axs[-1].set_ylabel('Height (m)')

    # remove x axis labels
    flat_axs[-1].set_xticks([])
    flat_axs[-1].set_xticklabels([])

    # fine tuning control of colorbar size and placement
    axins = inset_axes(ax, width = '2%', height = '75%', loc = 'center right', bbox_transform = fig.transFigure, borderpad = 5)

    # create a scalar mappable as a standin for contour so the colorbar remains standardized across different maps
    sm = mpl.cm.ScalarMappable(norm = norm, cmap = cmap.cmap('MPL_BrBG'))
    sm.set_array([])
    
    #  colorbar function passed using the scalar mappable 
    cbar = fig.colorbar(sm, cax = axins, orientation = 'vertical', extend = 'both', ticks = levels, boundaries = levels, aspect = 50)
    cbar.ax.tick_params(labelsize = 20) # size of numbers on colorbar
    cbar.set_label(cbar_label, size = 25, weight = 'bold', labelpad = 15)
    # cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter(cbar_label))

    if save:
        # all PNGs stored to anomaly_maps directory but ignored in Git
        save_path = current_file_directory / 'figures' / f'{save_name}.png'
        plt.savefig(save_path, dpi = 400)
        
    plt.show()

    return fig


six_panel_transect(5, 'noise', 'WVT (kg/m/s)', save = False)


# %%

def main():
    for month in range(5, 8):
        for experiment in ['noise', 'MODISImproved']:
            # six_panel_fig('RAINNC', month, 'd01', experiment, 'mm / 6 hours', save_name = f'pr_{month}_{experiment}_d01', save = True)
            # six_panel_fig('RAINNC', month, 'd02', experiment, 'mm / 6 hours', save_name = f'pr_{month}_{experiment}_d02', save = True)
            six_panel_transect(month, experiment, 'WVT (kg/m/s)', save_name = f'wvt_{month}_{experiment}', save = True)


if __name__ == '__main__':
    main()


# %%

experiment = 'noise'
month = 5

# define path for data call
path = current_file_directory / 'preprocessor' / 'wvt_analysis.nc'
ds = xr.open_dataset(path)
# calculate yearly anomalies
anom_data = ds['exp'].sel(experiment = experiment, month = month) - ds['ctl'].sel(experiment = experiment, month = month)
# pull five year average anomaly data
anom = ds['anom'].sel(experiment = experiment, month = month)
# %%
