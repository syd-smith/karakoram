#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 12:25:36 2026

@author: u1301408
"""

import cartopy.crs as ccrs
import glob
import matplotlib.animation as animation
import matplotlib.patches as patches
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from wrf import (getvar, ALL_TIMES)

ctl_path = '/uufs/chpc.utah.edu/common/home/strong-group7/husile/karakoram/model_result/wrfout_ctl/2020/'
exp_path = '/uufs/chpc.utah.edu/common/home/strong-group7/husile/karakoram/model_result/wrfout_MODISImproved/2020/'


def pr_animation(fpath, save_name):

    collect_files = sorted(glob.glob(fpath + 'wrfout_d02*06-30*')) + sorted(glob.glob(fpath + 'wrfout_d02*07*'))
    open_files = [Dataset(file) for file in collect_files]
    
    # precipitable water -> kg m-2
    accumulated_rain  =  getvar(open_files, 'RAINNC', timeidx = ALL_TIMES)
    hourly_pr = accumulated_rain.diff(dim = 'Time')
    
    open_bounds = Dataset('/uufs/chpc.utah.edu/common/home/strong-group7/husile/karakoram/model_result/wrfout_MODISImproved/2016/wrfout_d03_2016-06-08_00:00:00')
    area_bounds = getvar(open_bounds, 'RAINNC', timeidx=ALL_TIMES)
    min_long = float(area_bounds['XLONG'].values.min())
    max_long = float(area_bounds['XLONG'].values.max())
    min_lat = float(area_bounds['XLAT'].values.min())
    max_lat = float(area_bounds['XLAT'].values.max())
    
    fig = plt.figure()
    ax = plt.axes(projection = ccrs.PlateCarree())
    ax.coastlines()
    levels = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
    cmap = plt.get_cmap('GnBu')
    
    width = max_long - min_long
    height = max_lat - min_lat

    rect = patches.Rectangle((min_long, min_lat), width, height, 
                             linewidth=2, edgecolor='red', facecolor='none', 
                             label='D03 Bounds', zorder=5)

    ax.add_patch(rect)
    
    ims = []
    for i in range(len(hourly_pr.Time)):
        # Grab the data for one hour
        data = hourly_pr.isel(Time=i)
        
        # Plot the contour
        im = ax.contourf(data.XLONG, data.XLAT, data, 
                         levels=levels, cmap=cmap, extend='max',
                         transform=ccrs.PlateCarree())
        
        # Add a timestamp label that changes every frame
        t_str = str(data.Time.values)[:16]
        txt = ax.text(0.02, 0.95, t_str, transform=ax.transAxes, 
                      fontsize=12, fontweight='bold', bbox=dict(facecolor='white', alpha=0.7))
        
        # List of items to "draw" for this frame
        # (The first frame creates the colorbar, but we only need one)
        if i == 0:
            plt.colorbar(im, ax=ax, label='Precipitation (mm/hr)', shrink=0.6)
            
        current_frame = [im, txt]
        ims.append(current_frame)  
    
    
    ani = animation.ArtistAnimation(fig, ims, interval=100, blit=True, repeat_delay=1000)
    ani.save(save_name, writer='ffmpeg', fps=10)
    plt.show()
    
    return ani

control = pr_animation(ctl_path, 'ctl_d02_precip_july_2020.mp4')
experiment = pr_animation(exp_path, 'exp_d02_precip_july_2020.mp4')

