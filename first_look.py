#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 28 17:42:51 2026

@author: u1301408
"""

import sys
sys.path.append("/uufs/chpc.utah.edu/common/home/u0660911/Documents/python/python_mydefs") 
from court_stats import my_block_bootstrap_diffmeans_test
from setuptools import command
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

import xarray as xr
import numpy as np
import glob 

from matplotlib.colors import LinearSegmentedColormap, CenteredNorm
from netCDF4 import Dataset
from wrf import (getvar, to_np, ALL_TIMES, interplevel, latlon_coords)


from scipy.stats import permutation_test
import scipy.stats as stats


# read data
# paths 
bdir = '/uufs/chpc.utah.edu/common/home/strong-group7/court/proposals/glacierwx/'
# ctl_dir = '/scratch/general/vast/u1215181/wrf_k2_ctl/'
# exp_dir = '/scratch/general/vast/u1215181/wrf_k2_MODISImproved/
ctl_dir = '/uufs/chpc.utah.edu/common/home/strong-group7/husile/karakoram/model_result/wrfout_ctl/2016/'
exp_dir = '/uufs/chpc.utah.edu/common/home/strong-group7/husile/karakoram/model_result/wrfout_MODISImproved/2016/'

# D03 month to analyze
ctl_files_d02 =  sorted(glob.glob(ctl_dir + "wrfout_d02*2016-06-30*")) + sorted(glob.glob(ctl_dir + "wrfout_d02*2016-07*"))
exp_files_d02 =  sorted(glob.glob(exp_dir + "wrfout_d02*2016-06-30*")) + sorted(glob.glob(exp_dir + "wrfout_d02*2016-07*"))


# D02 -----------------------------------
ctl_list_d02 = [Dataset(x) for x in ctl_files_d02]
exp_list_d02 = [Dataset(x) for x in exp_files_d02]

ctl_pr_d02 = getvar(ctl_list_d02, 'RAINNC', timeidx=ALL_TIMES) + getvar(ctl_list_d02, 'RAINC', timeidx=ALL_TIMES)
exp_pr_d02 = getvar(exp_list_d02, 'RAINNC', timeidx=ALL_TIMES) + getvar(exp_list_d02, 'RAINC', timeidx=ALL_TIMES)


def backward_diff_time(dataarray):
    forward_diff = dataarray.diff(dim='Time')
    first_timestep = dataarray.isel(Time=0) * 0
    backward_diff = xr.zeros_like(dataarray)
    backward_diff.isel(Time=0)[:] = first_timestep
    backward_diff.isel(Time=slice(1, None))[:] = forward_diff
    return backward_diff
ctl_pr = backward_diff_time(ctl_pr_d02)
exp_pr = backward_diff_time(exp_pr_d02)

ctl_pr = ctl_pr.sel(Time=ctl_pr.Time.dt.month.isin([6, 7]))
exp_pr = exp_pr.sel(Time=exp_pr.Time.dt.month.isin([6, 7]))


delt = exp_pr - ctl_pr
n = delt.shape[0]
nprimes = np.full((delt.shape[1],delt.shape[2]),np.nan)
for i in np.arange(0,delt.shape[1]):
    for j in np.arange(0,delt.shape[2]):
        r1_delt = np.corrcoef(delt[:-1,i,j].data, delt[1:,i,j].data)[0, 1]
        nprimes[i,j] = n*(1-r1_delt)/(1+r1_delt)

s2delt = np.var(delt,axis=0,ddof=1)
z = np.mean(delt,axis=0)/np.sqrt(s2delt/nprimes)
p_value = 2 * (1 - stats.norm.cdf(abs(z)))


open_bounds = Dataset('/uufs/chpc.utah.edu/common/home/strong-group7/husile/karakoram/model_result/wrfout_MODISImproved/2016/wrfout_d03_2016-06-08_00:00:00')
area_bounds = getvar(open_bounds, 'RAINNC', timeidx=ALL_TIMES)
min_long = float(area_bounds['XLONG'].values.min())
max_long = float(area_bounds['XLONG'].values.max())
min_lat = float(area_bounds['XLAT'].values.min())
max_lat = float(area_bounds['XLAT'].values.max())


# Plotting
plt.figure(figsize=(10, 8))
ax = plt.gca()
delt_bar = np.mean(delt,axis=0)*31 # monthly mean 
delt_bar.plot.contourf(x='XLONG', y='XLAT', levels=np.arange(-7,7.5, 0.5) , cmap='BrBG')

width = max_long - min_long
height = max_lat - min_lat

rect = patches.Rectangle((min_long, min_lat), width, height, 
                         linewidth=2, edgecolor='red', facecolor='none', 
                         label='D03 Bounds', zorder=5)

ax.add_patch(rect)

sig_lats, sig_lons = np.where(p_value <= 0.1)
plt.scatter(exp_pr_d02.XLONG.values[sig_lats, sig_lons], exp_pr_d02.XLAT.values[sig_lats, sig_lons], color='black', marker='.', s=8)

plt.xlabel('Longitude')
plt.ylabel('Latitude')        




