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
def WVT(year, qvapor, ua, va, ht, experiment):

    # all tests done with ua and va values of 10
    # westerly wind -> 5.68 good
    # easterly wind -> -5.68 good
    # southerly wind -> 8.23. good
    # northerly wind -> -8.23 good
    # southwesterly wind -> 13.91 good 
        # also good becuase souwesterly wind should have a higher value than southerly or westerly wind respectively because it is more head on with transect
    # northeasterly wind -> -13.91 good 

    # end points of the vertical cross section
    start = CoordPair(lat = 35, lon = 72.5)
    stop = CoordPair(lat = 33, lon = 76)
    # adjust lon to be same units as lat since lon spacing changes as you move away from the equator
    # use angle of average latitude as theta for cos
    # convert degrees to radians
    avg_lat = (start.lat + stop.lat) / 2
    lat_rad = np.radians(avg_lat)
            
    # rise and run of cross section
    dlat = stop.lat - start.lat
    dlon = (stop.lon - start.lon) * np.cos(lat_rad) # adjusted by spacing of lat

    # find angle of cross section
    theta = np.arctan2(dlat, dlon)

    ## OLD way of calculating perpendicular wind that only returns positive values because the range of arccos is 0 to pi
    # # create a unit vector
    # t_cos = -np.cos(theta) # x or u value
    # t_sin = np.sin(theta) # y or v value

    # # get total magnitude of wind for dot product
    # M = np.sqrt(ua**2 + va**2)

    # # take the dot product of the wind direction and unit vector to get wind perpendicular to transect
    # dot = (ua/M) * t_cos + (va/M) * t_sin
    # # find alpha -> the angle of the wind to the transect
    # alpha = np.arccos(dot)

    # # get wind perpendicular to transect based on alpha
    # perpendicular_wind = M * np.sin(alpha)

    ## NEW wind calculations that preserves sign of wind with respect to trarnsect
    # if the unit vector originally is (-cos(theta), sin(theta)) it must be rotated 90 degrees clockwise to get the wind perpendicular to the transect
    # rotating clockwise makes westerly and southerly wind traveling perpendicular to the transect positive
    xofvector = np.abs(np.sin(theta))
    yofvector = np.abs(np.cos(theta))
    # taking the absolute value of the unit vectors makes it so the direction of the wind is preserved

    # take the dot product of the rotated unit vector and the wind
    perpendicular_wind = xofvector * ua[year] + yofvector * va[year]

    # convert qvapor to specific humidity
    q = qvapor[year] /(1 + qvapor[year]) 

    # WVT = specific humidity * vector of wind perpendicular to transect
    WVT = q * (perpendicular_wind)

    # pull openfiles data to feed to vertcross
    path = current_file_directory / f'wrfout_{experiment}/2016/'
    collect_files = sorted(glob.glob(str(path / f'wrfout_d01_2016-07*')))
    open_files = [Dataset(file) for file in collect_files]

    cross_section = vertcross(WVT, ht[year], wrfin = open_files, start_point = start, end_point = stop, latlon = True, autolevels = 100)

    return cross_section


for year, idx in zip(range(2016, 2021), range(0,4)):

    # calulate WVT for given year
    wvt = WVT(idx, noise_qvapor, noise_ua, noise_va, noise_ht, 'noise')

    # plot WVT transect
    fig = plt.figure(figsize = (12, 6))
    ax = plt.axes()

    norm = mcolors.TwoSlopeNorm(vmin = -0.01, vcenter = 0, vmax = 0.05)
    zvals = wvt['vertical'][0:39]
    x_vals = np.arange(wvt.shape[1])

    # Plot the cross-section data
    # to_np converts the xarray back to a numpy array for plotting
    wvt_contours = ax.contourf(x_vals, zvals, to_np(wvt.sel(vertical = slice(0, 12000))), cmap = cmap.cmap("MPL_BrBG"), norm = norm)

    # add y values that reflect the height in meters
    yvals = np.arange(0, 12000, 1000)
    ax.set_yticks(yvals)
    ax.set_yticklabels(yvals)
    ax.set_ylabel('Height (m)')

    # remove x axis labels
    ax.set_xticks([])
    ax.set_xticklabels([])

    # add colorbar
    plt.colorbar(wvt_contours, ax = ax, label = 'WVT (kg/m/s)')

    # add title
    plt.title(f'{year} WVT from 35N, 72.5E to 33N, 76E - Noise')
    plt.show()
    # plt.savefig()


# calculate five year average of transect for each experiment
# take anomalies

# %%
