#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 18:51:52 2023
Updated in Oct 2024
@author: Milan Salek

This script creates a plot that shows the main components of the 
radiation fluxes as measured by CERES project

The CERES data are from the web site:
https://ceres.larc.nasa.gov/data/#ebaftoa-level-3

Order EBAF-TOA, Level 3b data 
(Observed TOA all-sky and clear-sky fluxes; CERES-MODIS cloud properties. 
Clear-sky for cloud free areas of 1°x1° region.)

Some sources of algorithms for beginners in xarray and NetCDF:
https://earth-env-data-science.github.io/assignments/basic_xarray.html
https://fabienmaussion.info/climate_system/week_02/01_Lesson_NetCDF_Data.html

The Gistemp data are from https://data.giss.nasa.gov/gistemp/graphs_v4/graph_data/Monthly_Mean_Global_Surface_Temperature/graph.txt
You have to adjust/edit the dataset to the correct time range according to the CERES dataset
starting by 2000.21 (March 2000) and ending by the same month as the CERES data
See also # Defining the start and end date below!

The global average is computed by the same algorithm as the global average of
weighted temperature:
https://docs.xarray.dev/en/stable/examples/area_weighted_temperature.html

All the resulting values are 12-month rolling (moving) average!
"""

import xarray as xr
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
plt.rcParams['figure.figsize'] = (12, 6)
import netCDF4 as nc
from scipy import signal
import sys
from datetime import  datetime

# Defining the start and end date
start_date = "2000-03-01"
end_date = "2024-07-31"
avg_temp_Earth = 14.5  # rather arbitrary, not too important 

fpath_ceres = "CERES_EBAF-TOA_Ed4.2_Subset_200003-202407.nc"
ds = xr.open_dataset(fpath_ceres, engine='netcdf4')
print(ds.time)

# Reading the time series downloaded from GISTEMP NASA:
# https://data.giss.nasa.gov/gistemp/graphs_v4/graph_data/Monthly_Mean_Global_Surface_Temperature/graph.txt
# The file gistemp_2024-04_reduced.csv must contain the same time span as
# the CERES_EBAF*.nc file
gistemp_monthly_temperatures = []
with open("gistemp_2024-07_reduced.csv") as gistemp_p:
    next(gistemp_p)   # Skipping first line
    for line in gistemp_p:
#        line = gistemp_p.readline()
        array = line.split(',')
        monthly_temp_anomaly = float(array[2])
        monthly_temp = monthly_temp_anomaly + avg_temp_Earth        
        gistemp_monthly_temperatures.append(monthly_temp_anomaly)

time = pd.date_range(start=start_date, end=end_date, freq='MS')  # 'MS' for month start
time = time + pd.DateOffset(days=14)  # adjusting the time to the middle of the months
print (time)

# Create an xarray DataArray
gistemp_array = xr.DataArray(
    gistemp_monthly_temperatures,
    coords=[time],
    dims=["time"],
    name="monthly_time_series_data"
)

# Smoothing of th Gistemp data by 12-month moving average
gistemp_12months_smoothed = \
    gistemp_array.rolling(time = 12, 
            center = True).mean()    

# End of Gistemp data preparation
#------------------------------------------------
# CERES data preparation:
sw = ds.toa_sw_all_mon  # Short-wave outgoing radiation
solar_incmg = ds.solar_mon   # Solar incoming radiation
lw = ds.toa_lw_all_mon # Long-wave outgoing radiation
toa_net = ds.toa_net_all_mon

# Computation of weights that are proportional to cos(lat)
def cosine_function(x):
    return np.cos(np.radians(x))

result = xr.apply_ufunc(cosine_function, ds['lat'], keep_attrs=True)
ds['cosine_lat'] = result

# Computing the weighted average
weights = np.cos(np.deg2rad(ds.lat))
sw_weighted = sw.weighted(weights)
lw_weighted = lw.weighted(weights)
toa_net_weighted = toa_net.weighted(weights)
solar_incmg_weighted = solar_incmg.weighted(weights)

sw_weighted_mean = sw_weighted.mean(("lon", "lat")) # Computing of the areal (global) mean (short-wave outhoing) 
lw_weighted_mean = lw_weighted.mean(("lon", "lat")) # Computing of the areal (global) mean (long-wave outgoing) 
toa_net_weighted_mean = toa_net_weighted.mean(("lon", "lat")) # Computing of the areal (global) mean (radiation balance, difference) 
solar_incmg_weighted_mean = solar_incmg_weighted.mean(("lon", "lat")) # Computing of the areal (global) mean (solar incoming) 

# 12-month moving average of the global mean component of the radiation values: 
sw_weighted_mean_smoothed = \
    sw_weighted_mean.rolling(time = 12, 
                             center = True).mean()    
   
lw_weighted_mean_smoothed = \
    lw_weighted_mean.rolling(time = 12, 
                             center = True).mean()    

toa_net_weighted_mean_smoothed = \
    toa_net_weighted_mean.rolling(time = 12, 
                             center = True).mean()

solar_incmg_weighted_mean_smoothed = \
    solar_incmg_weighted_mean.rolling(time = 12, 
                             center = True).mean()    
    
plt.rcParams['figure.figsize'] = [10, 10]
fig, axs = plt.subplots(4, sharex=True)

ylim_range = 2.6
ylim_offset = 0.2

# English version:

axs[0].plot(sw_weighted_mean.time, sw_weighted_mean_smoothed, color="blue", linewidth=3)
axs[0].set_ylim(float(sw_weighted_mean_smoothed.min()) - ylim_offset, 
            float(sw_weighted_mean_smoothed.min() + ylim_range))
axs[0].set_ylabel('Short-wave Outgoing')
axs[0].set_title('Short-wave outgoing (reflected) radiation, ToA, global mean (W/m^2)', color="blue")
axs[0].set_xlim(datetime(2000,1,1), datetime(2025,1,1))
axs[0].grid('on', which='major', axis='y', linestyle=':')
axs[0].grid('on', which='major', axis='x', linestyle=':')
axs[0].text(datetime(2002, 2, 1), 100.8, '''Radiation balance components from CERES data (NASA/NOAA) 
  at the top of the atmosphere (ToA) and GISTEMP temperatures.
             Twelve-month running means of global quantities''',
            style='italic',
            bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 10}, fontsize=14)
axs[0].text(datetime(2016, 1, 1), 99.5, 'Data NASA/NOAA, author: M. Šálek', style='italic',
        bbox={'facecolor': 'grey', 'alpha': 0.1, 'pad': 10}, fontsize=8)


axs[1].plot(lw_weighted_mean.time, lw_weighted_mean_smoothed, color="red",linewidth=3)
axs[1].set_ylim(float(lw_weighted_mean_smoothed.min())  - ylim_offset, 
            float(lw_weighted_mean_smoothed.min() + ylim_range))
axs[1].set_ylabel('Long-wave Outgoing')
axs[1].set_title('Long-wave outgoing radiation, ToA, global mean  (W/m^2)', color = 'red')
axs[1].grid('on', which='major', axis='y', linestyle=':')
axs[1].grid('on', which='major', axis='x', linestyle=':')
axs[1].text(datetime(2016, 1, 1), 240, 'Data NASA/NOAA, author: M. Šálek', style='italic',
        bbox={'facecolor': 'grey', 'alpha': 0.1, 'pad': 10}, fontsize=8)


axs[2].plot(toa_net_weighted_mean.time, toa_net_weighted_mean_smoothed,
            color="black", linewidth=3)

# Adding Gistemp temperature anomaly
axs[2].plot(toa_net_weighted_mean.time, gistemp_12months_smoothed,
            color="red", linewidth=2)
axs[2].set_ylim(float(toa_net_weighted_mean_smoothed.min())  - ylim_offset, 
            float(toa_net_weighted_mean_smoothed.min() + ylim_range))
axs[2].set_ylabel('TOA-NET')
axs[2].set_ylabel('ToA net and gistemp T')
axs[2].set_title('Radiation balance (net flux), ToA, global mean (W/m^2, black) and GISTEMP anomaly (°C, red)')
axs[2].grid('on', which='major', axis='y', linestyle=':')
axs[2].grid('on', which='major', axis='x', linestyle=':')
axs[2].text(datetime(2016, 1, 1), 0.3, 'Data NASA/NOAA, author: M. Šálek', style='italic',
        bbox={'facecolor': 'grey', 'alpha': 0.1, 'pad': 10}, fontsize=8)

axs[3].plot(solar_incmg_weighted_mean.time, solar_incmg_weighted_mean_smoothed.data, color="blue", linewidth=3)
axs[3].set_xlabel('Time')
axs[3].set_ylim(float(solar_incmg_weighted_mean_smoothed.min())  - ylim_offset, 
            float(solar_incmg_weighted_mean_smoothed.min() + ylim_range))
axs[3].set_ylabel('SW')
axs[3].set_title('Solar incoming radiation (W/m^2)', color="blue")
axs[3].set_xlim(datetime(2000,1,1), datetime(2025,1,1))
axs[3].grid('on', which='major', axis='y', linestyle=':')
axs[3].grid('on', which='major', axis='x', linestyle=':')
axs[3].text(datetime(2016, 1, 1), 341.7, 'Data NASA/NOAA, author: M. Šálek', style='italic',
        bbox={'facecolor': 'grey', 'alpha': 0.1, 'pad': 10}, fontsize=8)

fig.savefig("CERES_radiation_fluxes_globe_ENG.png")

sys.exit()

# Czech version:

axs[0].plot(sw_weighted_mean.time, sw_weighted_mean_smoothed, color="blue", linewidth=3)
#axs[0].set_xlabel('Time')
axs[0].set_ylim(float(sw_weighted_mean_smoothed.min()) - ylim_offset, 
            float(sw_weighted_mean_smoothed.min() + ylim_range))
axs[0].set_ylabel('SW')
axs[0].set_ylabel('Krátkovlnné odraž. záření')
#axs[0].set_title('''Data globálního záření projektu CERES (NASA/NOAA) \na teplotní řady GISTEMP, dvanáctiměsíční klouzavé průměry\n\n, 
axs[0].set_title('''Krátkovlnné odražené záření na hor. hranici atmosf (W/m^2)''')
axs[0].set_xlim(datetime(2000,1,1), datetime(2025,1,1))
axs[0].grid('on', which='major', axis='y', linestyle=':')
axs[0].grid('on', which='major', axis='x', linestyle=':')
axs[0].text(datetime(2000, 1, 1), 100.8, '''Měření bilance záření přístroji CERES (NASA/NOAA) a teplotní řada GISTEMP
                    Dvanáctiměsíční klouzavé průměry pro celou Zemi''',
            style='italic',
            bbox={'facecolor': 'grey', 'alpha': 0.5, 'pad': 10}, fontsize=14)
axs[0].text(datetime(2016, 1, 1), 99.5, 'Data NASA/NOAA, author: M. Šálek', style='italic',
        bbox={'facecolor': 'grey', 'alpha': 0.1, 'pad': 10}, fontsize=8)


axs[1].plot(lw_weighted_mean.time, lw_weighted_mean_smoothed, color="red",linewidth=3)
#axs[1].set_xlabel('Time')
axs[1].set_ylim(float(lw_weighted_mean_smoothed.min())  - ylim_offset, 
            float(lw_weighted_mean_smoothed.min() + ylim_range))
axs[1].set_ylabel('LW')
axs[1].set_ylabel('Dlouhovlnnné vyzařování')
axs[1].set_title('Long - wave upwelling, global mean', color = 'red')
axs[1].set_title('CERES, dlouhovlnné vyzařování na horní hranici atmosféry (W/m^2)', color = 'red')
axs[1].text(datetime(2000, 6, 1), 241.88, 'Podle hypotézy zvyšujícího se skleníkového efektu má tato křivka klesat', style='italic',
        bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 10}, fontsize=14)
axs[1].grid('on', which='major', axis='y', linestyle=':')
axs[1].grid('on', which='major', axis='x', linestyle=':')
axs[1].text(datetime(2016, 1, 1), 240, 'Data NASA/NOAA, author: M. Šálek', style='italic',
        bbox={'facecolor': 'grey', 'alpha': 0.1, 'pad': 10}, fontsize=8)


axs[2].plot(toa_net_weighted_mean.time, toa_net_weighted_mean_smoothed,
            color="black", linewidth=3)

# Adding Gistemp temperature anomaly
axs[2].plot(toa_net_weighted_mean.time, gistemp_12months_smoothed,
            color="red", linewidth=2)
#axs[2].set_xlabel('Time')
axs[2].set_ylim(float(toa_net_weighted_mean_smoothed.min())  - ylim_offset, 
            float(toa_net_weighted_mean_smoothed.min() + ylim_range))
axs[2].set_ylabel('TOA-NET')
axs[2].set_ylabel('Bilance záření, teplota')
axs[2].set_title('Net flux, TOA, global mean')
axs[2].set_title('Bilance záření na hor. hran. atmosf. (W/m^2, černě) a odchylka teploty dle GISTEMP (°C, červeně)')
axs[2].grid('on', which='major', axis='y', linestyle=':')
axs[2].grid('on', which='major', axis='x', linestyle=':')
axs[2].text(datetime(2016, 1, 1), 0.3, 'Data NASA/NOAA, author: M. Šálek', style='italic',
        bbox={'facecolor': 'grey', 'alpha': 0.1, 'pad': 10}, fontsize=8)

axs[3].plot(solar_incmg_weighted_mean.time, solar_incmg_weighted_mean_smoothed.data, color="blue", linewidth=3)
axs[3].set_xlabel('Time')
axs[3].set_ylim(float(solar_incmg_weighted_mean_smoothed.min())  - ylim_offset, 
            float(solar_incmg_weighted_mean_smoothed.min() + ylim_range))
axs[3].set_ylabel('SW')
axs[3].set_ylabel('Přijímané krátkovlnné záření ')
axs[3].set_title('CERES measurement, solar incoming', color="blue")
axs[3].set_title('CERES, přijímané krátkovlnné záření od Slunce (W/m^2)', color="blue")
axs[3].set_xlim(datetime(2000,1,1), datetime(2025,1,1))
axs[3].grid('on', which='major', axis='y', linestyle=':')
axs[3].grid('on', which='major', axis='x', linestyle=':')
axs[3].text(datetime(2016, 1, 1), 341.7, 'Data NASA/NOAA, author: M. Šálek', style='italic',
        bbox={'facecolor': 'grey', 'alpha': 0.1, 'pad': 10}, fontsize=8)

fig.savefig("CERES_radiation_fluxes_globe_CZ.png")




