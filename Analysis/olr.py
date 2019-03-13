from netCDF4 import Dataset
import wrf
import xarray as xr
import numpy as np

import matplotlib
import matplotlib.cm as mpl_cm

import matplotlib.pyplot as plt

import cartopy.crs as crs
import cartopy.feature as cfe

###################################
## Load in WRF data
###################################
## 2_Nisg80_ThompsonDefault/
## 3_Nisg80_ThompsonAeroClim/
## 4_Nisg80_Thompson_naCCN0408_naCCN1100/
## 5_Archer_Default_AeroClim/
## 6_Archer_INITpl100e6/
## 7_Archer_INITpl100/
## 8_Archer_INITpl100_DRIVERpl100/
## 9_Archer_DRIVER_NWFA1D_100e6/
## 10_Archer_DRIVER_NWFA1D_100/
## 11_Archer_DRIVER_NWFA1D_100e3/
## 12_Archer_DRIVER_NWFA1D_x2/	# FAILED JOB - ERROR
## 13_Archer_DRIVER_NWFA1D_x05/
## 14_Archer_DRIVER_NWFA1D_150e3/
## 15_Archer_DRIVER_NWFA1D_150e3_K1/
## 16_Archer_DRIVER_NWFA1D_100e3_K1/

file_dir1 = '3_Nisg80_ThompsonAeroClim/'
file_dir2 = '10_Archer_DRIVER_NWFA1D_100/'

root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MAC/WRF_V4.0.1/RUNS/'

time_index = 32

nc1 = Dataset(root_dir+file_dir1+'wrfout_d02_2015-11-27_00:00:00')
olr1 = wrf.getvar(nc1, 'OLR', timeidx=time_index)

nc2 = Dataset(root_dir+file_dir2+'wrfout_d02_2015-11-27_00:00:00')
olr2 = wrf.getvar(nc2, 'OLR', timeidx=time_index)

## Quick Plot to check all is well
# olr.plot()

## Get the latitude and longitude points
lats, lons = wrf.latlon_coords(olr1)

## Get the cartopy mapping object
cart_proj = wrf.get_cartopy(olr1)

###################################
###################################
## WRF data processing
###################################
###################################

###################################
#####	FILE #1
###################################
theta1 = wrf.getvar(nc1, 'T', timeidx=time_index) + 300 # potential temperature in K
theta1.name = 'Potential temperature, K'

pressure1 = wrf.getvar(nc1, 'P', timeidx=time_index) + wrf.getvar(nc1, 'PB', timeidx=time_index)   # pressure in Pa
pressure1.name = 'Air pressure, Pa'

tempvar = float(287.05)/float(1005)
tempvar0 = (pressure1/100000)**tempvar
temperature1 = tempvar0 * theta1
temperature1.name = 'Air Temperature, K'

###################################
#####	FILE #2
###################################
theta2 = wrf.getvar(nc2, 'T', timeidx=time_index) + 300 # potential temperature in K
theta2.name = 'Potential temperature, K'

pressure2 = wrf.getvar(nc2, 'P', timeidx=time_index) + wrf.getvar(nc2, 'PB', timeidx=time_index)   # pressure in Pa
pressure2.name = 'Air pressure, Pa'

tempvar = float(287.05)/float(1005)
tempvar0 = (pressure2/100000)**tempvar
temperature2 = tempvar0 * theta2
temperature2.name = 'Air Temperature, K'

###################################
## MAP
###################################

data1 = wrf.to_np(olr1[:,:])
data2 = wrf.to_np(olr2[:,:])

# Create a figure
fig = plt.figure(figsize=(8,4))

# Set the GeoAxes to the projection used by WRF
ax = fig.add_axes([0.1,0.1,0.4,0.8], projection=cart_proj)	# left, bottom, width, height
# ax = plt.axes(projection=cart_proj)

# Add coastlines
ax.coastlines('50m', linewidth=0.8)
ax.add_feature(cfe.NaturalEarthFeature('physical', 'antarctic_ice_shelves_lines',
                                       '50m', linewidth=1.0, edgecolor='k', facecolor='none') )

# Plot contours
plt.contourf(wrf.to_np(lons), wrf.to_np(lats), data1, 10,
                transform=crs.PlateCarree(), cmap = mpl_cm.Reds)

# Add a color bar
cbar = plt.colorbar(ax=ax, shrink=.62)
cbar.set_label(olr1.name[-5:])

# Set the map limits.  Not really necessary, but used for demonstration.
# ax.set_xlim(wrf.cartopy_xlim(olr1))
# ax.set_ylim(wrf.cartopy_ylim(olr1))

# Add the gridlines
ax.gridlines(color="black", linestyle="dotted")

plt.title(olr1.name+'\n'+str(olr1.Time.values))


# Set the GeoAxes to the projection used by WRF
ax = fig.add_axes([0.55,0.1,0.4,0.8], projection=cart_proj)	# left, bottom, width, height
# ax = plt.axes(projection=cart_proj)

# Add coastlines
ax.coastlines('50m', linewidth=0.8)
ax.add_feature(cfe.NaturalEarthFeature('physical', 'antarctic_ice_shelves_lines',
                                       '50m', linewidth=1.0, edgecolor='k', facecolor='none') )

# Plot contours
plt.contourf(wrf.to_np(lons), wrf.to_np(lats), data2, 10,
                transform=crs.PlateCarree(), cmap = mpl_cm.Reds)

# Add a color bar
cbar = plt.colorbar(ax=ax, shrink=.62)
cbar.set_label(olr2.name[-5:])

# Set the map limits.  Not really necessary, but used for demonstration.
# ax.set_xlim(wrf.cartopy_xlim(olr2))
# ax.set_ylim(wrf.cartopy_ylim(olr2))

# Add the gridlines
ax.gridlines(color="black", linestyle="dotted")

plt.title(olr2.name+'\n'+str(olr2.Time.values))

plt.show()

