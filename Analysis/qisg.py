from netCDF4 import Dataset
import wrf
import xarray as xr
import numpy as np

import matplotlib
import matplotlib.cm as mpl_cm
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import cartopy.crs as crs
import cartopy.feature as cfe

###################################
## Load in WRF data
###################################
## 2_Nisg80_ThompsonDefault/
file_dir1 = '2_Nisg80_ThompsonDefault/'
file_dir2 = '3_Nisg80_ThompsonAeroClim/'

root_dir = '/data/mac/giyoung/MAC_WRFThompson/'

time_index = 32

nc1 = Dataset(root_dir+file_dir1+'wrfout_d02_2015-11-27_00:00:00')
qice1 = wrf.getvar(nc1, 'QICE', timeidx=time_index)
qsnow1 = wrf.getvar(nc1, 'QSNOW', timeidx=time_index)
qgraup1 = wrf.getvar(nc1, 'QGRAUP', timeidx=time_index)
qisg1 = qice1 + qsnow1 + qgraup1

nc2 = Dataset(root_dir+file_dir2+'wrfout_d02_2015-11-27_00:00:00')
qice2 = wrf.getvar(nc2, 'QICE', timeidx=time_index)
qsnow2 = wrf.getvar(nc2, 'QSNOW', timeidx=time_index)
qgraup2 = wrf.getvar(nc2, 'QGRAUP', timeidx=time_index)
qisg2 = qice2 + qsnow2 + qgraup2

## Quick Plot to check all is well
# qisg.plot()

## Get the latitude and longitude points
lats, lons = wrf.latlon_coords(qice1)

## Get the cartopy mapping object
cart_proj = wrf.get_cartopy(qice1)

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

rho1 = pressure1/(float(287.05) * temperature1)
rho1.name = 'Air density, kg m-3'

qisg1 = qisg1 * float(1e3)
qisg1.name = 'Qisg, gkg-1'

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

rho2 = pressure2/(float(287.05) * temperature2)
rho2.name = 'Air density, kg m-3'

qisg2 = qisg2 * float(1e3)
qisg2.name = 'Qisg, gkg-1'

###################################
## MAP
###################################

data1 = wrf.to_np(qisg1[16,:,:])
data2 = wrf.to_np(qisg2[16,:,:])

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
cbar.set_label(qisg1.name[-5:])

# Set the map limits.  Not really necessary, but used for demonstration.
# ax.set_xlim(wrf.cartopy_xlim(qisg1))
# ax.set_ylim(wrf.cartopy_ylim(qisg1))

# Add the gridlines
ax.gridlines(color="black", linestyle="dotted")

plt.title(qisg1.name+'\n'+str(qisg1.Time.values))


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
cbar.set_label(qisg2.name[-5:])

# Set the map limits.  Not really necessary, but used for demonstration.
# ax.set_xlim(wrf.cartopy_xlim(qisg2))
# ax.set_ylim(wrf.cartopy_ylim(qisg2))

# Add the gridlines
ax.gridlines(color="black", linestyle="dotted")

plt.title(qisg2.name+'\n'+str(qisg2.Time.values))

plt.show()

###################################
## VERTICAL PROFILE @ HALLEY
###################################

# Extract the model height 
z1 = wrf.getvar(nc1, "z")
z2 = wrf.getvar(nc2, "z")

plt.plot(np.squeeze(qisg1[:,137,183]),z1[:,137,183],label = 'Default')
plt.plot(np.squeeze(qisg2[:,137,183]),z2[:,137,183],label = 'AeroClim')
plt.ylim([0,2000])
plt.title(qisg1.name+'\n'+str(qisg1.Time.values))
plt.ylabel(z1.description)
plt.show()