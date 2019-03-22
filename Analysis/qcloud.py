from netCDF4 import Dataset
import wrf
import xarray as xr
import numpy as np

import matplotlib
import matplotlib.cm as mpl_cm
#from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import cartopy.crs as crs
import cartopy.feature as cfe

###################################
## Load in WRF data
###################################
###################################
## Morrison (GRL 2019)
###################################
## 30_DeMott_WATSAT_HM_noThresh_eta70_MYNN
## 31_DeMott_WATSAT_eta70_MYNN
## 32_DeMott_WATSAT_H2xHM_eta70_MYNN
## 33_DeMott_WATSAT_0xHM_eta70_MYNN
## 36_DeMott_WATSAT_2xHM_noThresh_eta70_MYNN
## 41_DeMott_WATSAT_2xHM_noThresh_WarmPIP_eta70_MYNN
## 42_ThompsonMP28
## 56_DeMott_WATSAT_10xHM_noThresh_eta70_MYNN
## 57_DeMott_WATSAT_5xHM_noThresh_eta70_MYNN

###################################
## Thompson AEROCLIM
###################################
## 2_Nisg80_ThompsonDefault/
## 3_Nisg80_ThompsonAeroClim/
## 4_Nisg80_Thompson_naCCN0408_naCCN1100/
## 5_Archer_Default_AeroClim
## 6_Archer_NWFApl100_AeroClim
## 7_Archer_INITpl100/
## 8_Archer_INITpl100_DRIVERpl100/
## 9_Archer_DRIVER_NWFA1D_100e6/
## 10_Archer_DRIVER_NWFA1D_100/
## 11_Archer_DRIVER_NWFA1D_100e3/       ## WRONG NAMELIST
## 12_Archer_DRIVER_NWFA1D_x2/       ## WRONG NAMELIST
## 13_Archer_DRIVER_NWFA1D_x05/       ## WRONG NAMELIST
## 14_Archer_DRIVER_NWFA1D_150e3/       ## WRONG NAMELIST
## 15_Archer_DRIVER_NWFA1D_150e3_K1/       ## WRONG NAMELIST
## 16_Archer_DRIVER_NWFA1D_100e3_K1/       ## WRONG NAMELIST
## 17_Archer_initialise_real_qnwfanow_x2/       ## WRONG NAMELIST
## 18_Archer_initialise_real_qnwfanow_K1_100e6/
## 19_Archer_initialise_real_qnwfanow_K1_x2/

file_dir1 = '42_ThompsonMP28/'
file_dir2 = '2_Nisg80_ThompsonDefault/'

root_dir1 = '/gws/nopw/j04/ncas_weather/gyoung/MAC/PWRF_V3.6.1/RUNS/'
root_dir2 = '/gws/nopw/j04/ncas_weather/gyoung/MAC/WRF_V4.0.1/RUNS/'

time_index = 32

nc1 = Dataset(root_dir1+file_dir1+'wrfout_d02_2015-11-27_00:00:00')
qcloud1 = wrf.getvar(nc1, 'QCLOUD', timeidx=time_index)

nc2 = Dataset(root_dir2+file_dir2+'wrfout_d02_2015-11-27_00:00:00')
qcloud2 = wrf.getvar(nc2, 'QCLOUD', timeidx=time_index)

## Quick Plot to check all is well
# qcloud.plot()

## Get the latitude and longitude points
lats, lons = wrf.latlon_coords(qcloud1)

## Get the cartopy mapping object
cart_proj = wrf.get_cartopy(qcloud1)

#### 	Define near-aircraft cloud box
xlon = wrf.getvar(nc2, 'XLONG', timeidx=time_index)
box = np.where(np.logical_and(xlon >=-29.5, xlon<=-26.5))

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

qcloud1 = qcloud1 * float(1e3)
qcloud1.name = 'Cloud LWC, g/kg'

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

qcloud2 = qcloud2 * float(1e3)
qcloud2.name = 'Cloud LWC, g/kg'

###################################
## MAP
###################################

data1 = wrf.to_np(qcloud1[16,:,:])
data2 = wrf.to_np(qcloud2[16,:,:])

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
cbar.set_label(qcloud1.name[-5:])

# Set the map limits.  Not really necessary, but used for demonstration.
# ax.set_xlim(wrf.cartopy_xlim(qcloud1))
# ax.set_ylim(wrf.cartopy_ylim(qcloud1))

# Add the gridlines
ax.gridlines(color="black", linestyle="dotted")

plt.title(qcloud1.name+'\n'+str(qcloud1.Time.values))


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
cbar.set_label(qcloud2.name[-5:])

# Set the map limits.  Not really necessary, but used for demonstration.
# ax.set_xlim(wrf.cartopy_xlim(qcloud2))
# ax.set_ylim(wrf.cartopy_ylim(qcloud2))

# Add the gridlines
ax.gridlines(color="black", linestyle="dotted")

plt.title(qcloud2.name+'\n'+str(qcloud2.Time.values))

plt.show()


###################################
## PROFILE OVER NEST SUBSET
###################################

# Extract the model height
z1 = wrf.getvar(nc1, "z")
z2 = wrf.getvar(nc2, "z")

datax1 = np.nanmean(np.nanmean(qcloud1[:,190:340,np.unique(box[1])],1),1)
datax2 = np.nanmean(np.nanmean(qcloud2[:,190:340,np.unique(box[1])],1),1)
datay1 = np.nanmean(np.nanmean(z1[:,190:340,np.unique(box[1])],1),1)
datay2 = np.nanmean(np.nanmean(z1[:,190:340,np.unique(box[1])],1),1)

plt.plot(datax1,datay1)
plt.plot(datax2,datay2)
plt.ylim([0,2000])
plt.xlabel(qcloud1.name)
plt.ylabel(z1.description)
plt.show()

###################################
## VERTICAL PROFILE @ HALLEY
###################################

##### HALLEY POSITION IN MODEL - NEAREST GRID POINT (LAT/LON)
### D01 = 118,  71 -> Z1[:,71,118]
### D02 = 183, 137 -> Z2[:,137,183]

plt.plot(np.squeeze(qcloud1[:,137,183]),z1[:,137,183],label = file_dir1[0:2])
plt.plot(np.squeeze(qcloud2[:,137,183]),z2[:,137,183],label = file_dir1[0:2])
plt.ylim([0,2000])
plt.title(qcloud1.name+'\n'+str(qcloud1.Time.values))
plt.ylabel(z1.description)
plt.show()
