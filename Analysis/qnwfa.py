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
## 11_Archer_DRIVER_NWFA1D_100e3/       ## WRONG NAMELIST
## 12_Archer_DRIVER_NWFA1D_x2/       ## WRONG NAMELIST -- failed job
## 13_Archer_DRIVER_NWFA1D_x05/       ## WRONG NAMELIST
## 14_Archer_DRIVER_NWFA1D_150e3/       ## WRONG NAMELIST
## 15_Archer_DRIVER_NWFA1D_150e3_K1/       ## WRONG NAMELIST
## 16_Archer_DRIVER_NWFA1D_100e3_K1/       ## WRONG NAMELIST
## 17_Archer_initialise_real_qnwfanow_x2/       ## WRONG NAMELIST
## 18_Archer_initialise_real_qnwfanow_K1_100e6/
## 19_Archer_initialise_real_qnwfanow_K1_x2/
## 20_Archer_initialise_real_qnwfanow_K1_x10/
## 21_Archer_initialise_real_qnwfanow_i1j_x10/

file_dir1 = '5_Archer_Default_AeroClim/'
file_dir2 = '5_Archer_Default_AeroClim/'

root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MAC/WRF_V4.0.1/RUNS/'

time_index = 32

nc1 = Dataset(root_dir+file_dir1+'wrfout_d02_2015-11-27_00:00:00')
qnwfa1 = wrf.getvar(nc1, 'QNWFA', timeidx=time_index)

nc2 = Dataset(root_dir+file_dir2+'wrfout_d02_2015-11-27_00:00:00')
qnwfa2 = wrf.getvar(nc2, 'QNWFA', timeidx=time_index)

## Quick Plot to check all is well
# qnwfa.plot()

## Get the latitude and longitude points
lats, lons = wrf.latlon_coords(qnwfa1)

## Get the cartopy mapping object
cart_proj = wrf.get_cartopy(qnwfa1)

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

qnwfa1 = (qnwfa1 * rho1) / float(1e6)
qnwfa1.name = 'water-friendly aerosol number con, cm-3'

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

qnwfa2 = (qnwfa2 * rho2) / float(1e6)
qnwfa2.name = 'water-friendly aerosol number con, cm-3'

###################################
## MAP
###################################

data1 = wrf.to_np(qnwfa1[0,:,:])
data2 = wrf.to_np(qnwfa2[0,:,:])

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
cbar.set_label(qnwfa1.name[-5:])

# Set the map limits.  Not really necessary, but used for demonstration.
# ax.set_xlim(wrf.cartopy_xlim(qnwfa1))
# ax.set_ylim(wrf.cartopy_ylim(qnwfa1))

# Add the gridlines
ax.gridlines(color="black", linestyle="dotted")

plt.title(qnwfa1.name+'\n'+str(qnwfa1.Time.values))


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
cbar.set_label(qnwfa2.name[-5:])

# Set the map limits.  Not really necessary, but used for demonstration.
# ax.set_xlim(wrf.cartopy_xlim(qnwfa2))
# ax.set_ylim(wrf.cartopy_ylim(qnwfa2))

# Add the gridlines
ax.gridlines(color="black", linestyle="dotted")

plt.title(qnwfa2.name+'\n'+str(qnwfa2.Time.values))

plt.show()

###################################
## VERTICAL CROSS-SECTION
###################################

# Extract the model height
z1 = wrf.getvar(nc1, "z")
z2 = wrf.getvar(nc2, "z")

##### HALLEY POSITION IN MODEL - NEAREST GRID POINT (LAT/LON)
### D01 = 118,  71 -> Z1[:,71,118]
### D02 = 183, 137 -> Z2[:,137,183]

plt.plot(np.squeeze(qnwfa1[:,137,183]),z1[:,137,183],label = file_dir1[0:2])
plt.plot(np.squeeze(qnwfa2[:,137,183]),z2[:,137,183],label = file_dir2[0:2])
plt.ylim([0,2000])
plt.title(qnwfa1.name+'\n'+str(qnwfa1.Time.values))
plt.legend()
plt.ylabel(z1.description)
plt.show()

# # Create the start point and end point for the cross section
# start_point = wrf.CoordPair(lat=-74.0, lon=-27.0)
# end_point = wrf.CoordPair(lat=-75.0, lon=-27.0)

# # Compute the vertical cross-section interpolation.  Also, include the
# # lat/lon points along the cross-section.
# height_range = np.arange(0,26)

# qnwfa_cross1 = wrf.vertcross(qnwfa1, z1, wrfin=nc1, start_point=start_point,
#                        end_point=end_point, latlon=True, meta=True)

# # Create the figure
# fig = plt.figure(figsize=(7,6.5))
# ax = plt.axes()

# # Make the contour plot
# qnwfa_contours1 = ax.contourf(wrf.to_np(qnwfa_cross1), cmap=mpl_cm.viridis)
# # ax.set_ylim([0,25])

# # Add the color bar
# plt.colorbar(qnwfa_contours1, ax=ax)

# # Set the x-ticks to use latitude and longitude labels.
# coord_pairs = wrf.to_np(qnwfa_cross1.coords["xy_loc"])
# x_ticks = np.arange(coord_pairs.shape[0])
# x_labels = [pair.latlon_str(fmt="{:.2f}, {:.2f}")
#             for pair in wrf.to_np(coord_pairs)]
# ax.set_xticks(x_ticks[::20])
# ax.set_xticklabels(x_labels[::20], rotation=45, fontsize=8)

# # Set the y-ticks to be height.
# vert_vals = wrf.to_np(qnwfa_cross1.coords["vertical"])
# v_ticks = np.arange(vert_vals.shape[0])
# ax.set_yticks(v_ticks[::20])
# ax.set_yticklabels(vert_vals[::20], fontsize=8)

# # Set the x-axis and  y-axis labels
# ax.set_xlabel("Latitude, Longitude", fontsize=12)
# ax.set_ylabel("Height (m)", fontsize=12)

# plt.title(qnwfa1.name+'\n'+str(qnwfa1.Time.values))

# plt.show()
