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
file_dir = '3_Nisg80_ThompsonAeroClim/'

root_dir = '/data/mac/giyoung/MAC_WRFThompson/'
nc = Dataset(root_dir+file_dir+'wrfout_d01_2015-11-27_00:00:00')
qnwfa = wrf.getvar(nc, 'QNWFA', timeidx=32)

## Quick Plot to check all is well
# qnwfa.plot()

## Get the latitude and longitude points
lats, lons = wrf.latlon_coords(qnwfa)

## Get the cartopy mapping object
cart_proj = wrf.get_cartopy(qnwfa)


###################################
## MAP
###################################

data = wrf.to_np(qnwfa[16,:,:])

# Create a figure
fig = plt.figure(figsize=(6,5))
# Set the GeoAxes to the projection used by WRF
ax = plt.axes(projection=cart_proj)

# Add coastlines
ax.coastlines('50m', linewidth=0.8)
ax.add_feature(cfe.NaturalEarthFeature('physical', 'antarctic_ice_shelves_lines', 
                                       '50m', linewidth=1.0, edgecolor='k', facecolor='none') )

# Plot contours
plt.contourf(wrf.to_np(lons), wrf.to_np(lats), data, 10, 
                transform=crs.PlateCarree(), cmap = mpl_cm.Reds)

# Add a color bar
cbar = plt.colorbar(ax=ax, shrink=.62)
cbar.set_label(qnwfa.units)

# Set the map limits.  Not really necessary, but used for demonstration.
ax.set_xlim(wrf.cartopy_xlim(qnwfa))
ax.set_ylim(wrf.cartopy_ylim(qnwfa))

# Add the gridlines
ax.gridlines(color="black", linestyle="dotted")

plt.title(qnwfa.description+'\n'+str(qnwfa.Time.values))

plt.show()

###################################
## VERTICAL CROSS-SECTION
###################################

# Extract the model height and wind speed
z = wrf.getvar(nc, "z")

# Create the start point and end point for the cross section
start_point = wrf.CoordPair(lat=-74.0, lon=-27.0)
end_point = wrf.CoordPair(lat=-75.0, lon=-27.0)

# Compute the vertical cross-section interpolation.  Also, include the
# lat/lon points along the cross-section.
height_range = np.arange(0,26)

qnwfa_cross = wrf.vertcross(qnwfa, z, wrfin=nc, start_point=start_point,
                       end_point=end_point, latlon=True, meta=True)

# Create the figure
fig = plt.figure(figsize=(7,6.5))
ax = plt.axes()

# Make the contour plot
qnwfa_contours = ax.contourf(wrf.to_np(qnwfa_cross), cmap=mpl_cm.viridis)
# ax.set_ylim([0,25])

# Add the color bar
plt.colorbar(qnwfa_contours, ax=ax)

# Set the x-ticks to use latitude and longitude labels.
coord_pairs = wrf.to_np(qnwfa_cross.coords["xy_loc"])
x_ticks = np.arange(coord_pairs.shape[0])
x_labels = [pair.latlon_str(fmt="{:.2f}, {:.2f}")
            for pair in wrf.to_np(coord_pairs)]
ax.set_xticks(x_ticks[::20])
ax.set_xticklabels(x_labels[::20], rotation=45, fontsize=8)

# Set the y-ticks to be height.
vert_vals = wrf.to_np(qnwfa_cross.coords["vertical"])
v_ticks = np.arange(vert_vals.shape[0])
ax.set_yticks(v_ticks[::20])
ax.set_yticklabels(vert_vals[::20], fontsize=8)

# Set the x-axis and  y-axis labels
ax.set_xlabel("Latitude, Longitude", fontsize=12)
ax.set_ylabel("Height (m)", fontsize=12)

plt.title(qnwfa.description+'\n'+str(qnwfa.Time.values))

plt.show()