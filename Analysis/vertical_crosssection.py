import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature

from netCDF4 import Dataset
from wrf import to_np, getvar, CoordPair, vertcross

# Open the NetCDF file
root_dir = '/data/mac/giyoung/MAC_WRFThompson/'
nc = Dataset(root_dir+'2_Nisg80_ThompsonDefault/wrfout_d02_2015-11-27_00:00:00')

# Extract the model height and wind speed
z = getvar(nc, "z")
wspd =  getvar(nc, "uvmet_wspd_wdir", units="kt")[0,:]

# Create the start point and end point for the cross section
start_point = CoordPair(lat=-74.0, lon=-27.0)
end_point = CoordPair(lat=-75.0, lon=-27.0)

# Compute the vertical cross-section interpolation.  Also, include the
# lat/lon points along the cross-section.
wspd_cross = vertcross(wspd, z, wrfin=nc, start_point=start_point,
                       end_point=end_point, latlon=True, meta=True)

# Create the figure
fig = plt.figure(figsize=(7,6.5))
ax = plt.axes()

# Make the contour plot
wspd_contours = ax.contourf(to_np(wspd_cross), cmap=get_cmap("viridis"))

# Add the color bar
plt.colorbar(wspd_contours, ax=ax)

# Set the x-ticks to use latitude and longitude labels.
coord_pairs = to_np(wspd_cross.coords["xy_loc"])
x_ticks = np.arange(coord_pairs.shape[0])
x_labels = [pair.latlon_str(fmt="{:.2f}, {:.2f}")
            for pair in to_np(coord_pairs)]
ax.set_xticks(x_ticks[::20])
ax.set_xticklabels(x_labels[::20], rotation=45, fontsize=8)

# Set the y-ticks to be height.
# ax.set_ylim([0,3000])
vert_vals = to_np(wspd_cross.coords["vertical"])
v_ticks = np.arange(vert_vals.shape[0])
ax.set_yticks(v_ticks[::20])
ax.set_yticklabels(vert_vals[::20], fontsize=8)

# Set the x-axis and  y-axis labels
ax.set_xlabel("Latitude, Longitude", fontsize=12)
ax.set_ylabel("Height (m)", fontsize=12)

plt.title("Vertical Cross Section of Wind Speed (kt)")

plt.show()