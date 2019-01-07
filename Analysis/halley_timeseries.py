##-----------------------------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------------------------
# The first line in a time-series output looks like this:
# 	
# 	NZCM McMurdo               2  7 mcm   (-77.850, 166.710) ( 153, 207) (-77.768, 166.500)   81.8 meters
# 	
# 	Those are name of the station, grid ID, time-series ID, station lat/lon, grid indices (nearest grid point to
# 	the station location), grid lat/lon, elevation.
# 	
# 	The variables from the time series output are:
# 	
# 	id, ts_hour, id_tsloc, ix, iy, t, q, u, v, psfc, glw, gsw, hfx, lh, tsk, tslb(1), rainc, rainnc, clw
# 	
# 	id:             grid ID
# 	ts_hour:        forecast time in hours
# 	id_tsloc:       time series ID
# 	ix,iy:          grid location (nearest grid to the station)
# 	t:              2 m Temperature (K)
# 	q:              2 m vapor mixing ratio (kg/kg)
# 	u:              10 m U wind (earth-relative)
# 	v:              10 m V wind (earth-relative)
# 	psfc:           surface pressure (Pa)
# 	glw:            downward longwave radiation flux at the ground (W/m^2, downward is positive)
# 	gsw:            net shortwave radiation flux at the ground (W/m^2, downward is positive)
# 	hfx:            surface sensible heat flux (W/m^2, upward is positive)
# 	lh:             surface latent heat flux (W/m^2, upward is positive)
# 	tsk:            skin temperature (K)
# 	tslb(1):        top soil layer temperature (K)
# 	rainc:          rainfall from a cumulus scheme (mm)
# 	rainnc:         rainfall from an explicit scheme (mm)
# 	clw:            total column-integrated water vapor and cloud variables
##-----------------------------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------------------------



import numpy as np

import matplotlib
import matplotlib.cm as mpl_cm
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import cartopy.crs as crs
import cartopy.feature as cfe

import pandas as pd

###################################
## Load in data
###################################

root_dir = '/data/mac/giyoung/MAC_WRFThompson/'
filename = "".join(root_dir+'2_Nisg80_ThompsonDefault/hal.d01.PH')
readtab = pd.read_table(filename,',')
df = pd.DataFrame(readtab)