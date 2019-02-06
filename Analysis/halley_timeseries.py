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

## 2_Nisg80_ThompsonDefault/
## 3_Nisg80_ThompsonAeroClim/
## 4_Nisg80_Thompson_naCCN0408_naCCN1100/
## 5_Archer_Default_AeroClim/
## 6_Archer_NWFApl100_AeroClim/

file_dir1 = '5_Archer_Default_AeroClim/'
file_dir2 = '6_Archer_NWFApl100_AeroClim/'

index = 'gsw'

root_dir = '/data/mac/giyoung/MAC_WRFThompson/'

###################################
## d01
###################################
filename1 = "".join(root_dir+file_dir1+'hal.d02.TS')
file1 = np.loadtxt(filename1,skiprows=1)
df1 = pd.DataFrame(file1)

###################################
## d02
###################################
filename2 = "".join(root_dir+file_dir2+'hal.d02.TS')
file2 = np.loadtxt(filename2,skiprows=1)
df2 = pd.DataFrame(file2)

###################################
## Quick check
###################################

df1.iloc[0,:] # prints first row of data
df1.head()

# data1 = df.values
# data1[0,0] ## first value
# df.columns ## headers

###################################
## Set column names
###################################

df1.columns = ['id','ts_hour','id_tsloc','ix','iy','t','q','u','v','psfc','glw','gsw','hfx','lh','tsk','tslb(1)','rainc','rainnc','clw']
df2.columns = ['id','ts_hour','id_tsloc','ix','iy','t','q','u','v','psfc','glw','gsw','hfx','lh','tsk','tslb(1)','rainc','rainnc','clw']

###################################
## Ignore 1st 24 hours (spin up)
###################################

plt.plot(df1.loc[np.size(df1.values[:,0])/float(2)-1:,'ts_hour']-24,df1.loc[np.size(df1.values[:,0])/float(2)-1:,index],label='d01')
plt.plot(df2.loc[np.size(df2.values[:,0])/float(2)-1:,'ts_hour']-24,df2.loc[np.size(df2.values[:,0])/float(2)-1:,index],label='d02')
plt.xlabel('Time, h [27-Nov-2018]')
plt.ylabel(index)
plt.xlim([0,24])
plt.legend()
plt.show()
