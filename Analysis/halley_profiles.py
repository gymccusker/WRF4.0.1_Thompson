##-----------------------------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------------------------
# Given a tslist file, for each location inside a model domain (either coarse
# or nested), below files are created:

# 1. pfx*.dNN.TS containing the regular time series output of surface variables.
# 2. pfx.dNN.UU containing a vertical profile of u wind component for each time step
# 3. pfx.dNN.VV containing a vertical profile of v wind component for each time step
# 4. pfx.dNN.TH containing a vertical profile of potential temperature for time step
# 5. pfx.dNN.PH containing a vertical profile of geopotential height in time step
# 6. pfx.dNN.QV containing a vertical profile of water vapor mixing ratio time step
##-----------------------------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------------------------

import numpy as np

import matplotlib
import matplotlib.cm as mpl_cm
# from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import cartopy.crs as crs
import cartopy.feature as cfe

import pandas as pd

from netCDF4 import Dataset
from wrf import to_np, getvar, CoordPair, vertcross

###################################
## Load in data
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

file_dir1 = '5_Archer_Default_AeroClim/'
file_dir2 = '11_Archer_DRIVER_NWFA1D_100e3/'

param = 'TH'

# root_dir = '/data/mac/giyoung/MAC_WRFThompson/' # BAS SCIHUB
root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MAC/WRF_V4.0.1/RUNS/'
obs_dir = '/gws/nopw/j04/ncas_weather/gyoung/MAC/FlightData/Halley/'

###################################
## d01
###################################
filename1_1 = "".join(root_dir+file_dir1+'hal.d01.'+param)
file1_1 = np.loadtxt(filename1_1,skiprows=1)
df1_1 = pd.DataFrame(file1_1, dtype='float')

###################################
## d02
###################################
filename2_1 = "".join(root_dir+file_dir1+'hal.d02.'+param)
file2_1 = np.loadtxt(filename2_1,skiprows=1)
df2_1 = pd.DataFrame(file2_1, dtype='float')

###################################
## d01
###################################
filename1_2 = "".join(root_dir+file_dir2+'hal.d01.'+param)
file1_2 = np.loadtxt(filename1_2,skiprows=1)
df1_2 = pd.DataFrame(file1_2, dtype='float')

###################################
## d02
###################################
filename2_2 = "".join(root_dir+file_dir2+'hal.d02.'+param)
file2_2 = np.loadtxt(filename2_2,skiprows=1)
df2_2 = pd.DataFrame(file2_2, dtype='float')

###################################
## Obs
###################################
filenameObs = "".join(obs_dir+'Halley_Sonde_Data_27-Nov-2015.txt')
dfObs = pd.read_table(filenameObs)

###################################
## WRF (vertical levels)
###################################
# Open the NetCDF file
nc1_1 = Dataset(root_dir+file_dir+'wrfout_d01_2015-11-27_00:00:00')

# Extract the model height and wind speed
Z1_1 = getvar(nc1_1, "z")
nc1_1.close()

# Open the NetCDF file
nc2_1 = Dataset(root_dir+file_dir+'wrfout_d02_2015-11-27_00:00:00')

# Extract the model height and wind speed
Z2_1 = getvar(nc2_1, "z")
nc2_1.close()


###################################
## WRF (vertical levels)
###################################
# Open the NetCDF file
nc1_2 = Dataset(root_dir+file_dir2+'wrfout_d01_2015-11-27_00:00:00')

# Extract the model height and wind speed
Z1_2 = getvar(nc1_2, "z")
nc1_2.close()

# Open the NetCDF file
nc2_2 = Dataset(root_dir+file_dir2+'wrfout_d02_2015-11-27_00:00:00')

# Extract the model height and wind speed
Z2_2 = getvar(nc2_2, "z")
nc2_2.close()

###################################
## Quick check
###################################

# df1_1.iloc[0,:] # prints first row of data
# df1_1.head()

###################################
## Set column names
###################################

df1_1.columns = ['ts_hour','z1','z2','z3','z4','z5','z6','z7','z8','z9','z10','z11','z12','z13','z14','z15']
df2_1.columns = ['ts_hour','z1','z2','z3','z4','z5','z6','z7','z8','z9','z10','z11','z12','z13','z14','z15']
dfObs.columns = ['Obtime','pres','z','temp','dewpt','U','V','RH']

###################################
## Obs quality control
###################################

dfObs.loc[:,'temp'][dfObs.loc[:,'temp'] == 'null'] = np.nan
dfObs.loc[:,'dewpt'][dfObs.loc[:,'dewpt'] == 'null'] = np.nan
dfObs.loc[:,'U'][dfObs.loc[:,'U'] == 'null'] = np.nan
dfObs.loc[:,'V'][dfObs.loc[:,'V'] == 'null'] = np.nan
dfObs.loc[:,'RH'][dfObs.loc[:,'RH'] == 'null'] = np.nan

###################################
## Convert temp to theta
###################################

kap = float(287.05)/float(1005)
tempvar = (pd.to_numeric(dfObs.loc[:,'pres'])/1000)**kap
theta = (pd.to_numeric(dfObs.loc[:,'temp'])+273.16)/tempvar

###################################
## Ignore 1st 24 hours (spin up)
###################################

time = np.where(df1_1.loc[:,'ts_hour']==(24.0+12.0))
ztemp = np.arange(0,15)

##### HALLEY POSITION IN MODEL - NEAREST GRID POINT (LAT/LON)
### D01 = 118,  71 -> Z1[:,71,118]
### D02 = 183, 137 -> Z2[:,137,183]

plt.plot(np.squeeze(df1_1.values[time,1:]),Z1_1[0:15,71,118],label = 'd01')
plt.plot(np.squeeze(df2_1.values[time,1:]),Z2_1[0:15,137,183],label = 'd02')
plt.plot(np.squeeze(df1_2.values[time,1:]),Z1_2[0:15,71,118],'--',label = 'd01')
plt.plot(np.squeeze(df2_2.values[time,1:]),Z2_2[0:15,137,183],'--',label = 'd02')
plt.plot(theta,dfObs.loc[:,'z'],'k',label = 'Obs')
plt.xlabel(param)
plt.ylabel('Z [m]')
plt.ylim([0,1200])
plt.xlim([270,290])
plt.legend()
plt.show()
