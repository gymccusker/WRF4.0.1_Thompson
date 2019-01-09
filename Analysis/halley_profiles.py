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
from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import cartopy.crs as crs
import cartopy.feature as cfe

import pandas as pd

from netCDF4 import Dataset
from wrf import to_np, getvar, CoordPair, vertcross

###################################
## Load in data
###################################

file_dir = '2_Nisg80_ThompsonDefault/'
param = 'TH'

root_dir = '/data/mac/giyoung/MAC_WRFThompson/'
obs_dir = '/data/scihub-users/giyoung/MAC/'

###################################
## d01
###################################
filename1 = "".join(root_dir+file_dir+'hal.d01.'+param)
file1 = np.loadtxt(filename1,skiprows=1)
df1 = pd.DataFrame(file1, dtype='float')

###################################
## d02
###################################
filename2 = "".join(root_dir+file_dir+'hal.d02.'+param)
file2 = np.loadtxt(filename2,skiprows=1)
df2 = pd.DataFrame(file2, dtype='float')

###################################
## Obs
###################################
filenameObs = "".join(obs_dir+'Halley_Sonde_Data_27-Nov-2015.txt')
dfObs = pd.read_table(filenameObs)

###################################
## WRF (vertical levels)
###################################
# Open the NetCDF file
nc1 = Dataset(root_dir+file_dir+'wrfout_d01_2015-11-27_00:00:00')

# Extract the model height and wind speed
Z1 = getvar(nc1, "z")
nc1.close()

# Open the NetCDF file
nc2 = Dataset(root_dir+file_dir+'wrfout_d02_2015-11-27_00:00:00')

# Extract the model height and wind speed
Z2 = getvar(nc2, "z")
nc2.close()

###################################
## Quick check
###################################
 
# df1.iloc[0,:] # prints first row of data
# df1.head()

###################################
## Set column names
###################################

df1.columns = ['ts_hour','z1','z2','z3','z4','z5','z6','z7','z8','z9','z10','z11','z12','z13','z14','z15']
df2.columns = ['ts_hour','z1','z2','z3','z4','z5','z6','z7','z8','z9','z10','z11','z12','z13','z14','z15']
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

# data1['theta'] = nc1.variables['T'][time_sci]+300 # potential temperature in K
# data1['p'] = (nc1.variables['P'][time_sci[0]]+nc1.variables['PB'][time_sci[0]])   # pressure in Pa
# tempvar = constants.R/float(1005)
# tempvar0 = (data1['p']/100000)**tempvar       
# data1['Tk'] = tempvar0*data1['theta']
# data1['rho'] = data1['p']/(constants.R*data1['Tk'])

kap = float(287.05)/float(1005)
tempvar = (pd.to_numeric(dfObs.loc[:,'pres'])/1000)**kap
theta = (pd.to_numeric(dfObs.loc[:,'temp'])+273.16)/tempvar

###################################
## Ignore 1st 24 hours (spin up)
###################################

time = np.where(df1.loc[:,'ts_hour']==(24.0+12.0))
ztemp = np.arange(0,15)

##### HALLEY POSITION IN MODEL - NEAREST GRID POINT (LAT/LON)
### D01 = 118,  71 -> Z1[:,71,118]
### D02 = 183, 137 -> Z2[:,137,183]

plt.plot(np.squeeze(df1.values[time,1:]),Z1[0:15,71,118],label = 'd01')
plt.plot(np.squeeze(df2.values[time,1:]),Z2[0:15,137,183],label = 'd02')
plt.plot(theta,dfObs.loc[:,'z'],'k',label = 'Obs')
plt.xlabel(param)
plt.ylabel('Z [m]')
plt.ylim([0,2000])
plt.xlim([270,290])
plt.legend()
plt.show()