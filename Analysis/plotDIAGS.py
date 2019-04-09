import time
import datetime
import numpy as np

###################################
## Load in WRF data
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

###################################
###################################
## Script to plot out WRF diagnostics using T-E uphys scheme
###################################
###################################

def TDtrans(data, nc, var, time_index):

    ###################################
    ###################################
    ## Transform units from kg/kg to g/m3
    ###################################
    ###################################

    import wrf

    theta = wrf.getvar(nc, 'T', timeidx=time_index) + 300 # potential temperature in K
    theta.name = 'Potential temperature, K'

    pressure = wrf.getvar(nc, 'P', timeidx=time_index) + wrf.getvar(nc, 'PB', timeidx=time_index)   # pressure in Pa
    pressure.name = 'Air pressure, Pa'

    tempvar = float(287.05)/float(1005)
    tempvar0 = (pressure/100000)**tempvar
    temperature = tempvar0 * theta
    temperature.name = 'Air Temperature, K'

    rho = pressure/(float(287.05) * temperature)
    rho.name = 'Air density, kg m-3'

    if var == 'QNCLOUD': data = (data * rho) / float(1e6)
    if var == 'QNWFA': data = (data * rho) / float(1e6)

    return data

def defName(data, var):

    if var == 'QNCLOUD': data.name = 'Cloud droplet number conc, cm-3'
    if var == 'QNWFA': data.name = 'water-friendly aerosol number con, cm-3'

def chooseData(nc1, nc2, var, time_index):

    ###################################
    ###################################
    ## Read in data from netCDF using WRF-Python
    ###################################
    ###################################

    import wrf

    data1 = wrf.getvar(nc1, var, timeidx=time_index)
    data2 = wrf.getvar(nc2, var, timeidx=time_index)
    # data3 = wrf.getvar(nc3, var, timeidx=time_index)

    data1 = TDtrans(data1, nc1, var, time_index)
    data2 = TDtrans(data2, nc2, var, time_index)
    # data3 = TDtrans(data3, nc3, var, time_index)

    return data1, data2

def plotmap(nc1, nc2, time_index, z_index):

    import matplotlib
    import matplotlib.cm as mpl_cm
    import matplotlib.pyplot as plt
    import cartopy.crs as crs
    import cartopy.feature as cfe
    import wrf

    ###################################
    ###################################
    ## Plot map (with cartopy)
    ###################################
    ###################################

    ## Load in chosen data
    data1, data2 = chooseData(nc1, nc2, var, time_index)

    data1 = defName(data1, var)
    data2 = defName(data2, var)

    ## Get the latitude and longitude points
    lats, lons = wrf.latlon_coords(data1)

    ## Get the cartopy mapping object
    cart_proj = wrf.get_cartopy(data1)

    data1 = wrf.to_np(data1[z_index,:,:])
    data2 = wrf.to_np(data2[z_index,:,:])
    # data3 = wrf.to_np(data3[z_index,:,:])

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
    cbar.set_label(data1.name[-5:])

    # Set the map limits.  Not really necessary, but used for demonstration.
    # ax.set_xlim(wrf.cartopy_xlim(qnwfa1))
    # ax.set_ylim(wrf.cartopy_ylim(qnwfa1))

    # Add the gridlines
    ax.gridlines(color="black", linestyle="dotted")

    plt.title(data1.name+'\n'+str(data1.Time.values))

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
    cbar.set_label(data2.name[-5:])

    # Set the map limits.  Not really necessary, but used for demonstration.
    # ax.set_xlim(wrf.cartopy_xlim(qnwfa2))
    # ax.set_ylim(wrf.cartopy_ylim(qnwfa2))

    # Add the gridlines
    ax.gridlines(color="black", linestyle="dotted")

    plt.title(data2.name+'\n'+str(data2.Time.values))

    plt.show()

def plotProfile(nc1, nc2, var, time_index, file_dir1, file_dir2):

    ###################################
    ###################################
    ## Plot profile of diags at Halley
    ###################################
    ###################################

    import wrf
    import matplotlib
    import matplotlib.cm as mpl_cm
    import matplotlib.pyplot as plt

    ## Load in chosen data
    data1, data2 = chooseData(nc1, nc2, var, time_index)

    data1 = defName(data1, var)
    data2 = defName(data2, var)

    # Extract the model height
    z1 = wrf.getvar(nc1, "z")
    z2 = wrf.getvar(nc2, "z")

    ##### HALLEY POSITION IN MODEL - NEAREST GRID POINT (LAT/LON)
    ### D01 = 118,  71 -> Z1[:,71,118]
    ### D02 = 183, 137 -> Z2[:,137,183]

    plt.plot(np.squeeze(data1[:,137,183]),z1[:,137,183],label = file_dir1[0:2])
    plt.plot(np.squeeze(data2[:,137,183]),z2[:,137,183],label = file_dir2[0:2])
    plt.ylim([0,2000])
    plt.title(data1.name+'\n'+str(data1.Time.values))
    plt.legend()
    plt.ylabel(z1.description)
    plt.show()

def plotSubset(nc1, nc2, var, time_index, file_dir1, file_dir2):

    ###################################
    ## PROFILE OVER NEST SUBSET
    ###################################

    import wrf
    import matplotlib
    import matplotlib.cm as mpl_cm
    import matplotlib.pyplot as plt

    ## Load in chosen data
    data1, data2 = chooseData(nc1, nc2, var, time_index)

    data1 = defName(data1, var)
    data2 = defName(data2, var)

    #### 	Define near-aircraft cloud box
    xlon = wrf.getvar(nc2, 'XLONG', timeidx=time_index)
    box = np.where(np.logical_and(xlon >=-29.5, xlon<=-26.5))

    # Extract the model height
    z1 = wrf.getvar(nc1, "z")
    z2 = wrf.getvar(nc2, "z")

    datax1 = np.nanmean(np.nanmean(data1[:,190:340,np.unique(box[1])],1),1)
    datax2 = np.nanmean(np.nanmean(data2[:,190:340,np.unique(box[1])],1),1)
    datay1 = np.nanmean(np.nanmean(z1[:,190:340,np.unique(box[1])],1),1)
    datay2 = np.nanmean(np.nanmean(z1[:,190:340,np.unique(box[1])],1),1)

    plt.plot(datax1,datay1)
    plt.plot(datax2,datay2)
    plt.ylim([0,2000])
    plt.xlabel(data1.name)
    plt.ylabel(z1.description)
    plt.show()

def main():

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
    ## 18_Archer_initialise_real_qnwfanow_K1_100e6/             ## GIVES THE SAME AS #5...?
    ## 19_Archer_initialise_real_qnwfanow_K1_x2/
    ## 20_Archer_initialise_real_qnwfanow_K1_x10/
    ## 21_Archer_initialise_real_qnwfanow_i1j_x10/
    ## 22_Archer_initialise_real_qnwfanow_x2_17redo/            ## GIVES THE SAME AS #5...?

    from netCDF4 import Dataset
    import wrf
    import xarray as xr
    import numpy as np

    START_TIME = time.time()
    print ''
    print 'Start: ' + time.strftime("%c")
    print ''

    file_dir1 = '2_Nisg80_ThompsonDefault/'
    file_dir2 = '3_Nisg80_ThompsonAeroClim/'
    # file_dir3 = '4_Nisg80_Thompson_naCCN0408_naCCN1100/'

    # root_dir = '/gws/nopw/j04/ncas_weather/gyoung/MAC/WRF_V4.0.1/RUNS/'
    root_dir = '/data/mac/giyoung/MAC_WRFThompson/'

    ### Read in netCDF files
    nc1 = Dataset(root_dir+file_dir1+'wrfout_d02_2015-11-27_00:00:00')
    nc2 = Dataset(root_dir+file_dir2+'wrfout_d02_2015-11-27_00:00:00')
    # nc3 = Dataset(root_dir+file_dir3+'wrfout_d02_2015-11-27_00:00:00')

    ## Choose time_index to plot (32 = 16h)
    time_index = 32

    ## Choose model level index (for plotting map)
    z_index = 16

    ## Choose variable to plot
    var = 'QNCLOUD'

    ## Plot map (cartopy)
    # map = plotmap(nc1, nc2, nc3, time_index, z_index)

    ## Plot vertical profile at Halley
    profile = plotProfile(nc1, nc2, var, time_index, file_dir1, file_dir2)

    ## Plot average diagnostics over nest subset
    subset = plotSubset(nc1, nc2, var, time_index, file_dir1, file_dir2)

    END_TIME = time.time()
    print ''
    print 'End: ' + time.strftime("%c")
    print ''

if __name__ == '__main__':

    main()
