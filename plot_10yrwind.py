import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import cartopy
import cartopy.crs as ccrs
import netCDF4 as nc

year = 2010
p_level = np.array((1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175,
                    200, 225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700,
                    750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000), dtype = np.float64)
level = 16 #250 mb; 21 for 500, 26 for 750, 30 for 850
matplotlib.rcParams.update({'font.size': 18})
upath = "./UMean_{0}-{1}".format(year)
ufile = nc.Dataset(upath)
vpath = "./VMean_{0}-{1}.nc".format(year)
vfile = nc.Dataset(vpath)
lat = np.array(ufile.variables['latitude'])
lon = np.array(ufile.variables['longitude'])

print("Variables imported")

(lon_min, lon_max, lat_min, lat_max) = (89, 121, -11, 11)
mean_u = ufile.variables['u wind'][:,:,:,:]
mean_v = vfile.variables['v wind'][:,:,:,:]
fig, axes = plt.subplots(nrows=3, ncols=4, dpi=300, figsize=(20,11), subplot_kw={'projection': ccrs.PlateCarree()}, sharex=True, sharey=True)

plt.setp(axes[-1, :], xlabel='Longitude (°E)')
plt.setp(axes[:, 0], ylabel='Latitude (°N)')
for hour in range(0, 24, 2):
    wind_u = mean_u[hour,level,:,:]
    wind_v = mean_v[hour,level,:,:]
    print("Variables read")
    if hour in range(0, 17):
        WIB = hour + 7
    else:
        WIB = hour + 7 - 24
    if WIB in range(0, 8):
        ycoord = 2
    elif WIB in range(9, 16):
        ycoord = 0
    else:
        ycoord = 1
    if WIB in [1, 9, 17]:
        xcoord = 0
    elif WIB in [3, 11, 19]:
        xcoord = 1
    elif WIB in [5, 13, 21]:
        xcoord = 2
    else:
        xcoord = 3
    ax = axes[ycoord, xcoord]
    ax.set_global()
    ax.coastlines()
    ax.set_extent([90, 120, -10, 10], crs=ccrs.PlateCarree())
    ax.coastlines()
    ax.set_xticks(range(90, 125, 5), crs=ccrs.PlateCarree())
    ax.set_yticks(range(-10, 15, 5), crs=ccrs.PlateCarree())
    print("Outline of graph created")
    Q = ax.quiver(lon[::5], lat[::5], wind_u[::5, ::5], wind_v[::5, ::5], scale_units='xy', scale=3)
    print("Quiver created")
    ax.set_title('{0}00 WIB'.format(WIB))

    qk = plt.quiverkey(Q, 
                  1, 1.04,                  
                   5,str(5)+' m/s',   
                   labelpos='E',                
                   coordinates='axes',
                   )
plt.suptitle('Mean Wind Patterns ({0} hPa)'.format(int(p_level[level])), y=0.95, fontsize=25)
plt.show()
