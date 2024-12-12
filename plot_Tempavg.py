import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import cartopy
import cartopy.crs as ccrs
import netCDF4 as nc

readpath = './TMean_{0}-{1}.nc'.format(yearlist[0],yearlist[-1])

f = nc.Dataset(readpath)

#plot lat, lon, precip
lat = f.variables['latitude']
lon = f.variables['longitude']
T = f.variables['T']
level = 30 #850 hPa

#set up axes
fig, axes = plt.subplots(nrows=3, ncols=4, dpi=300, figsize=(20, 10), subplot_kw={'projection': ccrs.PlateCarree()})
test1 = plt.setp(axes[-1, :], xlabel='Longitude (°E)')
test2 = plt.setp(axes[:, 0], ylabel='Latitude (°N)')

#"WIB" is converting from UTC to local time (WIB)
for hour in range(0, 24, 2):
    Temp = T[hour, level, :, :]
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
    ax.gridlines(color='grey', linestyle='--')
    ax.set_extent((90, 120, -10, 10), crs=ccrs.PlateCarree())
    ax.set_xticks([90, 95, 100, 105, 110, 115, 120])
    ax.set_yticks([-10, -5, 0, 5, 10]) 
    cs = ax.pcolormesh(lon, lat, Temp, transform=ccrs.PlateCarree(), cmap='jet')
    ax.set_title('{0}00 WIB'.format(WIB))
cb = plt.colorbar(cs, orientation='vertical', shrink=0.75, ax=axes)
cb.set_label('Temperature (K)', fontsize=15)
plt.suptitle('Temperature Daily Average 2010-2019', x=0.45, y=0.95, fontsize=25)
plt.show()
