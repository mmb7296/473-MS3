import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import cartopy
import cartopy.crs as ccrs
import netCDF4 as nc
from datetime import datetime, timedelta

time = datetime(2018,3,1,0,0)
endmo = datetime(2018,3,31,23)
readpath = './Hourly_M_2018_pcp_sum.nc'

f = nc.Dataset(readpath)

#plot lat, lon, precip
lat = f.variables['latitude']
lon = f.variables['longitude']
DCPCP = f.variables['pcp_total_rate']

#set up axes
fig, ax = plt.subplots(nrows=1, ncols=1, dpi=300, figsize=(20, 10), subplot_kw={'projection': ccrs.PlateCarree()})
hour = 12
valid = time + timedelta(hours = hour)
DCPrecip = DCPCP[hour, :, :]
ax.set_global()
ax.coastlines()
ax.gridlines(color='grey', linestyle='--')
ax.set_extent((90, 120, -10, 10), crs=ccrs.PlateCarree())
ax.set_xticks([90, 95, 100, 105, 110, 115, 120])
ax.set_yticks([-10, -5, 0, 5, 10]) 
ax.set_xlabel('Longitude (°E)')
ax.set_ylabel('Latitude (°N')
cs = ax.pcolormesh(lon, lat, DCPrecip, transform=ccrs.PlateCarree(), vmax=.75, vmin=0, cmap='Blues')

cb = plt.colorbar(cs, orientation='vertical', ax=ax)
cb.set_label('Total Precipitation (mm/h)', fontsize=15)
plt.title('Total Precipitation: {0}'.format(valid.strftime('%HZ %B %d %Y')), fontsize=25)
plt.show()
