import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import cartopy
import cartopy.crs as ccrs
import netCDF4 as nc
from datetime import datetime, timedelta

time = datetime(2018,3,1,0,0)
endmo = datetime(2018,3,31,23)
ufile = nc.Dataset('./473/e5.oper.an.sfc.128_165_10u.ll025sc.{0}_{1}.nc'.format(time.strftime('%Y%m%d%H'), endmo.strftime('%Y%m%d%H')))
vfile = nc.Dataset('./473/e5.oper.an.sfc.128_166_10v.ll025sc.{0}_{1}.nc'.format(time.strftime('%Y%m%d%H'), endmo.strftime('%Y%m%d%H')))
lat = ufile.variables['latitude']
lon = ufile.variables['longitude']
wind_u = ufile.variables['VAR_10U']
wind_v = vfile.variables['VAR_10V']

plt.figure(figsize=(10, 10))

ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent([90, 120, -10, 10], crs=ccrs.PlateCarree())
ax.coastlines()
ax.set_xlabel('Longitude (°E)')
ax.set_ylabel('Latitude (°N')
ax.set_xticks(range(90, 125, 5), crs=ccrs.PlateCarree())
ax.set_yticks(range(-10, 15, 5), crs=ccrs.PlateCarree())
hour = 12
valid = time + timedelta(hours = hour)
Q = plt.quiver(lon[::6],lat[::6],wind_u[hour,::6,::6],wind_v[hour,::6,::6], scale_units='xy', scale=3, width=0.0015)
ax.set_title('10m Winds: {0}'.format(valid.strftime('%HZ %B %d %Y')), fontsize=16)

ax.quiverkey(Q, X=0.90, Y=1.015, U=2, label='Quiver key, length = 2')

plt.show()
