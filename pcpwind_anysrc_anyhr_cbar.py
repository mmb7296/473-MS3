import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import cartopy
import cartopy.crs as ccrs
import netCDF4 as nc
from datetime import datetime, timedelta

time = datetime(2018,3,1,0,0)
endmo = datetime(2018,3,31,23)
molist = ['J','F','M','A','Ma','June','July','Aug','S','O','N','D']
monthstr = molist[time.month-1]
readpath = './Hourly_{0}_{1}_pcp_sum.nc'.format(monthstr, time.year) 

#the precipitation datasets are named in a way that doesn't allow use of datetime for the month names.
#I've included a list of abbreviations, then left it so the prospective user could choose one when picking
#the time.

pcpsources = ['Total', 'MCS', 'Deep Convection', 'Congestus']
pcpsource = pcpsources[0] #Can change to retrieve different precipitation source


f = nc.Dataset(readpath)

lat = f.variables['latitude']
lon = f.variables['longitude']

if pcpsource == 'Total':
    pcp = f.variables['pcp_total_rate']
elif pcpsource == 'MCS':
    pcp = f.variables['pcp_mcs_rate']
elif pcpsource == 'Deep Convection':
    pcp = f.variables['pcp_deepconvec_rate']
else:
    pcp = f.variables['pcp_congestus_rate']
    
ufile = nc.Dataset('./473/e5.oper.an.sfc.128_165_10u.ll025sc.{0}_{1}.nc'.format(time.strftime('%Y%m%d%H'), endmo.strftime('%Y%m%d%H')))
vfile = nc.Dataset('./473/e5.oper.an.sfc.128_166_10v.ll025sc.{0}_{1}.nc'.format(time.strftime('%Y%m%d%H'), endmo.strftime('%Y%m%d%H')))
wind_u = ufile.variables['VAR_10U']
wind_v = vfile.variables['VAR_10V']


#set up axes
fig, ax = plt.subplots(nrows=1, ncols=1, dpi=300, figsize=(20, 10), subplot_kw={'projection': ccrs.PlateCarree()})
hour = 12
valid = time + timedelta(hours = hour)
precip = pcp[hour, :, :]
ax.set_global()
ax.coastlines()
ax.gridlines(color='grey', linestyle='--')
ax.set_extent((90, 120, -10, 10), crs=ccrs.PlateCarree())
ax.set_xticks([90, 95, 100, 105, 110, 115, 120])
ax.set_yticks([-10, -5, 0, 5, 10]) 
#set up unique colorbar
cbarcolors = ['white', 'azure', 'paleturquoise', 'aquamarine', 'mediumaquamarine', 'darkturquoise', 'teal', 'darkslategrey']
custommap = mpl.colors.ListedColormap(cbarcolors)
#plots precipitation and wind together
cs = ax.pcolormesh(lon, lat, precip, transform=ccrs.PlateCarree(), vmin=0, vmax=0.8, cmap=custommap)
lat = ufile.variables['latitude']
lon = ufile.variables['longitude']
Q = plt.quiver(lon[::6],lat[::6],wind_u[hour,::6,::6],wind_v[hour,::6,::6], scale_units='xy', scale=3, width=0.0015)
ax.set_title('10m Winds: {0}'.format(valid.strftime('%HZ %B %d %Y')), fontsize=16)
ax.quiverkey(Q, X=1.075, Y=1.015, U=2, label='Quiver key, length = 2')
cb = plt.colorbar(cs, orientation='vertical', ax=ax)
cb.set_label('{0} Precipitation (mm/h)'.format(pcpsource), fontsize=15)
plt.title('{0} Precipitation and 10m Winds: {1}'.format(pcpsource, valid.strftime('%HZ %B %d %Y')), fontsize=25)
plt.show()
