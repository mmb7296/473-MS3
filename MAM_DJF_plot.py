from datetime import datetime
from datetime import timedelta
import xarray

pcpdict = {}
yearlist = [2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019]
seasonstring = 'MAM'
#seasonstring = 'DJF'

if seasonstring =='MAM':
    for year in yearlist:
        marchpath = 'Hourly_M_{0}_pcp_sum.nc'.format(year)
        aprilpath = 'Hourly_A_{0}_pcp_sum.nc'.format(year)
        maypath = 'Hourly_Ma_{0}_pcp_sum.nc'.format(year)
        
        marchfile = xarray.open_dataset(marchpath)
        aprilfile = xarray.open_dataset(aprilpath)
        mayfile = xarray.open_dataset(maypath)
        
        pcpdict["{0}".format(year)] = xarray.concat([marchfile,aprilfile,mayfile], dim="hour")
else:
    for year in yearlist:
        janpath = 'Hourly_J_{0}_pcp_sum.nc'.format(year)
        febpath = 'Hourly_F_{0}_pcp_sum.nc'.format(year)
        decpath = 'Hourly_D_{0}_pcp_sum.nc'.format(year)
        
        janfile = xarray.open_dataset(janpath)
        febfile = xarray.open_dataset(febpath)
        decfile = xarray.open_dataset(decpath)
        
        pcpdict["{0}".format(year)] = xarray.concat([janfile,febfile,decfile], dim="hour")
merge = xarray.concat([pcpdict["2010"], pcpdict["2011"], pcpdict["2012"], pcpdict["2013"], pcpdict["2014"], pcpdict["2015"], pcpdict["2016"], pcpdict["2017"], pcpdict["2018"], pcpdict["2019"]], dim="hour")
startime = datetime(yearlist[0],1,1,0,0)
endtime = datetime(yearlist[-1],12,31,23,59)
timedelta = timedelta(hours=1)
time = startime
hrlist = []
while (time <= endtime):
    hrlist.append(time.hour)
    time += timedelta

merge = merge.assign(hrlist=(['hour'],hrlist))
avg = merge.groupby("hrlist").mean()

avg.to_netcdf("./{0}_{1}-{2}_pcp_sum.nc").format(seasonstring,yarlist[0],yearlist[-1])

#Then plot:
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import cartopy
import cartopy.crs as ccrs
import netCDF4 as nc

start_year = 2010
end_year = 2019
seasonstring = 'MAM'
#seasonstring = 'DJF'

readpath = './{0}_{1}-{2}_pcp_sum.nc'.format(seasonstring, start_year, end_year)

f = nc.Dataset(readpath)

#plot lat, lon, precip
lat = f.variables['latitude']
lon = f.variables['longitude']
TP = f.variables['pcp_total_rate']

#set up axes
fig, axes = plt.subplots(nrows=3, ncols=4, dpi=300, figsize=(20, 10), subplot_kw={'projection': ccrs.PlateCarree()})
test1 = plt.setp(axes[-1, :], xlabel='Longitude (°E)')
test2 = plt.setp(axes[:, 0], ylabel='Latitude (°N)')
pcp_cus_colors = ['white', 'cyan', 'deepskyblue', 'green', 'yellow', 'orange', 'crimson', 'mediumorchid']
pcp_cus_cmap = mpl.colors.ListedColormap(pcp_cus_colors)

#"WIB" is converting from UTC to local time (WIB)
for hour in range(0, 24, 2):
    TotalPrecip = TP[hour, :, :]
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
    cs = ax.pcolormesh(lon, lat, TotalPrecip, transform=ccrs.PlateCarree(), vmax=2, vmin=0, cmap=pcp_cus_cmap)
    ax.set_title('{0}00 WIB'.format(WIB))
cb = plt.colorbar(cs, orientation='vertical', shrink=0.75, ax=axes)
cb.set_label('Total Precipitation (mm/h)', fontsize=15)
plt.suptitle('Precipitation Daily Average 2010-2019', x=0.45, y=0.95, fontsize=25)
plt.show()
