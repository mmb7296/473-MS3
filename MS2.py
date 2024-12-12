#Average precipitation, U and V winds over 24 hours
#Precipitation:

from datetime import datetime
from datetime import timedelta
import xarray

pcpdict = {}
yearlist = [2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019]

for year in yearlist:
    janpath = 'Hourly_J_{0}_pcp_sum.nc'.format(year)
    febpath = 'Hourly_F_{0}_pcp_sum.nc'.format(year)
    marchpath = 'Hourly_M_{0}_pcp_sum.nc'.format(year)
    aprilpath = 'Hourly_A_{0}_pcp_sum.nc'.format(year)
    maypath = 'Hourly_Ma_{0}_pcp_sum.nc'.format(year)
    junepath = 'Hourly_June_{0}_pcp_sum.nc'.format(year)
    julypath = 'Hourly_July_{0}_pcp_sum.nc'.format(year)
    augustpath = 'Hourly_Aug_{0}_pcp_sum.nc'.format(year)
    septpath = 'Hourly_S_{0}_pcp_sum.nc'.format(year)
    octpath = 'Hourly_O_{0}_pcp_sum.nc'.format(year)
    novpath = 'Hourly_N_{0}_pcp_sum.nc'.format(year)
    decpath = 'Hourly_D_{0}_pcp_sum.nc'.format(year)
    
    janfile = xarray.open_dataset(janpath)
    febfile = xarray.open_dataset(febpath)
    marchfile = xarray.open_dataset(marchpath)
    aprilfile = xarray.open_dataset(aprilpath)
    mayfile = xarray.open_dataset(maypath)
    junefile = xarray.open_dataset(junepath)
    julyfile = xarray.open_dataset(julypath)
    augfile = xarray.open_dataset(augustpath)
    septfile = xarray.open_dataset(septpath)
    octfile = xarray.open_dataset(octpath)
    novfile = xarray.open_dataset(novpath)
    decfile = xarray.open_dataset(decpath)
    
    pcpdict["{0}".format(year)] = xarray.concat([janfile,febfile,marchfile,aprilfile,mayfile,junefile,julyfile,augfile,septfile,octfile,novfile,decfile], dim="hour")

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

avg.to_netcdf("./{0}-{1}_pcp_sum.nc").format(yearlist[0],yearlist[-1])







#U and V-winds
#First, taking the winds over 10 years:

import netCDF4 as nc
import numpy as np
node_year = 2010 #This coce, as well as the V-wind code, is run for each year from 2010 to 2019
                 #and produces a new file each time.

startime = datetime(node_year,1,1,0,0)
endtime = datetime(node_year,12,31,23,59)

UVtime = startime
U_path = '/scratch/07687/samlg/ERA5_3D/{1}/e5.oper.an.pl.128_131_u.ll025uv.{0}00_{0}23.nc'.format(UVtime.strftime('%Y%m%d'), UVtime.strftime('%Y%m'))
#Source of file in TACC supercomputer

#U
Ufile = nc.Dataset(U_path)
latitude = Ufile.variables['latitude'][310:430]
longitude = Ufile.variables['longitude'][350:500]
level = Ufile.variables['level'][:]
time = Ufile.variables['time'][:]
time = len(time)
level = len(level)
lenlatitude = len(latitude)
lenlongitude = len(longitude)
print("Variables taken in")

Ufile.close()

i = 0
U_wind_total = np.zeros(92, dtype=object)
while (UVtime <= endtime):
    #variables that repeat
    for i in range(0, 92):
        U_path = '/scratch/07687/samlg/ERA5_3D/{1}/e5.oper.an.pl.128_131_u.ll025uv.{0}00_{0}23.nc'.format(UVtime.strftime('%Y%m%d'), UVtime.strftime('%Y%m'))
        V_path = '/scratch/07687/samlg/ERA5_3D/{1}/e5.oper.an.pl.128_132_v.ll025uv.{0}00_{0}23.nc'.format(UVtime.strftime('%Y%m%d'), UVtime.strftime('%Y%m'))
        Ufile = nc.Dataset(U_path)
        U_wind = Ufile.variables['U'][:,:,310:430,350:500]
        U_wind_total[i] = U_wind
        Ufile.close()
        print("Day {0} U wind complete".format(i+1))
        delta = datetime.timedelta(days=1)
        UVtime += delta
            
meanU = np.mean(U_wind_total, axis=0)
print("Mean complete")
    
#send this info back to netcdf4
Usavepath = '/work2/10049/mbennett1/stampede3/UMEAN_R{0}.nc'.format(node_year)
Uncfile = nc.Dataset(Usavepath, 'w', format = 'NETCDF4')
Uncfile.createDimension('time', time)
Uncfile.createDimension('level', level)
Uncfile.createDimension('lat', lenlatitude)
Uncfile.createDimension('lon', lenlongitude)
Uncfile.createVariable('u wind', np.float32, ('time', 'level', 'lat', 'lon'))
Uncfile.createVariable('latitude', np.float32, ('lat'))
Uncfile.createVariable('longitude', np.float32, ('lon'))
Uncfile.variables['u wind'][:,:,:,:] = meanU[:,:,:,:]
Uncfile.variables['latitude'][:] = latitude[:]
Uncfile.variables['longitude'][:] = longitude[:]
Uncfile.close()

UVtime = startime
V_path = '/scratch/07687/samlg/ERA5_3D/{1}/e5.oper.an.pl.128_132_v.ll025uv.{0}00_{0}23.nc'.format(UVtime.strftime('%Y%m%d'), UVtime.strftime('%Y%m'))

#V
Vfile = nc.Dataset(V_path)
latitude = Vfile.variables['latitude'][310:430]
longitude = Vfile.variables['longitude'][350:500]
level = Vfile.variables['level'][:]
time = Vfile.variables['time'][:]
time = len(time)
level = len(level)
lenlatitude = len(latitude)
lenlongitude = len(longitude)
print("Variables taken in")

Vfile.close()

i = 0
V_wind_total = np.zeros(92, dtype=object)
while (UVtime <= endtime):
    #variables that repeat
    for i in range(0, 92):
        U_path = '/scratch/07687/samlg/ERA5_3D/{1}/e5.oper.an.pl.128_131_u.ll025uv.{0}00_{0}23.nc'.format(UVtime.strftime('%Y%m%d'), UVtime.strftime('%Y%m'))
        V_path = '/scratch/07687/samlg/ERA5_3D/{1}/e5.oper.an.pl.128_132_v.ll025uv.{0}00_{0}23.nc'.format(UVtime.strftime('%Y%m%d'), UVtime.strftime('%Y%m'))
        Vfile = nc.Dataset(V_path)
        V_wind = Vfile.variables['V'][:,:,310:430,350:500]
        V_wind_total[i] = V_wind
        print("Day {0} V wind complete".format(i+1))
        Vfile.close()
        delta = datetime.timedelta(days=1)
        UVtime += delta
            
meanV = np.mean(V_wind_total, axis=0)
print("Mean complete")
    
#send this info back to netcdf4
Vsavepath = '/work2/10049/mbennett1/stampede3/VMEAN_R{0}.nc'.format(node_year)
Vncfile = nc.Dataset(Vsavepath, 'w', format = 'NETCDF4')
Vncfile.createDimension('time', time)
Vncfile.createDimension('level', level)
Vncfile.createDimension('lat', lenlatitude)
Vncfile.createDimension('lon', lenlongitude)
Vncfile.createVariable('v wind', np.float32, ('time', 'level', 'lat', 'lon'))
Vncfile.createVariable('latitude', np.float32, ('lat'))
Vncfile.createVariable('longitude', np.float32, ('lon'))
Vncfile.variables['v wind'][:,:,:,:] = meanV[:,:,:,:]
Vncfile.variables['latitude'][:] = latitude[:]
Vncfile.variables['longitude'][:] = longitude[:]
Vncfile.close()

#Then, averaging together:

U = {}
yearlist = [2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019]
for year in yearlist:
    Upath = './UVWind/UMEAN_R{0}.nc'.format(year)
    Ufile = xarray.open_dataset(Upath)
    
    U["{0}".format(year)] = Upath

Umerge = xarray.concat([U["2010"], U["2011"], U["2012"], U["2013"], U["2014"], U["2015"], U["2016"], U["2017"], U["2018"], U["2019"]], dim="time")
startime = datetime(yearlist[0],1,1,0,0)
endtime = datetime(yearlist[-1],12,31,23,59)
timedelta = timedelta(hours=1)
time = startime
hrlist = []
while (time <= endtime):
    hrlist.append(time.hour)
    time += timedelta
Umerge = Umerge.assign(hrlist=(['hour'],hrlist))
Uavg = Umerge.groupby("hrlist").mean()

Uavg.to_netcdf("./UMean_{0}-{1}.nc").format(yearlist[0],yearlist[-1])

V = {}
for year in yearlist:
    Vpath = './VVWind/VMEAN_R{0}.nc'.format(year)
    Vfile = xarray.open_dataset(Vpath)
    
    V["{0}".format(year)] = Vpath

Vmerge = xarray.concat([V["2010"], V["2011"], V["2012"], V["2013"], V["2014"], V["2015"], V["2016"], V["2017"], V["2018"], V["2019"]], dim="time")
startime = datetime(yearlist[0],1,1,0,0)
endtime = datetime(yearlist[-1],12,31,23,59)
timedelta = timedelta(hours=1)
time = startime
hrlist = []
while (time <= endtime):
    hrlist.append(time.hour)
    time += timedelta
Vmerge = Vmerge.assign(hrlist=(['hour'],hrlist))
Vavg = Vmerge.groupby("hrlist").mean()

Vavg.to_netcdf("./VMean_{0}-{1}.nc").format(yearlist[0],yearlist[-1])



#Plot average precipitation over 24 hours
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import cartopy
import cartopy.crs as ccrs
import netCDF4 as nc

start_year = 2010
end_year = 2019

readpath = './{0}-{1}_pcp_sum.nc'.format(start_year, end_year)

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

#Plot average winds over 24 hours at key levels in the atmosphere
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













#**Plot other variables' averages over 24 hours**
#Saving temperature and cloud cover is very similar to saving the u and v winds, given the similarities in ERA5 formatting:
#Temp
node_year = 2010 #can alter freely

startime = dt(node_year, 1, 1, 0, 0)
endtime = dt(node_year, 12, 31, 23, 59)

Ttime = startime
T_path = '/scratch/07687/samlg/ERA5_3D/{1}/e5.oper.an.pl.128_130_t.ll025sc.{0}00_{0}23.nc'.format(Ttime.strftime('%Y%m%d'), Ttime.strftime('%Y%m'))

#U
Tfile = nc.Dataset(T_path)
latitude = Tfile.variables['latitude'][310:430]
longitude = Tfile.variables['longitude'][350:500]
level = Tfile.variables['level'][:]
time = Tfile.variables['time'][:]
time = len(time)
level = len(level)
lenlatitude = len(latitude)
lenlongitude = len(longitude)
print("Variables taken in")

Tfile.close()

i = 0
T_total = np.zeros(92, dtype=object)
while (Ttime <= endtime):
    #variables that repeat
    for i in range(0, 92):
        T_path = '/scratch/07687/samlg/ERA5_3D/{1}/e5.oper.an.pl.128_130_t.ll025sc.{0}00_{0}23.nc'.format(Ttime.strftime('%Y%m%d'), Ttime.strftime('%Y%m'))
        Tfile = nc.Dataset(T_path)
        Temp = Tfile.variables['T'][:,:,310:430,350:500]
        T_total[i] = Temp
        Tfile.close()
        print("Day {0} T complete".format(i+1))
        delta = datetime.timedelta(days=1)
        Ttime += delta

meanT = np.mean(T_total, axis=0)
print("Mean complete")

#send this info back to netcdf4
Tsavepath = './TempMEAN_{0}.nc'.format(node_year)
Tncfile = nc.Dataset(Tsavepath, 'w', format = 'NETCDF4')
Tncfile.createDimension('time', time)
Tncfile.createDimension('level', level)
Tncfile.createDimension('lat', lenlatitude)
Tncfile.createDimension('lon', lenlongitude)
Tncfile.createVariable('temp', np.float32, ('time', 'level', 'lat', 'lon'))
Tncfile.createVariable('latitude', np.float32, ('lat'))
Tncfile.createVariable('longitude', np.float32, ('lon'))

#Cloud cover

node_year = 2010

startime = dt(node_year, 1, 1, 0, 0)
endtime = dt(node_year, 12, 31, 23, 59)

CCtime = startime
CC_path = '/scratch/07687/samlg/ERA5_3D/{1}/e5.oper.an.pl.128_248_cc.ll025sc.{0}00_{0}23.nc'.format(CCtime.strftime('%Y%m%d'), CCtime.strftime('%Y%m'))

#U
CCfile = nc.Dataset(CC_path)
latitude = CCfile.variables['latitude'][310:430]
longitude = CCfile.variables['longitude'][350:500]
level = CCfile.variables['level'][:]
time = CCfile.variables['time'][:]
time = len(time)
level = len(level)
lenlatitude = len(latitude)
lenlongitude = len(longitude)
print("Variables taken in")

CCfile.close()

i = 0
CC_total = np.zeros(92, dtype=object)
while (CCtime <= endtime):
    #variables that repeat
    for i in range(0, 92):
        CC_path = '/scratch/07687/samlg/ERA5_3D/{1}/e5.oper.an.pl.128_248_cc.ll025sc{0}00_{0}23.nc'.format(CCtime.strftime('%Y%m%d'), CCtime.strftime('%Y%m'))
        CCfile = nc.Dataset(CC_path)
        CC = CCfile.variables['CC'][:,:,310:430,350:500]
        CC_total[i] = CC
        CCfile.close()
        print("Day {0} CC complete".format(i+1))
        delta = datetime.timedelta(days=1)
        CCtime += delta

meanCC = np.mean(CC_total, axis=0)
print("Mean complete")

#send this info back to netcdf4
CCsavepath = './CCMEAN_{0}.nc'.format(node_year)
CCncfile = nc.Dataset(CCsavepath, 'w', format = 'NETCDF4')
CCncfile.createDimension('time', time)
CCncfile.createDimension('level', level)
CCncfile.createDimension('lat', lenlatitude)
CCncfile.createDimension('lon', lenlongitude)
CCncfile.createVariable('CC', np.float32, ('time', 'level', 'lat', 'lon'))
CCncfile.createVariable('latitude', np.float32, ('lat'))
CCncfile.createVariable('longitude', np.float32, ('lon'))
#likewise, averaging is similar
T = {}
yearlist = [2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019]
for year in yearlist:
    Tpath = './TempMEAN_{0}.nc'.format(year)
    Tfile = xarray.open_dataset(Tpath)
    
    T["{0}".format(year)] = Tpath

Tmerge = xarray.concat([T["2010"], T["2011"], T["2012"], T["2013"], T["2014"], T["2015"], T["2016"], T["2017"], T["2018"], T["2019"]], dim="time")
startime = datetime(yearlist[0],1,1,0,0)
endtime = datetime(yearlist[-1],12,31,23,59)
timedelta = timedelta(hours=1)
time = startime
hrlist = []
while (time <= endtime):
    hrlist.append(time.hour)
    time += timedelta
Tmerge = Tmerge.assign(hrlist=(['hour'],hrlist))
Tavg = Tmerge.groupby("hrlist").mean()

Tavg.to_netcdf("./TMean_{0}-{1}.nc").format(yearlist[0],yearlist[-1])

CC = {}
for year in yearlist:
    CCpath = './CCMEAN_{0}.nc'.format(year)
    CCfile = xarray.open_dataset(CCpath)
    
    CC["{0}".format(year)] = CCpath

CCmerge = xarray.concat([CC["2010"], CC["2011"], CC["2012"], CC["2013"], CC["2014"], CC["2015"], CC["2016"], CC["2017"], CC["2018"], CC["2019"]], dim="time")
startime = datetime(yearlist[0],1,1,0,0)
endtime = datetime(yearlist[-1],12,31,23,59)
timedelta = timedelta(hours=1)
time = startime
hrlist = []
while (time <= endtime):
    hrlist.append(time.hour)
    time += timedelta
CCmerge = CCmerge.assign(hrlist=(['hour'],hrlist))
CCavg = CCmerge.groupby("hrlist").mean()

CCavg.to_netcdf("./CCMean_{0}-{1}.nc").format(yearlist[0],yearlist[-1])

#Plotting
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

readpath = './CCMean_{0}-{1}.nc'.format(yearlist[0],yearlist[-1])

f = nc.Dataset(readpath)

#plot lat, lon, precip
lat = f.variables['latitude']
lon = f.variables['longitude']
CC = f.variables['CC']
level = 30 #650 hPA

#set up axes
fig, axes = plt.subplots(nrows=3, ncols=4, dpi=300, figsize=(20, 10), subplot_kw={'projection': ccrs.PlateCarree()})
test1 = plt.setp(axes[-1, :], xlabel='Longitude (°E)')
test2 = plt.setp(axes[:, 0], ylabel='Latitude (°N)')
cloud_colors = ['white', 'cyan', 'deepskyblue', 'green', 'yellow', 'orange', 'crimson', 'mediumorchid']
cloud_map = mpl.colors.ListedColormap(cloud_colors)

#"WIB" is converting from UTC to local time (WIB)
for hour in range(0, 24, 2):
    cloudcover = CC[hour, level :, :]
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
    cs = ax.pcolormesh(lon, lat, cloudcover, transform=ccrs.PlateCarree(), cmap=cloud_map)
    ax.set_title('{0}00 WIB'.format(WIB))
cb = plt.colorbar(cs, orientation='vertical', shrink=0.75, ax=axes)
cb.set_label('Cloud Cover (fraction)', fontsize=15)
plt.suptitle('Cloud Cover Daily Average 2010-2019', x=0.45, y=0.95, fontsize=25)
plt.show()











#getting MAM and DJF averages

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
