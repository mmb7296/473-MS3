import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import cartopy
import cartopy.crs as ccrs
import netCDF4 as nc

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
