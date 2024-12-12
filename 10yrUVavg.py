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
