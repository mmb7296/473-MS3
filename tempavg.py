import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import cartopy
import cartopy.crs as ccrs
import netCDF4 as nc

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
