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
