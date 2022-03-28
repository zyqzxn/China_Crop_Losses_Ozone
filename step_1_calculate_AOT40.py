##Create a new file to store the averaged data #
import netCDF4 as nc
from netCDF4 import Dataset
import numpy as np
import os
#create a file
newfile = '/Users/DyLi/MP/MP_Data/AOT40_Daily/aot40_2017_05.nc'
#Open the file as writable
fp_new = nc.Dataset(newfile,'w', format = 'NETCDF4')
#Create dimensions
fp_new. createDimension('time',N o n e )
fp_new. createDimension('lat',size= lats. size) fp_new. createDimension('lon',size= lons. size)
#Create variables
fp_new. createVariable("time","f8",("time",))
lat = fp_new.createVariable("lat","f8",("lat",)) lat.standard_name = 'latitude'
lat.long_name = 'latitude'
lat.units = 'degrees_north'
lat.axis = 'Y'
lon = fp_new.createVariable("lon","f8",("lon",)) lon.standard_name = 'longitude'
lon.long_name = 'longitude'
lon.units = 'degrees_east'
lon.axis = 'X'
aot40 = fp_new.createVariable("aot40","f8",("time","lat","lon")) aot40.standard_name = 'aot40'
aot40.long_name = 'O3 concentration in AOT40'
aot40.units = 'ppm'
aot40.comment = "center of grid cell"
fp_new['lat'][:] = lats[:] fp_new['lon'][:] = lons[:]
fp_new['aot40'][:] = AOT40_yr[:]
#fill days in
import datetime
day = np.arange(len(fp_new['time']))
fp_new['time'].units = 'days since 2017-1-1 00:00:00.0'
fp_new['time'].calendar = 'gregorian'
#dates = nc.num2date(day,units= fp_new['time'].units,calendar=times.calendar) #fp_new['time'][:] = nc.date2num(dates, units= fp_new['time'].units,calendar=times.cal 
fp_new['time'][:]= day[:]

#Adding attributes
fp_new.description = "Daily AOT40 in ppm, 2017" #Close the file
fp_new. close()

# Check the file info
file = '/Users/DyLi/MP/MP_Data/AOT40_Daily/aot40_2017.nc' fp = nc.Dataset(file, 'r')
fp['aot40']

# Check up metadata of the new file
%%bash
cdo -showname /Users/DyLi/MP/MP_Data/AOT40_Daily/aot40_2017.nc cdo -griddes /Users/DyLi/MP/MP_Data/AOT40_Daily/aot40_2017.nc
