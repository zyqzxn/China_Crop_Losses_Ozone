## read all files out
import glob
dir = '/Volumes/HD_Li/MP/MP_Data/O3/O3_2013/FGEOS_C4MOZ_L40CN_CEDS_MEIC.cam.h2' files = sorted(glob.glob(dir + '*.nc'))
files
import netCDF4 as nc
import numpy as np
fp = nc.Dataset('/Volumes/HD_Li/MP/MP_Data/O3/O3_2013/FGEOS_C4MOZ_L40CN_CEDS_MEIC.cam data = fp['O3_1'][0:12,55,:,:]*(10**9)
aot = np.where(data<40, np.nan, data-40)
aotsum = np.sum(aot, axis = 0) aotsum
### Calculation                .
## read all files out
import glob
dir = '/Volumes/HD_Li/MP/MP_Data/O3/O3_2013/FGEOS_C4MOZ_L40CN_CEDS_MEIC.cam.h2' files = sorted(glob.glob(dir + '*.nc'))
files
import netCDF4 as nc
import numpy as np
fp = nc.Dataset('/Volumes/HD_Li/MP/MP_Data/O3/O3_2013/FGEOS_C4MOZ_L40CN_CEDS_MEIC.cam data = fp['O3_1'][0:12,55,:,:]*(10**9)
aot = np.where(data<40, np.nan, data-40)
aotsum = np.sum(aot, axis = 0) aotsum
#define a ozone matrix function:
import numpy as np def AOT40(data):
sumthing = []
for n in range(0,12):
    aot = data[n,:,:]#daytime_data[n,:,:]
    aot1 = aot - 40
    #print(aot1)
    # replace negative values as NAN or 0; 
    aot1_temp = np.where(aot1<0, np.nan, aot1) 
    sumthing. append(aot1_temp)
sum1 = np.array(sumthing) 
sum2 = np.sum(sum1,axis = 0) 
return sum2
                
#calcAOT40 and store in an array
import netCDF4 as nc
import numpy as np
#create an array to store the daily AOT40 X = []
for i in range(180,182):#len(files)):
    file = files[i]
    fp = nc.Dataset(file, 'r')
    # Extract data from NetCDF file
    lats = fp['lat'][:]
    lons = fp['lon'][:]
    times = fp['time']
    level = fp['lev'][55]
    data = fp['O3_1'][:,55,:,:] * (10**9) #convert mol/mol to ppbv daytime_data = data[0:12, :, :] #local daytime hours (8:00–19:59) # local daytime; Green time: 8 hours difference #print(daytime_data.shape)
    #print(data.shape)
    daily_AOT40 = AOT40(daytime_data)#print(daily_AOT40)
    X. append(daily_AOT40)
               
#####before writing in new netCDF file,check the metadata info of the original files,
file = '/Volumes/HD_Li/MP/MP_Data/O3/O3_2011/FGEOS_C4MOZ_L40CN_CEDS_MEIC.cam.h2.2013- fp = nc.Dataset(file, 'r')
fp['lat']
 
#####regrid the files into new resolutions               
%%bash
#cdo -griddes /Volumes/HD_Li/MP/MP_Data/Ozone/2017/FGEOS_C4MOZ_L40CN_CEDS_MEIC.cam.h2.
cdo -griddes /Volumes/HD_Li/MP/MP_Data/O3/O3_2010/FGEOS_C4MOZ_L40
