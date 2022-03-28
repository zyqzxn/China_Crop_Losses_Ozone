#check the CP data info
%%bash
cdo -showname /Users/DyLi/MP/MP_Data/Crop/WHEAT_production_total_rain+irrig_doublePrec 
cdo -griddes /Users/DyLi/MP/MP_Data/Crop/WHEAT_production_total_rain+irrig_doublePreci

#regrid
%%writefile /Users/DyLi/MP/MP_Data/AOT40_Daily/grid.cdo
gridtype= lonlat         
xsize = 70
ysize = 45
xname= lon
xlongname= "longitude"
xunits= "degrees_east" 
yname= lat
ylongname = "latitude"
yunits= "degrees_north" 
xfirst= 70.5
xinc=1
yfirst = 10.5
yinc=1

#check
%%bash
cdo -selname,aot40 -remapnn,/Users/DyLi/MP/MP_Data/AOT40_Daily/grid.cdo /Users/DyLi/

#check - before regridding
file = '/Users/DyLi/MP/MP_Data/AOT40_Daily/aot40_2017.nc' fp = nc.Dataset(file, 'r')
fp['aot40']. shape

#check - after regridding
file2 = '/Users/DyLi/MP/MP_Data/AOT40_Daily/aot40_2017_upscale.nc' fp2 = nc.Dataset(file2, 'r')
fp2['aot40'].shape
