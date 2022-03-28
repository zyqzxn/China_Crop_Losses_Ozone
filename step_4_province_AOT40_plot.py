import netCDF4 as nc import numpy as np
import xarray as xray
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap from mpl_toolkits.basemap import maskoceans
import geopandas as gdp

## read all files out
import glob
#dir = '/Users/DyLi/MP/MP_Data/AOT40_Daily/FGEOS_C4MOZ_L40CN_CEDS_MEIC.cam.h2' #files = sorted(glob.glob(dir + '*.nc'))
dir = "/Users/DyLi/MP/MP_Data/MP_Agri_Ozone/FGEOS_C4MOZ_L40CN_CEDS_MEIC.cam.h2." files = sorted(glob.glob(dir + '*.zyq.nc'))
files

from rasterio import features
from affine import Affine
def transform_from_latlon(lat, lon):
lat = np.asarray(lat)
lon = np.asarray(lon)
trans = Affine.translation(lon[0], lat[0])
scale = Affine.scale(lon[1] - lon[0], lat[1] - lat[0]) return trans * scale
def rasterize(shapes, coords, latitude='latitude', longitude='longitude', fill=np.nan, **kwargs):
"""Rasterize a list of (geometry, fill_value) tuples onto the given xray coordinates. This only works for 1d latitude and longitude arrays.
"""
transform = transform_from_latlon(coords[latitude], coords[longitude]) out_shape = (len(coords[latitude]), len(coords[longitude]))
raster = features.rasterize(shapes, out_shape=out_shape,
fill=fill, transform=transform,
dtype=float, **kwargs)
spatial_coords = {latitude: coords[latitude], longitude: coords[longitude]}
return xray.DataArray(raster, coords=spatial_coords, dims=(latitude, longitude))

prov = gdp.read_file('/Users/DyLi/MP/MP_Data/Map/gadm36_CHN_shp/gadm36_CHN_1.shp') 
CHN = '/Users/DyLi/MP/MP_Data/Map/'
gdf2 = gdp.read_file('/Users/DyLi/MP/MP_Data/Map/gadm36_HKG_shp/gadm36_HKG_0.shp') 
gdf3 = gdp.read_file("/Users/DyLi/MP/MP_Data/Map/gadm36_MAC_shp/gadm36_MAC_0.shp") 
gdf4 = gdp.read_file('/Users/DyLi/MP/MP_Data/Map/gadm36_TWN_shp/gadm36_TWN_0.shp')

provs = prov.append([gdf2,gdf3,gdf4],ignore_index=True,sort=True) 
provs.loc[31, "NAME_1"] = "HongKong"
provs.loc[32, "NAME_1"] = "Macao"
provs.loc[33, "NAME_1"] = "Taiwan"
#provs = geopandas.read_file('/Users/DyLi/MP/MP_Data/Map/gadm36_CHN_shp/gadm36_CHN_1.s
provs_ids = {k: i for i, k in enumerate(provs.NAME_1)}
#shapes = [(shape, n) for n, shape in enumerate(provs.geometry)]
shapes = zip(provs.geometry, range(len(provs)))
#ds = xray.open_dataset('/Users/DyLi/MP/MP_Data/AOT40_Daily/aot40_2010_05.nc')
ds = xray.open_dataset('/Users/DyLi/MP/MP_Data/MP_Agri_Ozone/FGEOS_C4MOZ_L40CN_CEDS_M ds['provs'] = rasterize(shapes, ds.coords, longitude='lon', latitude='lat')

####Calc AOT40 by province
import pandas as pd
aotSum = pd.DataFrame(provs_ids.items(),columns=['Province', 'Order']) #dfWheat
#dfWheat["2010_AOT40_ppm"]= 0
for y in range (8):
    ds1 = xray.open_dataset(files[y]) for p in range(34):
    #a = ds1.aot40.where(ds.provs == p) #extract data of province #p
    a = ds1.AOT40_O3.where(ds.provs == p)
    # sum of a whole province aot level during Wheat growing season b = np.sum(a[:,:,:],axis = 0)
    #calculate the pixel-average AOT40 within a province
    c = np.sum(b)/np.count_nonzero(b)
    #calculate the pixel-average AOT40 within a province #dfWheat["2010_AOT40_ppm"][p]= b
    aotSum.loc[p, "201%d_AOT40_ppb"%y] = c
    #dfWheat["2010_RYL"][p]= b
    
#save to excel
aotSum. to_excel("aotSum.xlsx")

#####plot national AOT40
coun1 = gdp.read_file('/Users/DyLi/MP/MP_Data/Map/gadm36_CHN_shp/gadm36_CHN_0.shp') 
coun = coun1.append([gdf2,gdf3,gdf4],ignore_index=True,sort = True)
shapes = zip(coun.geometry, range(len(coun)))
ds = xray.open_dataset('/Users/DyLi/MP/MP_Data/AOT40_Daily/aot40_2010_05.nc') 
ds['country'] = rasterize(shapes, ds.coords, longitude='lon', latitude='lat')
mask = np.logical_or(ds["country"] == 0, ds["country"] == 3)

fig = plt.figure(figsize=(15,25))

vmin= 0
vmax= 64
cmap= plt. get_cmap('jet')
cmap. set_bad(color= 'w')
norm = mpl.colors.Normalize(vmin = vmin, vmax = vmax) levels = np.linspace(vmin,vmax,9)

for i in range(len(files)): file = files[i]
    fp = xray.open_dataset(file) lons = fp['lon'][5:135]
    lats = fp['lat'][10:]
    lon_0 = lons.mean() lat_0 = lats.mean()
    a = fp.AOT40_O3.where(mask)#extract data of province
    # sum of a whole province aot level during Wheat growing season
    data= np.sum(a,axis=0)
    data1 = np.array(data).reshape((90,140))
    data1[data1==0] = np.nan
    data2 = data1[10:,5:135]/(10**3) # sum of a provincal aot level during Wheat grow

    m = Basemap(projection='cyl',
                lat_0=lat_0, lon_0=lon_0,
                llcrnrlat= np. min(lats. values),urcrnrlat= np. max(lats. values),\ llcrnrlon= np. min(lons. values),urcrnrlon= np. max(lons. values),\ rsphere= 6371200,area_thresh= 10000)

    lon, lat = np.meshgrid(lons, lats) xi, yi = m(lon, lat)
    plt.subplot(4,2, i+1)
    cs2 = m.contourf(xi,yi,data2,cmap=cmap, norm = norm, levels = levels)
    # Add Grid Lines
    m.drawparallels(np.arange(-80., 81., 20.), labels=[1,0,0,0], fontsize=10) m.drawmeridians(np.arange(-180., 181., 20.), labels=[0,0,0,1], fontsize=10)
    # Add Province and Country Boundaries
    m. readshapefile(dir1,'Province',drawbounds= T r u e ) m. readshapefile(dir2,'Province',drawbounds= T r u e ) m. readshapefile(dir3,'Province',drawbounds= T r u e ) m. readshapefile(dir4,'Province',drawbounds= T r u e )
    # Add Colorbar
    cbar = m.colorbar(cs2,location='bottom', pad="10%")
    # Add Title
    plt.title('AOT40(ppm h), 201%d'%i)
    
plt. savefig('AOT40')
plt. tight_layout()
plt. show()


