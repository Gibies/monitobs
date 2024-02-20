import xarray 
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
#pyplot.switch_backend('agg')
from mpl_toolkits.basemap import Basemap
import cartopy
import cartopy.crs as ccrs

data1="/home/meena/q_20230610_00_ctl"
data2="/home/meena/q_20230610_00_exp"

a = xarray.open_dataset(data1).sel(hybrid_ht=0,longitude=slice(60,100),latitude=slice(0,40))
b = xarray.open_dataset(data2).sel(hybrid_ht=0,longitude=slice(60,100),latitude=slice(0,40))
var1 = a['q']
var2 = b['q']
diff = var2 - var1
print(diff)
fig = pyplot.figure(figsize=[12,5])
tmpdir='/home/umreg/ShortJobs/SWW/tmp/'
m = Basemap(projection='cyl',llcrnrlat=-85,urcrnrlat=85,llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='l')
ax = pyplot.axes(projection=ccrs.PlateCarree(central_longitude=0))
diff.plot(ax=ax,vmin=-0.02, vmax=0.02, cmap='Blues', transform=ccrs.PlateCarree())
#ax.coastlines()
#ax.gridlines()
m.drawcoastlines()
m.drawcountries()
m.readshapefile(tmpdir+'IND_WHOLE', 'IND_WHOLE', color='black',linewidth=0.5,  drawbounds=True)
pyplot.show()
pyplot.savefig('/scratch/${USER}/trial.png')

#diff.plot()
#plt.show()

