import xarray 
import matplotlib.pyplot as plt
import matplotlib as mpl
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
fig = plt.figure(figsize=[12,5])
ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
diff.plot(ax=ax,vmin=-0.02, vmax=0.02, cmap='Blues', transform=ccrs.PlateCarree())
ax.coastlines()
ax.gridlines()
plt.show()

#diff.plot()
#plt.show()

