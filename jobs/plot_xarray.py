import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy
import cartopy.crs as ccrs

a = xr.open_dataset('/home/meena/q_20230610_00_ctl')
b = xr.open_dataset('/home/meena/q_20230610_00_exp')
#print(a)
#print(b)
a1=a.sel(hybrid_ht=0,longitude=slice(60,100),latitude=slice(0,40))
b1=b.sel(hybrid_ht=0,longitude=slice(60,100),latitude=slice(0,40))
#a1=a.sel(hybrid_ht=0)
#b1=b.sel(hybrid_ht=0)
var1 = a1['q']
var2 = b1['q']
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

