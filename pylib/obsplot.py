#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 12:56:36 2019

@author: gibies
"""
from __future__ import print_function
import os,sys
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
PKGHOME=os.path.dirname(CURR_PATH)
OBSLIB=os.environ.get('OBSLIB',PKGHOME+"/pylib")
sys.path.append(OBSLIB)
OBSDIC=os.environ.get('OBSDIC',PKGHOME+"/pydic")
sys.path.append(OBSDIC)
OBSNML=os.environ.get('OBSNML',PKGHOME+"/nml")
sys.path.append(OBSNML)
SUBTYPNML=OBSNML+"/obs_subtype.nml"
import obslib
import domaindic
import numpy
import pandas
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
pyplot.switch_backend('agg')
import matplotlib.colors as colors
import matplotlib.cm as mplcm
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
#import Ngl

diaglev=int(os.environ.get('GEN_MODE',0))
def errprint(*args, **kwargs):
    if diaglev > 0: print(*args, file=sys.stderr, **kwargs)

def truncate_colormap(cpallet, minval=0.0, maxval=1.0, n=100):
    cmap=mplcm.get_cmap(cpallet)
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(numpy.linspace(minval, maxval, n)))
    return new_cmap

def get_keyrange(dataset,keyfield):
    lblset=set()
    for data in dataset["data"]:
       for keyfieldvalue in data[keyfield]:
          lblset.add(keyfieldvalue)
    lblist=list(lblset)
    return(lblist)

def plot_keyfield(dataset,plotfile,tagmark="",lblst=[],text="",textpos=(0.25, -0.20)):
    colors = ["b","g","y","r"]
    cmap= matplotlib.colors.ListedColormap(colors)
    clevs=range(0,24,6)
    fig = pyplot.figure()
    #colors = (0,0,0)
    area = 1.0      #numpy.pi*
    alpha=0.5
    parallels = numpy.arange(-80.,90,20.)
    meridians = numpy.arange(-180.,180.,30.)
    plot1=pyplot.subplot(212)
    fig= plot_cyl(dataset,fig,plot1,colors,area,alpha,parallels,meridians,tagmark=tagmark,lblst=lblst,text=text,textpos=textpos)
    #######
    pyplot.savefig(plotfile,bbox_inches='tight',dpi=200)
    return(fig)

def plot_latlon(data,plotfile,tagmark="",lblst=[],text="",textpos=(0.25, -0.20),fltrkey="subtype",subtypenml=SUBTYPNML,title=""):
    print(data)
    if fltrkey in data:
    	keylist=data[fltrkey].unique()
    	print(keylist)
    	datalist=[None]*len(keylist)
    	print(datalist)
    	for i,key in enumerate(keylist):
	   print(i,key)
	   datalist[i]=data[data[fltrkey] == key]
	   print(datalist[i])
    else:
	fltrkey=None
	keylist=["data"]
	datalist=[data]
    colors = ["b","g","y","r","violet","pink","purple","magenta"]
    cmap= matplotlib.colors.ListedColormap(colors)
    clevs=range(0,24,6)
    fig = pyplot.figure()
    #colors = (0,0,0)
    area = 1.0      #numpy.pi*
    alpha=0.5
    parallels = numpy.arange(-80.,90,20.)
    meridians = numpy.arange(-180.,180.,30.)
    #parallels = numpy.arange(-80.,90,20.)
    #meridians = numpy.arange(-180.,180.,30.)
    #######
    #plot1=pyplot.subplot(121)
    #fig = plot_ortho(data,fig,plot1,colors,area,alpha,parallels,meridians,polelat=10,polelon=70)
    #######
    #plot2=pyplot.subplot(122)
    #fig = plot_ortho(data,fig,plot1,colors,area,alpha,parallels,meridians,polelat=10,polelon=250)
    #######
    #plot1=pyplot.subplot(211)
    plot3=pyplot.subplot(212)
    if fltrkey == "subtype":
	lblst=[None]*len(keylist)
	for i,key in enumerate(keylist):
	   lblst[i]=obslib.get_subtype_name(subtypenml,key)
    else:
	lblst=keylist
    fig= plot_cyl(datalist,fig,plot3,colors,area,alpha,parallels,meridians,tagmark=tagmark,lblst=lblst,text=text,textpos=textpos,title=title)
    #######
    pyplot.savefig(plotfile,bbox_inches='tight',dpi=300)
    return(fig)

def plot_location(data,plotfile,tagmark="",lblst=[],text="",textpos=(0.25, -0.20)):
    colors = ["b","g","y","r","violet","pink","purple","magenta"]
    cmap= matplotlib.colors.ListedColormap(colors)
    clevs=range(0,24,6)
    fig = pyplot.figure()
    #colors = (0,0,0)
    area = 1.0      #numpy.pi*
    alpha=0.5
    parallels = numpy.arange(-80.,90,20.)
    meridians = numpy.arange(-180.,180.,30.)
    #######
    #plot1=pyplot.subplot(121)
    #fig = plot_ortho(data,fig,plot1,colors,area,alpha,parallels,meridians,polelat=10,polelon=70)
    #######
    #plot2=pyplot.subplot(122)
    #fig = plot_ortho(data,fig,plot1,colors,area,alpha,parallels,meridians,polelat=10,polelon=250)
    #######
    #plot1=pyplot.subplot(211)
    plot3=pyplot.subplot(212)
    fig= plot_cyl(data,fig,plot3,colors,area,alpha,parallels,meridians,tagmark=tagmark,lblst=lblst,text=text,textpos=textpos)
    #######
    pyplot.savefig(plotfile,bbox_inches='tight',dpi=200)
    return(fig)

def plot_merc(data,figure,plot,colors,area,alpha,parallels,meridians):
    plot = Basemap(projection='merc',llcrnrlat=-85,urcrnrlat=85,llcrnrlon=-180,urcrnrlon=180,lat_ts=20,resolution='c')
    plot.drawlsmask(land_color='wheat',ocean_color='lightblue',lakes=True)
    #map.bluemarble(scale=0.5);
    plot.drawcoastlines()
    #map.drawcountries()
    plot.drawparallels(parallels,labels=[True,False,True,False],fontsize=5)
    plot.drawmeridians(meridians,labels=[False,False,False,True],fontsize=5,rotation=45)
    x, y = plot(list(data.Longitude.values), list(data.Latitude.values))
    plot = pyplot.scatter(x,y,s=area,c=colors,alpha=alpha)
    return(figure)

def data_check(data):
        #x=numpy.array(data.Longitude.values)
        #y=numpy.array(data.Latitude.values)
	#xcheck=(x*10%10)
	#ycheck=(y*10%10)
	#data=data[obslib.compare(xcheck,ycheck)]
	data=data[data.Longitude.values != data.Latitude.values]
	print(data)
	return(data)	

def plot_cyl(datalist,figure,plot,colors,area,alpha,parallels,meridians,tagmark="",lblst=[],text="",textpos=(0.25, -0.20),display_count=True,title=""):
    plot = Basemap(projection='cyl', resolution='c', llcrnrlat= -90.,urcrnrlat= 90.,llcrnrlon=-180.,urcrnrlon=180.)
    #plot = Basemap(projection='cyl', resolution='c', llcrnrlat= 0.,urcrnrlat= 40.,llcrnrlon= 60.,urcrnrlon=95.)
    plot.drawlsmask(land_color='wheat',ocean_color='lightblue',lakes=True)
    #map.bluemarble(scale=0.5);
    plot.drawcoastlines()
    plot.drawparallels(parallels,labels=[True,False,True,False],fontsize=5)
    plot.drawmeridians(meridians,labels=[False,False,False,True],fontsize=5,rotation=45)
    count=len(datalist)
    print(count)
    if count < 5 :
       lblxpos=numpy.linspace(0.05, 0.95, num = count, endpoint = False)
       lblypos=[1.05]*count
    else:
	tmparr=list(numpy.linspace(0.05, 0.95, num = int((count//2)+(count%2)), endpoint = False))
	lblxpos=tmparr+tmparr
	lblypos=([1.05]*((count//2)+(count%2)))+([1.10]*int(count//2))
    print(lblxpos,lblypos)
    for idx in range(0,count,1):
        data = datalist[idx]
	data = data_check(data)
        x=numpy.array(data.Longitude.values)
        y=numpy.array(data.Latitude.values)
        if tagmark == "(b)" : print(x,y)
	print(idx,colors)
        plot = pyplot.scatter(x,y,s=area,c=colors[idx],alpha=alpha)
	print(lblst[idx],len(data))
        if display_count==True:
        	lbltxt=str(lblst[idx])+": "+str(len(data))
        	plot = pyplot.annotate(lbltxt,color=colors[idx],fontsize=5, xy=(lblxpos[idx], lblypos[idx]), xycoords='axes fraction')
    plot = pyplot.annotate(tagmark,fontsize=5, xy=(0.01, 1.05), xycoords='axes fraction')
    plot = pyplot.annotate(text,fontsize=5, xy=textpos, xycoords='axes fraction')
    plot = pyplot.annotate(title,fontsize=10, xy=(0.25,1.15), xycoords='axes fraction')
    return(figure)

def plot_ortho(data,figure,plot,colors,area,alpha,parallels,meridians,polelat,polelon):
    plot=Basemap(lat_0=polelat, lon_0=polelon, projection='ortho',resolution='l' )
    plot.drawlsmask(land_color='wheat',ocean_color='lightblue',lakes=True)
    plot.drawcoastlines()
    #plot.readshapefile('/home/gibies/maskfiles/AllIndia/AllIndia', 'AllIndia')
    plot.drawparallels(parallels,labels=[True,False,True,False],fontsize=5)
    plot.drawmeridians(meridians,labels=[False,False,False,True],fontsize=5)
    x, y = plot(list(data.Longitude.values), list(data.Latitude.values))
    plot = pyplot.scatter(x,y,s=area,c=colors,alpha=alpha)
    return(figure)

    #m.readshapefile('/home/gibies/maskfiles/AllIndia/AllIndia', 'AllIndia')
    

def plot_shaded(plotfile,data,title,clevs=range(0,10,1),cpallet="jet",extend="max"):
    cmap=mplcm.get_cmap(cpallet, len(clevs) - 1)
    glat=data.index
    glon=data.columns
    x,y=numpy.meshgrid(glon,glat)
    fig = pyplot.figure()
    map = Basemap(projection='cyl', resolution='c', llcrnrlat= -90.,urcrnrlat= 90.,llcrnrlon=-180.,urcrnrlon=180.)
    #map.bluemarble(scale=0.5);
    map.drawcoastlines()
    map.drawparallels(numpy.arange( -90., 90.,30.),labels=[1,0,0,0],fontsize=5)
    map.drawmeridians(numpy.arange(-180.,180.,30.),labels=[0,0,0,1],fontsize=5,rotation=45)
    #map.drawcountries()
    plot = map.contourf(x,y,data,clevs, cmap=cmap,extend=extend)
    pyplot.title(title)
    #ax = fig.add_subplot(311)
    cbar=fig.colorbar(plot,orientation='vertical')
    cbar.ax.tick_params(labelsize=5)
    pyplot.savefig(plotfile,bbox_inches='tight',dpi=200)

#def plot_button(plotfile,data,title,clevs,cpallet,extend="max"):
#    cmap=mplcm.get_cmap(cpallet, len(clevs) - 1)
#    glat=data.index
#    glon=data.columns
#    x,y=numpy.meshgrid(glon,glat)
#    fig = pyplot.figure()
#    map = Basemap(projection='cyl', resolution='c', llcrnrlat= -90.,urcrnrlat= 90.,llcrnrlon=-180.,urcrnrlon=180.)
#    #map.bluemarble(scale=0.5);
#    map.drawcoastlines()
#    map.drawparallels(numpy.arange( -90., 90.,30.),labels=[1,0,0,0],fontsize=5)
#    map.drawmeridians(numpy.arange(-180.,180.,30.),labels=[0,0,0,1],fontsize=5,rotation=45)
#    #map.drawcountries()
#    plot = map.contour(x,y,data,clevs, cmap=cmap,extend=extend)
#    pyplot.title(title)
#    #ax = fig.add_subplot(311)
#    cbar=fig.colorbar(plot,orientation='vertical')
#    cbar.ax.tick_params(labelsize=5)
#    pyplot.savefig(plotfile,bbox_inches='tight',dpi=200)

def plot_button(plotfile,data,title,clevs=range(0,10,1),cpallet="jet",extend="max"):
    cmap=mplcm.get_cmap(cpallet, len(clevs) - 1)
    glat=data.index
    glon=data.columns
    x,y=numpy.meshgrid(glon,glat)
    fig = pyplot.figure()
    map = Basemap(projection='cyl', resolution='c', llcrnrlat= -90.,urcrnrlat= 90.,llcrnrlon=-180.,urcrnrlon=180.)
    #map.bluemarble(scale=0.5);
    map.drawcoastlines()
    map.drawparallels(numpy.arange( -90., 90.,30.),labels=[1,0,0,0],fontsize=5)
    map.drawmeridians(numpy.arange(-180.,180.,30.),labels=[0,0,0,1],fontsize=5,rotation=45)
    #map.drawcountries()
    plot = map.scatter(x,y,data,clevs, cmap=cmap,extend=extend)
    pyplot.title(title)
    #ax = fig.add_subplot(311)
    cbar=fig.colorbar(plot,orientation='vertical')
    cbar.ax.tick_params(labelsize=5)
    pyplot.savefig(plotfile,bbox_inches='tight',dpi=200)
    return(fig)
    
def scatterplot(plotpath,odbnmlfile,data,cylcdatestr,obstype,varname,long_name,fieldname="Obsvalue",fill=True,extend="max"):
    data=obslib.odb_renamefield(data,odbnmlfile)
    #clevs=obslib.clevgen(long_name,fieldname)
    units=obslib.dataunit(long_name)
    plot_title=obstype.replace("_"," ")+"\n"+long_name.replace("_"," ")+" ("+units+") "+"scatter"+"\n"+cylcdatestr
    plotfile=plotpath+"/"+obstype+"_"+str(varname)+"_"+"scatter"+"_"+cylcdatestr+".png"
    fig = pyplot.figure()
    map = Basemap(projection='cyl', resolution='c', llcrnrlat= -90.,urcrnrlat= 90.,llcrnrlon=-180.,urcrnrlon=180.)
    #map.bluemarble(scale=0.5);
    map.drawcoastlines()
    map.drawparallels(numpy.arange( -90., 90.,30.),labels=[1,0,0,0],fontsize=5)
    map.drawmeridians(numpy.arange(-180.,180.,30.),labels=[0,0,0,1],fontsize=5,rotation=45)
    #map.drawcountries()
    colors = (0,0,0)
    area = numpy.pi*1.0
    alpha=0.5
    #print(list(data[str(varname)]))
    x, y = map(list(data.Longitude.values), list(data.Latitude.values))
    plot = pyplot.scatter(x,y,s=area,c=colors,alpha=alpha)
    pyplot.title(plot_title)
    #fig = pyplot.figure()
    #cbar=fig.colorbar(plot,orientation='vertical')
    #cbar.ax.tick_params(labelsize=5)
    pyplot.savefig(plotfile,bbox_inches='tight',dpi=200)
    return(fig)

def gridplot(plotpath,odbnmlfile,data,cylcdatestr,obstype,varname,long_name,fieldname="Obsvalue",gridopt="mean",cpallet="jet",fill=False,extend="max"):
    data=obslib.odb_renamefield(data,odbnmlfile)
    clevs=obslib.clevgen(long_name,fieldname)
    units=obslib.dataunit(long_name)
    plot_title=obstype.replace("_"," ")+"\n"+long_name.replace("_"," ")+" ("+units+") "+fieldname+"\n"+" (1x1 gridded "+gridopt+") "+"\n"+cylcdatestr
    plotfile=plotpath+"/"+obstype+"_"+str(varname)+"_"+fieldname+"_"+cylcdatestr+".png"
    if gridopt is "count" : gridded_data=obslib.gridded_count_1x1deg(data,varname)
    elif gridopt is "sum" : gridded_data=obslib.gridded_sum_1x1deg(data,varname,fieldname)
    elif gridopt is "mean" : gridded_data=obslib.gridded_mean_1x1deg(data,varname,fieldname)
    if fill:
        plot_shaded(plotfile,gridded_data,plot_title,clevs,cpallet,extend)
    else:
        plot_button(plotfile,gridded_data,plot_title,clevs,cpallet,extend)
    return(fig)

def plot_density(plotpath,nmlfile,data,cylcdatestr,obstype,element,long_name,fill=False,extend="max"):
    fieldname="ObsDensity"
    gridopt="count"
    cpallet=truncate_colormap("gist_ncar_r", 0.05, 0.4)
    gridplot(plotpath,nmlfile,data,cylcdatestr,obstype,element,long_name,fieldname,gridopt,cpallet,fill,extend)

def plot_gridmean(plotpath,nmlfile,data,cylcdatestr,obstype,element,long_name,fill=False,extend="both"):
    fieldname="Obsvalue"
    gridopt="mean"
    cpallet=truncate_colormap("gist_ncar", 0.1, 0.8)
    gridplot(plotpath,nmlfile,data,cylcdatestr,obstype,element,long_name,fieldname,gridopt,cpallet,fill,extend)

def plot_depart_firstguess(plotpath,nmlfile,data,cylcdatestr,obstype,element,long_name,fill=False,extend="both"):
    fieldname="FGDep"
    gridopt="mean"
    cpallet=truncate_colormap("gist_ncar", 0.1, 0.8)
    gridplot(plotpath,nmlfile,data,cylcdatestr,obstype,element,long_name,fieldname,gridopt,cpallet,fill,extend)

def plot_depart_anal(plotpath,nmlfile,data,cylcdatestr,obstype,element,long_name,fill=False,extend="both"):
    fieldname="AnalDep"
    gridopt="mean"
    cpallet=truncate_colormap("gist_ncar", 0.1, 0.8)
    gridplot(plotpath,nmlfile,data,cylcdatestr,obstype,element,long_name,fieldname,gridopt,cpallet,fill,extend)

def test_ngl():
    wkres = Ngl.Resources()
    rval=subprocess.call("export PYNGL_COLORMAPS=$HOME", shell=True)
    wkres.wkColorMap = "gibies_colourmap_20150117"
    wks = Ngl.open_wks("png","test",wkres)
    Ngl.draw_colormap(wks)
    Ngl.end()


def globalview(plotfile,data,title,clevs=range(0,10,1),cpallet="jet",extend="max",plotmode="shaded"):
    cmap=mplcm.get_cmap(cpallet, len(clevs) - 1)
    glat=data.index
    glon=data.columns
    x,y=numpy.meshgrid(glon,glat)
    fig = pyplot.figure()
    map = Basemap(projection='cyl', resolution='c', llcrnrlat= -90.,urcrnrlat= 90.,llcrnrlon=-180.,urcrnrlon=180.)
    #map.bluemarble(scale=0.5);
    map.drawcoastlines()
    map.drawparallels(numpy.arange( -90., 90.,30.),labels=[1,0,0,0],fontsize=5)
    map.drawmeridians(numpy.arange(-180.,180.,30.),labels=[0,0,0,1],fontsize=5,rotation=45)
    #map.drawcountries()
    if plotmode is "shaded": plot = map.contourf(x,y,data,clevs, cmap=cmap,extend=extend)
    if plotmode is "button": plot = map.scatter(x,y,data,clevs, cmap=cmap,extend=extend)
    pyplot.title(title)
    #ax = fig.add_subplot(311)
    cbar=fig.colorbar(plot,orientation='horizontal')
    cbar.ax.tick_params(labelsize=5)
    pyplot.savefig(plotfile,bbox_inches='tight',dpi=200)


def plot_field(data,plotfile,domain="global",varname="hloswind",textout=False):
	hgtmin=domaindic.dom_hgtmin[domain]
	hgtmax=domaindic.dom_hgtmax[domain]
	plotfile1=plotfile+"_"+domain+"_obs"
	title1="Observation (magnitude) at altitude ["+str(int(hgtmin))+"m to "+str(int(hgtmax))+"m ]"
	gridded_data_obs=obslib.gridded_rms_1x1deg(data,varname=varname,fieldname="observed")
	plot=globalview(plotfile1,gridded_data_obs,title=title1,clevs=[0,1,2,3,4,5,10,15,20,25,35,45],extend="max",plotmode="shaded")
	if textout: obslib.obs_frame_ascii(gridded_data,plotfile1)
	plotfile2=plotfile+"_"+domain+"_bkg"
	title2="Background (magnitude) at altitude ["+str(int(hgtmin))+"m to "+str(int(hgtmax))+"m ]"
	gridded_data_bkg=obslib.gridded_rms_1x1deg(data,varname=varname,fieldname="background")
	plot=globalview(plotfile2,gridded_data_bkg,title=title2,clevs=[0,1,2,3,4,5,10,15,20,25,35,45],extend="max",plotmode="shaded")
	if textout: obslib.obs_frame_ascii(gridded_data,plotfile2)
	plotfile3=plotfile+"_"+domain+"_omb"
	title3="OBS minus BKG (rms) at altitude ["+str(int(hgtmin))+"m to "+str(int(hgtmax))+"m ]"
	gridded_data_omb=obslib.gridded_rms_1x1deg(data,varname=varname,fieldname="obs_minus_bkg")
	plot=globalview(plotfile3,gridded_data_omb,title=title3,clevs=[0,1,2,3,4,5,10],extend="max",plotmode="shaded")
	if textout: obslib.obs_frame_ascii(gridded_data,plotfile3)
	return(gridded_data_omb)
