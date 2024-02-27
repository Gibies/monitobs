#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 12:56:36 2019

@author: gibies
"""

#############################################################################################################################
###  Merged from imdaanowcast package
#############################################################################################################################
from __future__ import print_function
import os,sys
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
PKGHOME=os.environ.get('PKGHOME',os.path.dirname(CURR_PATH))
PKGNAME=os.path.basename(PKGHOME)
LIB=os.environ.get('LIB',PKGHOME+"/pylib")
sys.path.append(LIB)
DIC=os.environ.get('DIC',PKGHOME+"/pydic")
sys.path.append(DIC)
NML=os.environ.get('NML',PKGHOME+"/nml")
sys.path.append(NML)
PALETTE=os.environ.get('PYNGL_COLORMAPS',PKGHOME+"/palette")
SUBTYPNML=NML+"/obs_subtype.nml"

diaglev=int(os.environ.get('GEN_MODE',0))
def errprint(*args, **kwargs):
    if diaglev > 0: print(*args, file=sys.stderr, **kwargs)

import datadic
import subprocess

import obslib
import pplib
import vardic
import glob,datetime
import Nio, Ngl
#import netCDF4
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
pyplot.switch_backend('agg')
import matplotlib.colors as colors
import matplotlib.cm as mplcm
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas
import numpy
import domaindic
import math

#import geocat.datafiles as gdf
#from geocat.viz import cmaps as gvcmaps
#from geocat.viz import util as gvutil

import cartopy
import cartopy.crs as ccrs
import xarray
import iris

cmapfile=os.environ.get('CMAP',PALETTE+"/gibies_colourmap_20150117.rgb")

diaglev=int(os.environ.get('GEN_MODE',0))
def errprint(*args, **kwargs):
    if diaglev > 0: print(*args, file=sys.stderr, **kwargs)

def load_cubes(infile):
    file_cubes=iris.load(infile)
    return(file_cubes)

def umff_to_grib2(infile,outfile=None):
    if outfile is None: outfile=infile.split(".",1)[0]+".grib2"
    file_cubes=iris.load(infile)
    iris.save(file_cubes,outfile)

def read_grib(gribfile):
	filtered_messages = []
	for msg in iris_grib.message.GribMessage.messages_from_filename(gribfile):
	    if msg.sections[4]['productDefinitionTemplateNumber'] is not None:
		print(msg.sections[4]['productDefinitionTemplateNumber'])
		print(msg.sections[4]['parameterNumber'])
		print(msg.sections[4]['typeOfFirstFixedSurface'])
		print(msg.sections[1])
		print(msg.sections[3])
		print(msg.sections[4])
		print(msg.sections[5])
		print(msg.sections[6])
		print(msg.sections[7])
		print(msg.sections[8])
		filtered_messages.append(msg)
	cubes_messages = iris_grib.load_pairs_from_fields(filtered_messages)
	print(cubes_messages)
	for cube, msg in cubes_messages:
	   prod_stat = msg.sections[1]['productionStatusOfProcessedData']
	   cube.attributes['productionStatusOfProcessedData'] = prod_stat
	print(cube.attributes['productionStatusOfProcessedData'])
		#msg.sections[4]['productDefinitionTemplateNumber']; 
		#msg.sections[4]['parameterNumber']
	#cubes = iris.load(gribfile)
	#return(cubes)

def cube_list(infile):
	cubes = read_grib(infile)
	cubes = list(cubes)
	um2grb2.tweaked_messages(cubes)
	return(cubes)

def save_grib(gribfile,cubes):
	iris.save(cubes, gribfile)

def get_keyrange(dataset,keyfield):
    lblset=set()
    for data in dataset["data"]:
       for keyfieldvalue in data[keyfield]:
          lblset.add(keyfieldvalue)
    lblist=list(lblset)
    return(lblist)

def data_check(data):
        #x=numpy.array(data.Longitude.values)
        #y=numpy.array(data.Latitude.values)
	#xcheck=(x*10%10)
	#ycheck=(y*10%10)
	#data=data[obslib.compare(xcheck,ycheck)]
	data=data[data.Longitude.values != data.Latitude.values]
	print(data)
	return(data)	

def reshape_frame(data,data_field,series_field="date",height_field="height"):
    data=data[[series_field,height_field,data_field]]
    indx=0
    timelist=list(data[series_field].unique())
    hoff_data=pandas.DataFrame(columns=list([series_field])+list(data.height.unique()))
    for i,time in enumerate(timelist):
	time_data=obslib.frame_select_filter(data,series_field,time)[data_field].values
	hoff_data=obslib.append(hoff_data,[time]+list(time_data))
    return(hoff_data)

def extract_series_height(data,series,height,data_field,series_field="date",height_field="height_cog",thickness=1000.0,intervel=1.0):
	height_relax=thickness/2
	series_relax=intervel/2
        if series_field in ["Latitude"]:
		data=obslib.frame_window_filter(data,series_field,(series-series_relax),(series+series_relax))
	else:
		data=obslib.frame_select_filter(data,series_field,series)
        data=obslib.frame_window_filter(data,height_field,(height-height_relax),(height+height_relax))
        obscount=len(data)
        stat_table=pandas.DataFrame(columns=[series_field,"height",data_field,obscount])
        if obscount > 0:
	   datamean=obslib.mean_of(data[data_field].values)
	   #print(obscount,datamean,series,height)
        else:
 	   datamean=numpy.nan
        obslib.append(stat_table,[series,height,datamean,obscount])
        return(stat_table)

def qcheck(data,checkfield="WindErr",cutofflim=10,flagfield="AladinConfFlag"):
	cutoffmin=-1*cutofflim
	cutoffmax=cutofflim
	if flagfield in data.columns: data=obslib.frame_select_filter(data,flagfield,0)
	data=obslib.frame_window_filter(data,checkfield,cutoffmin,cutoffmax)
	return(data)

def get_layermean(full_data,data_field,series_field="date",height_field="height_cog",height_min=1000,height_max=20000,thickness=1000):
    if series_field in ["Latitude"]:
	serieslist=numpy.arange(-90.0,91.0,1.0)
    else:
	serieslist=list(full_data[series_field].unique())
    tindx=0
    for i,series in enumerate(serieslist):
      hindx=0
      for height in range(int(height_min),int(height_max),int(thickness)):
	new_data=extract_series_height(full_data,series,height,data_field=data_field,series_field=series_field,height_field=height_field,thickness=thickness)
        if hindx == 0:
           stats=new_data
        else:
           stats=stats.append(new_data)
        hindx+=1
      series_data=reshape_frame(stats,data_field,series_field,"height")
      if tindx == 0:
	data=series_data
      else:
	data=data.append(series_data)
      tindx+=1
    return(data)

def filter(dataset,filter_field="AzimuthCOG",filter_min=0,filter_max=180):
    datalist=[]
    for data in dataset["data"]:
    	data=obslib.frame_window_filter(data,filter_field,filter_min,filter_max)
	datalist.append(data)
    dataset.update({"data":datalist})
    return(dataset)

def select(dataset,select_field="ChannelNumber",select_value=0):
    datalist=[]
    for data in dataset["data"]:
    	data=obslib.frame_select_filter(data,select_field,select_value)
	datalist.append(data)
    dataset.update({"data":datalist})
    return(dataset)

def transpo(data,index_name,column_name):
    data=data.rename(columns={str(index_name):str(column_name)})
    data=data.set_index(column_name)
    data=data.transpose()
    return(data)

def get_hofmuller(dataset,series_field,height_field,data_field,batchidx=0,height_min=1000,height_max=20000,thickness=1000,qcheck_enable=True):
    data=dataset["data"][batchidx]
    #field_val_range=daview.get_keyrange({"data":[full_data]},field)
    if qcheck_enable: data=qcheck(data)
    data=get_layermean(data,data_field=data_field,series_field=series_field,height_field=height_field,height_min=height_min,height_max=height_max,thickness=thickness)
    data=transpo(data,str(series_field),str("height"))
    print(data)
    outfile=dataset["plotfile"]+".txt"
    obslib.obs_frame_ascii(data,outfile)
    print(outfile)
    dataset.update({"data":data})
    if "cnlev" not in dataset: dataset.update({"cnlev":get_cnlev(data)})
    if "clrindx" not in dataset: dataset.update({"clrindx":range(172,225,1)})
    plot_raster_fill(dataset)
    return(data)


def get_keyrange(dataset,keyfield):
    lblset=set()
    for data in dataset["data"]:
	for keyfieldvalue in data[keyfield]:
	   lblset.add(keyfieldvalue)
    lblist=list(lblset)
    return(lblist)

def mpl_plot_keyfield(dataset,plotfile,tagmark="",lblst=[],text="",textpos=(0.25, -0.20)):
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

def get_cnlev(data):
	minval=numpy.nanmin(data.values)
	maxval=numpy.nanmax(data.values)
	print(minval,maxval)
	cnlev=numpy.arange(minval,maxval,10)
	return(cnlev)

def satpassloc(df):
	tmXBValues,indx = numpy.unique(df["Time"],return_index=True)
	tmXBLabels = ["("+"{:7.2f}".format(df["Latitude"][idx])+","+"{:7.2f}".format(df["Longitude"][idx])+")" for idx in indx]
	return(tmXBValues,tmXBLabels)


def ngl_plot_raster_fill(dataset={},data_hoff=None,cnlev=None,clrindx=None,title="",lstr="",rstr="",foottext="",plot_btmlev=0.0,plot_toplev=20000.0,plotfile="test",wks_type="png",cmapfile=cmapfile,fillval=-999.99):
	if "cnlev" not in dataset: dataset.update({"cnlev":[-25,-10,-5,-2,-1,-0.5,0.0,0.5,1,2,5,10,25]})
	if "cnlev" in dataset: cnlev=dataset["cnlev"]
	if "clrindx" not in dataset: dataset.update({"clrindx":range(172,204,1)+range(218,255,1)})
	if "clrindx" in dataset: clrindx=dataset["clrindx"]
	if "data" in dataset: data_hoff=dataset["data"]
	if "title" in dataset: title=dataset["title"]
	if "lstr" in dataset: lstr=dataset["lstr"]
	if "rstr" in dataset: rstr=dataset["rstr"]
	if "foottext" in dataset: foottext=dataset["foottext"]
	if "plot_btmlev" in dataset: plot_btmlev=dataset["plot_btmlev"]
	if "plot_toplev" in dataset: plot_toplev=dataset["plot_toplev"]
	if "plotfile" in dataset: plotfile=dataset["plotfile"]
	if "wks_type" in dataset: wks_type=dataset["wks_type"]
	if "cmapfile" in dataset: cmapfile=dataset["cmapfile"]
	if "fillval" in dataset: fillval=dataset["fillval"]
	print(plotfile+"."+wks_type)
	#if "height" not in data_hoff.columns:
 	#   col0=data_hoff.columns[0]
	#   data_hoff.rename(columns = {col0:"height"}, inplace = True)
	#   data_hoff=data_hoff.set_index("height")
	print(data_hoff)
	#data_hoff=data_hoff.truncate(plot_btmlev,plot_toplev)
	data_hoff.fillna(fillval,inplace=True)
	height=list(data_hoff.index)
	series=list(data_hoff.columns)
	tindx=list(range(0,len(series),1))
	hindx=list(range(0,len(height),1))

	if "tmXBValues" in dataset: 
		tmxbvals=dataset["tmXBValues"]
	else:
		tmxbvals=range(0,len(series))

	if "tmXBLabels" in dataset:
		tmxblbls=dataset["tmXBLabels"]
	else:
		tmxblbls=[str(int(series[lb]))+"" for lb in range(0,len(series))]

	if "tmYLValues" in dataset:
		tmylvals=dataset["tmYLValues"]
	else:
		tmylvals=height[::2]
		#tmylvals=[2000,4000,6000,8000,10000,12000,14000,16000,18000,20000]
	
	if "tmYLLabels" in dataset:
		tmyllbls=dataset["tmYLLabels"]
	else:
		tmyllbls=[str(int(hgt/1000))+"km" for hgt in tmylvals]
		#tmyllbls=["2km","4km","6km","8km","10km","12km","14km","16km","18km","20km"]


	print(tmylvals,tmyllbls)
	
	print(min(data_hoff),max(data_hoff))

	print(height)
	print(tmylvals)
	print(series)
	print(tmxbvals)

	cmap = Ngl.read_colormap_file(cmapfile)[clrindx,:]
	wks = Ngl.open_wks(wks_type,plotfile)
	resources = Ngl.Resources()

	resources.tiMainString   = "~F25~"+title  # Main title.
	resources.tiLeftString   = "~F25~"+title  # Main title.
	resources.tiRightString   = "~F25~"+title  # Main title.

	if "sfXArray" in dataset: resources.sfXArray = dataset["sfXArray"]
	resources.sfXCStartV = tmxbvals[0]  # Indicate start and end of left
	resources.sfXCEndV   = tmxbvals[-1]   # Y axes values.
	resources.tmXBMode      = "Explicit"   # Define your own tick mark labels.
	resources.tmXBLabelFont = "times-roman"  # Change font of labels.
	resources.tmXBLabelFontHeightF = 0.015 # Change font height of labels.
	resources.tmXBMinorOn   = False        # No minor tick marks.
	resources.tmXBValues    = tmxbvals # Location to put tick mark labels
	resources.tmXBLabels    = tmxblbls
	resources.tmXBLabelAngleF = 45
	resources.tmXBLabelJust = "TopRight"
	if "trXAxisType" in dataset: resources.trXAxisType = dataset["trXAxisType"]
	if "trXCoordPoints" in dataset: resources.trXCoordPoints = dataset["trXCoordPoints"]
	if "trXReverse" in dataset:
		resources.trXReverse  = dataset["trXReverse"]
	else:
		resources.trXReverse  = False    # Reverse the X values.
	print(resources.trXReverse)

	#resources.tiYAxisString  = "~F25~Height"  # Y axes label.
	resources.sfYCStartV = tmylvals[0]  # Indicate start and end of left
	resources.sfYCEndV   = tmylvals[-1]   # Y axes values.
	resources.tmYLMode      = "Explicit" # Define own tick mark labels.
	resources.tmYLLabelFont = "times-roman"  # Change the font.
	resources.tmYLMinorOn   = True        # No minor tick marks.
	resources.tmYLValues    = tmylvals
	resources.tmYLLabels    = tmyllbls
	if "trYReverse" in dataset:
                resources.trYReverse  = dataset["trYReverse"]
	else:
		resources.trYReverse  = False    # Reverse the Y values.
	resources.trYLog      = False    # Use log scale.

	resources.lbLabelBarOn 	= True
	resources.pmLabelBarDisplayMode = "Always"    # Turn off label bar.


	resources.cnFillOn      = True  # Turn on contour level fill.
	resources.cnFillMode	= "RasterFill"
	resources.cnLinesOn	= False
	resources.cnLineLabelsOn	= False
	resources.cnInfoLabelOn	= False
	resources.cnLevelSelectionMode = "ExplicitLevels"
	resources.cnLevels	= cnlev #[0,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500,550,600,650,700,750]
	resources.cnFillPalette = cmap
	resources.cnLineLabelAngleF = 0. # Draw contour line labels right-side up.
    	resources.sfMissingValueV = fillval

	resources.nglDraw  = False  # Don't draw the plot or advance the
	resources.nglFrame = False  # frame in the call to Ngl.contour.

	resources.nglMaximize = False
	print(data_hoff)
	contour = Ngl.contour(wks, data_hoff.values, resources)  # Create a contour plot.

	Ngl.draw(contour)  # Draw the contour plot.

	txres               = Ngl.Resources()    # Annotate plot with some text.
	txres.txFontHeightF = 0.015
	txres.txAngleF      = 0.
	Ngl.text_ndc(wks,"~F25~"+foottext,.5,.05,txres)

	Ngl.frame(wks) # Advance the frame.


def ngl_plot_hoff(data_hoff,cnlev,clrindx,title="",lstr="",rstr="",foottext="",plot_btmlev=0.0,plot_toplev=20000.0,plotfile="test",wks_type="png",cmapfile=cmapfile,fillval=-999.99):
	data_hoff=data_hoff.truncate(plot_btmlev,plot_toplev)
	data_hoff.fillna(fillval,inplace=True)
	height=list(data_hoff.index)
	time=list(data_hoff.columns)
	tindx=list(range(0,len(time),1))
	hindx=list(range(0,len(height),1))

	cmap = Ngl.read_colormap_file(cmapfile)[clrindx,:]
	wks = Ngl.open_wks(wks_type,plotfile)
	resources = Ngl.Resources()

	resources.tiMainString   = "~F25~"+title  # Main title.
	resources.tiLeftString   = "~F25~"+title  # Main title.
	resources.tiRightString   = "~F25~"+title  # Main title.

	resources.tmXBMode      = "Explicit"   # Define your own tick mark labels.
	resources.tmXBLabelFont = "times-roman"  # Change font of labels.
	resources.tmXBLabelFontHeightF = 0.015 # Change font height of labels.
	resources.tmXBMinorOn   = False        # No minor tick marks.
	resources.tmXBValues    = tindx[3::40] # Location to put tick mark labels
	resources.tmXBLabels    = time[3::40]
	resources.tmXBLabelAngleF = 20
	resources.tmXBLabelJust = "TopRight"
	resources.trXReverse  = True    # Reverse the X values.

	#resources.tiYAxisString  = "~F25~Height"  # Y axes label.
	resources.sfYCStartV = height[0]  # Indicate start and end of left
	resources.sfYCEndV   = height[-1]   # Y axes values.
	resources.tmYLMode      = "Explicit" # Define own tick mark labels.
	resources.tmYLLabelFont = "times-roman"  # Change the font.
	resources.tmYLMinorOn   = True        # No minor tick marks.
	resources.tmYLValues    = [2000,4000,6000,8000,10000,12000,14000,16000,18000,20000]
	resources.tmYLLabels    = ["2km","4km","6km","8km","10km","12km","14km","16km","18km","20km"]
	resources.trYReverse  = False    # Reverse the Y values.
	resources.trYLog      = False    # Use log scale.

	resources.lbLabelBarOn 	= True
	resources.pmLabelBarDisplayMode = "Always"    # Turn off label bar.

	resources.cnFillOn          = True  # Turn on contour level fill.
	resources.cnLinesOn	= False
	resources.cnLineLabelsOn	= False
	resources.cnInfoLabelOn	= False
	resources.cnLevelSelectionMode = "ExplicitLevels"
	resources.cnLevels	= cnlev #[0,10,20,30,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500,550,600,650,700,750]
	resources.cnFillPalette = cmap
	resources.cnLineLabelAngleF = 0. # Draw contour line labels right-side up.
    	resources.sfMissingValueV = fillval

	resources.nglDraw  = False  # Don't draw the plot or advance the
	resources.nglFrame = False  # frame in the call to Ngl.contour.

	resources.nglMaximize = False
	contour = Ngl.contour(wks, data_hoff.values, resources)  # Create a contour plot.

	Ngl.draw(contour)  # Draw the contour plot.

	txres               = Ngl.Resources()    # Annotate plot with some text.
	txres.txFontHeightF = 0.015
	txres.txAngleF      = 0.
	Ngl.text_ndc(wks,"~F25~"+foottext,.5,.05,txres)

	Ngl.frame(wks) # Advance the frame.

def get_fillval(dataset,fillvalfix=-999.99):
	try:fillval = dataset["_FillValue"]
	except:fillval = fillvalfix
	return(fillval)

def set_fillval(dataset):
	fillval=get_fillval(dataset)
    	dataset.update({"_FillValue":fillval})
	return(dataset)

def get_mask(dataset):
     if "data" in dataset:
	data=dataset["data"]
     else:
	if "zonal_data" in dataset:
	    data=dataset["zonal_data"]
	else:
	    if "merid_data" in dataset:
		data=dataset["merid_data"]
     mask=numpy.where(numpy.isnan(data))
     dataset.update({"mask":mask})
     return(dataset)

def mask_data(dataset,itemlist):
     for item in itemlist:
	if item in dataset:
	   data=dataset[item]
	   mask=dataset["mask"]
	   fillval=dataset["_FillValue"]	
	   data[mask]=fillval
	   dataset.update({item:data})
     return(dataset)

def fillval_mask(dataset):
	dataset=set_fillval(dataset)
	dataset=get_mask(dataset)
	dataset=mask_data(dataset,["data","zonal_data","merid_data"])
    	#data.fillna(fillval,inplace=True)
	#data=numpy.nan_to_num(data,copy=True,nan=fillval,posinf=None,neginf=None)
	return(dataset)

def get_2d_lat(dataset):
    shape=dataset["data"].shape
    datadir=dataset["datadir"]
    fname = dataset["filename"]
    datafile=datadir + "/" + fname
    infile = Nio.open_file(datafile,"r")
    latnam = dataset["latvar"]
    opt = dataset["dimopt"]
    if opt == 1:
	lat = infile.variables[latnam][:]
	lat2d=pplib.get_2dlat(lat,shape)
    if opt == 2:
    	lat2d = infile.variables[latnam][:,:]
    dataset.update({"lat2d":lat2d})
    return(dataset)

def get_2d_lon(dataset):
    shape=dataset["data"].shape
    datadir=dataset["datadir"]
    fname = dataset["filename"]
    datafile=datadir + "/" + fname
    infile = Nio.open_file(datafile,"r")
    lonnam = dataset["lonvar"]
    opt = dataset["dimopt"]
    if opt == 1:
	lon = infile.variables[lonnam][:]
	lon2d=pplib.get_2dlon(lon,shape)
    if opt == 2:
    	lon2d = infile.variables[lonnam][:,:]
    dataset.update({"lon2d":lon2d})
    return(dataset)

def get_data(dataset):
    datadir=dataset["datadir"]
    fname = dataset["filename"]
    vname = dataset["fvarname"]
    if "time_index" in dataset: 
	timindx = dataset["time_index"]
    else:
	timindx = None
    if "lev_index" in dataset: 
	levindx = dataset["lev_index"]
    else:
	levindx = None
    datafile=datadir + "/" + fname
    infile = Nio.open_file(datafile,"r")
    if timindx is not None and levindx is not None: dataset.update({"data":infile.variables[vname][timindx,levindx,:,:]})
    if timindx is None and levindx is not None: dataset.update({"data":infile.variables[vname][levindx,:,:]})
    if timindx is not None and levindx is None: dataset.update({"data":infile.variables[vname][timindx,:,:]})
    if timindx is None and levindx is None: dataset.update({"data":infile.variables[vname][:,:]})
    dataset.update(infile.variables[vname].attributes)
    dataset=fillval_mask(dataset)
    dataset = get_2d_lat(dataset)
    dataset = get_2d_lon(dataset)
    return(dataset)

def get_vector_magnitude(dataset):
    uwnd=dataset["zonal_data"]
    vwnd=dataset["merid_data"]
    usqr=uwnd.__pow__(2)
    vsqr=vwnd.__pow__(2)
    ssqr=usqr.__add__(vsqr)
    speed=numpy.sqrt(ssqr)
    dataset.update({"data":speed})
    return(dataset)

def get_vector_data(dataset):
    datadir=dataset["datadir"]
    fname = dataset["filename"]
    zonal_vname = dataset["fvarname"][0]
    merid_vname = dataset["fvarname"][1]
    if "time_index" in dataset: 
	timindx = dataset["time_index"]
    else:
	timindx = 0
    if "lev_index" in dataset: 
	levindx = dataset["lev_index"]
    else:
	levindx = None
    datafile=datadir + "/" + fname
    infile = Nio.open_file(datafile,"r")
    if timindx is not None and levindx is not None: 
	dataset.update({"zonal_data":infile.variables[zonal_vname][timindx,levindx,:,:]})
	dataset.update({"merid_data":infile.variables[merid_vname][timindx,levindx,:,:]})
    if timindx is None and levindx is not None: 
	dataset.update({"zonal_data":infile.variables[zonal_vname][levindx,:,:]})
	dataset.update({"merid_data":infile.variables[merid_vname][levindx,:,:]})
    if timindx is not None and levindx is None: 
	dataset.update({"zonal_data":infile.variables[zonal_vname][timindx,:,:]})
	dataset.update({"merid_data":infile.variables[merid_vname][timindx,:,:]})
    if timindx is None and levindx is None: 
	dataset.update({"zonal_data":infile.variables[zonal_vname][:,:]})
	dataset.update({"merid_data":infile.variables[merid_vname][:,:]})
    zonal_attributes=infile.variables[zonal_vname].attributes
    merid_attributes=infile.variables[merid_vname].attributes
    dataset.update({"zonal_attrib":zonal_attributes})
    dataset.update({"merid_attrib":merid_attributes})
    dataset=get_vector_magnitude(dataset)
    dataset=fillval_mask(dataset)
    if "grid_type" in dataset:
       grid_type = dataset["grid_type"]
    else:
       grid_type = None
    if grid_type is None:
       dataset = get_2d_lat(dataset)
       dataset = get_2d_lon(dataset)
    
    return(dataset)

def get_colrindx(colrdic):
    if "colrindx" not in colrdic or colrdic["colrindx"] is None:
    	try:colrmin=colrdic["colrmin"]
	except:colrmin=0
    	try:colrmax=colrdic["colrmax"]
	except:colrmax=170
	colrindx=numpy.arange(colrmin,colrmax,1)
    else:
    	colrindx=colrdic["colrindx"]
    return(colrindx)

def ngl_print_rgb(colrdic):
    colrindx=get_colrindx(colrdic)
    cmap = Ngl.read_colormap_file(cmapfile)[colrindx,:]
    print(cmap)

def ngl_orth_plot(wks,data,res):
	res.vpWidthF               =  0.8               #-- width of plot
	res.vpHeightF              =  0.8               #-- height of plot
	res.mpFillOn               =  True              #-- map fill on
	res.mpOceanFillColor       = "Transparent"      #-- default: dark blue
	res.mpLandFillColor        = "Gray90"           #-- default: dark red
	res.mpInlandWaterFillColor = "Transparent"      #-- default: white
	res.mpProjection           = "Orthographic"     #-- projection type
	res.tiMainString           = "Orthographic projection" #-- title
	#-- create the plot
	map = Ngl.map(wks,res)
	return(map)

def ngl_cyl_plot(wks,dataset,res):
    res.mpDataBaseVersion      = "MediumRes"
    res.mpFillOn              = True
    res.mpFillColors = [0,-1,-1,-1]
    res.mpFillAreaSpecifiers  = ["land"]
    res.mpSpecifiedFillColors = ["gray65"]
    res.mpGridLatSpacingF =  10.                  #-- grid lat spacing
    res.mpGridLonSpacingF =  10.                  #-- grid lon spacing
    res.mpLimitMode       = "LatLon"              #-- must be set using minLatF/maxLatF/minLonF/maxLonF
    res.mpMinLatF    = -90
    res.mpMaxLatF    = 90
    res.mpMinLonF    = 0
    res.mpMaxLonF    = 360
    res.nglDraw               = False               #-- don't draw individual plots
    res.nglFrame              = False               #-- don't advance frame

    if dataset["plot_type"] in ["contour","vector_scalar",]:
	print(res.sfXArray.shape)
	print(res.sfYArray.shape)
	print(res.cnLevels)
    						#plot1=Ngl.gsn_csm_contour_map_ce(wks,data,res) is not available
    if dataset["plot_type"] is "vector_scalar":
	uwnd=dataset["zonal_data"]
	vwnd=dataset["merid_data"]
    	data=dataset["data"]
	print(uwnd.shape)
	print(vwnd.shape)
	print(data.shape)
	res.vcMonoLineArrowColor = True  
    	plot1=Ngl.vector_scalar_map(wks, uwnd, vwnd, data, res)
    if dataset["plot_type"] in ["vector",]:
	uwnd=dataset["zonal_data"]
	vwnd=dataset["merid_data"]
	print(uwnd.shape)
	print(vwnd.shape)
	res.vcMonoLineArrowColor = False  # Draw vectors in color.
    	plot1=Ngl.vector_map(wks, uwnd, vwnd, res)
    if dataset["plot_type"] in ["contour",]:
	data=dataset["data"]
	print(data.shape)
	plot1=Ngl.contour_map(wks,data,res)
    return(plot1)
	

def ngl_get_subplot(dataset,wks):
    data=dataset["data"]
    lat2d=dataset["lat2d"]
    lon2d=dataset["lon2d"]
    if "cnlev" in dataset:
	cnlev= dataset["cnlev"]
    else:
	cnlev= list(numpy.arange(0.,32.,2.))
	print(cnlev)
    if "title" in dataset: 
	title=dataset["title"]
    else:
	title=""
    colrindx=get_colrindx(dataset)
    cmap = Ngl.read_colormap_file(cmapfile)[colrindx,:]
    res = Ngl.Resources()
    if "draworder" in dataset:
    	res.mpFillDrawOrder   = dataset["draworder"]
    
    dataset["plot_type"]="contour"

    if dataset["plot_type"] in ["contour","vector_scalar",]:
    	res.cnFillOn = True
   	res.cnLineOn = False
    	res.cnLineLabelsOn        = False
    	res.cnInfoLabelOn         = False
    	res.cnLevelSelectionMode = "ExplicitLevels"
    	res.cnLevels             = cnlev
    	res.cnFillPalette               = cmap

    if dataset["plot_type"] in ["vector","vector_scalar",]:
    	res.vcLevelPalette       = cmap
	res.vcMinFracLengthF     = 0.1   # Increase length of
	res.vcMinMagnitudeF      = 0.001  # vectors.
	res.vcRefLengthF         = 0.045
	res.vcRefMagnitudeF      = 2.0

    res.lbOrientation  = "Horizontal"
    res.sfXArray     = lon2d
    res.sfYArray     = lat2d
    res.sfMissingValueV = get_fillval(dataset)
    res.gsnAddCyclic = False  
    res.tiMainString = title
    if "gsnLeftString" in dataset:
    	res.gsnLeftString = dataset["gsnLeftString"]
    else:
	res.gsnLeftString = ""
    plot1=cyl_plot(wks,dataset,res)
    return(plot1)


def ngl_plot_data(dicset_list,plotfile,wks_type):
    wkres = Ngl.Resources()
    wks = Ngl.open_wks(wks_type,plotfile,wkres)
    plot = []
    for count,dataset in enumerate(dicset_list):
	if dataset["plot_type"] in ["contour"]:
		dataset=get_data(dataset)
	else:
		dataset=get_vector_data(dataset)
	plot1=get_subplot(dataset,wks)
	plot.append(plot1)
    count=len(plot)
    panelres                  =  Ngl.Resources()
    panelres.nglPanelLabelBar =  False           #-- common labelbar
    panelres.nglPanelYWhiteSpacePercent =  0        #-- reduce space between the panel plots
    panelres.nglPanelXWhiteSpacePercent =  0        #-- reduce space between the panel plots
    panelres.nglPanelTop      =  0.95               #-- top position of panel
    #-- create the panel
    Ngl.panel(wks,plot,[count,1],panelres)
    Ngl.end()
    return(plot)

def get_dataset_bias(dataset,dataref):
	dataset_bias=dataset.copy()
	data=dataset["data"]
	stnd=dataref["data"]
	bias=data-stnd
	dataset_bias.update({"title":dataset["title"]+" deviation"})
    	dataset_bias.update({"colrmin":172})
    	dataset_bias.update({"colrmax":255})
	dataset_bias.update({"colrindx":vardic.get_dicval(dataset["plotname"]+"_bias","colrindx")})		#numpy.array([172,179,186,193,200,205,226,233,240,247,254])
    	dataset_bias.update({"cnlev":vardic.get_dicval(dataset["plotname"]+"_bias","cnlev")})		#numpy.linspace(-1,1,9)
	dataset_bias.update({"data":bias})
	return(dataset_bias)


def ngl_plot_bias(dicset_list,plotfile,wks_type):
    wkres = Ngl.Resources()
    wks = Ngl.open_wks(wks_type,plotfile,wkres)
    plot = []
    for count,dicset1 in enumerate(dicset_list):
	dataset=get_data(dicset1)
	if count is 0:
		dataref=dataset
	else:
		dataset=get_dataset_bias(dataset,dataref)
	print(count)
	print(dataset)
	plot1=get_subplot(dataset,wks)
	plot.append(plot1)
    count=len(plot)
    panelres                  =  Ngl.Resources()
    panelres.nglPanelLabelBar =  False           #-- common labelbar
    panelres.nglPanelYWhiteSpacePercent =  0        #-- reduce space between the panel plots
    panelres.nglPanelXWhiteSpacePercent =  0        #-- reduce space between the panel plots
    panelres.nglPanelTop      =  0.95               #-- top position of panel
    #-- create the panel
    Ngl.panel(wks,plot,[count,1],panelres)
    Ngl.end()
    return(plot)

def correct_cnlev(dataset,lev):
    if dataset["fvarname"] in ["air_temperature",]:
      if int(lev) >= 850:
	dataset.update({"cnlev" : numpy.arange(260.,305.,5.)})
      if int(lev) < 850 and int(lev) >= 500:
	dataset.update({"cnlev" : numpy.arange(230.,275.,5.)})
      if int(lev) < 500 and int(lev) >= 200:
	dataset.update({"cnlev" : numpy.arange(200.,235.,5.)})
    print(lev,dataset["cnlev"])
    return(dataset)

def convert_theta_to_airt(dataset,lev):
    if dataset["fvarname"] in ["air_potential_temperature",]:
	mulfactor=math.pow((float(lev)*0.001),0.286)
	data=dataset["data"]
	data=data*mulfactor
	mask=dataset["mask"]
	fillval=dataset["_FillValue"]
	data[mask]=fillval
	dataset.update({"data":data})
	dataset.update({"fvarname":"air_temperature"})
	dataset=correct_cnlev(dataset,lev)
    else:
	print("File variable name is "+str(dataset["fvarname"]))
    return(dataset)

def ngl_plot_stdlev(dicset_list,plotfile,wks_type,lev_list=vardic.stdlev):
    dicset=dicset_list[0].copy()
    wkres = Ngl.Resources()
    wks = Ngl.open_wks(wks_type,plotfile,wkres)
    plot = []
    for count,lev in enumerate(lev_list):
	dataset=dicset.copy()
	dataset.update({"title":""})
	dataset.update({"gsnLeftString":"Level "+str(lev)+" hPa"})
	fnampreix=dataset["filename"].split('.')[0]
	dataset.update({"filename":fnampreix+".nc"})
	lev_index=obslib.str_list_index(lev,vardic.lev_list)
	print(lev,lev_index, vardic.lev_list)
	dataset.update({"lev_index":lev_index})
	if dataset["plot_type"] in ["contour"]:
		dataset=get_data(dataset)
	else:
		dataset=get_vector_data(dataset)
	dataset=convert_theta_to_airt(dataset,lev)
	plot1=get_subplot(dataset,wks)
	plot.append(plot1)
    count=len(plot)
    panelres                  =  Ngl.Resources()
    panelres.nglPanelLabelBar =  False           #-- common labelbar
    panelres.nglPanelYWhiteSpacePercent =  0        #-- reduce space between the panel plots
    panelres.nglPanelXWhiteSpacePercent =  0        #-- reduce space between the panel plots
    panelres.nglPanelTop      =  0.95               #-- top position of panel
    #-- create the panel
    Ngl.panel(wks,plot,[count,1],panelres)
    Ngl.end()
    return(plot)




def ngl_plot_databias(dicset_list,plotfile,wks_type):
    wkres = Ngl.Resources()
    #wkres.wkColorMap = "default"
    wks = Ngl.open_wks(wks_type,plotfile,wkres)
    plot = []
    if len(dicset_list) == 2:
	dicset1=dicset_list[0]
	dicset2=dicset_list[1]
	dataset1=get_data(dicset1)
	plot1=get_subplot(dataset1,wks)
	plot.append(plot1)
	dataset2=get_data(dicset2)
	plot2=get_subplot(dataset2,wks)
	plot.append(plot2)
	dataset3=get_dataset_bias(dataset2,dataset1)
	plot3=get_subplot(dataset3,wks)
	plot.append(plot3)
    count=len(plot)
    panelres                  =  Ngl.Resources()
    panelres.nglPanelLabelBar =  False           #-- common labelbar
    panelres.nglPanelYWhiteSpacePercent =  0        #-- reduce space between the panel plots
    panelres.nglPanelXWhiteSpacePercent =  0        #-- reduce space between the panel plots
    panelres.nglPanelTop      =  0.95               #-- top position of panel
    #-- create the panel
    Ngl.panel(wks,plot,[count,1],panelres)
    Ngl.end()
    return(plot)

def ngl_plot_grid_data(data,plotfile,wks_type):
	wkres = Ngl.Resources()
	wkres.wkColorMap = "default"
	wks = Ngl.open_wks(wks_type,plotfile,wkres)
	res = Ngl.Resources()
	plot = Ngl.contour_map(wks,data,res)
	return(plot)

def mpl_truncate_colormap(cpallet, minval=0.0, maxval=1.0, n=100):
    cmap=mplcm.get_cmap(cpallet)
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(numpy.linspace(minval, maxval, n)))
    return new_cmap

def mpl_plot_location(data,plotfile,plotmode="cyl"):
    fig = pyplot.figure()
    colors = (0,0,0)
    area = 1.0      #numpy.pi*
    alpha=0.5
    parallels = numpy.arange(-80.,90,20.)
    meridians = numpy.arange(-180.,180.,30.)
    if plotmode is "ortho": 
	plot1=pyplot.subplot(121)
	fig = plot_ortho(data,fig,plot1,colors,area,alpha,parallels,meridians,polelat=10,polelon=70)
	plot2=pyplot.subplot(122)
	fig = plot_ortho(data,fig,plot2,colors,area,alpha,parallels,meridians,polelat=10,polelon=250)
    if plotmode is "merc":
	plot=pyplot.subplot(211)
    	fig= plot_merc(data,fig,plot,colors,area,alpha,parallels,meridians)
    if plotmode is "cyl":
	plot=pyplot.subplot(211)
    	fig= plot_cyl(data,fig,plot,colors,area,alpha,parallels,meridians)
    pyplot.savefig(plotfile,bbox_inches='tight',dpi=200)
    return(fig)

def mpl_plot_merc(data,figure,plot,colors,area,alpha,parallels,meridians):
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

def mpl_plot_cyl(data,figure,plot,colors,area,alpha,parallels,meridians):
    plot = Basemap(projection='cyl', resolution='c', llcrnrlat= -90.,urcrnrlat= 90.,llcrnrlon=-180.,urcrnrlon=180.)
    plot.drawlsmask(land_color='wheat',ocean_color='lightblue',lakes=True)
    #map.bluemarble(scale=0.5);
    plot.drawcoastlines()
    #map.drawcountries()
    plot.drawparallels(parallels,labels=[True,False,True,False],fontsize=5)
    plot.drawmeridians(meridians,labels=[False,False,False,True],fontsize=5,rotation=45)
    x, y = plot(list(data.Longitude.values), list(data.Latitude.values))
    plot = pyplot.scatter(x,y,s=area,c=colors,alpha=alpha)
    return(figure)

def mpl_plot_ortho(data,figure,plot,colors,area,alpha,parallels,meridians,polelat,polelon):
    plot=Basemap(lat_0=polelat, lon_0=polelon, projection='ortho',resolution='l' )
    plot.drawlsmask(land_color='wheat',ocean_color='lightblue',lakes=True)
    plot.drawcoastlines()
    #plot.readshapefile('/home/gibies/maskfiles/AllIndia/AllIndia', 'AllIndia')
    plot.drawparallels(parallels,labels=[True,False,True,False],fontsize=5)
    plot.drawmeridians(meridians,labels=[False,False,False,True],fontsize=5)
    x, y = plot(list(data.Longitude.values), list(data.Latitude.values))
    plot = pyplot.scatter(x,y,s=area,c=colors,alpha=alpha)
    return(figure)

def mpl_globalview(plotfile,data,title,clevs=range(0,10,1),cpallet="jet",extend="max",plotmode="shaded"):
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

def mpl_plot_colourbutton():
   lon=data[:,1]
   lat=data[:,0]
   val=data[:,2]

   cmap=daview.truncate_colormap('jet',minval=0.35, maxval=1.0, n=100)
   cticks=[1,2,5,10,25,50,100,200]
   colrs=cmap              #["lightblue","lightgreen","yellow","orange","red"]
   cnorm=matplotlib.colors.LogNorm(1,250)
   ##########################
   m1 = Basemap(projection='cyl', resolution='l', llcrnrlat=-90, llcrnrlon=0, urcrnrlat=90, urcrnrlon=360)
   m1.drawcoastlines()
   m1.drawcoastlines(color='white', linewidth=0.2)
   m1.drawparallels(numpy.arange(-90.,90.,30.),labels=[1,0,0,0],fontsize=7)
   m1.drawmeridians(numpy.arange(0.,360.,60.),labels=[0,0,0,1],fontsize=7)
   lons,lats=m1(lon,lat)
   sc1=m1.scatter(lons, lats,c=val,cmap=colrs,marker = 'o',s=5,norm=cnorm)
   del val
   #plt.clim(0,150)
   #cbar=plt.colorbar(sc1)
   cbar=plt.colorbar(sc1,ticks=cticks,orientation="vertical",extend="max",shrink=0.9, pad=0.01)
   ticklabs = cbar.ax.get_yticklabels()
   cbar.ax.set_yticklabels(cticks, fontsize=7)
   sc1 = plt.annotate("(a) Reception count (0.5 deg grid)",fontsize=7, xy=(0.01, 1.05), xycoords='axes fraction')
   
   #fig2, ax2 = plt.subplots(212,figsize=(10,8))
   plot2=plt.subplot(212)
   return(sc1)


#############################################################################################################################
###  Merged from imdaanowcast package
#############################################################################################################################

def ngl_resof(plotinfo):
    if "res" in plotinfo:
        res=plotinfo["res"]
    else:
        res = Ngl.Resources()
    return(res)

def ngl_gen_wks(plotinfo):
    plotfile=plotinfo["plotfile"]
    wks_type=plotinfo["wks_type"]
    wkres = Ngl.Resources()
    wkres.wkColorMap = "gibies_colourmap_20150117"
    wks = Ngl.open_wks(wks_type,plotfile,wkres)
    plotinfo.update({"wks":wks})
    return(plotinfo)


def ngl_get_colrindx(colrdic):
    if "colrindx" not in colrdic or colrdic["colrindx"] is None:
    	try:colrmin=colrdic["colrmin"]
	except:colrmin=0
    	try:colrmax=colrdic["colrmax"]
	except:colrmax=170
	colrindx=numpy.arange(colrmin,colrmax,1)
    else:
    	colrindx=colrdic["colrindx"]
    return(colrindx)

def ngl_get_cmap(plotinfo):
    colrindx=get_colrindx(plotinfo)
    cmap = Ngl.read_colormap_file(cmapfile)[colrindx,:]
    if "revcmap" in plotinfo and plotinfo["revcmap"]: cmap=cmap[::-1,:]
    plotinfo.update({"cmap":cmap})
    return(plotinfo)

def print_rgb(colrdic):
    if "cmap" not in colrdic: colrdic=get_cmap(colrdic)
    cmap=colrdic["cmap"]
    print(cmap)

def mpl_truncate_colormap(cpallet, minval=0.0, maxval=1.0, n=100):
    cmap=mplcm.get_cmap(cpallet)
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(numpy.linspace(minval, maxval, n)))
    return new_cmap

def check_plot(plotinfo):
    outfile=os.path.abspath(plotinfo["plotfile"]+"."+plotinfo["wks_type"])
    outfile=subprocess.call("ls "+outfile, shell=True)
    return(outfile)

def ngl_test_ngl(plotinfo={}):
    if "plotfile" not in plotinfo: plotinfo={"plotfile":"palette_check"}
    if "wks_type" not in plotinfo: plotinfo.update({"wks_type":"png"})
    outfile=os.path.abspath(plotinfo["plotfile"]+"."+plotinfo["wks_type"])
    if "wks" not in plotinfo: plotinfo=gen_wks(plotinfo)
    os.environ["NCARG_COLORMAPS"]=PALETTE
    os.environ["PYNGL_COLORMAPS"]=PALETTE
    cmapfile=subprocess.call("ls "+PALETTE+"/gibies_colourmap_20150117.rgb", shell=True)
    
    plot=Ngl.draw_colormap(plotinfo["wks"])
    Ngl.end()
    outfile=check_plot(plotinfo)
    print(cmapfile,outfile)
    return(outfile)

def ngl_write_tlstring(plotinfo):
    vpx = Ngl.get_float(plotinfo["plot"],"vpXF")
    vpy = Ngl.get_float(plotinfo["plot"],"vpYF")
    vpw = Ngl.get_float(plotinfo["plot"],"vpWidthF")
    vph = Ngl.get_float(plotinfo["plot"],"vpHeightF")
    if not "tlsres" in plotinfo:
	plotinfo.update({"tlsres":Ngl.Resources()})
	plotinfo["tlsres"].txFontHeightF = 0.014
	plotinfo["tlsres"].txJust        = "CenterLeft"
    if not "tlstring" in plotinfo: 
 	print("No string for Top Left corner")
        plotinfo.update({"tlstring":"____"})
    Ngl.text_ndc(plotinfo["wks"], plotinfo["tlstring"], vpx+0.05*vpw, vpy, plotinfo["tlsres"])
    return(plotinfo)

def ngl_write_trstring(plotinfo):
    vpx = Ngl.get_float(plotinfo["plot"],"vpXF")
    vpy = Ngl.get_float(plotinfo["plot"],"vpYF")
    vpw = Ngl.get_float(plotinfo["plot"],"vpWidthF")
    vph = Ngl.get_float(plotinfo["plot"],"vpHeightF")
    if not "trsres" in plotinfo:
	plotinfo.update({"trsres":Ngl.Resources()})
	plotinfo["trsres"].txFontHeightF = 0.014
	plotinfo["trsres"].txJust        = "CenterRight"
    if not "trstring" in plotinfo: 
 	print("No string for Top Right corner")
        plotinfo.update({"trstring":"____"})
    Ngl.text_ndc(plotinfo["wks"], plotinfo["trstring"], vpx+0.95*vpw, vpy, plotinfo["trsres"])
    return(plotinfo)

def ngl_write_blstring(plotinfo):
    vpx = Ngl.get_float(plotinfo["plot"],"vpXF")
    vpy = Ngl.get_float(plotinfo["plot"],"vpYF")
    vpw = Ngl.get_float(plotinfo["plot"],"vpWidthF")
    vph = Ngl.get_float(plotinfo["plot"],"vpHeightF")
    if not "blsres" in plotinfo:
	plotinfo.update({"blsres":Ngl.Resources()})
	plotinfo["blsres"].txFontHeightF = 0.014
	plotinfo["blsres"].txJust        = "CenterLeft"
    if not "blstring" in plotinfo: 
 	print("No string for Bottom Left corner")
        plotinfo.update({"blstring":"____"})
    Ngl.text_ndc(plotinfo["wks"], plotinfo["blstring"], vpx+0.05*vpw, vpy-1.2*vph, plotinfo["blsres"])
    return(plotinfo)

def ngl_write_brstring(plotinfo):
    vpx = Ngl.get_float(plotinfo["plot"],"vpXF")
    vpy = Ngl.get_float(plotinfo["plot"],"vpYF")
    vpw = Ngl.get_float(plotinfo["plot"],"vpWidthF")
    vph = Ngl.get_float(plotinfo["plot"],"vpHeightF")
    if not "brsres" in plotinfo:
	plotinfo.update({"brsres":Ngl.Resources()})
	plotinfo["brsres"].txFontHeightF = 0.014
	plotinfo["brsres"].txJust        = "CenterRight"
    if not "brstring" in plotinfo: 
 	print("No string for Bottom Right corner")
        plotinfo.update({"brstring":"____"})
    Ngl.text_ndc(plotinfo["wks"], plotinfo["brstring"], vpx+0.95*vpw, vpy-1.2*vph, plotinfo["brsres"])
    return(plotinfo)

def write_info_text(plotinfo):
    if "tlstring" in plotinfo : plotinfo=write_tlstring(plotinfo)
    if "trstring" in plotinfo : plotinfo=write_trstring(plotinfo)
    if "blstring" in plotinfo : plotinfo=write_blstring(plotinfo)
    if "brstring" in plotinfo : plotinfo=write_brstring(plotinfo)
    return(plotinfo)

def ngl_add_info_string(plotinfo):
    #-- add units and copyrights to wks
    plot=plotinfo["plot"]
    wks=plotinfo["wks"]
    vpx = Ngl.get_float(plot,"vpXF")             #-- retrieve value of res.vpXF from plot
    vpy = Ngl.get_float(plot,"vpYF")             #-- retrieve value of res.vpYF from plot
    vpw = Ngl.get_float(plot,"vpWidthF")         #-- retrieve value of res.vpWidthF from plot
    vph = Ngl.get_float(plot,"vpHeightF")        #-- retrieve value of res.vpHeightF from plot

    print(plot,vpx,vpy,vpw,vph)

    units = "test_unit"
    txres = Ngl.Resources()
    txres.txFontHeightF = 0.014                 #-- font size for left, center and right string
    txres.txJust        = "CenterCenter"        #-- text justification
    Ngl.text_ndc(wks, "["+units+"]", vpx-0.04, vpy, txres)  #-- add text to wks

    txres.txFontHeightF = 0.012                 #-- font size for left, center and right string
    Ngl.text_ndc(wks,"footnote", vpx+vpw-0.07, vpy-vph, txres) #-- plot copyright info
    #plotinfo.update({"wks":wks})
    return(plotinfo)


def open_file(filename,mode="r",opt=None):
    f = Nio.open_file(filename,mode=mode,options=opt)
    return(f)

def get_var_data(infile,varname,slicestring=":"):
#    infile=open_file(filename)
    var = infile.variables[varname]
    print(var.dimensions)
    data=var[slicestring]
    return(data)

def time_slice_data(datain,dimpos=0,dim_indx=0):
    if dimpos is 0 : data=datain[dim_indx]
    return(data)

def level_slice_data(datain,dimpos=0,dim_indx=0):
    if dimpos is 0 : data=datain[dim_indx]
    return(data)


def get_fillval(dataset,fillvalfix=-999.99):
	try:fillval = dataset["_FillValue"]
	except:fillval = fillvalfix
	return(fillval)

def set_fillval(dataset):
	fillval=get_fillval(dataset)
    	dataset.update({"_FillValue":fillval})
	return(dataset)

def get_mask(dataset):
     if "data" in dataset:
	data=dataset["data"]
     else:
	if "zonal_data" in dataset:
	    data=dataset["zonal_data"]
	else:
	    if "merid_data" in dataset:
		data=dataset["merid_data"]
     mask=numpy.where(numpy.isnan(data))
     dataset.update({"mask":mask})
     return(dataset)

def mask_data(dataset,itemlist):
     for item in itemlist:
	if item in dataset:
	   data=dataset[item]
	   mask=dataset["mask"]
	   fillval=dataset["_FillValue"]	
	   data[mask]=fillval
	   dataset.update({item:data})
     return(dataset)

def fillval_mask(dataset):
	dataset=set_fillval(dataset)
	dataset=get_mask(dataset)
	dataset=mask_data(dataset,["data","zonal_data","merid_data"])
    	#data.fillna(fillval,inplace=True)
	#data=numpy.nan_to_num(data,copy=True,nan=fillval,posinf=None,neginf=None)
	return(dataset)

def get_2d_lat(datainfo,dataset):
    shape=dataset["data"].shape
    datadir=datainfo["datadir"]
    fname = datainfo["filename"]
    datafile=datadir + "/" + fname
    infile = Nio.open_file(datafile,"r")
    latnam = datainfo["latvar"]
    opt = datainfo["dimopt"]
    if opt == 1:
	lat = infile.variables[latnam][:]
	lat2d=pplib.get_2dlat(lat,shape)
    if opt == 2:
    	lat2d = infile.variables[latnam][:,:]
    dataset.update({"lat2d":lat2d})
    return(dataset)

def get_2d_lon(datainfo,dataset):
    shape=dataset["data"].shape
    datadir=datainfo["datadir"]
    fname = datainfo["filename"]
    datafile=datadir + "/" + fname
    infile = Nio.open_file(datafile,"r")
    lonnam = datainfo["lonvar"]
    opt = datainfo["dimopt"]
    if opt == 1:
	lon = infile.variables[lonnam][:]
	lon2d=pplib.get_2dlon(lon,shape)
    if opt == 2:
    	lon2d = infile.variables[lonnam][:,:]
    dataset.update({"lon2d":lon2d})
    return(dataset)

def get_vardims(datainfo,dataset):
    shape=dataset["data"].shape
    datadir=datainfo["datadir"]
    fname = datainfo["filename"]
    datafile=datadir + "/" + fname
    infile = Nio.open_file(datafile,"r")
    dimlist=datainfo["dimnames"]
    if len(shape) is not len(dimlist):
       err_string="Data shape missmatch with number of dimention: "+str(len(shape))+" not equal to "+str(len(dimlist))+" "
       print(err_string)
    else:
       for dimname in dimlist:
           dataset.update({dimname:infile.variables[dimname][:]})
    dataset.update({"dimnames":dimlist})
    return(dataset)

def get_slice_string(datainfo):
    keylist=["timeslice","levslice","latslice","lonslice"]
    slicestring=""
    for key in keylist:
        if key in datainfo: slicestring=slicestring+datainfo[key]+" "
    if slicestring is "" : slicestring=":"
    datainfo.update({"slicestring":slicestring})
    print(datainfo["slicestring"])
    return(datainfo)

def get_data(datainfo,vname=None):
    dataset={}
    datadir=datainfo["datadir"]
    fname = datainfo["filename"]
    if vname is None : vname = datainfo["fvarname"]
    if "grid_type" in datainfo:
       grid_type = datainfo["grid_type"]
    else:
       grid_type = None
    if "time_index" in datainfo: 
	timindx = datainfo["time_index"]
    else:
	timindx = None
    if "lev_index" in datainfo: 
	levindx = datainfo["lev_index"]
    else:
	levindx = None
    datafile=datadir + "/" + fname
    infile = Nio.open_file(datafile,"r")
    datainfo.update({"timedimname":infile.variables[vname].dimensions[datainfo["timepos"]]})
    datainfo.update({"latdimname":infile.variables[vname].dimensions[datainfo["latpos"]]})
    datainfo.update({"londimname":infile.variables[vname].dimensions[datainfo["lonpos"]]})
    datainfo.update({"dimnames":[datainfo["latdimname"],datainfo["londimname"]]})
    datainfo=get_slice_string(datainfo)
    slicestring=datainfo["slicestring"]
    data=get_var_data(infile,vname,slicestring)
    print(data)
    dataset.update({"data":data})
    dataset.update(infile.variables[vname].attributes)
    dataset=fillval_mask(dataset)
    dataset=get_vardims(datainfo,dataset)
    if grid_type is None:
       dataset = get_2d_lat(datainfo,dataset)
       dataset = get_2d_lon(datainfo,dataset)
    return(dataset)

def get_vector_magnitude(dataset):
    uwnd=dataset["zonal_data"]["data"]
    vwnd=dataset["merid_data"]["data"]
    usqr=uwnd.__pow__(2)
    vsqr=vwnd.__pow__(2)
    ssqr=usqr.__add__(vsqr)
    speed=numpy.sqrt(ssqr)
    dataset.update({"data":speed})
    return(dataset)

def get_vector_data(datainfo):
    dataset={}
    datadir=datainfo["datadir"]
    fname = datainfo["filename"]
    zonal_vname = datainfo["fvarname"][0]
    zonal_data=get_data(datainfo,vname=zonal_vname)
    merid_vname = datainfo["fvarname"][1]
    merid_data=get_data(datainfo,vname=merid_vname)
    dataset.update({"zonal_data":zonal_data})
    dataset.update({"merid_data":merid_data})
    for key in list(zonal_data.keys()):
        if key is not "data" : dataset.update({key:dataset["zonal_data"][key]})
    return(dataset)

def data_filename(datafield,datadir,init_date,init_hour,fcst_hour,datainfo=None,fname=None):
    if datainfo is None: datainfo=datadic.get_data_info(datafield)
    datainfo.update({"datadir":datadir})
    datainfo.update({"grid_type":"regular"})
    if fname is None:
    	fnprefix="ncum_imdaa_nearrealtime_12km_"
    	fnsufix="_"+str(init_date)+"_"+str(init_hour).zfill(2)+"z_hrofday.nc"
    	datainfo.update({"filename":fnprefix+str(datainfo["fnamekey"])+fnsufix})
    else:
	datainfo.update({"filename":fname})
    datainfo.update({"timeslice":"time|"+str(fcst_hour).zfill(2)})
    print(datainfo)
    return(datainfo)

def get_field_list():
    field_list=datadic.field_list
    return(field_list)


def umff_to_grib2(infile,outfile=None):
    if outfile is None: outfile=infile.split(".",1)[0]+".grib2"
    file_cubes=iris.load(infile)
    iris.save(file_cubes,outfile)

#############################################################################################################################
### Matlab Plotting Library based functions
#############################################################################################################################
def mpl_plot_keyfield(dataset,plotfile,tagmark="",lblst=[],text="",textpos=(0.25, -0.20)):
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

def mpl_plot_latlon(data,plotfile,tagmark="",lblst=[],text="",textpos=(0.25, -0.20),fltrkey="subtype",subtypenml=SUBTYPNML,title=""):
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
    fig= mpl_plot_cyl(datalist,fig,plot3,colors,area,alpha,parallels,meridians,tagmark=tagmark,lblst=lblst,text=text,textpos=textpos,title=title)
    #######
    pyplot.savefig(plotfile,bbox_inches='tight',dpi=300)
    return(fig)

def mpl_plot_location(data,plotfile,tagmark="",lblst=[],text="",textpos=(0.25, -0.20)):
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
    #fig = mpl_plot_ortho(data,fig,plot1,colors,area,alpha,parallels,meridians,polelat=10,polelon=70)
    #######
    #plot2=pyplot.subplot(122)
    #fig = mpl_plot_ortho(data,fig,plot1,colors,area,alpha,parallels,meridians,polelat=10,polelon=250)
    #######
    #plot1=pyplot.subplot(211)
    plot3=pyplot.subplot(212)
    fig= mpl_plot_cyl(data,fig,plot3,colors,area,alpha,parallels,meridians,tagmark=tagmark,lblst=lblst,text=text,textpos=textpos)
    #######
    pyplot.savefig(plotfile,bbox_inches='tight',dpi=200)
    return(fig)

def mpl_plot_merc(data,figure,plot,colors,area,alpha,parallels,meridians):
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

def mpl_plot_cyl(datalist,figure,plot,colors,area,alpha,parallels,meridians,tagmark="",lblst=[],text="",textpos=(0.25, -0.20),display_count=True,title=""):
    plot = Basemap(projection='cyl', resolution='c', llcrnrlat= -90.,urcrnrlat= 90.,llcrnrlon= -180.,urcrnrlon=180.)
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

def mpl_plot_ortho(data,figure,plot,colors,area,alpha,parallels,meridians,polelat,polelon):
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
    

def mpl_plot_shaded(plotfile,data,title,clevs=range(0,10,1),cpallet="jet",extend="max"):
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

def mpl_plot_button(plotfile,data,title,clevs=range(0,10,1),cpallet="jet",extend="max"):
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
    
def mpl_scatterplot(plotpath,odbnmlfile,data,cylcdatestr,obstype,varname,long_name,fieldname="Obsvalue",fill=True,extend="max"):
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

def mpl_gridplot(plotpath,odbnmlfile,data,cylcdatestr,obstype,varname,long_name,fieldname="Obsvalue",gridopt="mean",cpallet="jet",fill=False,extend="max"):
    data=obslib.odb_renamefield(data,odbnmlfile)
    clevs=obslib.clevgen(long_name,fieldname)
    units=obslib.dataunit(long_name)
    plot_title=obstype.replace("_"," ")+"\n"+long_name.replace("_"," ")+" ("+units+") "+fieldname+"\n"+" (1x1 gridded "+gridopt+") "+"\n"+cylcdatestr
    plotfile=plotpath+"/"+obstype+"_"+str(varname)+"_"+fieldname+"_"+cylcdatestr+".png"
    if gridopt is "count" : gridded_data=obslib.gridded_count_1x1deg(data,varname)
    elif gridopt is "sum" : gridded_data=obslib.gridded_sum_1x1deg(data,varname,fieldname)
    elif gridopt is "mean" : gridded_data=obslib.gridded_mean_1x1deg(data,varname,fieldname)
    if fill:
        mpl_plot_shaded(plotfile,gridded_data,plot_title,clevs,cpallet,extend)
    else:
        mpl_plot_button(plotfile,gridded_data,plot_title,clevs,cpallet,extend)
    return(fig)

def mpl_plot_density(plotpath,nmlfile,data,cylcdatestr,obstype,element,long_name,fill=False,extend="max"):
    fieldname="ObsDensity"
    gridopt="count"
    cpallet=mpl_truncate_colormap("gist_ncar_r", 0.05, 0.4)
    mpl_gridplot(plotpath,nmlfile,data,cylcdatestr,obstype,element,long_name,fieldname,gridopt,cpallet,fill,extend)

def mpl_plot_gridmean(plotpath,nmlfile,data,cylcdatestr,obstype,element,long_name,fill=False,extend="both"):
    fieldname="Obsvalue"
    gridopt="mean"
    cpallet=mpl_truncate_colormap("gist_ncar", 0.1, 0.8)
    mpl_gridplot(plotpath,nmlfile,data,cylcdatestr,obstype,element,long_name,fieldname,gridopt,cpallet,fill,extend)

def mpl_plot_depart_firstguess(plotpath,nmlfile,data,cylcdatestr,obstype,element,long_name,fill=False,extend="both"):
    fieldname="FGDep"
    gridopt="mean"
    cpallet=mpl_truncate_colormap("gist_ncar", 0.1, 0.8)
    mpl_gridplot(plotpath,nmlfile,data,cylcdatestr,obstype,element,long_name,fieldname,gridopt,cpallet,fill,extend)

def mpl_plot_depart_anal(plotpath,nmlfile,data,cylcdatestr,obstype,element,long_name,fill=False,extend="both"):
    fieldname="AnalDep"
    gridopt="mean"
    cpallet=mpl_truncate_colormap("gist_ncar", 0.1, 0.8)
    mpl_gridplot(plotpath,nmlfile,data,cylcdatestr,obstype,element,long_name,fieldname,gridopt,cpallet,fill,extend)

def mpl_globalview(plotfile,data,title,clevs=range(0,10,1),cpallet="jet",extend="max",plotmode="shaded"):
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


def mpl_plot_field(data,plotfile,domain="global",varname="hloswind",textout=False,clevs=None,dlevs=None):
	if clevs is None : clevs = [0,1,2,3,4,5,10,15,20,25,35,45]
	if dlevs is None : dlevs = [0,1,2,3,4,5,10]
	hgtmin=domaindic.dom_hgtmin[domain]
	hgtmax=domaindic.dom_hgtmax[domain]
	plotfile1=plotfile+"_"+domain+"_obs"
	title1="Observation (magnitude) at altitude ["+str(int(hgtmin))+"m to "+str(int(hgtmax))+"m ]"
	gridded_data_obs=obslib.gridded_rms_1x1deg(data,varname=varname,fieldname="observed")
	plot=mpl_globalview(plotfile1,gridded_data_obs,title=title1,clevs=clevs,extend="max",plotmode="shaded")
	if textout: obslib.obs_frame_ascii(gridded_data,plotfile1)
	plotfile2=plotfile+"_"+domain+"_bkg"
	title2="Background (magnitude) at altitude ["+str(int(hgtmin))+"m to "+str(int(hgtmax))+"m ]"
	gridded_data_bkg=obslib.gridded_rms_1x1deg(data,varname=varname,fieldname="background")
	plot=mpl_globalview(plotfile2,gridded_data_bkg,title=title2,clevs=clevs,extend="max",plotmode="shaded")
	if textout: obslib.obs_frame_ascii(gridded_data,plotfile2)
	plotfile3=plotfile+"_"+domain+"_omb"
	title3="OBS minus BKG (rms) at altitude ["+str(int(hgtmin))+"m to "+str(int(hgtmax))+"m ]"
	gridded_data_omb=obslib.gridded_rms_1x1deg(data,varname=varname,fieldname="obs_minus_bkg")
	plot=mpl_globalview(plotfile3,gridded_data_omb,title=title3,clevs=dlevs,extend="max",plotmode="shaded")
	if textout: obslib.obs_frame_ascii(gridded_data,plotfile3)
	return(gridded_data_omb)

def mpl_plot_indian(dataset,plot):
	data1=dataset["data1"]
	data2=dataset["data2"]
	a = xar_slice_indian(data1)
	b = xar_slice_indian(data2)
	#a = xarray.open_dataset(data1).sel(hybrid_ht=0,longitude=slice(60,100),latitude=slice(0,40))
	#b = xarray.open_dataset(data2).sel(hybrid_ht=0,longitude=slice(60,100),latitude=slice(0,40))
	var1 = a['q']
	var2 = b['q']
	diff = var2 - var1
	print(diff)
	fig = pyplot.figure(figsize=[12,5])
	ax = pyplot.axes(projection=ccrs.PlateCarree(central_longitude=0))
	diff.plot(ax=ax,vmin=-0.02, vmax=0.02, cmap='Blues', transform=ccrs.PlateCarree())
	ax.coastlines()
	ax.gridlines()
	pyplot.show()
	return(plot)

def mpl_plot_integrated_gl(dataset,plot):
	q_ctl=dataset["q_ctl"]
        rho_ctl=dataset["rho_ctl"]
	q_ctl = xarray.open_dataset(q_ctl)
	rho_ctl = xarray.open_dataset(rho_ctl)
	#thickness = xarray.DataArray(data=numpy.zeros(len(q_ctl.hybrid_ht)),dims=["hybrid_ht"],coords={"hybrid_ht": q_ctl.hybrid_ht},name="thickness")
	#thickness = thickness.isel(hybrid_ht=slice(None, -1))
	#for i in range(1, len(q_ctl.hybrid_ht) - 1):
    	#	thickness[i] = ((q_ctl.hybrid_ht[i] - q_ctl.hybrid_ht[i - 1]) / 2) + ((q_ctl.hybrid_ht[i + 1] - q_ctl.hybrid_ht[i]) / 2)
	#thickness[0] = (q_ctl.hybrid_ht[1] - q_ctl.hybrid_ht[0]) / 2
 	thickness=xar_layer_thickness(q_ctl)
	weighted_q = (q_ctl.q * thickness* rho_ctl['unspecified'].values)/82369
	integrated_weighted_q = weighted_q.sum('hybrid_ht')
	integrated_weighted_q.to_netcdf('integrated_weighted_q_ctl.nc')
	m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='c')
	fig = pyplot.figure(figsize=[12,5])
	ax = pyplot.axes(projection=ccrs.PlateCarree(central_longitude=0))
	integrated_weighted_q.plot(ax=ax,cmap='Blues', transform=ccrs.PlateCarree())
	m.drawcoastlines()
	pyplot.show()
	return(plot)

def mpl_truncate_colormap(cpallet, minval=0.0, maxval=1.0, n=100):
    cmap=mplcm.get_cmap(cpallet)
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(numpy.linspace(minval, maxval, n)))
    return new_cmap

def mpl_plot_ose_scalar(plotdic):
	data_ctl=plotdic["data_ctl"]
	data_exp=plotdic["data_exp"]
	plotfile=plotdic["plotfile"]
	axlbl_y_ctl=plotdic["ctlname"]
	axlbl_y_exp=plotdic["expname"]

	integrated_q_diff = data_exp - data_ctl

	fig, axes = pyplot.subplots(nrows=3, ncols=1,figsize=[8,20], subplot_kw={'projection': ccrs.PlateCarree(central_longitude=0)})
	plot=[None]*3
	axlbly=[None]*3

	plot[0]=data_ctl.plot(ax=axes[0], cmap='Blues', transform=ccrs.PlateCarree())
	m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='c',ax=axes[0])
	m.drawcoastlines()
	axlbly[0]=axes[0].text(-0.1, 0.5, axlbl_y_ctl, va='center', ha='center', rotation='vertical', transform=axes[0].transAxes)
	axins = inset_axes(axes[0], width = "5%", height = "100%", loc = 'lower left', bbox_to_anchor = (1.09, 0., 1, 1), bbox_transform = axes[0].transAxes, borderpad = 0)
	fig.colorbar(plot[0], cax = axins)	

	plot[1]=data_exp.plot(ax=axes[1], cmap='Blues', transform=ccrs.PlateCarree())
	m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='c',ax=axes[1])
	m.drawcoastlines()
	axlbly[1]=axes[1].text(-0.1, 0.5, axlbl_y_ctl, va='center', ha='center', rotation='vertical', transform=axes[1].transAxes)
	axins = inset_axes(axes[1], width = "5%", height = "100%", loc = 'lower left', bbox_to_anchor = (1.09, 0., 1, 1), bbox_transform = axes[1].transAxes, borderpad = 0)
	fig.colorbar(plot[1], cax = axins)	

	plot[2]=integrated_q_diff.plot(ax=axes[2],vmin=-6,vmax=6, cmap='RdBu_r', transform=ccrs.PlateCarree())
	m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,llcrnrlon=-180,urcrnrlon=180,resolution='c',ax=axes[2])
	m.drawcoastlines()
	axlbly[2]=axes[2].text(-0.1, 0.5, 'EXP-CTL', va='center', ha='center', rotation='vertical', transform=axes[2].transAxes)
	axins = inset_axes(axes[2], width = "5%", height = "100%", loc = 'lower left', bbox_to_anchor = (1.09, 0., 1, 1), bbox_transform = axes[2].transAxes, borderpad = 0)
	fig.colorbar(plot[2], cax = axins)	

	#im = axes[0].imshow(data_ctl)	

	pyplot.tight_layout(pad=10)
	pyplot.savefig(plotfile)
	return(plotfile)



#############################################################################################################################
### NCAR Graghics Library based functions
#############################################################################################################################


def ngl_plot_filename(datafield,plotdir,init_date,fcst_hour,prefix="",sufix=""):
    plotinfo={}
    sufix="_"+str(fcst_hour).zfill(2)+"Z"
    plotinfo.update({"plotfile":plotdir+"/"+prefix+datafield+sufix})
    plotinfo.update({"wks_type":"png"})
    plotinfo.update({"trstring":str(init_date)+"_"+str(fcst_hour).zfill(2)+"Z"})
    plotinfo=ngl_gen_wks(plotinfo)
    return(plotinfo)


def ngl_test_ngl():
    wkres = Ngl.Resources()
    rval=subprocess.call("export PYNGL_COLORMAPS=$HOME", shell=True)
    wkres.wkColorMap = "gibies_colourmap_20150117"
    wks = Ngl.open_wks("png","test",wkres)
    Ngl.draw_colormap(wks)
    Ngl.end()


def ngl_plot_all(datadir,plotdir,init_date,init_hour,fcst_lead,field_list=None,fname=None):
    hoursofday=numpy.arange(init_hour,init_hour+fcst_lead,1)
    if field_list is None : field_list=get_field_list()
    for datafield in field_list:
        for fcst_hour in hoursofday:
            print(fcst_hour)
            print(datafield)
            plotfile=plot_data([datafield],datadir,plotdir,init_date,init_hour,fcst_hour,fname)
    Ngl.end()
    return(plotfile)

def ngl_vcstyle(plotinfo):
    res=resof(plotinfo)
    if "vcstyle" in plotinfo : 
       vcstyle=plotinfo["vcstyle"]
    else : 
       vcstyle="CurlyVector"
    if vcstyle is "CurlyVector":
        res.vcLineArrowThicknessF = 3.0
	res.vcMonoLineArrowColor = False  
    if vcstyle is "FillArrow" :
        res.vcFillArrowsOn       = True
        res.vcMonoFillArrowFillColor  =  False
        res.vcFillArrowWidthF         = 0.1
    res.vcGlyphStyle = vcstyle
    plotinfo.update({"res":res})
    return(plotinfo) 

def ngl_vcmag(plotinfo):
    res=resof(plotinfo)
    if plotinfo["plot_type"] in ["vector","vector_scalar",]:
    	if "cmap" in plotinfo : res.vcLevelPalette = plotinfo["cmap"]
	res.vcMinFracLengthF     = 0.008   # Increase length of vector
	res.vcRefLengthF         = 0.045
	if "vcmag" in plotinfo : 
           res.vcRefMagnitudeF      = plotinfo["vcmag"]
        else :
           res.vcRefMagnitudeF      = 2.0
    plotinfo.update({"res":res})
    return(plotinfo)

def ngl_plot_data(field_list,datadir,plotdir,init_date,init_hour,fcst_hour,fname=None):
    plotinfo=plot_filename(field_list[0],plotdir,init_date,fcst_hour,prefix="ncumda_g12_imdaa_",sufix="_"+str(fcst_hour).zfill(2)+"Z")
    print(plotinfo)
    plot = []
    for count,datafield in enumerate(field_list):
        plotinfo=datadic.get_plot_info(datafield,plotinfo)
        datainfo=data_filename(datafield,datadir,init_date,init_hour,fcst_hour,fname=fname)
        plotinfo=get_cmap(plotinfo)
        res=resof(plotinfo)
        plotinfo.update({"res":res})
	if plotinfo["plot_type"] in ["contour"]:
		dataset=get_data(datainfo)
	else:
		dataset=get_vector_data(datainfo)
	plot1=get_subplot(dataset,plotinfo)
	plot.append(plot1)
    count=len(plot)
    panelres                  =  Ngl.Resources()
    panelres.nglPanelLabelBar =  False           #-- common labelbar
    panelres.nglPanelYWhiteSpacePercent =  0        #-- reduce space between the panel plots
    panelres.nglPanelXWhiteSpacePercent =  0        #-- reduce space between the panel plots
    panelres.nglPanelTop      =  0.95               #-- top position of panel
    #-- create the panel
    Ngl.panel(plotinfo["wks"],plot,[count,1],panelres)
    outfile=check_plot(plotinfo)
    del plotinfo, datainfo, plot
    print(outfile)
    return(outfile)

def ngl_get_subplot(dataset,plotinfo):
    wks=plotinfo["wks"]
    res=resof(plotinfo)
#   lat2d=dataset["lat2d"]
#   lon2d=dataset["lon2d"]
    if "title" in plotinfo: 
	title=plotinfo["title"]
    else:
	title=""
    if "draworder" in plotinfo:
    	res.mpFillDrawOrder   = plotinfo["draworder"]

    if plotinfo["plot_type"] in ["contour","vector_scalar",]:
    	res.cnFillOn = True
        if "cnline" in plotinfo and plotinfo["cnline"] is not None :
		res.cnLevelFlags = plotinfo["cnline"]
        else:
   		res.cnLineOn = False
		res.cnLevelFlags = "NoLine"
    		res.cnLineLabelsOn        = False
    		res.cnInfoLabelOn         = False
    	res.cnLevelSelectionMode = "ExplicitLevels"
    	if "cnlev" in plotinfo : res.cnLevels 	= plotinfo["cnlev"]
        if "cmap" in plotinfo :	res.cnFillPalette = plotinfo["cmap"]


    res.lbOrientation  = "Horizontal"
#    res.sfXArray     = lon2d
#    res.sfYArray     = lat2d
    res.sfMissingValueV = get_fillval(dataset)
    res.gsnAddCyclic = False  
    res.tiMainString = title
    if "gsnLeftString" in plotinfo:
    	res.gsnLeftString = plotinfo["gsnLeftString"]
    else:
	res.gsnLeftString = ""
    plotinfo.update({"wks":wks})
    plotinfo.update({"res":res})
    plot1=cyl_plot(dataset,plotinfo)
    return(plot1)

def ngl_cyl_plot(dataset,plotinfo):
    wks=plotinfo["wks"]
    res=resof(plotinfo)
    res.nglDraw               = False               #-- don't draw individual plots
    res.nglFrame              = False               #-- don't advance frame
    res.mpDataBaseVersion      = "MediumRes"
    res.mpFillOn              = True
    res.mpFillColors = [0,-1,-1,-1]
    res.mpFillAreaSpecifiers  = ["land"]
    res.mpSpecifiedFillColors = ["gray65"]
    res.mpGridLatSpacingF =  10.                  #-- grid lat spacing
    res.mpGridLonSpacingF =  10.                  #-- grid lon spacing
    res.mpLimitMode       = "LatLon"              #-- must be set using minLatF/maxLatF/minLonF/maxLonF
    print(dataset)   
    if "latitude" in dataset["dimnames"]:
       res.mpMinLatF = min(dataset["latitude"])
       res.mpMaxLatF = max(dataset["latitude"])
       if plotinfo["plot_type"] in ["contour","vector_scalar"]:
          res.sfYCStartV = float(min(dataset["latitude"]))
          res.sfYCEndV   = float(max(dataset["latitude"]))
       if plotinfo["plot_type"] in ["vector","vector_scalar"]:
          res.vfYCStartV = float(min(dataset["latitude"]))
          res.vfYCEndV   = float(max(dataset["latitude"]))
    else:
       res.mpMinLatF    = -90
       res.mpMaxLatF    = 90
       
    if "longitude" in dataset["dimnames"]:
       res.mpMinLonF = min(dataset["longitude"])
       res.mpMaxLonF = max(dataset["longitude"])
       if plotinfo["plot_type"] in ["contour","vector_scalar"]:
          res.sfXCStartV = float(min(dataset["longitude"]))   # Define where contour plot
          res.sfXCEndV   = float(max(dataset["longitude"]))   # should lie on the map plot.
       if plotinfo["plot_type"] in ["vector","vector_scalar"]:
          res.vfXCStartV = float(min(dataset["longitude"]))   # Define where contour plot
          res.vfXCEndV   = float(max(dataset["longitude"]))   # should lie on the map plot.
    else:
       res.mpMinLonF    = 0
       res.mpMaxLonF    = 360
      
    if plotinfo["plot_type"] in ["vector","vector_scalar",]:
        plotinfo=vcmag(plotinfo)
        plotinfo=vcstyle(plotinfo)
        res=resof(plotinfo)
        res.vcMinDistanceF = 0.01
	uwnd=dataset["zonal_data"]["data"]
	vwnd=dataset["merid_data"]["data"]
	print(uwnd.shape)
	print(vwnd.shape)
    if plotinfo["plot_type"] in ["vector_scalar",]:
    	data=get_vector_magnitude(dataset)["data"]
	print(data.shape)
    	plotinfo["plot"]=Ngl.vector_scalar_map(wks, uwnd, vwnd, data, res)
    if plotinfo["plot_type"] in ["vector",]:
	res.vcMonoLineArrowColor = False  # Draw vectors in color.
    	plotinfo["plot"]=Ngl.vector_map(wks, uwnd, vwnd, res)
    if plotinfo["plot_type"] in ["contour",]:
	data=dataset["data"]
	print(data.shape)
	plotinfo["plot"]=Ngl.contour_map(wks,data,res)
    plotinfo=write_info_text(plotinfo)
    return(plotinfo["plot"])
	

#############################################################################################################################
### XARRAY based functions
#############################################################################################################################


def get_dataset(filename,time=0,lev=0,lon_min=0,lon_max=-1,lon_skip=None,lat_min=0,lat_max=-1,lat_skip=None):
	file_in = xarray.open_dataset(gdf.get(filename))
	dataset = file_in.isel(time=time, lev=lev, lon=slice(lon_min, lon_max, lon_skip), lat=slice(lat_min, lat_max, lat_skip))
	return(datasset)

def xar_slice_indian(data):
	ind_data=xarray.open_dataset(data).sel(hybrid_ht=0,longitude=slice(60,100),latitude=slice(0,40))
	return(ind_data)

def xar_layer_thickness(q_ctl):
	thickness = xarray.DataArray(data=numpy.zeros(len(q_ctl.hybrid_ht)),dims=["hybrid_ht"],coords={"hybrid_ht": q_ctl.hybrid_ht},name="thickness")
	thickness = thickness.isel(hybrid_ht=slice(None, -1))
	for i in range(1, len(q_ctl.hybrid_ht) - 1):
	    thickness[i] = ((q_ctl.hybrid_ht[i] - q_ctl.hybrid_ht[i - 1]) / 2) + ((q_ctl.hybrid_ht[i + 1] - q_ctl.hybrid_ht[i]) / 2)
	thickness[0] = (q_ctl.hybrid_ht[1] - q_ctl.hybrid_ht[0]) / 2
	return(thickness)

def xar_quot_rsqure(rho_ctl):
	R_ctl = rho_ctl.hybrid_ht + 6371*(10**3)
	R_square_ctl = R_ctl**2
	rho_ctl['density'] = rho_ctl['unspecified'] / R_square_ctl
	return(rho_ctl)

def xar_qrhodh(q_ctl,rho_ctl):
	if "density" not in rho_ctl.data_vars:
		rho_ctl=xar_quot_rsqure(rho_ctl)
	thickness=xar_layer_thickness(q_ctl)
	weighted_q_ctl = q_ctl.q * thickness*rho_ctl['density'].values
	return(weighted_q_ctl)

def xar_ipw(q_ctl,rho_ctl):
	weighted_q_ctl = xar_qrhodh(q_ctl,rho_ctl)
	data_ctl = weighted_q_ctl.sum('hybrid_ht')
	return(data_ctl)

def xar_orthogrid_wind(u_wind_ctl,v_wind_ctl):
	v_wind_ctl_interp = v_wind_ctl.interp(latitude=u_wind_ctl.latitude, longitude=u_wind_ctl.longitude)
	return(v_wind_ctl_interp)

def xar_qtransdh(q_ctl,rho_ctl,u_wind_ctl):
	if "density" not in rho_ctl.data_vars:
		rho_ctl=xar_quot_rsqure(rho_ctl)
	thickness=xar_layer_thickness(q_ctl)
	weighted_q_u_ctl = q_ctl.q * thickness*rho_ctl['density'].values*u_wind_ctl['u'].values
	return(weighted_q_u_ctl)

