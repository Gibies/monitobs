#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 30 04:09:45 2021

@author: gibies
"""
import os,sys
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
PKGHOME=os.environ.get('PKGHOME',os.path.dirname(CURR_PATH))
PKGNAME=os.path.basename(PKGHOME)
LIB=os.environ.get('LIB',PKGHOME+"/pylib")
sys.path.append(LIB)
import daview
import datetime
import pandas
import numpy
import re

fnamekey={
"default":"isobaric",
"wet_soil":"soil",
"t_sfc":"screen",
"temp_sfc":"screen",
"temp_skin":"pbl",
"hght_pbl":"pbl",
"roughness":"pbl",
"pres_sfc":"screen",
"pres_msl":"msl",
"rhum_sfc":"screen",
"vis_sfc":"screen",
"cld_high":"cloud",
"cld_med":"cloud",
"cld_low":"cloud",
"cld_vrlo":"cloud",
"cld_tot":"cloud",
"uwnd_850hpa":"wind",
"wind_10m":"pole",
"wind_50m":"tower",
"wind_850hpa":"wind",
"wind_500hpa":"wind",
"wind_200hpa":"wind",
}

fvarname={
"default":["u","v"],
"wet_soil":"sm",
"msl_pressure":"p",
"temp":"temp",
"temp_sfc":"temp",
"temp_skin":"temp",
"hght_pbl":"blht",
"roughness":"ROHLENABL",
"pres":"p",
"pres_sfc":"p",
"pres_msl":"p",
"rhum":"rh",
"rhum_sfc":"rh",
"vis_sfc":"VIS",
"cld_high":"HCDC",
"cld_med":"MCDC",
"cld_low":"LCDC",
"cld_vrlo":"VLCDC",
"cld_tot":"TCDCRO",
"uwnd_850hpa":"u",
"wind_10m":["u","v"],
"wind_50m":["u","v"],
"wind_850hpa":["u","v"],
}

timeslice={
"default":"time|i0",
"wet_soil":"time|i0",
"cloud":"time|i0",
"temp_sfc":"time|i0",
"temp_skin":"time|i0",
"pres_sfc":"time|i0",
"pres_msl":"time|i0",
"rhum_sfc":"time|i0",
"vis_sfc":"time|i0",
"cld_high":"time|i0",
"cld_med":"time|i0",
"cld_low":"time|i0",
"cld_vrlo":"time|i0",
"cld_tot":"time|i0",
"uwnd_850hpa":"time|i0",
"wind_10m":"time|i0",
"wind_50m":"time|i0",
"wind_850hpa":"time|i0",
}

levslice={
"default":"p|850",
"850hpa":"p|850",
"500hpa":"p|500",
"200hpa":"p|200",
"wet_soil":"level6|i0",
"cloud":"unspecified|i0",
"temp_sfc":"ht|i0",
"temp_skin":"surface|i0",
"pres_sfc":"surface|i0",
"pres_msl":"msl|i0",
"rhum_sfc":"ht|i0",
"vis_sfc":"ht|i0",
"cld_high":"unspecified|i0",
"cld_med":"unspecified|i0",
"cld_low":"unspecified|i0",
"cld_vrlo":"unspecified|i0",
"cld_tot":"unspecified|i0",
"uwnd_850hpa":"p|850",
"wind_10m":"ht|i0",
"wind_50m":"ht|i0",
"wind_850hpa":"p|850",
"wind_500hpa":"p|500",
"wind_200hpa":"p|200",
}


latslice={
"default":"latitude|:",
"wet_soil":"latitude|:",
"cloud":"latitude|:",
"temp_sfc":"latitude|:",
"temp_skin":"latitude|:",
"pres_sfc":"latitude|:",
"pres_msl":"latitude|:",
"rhum_sfc":"latitude|:",
"vis_sfc":"latitude|:",
"cld_high":"latitude|:",
"cld_med":"latitude|:",
"cld_low":"latitude|:",
"cld_vrlo":"latitude|:",
"cld_tot":"latitude|:",
"uwnd_850hpa":"{latitude|5:15}",
"wind_10m":"latitude|:",
"wind_50m":"latitude|:",
"wind_850hpa":"latitude|:",
}

lonslice={
"default":"longitude|:",
"wet_soil":"longitude|:",
"cloud":"longitude|:",
"temp_sfc":"longitude|:",
"temp_skin":"longitude|:",
"pres_sfc":"longitude|:",
"pres_msl":"longitude|:",
"rhum_sfc":"longitude|:",
"vis_sfc":"longitude|:",
"cld_high":"longitude|:",
"cld_med":"longitude|:",
"cld_low":"longitude|:",
"cld_vrlo":"longitude|:",
"cld_tot":"longitude|:",
"uwnd_850hpa":"{longitude|60:80}",
"wind_10m":"longitude|:",
"wind_50m":"longitude|:",
"wind_850hpa":"longitude|:",
}

timepos={
"default":0,
"wet_soil":0,
"cloud":0,
"msl_pressure":0,
"temp_sfc":0,
"temp_skin":0,
"pres_sfc":0,
"pres_msl":0,
"rhum_sfc":0,
"vis_sfc":0,
"cld_high":0,
"cld_med":0,
"cld_low":0,
"cld_vrlo":0,
"cld_tot":0,
"uwnd_850hpa":0,
"wind_10m":0,
"wind_50m":0,
"wind_850hpa":0,
}

latpos={
"default":2,
"wet_soil":2,
"cloud":2,
"msl_pressure":2,
"temp_sfc":2,
"temp_skin":2,
"pres_sfc":2,
"pres_msl":2,
"rhum_sfc":2,
"vis_sfc":2,
"cld_high":2,
"cld_med":2,
"cld_low":2,
"cld_vrlo":2,
"cld_tot":2,
"uwnd_850hpa":2,
"wind_10m":2,
"wind_50m":2,
"wind_850hpa":2,
}

lonpos={
"default":3,
"wet_soil":3,
"cloud":3,
"msl_pressure":3,
"temp_sfc":3,
"temp_skin":3,
"pres_sfc":3,
"pres_msl":3,
"rhum_sfc":3,
"vis_sfc":3,
"cld_high":3,
"cld_med":3,
"cld_low":3,
"cld_vrlo":3,
"cld_tot":3,
"uwnd_850hpa":3,
"wind_10m":3,
"wind_50m":3,
"wind_850hpa":3,
}

cnlev = {
"cld_tot": numpy.arange(5., 10., 1.),
"sst" : numpy.arange(0.,32.,2.),
"airt850" : numpy.arange(280.,300.,1.),
"vis_sfc": numpy.arange(0.,1000.,50.),
"rhum": numpy.arange(50.,100.,5.),
"rhum_sfc": numpy.arange(50.,100.,5.),
"temp": numpy.arange(270.,320.,2.),
"temp_sfc": numpy.arange(270.,320.,2.),
"uwnd_850hpa": numpy.arange(5.,30.,2.),
"temp_850hpa": numpy.arange(260.,300.,2.),
"temp_500hpa": numpy.arange(240.,280.,2.),
"temp_200hpa": numpy.arange(200.,230.,2.),
#"temp_skin": numpy.arange(270.,310.,2.),
"temp_skin": numpy.arange(270.,320.,2.),
#"pres_sfc": numpy.arange(70000.,102000.,500.),
"pres_sfc": numpy.array([70000.,80000.,90000.,95000.,97000.,98000.,99000.,99500.,99700.,99800.,99900.,100000.,100500.,100600.,100700.,100800.,100900.,101000.,101100.,101200.,101300.,101400.,101500.,102000.,]),
"pres_msl": numpy.array([100000.,100500.,100600.,100700.,100800.,100900.,101000.,101100.,101200.,101300.,101400.,101500.,101600.,101700.,101800.,101900.,102000.,102500.,103000.,]),
#"pres_msl": numpy.arange(100100.,102000.,100.),
	}

colrmin = {
"cld": 0,
"sst" : 30,
"airt850" : 30,
"wet_soil":40,
"cld_high": 0,
"cld_med": 0,
"cld_low": 0,
"cld_vrlo": 0,
"cld_tot": 0,
	}

colrmax = {
"cld": 50,
"sst" : 170,
"airt850" : 170,
"wet_soil":140,
"cld_high": 50,
"cld_med": 50,
"cld_low": 50,
"cld_vrlo": 50,
"cld_tot": 50,
"rhum": 50,
	}

revcmap = {
"cld": True,
"wet_soil":True,
"cld_high": True,
"cld_med": True,
"cld_low": True,
"cld_vrlo": True,
"cld_tot": True,
"vis_sfc": True,
"pres_msl": True,
}

cnline = {
#Options: NoLine (default) LineOnly LabelOnly LineAndLabel

}

vcmag = {
"wind":20,
"wind_10m":20,
"wind_50m":20,
"wind_850hpa":20,
"wind_500hpa":20,
"wind_200hpa":20,
}

draworder = {
	"sst" : "Postdraw",
	"airt850" : "Predraw",
	}

writenc = {
	"sst" : False,
	"airt850" : True,
	}

grpkeylist=["wind","temp","rhum","pres","cld","sfc","10m","50m","850hpa","500hpa","200hpa"]

data_info_key_list=["fnamekey","fvarname","timeslice","levslice","latslice","lonslice","timepos","latpos","lonpos",]

plot_info_key_list=["cnlev","colrmin","colrmax","revcmap","cnline","draworder","vcmag"]


field_list=[
"cld_high",
"cld_med",
"cld_low",
"cld_vrlo",
"cld_tot",
"wet_soil",
"pres_sfc",
"temp_sfc",
"temp_skin",
"rhum_sfc",
"vis_sfc",
"pres_msl",
"wind_10m",
"wind_50m",
"uwnd_850hpa",
"wind_850hpa",
"wind_500hpa",
"wind_200hpa",
"temp_850hpa",
"temp_500hpa",
"temp_200hpa",
"rhum_850hpa",
"rhum_500hpa",
"rhum_200hpa",
]

def update(diction,item,field):
    if field in globals()[item] :
       diction.update({item:globals()[item][field]})
    if item not in diction :
       for grpkey in grpkeylist:
           if re.search(grpkey,field):
               if grpkey in globals()[item]:
                  diction.update({item:globals()[item][grpkey]})
    if item not in diction :
       if "default" in globals()[item] :
           diction.update({item:globals()[item]["default"]})
    return(diction)

def get_data_info(datafield):
    datainfo={}
    for dik in data_info_key_list:
        datainfo=update(datainfo,dik,datafield)
    return(datainfo)


def get_plot_info(datafield,plotinfo):
    plotinfo.update({"tlstring":datafield})
    for pik in plot_info_key_list:
        plotinfo=update(plotinfo,pik,datafield)
    vckeylist = [ "wind", ]
    for vckey in vckeylist :
       if re.search(vckey,datafield) :
          plotinfo.update({"plot_type":"vector_scalar"})
    if "plot_type" not in plotinfo :
       plotinfo.update({"plot_type":"contour"})
    return(plotinfo)


