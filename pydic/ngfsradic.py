#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 05:57:45 2020

@author: gibies
"""
import sys
import os
MONITOBS=os.environ.get('MONITOBS',"../")
OBSLIB=os.environ.get('OBSLIB',MONITOBS+"/pylib")
sys.path.append(OBSLIB)
import obslib
import obsdic
import datetime
import pandas 
import netCDF4


SYMOBS=os.environ.get('SYMOBS',True)
diaglev=int(os.environ.get('GEN_MODE',0))

ngfsrapath="/home/gibies/data/ngfsra"

def filename(element,year,month,day=None):
    fileprefix="ngfs_reanl_"
    datafile = {
        "gph" : "%s/%s"%(ngfsrapath,fileprefix+"HR_HGT-prl_"+str(year)+str(month)+"0100-"+str(year)+str(month)+"3118.nc"),
        "tmp" : "%s/%s"%(ngfsrapath,fileprefix+"HR_TMP-prl_"+str(year)+str(month)+"0100-"+str(year)+str(month)+"3118.nc"),
        "shum" : "%s/%s"%(ngfsrapath,fileprefix+"HR_SPFH-2m-agl_"+str(year)+str(month)+"0100-"+str(year)+str(month)+"3118.nc"),
        "rhum" : "%s/%s"%(ngfsrapath,fileprefix+"HR_RH-prl_"+str(year)+str(month)+"0100-"+str(year)+str(month)+"3118.nc"),
        "uwnd" : "%s/%s"%(ngfsrapath,fileprefix+"HR_UGRD-prl_"+str(year)+str(month)+"0100-"+str(year)+str(month)+"3118.nc"),
        "vwnd" : "%s/%s"%(ngfsrapath,fileprefix+"HR_VGRD-prl_"+str(year)+str(month)+"0100-"+str(year)+str(month)+"3118.nc"),
        "tsfc" : "%s/%s"%(ngfsrapath,fileprefix+"HR_TMP-2m-agl_"+str(year)+str(month)+"0100-"+str(year)+str(month)+"3118.nc"),
        "rhsfc" : "%s/%s"%(ngfsrapath,fileprefix+"HR_RH-2m-agl_"+str(year)+str(month)+"0100-"+str(year)+str(month)+"3118.nc"),
        "usfc" : "%s/%s"%(ngfsrapath,fileprefix+"HR_UGRD-2m-agl_"+str(year)+str(month)+"0100-"+str(year)+str(month)+"3118.nc"),
        "vsfc" : "%s/%s"%(ngfsrapath,fileprefix+"HR_VGRD-2m-agl_"+str(year)+str(month)+"0100-"+str(year)+str(month)+"3118.nc"),
        "psfc" : "%s/%s"%(ngfsrapath,fileprefix+"HR_PRES-sfc_"+str(year)+str(month)+"0100-"+str(year)+str(month)+"3118.nc"),
        "slp" : "%s/%s"%(ngfsrapath,fileprefix+"HR_PRES-msl_"+str(year)+str(month)+"0100-"+str(year)+str(month)+"3118.nc"),
        "t2m" : "%s/%s"%(ngfsrapath,fileprefix+"HR_TMP-2m-agl_"+str(year)+str(month)+"0100-"+str(year)+str(month)+"3118.nc"),
        "rh2m" : "%s/%s"%(ngfsrapath,fileprefix+"HR_RH-2m-agl_"+str(year)+str(month)+"0100-"+str(year)+str(month)+"3118.nc"),
        "sh2m" : "%s/%s"%(ngfsrapath,fileprefix+"HR_SPFH-2m-agl_"+str(year)+str(month)+"0100-"+str(year)+str(month)+"3118.nc"),
        "u2m" : "%s/%s"%(ngfsrapath,fileprefix+"HR_GUST-sfc_"+str(year)+str(month)+"0100-"+str(year)+str(month)+"3118.nc"),
        "v2m" : "%s/%s"%(ngfsrapath,fileprefix+"HR_GUST-sfc_"+str(year)+str(month)+"0100-"+str(year)+str(month)+"3118.nc"),
        "u10m" : "%s/%s"%(ngfsrapath,fileprefix+"HR_UGRD-10m-agl_"+str(year)+str(month)+"0100-"+str(year)+str(month)+"3118.nc"),
        "v10m" : "%s/%s"%(ngfsrapath,fileprefix+"HR_VGRD-10m-agl_"+str(year)+str(month)+"0100-"+str(year)+str(month)+"3118.nc"),
        }
    return(datafile[element])

datavar = {
        "time"	: "time",
        "lev"	: "level",
        "lat"	: "latitude",
        "lon"	: "longitude",
        "gph"	: "HGT",
        "pres"	: "pres",
        "tmp"	: "TMP",
        "rhum"	: "RH",
        "shum"	: "SPFH",
        "uwnd"	: "UGRD",
        "vwnd"	: "VGRD",
        "slp"	: "PRES_meansealevel",
	"psfc"	: "PRES_surface",
	"t2m"	: "TMP_2maboveground",
	"rh2m"	: "RH_2maboveground",
	"sh2m"	: "SPFH_2maboveground",
        "u2m"	: "UGRD_2maboveground",
        "v2m"	: "VGRD_2maboveground",
        "u10m"	: "UGRD_10maboveground",
        "v10m"	: "VGRD_10maboveground",
        }

def getdata(var="time",element="gph",year=None,month=None,day=None):
    print(var,element)
    with netCDF4.Dataset(filename(element,year,month,day), 'r',  format='NETCDF4_CLASSIC') as fileptr : 
        data = fileptr.variables[datavar[var]][:]
    return(data)
    
def getunits(var="time",element="gph",year=None,month=None,day=None):
    with netCDF4.Dataset(filename(element,year,month,day), 'r',  format='NETCDF4_CLASSIC') as fileptr : 
        units = fileptr.variables[datavar[var]].units
    return(units)



