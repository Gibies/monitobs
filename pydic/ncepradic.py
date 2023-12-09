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
import pandas as pd
import netCDF4


SYMOBS=os.environ.get('SYMOBS',True)
diaglev=int(os.environ.get('GEN_MODE',0))

netcdfpath="/home/gibies/DATA/NCEP_RA_6hr"

def filename(Year,element):
    datafile = {
        "gph" : "%s/%s"%(netcdfpath,"hgt."+str(Year)+".nc"),
        "tmp" : "%s/%s"%(netcdfpath,"air."+str(Year)+".nc"),
        "rhum" : "%s/%s"%(netcdfpath,"rhum."+str(Year)+".nc"),
        "uwnd" : "%s/%s"%(netcdfpath,"uwnd."+str(Year)+".nc"),
        "vwnd" : "%s/%s"%(netcdfpath,"vwnd."+str(Year)+".nc"),
        "tsfc" : "%s/%s"%(netcdfpath,"air.2m.gauss."+str(Year)+".nc"),
        "rhsfc" : "%s/%s"%(netcdfpath,"rhum.2m.gauss."+str(Year)+".nc"),
        "usfc" : "%s/%s"%(netcdfpath,"uwnd.2m.gauss."+str(Year)+".nc"),
        "vsfc" : "%s/%s"%(netcdfpath,"vwnd.2m.gauss."+str(Year)+".nc"),
        "psfc" : "%s/%s"%(netcdfpath,"pres.sfc."+str(Year)+".nc"),
        "slp" : "%s/%s"%(netcdfpath,"slp."+str(Year)+".nc"),
        "t2m" : "%s/%s"%(netcdfpath,"air.2m.gauss."+str(Year)+".nc"),
        "rh2m" : "%s/%s"%(netcdfpath,"rhum.2m.gauss."+str(Year)+".nc"),
        "sh2m" : "%s/%s"%(netcdfpath,"shum.2m.gauss."+str(Year)+".nc"),
        "u2m" : "%s/%s"%(netcdfpath,"uwnd.2m.gauss."+str(Year)+".nc"),
        "v2m" : "%s/%s"%(netcdfpath,"vwnd.2m.gauss."+str(Year)+".nc"),
        "u10m" : "%s/%s"%(netcdfpath,"uwnd.10m.gauss."+str(Year)+".nc"),
        "v10m" : "%s/%s"%(netcdfpath,"vwnd.10m.gauss."+str(Year)+".nc"),
        }
    return(datafile[element])

datavar = {
        "time" : "time",
        "lev" : "level",
        "lat" : "lat",
        "lon" : "lon",
        "gph" : "hgt",
        "slp" : "slp",
        "pres" : "pres",
        "tmp" : "air",
        "rhum" : "rhum",
        "shum" : "shum",
        "uwnd" : "uwnd",
        "vwnd" : "vwnd",
        }

def getdata(Year,var="time",element="gph"):
    with netCDF4.Dataset(filename(Year,element), 'r',  format='NETCDF4_CLASSIC') as fileptr : 
        data = fileptr.variables[datavar[var]][:]
    return(data)
    
def getunits(Year,var="time",element="gph"):
    with netCDF4.Dataset(filename(Year,element), 'r',  format='NETCDF4_CLASSIC') as fileptr : 
        units = fileptr.variables[datavar[var]].units
    return(units)


