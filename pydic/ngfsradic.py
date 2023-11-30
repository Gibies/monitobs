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

ngfsrapath="/home/gibies/data/ngfsra"

def filename(Year,element):
    fileprefix="ngfs_reanl_"
    datafile = {
        "gph" : "%s/%s"%(ngfsrapath,fileprefix+"HR_HGT-prl_"+str(Year)+"120100-"+str(Year)+"123118.nc"),
        "tmp" : "%s/%s"%(ngfsrapath,fileprefix+"HR_TMP-prl_"+str(Year)+"120100-"+str(Year)+"123118.nc"),
        "shum" : "%s/%s"%(ngfsrapath,fileprefix+"HR_SPFH-2m-agl_"+str(Year)+"120100-"+str(Year)+"123118.nc"),
        "rhum" : "%s/%s"%(ngfsrapath,fileprefix+"HR_RH-prl_"+str(Year)+"120100-"+str(Year)+"123118.nc"),
        "uwnd" : "%s/%s"%(ngfsrapath,fileprefix+"HR_UGRD-prl_"+str(Year)+"120100-"+str(Year)+"123118.nc"),
        "vwnd" : "%s/%s"%(ngfsrapath,fileprefix+"HR_VGRD-prl_"+str(Year)+"120100-"+str(Year)+"123118.nc"),
        "tsfc" : "%s/%s"%(ngfsrapath,fileprefix+"HR_TMP-2m-agl_"+str(Year)+"120100-"+str(Year)+"123118.nc"),
        "rhsfc" : "%s/%s"%(ngfsrapath,fileprefix+"HR_RH-2m-agl_"+str(Year)+"120100-"+str(Year)+"123118.nc"),
        "usfc" : "%s/%s"%(ngfsrapath,fileprefix+"HR_UGRD-2m-agl_"+str(Year)+"120100-"+str(Year)+"123118.nc"),
        "vsfc" : "%s/%s"%(ngfsrapath,fileprefix+"HR_VGRD-2m-agl_"+str(Year)+"120100-"+str(Year)+"123118.nc"),
        "psfc" : "%s/%s"%(ngfsrapath,fileprefix+"HR_PRES-sfc_"+str(Year)+"120100-"+str(Year)+"123118.nc"),
        "slp" : "%s/%s"%(ngfsrapath,fileprefix+"HR_PRES-msl_"+str(Year)+"120100-"+str(Year)+"123118.nc"),
        "t2m" : "%s/%s"%(ngfsrapath,fileprefix+"HR_TMP-2m-agl_"+str(Year)+"120100-"+str(Year)+"123118.nc"),
        "rh2m" : "%s/%s"%(ngfsrapath,fileprefix+"HR_RH-2m-agl_"+str(Year)+"120100-"+str(Year)+"123118.nc"),
        "sh2m" : "%s/%s"%(ngfsrapath,fileprefix+"HR_SPFH-2m-agl_"+str(Year)+"120100-"+str(Year)+"123118.nc"),
        "u2m" : "%s/%s"%(ngfsrapath,fileprefix+"HR_GUST-sfc_"+str(Year)+"120100-"+str(Year)+"123118.nc"),
        "v2m" : "%s/%s"%(ngfsrapath,fileprefix+"HR_GUST-sfc_"+str(Year)+"120100-"+str(Year)+"123118.nc"),
        "u10m" : "%s/%s"%(ngfsrapath,fileprefix+"HR_UGRD-10m-agl_"+str(Year)+"120100-"+str(Year)+"123118.nc"),
        "v10m" : "%s/%s"%(ngfsrapath,fileprefix+"HR_VGRD-10m-agl_"+str(Year)+"120100-"+str(Year)+"123118.nc"),
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

###Usage filevar=filevartable(Year)
#def filevartable(Year):
#    datafile=filename(Year,"gph")
#    filevar=pd.DataFrame([[datafile["gph"],"time"]],index=["time"],columns=["file","fvname"])
#    filevar = filevar.append(pd.DataFrame([["%s/%s"%(ngfsrapath,"hgt."+str(Year)+".nc"),"level"]],index=["lev"],columns=["file","fvname"]))
#    filevar = filevar.append(pd.DataFrame([["%s/%s"%(ngfsrapath,"hgt."+str(Year)+".nc"),"lat"]],index=["lat"],columns=["file","fvname"]))
#    filevar = filevar.append(pd.DataFrame([["%s/%s"%(ngfsrapath,"hgt."+str(Year)+".nc"),"lon"]],index=["lon"],columns=["file","fvname"]))
#    filevar = filevar.append(pd.DataFrame([["%s/%s"%(ngfsrapath,"hgt."+str(Year)+".nc"),"hgt"]],index=["gph"],columns=["file","fvname"]))
#    filevar = filevar.append(pd.DataFrame([["%s/%s"%(ngfsrapath,"air."+str(Year)+".nc"),"air"]],index=["tmp"],columns=["file","fvname"]))
#    filevar = filevar.append(pd.DataFrame([["%s/%s"%(ngfsrapath,"rhum."+str(Year)+".nc"),"rhum"]],index=["rhum"],columns=["file","fvname"]))
#    filevar = filevar.append(pd.DataFrame([["%s/%s"%(ngfsrapath,"uwnd."+str(Year)+".nc"),"uwnd"]],index=["uwnd"],columns=["file","fvname"]))
#    filevar = filevar.append(pd.DataFrame([["%s/%s"%(ngfsrapath,"vwnd."+str(Year)+".nc"),"vwnd"]],index=["vwnd"],columns=["file","fvname"]))
#    return(filevar)
#    

def getdata(Year,var="time",element="gph"):
    print(var,element)
    with netCDF4.Dataset(filename(Year,element), 'r',  format='NETCDF4_CLASSIC') as fileptr : 
        data = fileptr.variables[datavar[var]][:]
    return(data)
    
def getunits(Year,var="time",element="gph"):
    with netCDF4.Dataset(filename(Year,element), 'r',  format='NETCDF4_CLASSIC') as fileptr : 
        units = fileptr.variables[datavar[var]].units
    return(units)

#gphdata = {
#            "gph" : getdata(DT,"gph"),
#            "time" : getdata(DT,"time"), # with netCDF4.Dataset(filevar.file["gph"], 'r',  format='NETCDF4_CLASSIC') as gphfile : gphfile.variables[filevar.fvname["time"]][:],
#            "time_units" : getunits(DT,"time"), # with netCDF4.Dataset(filevar.file["gph"], 'r',  format='NETCDF4_CLASSIC') as gphfile : gphfile.variables[filevar.fvname["time"]].units,
#            #"lev" : with netCDF4.Dataset(filevar.file["gph"], 'r',  format='NETCDF4_CLASSIC') as gphfile : gphfile.variables[filevar.fvname["lev"]][:],
#            #"lat" : with netCDF4.Dataset(filevar.file["gph"], 'r',  format='NETCDF4_CLASSIC') as gphfile : gphfile.variables[filevar.fvname["lat"]][:],
#            #"lon" : with netCDF4.Dataset(filevar.file["gph"], 'r',  format='NETCDF4_CLASSIC') as gphfile : gphfile.variables[filevar.fvname["lon"]][:]
#            }

