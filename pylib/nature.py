#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 17:22:59 2021

@author: gibies
"""
from __future__ import print_function
import sys
import os
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
PKGHOME=os.path.dirname(CURR_PATH)
OBSLIB=os.environ.get('OBSLIB',PKGHOME+"/pylib")
sys.path.append(OBSLIB)
OBSDIC=os.environ.get('OBSDIC',PKGHOME+"/pydic")
sys.path.append(OBSDIC)
OBSNML=os.environ.get('OBSNML',PKGHOME+"/nml")
sys.path.append(OBSNML)
obs_index_nml="obs_index_nml"
nmlfile="%s/%s" % (OBSNML,obs_index_nml)
import obslib
import obsdic
import netCDF4
import ngfsradic as datadic
#import ncepradic as datadic
#import imdaadic as datadic


def getdata(var="time",element=None,year=None,month=None,day=None):
    if element is None: element=var
    print(var,element)
    with netCDF4.Dataset(datadic.filename(element,year,month,day), 'r',  format='NETCDF4_CLASSIC') as fileptr : 
        data = fileptr.variables[datadic.datavar[var]][:]
    return(data)
    
def getunits(var="time",element=None,year=None,month=None,day=None):
    if element is None: element=var
    with netCDF4.Dataset(datadic.filename(element,year,month,day), 'r',  format='NETCDF4_CLASSIC') as fileptr : 
        units = fileptr.variables[datadic.datavar[var]].units
    return(units)

def getfiledim(Tnow,element=None):
    if element is None: element="tmp"
    Year=Tnow.year
    Month=Tnow.month
    Day=Tnow.day
    time=getdata("time",element,year=Year,month=Month,day=Day)
    time_units=getunits("time",element,year=Year,month=Month,day=Day)
    lat=getdata("lat",element,year=Year,month=Month,day=Day)
    lon=getdata("lon",element,year=Year,month=Month,day=Day)
    lev=getdata("lev",element,year=Year,month=Month,day=Day)
    filedim={
            "time":time,
            "time_units":time_units,
            "lat":lat,
            "lon":lon,
            "lev": lev,
            }
    return(Year,filedim)



#def getdata(Year,var="time",element=None,year=None,month=None,day=None):
#	if element is None: element=var
#	data=ngfsradic.getdata(Year,var,element)
#	data=ncepradic.getdata(Year,var,element)
#	data=imdaadic.getdata(Year,var,element)
#	return(data)
#
#def getunits(Year,var="time",element=None):
#    if element is None: element=var
#    units=ngfsradic.getunits(Year,var,element)
#    return(units)
