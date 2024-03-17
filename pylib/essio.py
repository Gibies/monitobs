#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 12:56:36 2019

@author: gibies
"""
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
import xarray
import iris
from iris.util import new_axis


#############################################################################################################################
### 
#############################################################################################################################





















#############################################################################################################################
### IRIS based functions
#############################################################################################################################

def iri_load_cubes(infile,cnst=None,callback=None,stashcode=None,option=0,dims=None):
    opt=str(option)
    if stashcode is not None: cnst=iris.AttributeConstraint(STASH=stashcode)
    switcher = {
       "0" :lambda: iris.load(infile,constraints=cnst, callback=callback),
       "1" :lambda: iris.load_cubes(infile,constraints=cnst, callback=callback),
       "2" :lambda: iris.load_cube(infile,constraint=cnst, callback=callback),
       "3" :lambda: iris.load_raw(infile,constraints=cnst, callback=callback),
    }
    func = switcher.get(opt, lambda: 'Invalid option')
    cubes = func()
    cubedims=[coord.name() for coord in cubes.dim_coords]
    cubeauxc=[coord.name() for coord in cubes.aux_coords]
    if dims is not None:
	for dimnam in dims:
	   if dimnam in cubeauxc:
		if len(cubes.coord(dimnam).points) is 1:
		   cubes=new_axis(cubes,dimnam)
    return(cubes)

def iri_to_nc(infile,varnames,outfile,callback=None,stashcode=None,option=2,dims=None,coords=None):
	cube=iri_load_cubes(infile,cnst=varnames,callback=callback,stashcode=stashcode,option=option,dims=dims)
	nc_file=iris.save(cube,outfile)
	return(nc_file)	

#############################################################################################################################
### IRIS and XARRAY combination based functions
#############################################################################################################################

def irx_cube_array(cube,varnames,dims=None,coords=None):
	cubedims=[coord.name() for coord in cube.dim_coords]
	cubeauxc=[coord.name() for coord in cube.aux_coords]
	if dims is None: dims=cubedims	
	datset=xarray.Dataset()
	if coords is None: 
		coords=datset.coords
		for dimnam in dims:
			coords.update({dimnam:cube.coord(dimnam).points,})
			unit=cube.coord(dimnam).units
			datset[dimnam].attrs['units'] = unit
	for var in varnames:
		data1=cube.data
		units=cube.units
		datset[var]=xarray.DataArray(data=data1,dims=dims,coords=coords,name=var)
		datset[var].attrs['units'] = units
	return(datset)

def irx_load_cubray(infile,varnames,callback=None,stashcode=None,option=2,dims=None,coords=None):
	cube=iri_load_cubes(infile,cnst=varnames,callback=callback,stashcode=stashcode,option=option,dims=dims)
	datset=irx_cube_array(cube,varnames,dims=dims,coords=coords)
	return(datset)


#############################################################################################################################
### IRIS, XARRAY and NIO combination based functions
#############################################################################################################################

def ixn_extract(infile,varnames,callback=None,stashcode=None,option=2,dims=None,coords=None,outfile=None,):
	datset=irx_load_cubray(infile,varnames,callback=callback,stashcode=stashcode,option=option,dims=dims,coords=coords)
	var_lst_str=obslib.underscore(varnames)
	if outfile is None: outfile=infile.split(".")[0]+"_"+var_lst_str+".nc"
	if dims is None: dims=datset.dims
	if coords is None: coords=datset.coords
	#void=nix_write(datset,outfile,dims,varnames)
	datset_new=nix_read(outfile,dims,varnames)
	#datset_new=xar_extract(outfile,dims,varnames)
	return(None)

#############################################################################################################################
### NCAR Input Output Library (NIO) and XARRAY combination based functions
#############################################################################################################################


def nix_write_varattr(datset,fileptr,varnam,attrnam,attrtyp="str"):
	attrval=datset[varnam].attrs[attrnam]
	if attrtyp is "str": attrval=str(attrval)
	attrptr=setattr(fileptr.variables[varnam],attrnam,attrval)
	return(fileptr)

def nix_write_var(datset,fileptr,varnam,vartyp="d",varattlst=None):
	if varattlst is None: varattlst=["units"]
	data=datset[varnam]
	varptr = fileptr.create_variable(varnam,vartyp,data.dims)
	fileptr.variables[varnam].assign_value(data)
	for attrnam in varattlst:
		fileptr=nix_write_varattr(datset,fileptr,varnam,attrnam)
	return(fileptr)

def nix_write(datset,filenam,dimlist,varlist):
	fileptr=Nio.open_file(filenam, "rw")
	for dimnam in dimlist:
		dimptr=fileptr.create_dimension(dimnam,len(datset[dimnam]))
		fileptr=nix_write_var(datset,fileptr,dimnam,vartyp="d",varattlst=["units"])
	for varnam in varlist:
		fileptr=nix_write_var(datset,fileptr,varnam,vartyp="d",varattlst=["units"])
	fileptr.close()
	return(None)

def nix_read(filenam,dimlist,varlist):
	fileptr=Nio.open_file(filenam, "r")
	datset=xarray.Dataset()
	for varnam in varlist:
		datset=nix_read_var(fileptr,varnam,varattlst=["units"],datset=datset)
	for dimnam in dimlist:
		datset=nix_read_var(fileptr,dimnam,varattlst=["units"],datset=datset)
	xar_print(datset,dimlist,varlist)
	return(datset)

def nix_read_var(fileptr,varnam,varattlst=None,datset=None):
	if varattlst is None: varattlst=["units"]
	if datset is None: datset=xarray.Dataset()
	var=fileptr.variables[varnam]
	type = var.typecode()
	numDims = var.rank
	dimSizes = var.shape
	dims = var.dimensions
	data=var.get_value()
	datset[varnam]=xarray.DataArray(data,name=varnam,dims=dims)
	for attrnam in varattlst:
		datset=nix_read_varattr(fileptr,varnam,attrnam,datset=datset)
	return(datset)

def nix_read_varattr(fileptr,varnam,attrnam,datset=None):
	if datset is None: datset=xarray.Dataset()
	attrval=fileptr.variables[varnam].attributes[attrnam]
	datset.variables[varnam].attrs[attrnam]=attrval
	return(datset)

#############################################################################################################################
### XARRAY based functions
#############################################################################################################################

def xar_print(datset,dimlist,varlist):
	for dimnam in dimlist:
		print(datset.variables[dimnam])
		print(datset.variables[dimnam].attrs)
	for varnam in varlist:
		print(datset.variables[varnam].attrs)
	datset.close()
	return(None)

def xar_extract(filenam,dimlist,varlist):
	datset=xarray.open_dataset(filenam)
	print("Dataset is loaded")
	xar_print(datset,dimlist,varlist)
	return(datset)

#############################################################################################################################
### Local functions
#############################################################################################################################

def datset_extract(infile,varnames,callback=None,stashcode=None,option=2,dims=None,coords=None,outfile=None,):
	datset=ixn_extract(infile,varnames,callback=callback,stashcode=stashcode,option=option,dims=dims,coords=coords,outfile=outfile)
	return(datset)
