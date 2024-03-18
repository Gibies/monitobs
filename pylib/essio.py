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

def iri_load_cubes(infile,cnst=None,callback=None,stashcode=None,option=0,dimlst=None):
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
    cubedimlst=[coord.name() for coord in cubes.dim_coords]
    cubeauxc=[coord.name() for coord in cubes.aux_coords]
    if dimlst is not None:
	for dimnam in dimlst:
	   if dimnam in cubeauxc:
		if len(cubes.coord(dimnam).points) is 1:
		   cubes=new_axis(cubes,dimnam)
    return(cubes)

def iri_to_nc(infile,varlst,outfile,callback=None,stashcode=None,option=2,dimlst=None,coords=None):
	cube=iri_load_cubes(infile,cnst=varlst,callback=callback,stashcode=stashcode,option=option,dimlst=dimlst)
	nc_file=iris.save(cube,outfile)
	return(nc_file)	

#############################################################################################################################
### IRIS and XARRAY combination based functions
#############################################################################################################################

def irx_cube_array(cube,varlst,dimlst=None,coords=None):
	cubedimlst=[coord.name() for coord in cube.dim_coords]
	cubeauxc=[coord.name() for coord in cube.aux_coords]
	if dimlst is None: dimlst=cubedimlst	
	datset=xarray.Dataset()
	if coords is None: 
		coords=datset.coords
		for dimnam in dimlst:
			coords.update({dimnam:cube.coord(dimnam).points,})
			unit=cube.coord(dimnam).units
			datset[dimnam].attrs['units'] = unit
	for var in varlst:
		data1=cube.data
		units=cube.units
		datset[var]=xarray.DataArray(data=data1,dims=dimlst,coords=coords,name=var)
		datset[var].attrs['units'] = units
	return(datset)

def irx_load_cubray(infile,varlst,callback=None,stashcode=None,option=2,dimlst=None,coords=None):
	cube=iri_load_cubes(infile,cnst=varlst,callback=callback,stashcode=stashcode,option=option,dimlst=dimlst)
	datset=irx_cube_array(cube,varlst,dimlst=dimlst,coords=coords)
	return(datset)


#############################################################################################################################
### IRIS, XARRAY and NIO combination based functions
#############################################################################################################################

def ixn_extract(infile,varlst,callback=None,stashcode=None,option=2,dimlst=None,coords=None,outfile=None,):
	datset=irx_load_cubray(infile,varlst,callback=callback,stashcode=stashcode,option=option,dimlst=dimlst,coords=coords)
	if varlst is None: varlst=[var.name for var in dataset.data_vars]
	var_lst_str=obslib.underscore(varlst)
	if outfile is None: outfile=infile.split(".")[0]+"_"+var_lst_str+".nc"
	if dimlst is None: dimlst=datset.dimlst
	if coords is None: coords=datset.coords
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
	varptr = fileptr.create_variable(varnam,vartyp,data.dimlst)
	fileptr.variables[varnam].assign_value(data)
	for attrnam in varattlst:
		fileptr=nix_write_varattr(datset,fileptr,varnam,attrnam)
	return(fileptr)

def nix_write(datset,filenam,dimlst=None,varlst=None):
	if dimlst is None: dimlst=xar_dimlst(datset)
	if varlst is None: varlst=xar_varlst(datset)
	fileptr=Nio.open_file(filenam, "rw")
	for dimnam in dimlst:
		dimptr=fileptr.create_dimension(dimnam,len(datset[dimnam]))
		fileptr=nix_write_var(datset,fileptr,dimnam,vartyp="d",varattlst=["units"])
	for varnam in varlst:
		fileptr=nix_write_var(datset,fileptr,varnam,vartyp="d",varattlst=["units"])
	fileptr.close()
	return(None)

def nix_read(filenam,dimlst,varlst):
	fileptr=Nio.open_file(filenam, "r")
	datset=xarray.Dataset()
	for varnam in varlst:
		datset=nix_read_var(fileptr,varnam,varattlst=["units"],datset=datset)
	for dimnam in dimlst:
		datset=nix_read_var(fileptr,dimnam,varattlst=["units"],datset=datset)
	return(datset)

def nix_read_var(fileptr,varnam,varattlst=None,datset=None):
	if varattlst is None: varattlst=["units"]
	if datset is None: datset=xarray.Dataset()
	var=fileptr.variables[varnam]
	type = var.typecode()
	numDims = var.rank
	dimSizes = var.shape
	dimlst = var.dimensions
	data=var.get_value()
	datset[varnam]=xarray.DataArray(data,name=varnam,dims=dimlst)
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

def xar_dimlst(datset):
	dimlst=datset.dims
	return(dimlst)

def xar_varlst(datset):
	varlst=[var for var in datset.data_vars]
	return(varlst)

def xar_print(datset,dimlst=None,varlst=None):
	if dimlst is None: dimlst=xar_dimlst(datset)
	if varlst is None: varlst=xar_varlst(datset)
	for dimnam in dimlst:
		print(datset.variables[dimnam])
		print(datset.variables[dimnam].attrs)
	for varnam in varlst:
		print(datset.variables[varnam].attrs)
	datset.close()
	return(None)

def xar_extract(filenam,dimlst=None,varlst=None):
	datset=xarray.open_dataset(filenam)
	if dimlst is None: dimlst=xar_dimlst(datset)
	if varlst is None: varlst=xar_varlst(datset)
	print("Dataset is loaded")
	xar_print(datset,dimlst,varlst)
	return(datset)

#############################################################################################################################
### Local functions
#############################################################################################################################

def datset_extract(infile,varlst,callback=None,stashcode=None,option=2,dimlst=None,coords=None,outfile=None,):
	switcher = {
		"0" :lambda: ixn_extract(infile,varlst,callback=callback,stashcode=stashcode,option=option,dimlst=dimlst,coords=coords),
		"1" :lambda: ixn_extract(infile,varlst,callback=callback,stashcode=stashcode,option=option,dimlst=dimlst,coords=coords),
		"2" :lambda: ixn_extract(infile,varlst,callback=callback,stashcode=stashcode,option=option,dimlst=dimlst,coords=coords),
		"3" :lambda: ixn_extract(infile,varlst,callback=callback,stashcode=stashcode,option=option,dimlst=dimlst,coords=coords),
		"4" :lambda: nix_read(infile,dimlst,varlst),
    	}
	func = switcher.get(opt, lambda: 'Invalid option')
	datset = func()
	dimlst=xar_dimlst(datset)
	varlst=xar_varlst(datset)
	#datset=xar_extract(outfile,dimlst,varlst)
	#xar_print(datset,dimlst,varlst)
	print(varlst)
	#print(varlst[0].name)
	void=nix_write(datset,outfile,dimlst,varlst)
	datset=nix_read(outfile,dimlst,varlst)
	#datset=xar_extract(outfile,dimlst,varlst)
	xar_print(datset,dimlst,varlst)
	return(datset)
