#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 12:30:50 2020

@author: gibies
"""


import numpy
import iris
import stratify
import sys,os,glob,datetime,netCDF4
iris.FUTURE.netcdf_no_unlimited=True

from functools import partial
interpolator = partial(stratify.interpolate,
                       interpolation=stratify.INTERPOLATE_NEAREST,
                       extrapolation=stratify.EXTRAPOLATE_LINEAR)

CURR_PATH=os.path.dirname(os.path.abspath(__file__))
PKGHOME=os.path.dirname(CURR_PATH)
OBSLIB=os.environ.get('OBSLIB',PKGHOME+"/pylib")
sys.path.append(OBSLIB)
OBSDIC=os.environ.get('OBSDIC',PKGHOME+"/pydic")
sys.path.append(OBSDIC)
OBSNML=os.environ.get('OBSNML',PKGHOME+"/nml")
sys.path.append(OBSNML)

import vardic

def hPa_to_units(lev_list,units):
    new_list=[]
    for lev in lev_list:
	if units == "Pa":
		lev=(float(lev)*100)
	if units == "hPa":
		lev=(float(lev))
	new_list.append(lev)
    return(new_list)

def extract_field(infile,outpath,varlist,outfile_prefix="data"):
	cubes=iris.load(infile)
	data=[]
	for var in varlist:
		select=vardic.iris_select[var]
		data.append(cubes.extract(select)[0])
	iris.save(data, outpath+'/'+outfile_prefix+'.nc');
	return(data)


def extract_level(infile,outpath,varlist,levpos,outfile_prefix="data"):
        cubes=iris.load(infile)
        data=[]
        for var in varlist:
                select=vardic.iris_select[var]
                data.append(cubes.extract(select)[0][levpos,:,:])
        iris.save(data, outpath+'/'+outfile_prefix+'.nc');
        return(data)

def extract_meridian(infile,outpath,varlist,levpos,latmin,latmax,outfile_prefix="data"):
        cubes=iris.load(infile)
        data=[]
        for var in varlist:
                select=vardic.iris_select[var]
                data.append(cubes.extract(select)[0][levpos,:,:].intersection(latitude=(latmin,latmax)))
        iris.save(data, outpath+'/'+outfile_prefix+'.nc');
        return(data)

def extract_tropics(infile,outpath,varlist,levpos,latminpos,latmaxpos,outfile_prefix="data"):
        cubes=iris.load(infile)
        data=[]
        for var in varlist:
                select=vardic.iris_select[var]
                data.append(cubes.extract(select)[0][levpos,latminpos:latmaxpos,:])
        iris.save(data, outpath+'/'+outfile_prefix+'.nc');
        return(data)

def extract_zone(infile,outpath,varlist,levpos,latminpos,latmaxpos,lonminpos,lonmaxpos,outfile_prefix="data"):
        cubes=iris.load(infile)
        data=[]
        for var in varlist:
                select=vardic.iris_select[var]
                data.append(cubes.extract(select)[0][levpos,latminpos:latmaxpos,lonminpos:lonmaxpos])
        iris.save(data, outpath+'/'+outfile_prefix+'.nc');
        return(data)

def get_coord_indx(cube,coord_name):
	co =  cube.dim_coords
    	axidx = 0
    	for i in range(len(co)):
        	if co[i].standard_name == coord_name:
            		axidx = i
            		break
	return(axidx)

def get_var_index(cubes,varname):
	indx=0
	for i in range(len(cubes)):
		if cubes[i].standard_name == varname:
			indx = i
			break
	return(indx)
	
def get_var_cube(cubes,varname):
	indx=get_var_index(cubes,varname)
	cube=cubes[indx]
	return(cube)

def get_var_data(cubes,varname):
	cube=get_var_cube(cubes,varname)
	vardata=cube.data
	return(vardata)

def get_var_units(cubes,varname):
	cube=get_var_cube(cubes,varname)
	units=cube.units
	return(units)

def get_alt(cube):
	orography = cube.coord('surface_altitude').points
	model_top = cube.coord('level_height').points[-1]
	model_level_len = len(cube.coord('level_height').points)
	sigma = 1.1 - numpy.logspace(-1, numpy.log10(1.1), model_level_len)
	sigma = sigma[:, numpy.newaxis, numpy.newaxis]
	altitude = (orography * sigma) + (model_top * (1 - sigma))
	return(altitude)

def isotopic_level_select(data,alti,lev):
	axidx=get_coord_indx(data[0],'model_level_number')
	isobar=stratify.interpolate(lev,pres,data,axis=axidx)
	return(isobar)

def regrid(prescube,varcube):
	scheme=iris.analysis.Nearest()
	pres_newgrid=prescube.regrid(varcube,scheme)
	return(pres_newgrid)

def frame_cube(data,attr_source):
	unit = attr_source.units
	sname = attr_source.standard_name
	lname = attr_source.long_name
	attr = attr_source.attributes
	cm = attr_source.cell_methods[0] if attr_source.cell_methods else None
	cube=iris.cube.Cube(data=data,units=unit, standard_name=sname,long_name=lname, attributes=attr)
	return(cube)

def copy_coords(cube,coord_source,indx_list=None):
	coords=coord_source.dim_coords
	if indx_list is None:
		indx_list=range(0,len(coords),1)
	for indx in indx_list:
		cube.add_dim_coord(coords[indx],indx)
	return(cube)

def dim_confirm(prescube,varcube):
	lev_name=vardic.fvarname['mdl']
	axidx=get_coord_indx(varcube,lev_name)
	print(prescube)
	pres_newgrid=regrid(prescube,varcube)
	coords=varcube.dim_coords
	print(coords)
	target_levs=numpy.array(coords[0].points)
	source_levs=numpy.array(pres_newgrid.dim_coords[0].points)
	print(source_levs)
	source_levs=source_levs[:, numpy.newaxis, numpy.newaxis]
	source_levs=source_levs*numpy.ones((1,len(coords[1].points),len(coords[2].points)),dtype=numpy.int)
	print(source_levs)
	pres_newlev=stratify.interpolate(target_levs,source_levs,pres_newgrid.data,axis=axidx)
	pres_newlev = numpy.ma.masked_invalid(pres_newlev)
	pres=frame_cube(pres_newlev,prescube)
	pres=copy_coords(pres,varcube)
	return(pres)

def extract_where(varcube,prescube,level,lev_name):
	varcube=regrid(varcube,prescube)
	axidx=get_coord_indx(varcube,lev_name)
	prescube=dim_confirm(prescube,varcube)
	newdata=stratify.interpolate(level.points,prescube.data[1:,:,:],varcube.data[1:,:,:],axis=axidx)
	newcube=frame_cube(newdata,varcube)
	newcube.add_dim_coord(level,0)
	newcube=copy_coords(newcube,varcube,[1,2])
	return(newcube)

def isobaric_level_select(cubes,lev_list,varlist=vardic.varlist):
	pres_name=vardic.fvarname['pres']
	prescube=get_var_cube(cubes,pres_name)
	pres_unit=prescube.units
	lev_name=vardic.fvarname['mdl']
	lev_list=hPa_to_units(lev_list,pres_unit)
	levels=numpy.array(lev_list)
 	print(levels)	
	level=iris.coords.DimCoord(levels, standard_name=pres_name, units=pres_unit, attributes={'comments': 'Isobaric levels'})
	newcubes=[]
	for var in varlist:
		if var != "pres":
			fvarname=vardic.fvarname[var]
			print(fvarname)
			varcube=get_var_cube(cubes,fvarname)
			isobar=extract_where(varcube,prescube,level,lev_name)
			newcubes.append(isobar)
	return(newcubes)

def horizontal_crop(cubes,latmin=-90,latmax=90,lonmin=0,lonmax=360,varlist=vardic.varlist):
	newcubes=[]
	for var in varlist:
		varname=vardic.fvarname[var]
		varcube=get_var_cube(cubes,varname)
		print(latmin,latmax)
		varcube=varcube.intersection(latitude=(latmin,latmax))
		varcube=varcube.intersection(longitude=(lonmin,lonmax))
		newcubes.append(varcube)
	return(newcubes)

def print_cubes(cubes):
	for indx in range(0,len(cubes),1):
		print(cubes[indx])


def write_netcdf(infile,outpath,outfile_prefix="data",latmin=-90,latmax=90,lonmin=0,lonmax=360,lev_list=vardic.lev_list,varlist=vardic.umvarlist):
	#if datadic is not None:
	#    diclist=["outfile_prefix","latmin", "latmax", "lonmin", "lonmax", "lev_list", "varlist"]
	#    for item in diclist:
	#	if item in datadic:
	#		outfile_prefix = datadic["outfile_prefix"]
	print(infile)
	cubes=iris.load(infile)
	print(cubes)
	#alt=get_alt(cubes[0])
	print(latmin,latmax)
	cubes=horizontal_crop(cubes,latmin=latmin,latmax=latmax,lonmin=lonmin,lonmax=lonmax,varlist=varlist)
	print(len(lev_list))
	newcubes=[]
	print(varlist)
	if len(lev_list) > 0:
	    levcubes=isobaric_level_select(cubes,lev_list,varlist)
	    print_cubes(levcubes)
	    iris.save(levcubes, outpath+'/'+outfile_prefix+'_stdlev.nc');
	else:
	    iris.save(cubes, outpath+'/'+outfile_prefix+'.nc');
	return(cubes)
		
	#altitude = get_alt(data[0])
	#data=data.extract(iris.Constraint(air_pressure=(85000,50000)))
	#axidx=get_coord_indx(data[0],'model_level_number')
	#target_altitudes= numpy.array([800,])
	#varindx=get_var_index(data,"air_pressure")
	#vardata=data[varindx].data
	#newCube = stratify.interpolate(target_altitudes,vardata,vardata,axis=axidx)
	#newCube = numpy.ma.masked_invalid(newCube)
	#print(newCube)
