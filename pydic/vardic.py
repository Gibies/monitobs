#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 23 12:56:36 2019

@author: gibies
"""

import iris
import numpy

varlist=["sst", "sss", "ssh", "mld", "temp", "pres", "airt850", "airt500", "airt200"]
umvarlist=[ "pres", "airt", "uwnd", "vwnd", "wwnd", "sphum", "clud", ]
stdlev=["200", "500", "850"]
lev_list=[1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200]

plotname = { 
	"sst" : "sst",
	"temp"	: "airt",
	"airt850" : "airt",
	}

filename = { 
	"sst" : "foamlite.grid_T.nc",
	"airt"	: "data.nc",
	"airt850" : "data_850.nc",
	}

fvarname = { 
	"sst" : "sosstsst",
	"temp"	: "air_potential_temperature",
	"mdl"	: "model_level_number",
	"pres" 	: "air_pressure",
	"uwnd"	: "x_wind",
	"vwnd"	: "y_wind",
	"wwnd"	: "upward_air_velocity",
	"airt"	: "air_potential_temperature",
	"clud"	: "cloud_volume_fraction_in_atmosphere_layer",
	"sphum"	: "specific_humidity",
	"airt850"	: "air_potential_temperature",
	"airt500"	: "air_potential_temperature",
	"airt200"	: "air_potential_temperature",
	}

latvar = {
	"sst" : "nav_lat",
	"airt" : "latitude",
	"wind" : "latitude",
	"uwnd" : "latitude",
	"vwnd" : "latitude",
	}

lonvar = {
	"sst" : "nav_lon",
	"airt" : "longitude",
	"uwnd" : "longitude",
	"vwnd" : "longitude",
	"wind" : "longitude",
	}

dimopt = {
	"sst" : 2,
	"airt" : 1,
	"uwnd" : 1,
	"vwnd" : 1,
	"wind" : 1,
	}

iris_select={
"temp"	: iris.Constraint('air_potential_temperature'),
"pres"	: iris.Constraint('air_pressure'),
"uwnd"	: iris.Constraint('x_wind'),
"vwnd"	: iris.Constraint('y_wind'),
"wwnd"	: iris.Constraint('upward_air_velocity'),
}

cnlev = {
	"sst" : numpy.arange(0.,32.,2.),
	"airt850" : numpy.arange(280.,300.,1.),
	"uwnd" : numpy.arange(-10.,10.,2.),
	}

colrmin = {
	"sst" : 30,
	"airt850" : 30,
	}

colrmax = {
	"sst" : 170,
	"airt850" : 170,
	}

draworder = {
	"sst" : "Postdraw",
	"airt850" : "Predraw",
	}

writenc = {
	"sst" : False,
	"airt850" : True,
	}

sss_dic = {
	"shrtnam" : "sss",
	"filename" : "foamlite.grid_T.nc",
	"fvarname" : "sosaline",
	"latvar" : "nav_lat",
	"lonvar" : "nav_lon",
	"time_index" : 1,
	"dimopt" : 2,
	"cnlev" : numpy.arange(20.,40.,2.),
	"colrmin" : 30,
	"colrmax" : 170,
	"colrindx" :numpy.arange(30,170,1),
	"draworder" : "Postdraw",
	"writenc" : False,
	"plotname" : "sss",
	}

mld_dic = {
	"shrtnam" : "mld",
	"filename" : "foamlite.grid_T.nc",
	"fvarname" : "sokaraml",
	"latvar" : "nav_lat",
	"lonvar" : "nav_lon",
	"time_index" : 1,
	"dimopt" : 2,
	"cnlev" : numpy.arange(0.,500.,50.),
	"colrmin" : 30,
	"colrmax" : 170,
	"colrindx" :numpy.arange(30,170,1),
	"draworder" : "Postdraw",
	"writenc" : False,
	"plotname" : "mld",
	}

ssh_dic = {
	"shrtnam" : "ssh",
	"filename" : "foamlite.grid_T.nc",
	"fvarname" : "sossheig",
	"latvar" : "nav_lat",
	"lonvar" : "nav_lon",
	"time_index" : 1,
	"dimopt" : 2,
	"cnlev" : numpy.arange(0.,1.,0.1),
	"colrmin" : 30,
	"colrmax" : 170,
	"colrindx" :numpy.arange(30,170,1),
	"draworder" : "Postdraw",
	"writenc" : False,
	"plotname" : "ssh",
	}



sst_dic = {
	"shrtnam" : "sst",
	"filename" : "foamlite.grid_T.nc",
	"fvarname" : "sosstsst",
	"latvar" : "nav_lat",
	"lonvar" : "nav_lon",
	"time_index" : 1,
	"dimopt" : 2,
	"cnlev" : numpy.arange(0.,32.,2.),
	"colrmin" : 30,
	"colrmax" : 170,
	"colrindx" :numpy.arange(30,170,1),
	"draworder" : "Postdraw",
	"writenc" : False,
	"plotname" : "sst",
	}

pres_dic = {
	"shrtnam" : "pres",
	"latvar" : "latitude",
	"lonvar" : "longitude",
	"time_index" : 0,
	"dimopt" : 1,
	"colrmin" : 30,
	"colrmax" : 170,
	"draworder" : "Predraw",
	"writenc" : False,
	}

airt_dic = {
	"shrtnam" : "airt",
	"filename" : filename["airt"],
	"fvarname" : fvarname["airt"],
	"latvar" : "latitude",
	"lonvar" : "longitude",
	"time_index" : 0,
	"dimopt" : 1,
	"cnlev" : numpy.arange(280.,305.,5.),
	"colrmin" : 30,
	"colrmax" : 170,
	"colrindx" : numpy.arange(30,170,1),
	"draworder" : "Predraw",
	"writenc" : False,
	"plotname" : "airt_850",
	}

sst_bias={
	"colrindx"	: numpy.array([172,179,186,193,200,205,226,233,240,247,254]),
    	"cnlev"	: numpy.linspace(-1,1,9),
	}

sss_bias={
	"colrindx"	: numpy.array([172,179,186,193,200,205,226,233,240,247,254]),
    	"cnlev"	: numpy.linspace(-0.1,0.1,5),
	}

ssh_bias={
	"colrindx"	: numpy.array([172,179,186,193,200,205,226,233,240,247,254]),
    	"cnlev"	: numpy.linspace(-0.1,0.1,9),
	}

mld_bias={
	"colrindx"	: numpy.array([172,179,186,193,200,205,226,233,240,247,254]),
    	"cnlev"	: numpy.linspace(-50,50,5),
	}

airt_bias={
	#"colrindx"	: numpy.array([172,179,186,193,200,205,226,233,240,247,254]),
	"colrindx"	: numpy.arange(172,254,1),
    	"cnlev"	: numpy.linspace(-0.8,0.8,4),
	}

getdic = {
	"sst" : sst_dic,
	"sst_bias" : sst_bias,
	"sss" : sss_dic,
	"sss_bias" : sss_bias,
	"ssh" : ssh_dic,
	"ssh_bias" : ssh_bias,
	"mld" : mld_dic,
	"mld_bias" : mld_bias,
	"airt" : airt_dic,
	"airt_bias" : airt_bias,
	"airt_850_bias" : airt_bias,
	"airt850" : airt_dic,
	"airt500" : airt_dic,
	"airt200" : airt_dic,
	}

diclist=["sst","sst_bias","sss","sss_bias","ssh","ssh_bias","mld","mld_bias","airt",]

def get_dic(dicname):
    if dicname in diclist:
	dic=getdic[dicname].copy()
    else:
 	dic={}
	dic.update({"shrtnam" : dicname})
	if dicname in ["wind",]:
	     dic.update({"fvarname":(fvarname["uwnd"],fvarname["vwnd"])})
	else:
	     dic.update({"fvarname":fvarname[dicname]})
	dic.update({"latvar":latvar[dicname],"lonvar":lonvar[dicname]})
	dic.update({"dimopt":dimopt[dicname]})
    return(dic)

def get_dicval(dicname,valname):
	dic=getdic[dicname].copy()
	val=dic[valname]
	return(val)
