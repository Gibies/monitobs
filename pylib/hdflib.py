import os,sys
USER=os.environ.get('USER',"myhome")
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
PKGHOME=os.path.dirname(CURR_PATH)
OBSLIB=os.environ.get('OBSLIB',PKGHOME+"/pylib")
sys.path.append(OBSLIB)
OBSDIC=os.environ.get('OBSDIC',PKGHOME+"/pydic")
sys.path.append(OBSDIC)
OBSNML=os.environ.get('OBSNML',PKGHOME+"/nml")
sys.path.append(OBSNML)
print(OBSLIB)
import obslib
import pandas
import numpy
import math
import subprocess
from collections import OrderedDict
import h5py
import re

def get_grp_lst(filename):
    with h5py.File(filename, "r") as f:
	grplst=list(f.keys())
	return(grplst)

def get_var_lst(filename,grpindx=0,grpnam=None):
    with h5py.File(filename, "r") as f:
	if grpnam is None: grpnam = list(f.keys())[grpindx]
        data = list(f[grpnam])
        return(data) #List all the Datasets

def get_grp_attrs(filename,grpnam):
    with h5py.File(filename, "r") as f:
        attrs = dict(f['/'+grpnam].attrs.items())
	return(attrs)

def get_grp_attrs_val(filename,grpnam,key):
    attrs=get_grp_attrs(filename,grpnam)
    if key not in attrs:
	key=key.replace("_"," ")
    if key not in attrs:
	key=key.title()
    if key not in attrs:
	val=None
    else:
	val=attrs[key]
    return(val)

def scale_data(data,scale,formula):
	if scale is not None:
		scale=float(scale) 
	else:
		formula=""
	if "Scale * Value" in formula:
		data = numpy.multiply(data,scale)
	else:
		print("Formula not handled")
	return(data)

def shape_data(data1,shape1=None):
	if shape1 is None: 
		shape1=data1.shape
	dimcnt=len(shape1)
	if dimcnt is 2:
		tracklen=shape1[0]
		scanwdth=shape1[1]
	shpold=data1.shape
	if shpold[1] == tracklen: 
		data1=numpy.swapaxes(data1,0,1)
	if shpold[0] == 1 : 
		data=numpy.tile(data1,scanwdth)
	else:
		data=data1
	return(data)

def get_var_data(filename,grpnam,varnam,dims=None):
    with h5py.File(filename, "r") as f:
        dataptr = f.get(grpnam+"/"+varnam)
        data1 = numpy.array(dataptr)
    	#data1 = numpy.ma.masked_object(data1,int(32767))
    	#data1 = numpy.ma.masked_object(data1,int(65535))
        if str(data1.dtype) == "int16": 
		data1 = obslib.mask_array(data1, 32767)
	else: 
		print(data1.dtype)
        if str(data1.dtype) == "uint16": 
		data1 = obslib.mask_array(data1, 65535)
	else: 
		print(data1.dtype)
	scale=get_grp_attrs_val(filename,grpnam,varnam+" Scale")
	formula=get_grp_attrs_val(filename,grpnam,"Formula to derive value of a Parameter")
	data1=scale_data(data1,scale,formula)
	data=shape_data(data1,dims)
	if dims[0] == data.shape[0] and dims[1] == data.shape[1]: 
		data=data.flatten()
		print(data)
	else:
		print("Shape missmatch",dims,data.shape)
	return(data)

def frame_var_data(filename,grpnam,varnam,dims=(1720,144),data=None):
	data1=get_var_data(filename,grpnam,varnam,dims=dims)
	#print(data1)
	if data is None:
		data=pandas.DataFrame({varnam:data1})
	else:
		print(data1)
		data[varnam]=data1
	return(data)

def frame_data(filename,varnml,elist=None):
	varlstinfo=pandas.read_table(varnml)
	if elist is None: elist=varlstinfo.varname.values
	print(varlstinfo)
	data=pandas.DataFrame()
	for varnam in elist:
		grpnam=varlstinfo.query("varname == @varnam").grpname.values[0]
		print(filename,grpnam,varnam)
		data=frame_var_data(filename,grpnam,varnam,data=data)
	return(data)
