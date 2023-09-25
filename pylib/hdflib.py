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
import obslib
import pandas
import numpy
import math
import subprocess
from collections import OrderedDict
import h5py
import re

def hdf_keylst(filename,prefix=''):
	with h5py.File(filename, "r") as f:
	    if len(prefix) > 0:
		if isinstance(f[prefix], h5py.Dataset):
			grps=None
		else:
			grps=f[prefix].keys()
	    else:
		grps=f.keys()
	return(grps)

def hdf_path(filename,prefix='',grpnam=None,varnam=None):
	if grpnam in hdf_keylst(filename,prefix):
		fpath=prefix+'/'+grpnam
	else:
	    print("Group name : "+str(grpnam)+" is not in path : "+str(prefix))
	    fpath=prefix
	with h5py.File(filename, "r") as f:
	    if varnam is not None: 
		if varnam in f[fpath].keys():
			fpath=fpath+'/'+varnam
	    if fpath in f: return(fpath)

def hdf_items(filename,prefix='',grpnam=None):
	fpath=hdf_path(filename,prefix,grpnam)
	if not fpath is None : itms=hdf_keylst(filename,fpath)
	else: itms=None
	hdfinfo={
		"path":fpath,
		"items":itms,
		}
	return(hdfinfo)

def hdf_locate(Name,filename,prefix='',grps=None):
	filepath=None
	if isinstance(grps, str):
	    if grps in hdf_keylst(filename,prefix):
		fdic=hdf_items(filename,prefix,grps)   
		if not fdic["items"] is None: filepath=hdf_locate(Name,filename,fdic["path"],fdic["items"])
	if isinstance(grps, list):
	    if Name in grps:
		filepath=hdf_path(filename,prefix,Name)
	    else:
		for grp in grps:
		   if filepath is None:
			fdic=hdf_items(filename,prefix,grp)
			if not fdic["items"] is None: filepath=hdf_locate(Name,filename,fdic["path"],fdic["items"])
	if grps is None: 
		fdic={"path":'',"items":hdf_keylst(filename,prefix)}
		if not fdic["items"] is None: filepath=hdf_locate(Name,filename,fdic["path"],fdic["items"])
	if filepath is not None: return filepath

def hdf_get_data(filename,fpath):
    with h5py.File(filename, "r") as f:
        dataptr = f.get(fpath)
        data1 = numpy.array(dataptr)
    return(data1)

def hdf_var_dims(filename,varname):
	loc=hdf_locate(varname,filename)
	data=hdf_get_data(filename,loc)
	dims=data.shape
	return(dims)

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
		formula=""
	if formula is not "":
		print("Warnning: Formula not handled")
		print("Formula: "+formula)
	return(data)

def shape_data(data1,shape1=None):
	if shape1 is None: 
		shape1=data1.shape
	dimcnt=len(shape1)
	if dimcnt is 2:
		tracklen=shape1[0]
		scanwdth=shape1[1]
	shpold=data1.shape
	if len(shpold) == 1:
		data1 = numpy.expand_dims(data1,axis=1)
	shpold=data1.shape
	if shpold[1] == tracklen: 
		data1=numpy.swapaxes(data1,0,1)
	shpold=data1.shape
	if shpold[1] == 1 : 
		data1=numpy.tile(data1,scanwdth)
	data=data1
	return(data)

def get_var_data(filename,grpnam,varnam,dims=None,prefix=''):
    fpath=hdf_locate(varnam,filename,prefix,grpnam)
    print(fpath)
    with h5py.File(filename, "r") as f:
        dataptr = f.get(fpath)
        data1 = numpy.array(dataptr)
        if str(data1.dtype) == "int16": 
		data1 = obslib.mask_array(data1)	#, 32767)
	else: 
		print(data1.dtype)
        if str(data1.dtype) == "uint16": 
		data1 = obslib.mask_array(data1)	#, 65535)
	else: 
		print(data1.dtype)
	scale=get_grp_attrs_val(filename,grpnam,varnam+" Scale")
	formula=get_grp_attrs_val(filename,grpnam,"Formula to derive value of a Parameter")
	data1=scale_data(data1,scale,formula)
	data=shape_data(data1,dims)
	shpdata=data.shape
	if len(shpdata) == 3 and dims[0] == shpdata[0] and dims[1] == shpdata[1]:
	   data=data.reshape((shpdata[0]*shpdata[1],shpdata[2]))
	else:
		if len(shpdata) == 2 and dims[0] == shpdata[0] and dims[1] == shpdata[1]: 
			data=data.flatten()
		else:
			print("Shape missmatch",dims,shpdata)
			print("Warnning: "+str(shpdata))
	print(data.shape)
	return(data)

def frame_const_data(const,elenam,dims,data=None,fillval=numpy.nan):
	if data is None: data=pandas.DataFrame()
	if elenam in const:
		data1=[const[elenam]]*(dims[0]*dims[1])
	else :
		data1=[fillval]*(dims[0]*dims[1])
	data[elenam]=data1
	return(data)

def frame_var_data(filename,grpnam,varnam,elenam,dims,data=None,const={},fillval=numpy.nan):
	if grpnam in hdf_keylst(filename): 
		data1=get_var_data(filename,grpnam,varnam,dims=dims)
		#print(data1)
		print(data1.shape)
		if data is None: data=pandas.DataFrame()
		if len(data1.shape) == 1 : data[elenam]=data1
		if len(data1.shape) == 2 :
			chnls = range(1,(data1.shape[1]+1),1)
			print(chnls)
			for indx in chnls:
				data[elenam+"_"+str(indx)]=data1[:,(indx-1)]
	else : 
		data=frame_const_data(const,elenam,dims,data,fillval)
	return(data)

def frame_data(filename,varnml,elist=None,dtfmt='%Y-%jT%H:%M:%S.%f',const={},dims=None,refvar="Latitude",fillval=numpy.nan):
	filename=obslib.globlist(filename)[0]
	if dims is None: dims=hdf_var_dims(filename,refvar)
	varlstinfo=pandas.read_table(varnml)
	if elist is None: elist=varlstinfo.indx.values
	print(varlstinfo)
	data=pandas.DataFrame()
	for indx in elist:
		grpnam=varlstinfo.query("indx == @indx").grpname.values[0]
		varnam=varlstinfo.query("indx == @indx").varname.values[0]
		elenam=varlstinfo.query("indx == @indx").elename.values[0]
		print(filename,grpnam,varnam,elenam)
		data=frame_var_data(filename,grpnam,varnam,elenam,dims,data,const,fillval)
	data=obslib.pandas_dtfmt(data,dtfmt)
	return(data)

def datetime(filename,varnml,dtfmt='%Y-%jT%H:%M:%S.%f',indx=0):
	filename=obslib.globlist(filename)[0]
	varlstinfo=pandas.read_table(varnml)
	elist=[3,4,5,6,7,8]
	data=pandas.DataFrame()
	grpnam=varlstinfo.query("indx == @elist[0]").grpname.values[0]
	varnam=varlstinfo.query("indx == @elist[0]").varname.values[0]
	fpath=hdf_locate(varnam,filename)
    	with h5py.File(filename, "r") as f:
        	dataptr = f.get(fpath)
        	data1 = numpy.array(dataptr).flatten()[indx]
	datetime=obslib.datetimeframe(data1,dtfmt)
	return(datetime)
