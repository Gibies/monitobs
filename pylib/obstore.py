#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 15:45:39 2018

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
import obslib
#import obsdic
import obsheader
import fixheader
import numpy
import pandas
import struct
import datetime
import itertools
diaglev=int(os.environ.get('GEN_MODE',0))
MAXINDX=int(os.environ.get('MAXINDX',fixheader.MAXINDX))
HDRSIZE=int(os.environ.get('HDRSIZE',fixheader.HDRSIZE))
LUTSIZE=int(os.environ.get('LUTSIZE',fixheader.LUTSIZE))
HBpos=int(os.environ.get('HBpos',fixheader.HBpos))
HBlen=int(os.environ.get('HBlen',fixheader.HBlen))
HCpos=int(os.environ.get('HCpos',fixheader.HCpos))
HClen=int(os.environ.get('HClen',fixheader.HClen))
HDR20=fixheader.HDR20
HDRgam=fixheader.HDRgam
FIXHDR=fixheader.FIXHDR
TREF=19700101
HLFTW=3
NAN_VAL_INT=-32768
NAN_VAL=-1.07374182e+09
NML_OBS_INDX="%s/%s" % (OBSNML,"obs_index_nml")

def obstore_info_fix(obstore_info):
	if "timewindowhalf" not in obstore_info: obstore_info["timewindowhalf"]=HLFTW
	return(obstore_info)

def errprint(*args, **kwargs):
    if diaglev > 0: print(*args, file=sys.stderr, **kwargs)

def binfmt(sec_nam,size):
    return {
        'fix1': ">"+str(size)+"q",
        'fix2': ">"+str(size)+"d",
        'hdralp': ">"+str(size)+"d",
        'hdrbet': ">"+str(size)+"d",
        'hdrgam': ">"+str(size)+"d",
        'alpha': ">"+str(size)+"q",
        'beeta': ">"+str(size)+"q",
        'gamma': ">"+str(size)+"d",
        'ldc' : ">"+str(size)+"d",
        'rdc' : ">"+str(size)+"d",
        'cdc' : ">"+str(size)+"d",
        'lut' : ">"+str(size)+"q",
        'data': ">"+str(size)+"d"
    }.get(sec_nam, ">"+str(size)+"d")

def binfill(sec_nam,size,fillval=None,fillvalint=None):
    if fillval is None: fillval=NAN_VAL
    if fillvalint is None: fillvalint=NAN_VAL_INT
    #fillval=numpy.nan
    return {
	'hdr':	[fillval]*size,
        'fix1':	[fillvalint]*size,
        'fix2': [fillval]*size,
        'hdralp':[fillval]*size,
        'hdrbet':[fillval]*size,
        'hdrgam':[fillval]*size,
        'alpha':[fillvalint]*size,
        'beeta':[fillvalint]*size,
        'gamma':[fillval]*size,
        'ldc' :	[fillval]*size,
        'rdc' :	[fillval]*size,
        'cdc' :	[fillval]*size,
        'lut' :	[fillvalint]*size,
        'data':	[fillval]*size
    }.get(sec_nam, [fillval]*size)

def binpos(obsfile,sec_nam):
    return {
	'hdr': (1,340,1),
        'fix1': (1,305,1),
        'fix2': obstore_get_pos(obsfile,105,106), #(306,34,1),
        'hdralp': (1,256,1),
        'hdrbet': obstore_get_pos(obsfile,100,101,fmtkey="d"), #(257,49,1),
        'hdrgam': obstore_get_pos(obsfile,105,106,fmtkey="d"), #(306,34,1),
        'alpha': (1,256,1),
        'beeta': obstore_get_pos(obsfile,100,101), #(257,49,1),
        'gamma': obstore_get_pos(obsfile,105,106), #(306,34,1),
        'ldc' : obstore_get_pos(obsfile,110,111,112),
        'rdc' : obstore_get_pos(obsfile,115,116,117),
        'cdc' : obstore_get_pos(obsfile,120,121,122),
        'lut' : obstore_get_pos(obsfile,150,151,152),
        'data': obstore_get_pos(obsfile,160,161,162)
    }.get(sec_nam, (1,340,1))
    
def obstore_read_header(obsfile,pos,size,fmtkey="q"):
    val=obslib.binary_read(obsfile,pos,size,fmtkey)
    return(val)
    #obsfile.seek((pos-1)*8,0)
    #data = obsfile.read(size*8)
    #val = struct.unpack(">"+str(size)+"q", data)
    #if size is 1:
    #    return(val)[0]
    #else:
    #    return(val)

def bin_write(obsfile,pos,data_fmt,data):
    #print(data,data_fmt,pos)
    #obstore_overwrite_formatted(data,data_fmt,obsfile,pos)
    obslib.binary_write(data,pos,obsfile)
    #####Different syntax####
    #obsfile.seek((pos-1)*8,0)
    #frmtstr=struct.Struct(data_fmt)
    #obsfile.write(frmtstr.pack(data))

def obstore_write_int8byte(data,obsfile):
    frmtstr=struct.Struct(">%sq"%(len(data)))
    obsfile.write(frmtstr.pack(*data))
    #struct.pack(">%sq"%(len(data)),*data)

def obstore_read_real(obsfile,pos,size,fmtkey="d"):
    val=obslib.binary_read(obsfile,pos,size,fmtkey)
    return(val)
    #obsfile.seek((pos-1)*8,0)
    #data = obsfile.read(size*8)
    #val = struct.unpack(">"+str(size)+"d", data)
    #if size is 1:
    #    return(val)[0]
    #else:
    #    return(val)

def obstore_write_dec8byte(data,obsfile):
    frmtstr=struct.Struct(">%sd"%(len(data)))
    obsfile.write(frmtstr.pack(*data))
    #struct.pack(">%sq"%(len(data)),*data) 

def obstore_overwrite_dec8byte(data,obsfile,pos):
    obsfile.seek((pos-1)*8,0)
    frmtstr=struct.Struct(">%sd"%(len(data)))
    obsfile.write(frmtstr.pack(*data))

def obstore_read_formatted(obsfile,pos,size,data_fmt):
    fmtkey=data_fmt[-1]
    val=obslib.binary_read(obsfile,pos,size,fmtkey)
    return(val)
    #obsfile.seek((pos-1)*8,0)
    #data = obsfile.read(size*8)
    #val = struct.unpack(data_fmt, data)
    #if size is 1:
    #    return(val)[0]
    #else:
    #    return(val)

def obstore_write_formatted(data,data_fmt,obsfile):
    frmtstr=struct.Struct(data_fmt)
    obsfile.write(frmtstr.pack(*data))

def obstore_overwrite_formatted(data,data_fmt,obsfile,pos):
    obsfile.seek((pos-1)*8,0)
    frmtstr=struct.Struct(data_fmt)
    if type(data) is int:
        obsfile.write(frmtstr.pack(data))
    elif type(data) is float:
        obsfile.write(frmtstr.pack(data))
    else:
        obsfile.write(frmtstr.pack(*data))

def obstore_write_toheaderpos(obsfile,pos,size,data):
    data_fmt=">"+str(size)+"q"
    obstore_overwrite_formatted(data,data_fmt,obsfile,pos)

def obstore_write_obsgroup(obsgroup,obsfile):
    obstore_write_toheaderpos(obsfile,287,1,obsgroup)
    #obstore_write_toheaderpos(obsfile,287,1,obsgroup)

def obstore_write_batchcount(batchcount,obsfile):
    obstore_write_toheaderpos(obsfile,112,1,batchcount)
    obstore_write_toheaderpos(obsfile,117,1,batchcount)
    obstore_write_toheaderpos(obsfile,122,1,batchcount)
    obstore_write_toheaderpos(obsfile,152,1,batchcount)
    obstore_write_toheaderpos(obsfile,288,1,batchcount)
    obstore_write_toheaderpos(obsfile,305,1,batchcount)

def obstore_write_datapos(datapos,obsfile):
    obstore_write_toheaderpos(obsfile,160,1,datapos)

def obstore_write_datalen(datalen,obsfile):
    obstore_write_toheaderpos(obsfile,161,1,datalen)
    obstore_write_toheaderpos(obsfile,284,1,datalen)

def obstore_write_stime(stime,obsfile):
    obstore_write_toheaderpos(obsfile,270,1,stime)
    obstore_write_toheaderpos(obsfile,286,1,stime)

def obstore_write_etime(etime,obsfile):
    obstore_write_toheaderpos(obsfile,271,1,etime)

def obstore_write_year(year,obsfile,nmlfile):
    obstore_overwrite_formatted(year,">1q",obsfile,21)
    obstore_overwrite_formatted(year,">1q",obsfile,28)
    obstore_overwrite_formatted(year,">1q",obsfile,35)
    nsubtype=obstore_read_header(obsfile,152,1)
    for irow,subtype in enumerate(obstore_read_subtype(obsfile), start=1):
        obstore_write_subhead_segment(obsfile,"lut",irow,1,1,int(year))
        obstore_write_subhead_segment(obsfile,"lut",irow,7,1,int(year))
        #obstore_write_data_element(obsfile,nmlfile,irow,"Year",year)

def obstore_write_month(month,obsfile,nmlfile):
    obstore_overwrite_formatted(month,">1q",obsfile,22)
    obstore_overwrite_formatted(month,">1q",obsfile,29)
    obstore_overwrite_formatted(month,">1q",obsfile,36)
    nsubtype=obstore_read_header(obsfile,152,1)
    for irow,subtype in enumerate(obstore_read_subtype(obsfile), start=1):
        obstore_write_subhead_segment(obsfile,"lut",irow,2,1,int(month))
        obstore_write_subhead_segment(obsfile,"lut",irow,8,1,int(month))
        #obstore_write_data_element(obsfile,nmlfile,irow,"Month",month)

def obstore_write_day(day,obsfile,nmlfile):
    obstore_overwrite_formatted(day,">1q",obsfile,23)
    obstore_overwrite_formatted(day,">1q",obsfile,30)
    obstore_overwrite_formatted(day,">1q",obsfile,37)
    nsubtype=obstore_read_header(obsfile,152,1)
    for irow,subtype in enumerate(obstore_read_subtype(obsfile), start=1):
        obstore_write_subhead_segment(obsfile,"lut",irow,3,1,int(day))
        obstore_write_subhead_segment(obsfile,"lut",irow,9,1,int(day))
        #obstore_write_data_element(obsfile,nmlfile,irow,"Day",day)

def obstore_write_hour(hour,obsfile,nmlfile):
    obstore_overwrite_formatted(hour,">1q",obsfile,24)
    obstore_overwrite_formatted(hour,">1q",obsfile,31)
    obstore_overwrite_formatted(hour,">1q",obsfile,38)
    nsubtype=obstore_read_header(obsfile,152,1)
    for irow,subtype in enumerate(obstore_read_subtype(obsfile), start=1):
        obstore_write_subhead_segment(obsfile,"lut",irow,4,1,int(hour))
        obstore_write_subhead_segment(obsfile,"lut",irow,10,1,int(hour))
        #obstore_write_data_element(obsfile,nmlfile,irow,"Hour",hour)

def obstore_write_minute(minute,obsfile,nmlfile):
    obstore_overwrite_formatted(minute,">1q",obsfile,25)
    obstore_overwrite_formatted(minute,">1q",obsfile,32)
    obstore_overwrite_formatted(minute,">1q",obsfile,39)
    nsubtype=obstore_read_header(obsfile,152,1)
    for irow,subtype in enumerate(obstore_read_subtype(obsfile), start=1):
        obstore_write_subhead_segment(obsfile,"lut",irow,5,1,int(minute))
        obstore_write_subhead_segment(obsfile,"lut",irow,11,1,int(minute))
        #obstore_write_data_element(obsfile,nmlfile,irow,"Minute",minute)

def obstore_write_second(second,obsfile,nmlfile):
    obstore_overwrite_formatted(second,">1q",obsfile,26)
    obstore_overwrite_formatted(second,">1q",obsfile,33)
    obstore_overwrite_formatted(second,">1q",obsfile,40)

def obstore_write_julday(julday,obsfile,nmlfile):
    obstore_overwrite_formatted(julday,">1q",obsfile,27)
    obstore_overwrite_formatted(julday,">1q",obsfile,34)
    obstore_overwrite_formatted(julday,">1q",obsfile,41)
    nsubtype=obstore_read_header(obsfile,152,1)
    for irow,subtype in enumerate(obstore_read_subtype(obsfile), start=1):
        obstore_write_subhead_segment(obsfile,"lut",irow,6,1,int(julday))
        obstore_write_subhead_segment(obsfile,"lut",irow,12,1,int(julday))

def obstore_write_datetime(DT,obsfile,nmlfile):
    obstore_write_year(DT.year,obsfile,nmlfile)
    obstore_write_month(DT.month,obsfile,nmlfile)
    obstore_write_day(DT.day,obsfile,nmlfile)
    obstore_write_hour(DT.hour,obsfile,nmlfile)
    obstore_write_minute(DT.minute,obsfile,nmlfile)
    obstore_write_second(DT.second,obsfile,nmlfile)
    obstore_write_julday(DT.timetuple().tm_yday,obsfile,nmlfile)
    
def obstore_read_data_header(obsfile):
    pos = obstore_read_header(obsfile,150,1)
    tcol = obstore_read_header(obsfile,151,1)
    nrow = obstore_read_header(obsfile,152,1)
    pos_data = obstore_read_header(obsfile,160,1)
    size=nrow*tcol
    obsfile.seek((pos-1)*8,0)
    data = obsfile.read(nrow*tcol*8)
    data_header = numpy.asarray(struct.unpack(">"+str(size)+"q", data)).reshape(nrow,tcol)
    return(data_header)
        
def obstore_read_LDC(obsfile,pos,nrow,tcol):
    size=nrow*tcol
    obsfile.seek((pos-1)*8,0)
    data = obsfile.read(nrow*tcol*8)
    val = numpy.asarray(struct.unpack(">"+str(size)+"d", data)).reshape(nrow,tcol)
    if size is 1:
        return(val)[0]
    else:
        return(val)
    
def obstore_read_LUT(obsfile,pos,nrow,tcol):
    size=nrow*tcol
    obsfile.seek((pos-1)*8,0)
    data = obsfile.read(nrow*tcol*8)
    val = numpy.asarray(struct.unpack(">"+str(size)+"q", data)).reshape(nrow,tcol)
    if size is 1:
        return(val)[0]
    else:
        return(val)
        
def obstore_read_subhead(obsfile,pos,nrow,tcol,icol):
    size=nrow*tcol
    #print(pos,nrow,tcol,icol)
    obsfile.seek((pos-1)*8,0)
    data = obsfile.read(nrow*tcol*8)
    val = numpy.asarray(struct.unpack(">"+str(size)+"q", data)).reshape(nrow,tcol)
    if size is 1:
        return(val)[0]
    else:
        return(val)[:,(icol-1)]

def obstore_get_pos(obsfile,p1,p2,p3=None,fmtkey="q"):
    pos = obstore_read_header(obsfile,p1,1,fmtkey=fmtkey)
    tcol = obstore_read_header(obsfile,p2,1,fmtkey=fmtkey)
    if p3 is not None :
       nrow = obstore_read_header(obsfile,p3,1,fmtkey=fmtkey)
    else :
       nrow = 1
    return(pos,tcol,nrow)
    
def obstore_set_hdrpos(obsfile,hdrbetpos=HBpos,hdrbetlen=HBlen,hdrgampos=HCpos,hdrgamlen=HClen,HDR20=HDR20):
    obslib.binary_write(HDR20,1,obsfile)	#	bin_write(obsfile,1,"",HDR20)
    bin_write(obsfile,100,">1q",hdrbetpos)
    bin_write(obsfile,101,">1q",hdrbetlen)
    bin_write(obsfile,105,">1q",hdrgampos)
    bin_write(obsfile,106,">1q",hdrgamlen)
    hdrsize = int(hdrgampos) + int(hdrgamlen) - 1
    return(hdrsize)
    
def obstore_set_batchpos(obsfile,numbatch,batch_data_offset=0,batch_data_length=0,hdrsize=339,maxindx=MAXINDX,lutsize=LUTSIZE):
    #print(hdrsize)
    hdrsize=obstore_set_hdrpos(obsfile)
    ldc_begin=int(hdrsize)+1
    bin_write(obsfile,110,">1q",ldc_begin)
    ldc_ncols=int(maxindx)
    bin_write(obsfile,111,">1q",ldc_ncols)
    ldc_nrows=int(numbatch)
    bin_write(obsfile,112,">1q",ldc_nrows)
    ldc_end=int(ldc_begin)+int(ldc_ncols)*int(ldc_nrows)-1
    rdc_begin=int(ldc_end)+1
    bin_write(obsfile,115,">1q",rdc_begin)
    rdc_ncols=int(maxindx)
    bin_write(obsfile,116,">1q",rdc_ncols)
    rdc_nrows=int(numbatch)
    bin_write(obsfile,117,">1q",rdc_nrows)
    rdc_end=int(rdc_begin)+int(rdc_ncols)*int(rdc_nrows)-1
    cdc_begin=int(rdc_end)+1
    bin_write(obsfile,120,">1q",cdc_begin)
    cdc_ncols=int(maxindx)
    bin_write(obsfile,121,">1q",cdc_ncols)
    cdc_nrows=int(numbatch)
    bin_write(obsfile,122,">1q",cdc_nrows)
    cdc_end=int(cdc_begin)+int(cdc_ncols)*int(cdc_nrows)-1
    lut_begin=int(cdc_end)+1
    bin_write(obsfile,150,">1q",lut_begin)
    lutsize=int(lutsize)
    bin_write(obsfile,151,">1q",lutsize)
    lut_nrows=int(numbatch)
    bin_write(obsfile,152,">1q",lut_nrows)
    lut_end=int(lut_begin)+int(lutsize)*int(lut_nrows)-1
    data_begin=int(lut_end)+1
    bin_write(obsfile,160,">1q",data_begin)
    data_length=batch_data_offset+batch_data_length
    bin_write(obsfile,161,">1q",data_length)
    bin_write(obsfile,data_begin,">1q",1)
    data_end=int(data_begin)+int(data_length)
    obstore_write_datapos(data_begin,obsfile)
    obstore_write_batchcount(numbatch,obsfile)
    file_data_begin = obstore_read_header(obsfile,160,1)
    file_data_length = obstore_read_header(obsfile,161,1)
    file_data_end = int(file_data_begin)+int(file_data_length)-1
    #print(file_data_begin,file_data_length,file_data_end)
    return(data_begin)       #binpos(obsfile,"data")[0]

def obstore_erase_subhead(obsfile,sec_nam,irow):
    (pos,tcol,nrow)=binpos(obsfile,sec_nam)
    size=tcol
    pointer=pos+((irow-1)*tcol)
    obsfile.seek((pointer-1)*8,0)
    fmtstr=struct.Struct(binfmt(sec_nam,size))
    if size is 1:
        obsfile.write(fmtstr.pack(binfill(sec_nam,size)))
    else:
        obsfile.write(fmtstr.pack(*binfill(sec_nam,size)))

def obstore_write_subhead(obsfile,sec_nam,irow,data):
    (pos,tcol,nrow)=binpos(obsfile,sec_nam)
    size=tcol
    pointer=pos+((irow-1)*tcol)
    obsfile.seek((pointer-1)*8,0)
    fmtstr=struct.Struct(binfmt(sec_nam,size))
    if size is 1:
        obsfile.write(fmtstr.pack(data))
    else:
        obsfile.write(fmtstr.pack(*data))
        
def obstore_write_subhead_segment(obsfile,sec_nam,irow,icol,ilen,data):
    size=ilen
    (pos,tcol,nrow)=binpos(obsfile,sec_nam)
    pointer=pos+((irow-1)*tcol)+(icol-1)
    obsfile.seek((pointer-1)*8,0)
    fmtstr=struct.Struct(binfmt(sec_nam,size))
    if size is 1:
        obsfile.write(fmtstr.pack(data))
    else:
        obsfile.write(fmtstr.pack(*data))

def obstore_erase_subhead_segment(obsfile,sec_nam,irow,icol,ilen):
    size=ilen
    (pos,tcol,nrow)=binpos(obsfile,sec_nam)
    pointer=pos+((irow-1)*tcol)+(icol-1)
    print(sec_nam+ " at position "+ str(pointer) + "erased")
    obsfile.seek((pointer-1)*8,0)
    fmtstr=struct.Struct(binfmt(sec_nam,size))
    if size is 1:
        obsfile.write(fmtstr.pack(binfill(sec_nam,size)))
    else:
        obsfile.write(fmtstr.pack(*binfill(sec_nam,size)))

def obstore_read_bin8(obsfile,pointer,size,sec_nam=None):
    ptrpos=((pointer-1)*8)
    obsfile.seek(ptrpos,0)
    data=obsfile.read(size*8)
    if size is 1:
        return struct.unpack(binfmt(sec_nam,size), data)[0]
    else:
        return struct.unpack(binfmt(sec_nam,size), data)


def hdrseclst(ftype):
	hdrsec={
	"fixhdr": ["alpha","beeta","gamma"],
	"hdr": ["hdrldc","hdrrdc","hdrcdc","hdrlut"],
	"obstore": ["alpha","beeta","gamma","ldc","rdc","cdc","lut"]
	}
	return(hdrsec.get(ftype,["data"]))

    
def obstore_read_subhead_segment(obsfile,sec_nam,irow=1,icol=1,ilen=None,batchid=None,maxindx=MAXINDX,lutsize=LUTSIZE):
    sec_tag=sec_nam.replace("hdr","")
    if sec_nam == "hdralp" :
	pos = 1
	tcol = 256
	nrow = 1
    if sec_nam == "hdrbet" :
	pos = 257
	tcol = 49
	nrow = 1
    if sec_nam == "hdrgam" :
	pos = 307
	tcol = 34
	nrow = 1
    if sec_tag == "ldc" :
	pos = 1
	tcol = maxindx
	nrow = 1
    if sec_tag == "rdc" :
	pos = 1+int(maxindx)
	tcol = maxindx
	nrow = 1
    if sec_tag == "cdc" :
	pos = 1+int(maxindx)*2
	tcol = maxindx
	nrow = 1
    if sec_tag == "lut" :
	pos = 1+int(maxindx)*3
	tcol = lutsize
	nrow = 1
    if not sec_nam in ["hdralp","hdrbet","hdrgam","hdrldc","hdrrdc","hdrcdc","hdrlut"] :
    	(pos,tcol,nrow)=binpos(obsfile,sec_nam)
    if ilen is not None : size = ilen
    else : size = tcol
    if batchid is not None: irow=batchid
    if sec_nam in ["alpha", "beeta", "gamma", "hdralp","hdrbet","hdrgam"] : irow = 1
    pointer=pos+((irow-1)*tcol)+(icol-1)
    hdrinfo=obstore_read_bin8(obsfile,pointer,size,sec_nam=sec_nam)
    return(hdrinfo)
    
def obstore_read_LookUpTable(obsfile,batchid=1,lutsize=LUTSIZE):
    LUT=obstore_read_subhead_segment(obsfile,"lut",batchid,1,lutsize)
    return(LUT)
##########following code is depreciated due to bugg #######
#    (pos,tcol,nrow)=binpos(obsfile,'lut')
#    pointer=pos+((irow-1)*tcol)
#    obsfile.seek((pointer-1)*8,0)
#    data = obsfile.read(tcol*8)
#    size=tcol-64
#    return struct.unpack(">62q2d"+str(size)+"q", data)

def obstore_read_batch_header(obsfile,batchid=1,lutsize=LUTSIZE,maxindx=MAXINDX,ftype="obstore",sec_lst=None):
    #alpha=obstore_read_subhead_segment(obsfile,"alpha",1,1,256)
    #beeta=obstore_read_subhead_segment(obsfile,"beeta",1,1,49)
    #gamma=obstore_read_subhead_segment(obsfile,"gamma",1,1,34)
    #LDC=obstore_read_subhead_segment(obsfile,"ldc",batchid,1,maxindx)
    #RDC=obstore_read_subhead_segment(obsfile,"rdc",batchid,1,maxindx)
    #CDC=obstore_read_subhead_segment(obsfile,"cdc",batchid,1,maxindx)
    #LUT=obstore_read_subhead_segment(obsfile,"lut",batchid,1,lutsize)
    if ftype in ["fixhdr","obstore"]: maxindx=obslib.binary_read(obsfile,111,1)
    if maxindx is None : maxindx=MAXINDX
    if sec_lst is None: sec_lst=hdrseclst(ftype)
    hdr_info={}
    for indx in range(len(sec_lst)) :
	sec_nam=sec_lst[indx]
	sec_tag=sec_nam.replace("hdr","")
    	hdr_info[sec_tag]=obstore_read_subhead_segment(obsfile,sec_nam,batchid=batchid,maxindx=maxindx,lutsize=lutsize)
	#print(hdr_info[sec_tag])
    return(hdr_info)

def write_batchheader(DT,obsfile,nmlfile,subtype,elist,batchid=1,batch_data_offset=0,obscount=1,batchcount=1,hdrsize=339,maxindx=MAXINDX,lutsize=LUTSIZE,Tref=TREF,HlfTW=HLFTW):
    halftw=obslib.timeperiod(hh=int(HlfTW))
    stime=int(obslib.mins_sinceTref(DT-halftw,Tref=Tref))
    etime=int(obslib.mins_sinceTref(DT+halftw,Tref=Tref)-1)
    print(stime,etime)
    datapos=obstore_set_batchpos(obsfile,batchcount,batch_data_offset=batch_data_offset,hdrsize=hdrsize,maxindx=maxindx,lutsize=lutsize)    
    elist=obstore_write_elist(obsfile,nmlfile,batchid,subtype,elist,maxindx=maxindx)
    (obscount,reclen,batchlen)=obstore_write_lut(obsfile,batchid,subtype,elist,stime=stime,etime=etime,obscount=obscount,maxindx=maxindx)
    datalen=update_datalen(obsfile,batchid=batchid,batch_data_offset=batch_data_offset)
    #print(str(obscount)+ " observations of length "+str(reclen)+ " make a total length of "+str(batchlen)+ " for the batch "+str(batchid)+"." )
    #print("Data start position is "+str(datapos))
    #batch_data_begin=datapos+batch_data_offset
    #datalen=batch_data_offset+obstore_read_subhead_segment(obsfile,"lut",idx,30,1)
    #print("(Batch "+str(idx)+ ") length : "+str(batchlen))
    #print("Updated data length : "+str(datalen))
    #dataend=datapos+datalen
    #print("Data end position is "+str(dataend))
    #obstore_write_datalen(datalen,obsfile,nmlfile)
    obstore_write_datetime(DT,obsfile,nmlfile)
    obstore_write_stime(stime,obsfile)
    obstore_write_etime(etime,obsfile)
    return(datalen)

def update_datalen(obsfile,batchid=0,batch_data_offset=0,):
    datapos=obstore_read_header(obsfile,160,1)
    batchlen=obstore_read_subhead_segment(obsfile,"lut",batchid,30,1)
    datalen=batch_data_offset+batchlen
    obstore_write_datalen(datalen,obsfile)
    dataend=datapos+datalen
    print("(Batch "+str(batchid)+ ") length : "+str(batchlen))
    print("Updated data length : "+str(datalen))
    print("Data end position is "+str(dataend))
    return(datalen)

def batch_header_read(input_file,nmlfile,maxindx=MAXINDX,batchid=1,lutsize=LUTSIZE):
    with open(input_file, "rb") as obsfile:
        batchinfo=obstore_headerinfo(obsfile,nmlfile,maxindx)
        #print(obstore_read_LookUpTable(obsfile,batchid))
        #print(batchinfo.Batch_Data_Position[0],batchinfo.Batch_Data_Length[0],batchinfo.Batch_Data_End[0])
    return(batchinfo.elist_map[batchid])


def obstore_read_data(obsfile):
    (pos,data_len,void)=obstore_get_pos(obsfile,160,161,162)
    obsfile.seek((pos-1)*8,0)
    data = obsfile.read(data_len*8)
    val = struct.unpack(">"+str(data_len)+"d", data)
    return(val)

def read_data(input_file):
    with open(input_file, "rb") as obsfile:
	data=obstore_read_data(obsfile)
    return(data)
    
def obstore_erase_data(obsfile):
    pos = obstore_read_header(obsfile,160,1)
    data_len = obstore_read_header(obsfile,161,1)
    data_len+=100  ### Temperory testing purpose
    void = obstore_read_header(obsfile,162,1)
    obsfile.seek((pos-1)*8,0)
    frmtstr=struct.Struct(">"+str(data_len)+"d")
    data=[-1073741824.0] * data_len
    obsfile.write(frmtstr.pack(*data))

def obstore_close_data(obsfile):
    pos = obstore_read_header(obsfile,160,1)
    data_len = obstore_read_header(obsfile,161,1)
    data_end = int(pos)+int(data_len)
    obsfile.seek((data_end)*8,0)
    frmtstr=struct.Struct(">"+str(8)+"q")
    data=[0]*8
    obsfile.write(frmtstr.pack(*data))
    return(pos,data_len,data_end) 

def obstore_read_EOF(obsfile):
    pos = obstore_read_header(obsfile,160,1)
    data_len = obstore_read_header(obsfile,161,1)
    data_end = int(pos)+int(data_len)
    obsfile.seek((data_end)*8,0)
    data = obsfile.read(12)
    val = struct.unpack(">"+str(12)+"s", data)
    return(str(val[0]))
    
def obstore_get_tail_length(obsfile,seekpos):
    pos = obstore_read_header(obsfile,160,1)
    data_len = obstore_read_header(obsfile,161,1)
    data_end = int(pos)+int(data_len)
    tail_len= int(data_end)-int(seekpos)
    return(tail_len)

def obstore_check_data_length(obsfile,seekpos):
    pos = obstore_read_header(obsfile,160,1)
    data_len = obstore_read_header(obsfile,161,1)
    data_end = int(pos)+int(data_len)
    tail_len= int(data_end)-int(seekpos)
    if diaglev > 9 : errprint(seekpos,data_end,tail_len)
    if int(tail_len) < 0:
        new_len=int(seekpos)+1-int(pos)
        if diaglev > 9 : errprint("Correcting data_lenth, Old: %s, New: %s"%(data_len,new_len))
        obstore_overwrite_formatted(new_len,">1q",obsfile,161)
    return(obstore_read_header(obsfile,161,1))

def obstore_read_subtype(obsfile):
    (pos_lut,tcol_lut,nrow_lut)=obstore_get_pos(obsfile,150,151,152)
    subtypelist=obstore_read_subhead(obsfile,pos_lut,nrow_lut,tcol_lut,68)
    subtypelist=obslib.binsort(subtypelist,missing=-32768).flatten()
    return(subtypelist)

def obstore_write_subtype(obsfile,irow,subtype):
    obstore_write_subhead_segment(obsfile,"lut",irow,68,1,int(subtype))

def obstore_write_stime_lut(obsfile,irow,stime):
    obstore_write_subhead_segment(obsfile,"lut",irow,69,1,int(stime))

def obstore_write_etime_lut(obsfile,irow,etime):
    obstore_write_subhead_segment(obsfile,"lut",irow,70,1,int(etime))

def obstore_write_subtype_nrows(obsfile,irow,nrows,cscnt=NAN_VAL_INT):
    obstore_write_subhead_segment(obsfile,"lut",irow,19,1,int(nrows))
    obstore_write_subhead_segment(obsfile,"lut",irow,66,1,int(nrows))

def obstore_write_subtype_tcols(obsfile,irow,tcols,cspos=NAN_VAL_INT):
    obstore_write_subhead_segment(obsfile,"lut",irow,18,1,int(tcols))
    obstore_write_subhead_segment(obsfile,"lut",irow,67,1,int(tcols))

def obstore_write_subtype_ncells(obsfile,irow,ncells):
    obstore_write_subhead_segment(obsfile,"lut",irow,15,1,int(ncells))
    obstore_write_subhead_segment(obsfile,"lut",irow,30,1,int(ncells))

def obstore_write_subtype_filepos(obsfile,irow):
    pos_data = obstore_read_header(obsfile,160,1)
    filepos=pos_data-1
    for i in range(1,irow,1):
        filepos+=obstore_read_subhead_segment(obsfile,"lut",i,30,1)
    obstore_write_subhead_segment(obsfile,"lut",irow,29,1,int(filepos))

def obstore_write_subtype_bodypos(obsfile,irow):
    bodypos=1
    for i in range(1,irow,1):
        bodypos+=obstore_read_subhead_segment(obsfile,"lut",i,15,1)
    obstore_write_subhead_segment(obsfile,"lut",irow,40,1,int(bodypos))

def obstore_write_subtype_lut38(obsfile,irow,val=NAN_VAL_INT):
    val=1111
    obstore_write_subhead_segment(obsfile,"lut",irow,38,1,int(val))

def obstore_write_subtype_lut39(obsfile,irow,val=NAN_VAL_INT):
    val=1
    obstore_write_subhead_segment(obsfile,"lut",irow,39,1,int(val))

def obstore_write_subtype_lut42(obsfile,irow,val=NAN_VAL_INT):
    val=1
    obstore_write_subhead_segment(obsfile,"lut",irow,42,1,int(val))

def obstore_write_subtype_lut_zero(obsfile,irow,fillval=NAN_VAL_INT,lutsize=LUTSIZE):
    obstore_write_subhead_segment(obsfile,"lut",irow,13,2,[fillval,fillval])
    obstore_write_subhead_segment(obsfile,"lut",irow,16,2,[fillval,fillval])
    obstore_write_subhead_segment(obsfile,"lut",irow,20,2,[0,0])
    obstore_write_subhead_segment(obsfile,"lut",irow,22,1,int(fillval))
    obstore_write_subhead_segment(obsfile,"lut",irow,23,6,[0,0,0,0,0,0])
    obstore_write_subhead_segment(obsfile,"lut",irow,31,7,[0,0,0,0,0,0,0])
    obstore_write_subhead_segment(obsfile,"lut",irow,41,1,0)
    obstore_write_subhead_segment(obsfile,"lut",irow,43,20,[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
    obstore_write_subhead_segment(obsfile,"lut",irow,65,1,int(fillval))
    laseglen=int(lutsize-71)
    obstore_write_subhead_segment(obsfile,"lut",irow,71,laseglen,[fillval]*laseglen)

def obstore_write_subtype_lut63(obsfile,irow,val=NAN_VAL_INT):
    val=int(-4481081629233643520)
    ### operational obstore have fixed value ##
    obstore_write_subhead_segment(obsfile,"lut",irow,63,1,val)

def obstore_write_subtype_lut64(obsfile,irow,val=NAN_VAL_INT):
    val=int(4607182418800017408)
    ### operational obstore have fixed value ##
    obstore_write_subhead_segment(obsfile,"lut",irow,64,1,val)
    
def obstore_read_batchinfo(obsfile,indx):
    pos_data = obstore_read_header(obsfile,160,1)
    data_header=obstore_read_data_header(obsfile)
    #print(indx, pos_data, data_header[0,14], pos_data+data_header[0,14])
    for i in range(0,(indx-1),1):
        #print(data_header[i,14],obstore_read_subhead_segment(obsfile,"lut",i,15,1))
        pos_data=pos_data+data_header[i,14]
    obs_count=data_header[indx-1,18]
    obs_nele=data_header[indx-1,17]
    data_len=data_header[indx-1,14]
    data_end=pos_data+data_len-1
    return(indx,pos_data,obs_count,obs_nele,data_len,data_end)

def obstore_batchinfo(obsfile,nmlfile,indx,maxindx=MAXINDX):
    subtype=obstore_read_index_subtype(obsfile,indx)
    (indx,pos_data,obs_count,obs_nele,data_len,data_end)=obstore_read_batchinfo(obsfile,indx)
    elist=obstore_read_batch_elements(obsfile,indx,nmlfile,maxindx)
    batchinfo=pandas.DataFrame([[indx,subtype,pos_data,obs_count,obs_nele,data_len,data_end,elist]],index=[indx],columns=["Batch_Index","Batch_Subtype","Batch_Data_Position","Batch_Obs_Count","Batch_Col_Count","Batch_Data_Length","Batch_Data_End","elist_map"])
    return(batchinfo)

def print_batchinfo(batchinfo):
    for indx in batchinfo.Batch_Index.values:
        if indx not in [0]:
            print(indx)
            print(batchinfo[["Batch_Index","Batch_Subtype","Batch_Data_Position","Batch_Obs_Count","Batch_Col_Count","Batch_Data_Length","Batch_Data_End"]].loc[indx])
            print(batchinfo['elist_map'].loc[indx])

def obstore_headerinfo(obsfile,nmlfile,maxindx=MAXINDX):
    file_data_begin = obstore_read_header(obsfile,160,1)
    file_data_length = obstore_read_header(obsfile,161,1)
    file_data_end = int(file_data_begin)+int(file_data_length)-1
    subtype_list=obstore_read_subtype(obsfile)
    #print(file_data_begin,file_data_length,file_data_end)
    batchinfo=pandas.DataFrame([[0,"Full Data",file_data_begin,"----","----",file_data_length,file_data_end,subtype_list]],index=[0],columns=["Batch_Index","Batch_Subtype","Batch_Data_Position","Batch_Obs_Count","Batch_Col_Count","Batch_Data_Length","Batch_Data_End","elist_map"])
    for indx,subtype in enumerate(subtype_list,start=1):
        batchinfo=batchinfo.append(obstore_batchinfo(obsfile,nmlfile,indx,maxindx))
    return(batchinfo)

def obstore_write_header(obsfile,nmlfile,batchinfo,batch_count=None):
    if batch_count is None: batch_count=numpy.max(batchinfo.Batch_Index.values)
    file_data_pos=obstore_set_batchpos(obsfile,batch_count)
    file_data_begin = batchinfo.Batch_Data_Position[0]
    file_data_end = batchinfo.Batch_Data_End[batch_count]
    file_data_length = int(file_data_end) - int(file_data_begin)
    bin_write(obsfile,160,">1q",file_data_begin)
    bin_write(obsfile,161,">1q",file_data_length)
    bin_write(obsfile,284,">1q",file_data_length)
    bin_write(obsfile,288,">1q",batch_count)
    bin_write(obsfile,305,">1q",batch_count)
    for indx in batchinfo.Batch_Index.values:
        if indx not in [0]:
            obstore_write_batch_header(obsfile,nmlfile,batchinfo,indx)
    return(obstore_read_subtype(obsfile))

def obstore_write_lut(obsfile,irow,subtype,elist,stime=0,etime=0,obscount=0,maxindx=MAXINDX):
    reclen=int(sum(elist.LDC.values))
    batchlen=int(obscount*reclen)      ###int(obs_count*(sum(elist.LDC.values)))
    obstore_write_subtype_lut_zero(obsfile,irow)
    obstore_write_subtype(obsfile,irow,subtype)
    obstore_write_stime_lut(obsfile,irow,stime)
    obstore_write_etime_lut(obsfile,irow,etime)
    obstore_write_subtype_nrows(obsfile,irow,obscount)
    obstore_write_subtype_tcols(obsfile,irow,reclen)
    obstore_write_subtype_ncells(obsfile,irow,batchlen)
    obstore_write_subtype_filepos(obsfile,irow)
    obstore_write_subtype_bodypos(obsfile,irow)
    obstore_write_subtype_lut38(obsfile,irow)
    obstore_write_subtype_lut39(obsfile,irow)
    obstore_write_subtype_lut42(obsfile,irow)
    obstore_write_subtype_lut63(obsfile,irow)
    obstore_write_subtype_lut64(obsfile,irow)
    return(obscount,reclen,batchlen)

def obstore_write_batch_header(obsfile,nmlfile,batchinfo,irow,maxindx=MAXINDX):
    subtype=batchinfo.Batch_Subtype[irow]
    obscount=batchinfo.Batch_Obs_Count[irow]
    elist=batchinfo.elist_map[irow]
    stime=batchinfo.StartTime[irow]
    etime=batchinfo.EndTime[irow]
    obstore_write_elist(obsfile,nmlfile,irow,subtype,elist)
    ########LookUpTable Editing########
    obstore_erase_subhead_segment(obsfile,"lut",irow,1,LUTSIZE)
    obstore_write_lut(obsfile,irow,subtype,elist,stime=stime,etime=etime,obscount=obscount,maxindx=maxindx)
    return(obstore_read_batch_elements(obsfile,irow,nmlfile,maxindx))
    
def obstore_read_index_subtype(obsfile,indx):
    return(obstore_read_subhead_segment(obsfile,"lut",int(indx),68,1))

def obstore_read_subtype_index(obsfile,subtype):
    indx=numpy.where(obstore_read_subtype(obsfile) == int(subtype))[0][0] + 1
    return(indx)

def obstore_write_subtype_index(obsfile,nmlfile,irow,subtype,obs_count,elist,maxindx=MAXINDX):
    obstore_write_elist(obsfile,nmlfile,irow,subtype,elist)
    obstore_write_subtype(obsfile,irow,subtype)
    obstore_write_subtype_nrows(obsfile,irow,obs_count)
    obstore_write_subtype_tcols(obsfile,irow,int(sum(elist.LDC.values)))
    obstore_write_subtype_ncells(obsfile,irow,int(obs_count*(sum(elist.LDC.values))))
    obstore_write_subtype_filepos(obsfile,irow)
    obstore_write_subtype_bodypos(obsfile,irow)
    obstore_write_subtype_lut38(obsfile,irow)
    obstore_write_subtype_lut39(obsfile,irow)
    return(obstore_read_batch_elements(obsfile,irow,nmlfile,maxindx))

def obstore_read_station_record_position(obsfile,subtype,StnNo,nmlfile):
    WMOStnNo = obstore_read_data_element(obsfile,nmlfile,indx,"WMOStnNo")
    return(WMOStnNo.index[WMOStnNo['WMOStnNo'] == StnNo].tolist())
    
def obstore_read_data_subtype(obsfile,subtype):
    (indx,pos_data,obs_count,tcols,data_len,data_end)=obstore_read_batchinfo(obsfile,indx)
    obsfile.seek((pos_data-1)*8,0)
    data = obsfile.read(data_len*8)
    val = numpy.asarray(struct.unpack(">"+str(data_len)+"d", data)).reshape(obs_nele,obs_count)
    if data_len is 1:
        return(val[0])
    else:
        return(val[:,:])  #(val)[:,:]

def frame_batch_elist(hdrinfodic,nmlfile=NML_OBS_INDX,maxindx=MAXINDX):
    ldc=numpy.array(hdrinfodic["ldc"])
    rdc=numpy.array(hdrinfodic["rdc"])
    cdc=numpy.array(hdrinfodic["cdc"])
    cdc_sortlist=cdc
    #print(cdc)
    cdc_sortlist=obslib.binsort(cdc_sortlist)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=numpy.nan)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=-4481081629233643520)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=4607182418800017408)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=-1073741820.0)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=-1073741824.0)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=-1.07374182e+09)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=-3.2768e+04)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=-32768.0)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=-32768)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=0.)
    cdc_sortlist=obslib.binsort(cdc_sortlist,binmin=1.0)
    for i in cdc_sortlist:
      idx_arr=numpy.where(cdc == i)[0]
      if idx_arr.size > 0:
        obs_index=(idx_arr[0])+1
	#if rdc[obs_index-1] < 0 : rdc[obs_index-1] = 0
        LDC=int(ldc[obs_index-1])
        RDC=int(rdc[obs_index-1])
        CDC=int(cdc[obs_index-1])
	if obs_index > maxindx : print("Index not handled :"+obs_index)
        with open(nmlfile, "r") as nml:
            Element= pandas.read_table(nml, skiprows=None, header=0).query("eleindex == @obs_index").elename.reset_index(drop=True)[0]
        if i == cdc_sortlist[0]: 
            elist=pandas.DataFrame([[CDC,RDC,LDC,obs_index,Element]],index=[i],columns=["CDC","RDC","LDC","obs_index","Element"])
        else:
            elist=elist.append(pandas.DataFrame([[CDC,RDC,LDC,obs_index,Element]],index=[i],columns=["CDC","RDC","LDC","obs_index","Element"]))
    #print(elist)
    return(elist)

def obstore_read_batch_elements(obsfile,obsidx,nmlfile=NML_OBS_INDX,maxindx=MAXINDX):
    (pos_ldc,tcol_ldc,nrow_ldc)=obstore_get_pos(obsfile,110,111,112)
    (pos_rdc,tcol_rdc,nrow_rdc)=obstore_get_pos(obsfile,115,116,117)
    (pos_cdc,tcol_cdc,nrow_cdc)=obstore_get_pos(obsfile,120,121,122)
    (pos_lut,tcol_lut,nrow_lut)=obstore_get_pos(obsfile,150,151,152)
    ldc = numpy.array(obstore_read_subhead_segment(obsfile,"ldc",int(obsidx),1,maxindx)) 
    #obstore_read_LDC(obsfile,pos_ldc,nrow_ldc,tcol_ldc)[obsidx,:]
    rdc = numpy.array(obstore_read_subhead_segment(obsfile,"rdc",int(obsidx),1,maxindx))
    #obstore_read_LDC(obsfile,pos_rdc,nrow_rdc,tcol_rdc)[obsidx,:]
    cdc = numpy.array(obstore_read_subhead_segment(obsfile,"cdc",int(obsidx),1,maxindx)) 
    #obstore_read_LDC(obsfile,pos_cdc,nrow_cdc,tcol_cdc)[obsidx,:]
    cdc_sortlist=cdc
    #print(cdc)
    cdc_sortlist=obslib.binsort(cdc_sortlist)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=numpy.nan)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=-4481081629233643520)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=4607182418800017408)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=-1073741820.0)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=-1073741824.0)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=-1.07374182e+09)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=-3.2768e+04)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=-32768.0)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=-32768)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=0.)
    cdc_sortlist=obslib.binsort(cdc_sortlist,binmin=1.0)
    subtype=obstore_read_subtype(obsfile)[obsidx-1]
    for i in cdc_sortlist:
      idx_arr=numpy.where(cdc == i)[0]
      if idx_arr.size > 0:
        obs_index=(idx_arr[0])+1
	#if rdc[obs_index-1] < 0 : rdc[obs_index-1] = 0
        LDC=int(ldc[obs_index-1])
        RDC=int(rdc[obs_index-1])
        CDC=int(cdc[obs_index-1])
	if obs_index > maxindx : print("Index not handled :"+obs_index)
        with open(nmlfile, "r") as nml:
            Element= pandas.read_table(nml, skiprows=None, header=0).query("eleindex == @obs_index").elename.reset_index(drop=True)[0]
        if i == cdc_sortlist[0]: 
            elist=pandas.DataFrame([[CDC,RDC,LDC,obs_index,Element]],index=[i],columns=["CDC","RDC","LDC","obs_index","Element"])
        else:
            elist=elist.append(pandas.DataFrame([[CDC,RDC,LDC,obs_index,Element]],index=[i],columns=["CDC","RDC","LDC","obs_index","Element"]))
    #print(elist)
    return(elist)
    
def element_frame(idx,CDC,RDC,LDC,obs_index,nmlfile):
    #print(idx,CDC,RDC,LDC,obs_index)
    with open(nmlfile, "r") as nml:
        Element= pandas.read_table(nml, skiprows=None, header=0).query("eleindex == @obs_index").elename.reset_index(drop=True)[0]
    eframe=pandas.DataFrame([[CDC,RDC,LDC,obs_index,Element]],index=[idx],columns=["CDC","RDC","LDC","obs_index","Element"])
    return(eframe)

def obstore_create_element_table(nmlfile,obsidxlist):
    idx=0
    elenamlst=[]
    for i,obs_index in enumerate(obsidxlist,1):
        if obs_index not in elenamlst:
            idx+=1
            if obs_index not in obsidxlist[i:len(obsidxlist)]:
                LDC=1
                RDC=1
                CDC=i
                elenamlst.append(obs_index)
            else:
                j=obsidxlist[i:len(obsidxlist)].index(obs_index)+1
                k=obsidxlist[i:len(obsidxlist)].count(obs_index)+1
                LDC=k
                RDC=j
                CDC=i
                elenamlst.append(obs_index)
            if idx == 1: 
                elist=element_frame(idx,CDC,RDC,LDC,obs_index,nmlfile)
            else:
                elist=elist.append(element_frame(idx,CDC,RDC,LDC,obs_index,nmlfile))
    return(elist)

def obstore_read_subtype_elements(obsfile,subtype,nmlfile,maxindx=MAXINDX):
    obsidx = numpy.where(obstore_read_subtype(obsfile) == int(subtype))[0][0]
    return(obstore_read_batch_elements(obsfile,obsidx,nmlfile,maxindx))
    
def obstore_write_elist(obsfile,nmlfile,irow,subtype,elist,maxindx=MAXINDX):
    if diaglev > 0 : print(elist)
    (pos_ldc,tcol_ldc,nrow_ldc)=obstore_get_pos(obsfile,110,111,112)
    (pos_rdc,tcol_rdc,nrow_rdc)=obstore_get_pos(obsfile,115,116,117)
    (pos_cdc,tcol_cdc,nrow_cdc)=obstore_get_pos(obsfile,120,121,122)
    (pos_lut,tcol_lut,nrow_lut)=obstore_get_pos(obsfile,150,151,152)
    obstore_erase_subhead_segment(obsfile,"ldc",irow,1,maxindx)
    obstore_erase_subhead_segment(obsfile,"rdc",irow,1,maxindx)
    obstore_erase_subhead_segment(obsfile,"cdc",irow,1,maxindx)
    for i,icol in enumerate(elist.obs_index.values):
        obstore_write_subhead_segment(obsfile,"ldc",irow,icol,1,elist.LDC.values[i])
        obstore_write_subhead_segment(obsfile,"rdc",irow,icol,1,elist.RDC.values[i])
        obstore_write_subhead_segment(obsfile,"cdc",irow,icol,1,elist.CDC.values[i])
    return(obstore_read_batch_elements(obsfile,irow,nmlfile,maxindx))
    
def obstore_read_element(obsfile,elist,element,pos_data,record_pos,tcols,nmlfile):
    elem_name_list=elist.Element.values
    #print(elem_name_list)
    if element in elem_name_list:
        if diaglev > 9: print(element)
        (rdc,cdc,ldc) = elist.query("Element == @element")[["RDC","CDC","LDC"]].reset_index(drop=True).values[0]
        ele_pos_offset=int(cdc)-1
        if diaglev > 100 : errprint(elist, element, cdc, ldc, rdc)
        if element in ["CharData"]:
            ncols=1
            frmtstr=struct.Struct(">"+str(ldc)+"d")
            data_size=int(ldc)*8
            if diaglev > 100 : errprint("Reading CharData formated %s with size %s values"%(frmtstr,data_size))
        else:
            ncols=int(ldc)
            frmtstr=struct.Struct(">"+str(1)+"d")
            data_size=8
        if diaglev > 100: errprint("ncols ="+str(ncols))
        if ncols is 1 :
            levs=[1]
            data_element=pandas.DataFrame(index=record_pos,columns=[element])
        else:
            levs=range(1,ncols+1,1)
            data_element=pandas.DataFrame(index=record_pos,columns=[element+str(i) for i in range(1,ncols+1,1)])
        if diaglev > 9: errprint("ncols ="+str(ncols)+"levs = %s"%(levs)+"data_size ="+str(data_size))
        if element in ["CharData"]:
            for i in record_pos:
                seekpos=pos_data+((i-1)*tcols)+int(cdc-1)
                obsfile.seek((seekpos-1)*8,0)
                data = obsfile.read(data_size)
                data_element.xs(i)[0] = obslib.getstring(frmtstr.unpack(data))
                if diaglev > 50: errprint(data_element.xs(i)[0])
                if diaglev > 100: errprint("unpack = %s"%(str(frmtstr.unpack(data))))
        else:
            for i in record_pos:
                for j in levs:
                    seekpos=pos_data+((i-1)*tcols)+int(cdc-1)+((j-1)*rdc)
                    if diaglev > 100 : errprint("pos_data=%s,i=%s,tcols=%s,cdc=%s,j=%s,rdc =%s,seekpos = %s"%(pos_data,i,tcols,cdc,j,rdc,seekpos))
                    obsfile.seek((seekpos-1)*8,0)
                    if diaglev > 100: errprint("Remaining bytes = "+str(obstore_get_tail_length(obsfile,seekpos)))
                    data = obsfile.read(data_size)
                    #errprint(data_size)
                    errprint("Check point 6: data = %s"%(data))
                    #print(i,j,data_size,frmtstr)
                    data_element.xs(i)[j-1] = numpy.asarray(frmtstr.unpack(data))[0]
        return(data_element)
    else:
        errprint("Element %s is not available "%(element))
     

def obstore_read_element_level(obsfile,subtype,element,pos_data,record_pos,obs_nele,lev_pos,nmlfile):
    elist = obstore_read_elements(obsfile,int(subtype),nmlfile)
    (rdc,cdc,ldc) = elist.query("Element == @element")[["RDC","CDC","LDC"]].reset_index(drop=True).values[0]
    ele_pos_offset=int(cdc)-1
    ncols=int(ldc)
    data_element=pandas.DataFrame(index=record_pos,columns=[element])
    for a,i in enumerate(record_pos):
                lev_idx=lev_pos[a]-1
                seekpos=pos_data+((i-1)*obs_nele)+int(cdc-1)+((lev_pos[a])*rdc)
                obsfile.seek((seekpos-1)*8,0)
                data = obsfile.read(1*8)
                data_element.xs(i)[a] = numpy.asarray(struct.unpack(">"+str(1)+"d", data))[0]
    return(data_element) 

def obstore_write_element_level(obsfile,nmlfile,subtype,elist,element,pos_data,record_pos,lev_pos,data,maxindx=MAXINDX):
    tcols=sum(elist.CDC.values)
    (rdc,cdc,ldc) = elist.query("Element == @element")[["RDC","CDC","LDC"]].reset_index(drop=True).values[0]
    ele_pos_offset=int(cdc)-1
    ncols=int(ldc)
    for a,i in enumerate(record_pos):
                lev_idx=lev_pos[a]-1
                seekpos=pos_data+((i-1)*tcols)+int(cdc-1)+((lev_idx)*rdc)
                obsfile.seek((seekpos-1)*8,0)
                fmtstr=struct.Struct(binfmt("data",1))
                obsfile.write(fmtstr.pack(data[a]))
      
def obstore_read_data_element(obsfile,nmlfile,indx,element,maxindx=MAXINDX):
    elist = obstore_read_batch_elements(obsfile,indx,nmlfile,maxindx)
    (indx,pos_data,obs_count,tcols,data_len,data_end)=obstore_read_batchinfo(obsfile,indx)
    #print("indx,pos_data,obs_count,tcols,data_len,data_end")
    #print(indx,pos_data,obs_count,tcols,data_len,data_end)
    #print(elist)
    if element in elist.Element.values:
        if diaglev > 0 : errprint("Batch No: %s ; Element : %s"%(indx,element))
        record_pos=range(1,obs_count+1,1)
        return(obstore_read_element(obsfile,elist,element,pos_data,record_pos,tcols,nmlfile))
    else:
        if diaglev > 0 : errprint("Batch No: %s does not contain Element: %s "%(indx,element))

def obstore_write_data_element(obsfile,nmlfile,indx,element,data,maxindx=MAXINDX):
    elist = obstore_read_batch_elements(obsfile,indx,nmlfile,maxindx)
    (indx,pos_data,obs_count,tcols,data_len,data_end)=obstore_read_batchinfo(obsfile,indx)
    if element in elist.Element.values:
        record_pos=range(1,obs_count+1,1)
        return(obstore_write_element(obsfile,indx,elist,element,pos_data,record_pos,tcols,data,maxindx=maxindx))

def obstore_erase_data_element(obsfile,nmlfile,indx,element,maxindx=MAXINDX):
    elist = obstore_read_batch_elements(obsfile,indx,nmlfile,maxindx)
    (indx,pos_data,obs_count,tcols,data_len,data_end)=obstore_read_batchinfo(obsfile,indx)
    if element in elist.Element.values:
        record_pos=range(1,obs_count+1,1)
        return(obstore_erase_element(obsfile,indx,elist,element,pos_data,record_pos,tcols))

def obstore_read_data_record(obsfile,indx=1,record_pos=[1],nmlfile=NML_OBS_INDX):
    elist = obstore_read_batch_elements(obsfile,indx,nmlfile)
    subtype=obstore_read_index_subtype(obsfile,indx)
    (indx,pos_data,obs_count,tcols,data_len,data_end)=obstore_read_batchinfo(obsfile,indx)
    data_len=tcols #  since only one record is fetching at a time.
    if data_len is 1:
        data_element=pandas.DataFrame(index=record_pos,columns=[element])
    else:
        data_element=pandas.DataFrame(index=record_pos,columns=range(1,data_len+1,1))

    for i in record_pos:
        seekpos=pos_data+((i-1)*tcols)
        obsfile.seek((seekpos-1)*8,0)
        data = obsfile.read(data_len*8)
        data_element.xs(i)[:] = numpy.asarray(struct.unpack(">"+str(data_len)+"d", data))
    return(data_element)

def obstore_erase_data_record(obsfile,subtype,irow,icol,ilen):
    (indx,pos_data,obs_count,tcols,data_len,data_end)=obstore_read_batchinfo(obsfile,indx)
    size=ilen
    pointer=pos_data+((irow-1)*obs_nele)+(icol-1)
    obsfile.seek((pointer-1)*8,0)
    fmtstr=struct.Struct(binfmt("data",size))
    if size is 1:
        obsfile.write(fmtstr.pack(binfill("data",size)))
    else:
        obsfile.write(fmtstr.pack(*binfill("data",size)))

def obstore_write_data_record(obsfile,subtype,irow,icol,ilen,data):
    (indx,pos_data,obs_count,tcols,data_len,data_end)=obstore_read_batchinfo(obsfile,indx)
    size=ilen
    pointer=pos_data+((irow-1)*obs_nele)+(icol-1)
    obsfile.seek((pointer-1)*8,0)
    fmtstr=struct.Struct(binfmt("data",size))
    if size is 1:
        obsfile.write(fmtstr.pack(data))
    else:
        obsfile.write(fmtstr.pack(*data))

def obstore_read_data_record_element(obsfile,subtype,record_pos,element,nmlfile):
    elist = obstore_read_elements(obsfile,int(subtype),nmlfile)
    (indx,pos_data,obs_count,tcols,data_len,data_end)=obstore_read_batchinfo(obsfile,indx)
    return(obstore_read_element(obsfile,elist,element,pos_data,record_pos,obs_nele,nmlfile))

def obstore_read_data_station(obsfile,subtype,StnNo,nmlfile):
    elist = obstore_read_elements(obsfile,int(subtype),nmlfile)
    (indx,pos_data,obs_count,tcols,data_len,data_end)=obstore_read_batchinfo(obsfile,indx)
    data_len=obs_nele #  since only one record is fetching at a time.
    record_pos = obstore_read_station_record_position(obsfile,subtype,StnNo,nmlfile)
    if data_len is 1:
        data_element=pandas.DataFrame(index=record_pos,columns=[element])
    else:
        data_element=pandas.DataFrame(index=record_pos,columns=range(1,data_len+1,1))
    for i in record_pos:
        seekpos=pos_data+((i-1)*obs_nele)
        obsfile.seek((seekpos-1)*8,0)
        data = obsfile.read(data_len*8)
        data_element.xs(i)[:] = numpy.asarray(struct.unpack(">"+str(data_len)+"d", data))
    return(data_element)

def obstore_read_data_station_element(obsfile,subtype,StnNo,element,nmlfile):
    elist = obstore_read_elements(obsfile,int(subtype),nmlfile)
    (indx,pos_data,obs_count,tcols,data_len,data_end)=obstore_read_batchinfo(obsfile,indx)
    record_pos = obstore_read_station_record_position(obsfile,subtype,StnNo,nmlfile)
    return(obstore_read_element(obsfile,elist,element,pos_data,record_pos,obs_nele,nmlfile))
    
def obstore_read_level_position(obsfile,subtype,record,PLEV,nmlfile):
    try:
        lev_pos=numpy.where(numpy.asarray(obstore_read_data_record_element(obsfile,subtype,[record],"PlevelsA",nmlfile).values) == PLEV)[1][0]+1
        return(lev_pos)
    except:errprint("level %s not listed in record %s"%(PLEV,record))
    
def obstore_read_data_record_element_plevel(obsfile,subtype,record_pos,element,plev,nmlfile):
    elist = obstore_read_elements(obsfile,int(subtype),nmlfile)
    (indx,pos_data,obs_count,tcols,data_len,data_end)=obstore_read_batchinfo(obsfile,indx)
    if element in elist.Element.values:
        lev_pos = numpy.zeros(len(record_pos), dtype=numpy.int)
        for i,recpos in enumerate(record_pos):
            lev_pos[i] = obstore_read_level_position(obsfile,subtype,recpos,plev,nmlfile)
        return(obstore_read_element_level(obsfile,subtype,element,pos_data,record_pos,obs_nele,lev_pos,nmlfile))

def obstore_read_data_station_element_plevel(obsfile,subtype,StnNo,element,plev,nmlfile):
    elist = obstore_read_elements(obsfile,int(subtype),nmlfile)
    (indx,pos_data,obs_count,tcols,data_len,data_end)=obstore_read_batchinfo(obsfile,indx)
    if element in elist.Element.values:
        record_pos=obstore_read_station_record_position(obsfile,subtype,StnNo,nmlfile)
        lev_pos = numpy.zeros(len(record_pos), dtype=numpy.int)
        for i,recpos in enumerate(record_pos):
            lev_pos[i] = obstore_read_level_position(obsfile,subtype,recpos,plev,nmlfile)
        return(obstore_read_element_level(obsfile,subtype,element,pos_data,record_pos,obs_nele,lev_pos,nmlfile))

def obstore_read_data_element_plevel(obsfile,subtype,element,plev,nmlfile):
    elist = obstore_read_elements(obsfile,int(subtype),nmlfile)
    (indx,pos_data,obs_count,tcols,data_len,data_end)=obstore_read_batchinfo(obsfile,indx)
    if element in elist.Element.values:
        record_pos=range(1,obs_count+1,1)
        lev_pos = numpy.zeros(len(record_pos), dtype=numpy.int)
        for i,recpos in enumerate(record_pos):
            try:lev_pos[i] = obstore_read_level_position(obsfile,subtype,recpos,plev,nmlfile)
            except:errprint(obstore_read_data_record_element(obsfile,subtype,[recpos],"PlevelsA",nmlfile).values)
        return(obstore_read_element_level(obsfile,subtype,element,pos_data,record_pos,obs_nele,lev_pos,nmlfile))

def obstore_read_element_lvl(obsfile,subtype,element,pos_data,record_pos,obs_nele,lev_pos,nmlfile):
    elist = obstore_read_elements(obsfile,int(subtype),nmlfile)
    (rdc,cdc,ldc) = elist.query("Element == @element")[["RDC","CDC","LDC"]].reset_index(drop=True).values[0]
    ele_pos_offset=int(cdc)-1
    ncols=int(ldc)
    data_element=pandas.DataFrame(index=record_pos,columns=[element])
    for a,i in enumerate(record_pos):
                lev_idx=lev_pos-1
                seekpos=pos_data+((i-1)*obs_nele)+int(cdc-1)+((lev_pos)*rdc)
                obsfile.seek((seekpos-1)*8,0)
                data = obsfile.read(1*8)
                data_element.xs(i)[0] = numpy.asarray(struct.unpack(">"+str(1)+"d", data))[0]
    return(data_element) 
        
def obstore_read_data_record_element_lvl(obsfile,subtype,record_pos,element,lev_pos,nmlfile):
    elist = obstore_read_elements(obsfile,int(subtype),nmlfile)
    (indx,pos_data,obs_count,tcols,data_len,data_end)=obstore_read_batchinfo(obsfile,indx)
    if element in elist.Element.values:
        return(obstore_read_element_lvl(obsfile,subtype,element,pos_data,record_pos,obs_nele,lev_pos,nmlfile))

def obstore_read_data_station_element_lvl(obsfile,subtype,StnNo,element,lev_pos,nmlfile):
    elist = obstore_read_elements(obsfile,int(subtype),nmlfile)
    (indx,pos_data,obs_count,tcols,data_len,data_end)=obstore_read_batchinfo(obsfile,indx)
    if element in elist.Element.values:
        record_pos=obstore_read_station_record_position(obsfile,subtype,StnNo,nmlfile)
        return(obstore_read_element_lvl(obsfile,subtype,element,pos_data,record_pos,obs_nele,lev_pos,nmlfile))

def obstore_read_data_element_lvl(obsfile,subtype,element,lev_pos,nmlfile):
    elist = obstore_read_elements(obsfile,int(subtype),nmlfile)
    (indx,pos_data,obs_count,tcols,data_len,data_end)=obstore_read_batchinfo(obsfile,indx)
    if element in elist.Element.values:
        record_pos=range(1,obs_count+1,1)
        return(obstore_read_element_lvl(obsfile,subtype,element,pos_data,record_pos,obs_nele,lev_pos,nmlfile))

def obstore_copy_data_element(outfile,nmlfile,indx,element,infile,maxindx=MAXINDX):
    if element in obstore_read_batch_elements(infile,indx,nmlfile,maxindx).Element.values:
        data = obstore_read_data_element(infile,nmlfile,indx,element)
        obstore_write_data_element(outfile,nmlfile,indx,element,data)
        return data

def obstore_copy_subtype_index(outfile,nmlfile,pos,subtype,infile):
    (indx,pos_data,obs_count,tcols,data_len)=obstore_read_subtype_index(infile,subtype)
    elist = obstore_read_elements(infile,int(subtype),nmlfile)
    return(obstore_write_subtype_index(outfile,nmlfile,pos,subtype,obs_count,elist))
    
def obstore_copy_batchinfo(outfile,nmlfile,indx,subtype,infile,maxindx=MAXINDX):
    (indx,pos_data,obs_count,tcols,data_len,data_end)=obstore_read_batchinfo(infile,indx)
    elist = obstore_read_batch_elements(infile,indx,nmlfile,maxindx)
    return(obstore_write_subtype_index(outfile,nmlfile,indx,subtype,obs_count,elist))

def print_data_batch(obsfile,nmlfile,subtype,indx,maxindx=MAXINDX):
    if diaglev > 0 : print(subtype)
    elist=obstore_read_batch_elements(obsfile,indx,nmlfile,maxindx)  
    if diaglev > 0 : print(elist)
    elenams=elist.Element.values
    if diaglev > 0 : errprint(elenams)
    for element in elenams:
        if diaglev > 0 : errprint(element)
        print(obstore_read_data_element(obsfile,nmlfile,indx,element))

def getelenams(obsfile,nmlfile,subtype=None,indx=None,maxindx=MAXINDX):
    if subtype is None: subtype=obstore_read_index_subtype(obsfile,indx)
    if indx is None: indx=obstore_read_subtype_index(obsfile,subtype)
    elist=obstore_read_batch_elements(obsfile,indx,nmlfile,maxindx)  
    if diaglev > 0 : print(elist)
    elenams=elist.Element.values
    if diaglev > 0 : errprint(elenams)
    return(elenams)

def getldc(elenam,elist=pandas.DataFrame(),obsfile=None,indx=None,nmlfile=None,maxindx=MAXINDX):
    if elist.empty: elist=obstore_read_batch_elements(obsfile,indx,nmlfile,maxindx)
    element_list=elist.Element.values
    ldc=elist.LDC.values
    print(ldc)
    for i,element in enumerate(element_list,start=0):
        if element in [elenam]:
            print(element,ldc[i])
            return(ldc[i])

def frame_data_batch(obsfile,nmlfile,indx,elenams,maxindx=MAXINDX):
    print("inside frame_data_batch function of obstore")
    print(elenams)
    dataframelist=[None]*len(elenams)
    for i,element in enumerate(elenams):
        if diaglev > 0 : errprint(element)
        dataframelist[i]=obstore_read_data_element(obsfile,nmlfile,indx,element,maxindx=maxindx)
	#print(dataframelist[i])
    data=obslib.obsdfcat(dataframelist)
    print(data)
    return(data)

def erase_data_batch(obsfile,nmlfile,subtype,indx,maxindx=MAXINDX):
    if diaglev > 0 : errprint(subtype)
    elist=obstore_read_batch_elements(obsfile,indx,nmlfile,maxindx)  
    if diaglev > 0 : errprint(elist)
    elenams=elist.Element.values
    if diaglev > 0 : errprint(elenams)
    for element in elenams:
        if diaglev > 0 : errprint(element)
        print(obstore_erase_data_element(obsfile,nmlfile,indx,element))

def obstore_erase_element_backup(obsfile,batchindx,elist,element,pos_data,record_pos,tcols,fillval=NAN_VAL,filedata=None):
    (rdc,cdc,ldc) = elist.query("Element == @element")[["RDC","CDC","LDC"]].reset_index(drop=True).values[0]
    ele_pos_offset=int(cdc)-1
    ncols=int(ldc)
    fmtstr=struct.Struct(binfmt("data",1))
    data=fillval
    if ncols is 1 :
            levs=[1]
    else:
            levs=range(1,ncols+1,1)
    for i in record_pos:
            for j in levs:
                        errprint(data)
                        errprint(fmtstr.pack(data))
                        seekpos=int(pos_data+((i-1)*tcols)+int(cdc-1)+(j-1)*rdc)
                        errprint(seekpos)
                        obsfile.seek((seekpos-1)*8,0)
                        #obsfile.write(fmtstr.pack(data))
	    	    	filedata=file_data_modify([seekpos],filedata=filedata)
                        obstore_check_data_length(obsfile,seekpos)
    return(filedata)

def obstore_write_element_backup(obsfile,batchindx,elist,element,pos_data,record_pos,tcols,data,maxindx=MAXINDX,option=0,fillval=NAN_VAL,filedata=None):
    print("Write_element started")
	#### Need to implement indexing insted of record loop for accelerating the code execution.
    #obstore_erase_element(obsfile,batchindx,elist,element,pos_data,record_pos,tcols,fillval=fillval)
    #print(len(data))
    #print(len(record_pos))
    #print(data,record_pos)
    (rdc,cdc,ldc) = elist.query("Element == @element")[["RDC","CDC","LDC"]].reset_index(drop=True).values[0]
    ele_pos_offset=int(cdc)-1
    ncols=int(ldc)
    fmtstr=struct.Struct(binfmt("data",1))
    if diaglev > 5: 
        if element in ["CharData"]: print(data)
    for i,rec in enumerate(record_pos):
	print(rec)
        if type(data) is int or type(data) is float: # This option is for packing constant values like Year, Month, Day
            seekpos=int(pos_data+((i)*tcols)+int(cdc-1))
            obsfile.seek((seekpos-1)*8,0)
            #obsfile.write(fmtstr.pack(data))
	    filedata=file_data_append(data,seekpos,filedata=filedata)
        else:
            for j,lev in enumerate(data.columns):
                if type(lev) is str:
                    #print(lev,rec,data.xs(rec)[lev])
		    if "_" in lev : j=(int(lev.split("_")[1])-1)
                    data_one=data.xs(rec)[lev]
		    if rec == 0: print(data_one)
                else:
                    data_one=data.xs(rec)[int(lev)]
                if data_one is None:
                    if diaglev > 500: print(data_one, element, "Case A")
                    if diaglev > 0: errprint("%s at record %s, level %s is %s"%(element,rec,lev,data_one))
                elif type(data_one) is str:
                    if diaglev > 500: print(data_one, element, "Case B")
                    if diaglev > 50: errprint(data_one)
                    data_array=obslib.getascii(data_one)
                    seekpos=int(pos_data+((i)*tcols)+int(cdc-1))
                    obsfile.seek((seekpos-1)*8,0)
                    for k,data_one in enumerate(data_array):
                            packed_data=fmtstr.pack(data_one)
                            if diaglev > 50: errprint(k,data_one,packed_data)
                            #obsfile.write(packed_data)
	    		    filedata=file_data_append(data_one,seekpos,filedata=filedata)
                else:
                    if diaglev > 500: print(data_one, element, "Case C")
		    if rdc < 0 : rdc = 0
                    seekpos=int(pos_data+((i)*tcols)+int(cdc-1)+j*rdc)
                    obsfile.seek((seekpos-1)*8,0)
                    #obsfile.write(fmtstr.pack(data_one))
	    	    filedata=file_data_append(data_one,seekpos,filedata=filedata)
                    obstore_check_data_length(obsfile,seekpos)
    print("Write_element finished")
    return(filedata)

def obs_location(inpath,nmlpath,obstype,filevar=None,maxindx=MAXINDX):
    obstypedic=obsdic.obstype[obstype]
    filename=obstypedic["filename"]
    input_file="%s/%s" % (inpath,filename)
    obs_index_nml="obs_index_nml"
    nmlfile="%s/%s" % (nmlpath,obs_index_nml)
    obsgroup=obstypedic["obsgroup"]
    with open(input_file, "rb") as infile:
            obsgroup=obstore_read_header(infile,287,1)
            subtypegroup=obstore_read_subtype(infile)
            batchcount=len(subtypegroup)
            elistgroup=[None]*batchcount
            datagroup=[None]*batchcount
            for indx,subtype in enumerate(obstore_read_subtype(infile)[0:],start=1):
                elist=obstore_read_batch_elements(infile,indx,nmlfile,maxindx)
                elistgroup[indx-1]=elist
                count=0
                for element in obsdic.station_call_list:
                    elemdata = obstore_read_data_element(infile,nmlfile,indx,element)
                    if elemdata is not None:
                        count+=1
                        if count == 1:
                            location=elemdata
                        else:
                            location=location.join(elemdata)
                datagroup[indx-1]=location
    dataset={"obsgroup":obsgroup, "subtype":subtypegroup, "data":datagroup }
    return(dataset)

def obstore_read_file(inpath,obstype,nmlpath=OBSNML,filevar=None,maxindx=MAXINDX):
    keyinfofile=OBSNML+"/keys_"+obstype+".nml"
    infokeylist=["OBSTORE", "obsgroup","maxindx","hdrsize","lutsize"]
    obstypedic=obslib.get_key_dic(keyinfofile,infokeylist,infodic={})
    print(obstypedic)
    filename=obstypedic["OBSTORE"]+".obstore"
    input_file="%s/%s" % (inpath,filename)
    obs_index_nml="obs_index_nml"
    nmlfile="%s/%s" % (nmlpath,obs_index_nml)
    obsgroup=obstypedic["obsgroup"]
    with open(input_file, "rb") as infile:
            obsgroup=obstore_read_header(infile,287,1)
            subtypegroup=obstore_read_subtype(infile)
            batchcount=len(subtypegroup)
            elistgroup=[None]*batchcount
            datagroup=[None]*batchcount
            for indx,subtype in enumerate(obstore_read_subtype(infile)[0:],start=1):
                elist=obstore_read_batch_elements(infile,indx,nmlfile,maxindx)
                elistgroup[indx-1]=elist
                count=0
                elenams=elist.Element.values
                for element in elenams:
                    elemdata = obstore_read_data_element(infile,nmlfile,indx,element)
                    if elemdata is not None:
                        count+=1
                        if count == 1:
                            location=elemdata
                        else:
                            location=location.join(elemdata)
                datagroup[indx-1]=location
    dataset={"obsgroup":obsgroup, "subtype":subtypegroup, "data":datagroup }
    return(dataset)

#############202210#####################################################################################

def assign_subtype(data,nmlfile):
    subtype=obslib.get_key_info(nmlfile,"obsubtyp")
    print(subtype)
    data=data.assign(subtype=[int(subtype)]*len(data))
    return(data)

def uniform_batch(data,btchcnt):
	obscount=len(data)
	batch_obs_cnt=[int(obscount/btchcnt)]*btchcnt
	residue_count=obscount-int(obscount/btchcnt)*btchcnt
	for i in range(1,btchcnt+1):
		if residue_count >= i : batch_obs_cnt[i-1]=batch_obs_cnt[i-1]+1
		print(i, batch_obs_cnt[i-1])
        return(batch_obs_cnt)

def data_batch(data,btch_cnt_list):
	indxzero=0
	datagroup=[]
	for btchid,btch_cnt in enumerate(btch_cnt_list,start=1):
		indxfrst=indxzero
		indxlast=indxzero+btch_cnt
		data1=data.iloc[indxfrst:indxlast]
		datagroup=datagroup+[data1]
		indxzero=indxlast
	return(datagroup)

def dataselect(data,element,ldc):
    #print("Dataselect started")
    if ldc == 1:
        datanew=data[[element]]
    if ldc > 1:
	index=data.index
	datanew=obslib.DataFrame(index=index)
	for i in range(1,ldc+1,1):
	    elenew=element+"_"+str(i)
            if elenew in data.columns:
		datanew[elenew]=data[[elenew]]
            if element in data.columns:
                data1=data[[element]]
                data1.columns = [element+"_"+str(i)]
                if i is 1:
                    datanew=data1
                else:
                    datanew=datanew.join(data1)
            if diaglev > 999 : 
		print(datanew)
    #datanew=obslib.reset_index(datanew,index=data.index.copy())
    #print(data.xs(data.index[0]))
    #print(datanew.xs(datanew.index[0]))
    #print("Dataselect finished")
    return(datanew)

def obstore_file_data_write(obsfile,filedata,textfile=None,datfile=None,fillval=NAN_VAL):
    #print(filedata)
    pos=filedata.index[0]
    data=numpy.array(filedata["data"].values)
    try:data=data.astype(float)
    except ValueError:data=fillval
    if textfile is not None: obslib.obs_frame_ascii(filedata,textfile,0) 
    if datfile is not None :
	print(datfile)
	with file(datfile, "w") as datout:
		fmt0 = '%15.6f'
    	   	numpy.savetxt(datout, X=data, delimiter="\n",fmt=fmt0)
    if obsfile is not None :
    	status=obslib.binary_write(data,pos,obsfile)
    else:
	status=-1
    return(status)
    
def file_data_create(dpos,dlen=1,fmtflg="data",filedata=None,data=NAN_VAL,fillval=NAN_VAL,elenam="",index=None):
    #print("index inside file_data_create")
    #print(index)
    if index is None : index=range(dpos,(dpos+dlen+1))
    if filedata is None : filedata=obslib.DataFrame(columns=["pos","len","data","fmtflg","elenam"],index=index)
    #print("filedata length is "+str(len(filedata)))
    #print(filedata.index)
    filedata=filedata.fillna(fillval)
    filedata["pos"]=filedata.index
    filedata=filedata.set_index("pos")
    #print("filedata length is "+str(len(filedata)))
    #print(filedata.index)
    return(filedata)

def file_data_append(endpos,dpos=1,dlen=1,fmtflg="data",filedata=None,data=NAN_VAL,fillval=NAN_VAL,elenam=""):
    if type(dpos) is list : dpos=dpos[0]
    if filedata is None : filedata=file_data_create(dpos)
    newpos=filedata.index[-1]+1
    newlen=endpos-newpos+1
    newdf=file_data_create(newpos,newlen,fillval=fillval,elenam=elenam)
    filedata = filedata.append(newdf)
    return(filedata)

def check_length(dpos,data,elenam,message=""):
    print(message)
    if type(data) is not int and type(data) is not float:
    	if len(dpos) == len(data):
		print("data and index length are matching for element "+elenam)
    	else:
		print("length missmatch for element "+elenam)
		print(len(dpos))
		print(len(data))

def file_data_modify(dpos,dlen=1,fmtflg="data",filedata=None,data=NAN_VAL,diagflg=0,diagout=None,elenam="",fdindx=None):
    #index1=pandas.Series(dpos,index=range(1,len(dpos)+1))
    index=pandas.Series(dpos,index=range(1,len(dpos)+1))
    #print("index inside file_data_modify")
    #print(index)
    minlen=dpos[-1] - dpos[0] + 1
    #print("Minlen is "+str(minlen))
    if fdindx is None : fdindx=range(dpos[0],dpos[-1])
    if len(filedata) < minlen : filedata=file_data_create(dpos[0],minlen,index=fdindx)
    #print("filedata length is "+str(len(filedata)))
    #print(filedata.index)
    check_length(dpos,data,elenam)
    if type(data) is int or type(data) is float:
	#print("Case A: Constant value data")
	data = [data]*len(index)
	data_new=pandas.DataFrame({"data" :data },index=index)
	data_new.loc[index]["pos"] = index
	data_new["len"] = dlen
	data_new["fmtflg"] = fmtflg
	data_new["elenam"] = elenam
	filedata.loc[index] = numpy.nan
	filedata = filedata.combine_first(data_new.loc[index])
    elif type(data) is list :
	#print("Case B: List of size "+str(len(data)))
	data_new=pandas.DataFrame({"data" :data },index=index)
	data_new.loc[index]["pos"] = index
	data_new["len"] = dlen
	data_new["fmtflg"] = fmtflg
	data_new["elenam"] = elenam
	filedata.loc[index] = numpy.nan
	filedata = filedata.combine_first(data_new.loc[index])
    elif type(data) is not pandas.DataFrame and type(data) is pandas.Series :
	#print("Case C: Series of size "+str(data.size))
	data_new=pandas.DataFrame({"data" :data.values },index=index)
	data_new.loc[index]["pos"] = index
	data_new["len"] = dlen
	data_new["fmtflg"] = fmtflg
	data_new["elenam"] = elenam
	filedata.loc[index] = numpy.nan
	filedata = filedata.combine_first(data_new.loc[index])
    elif len(data.columns) == 1 :
	#print("Case D: Dataframe of "+str(len(data.columns))+" columns.")
	data_new=pandas.DataFrame({"data" : data[data.columns[0]].values },index=index)
	data_new.loc[index]["pos"] = index
	data_new["len"] = dlen
	data_new["fmtflg"] = fmtflg
	data_new["elenam"] = elenam
	filedata.loc[index] = numpy.nan
	filedata = filedata.combine_first(data_new.loc[index])
    else :
	print("Case E: Dataframe of "+str(len(data.columns))+" columns.")
	print(data)
	print(len(data.columns))
	print(filedata)
    if diagflg > 10 : 
		#print(index)
		print(dpos[0],dpos[-1],minlen)
		print(data)
		print(filedata.loc[index])
    return(filedata)

def obstore_erase_element(obsfile,batchindx,elist,element,pos_data,record_pos,tcols,fillval=NAN_VAL,filedata=None,diagflg=0,diagout=None):
    (rdc,cdc,ldc) = elist.query("Element == @element")[["RDC","CDC","LDC"]].reset_index(drop=True).values[0]
    ele_pos_offset=int(cdc)-1
    ncols=int(ldc)
    fmtstr=struct.Struct(binfmt("data",1))
    data=fillval
    if ncols is 1 :
            levs=[1]
    else:
            levs=range(1,ncols+1,1)
    for j,lev in enumerate(levs):
	if rdc < 0 : rdc = 0
	dpos=numpy.array([pos_data+((i)*tcols)+(cdc-1)+j*rdc for i in range(0,len(record_pos))])
	#print(dpos)
	filedata=file_data_modify(dpos,data=fillval,filedata=filedata,diagflg=diagflg,diagout=diagout,elenam=element)
    return(filedata)
                
def obstore_write_element(obsfile,batchindx,elist,element,pos_data,record_pos,tcols,data,maxindx=MAXINDX,option=0,fillval=NAN_VAL,filedata=None,diagflg=0,diagout=None):
    print("Write_element "+element+" to dataframe started")
    (rdc,cdc,ldc) = elist.query("Element == @element")[["RDC","CDC","LDC"]].reset_index(drop=True).values[0]
    ele_pos_offset=int(cdc)-1
    ncols=int(ldc)
    fmtstr=struct.Struct(binfmt("data",1))
    nele=tcols
    nrec=len(record_pos)
    dposend=nele*nrec+pos_data
    fdindx=range(pos_data,dposend)
    minlen=len(fdindx)
    if len(filedata) < minlen : filedata=file_data_create(pos_data,minlen,index=fdindx)
    #print("inside obstore_write_element function")
    #print(filedata.index)
    #print("inside obstore_write_element function")
    if diaglev > 5: 
        if element in ["CharData"]: print(data)
    if type(data) is int or type(data) is float:
	dpos=[pos_data+((i)*tcols)+(cdc-1) for i in range(0,len(record_pos))]
	#print(dpos)
	filedata=file_data_modify(dpos,data=data,filedata=filedata,diagflg=diagflg,diagout=diagout,elenam=element,fdindx=fdindx)
    else:
	for j,lev in enumerate(data.columns):
                if type(lev) is str:
		    if "_" in lev : j=(int(lev.split("_")[1])-1)
                    data_one=data[lev]
                else:
                    data_one=data[int(lev)]
                if data_one is None:
                    if diaglev > 500: print(data_one, element, "Case A")
                    if diaglev > 10: errprint("%s at record %s, level %s is %s"%(element,rec,lev,data_one))
                elif type(data_one) is str:
                    if diaglev > 500: print(data_one, element, "Case B")
                    if diaglev > 50: errprint(data_one)
                    data_array=obslib.getascii(data_one)
                    for k,data_one in enumerate(data_array):
                        packed_data=fmtstr.pack(data_one)
                        if diaglev > 50: errprint(k,data_one,packed_data)
			dpos=[pos_data+((i)*tcols)+(cdc-1) for i in range(0,len(record_pos))]
			filedata=file_data_modify(dpos,data=data_one,filedata=filedata,diagflg=diagflg,diagout=diagout,elenam=element,fdindx=fdindx)
                else:
                    	if diaglev > 500: print(data_one, element, "Case C")
		    	if rdc < 0 : rdc = 0
			dpos=numpy.array([pos_data+((i)*tcols)+(cdc-1)+j*rdc for i in range(0,len(record_pos))])
			filedata=file_data_modify(dpos,data=data_one,filedata=filedata,diagflg=diagflg,diagout=diagout,elenam=element,fdindx=fdindx)
    #print("Write_element to dataframe finished")
    return(filedata)

def obstore_write_data(obsfile,nmlfile,elist,data,batchindx=1,maxindx=MAXINDX,filedata=None,diagflg=0,diagout=None):
    elenams=elist.Element.values
    ldc=elist.LDC.values
    (indx,pos_data,obs_count,tcols,data_len,data_end)=obstore_read_batchinfo(obsfile,batchindx)
    if filedata is None : filedata=file_data_create(pos_data,data_len)
    elem_name_list=list(data.columns[:])
    for i in range(len(elem_name_list)):
	elem_name_list[i]=elem_name_list[i].split("_")[0]
    #print(elem_name_list[0])
    #print(elenams[0])
    diagout_name=diagout
    for idx,element in enumerate(elenams):
    	if diagflg > 0 : 
		diagout=diagout_name+"_"+str(batchindx)+"_"+str(idx+1)
        if element in elem_name_list:
            print(element+" is available")
            record_pos=data.index
            #range(1,len(data)+1,1)
            eledata=dataselect(data,element,ldc[idx])
            #print(ldc[idx], eledata, element)
    	    if diagflg > 99 : 
		print(eledata)
		obslib.obs_frame_ascii(eledata,diagout,0) 
            filedata=obstore_write_element(obsfile,batchindx,elist,element,pos_data,record_pos,tcols,eledata,maxindx=maxindx,filedata=filedata,diagflg=diagflg,diagout=diagout)
        else:
            filedata=obstore_erase_element(obsfile,batchindx,elist,element,pos_data,record_pos,tcols,fillval=NAN_VAL,filedata=filedata,diagflg=diagflg,diagout=diagout)
            print(element+" is not available")
    #print(filedata)
    return(filedata)

def obstore_write_batches(DT,obsfile,nmlfile,obsgroup,subtypegroup,elistgroup,datagroup,batchcount=1,hdrsize=HDRSIZE,maxindx=MAXINDX,lutsize=LUTSIZE,filedata=None,diagflg=0,diagout=None,HlfTW=HLFTW,Tref=TREF,callsignflag=False):
    #hdrsize=obsheader.write_obsheader(obsfile,nmlfile,obsgroup,maxindx=maxindx,callsignflag=callsignflag)
    print("File Header size is "+str(hdrsize))
    datapos=obstore_set_batchpos(obsfile,batchcount,batch_data_offset=0,hdrsize=hdrsize,maxindx=maxindx,lutsize=lutsize)    
    print("Data start position is "+str(datapos))
    datalen=0
    if filedata is None : filedata=file_data_create(datapos,datalen)
    for idx in range(0,batchcount,1):
        subtype=subtypegroup[idx]
        elist=elistgroup[idx]
        data=datagroup[idx]
        if data is not None:
            obscount=len(data)
	    batch_data_offset=datalen
            datalen=write_batchheader(DT,obsfile,nmlfile,subtype,elist,batchid=idx+1,batch_data_offset=batch_data_offset,obscount=obscount,batchcount=batchcount,hdrsize=hdrsize,maxindx=maxindx,lutsize=lutsize,HlfTW=HlfTW,Tref=Tref)
	    print(idx,obscount,batch_data_offset,datalen)
	    dataendpos=datapos+batch_data_offset+datalen
	    print(len(filedata.index),datalen)
	    #if len(filedata.index) < datalen : filedata=file_data_append(dataendpos,filedata=filedata)
            filedata=obstore_write_data(obsfile,nmlfile,elist,data,batchindx=(idx+1),maxindx=maxindx,filedata=filedata,diagflg=diagflg,diagout=diagout)
    obstore_write_obsgroup(obsgroup,obsfile)
    (datapos,datalen,dataend)=obstore_close_data(obsfile)
    print(datapos,datalen,dataend)
    #print(filedata)
    return(filedata)

def totobs(obstore_info):
    	outpath = obstore_info["outpath"]
    	batch_obs_cnt = obstore_info["count_list"]
	#############################################################################
	with file(outpath+"/totobs.dat", "w") as datout:
		fmt01 = '%15.0i'
    	   	numpy.savetxt(datout, X=batch_obs_cnt, delimiter="\n",fmt=fmt01)
	obscnt=numpy.sum(batch_obs_cnt)
	elist=elistgroup[0]
	elecnt=numpy.sum(elist.LDC)
	datlen=obscnt*elecnt
	with file(outpath+"/cnt.dat", "w") as datout:
		fmt01 = '%15.0i'
    	   	numpy.savetxt(datout, X=[datlen,elecnt], delimiter="\n",fmt=fmt01)
	############################################################################
	datapos=obstore_set_batchpos(obsfile,batchcount,batch_data_offset=0,hdrsize=hdrsize,maxindx=maxindx,lutsize=lutsize)  
	if filedata is None : filedata=file_data_create(datapos,datlen-1)
	return(filedata)
    
def obstore_create_file(obstore_info,diagflg=0,callsignflag=False,filedata=None,):
    obstore_info=obstore_info_fix(obstore_info)
    obstype = obstore_info["obstype"]
    obsgroup = int(obstore_info["obsgroup"])
    subtypegroup = obstore_info["subtypegroup"]
    elistgroup = obstore_info["elistgroup"]
    datagroup = obstore_info["datagroup"]
    DT = obstore_info["timeinfo"]
    if "timewindowhalf" in obstore_info: 
	HlfTW=obstore_info["timewindowhalf"]
    else:
	HlfTW=None
    if "timereference" in obstore_info: 
	Tref=obstore_info["timereference"]
    else:
	Tref=None
    outpath = obstore_info["outpath"]
    filename = obstore_info["filename"]
    batchcount = int(obstore_info["batchcount"])
    if "nmlpath" in obstore_info:
	nmlpath = obstore_info["nmlpath"]
    else:
	nmlpath = OBSNML
    if "maxindx" in obstore_info:
	maxindx = int(obstore_info["maxindx"])
    else:
	maxindx = MAXINDX
    if "hdrsize" in obstore_info:
	hdrsize = int(obstore_info["hdrsize"])
    else:
	hdrsize = HDRSIZE
    if "lutsize" in obstore_info:
	lutsize = int(obstore_info["lutsize"])
    else:
	lutsize = LUTSIZE
    datfile="%s/%s" % (outpath,filename.replace(".","_")+".dat")
    output_file="%s/%s" % (outpath,filename)
    obs_index_nml="obs_index_nml"
    nmlfile="%s/%s" % (nmlpath,obs_index_nml)
    if diagflg > 0:
    	textfile=outpath+"/"+filename.replace(".","_")+".txt"
    	diagout=outpath+"/"+filename.replace(".","_")+"_diag"
    else :
	textfile=None
	diagout=None
    print("Writting "+str(batchcount)+" batches of data to "+ output_file)
    obslib.mkdir(output_file.rsplit("/",1)[0])
    if "count_list" in obstore_info:
    	filedata=totobs(obstore_info)
    with open(output_file, "wb+") as obsfile:
	hdrsize=obsheader.write_obsheader(obsfile,nmlfile,obsgroup,maxindx=maxindx,callsignflag=callsignflag)
    	#obslib.binary_write(FIXHDR,1,obsfile)
    	#obslib.binary_write(HDR20,1,obsfile)
    	#obslib.binary_write(HDRgam,306,obsfile)
	#############################################################################
	#with file(outpath+"/totobs.dat", "w") as datout:
	#	fmt01 = '%15.0i'
    	#   	numpy.savetxt(datout, X=batch_obs_cnt, delimiter="\n",fmt=fmt01)
	#obscnt=numpy.sum(batch_obs_cnt)
	#elist=elistgroup[0]
	#elecnt=numpy.sum(elist.LDC)
	#datlen=obscnt*elecnt
	#with file(outpath+"/cnt.dat", "w") as datout:
	#	fmt01 = '%15.0i'
    	#   	numpy.savetxt(datout, X=[datlen,elecnt], delimiter="\n",fmt=fmt01)
	#datapos=obstore_set_batchpos(obsfile,batchcount,batch_data_offset=0,hdrsize=hdrsize,maxindx=maxindx,lutsize=lutsize) 
	#if filedata is None : filedata=file_data_create(datapos,datlen-1)
	############################################################################
        filedata=obstore_write_batches(DT,obsfile,nmlfile,obsgroup,subtypegroup,elistgroup,datagroup,batchcount=batchcount,hdrsize=hdrsize,maxindx=maxindx,lutsize=lutsize,filedata=filedata,diagflg=diagflg,diagout=diagout,HlfTW=HlfTW,Tref=Tref,callsignflag=callsignflag)
	statusflg=obstore_file_data_write(obsfile,filedata,textfile=textfile,datfile=datfile)
        #print(filedata)
    print("Writting to "+output_file+ " is completed")
    return(obstore_info)

def obstore_write(data,keyinfofile,outpath,btchcnt=None,obstore_info=None,cntmax=None,DT=None,Tref=None,HlfTW=None,diagflag=0,callsignflag=False,missing_value=-1073741824.00000):
	data = data[data[["Latitude", "Longitude"]].notnull().all(1)]
	data = data.replace(numpy.nan,missing_value)
	print(data)
	infokeylist=["obsgroup","maxindx","hdrsize","lutsize"]
	obstore_info=obslib.get_key_dic(keyinfofile,infokeylist,obstore_info)
	obstore_info["obstype"]=obslib.get_key_info(keyinfofile,"obsname")
	obstore_info["maxindx"]=int(obslib.get_key_info(keyinfofile,"maxindx"))
	outfile=obslib.get_key_info(keyinfofile,"OBSTORE")
	obstore_info["outpath"]=outpath
	obstore_info["filename"]=outfile+".obstore"
	if "nmlfile" not in obstore_info: obstore_info["nmlfile"]=NML_OBS_INDX
	if "subtype" not in data: data=assign_subtype(data,keyinfofile)
	if cntmax is not None : data=obslib.data_thinning(data,cntmax=cntmax,callsign="SatID",fillval=NAN_VAL)
	#print(list(data.columns))
	subtype_list=data.subtype.unique()
	#print(subtype_list)
	if btchcnt is None : btchcnt=len(subtype_list)
	if len(subtype_list) > 1 : btchcnt=len(subtype_list)
	if btchcnt > len(subtype_list) or btchcnt == 1 :
		batch_obs_cnt=uniform_batch(data,btchcnt)
		elistkey="elemlist_"+str(subtype_list[0])
		print(elistkey)
		elemlist=obslib.get_key_list_info(keyinfofile,elistkey)
		nmlfile=obstore_info["nmlfile"]
		elist=obstore_create_element_table(nmlfile,elemlist)
		print(elist)
	else:
	    batch_obs_cnt=[]
	    for btchindx,subtype in enumerate(subtype_list,start=1):
		#print(btchindx,subtype)
		data1=data[data.subtype==subtype]
		batch_obs_cnt=batch_obs_cnt+[len(data1)]
		elistkey="elemlist_"+str(subtype)
		elemlist=obslib.get_key_list_info(keyinfofile,elistkey)
		nmlfile=obstore_info["nmlfile"]
		elist=obstore_create_element_table(nmlfile,elemlist)
		print(elist)
	#print(batch_obs_cnt)
	datagroup=data_batch(data,batch_obs_cnt)
	btchcnt=len(batch_obs_cnt)
	obstore_info["batchcount"] = btchcnt
	obstore_info["count_list"] = batch_obs_cnt
        obstore_info["datagroup"] = datagroup
	subtypegroup=[subtype_list[0]]*btchcnt
	obstore_info["subtypegroup"]=subtypegroup
	obstore_info["elistgroup"]=[elist]*btchcnt
	if DT is not None :
		obstore_info["timeinfo"]=DT
	else :
		obstore_info["timeinfo"]=obslib.get_date_info(data)
	if Tref is not None :
		obstore_info["timereference"]=Tref
	else :
		obstore_info["timereference"]=TREF
	if HlfTW is not None :
		obstore_info["timewindowhalf"] = HlfTW
	else :
		obstore_info["timewindowhalf"] = HLFTW
	obstore_info=obstore_create_file(obstore_info,diagflag,callsignflag=callsignflag)
	return(obstore_info)

#############202307#####################################################################################
def obs_hdr_read(inputfile,maxindx=int(MAXINDX),batchid=1,debug=False):
	batchid=int(batchid)
	if "fixhdr" in inputfile:
		ftype="fixhdr"
	else: ftype=inputfile.split(".")[1]

	print(maxindx)
	with open(inputfile, "rb") as obsfile:
	    hdr_info=obstore_read_batch_header(obsfile,maxindx=maxindx,ftype=ftype,batchid=batchid)
	if not ftype in ["fixhdr"]:
		if debug: 
			for i,val in enumerate(hdr_info["cdc"],1): 
				print(i,val)
		hdr_info["elist"]=frame_batch_elist(hdr_info)
	hdrinfo={}
	for itm in ["alpha","beeta","gamma","elist","lut"]:
		if itm in hdr_info: hdrinfo[itm]=hdr_info[itm]
	return(hdrinfo)

def read_data_record(inputfile,maxindx=MAXINDX,batchid=1,recno=1):
	hdr_info=obs_hdr_read(inputfile,maxindx=maxindx,batchid=batchid)
	elemcnt=hdr_info["elist"].LDC.sum(axis=0)
	recpos=(recno-1)*elemcnt
	recend=(recno)*elemcnt
	data_rec=read_data(inputfile)[recpos:recend]
	return(data_rec)

def read_binfile_strip(inputfile,pos=1,size=1,fmtkey="q"):
    with open(inputfile, "rb") as obsfile:
	datastrip=obstore_read_header(obsfile,pos,size,fmtkey=fmtkey)
    return(datastrip)

