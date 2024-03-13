#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 21 15:45:39 2018

@author: gibies
"""
from __future__ import print_function
import sys
import os
diaglev=int(os.environ.get('GEN_MODE',0))
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
PKGHOME=os.path.dirname(CURR_PATH)
OBSLIB=os.environ.get('OBSLIB',PKGHOME+"/pylib")
sys.path.append(OBSLIB)
OBSDIC=os.environ.get('OBSDIC',PKGHOME+"/pydic")
sys.path.append(OBSDIC)
OBSNML=os.environ.get('OBSNML',PKGHOME+"/nml")
sys.path.append(OBSNML)
varobs_nml=OBSNML+"/varobs_nml"
import obslib
import obsdic
import obsmod
import obsheader
import numpy
import pandas
import struct
import datetime
import itertools

def errprint(*args, **kwargs):
    if diaglev > 0: print(*args, file=sys.stderr, **kwargs)

def header_info(obsfile):
    hbet = obslib.binary_file_segment_read(obsfile,"hbet")
    varobs_header={
    "gridcnt_x"  : hbet[5],
    "gridcnt_y"  : hbet[6],
    "gridcnt_z"  : hbet[7],
    "rec_len"    : hbet[28],
    "obs_grp"    : hbet[30],
    "cnt_meta"   : hbet[39],
    "cnt_data"   : hbet[40],
    "cnt_var"    : hbet[41],
    "cnt_lev"    : hbet[46],
    "cnt_batch"  : hbet[48],
    }
    return(varobs_header)

def read_record(obsfile,batchnum=1,recnum=1):
    (halp,hbet,hgam,ldc,rdc,cdc,lut,necklace)=obslib.binary_read_data_header(obsfile)
    print(hgam)
    hdr=header_info(obsfile)
    obslib.dic_print(hdr)
    #for (key,val) in hdr.items():
    #   print(key+" : "+str(val))
    data_rec_len=hdr["rec_len"]
    datapos=halp[159]
    dataoffset=data_rec_len*(recnum-1)
    recpos=datapos+dataoffset
    data_rec=obslib.binary_read_data(obsfile,recpos,data_rec_len)
    return(data_rec)

def rec_meta(obsfile,batchnum=1,recnum=1):
    hdr=header_info(obsfile)
    data=read_record(obsfile,batchnum,recnum)
    rec_meta=data[0:(hdr["cnt_meta"])]
    time=int(rec_meta[0])
    subtype=int(rec_meta[1])
    print("Time :"+str(time))
    print("Subtype :"+str(subtype))
    return(rec_meta)

def batch_info(obsfile,bthnum=1):
    bthptr=(bthnum-1)
    lut = obslib.binary_file_segment_read(obsfile,"lut")
    bth_hdr = {
    "bth_size"  : lut[bthptr,14],
    "bth_start" : lut[bthptr,28],
    "bth_i30"   : lut[bthptr,29],
    "bth_i39"   : lut[bthptr,38],
    "bth_i40"   : lut[bthptr,39],
    "bth_i65"   : lut[bthptr,64],
    "bth_cnt_obs" : lut[bthptr,65],
    "bth_i67"   : lut[bthptr,66],
    }
    return(bth_hdr)

def get_data(obsfile,bthnum=1,obsnum=1,varnum=1,colnum=1,levnum=1):
    hdr=header_info(obsfile)
    bhdr=batch_info(obsfile,bthnum)
    datapos=obslib.get_datapos(obsfile)

    obslib.dic_print(bhdr)
    batchpos=(bhdr["bth_start"]+1) 
    obspos=batchpos+hdr["rec_len"]*(obsnum-1)
    
    var_width=(hdr["cnt_data"]*hdr["cnt_lev"])
    varpos=obspos+hdr["cnt_meta"]+var_width*(varnum-1)
    colpos=varpos+hdr["cnt_lev"]*(colnum-1)
    levpos=colpos+(levnum-1)
    print("Read position : "+str(levpos))
    data=obslib.binary_read_data(obsfile,levpos,1)

    bthtailpos=(bhdr["bth_start"]+bhdr["bth_size"]+1)
    bthtaillen=(bhdr["bth_i30"]-bhdr["bth_size"])
    tail=obslib.binary_read_data(obsfile,bthtailpos,bthtaillen)
    print(bthtailpos,bthtaillen)
    print(tail[:])
    return(tail)

def read_batch_tail(obsfile,bthnum=1):
    hdr=header_info(obsfile)
    bhdr=batch_info(obsfile,bthnum)
    obslib.dic_print(bhdr)

    print(bhdr["bth_size"])
    tailpos=(bhdr["bth_start"]+bhdr["bth_size"]+1)
    print(tailpos)
    batchend=(bhdr["bth_start"]+bhdr["bth_i30"])
    print(batchend)
    taillen=(bhdr["bth_cnt_obs"]*bhdr["bth_i67"])
    print(taillen)
    data=obslib.binary_read_data(obsfile,tailpos,taillen)
    data=numpy.reshape(data,[bhdr["bth_cnt_obs"],bhdr["bth_i67"]])
    return(data)


def get_varlist(obsfile,nmlfile=varobs_nml):
    cdc_var_pos=obslib.read_cdc(obsfile,4)
    cdc_num_lev=obslib.read_cdc(obsfile,8)
    print(obslib.binary_file_segment_read(obsfile,"cdc"))
    elist=element_frame(cdc_var_pos,cdc_num_lev,nmlfile=nmlfile)
    return(elist)


def element_frame(cdc,ldc,nmlfile=varobs_nml):
    cdc_sortlist=cdc
    cdc_sortlist=obslib.binsort(cdc_sortlist)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=-1.07374182e+09)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=-3.2768e+04)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=-1073741820.0)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=-32768)
    cdc_sortlist=obslib.binsort(cdc_sortlist,missing=0.)
    for idx in cdc_sortlist:
       with open(nmlfile, "r") as nml:
          obs_index=(numpy.where(cdc == idx)[0][0])+1
          CDC=cdc[obs_index-1]
          LDC=ldc[obs_index-1]
          Element=pandas.read_table(nml, skiprows=None, header=0).query("eleindex == @obs_index").elename.reset_index(drop=True)[0]
          if idx == cdc_sortlist[0] : 
             eframe=pandas.DataFrame([[CDC,LDC,obs_index,Element]],index=[idx],columns=["CDC","LDC","obs_index","Element"])
          else : 
             eframe=eframe.append(pandas.DataFrame([[CDC,LDC,obs_index,Element]],index=[idx],columns=["CDC","LDC","obs_index","Element"]))
    return(eframe)

def get_varobs_batch_cnt(obsfile):
    batch_cnt=int(obslib.binary_read_record(obsfile,152,1,1,"q"))
    return(batch_cnt)

def get_varobs_obs_cnt(obsfile):
    obs_cnt=int(obslib.binary_read_record(obsfile,284,1,1,"q"))
    return(obs_cnt)

def get_batch_pos(obsfile,batchid=1):
    lut=obslib.read_lut(obsfile,recno=batchid)
    batchpos=lut[28]+1
    return(batchpos)

def get_batch_len(obsfile,batchid=1):
    lut=obslib.read_lut(obsfile,recno=batchid)
    return(lut[29])

def get_batch_data_len(obsfile,batchid=1):
    lut=obslib.read_lut(obsfile,recno=batchid)
    return(lut[14])

def get_batch_obs_cnt(obsfile,batchid=1):
    lut=obslib.read_lut(obsfile,recno=batchid)
    return(lut[65])

def get_latlon(obsfile,batchid=1):
    batchcnt=get_varobs_batch_cnt(obsfile)
    latlonstrt=get_batch_pos(obsfile,batchcnt-1)
    for idx in range(1,batchid):
	print(idx)
        latlonstrt=latlonstrt+2*get_batch_obs_cnt(obsfile,idx)
    batch_obs_cnt=get_batch_obs_cnt(obsfile,batchid)
    latlon=obslib.binary_read_record(obsfile,latlonstrt,batch_obs_cnt,2,"d")
    #obsmod.ascii_file_write(latlon,option=1)
    return(latlon)

def get_time(obsfile,batchid=1):
    batchcnt=get_varobs_batch_cnt(obsfile)
    timestrt=get_batch_pos(obsfile,batchcnt)
    for idx in range(1,batchid):
	print(idx)
        timestrt=timestrt+get_batch_obs_cnt(obsfile,idx)
    batch_obs_cnt=get_batch_obs_cnt(obsfile,batchid)
    time_arr=obslib.binary_read_record(obsfile,timestrt,batch_obs_cnt,1,"d")
    #obsmod.ascii_file_write(time_arr,option=1)
    return(time_arr)

def get_data(obsfile,bthnum=1):
    hdr=header_info(obsfile)
    bhdr=batch_info(obsfile,bthnum)
    datapos=obslib.get_datapos(obsfile)
    batchpos=(bhdr["bth_start"]+1)
    batchend=(bhdr["bth_start"]+bhdr["bth_i30"])
    batchlen=(batchend-batchpos)
    data=obslib.binary_read_data(obsfile,batchpos,batchlen)
    data=numpy.reshape(data[0:bhdr["bth_size"]],[bhdr["bth_cnt_obs"],hdr["rec_len"]])
    #obsmod.ascii_file_write(data,option=1)
    return(data)

def read_batch(obsfile,batchid=1):
    elist=get_varlist(obsfile)
    obsmod.ascii_file_write(elist)
    header=["time_sec","lat","lon","time_min","subtype"]
    header=header+["CS01","CS02","CS03","CS04","CS05","CS06","CS07","CS08","CS09","CS10","CS11","CS12","CS13","CS14","CS15","CS16"]
    header=header+["plev","flag"]
    for idx in range(0,len(elist)):
        varid=elist.Element.values[idx]
        numlev=int(elist.LDC.values[idx])
        for lev in range(1,numlev+1):
            header=header+["obs_val_"+varid+"_"+str(lev)]
        for lev in range(1,numlev+1):
            header=header+["obs_err_"+varid+"_"+str(lev)]
        for lev in range(1,numlev+1):
            header=header+["pge_"+varid+"_"+str(lev)]
    data=get_data(obsfile,batchid)
    latlon=get_latlon(obsfile,batchid)
    obstim=get_time(obsfile,batchid)
    batch_obs_cnt=get_batch_obs_cnt(obsfile,batchid)
    for idx in range(0,batch_obs_cnt):
        obs_rec=numpy.concatenate((obstim[idx],latlon[idx],data[idx]))
        if idx == 0 :
           df=pandas.DataFrame([obs_rec],index=[idx+1],columns=header)
        else :
           df=df.append(pandas.DataFrame([obs_rec],index=[idx+1],columns=header))
    data=df
    obsmod.ascii_file_write(data,option=1)
    return(data)

def read_latlon(infilename,btchid=1): 
    with open(infilename, "rb") as obsfile:
	a=get_latlon(obsfile,batchid=btchid)
