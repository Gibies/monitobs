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
import obslib
import obsdic
import obsmod
import obsheader
import numpy
import pandas
import struct
import datetime
import itertools


def header_info(obsfile):
    hbet = obslib.binary_file_segment_read(obsfile,"hbet")
    varobs_header={
    "gridcnt_x"  : hbet[5],
    "gridcnt_y"  : hbet[6],
    "plev_cnt"  : hbet[7],
    "wet_lev_cnt": hbet[8],
    "_hi10_"	 : hbet[9],
    "_hi12_"	 : hbet[11],
    "_hi13_"	 : hbet[12],
    "_hi17_"	 : hbet[16],
    "_hi24_"	 : hbet[23],
    "_hi25_"	 : hbet[24],
    "tot_obs_cnt": hbet[27],
    "rec_len"    : hbet[29],
    ###rec_len index is different for varobs and varcx
    "obs_grp"    : hbet[30],
    "cnt_meta"   : hbet[39],
    "cnt_data"   : hbet[40],
    "cnt_var"    : hbet[41],
    "cnt_lev"    : hbet[47],
    ###cnt_lev intex is different for varobs and  varcx
    "cnt_batch"  : hbet[48],
    }
    return(varobs_header)

def batch_info(obsfile,bthnum=1):
    bthptr=(bthnum-1)
    lut = obslib.binary_file_segment_read(obsfile,"lut")
    bth_hdr = {
    "bth_data_len" : lut[bthptr,14],
    "bth_strt_pos" : lut[bthptr,28],
    "bth_tot_len"  : lut[bthptr,29],
    "bth_i39"	   : lut[bthptr,38],
    "bth_obs_cnt"  : lut[bthptr,65],
    }
    return(bth_hdr)

def element_frame(cdc,ldc,nmltype):
    nmlfile=obsmod.get_nml(nmltype)
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
          if ldc is not None :
             LDC=ldc[obs_index-1]
          else :
             LDC=1
          Element=pandas.read_table(nml, skiprows=None, header=0).query("eleindex == @obs_index").elename.reset_index(drop=True)[0]
          if idx == cdc_sortlist[0] : 
             eframe=pandas.DataFrame([[CDC,LDC,obs_index,Element]],index=[idx],columns=["CDC","LDC","obs_index","Element"])
          else : 
             eframe=eframe.append(pandas.DataFrame([[CDC,LDC,obs_index,Element]],index=[idx],columns=["CDC","LDC","obs_index","Element"]))
    return(eframe)

def get_varlist(obsfile):
    nmltype="cx_surf"
    cdc_var_pos=obslib.read_cdc(obsfile,6)
    cdc_lev_cnt=None
    elist=element_frame(cdc_var_pos,cdc_lev_cnt,nmltype=nmltype)
    nmltype="cx_uair"
    cdc_var_pos=obslib.read_cdc(obsfile,7)
    cdc_lev_cnt=obslib.read_cdc(obsfile,8)
    elist=elist.append(element_frame(cdc_var_pos,cdc_lev_cnt,nmltype=nmltype))
    obsmod.ascii_file_write(elist)
    return(elist)

def get_data(obsfile,bthnum=1):
    hdr=header_info(obsfile)
    bhdr=batch_info(obsfile,bthnum)
    obsmod.ascii_file_write(bhdr,option=1)
    datapos=obslib.get_datapos(obsfile)
    batchpos=(bhdr["bth_strt_pos"]+1)
    batchend=(bhdr["bth_strt_pos"]+bhdr["bth_tot_len"])
    batchlen=(batchend-batchpos)
    data=obslib.binary_read_data(obsfile,batchpos,batchlen)
    data=numpy.reshape(data[0:bhdr["bth_data_len"]],[bhdr["bth_obs_cnt"],hdr["rec_len"]])
    #obsmod.ascii_file_write(data,option=1)
    return(data)

def get_latlon(obsfile,batchid=1):
    batchcnt=header_info(obsfile)["cnt_batch"]
    latlonstrt=batch_info(obsfile,batchcnt)["bth_strt_pos"]+1
    for idx in range(1,batchid):
	print(idx)
        latlonstrt=latlonstrt+2*get_batch_obs_cnt(obsfile,idx)
    batch_obs_cnt=batch_info(obsfile,batchid)["bth_obs_cnt"]
    latlon=obslib.binary_read_record(obsfile,latlonstrt,batch_obs_cnt,2,"d")
    #obsmod.ascii_file_write(latlon,option=1)
    return(latlon)

def read_batch(obsfile,batchid=1):
    elist=get_varlist(obsfile)
    header=[]
    for idx in range(0,len(elist)):
        varid=elist.Element.values[idx]
        numlev=int(elist.LDC.values[idx])
        for lev in range(1,numlev+1):
            header=header+["cx_val_"+varid+"_"+str(lev)]
    data=get_data(obsfile,batchid)
    #latlon=get_latlon(obsfile,batchid)
    #obstim=get_time(obsfile,batchid)
    batch_obs_cnt=batch_info(obsfile,batchid)["bth_obs_cnt"]
    print(batch_obs_cnt)
    for idx in range(0,batch_obs_cnt):
        obs_rec=data[idx]
        if idx == 0 :
           df=pandas.DataFrame([obs_rec],index=[idx+1],columns=header)
        else :
           df=df.append(pandas.DataFrame([obs_rec],index=[idx+1],columns=header))
    data=df
    obsmod.ascii_file_write(data,option=1)
    return(data)

