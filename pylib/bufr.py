#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 22 00:24:43 2020

@author: gibies
"""
from __future__ import print_function
import sys
import os
import pandas
import numpy
import string
import re
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
PKGHOME=os.path.dirname(CURR_PATH)
OBSLIB=os.environ.get('OBSLIB',PKGHOME+"/pylib")
sys.path.append(OBSLIB)
OBSDIC=os.environ.get('OBSDIC',PKGHOME+"/pydic")
sys.path.append(OBSDIC)
OBSNML=os.environ.get('OBSNML',PKGHOME+"/nml")
sys.path.append(OBSNML)
import bufrdic
import obstore
import obsmod
import obsdic

def bufr_element_frame(idx,bufr_pos,obs_index,nmlfile):
    with open(nmlfile, "r") as nml:
        Element= pandas.read_table(nml, skiprows=None, header=0).query("eleindex == @obs_index").elename.reset_index(drop=True)[0]
    eframe=pandas.DataFrame([[bufr_pos,obs_index,Element]],index=[idx],columns=["bufr_ele_pos","obs_index","Element"])
    return(eframe)

def bufr_create_element_table(nmlfile,obsidxlist):
    idx=0
    elenamlst=[]
    for i,obs_index in enumerate(obsidxlist,1):
        if obs_index not in elenamlst:
            idx+=1
            if obs_index is not None:
                bufr_pos=i
                elenamlst.append(obs_index)
                if idx == 1: 
                    elist=bufr_element_frame(idx,bufr_pos,obs_index,nmlfile)
                else:
                    elist=elist.append(bufr_element_frame(idx,bufr_pos,obs_index,nmlfile))
    return(elist)

def record_count(TEXT_FILE,HEADER_OFFSET,RECORD_LENGTH):
    num_lines = sum(1 for line in open(TEXT_FILE))
    return(int((num_lines-HEADER_OFFSET)/RECORD_LENGTH))

def read_bufr_dump(INFILE_LIST,HEADER_OFFSET,RECORD_LENGTH,bufr_ele_list):
    nmlfile=obsmod.obs_nml
    bufr_elist=bufr_create_element_table(nmlfile,bufr_ele_list)
    data=pandas.DataFrame(columns=bufr_elist.Element.values)
    batch_count=len(INFILE_LIST)
    datagroup=[None]*batch_count
    batch_id = 1
    for TEXT_FILE in INFILE_LIST:
        print(TEXT_FILE)
        rec_on_batch=record_count(TEXT_FILE,HEADER_OFFSET,RECORD_LENGTH)
        print(rec_on_batch)
        rec_count=0
        ele_count=-1
        value=""
        batch_offset=HEADER_OFFSET
        batch_size=RECORD_LENGTH*rec_on_batch
        with open(TEXT_FILE, "r") as fileptr:
            for i, line in enumerate(fileptr):
                if rec_count > rec_on_batch:
                    print(rec_count, rec_on_batch, rec_on_batch*RECORD_LENGTH)
                    print(batch_offset, i, rec_count, ele_count)
                    if data is not None:
                     	datagroup[batch_id-1]=data
                    batch_id+=1
                    batch_offset=i
                if batch_id > batch_count:
                    break
                if i == batch_offset:
                    rec_count=1
                    data=pandas.DataFrame(columns=bufr_elist.Element.values)
                    val_list=[]
                if rec_count >= 1 and rec_count <= rec_on_batch:
                    word=re.sub('['+string.punctuation+']','',line).split()
                    if len(word) >= 3 and word[0].isdigit():
                        ele_count=int(word[0])
                        strval=line[40:65]
                        try:
                            value=float(strval)
                        except ValueError:
                            value=numpy.nan
                    else:
                        ele_count = 0
                    if ele_count == 0:
                        if len(val_list) == len(bufr_elist.obs_index.values):
                            data.loc[rec_count]=val_list
                        val_list=[]
                        rec_count+=1
                    if ele_count > 0:
                        index=bufr_ele_list[int(ele_count)-1]
                        if index in bufr_elist.obs_index.values:
                            val_list.append(value)
    return(datagroup)

def write_to_obstore(INFILE_LIST,HEADER_OFFSET,RECORD_LENGTH,bufr_ele_list,outpath,nmlpath,DT,obstype,maxindx=608):
    batch_count=len(INFILE_LIST)
    obstypedic=obsdic.obstype[obstype]
    filename=obstypedic["filename"]
    output_file="%s/%s" % (outpath,filename)
    obs_index_nml="obs_index_nml"
    nmlfile="%s/%s" % (nmlpath,obs_index_nml)
    obsgroup=obstypedic["obsgroup"]
    subtype=obstypedic["subtype"][0]
    subtypegroup=[subtype]*batch_count
    obs_index=obstypedic[str(subtype)]
    elistgroup=[obstore.obstore_create_element_table(nmlfile,obs_index)]*batch_count
    print(elistgroup)
    datagroup=read_bufr_dump(INFILE_LIST,HEADER_OFFSET,RECORD_LENGTH,bufr_ele_list)
    print(datagroup)
    print(DT)
    with open(output_file, "wb+") as outfile:
        (datapos,datalen,dataend)=obstore.create_obstore(DT,outfile,nmlfile,obsgroup,subtypegroup,elistgroup,datagroup,batchcount=batch_count,header_offset=339,maxindx=maxindx,lut_ncols=LUTSIZE)
    print("Writting to "+output_file+ " is completed. Data position:"+str(datapos)+" Data length:"+str(datalen)+" Data end:"+str(dataend))
    obsmod.obs_frame(datagroup,subtypegroup,outpath,filename=obstype,option=1)
    return(datagroup)

