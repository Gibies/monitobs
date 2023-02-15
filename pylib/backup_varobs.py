#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 15:45:39 2020

@author: gibies

Revised on Thu Aug 13 13:08:00 2020
by Ruchika
for parallel processing
"""

import os,sys
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
CYLCROOT=os.path.dirname(os.path.dirname(os.path.dirname(CURR_PATH)))
CYLCPATH=os.environ.get('CYLCPATH',CYLCROOT)
MONITOBS=os.environ.get('MONITOBS',CYLCPATH+"/modules/monitobs")
LIB=os.environ.get('LIB',MONITOBS+"/pylib")
sys.path.append(LIB)
DIC=os.environ.get('DIC',MONITOBS+"/pydic")
sys.path.append(DIC)
NML=os.environ.get('NML',MONITOBS+"/nml")
sys.path.append(NML)
import textview
import pandas
import numpy
import domaindic
import obslib
#import daview
import multiprocessing as multip
import pymp 

def reshape_frame(data,CYLC_TASK_CYCLE_TIME,period,field):
    data=data[["cylc_dt","height",field]]
    dt=obslib.cylcdate_to_pydate(CYLC_TASK_CYCLE_TIME)
    tdel=obslib.timeperiod(hh=6)
    cycle_list=obslib.dtlist_backward(dt,tdel,period)[:]
    indx=0
    hoff_data=pandas.DataFrame(columns=list(["cylc_dt"])+list(data.height.unique()))
    for i,dt in enumerate(cycle_list):
	cylc_dt=obslib.cylcdate(dt)
	cycle_data=obslib.frame_select_filter(data,"cylc_dt",cylc_dt)[field].values
	print(cycle_data)
	hoff_data=obslib.append(hoff_data,[cylc_dt]+list(cycle_data))
    return(hoff_data)

def extract_time_height(data,dt,height):
        cylc_dt=obslib.cylcdate(dt)
        data=obslib.frame_select_filter(data,"date",cylc_dt)
        data=obslib.frame_window_filter(data,"height_cog",(height-500),(height+500))
        obscount=len(data)
        stat_table=pandas.DataFrame(columns=["cylc_dt","height","obs_count","obsmean","bkgmean","ombmean"])
        if obscount > 0:
           obsmean=obslib.mean_of(data.observed.values,absol=True)
           bkgmean=obslib.mean_of(data.background.values,absol=True)
           ombmean=obslib.mean_of(data.obs_minus_bkg.values)
        else:
           obsmean=numpy.nan
           bkgmean=numpy.nan
           ombmean=numpy.nan
        obslib.append(stat_table,[cylc_dt,height,obscount,obsmean,bkgmean,ombmean])
        return(stat_table)

def get_hofmuller(full_data,CYLC_TASK_CYCLE_TIME,period):
    dt=obslib.cylcdate_to_pydate(CYLC_TASK_CYCLE_TIME)
    tdel=obslib.timeperiod(hh=6)
    cycle_list=obslib.dtlist_backward(dt,tdel,period)[:]
    indx=0
    for i,dt in enumerate(cycle_list):
      for height in range(1000,20000,1000):
        if indx == 0:
           stats=extract_time_height(full_data,dt,height)
        else:
           stats=stats.append(extract_time_height(full_data,dt,height))
        indx+=1
    return(stats)

def get_stat(data,domain,varname="hloswind",absol=False,report=False,plot=True,plotfile=None):
	lonmin=domaindic.dom_lonmin[domain]
	lonmax=domaindic.dom_lonmax[domain]
	latmin=domaindic.dom_latmin[domain]
	latmax=domaindic.dom_latmax[domain]
	hgtmin=domaindic.dom_hgtmin[domain]
	hgtmax=domaindic.dom_hgtmax[domain]
	domstat=pandas.DataFrame(columns=["domain","lonmin","lonmax","latmin","latmax","hgtmin","hgtmax","obscount","obsmean","bkgmean","ombmean"])
        obscount=len(data)
        if obscount > 0:
                obsmean=obslib.mean_of(data.observed.values,absol=absol)
                bkgmean=obslib.mean_of(data.background.values,absol=absol)
                ombmean=obslib.mean_of(data.obs_minus_bkg.values)
        else:
                obsmean=numpy.nan
                bkgmean=numpy.nan
                ombmean=numpy.nan
        obslib.append(domstat,[domain,lonmin,lonmax,latmin,latmax,hgtmin,hgtmax,obscount,obsmean,bkgmean,ombmean])
        return(domstat)


def domain_stat(data,domain,varname="hloswind",absol=False,report=False,plot=True,plotfile=None):
	data=domain_filter(data,domain,varname="hloswind",absol=False,report=False,plot=True,plotfile=None)
	domstat=get_stat(data,domain,varname="hloswind",absol=False,report=False,plot=True,plotfile=None)
        if report:
                print("________________________________")
                print(domaindic.dom_name[domain])
                print("Latitude : "+str(latmin)+" to "+str(latmax))
                print("Height : "+str(hgtmin)+" to "+str(hgtmax))
                print("Observation count : "+str(obscount))
                print("Mean observation : "+str(obsmean))
                print("Mean background : "+str(bkgmean))
                print("Mean observation minus background : "+str(ombmean))
                print("________________________________")
        if plot:
               #obsplot.plot_field(data,plotfile,domain=domain,varname=varname)
		daview.plot_field(data,plotfile,domain=domain,varname=varname)
        return(domstat)

def domain_filter(data,domain,varname="hloswind",absol=False,report=False,plot=True,plotfile=None):
	lonmin=domaindic.dom_lonmin[domain]
	lonmax=domaindic.dom_lonmax[domain]
	latmin=domaindic.dom_latmin[domain]
	latmax=domaindic.dom_latmax[domain]
	hgtmin=domaindic.dom_hgtmin[domain]
	hgtmax=domaindic.dom_hgtmax[domain]
	data=obslib.frame_window_filter(data,"height_cog",hgtmin,hgtmax)
	data=obslib.frame_window_filter(data,"Latitude",latmin,latmax)
	data=obslib.frame_window_filter(data,"Longitude",lonmin,lonmax)
	return(data)

def read_index(infile):
	strng="Varfield Levels"
	block=textview.get_visual_block(infile,strng)
	varindx= textview.column_numeric(block[1:-1],-2)
	return(varindx)

def head_line_num(infile):
	header_line_count = textview.get_line_number(infile,"batch ob num  level  field")
	return(header_line_count)

def read_data(infile,offsetpos,data_len=1):
	header_line_count = head_line_num(infile)
	data=textview.get_lines(infile,(header_line_count+offsetpos),num_lines=data_len)
	return(data)

def read_table_header(infile):
	line = head_line_num(infile)
	header=textview.get_lines(infile,line)
	col_head=[]
	tmp_word=""
	i=-1
	for word in header[0].split():
		word=textview.alnum(word)
		if word == "ob": 
			tmp_word=word
		else:
			if tmp_word == "":
				col_head.append(word)
			else:
				col_head.append(tmp_word+"_"+word)
			tmp_word=""
	return(col_head)

def count_obs(infile):
	varindx=read_index(infile)
	fld_cnt=len(varindx)
	file_line_count=textview.get_file_length(infile)
	header_line_count=head_line_num(infile)
	data_len=file_line_count-header_line_count
	obscount=data_len/fld_cnt
	return(fld_cnt,obscount)

def list_modify(inlist,reject_list=[],action="trim"):
	newlist=[]
	for item in inlist:
		if item not in reject_list: newlist.append(item)
	return(newlist)

def frame_data(infile,filetype="varobs",obscnt=None,offset=0):
	(field_count,obs_count)=count_obs(infile)
	if obscnt is not None: obs_count=obscnt
	offsetpos=field_count*offset
	data_len=field_count*obs_count
	header=read_table_header(infile)
	if filetype == "modelobs" : header=list_modify(header,["Callsign"])
	data=read_data(infile,(offsetpos+1),data_len)
	index=numpy.arange(1, data_len+1)
	df=pandas.DataFrame(columns=header,index=index)
	for i in index:	
		data1=numpy.array(data[i-1].split())
		if data1.shape == df.loc[i].shape: df.loc[i]=data1
		else:print(filetype+" shape conflict at index :"+str(i)+ " as "+ str(data1.shape) + " and " + str(df.loc[i].shape))
	if "lat" in df: df=df.rename(columns={"lat": "Latitude"})
	if "lon" in df: df=df.rename(columns={"lon": "Longitude"})
	return(df)

def field_query(infile,field_id=103,field_name="height_cog",col_list=["batch","ob_num"],filetype="varobs",obscnt=None,offset=0):
	obs_df=frame_data(infile,filetype=filetype,obscnt=obscnt,offset=offset)
	field_df=obs_df[obs_df.field.values == str(field_id)]
	field_df=field_df.rename(columns={"ob_value": field_name})
	col_list=col_list+[field_name]
	field_df=field_df[col_list]
	return(field_df)

def varobs_read(obsin,obscnt=None,offset=0,obs_field="obs_data",field_id=98,col_list=["batch","ob_num","field","level","Latitude","Longitude"]):
	(field_count,obs_count)=count_obs(obsin)
	if obscnt is None: obscnt=obs_count
	field_data=field_query(obsin,field_id=field_id,field_name=obs_field,col_list=col_list,filetype="varobs",obscnt=obscnt,offset=offset)
	return(field_data)

def omb_calc(thread,datalist,obsin,bkgin,obscnt=None,offset=0,obs_field="hlosw_obs",bkg_field="hlosw_bkg",field_id=98):
     (field_count,obs_count)=count_obs(obsin)
     if obscnt is None: obscnt=obs_count
     index=thread.range(obscnt)
     num=thread.thread_num
     data=pandas.DataFrame()
     for indx in index:
	obscnt=1
	offset=indx
     	obs_list=["batch","ob_num","field","level","Latitude","Longitude","pge","ob_error"]
     	bkg_list=["batch","ob_num","field"]
     	obs_df=field_query(obsin,field_id=field_id,field_name=obs_field,col_list=obs_list,filetype="varobs",obscnt=obscnt,offset=offset)
     	bkg_df=field_query(bkgin,field_id=field_id,field_name=bkg_field,col_list=bkg_list,filetype="modelobs",obscnt=obscnt,offset=offset)
     	newdf=obs_df.merge(bkg_df,on=bkg_list)
     	newdf['obs_minus_bkg']=newdf[obs_field].astype(numpy.float)-newdf[bkg_field].astype(numpy.float)
     	height_cog=field_query(obsin,field_id=103,field_name="height_cog",filetype="varobs",obscnt=obscnt,offset=offset)
     	newdf=newdf.merge(height_cog[["batch","ob_num","height_cog"]],on=["batch","ob_num"])
     	#newdf=newdf.rename(columns={"lat": "Latitude", "lon": "Longitude"})
     	data=data.append(newdf,ignore_index=True)
	#datalist=datalist.append((thread.thread_num,data))
     return({num:data})

def sorted_rebuild(datalist):
     data=pandas.DataFrame()
     datadic={}
     for indx in range(len(datalist)):
	data1=datalist[indx]
     	datadic.update(data1)
     for key in sorted(datadic.keys()):
	data1=datadic[key]
	data=data.append(data1,ignore_index=True)
     return(data)

def obs_minus_bkg(obsin,bkgin,obscnt=None,offset=0,omp_num_threads=128):
     obs_field="hlosw_obs"
     bkg_field="hlosw_bkg"
     (field_count,obs_count)=count_obs(obsin)
     if obscnt is None: obscnt=obs_count
     data_len=field_count*obs_count
     datalist=pymp.shared.list()
     with pymp.Parallel(omp_num_threads) as thread:
	data1=omb_calc(thread,datalist,obsin,bkgin,obscnt,offset,obs_field,bkg_field)
	datalist.append(data1)
     data=sorted_rebuild(datalist)
     return(data)

