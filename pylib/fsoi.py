#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 03 15:45:39 2021

@author: gibies
"""
from __future__ import print_function
import sys,os
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
import numpy
import pandas
import shutil
import glob
   
stnlst_path=os.environ.get('StnLstDir',CYLCPATH+"/share/data/etc/stationlists/atmos")
varnml=os.environ.get('VarNml',OBSNML+"/varobs_nml")

diaglev=int(os.environ.get('GEN_MODE',0))
def errprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def read_fsoi(fsoi_info_string,outpath,workdir,yyyymmdd,hh,freq_cutoff=None,stnlst_path=stnlst_path,varnml=varnml,fltr=None):
    data=obslib.DataFrame()
    outfile=None
    clmnhdr=None
    colspecs=None
    infile_list=glob.glob(fsoi_info_string)
    if len(infile_list) == 0 : errprint("Error: No input file found at " +fsoi_info_string)
    else : print(infile_list)
    if len(infile_list) > 1 : errprint("Error: Multiple input files found at " +fsoi_info_string)
    checklist=workdir+"/checklist"
    if os.path.exists(checklist):
       os.remove(checklist)
    obslib.mkdir(workdir)
    with open(checklist,"w") as cklst:
        cklst.write("obstype\tfilelist\n")
    for infile in infile_list:
        #if clmnhdr is None : clmnhdr=["num", "obval", "inova", "sens", "lat", "lon", "obpres", "group_type", "subtype", "idt_ob", "obserr", "bkgerr", "obstype_callsign_chnlno"]
        if clmnhdr is None : clmnhdr=["num", "obval", "inova", "sens", "lat", "lon", "obpres", "group_type", "subtype", "idt_ob", "obserr", "bkgerr", "obsidtyinfo"]
	if colspecs is None : colspecs=[(1,8),(9,24),(25,40),(41,57),(58,64),(65,71),(72,82),(83,86),(87,92),(93,99),(100,114),(115,121),(122,150)]
        if checklist is None : checklist=workpath+"/checklist"
        infldr=infile.rsplit("/",1)[0]
        fname=infile.rsplit("/",1)[1]
        num_vals=obslib.get_numerics(fname)
        ftime=num_vals[0]
        ftstr=str(ftime)
        print(ftstr)
        if ftstr == str(yyyymmdd)+str(hh) :
           fsoi_data=pandas.read_fwf(infile,  colspecs=colspecs, skiprows=None, header=None, names=clmnhdr, converters={"obstype_callsign_chnlno":str})
	else:
	   print("No data found due to date missmatch")
        data=data.append(fsoi_data)
    if fltr is not None : data=data[data["obsidtyinfo"].str.contains(fltr)]
    return(data)

def fsoi_obstype_list(data):
    obstype_list=[]
    for str in data.obstype_callsign_chnlno.unique():
        str=str.split()[0]
        if str not in obstype_list:
           obstype_list=obstype_list+[str]
    return(obstype_list)


def fsoi_station_list(data,subtype):
	data=data[data.subtype==subtype]
	station_list=[]
	for str in data.obstype_callsign_chnlno.unique():
		strlst=str.split()
		if len(strlst) > 1:
			str=strlst[1]
		else:
			str=strlst[0]
		if str not in station_list:
			station_list=station_list+[str]
	return(station_list)

def fsoi_subtype_nml(data):
    obstype_list=fsoi_obstype_list(data)
    nmldata=None
    for obstype in obstype_list:
        subtype_list=data[data["obstype_callsign_chnlno"].str.contains(obstype)].subtype.unique()
        for subtype in subtype_list:
	    nmldata=obslib.frame_data([[subtype,obstype]],columns=["subtype","obstype"],data=nmldata)
    return(nmldata)

def fsoi_station_nml(data):
    subtype_nml=fsoi_subtype_nml(data)
    for subtype in subtype_nml.subtype:
	obstype=subtype_nml[subtype_nml.subtype==subtype].obstype
	station_list=fsoi_station_list(data,subtype)
	for station in station_list:
            nmldata=obslib.frame_data([[subtype,obstype,station]],columns=["subtype","obstype","station"],data=nmldata) 
	    print(subtype,obstype,station)
    return(nmldata)

def fsoi_subtype_list(fsoi_info_string,outpath,workdir,yyyymmdd,hh,freq_cutoff=None,stnlst_path=stnlst_path,varnml=varnml):
    data=read_fsoi(fsoi_info_string,outpath,workdir,yyyymmdd,hh,freq_cutoff=None,stnlst_path=stnlst_path,varnml=varnml)
    subtype_list=data.subtype.unique()
    return(subtype_list)

def fsoi_subtype_read(subtype,fsoi_info_string,outpath,workdir,yyyymmdd,hh,freq_cutoff=None,stnlst_path=stnlst_path,varnml=varnml):
    data=read_fsoi(fsoi_info_string,outpath,workdir,yyyymmdd,hh,freq_cutoff=None,stnlst_path=stnlst_path,varnml=varnml)
    subtype_data=obslib.frame_select_filter(data,"subtype",int(subtype))
    return(subtype_data)

def proactive_stnlist(fsoi_info_string,outpath,workdir,yyyymmdd,hh,freq_cutoff=None,stnlst_path=stnlst_path,varnml=varnml):
    outfile=None
    infile_list=glob.glob(fsoi_info_string)
    if len(infile_list) == 0 : errprint("Error: No input file found at " +fsoi_info_string)
    else : print(infile_list)
    checklist=workdir+"/checklist"
    if os.path.exists(checklist):
       os.remove(checklist)
    obslib.mkdir(workdir)
    with open(checklist,"w") as cklst:
        cklst.write("obstype\tfilelist\n")
    clmnhdr=["lat", "lon", "plev", "varno", "subtype", "time", "obstype", "callsign"]
    for infile in infile_list:
        (obstype,wrkfile_list)=get_stn_list(infile,workdir,yyyymmdd,hh,freq_cutoff,stnlst_path,varnml,checklist,clmnhdr=clmnhdr) 
        with open(checklist,"a") as cklst:
             cklst.write(obstype + "\t" + wrkfile_list + "\n")
    ckdb=pandas.read_table(checklist, skiprows=None, header=0)
    print(ckdb)
    for obstype in ckdb.obstype.values:
        wrkfile_list=workdir+"/"+str(obstype)+"/file_list"
        print(wrkfile_list)
        outfile=generate_stnlst_nl(wrkfile_list,outpath,workdir,obstype,stnlst_path)
    return(outfile)


def get_stn_list(infile,workpath,cylcdate,cylchour,freq_cutoff=None,stnlst_path=stnlst_path,varnml=varnml,checklist=None,clmnhdr=None):
    outfile=None
    if clmnhdr is None : clmnhdr=["num", "obval", "xinov", "sens1", "lat", "lon", "obpres", "group_type", "instype", "idt_ob", "obserr", "bkgerr", "obstype", "callsign", "chnlno"]
    if checklist is None : checklist=workpath+"/checklist"
    infldr=infile.rsplit("/",1)[0]
    fname=infile.rsplit("/",1)[1]
    num_vals=obslib.get_numerics(fname)
    ftime=num_vals[0]
    ftstr=str(ftime)
    if ftstr == cylcdate+cylchour :
      file_data=pandas.read_fwf(infile, skiprows=None, header=None, names=clmnhdr, converters={"callsign":str})
      obstype=get_obstype(file_data.obstype.values[0])
      workdir=workpath+"/"+str(obstype)
      wrkfile_list=workdir+"/file_list"
      ckdb=pandas.read_table(checklist, skiprows=None, header=0)
      if obstype not in ckdb.obstype.values and os.path.exists(wrkfile_list) : 
         print(wrkfile_list + " is to be removed")
         os.remove(wrkfile_list)
      wrkfile_list=subtype_filter(file_data,workdir,cylchour,obstype,varnml,freq_cutoff,wrkfile_list)
    else: errprint("Error: Timestamp of the input file name missmatch : "+ ftstr + " != " + yyyymmdd+hh)
    return(obstype,wrkfile_list)

def get_obstype(OBSTYPE):
    return{
    "SYNOP" : "Surface",
    "BUOY"  : "Surface",
    "SONDE" : "Sonde",
    "TEMP"  : "Sonde",
    }.get(OBSTYPE,"")

def subtype_filter(data,workpath,cylchour,obstype,varnml=varnml,freq_cutoff=None,file_list=None):
    subtype_list=data.subtype.unique()
    ctlbk=None
    ctldic=None
    for subtype in subtype_list:
        subtype_data=obslib.frame_select_filter(data,"subtype",int(subtype))

def subtype_filter(data,workpath,cylchour,obstype,varnml=varnml,freq_cutoff=None,file_list=None):
    subtype_list=data.subtype.unique()
    ctlbk=None
    ctldic=None
    for subtype in subtype_list:
        subtype_data=obslib.frame_select_filter(data,"subtype",int(subtype))
        workdir=workpath+"/"+str(subtype)
        ctldic=get_ctlblk_header("header_string","&Station ObsType= '"+str(subtype)+"',\n",None)
        ctldic=get_ctlblk_header("comment","\n ! FSOI based rejection for "+str(subtype),ctldic)
        file_list=time_filter(subtype_data,workdir,ctldic.copy(),cylchour,varnml,freq_cutoff,file_list)
    return(file_list)

def time_filter(data,workpath,ctlbk,cylchour,varnml=varnml,freq_cutoff=None,file_list=None):
    time_list=data.time.unique()
    for obstime in time_list:
        ctldic=ctlbk.copy()
        time_data=obslib.frame_select_filter(data,"time",int(obstime))
        workdir=workpath+"/"+str(obslib.obs_clock_hour(obstime,cylchour))
        ctldic=get_ctlblk_header("time_ctlstr",get_time_ctlstr(obstime,cylchour),ctldic)
        ctldic=get_ctlblk_header("comment"," Time "+str(obstime)+" ",ctldic)
        file_list=lev_filter(data,workdir,ctldic.copy(),varnml,freq_cutoff,file_list)
    return(file_list)

def lev_filter(data,workpath,ctlbk,varnml=varnml,freq_cutoff=None,file_list=None):
    lev_list=data.plev.unique()
    mfactor=100
    lev_list=obslib.unique_int(lev_list,mfactor,missing=-9999.9999)
    for obslev in lev_list:
        ctldic=ctlbk.copy()
        lev_btm=int(obslev+(mfactor/2))
        lev_top=int(obslev-(mfactor/2))
        lev_data=obslib.frame_window_filter(data,"plev",lev_top,lev_btm)
        if len(lev_data.index) > 0 :
           workdir=workpath+"/"+str(obslev)
           ctldic=get_ctlblk_header("lev_ctlstr",get_lev_ctlstr(lev_btm,lev_top),ctldic)
           ctldic=get_ctlblk_header("comment"," Lev "+str(obslev)+" ",ctldic)
           file_list=var_filter(lev_data,workdir,ctldic.copy(),varnml,freq_cutoff,file_list)
        else : print(" No data betveen "+str(lev_btm) + " and "+ str(lev_top))
    return(file_list)

def var_filter(data,workpath,ctlbk,varnml,freq_cutoff,file_list):
    var_list=data.varno.unique()
    for varno in var_list:
        ctldic=ctlbk.copy()
        varname=pandas.read_table(varnml, skiprows=None, header=0).query("eleindex == @varno").elename.values[0]
        workdir=workpath+"/"+str(varname)
        ctldic=get_ctlblk_header("var_ctlstr",get_ctlstr(varno),ctldic)
        ctldic=get_ctlblk_header("comment"," of "+str(varname)+" ",ctldic)
        var_data=obslib.frame_select_filter(data,"varno",int(varno))
        file_list=consistancy_check(var_data,workdir,ctldic.copy(),freq_cutoff,file_list)
    return(file_list)

def consistancy_check(data,workdir,ctlbk,freq_cutoff=None,file_list=[]):
    detri_stn_list=data.callsign.unique()
    stats=None
    for indx,stnid in enumerate(detri_stn_list):
        ctldic=ctlbk.copy()
        freq = numpy.count_nonzero(data.callsign == stnid )
        stats = obslib.frame_data([[stnid,freq]],["StnID","Detri_freq"],[indx],stats)
    if freq_cutoff is None: freq_cutoff=min(stats.Detri_freq)
    stn_list=stats.query("Detri_freq >= @freq_cutoff").StnID.values
    outfile=ops_station_namelist_segment(workdir,ctldic.copy(),stn_list)
    with open(file_list,"a") as lst:
         lst.write(outfile+"\n")
    return(file_list)

def ops_station_namelist_segment(workdir,ctlbk,stn_lst):
    outfile=workdir+"/station_id_detrimental.txt"
    filedir=outfile.rsplit("/",1)[0]
    obslib.mkdir(filedir)
    with open(outfile,"w") as op:
        stn_lst_str_arr=stnlist_array_to_string(stn_lst)
        for indx,stn_lst_str in enumerate(stn_lst_str_arr): 
            ctldic=ctlbk.copy()
            if len(stn_lst_str_arr) > 1 : ctldic=get_ctlblk_header("comment"," Batch "+str(indx+1)+"\n\n",ctldic)
            else : ctldic=get_ctlblk_header("comment","\n\n",ctldic)
            if "comment" in ctldic.keys() : op.write(ctldic["comment"])
            if "header_string" in ctldic.keys() : op.write(ctldic["header_string"])
            if "time_ctlstr" in ctldic.keys() : op.write(ctldic["time_ctlstr"])
            if "lev_ctlstr" in ctldic.keys() : op.write(ctldic["lev_ctlstr"])
            if "var_ctlstr" in ctldic.keys() : op.write(ctldic["var_ctlstr"])
            if len(stn_lst_str) == 7: op.write("Id= "+stn_lst_str+" /\n")
            else: op.write("Ids= "+stn_lst_str+" /\n")
            op.write("\n")
    return(outfile)

def stnlist_array_to_string(stn_lst):
    stn_lst_str_arr=[]
    if len(stn_lst) < 1000:
       stn_lst_str_arr=stn_lst_str_arr+[numpy.array2string(stn_lst, separator=',', threshold=100000)[1:-1]]
    else :
       stn_lst_str_arr=stn_lst_str_arr+[numpy.array2string(stn_lst[0:999], separator=',', threshold=100000)[1:-1]]
       stn_lst=stn_lst[1000:]
       stn_lst_str_arr=stn_lst_str_arr+stnlist_array_to_string(stn_lst)
    return(stn_lst_str_arr)

def get_ctlblk_header(handle,string,ctlblk=None):
    if ctlblk is None: ctlblk={}
    if handle in ctlblk.keys() : ctlblk[handle]=ctlblk[handle]+string
    else : ctlblk[handle]=string
    return(ctlblk)

def get_lev_ctlstr(lev_btm,lev_top):
    lev_ctlstr="PressureBottom="+str(lev_btm)+", PressureTop="+str(lev_top)+", \n"
    return(lev_ctlstr)

def get_ctlstr(varno):
    return{
    1	: "RejP='T',\n",
    2	: "RejT='T',\n",
    3	: "RejRH='T',\n",
    4	: "RejUV='T',\n",
    5	: "RejUV='T',\n",
    }.get(varno,"")

def get_time_ctlstr(obstime,cylchour):
    hour_of_day=obslib.obs_clock_hour(obstime,cylchour)
    hour_start=obslib.clock_24_hour(hour_of_day-1)
    hour_end=obslib.clock_24_hour(hour_of_day+1)
    time_ctlstr="HourStart="+str(hour_start)+", HourEnd="+str(hour_end)+", \n"
    return(time_ctlstr)


def generate_stnlst_nl(file_list,outpath,workdir,obstype,stnlst_path):
    obslib.mkdir(outpath)
    outfile=outpath+"/"+str(obstype)+"_stlist.nl"
    template_file=get_stlst_template(obstype,stnlst_path)
    if template_file is None : errprint("Station list template not available")
    infile_list=[]
    infile_list=infile_list+[template_file]
    with open(outfile,"w") as stnlst:
       with open(file_list,"r") as flp: 
          flist=flp.readlines() 
          for filesegment in flist:
             infile_list=infile_list+[filesegment[:-1]]
             #print(infile_list)
       for indx in range(0,len(infile_list)):
            subfile=infile_list[indx]
            print(indx)
            print(subfile)
            with open(subfile,"r") as fp:
                for line in fp:
                   stnlst.write(line)
                stnlst.write("\n")
    return(outfile)

