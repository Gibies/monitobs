#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 15:29:58 2019

@author: gibies
"""
from __future__ import print_function
import subprocess
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
obs_index_nml=OBSNML+"/obs_index_nml"
odb_index_nml=OBSNML+"/odb_index_nml"
varobs_nml=OBSNML+"/varobs_nml"
varcx_surf_nml=OBSNML+"/varcx_surf_nml"
varcx_uair_nml=OBSNML+"/varcx_uair_nml"
import obslib
import obsdic
import obstore
import fixheader
import ecbufr
#import ncbufr
#import varobs
#import varcx
#import fsoi
import symobs
import daview
import essio
import sqlobs
#import sqlodb
#import obsgui
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pyplot
import pandas
import numpy
import glob
#import odb

diaglev=int(os.environ.get('GEN_MODE',0))
MAXINDX=int(os.environ.get('MAXINDX',fixheader.MAXINDX))
HDRSIZE=int(os.environ.get('HDRSIZE',fixheader.HDRSIZE))
LUTSIZE=int(os.environ.get('LUTSIZE',fixheader.LUTSIZE))
HBpos=int(os.environ.get('HBpos',fixheader.HBpos))
HBlen=int(os.environ.get('HBlen',fixheader.HBlen))
HCpos=int(os.environ.get('HCpos',fixheader.HCpos))
HClen=int(os.environ.get('HClen',fixheader.HClen))
obs_nml=os.environ.get('obs_index_nml',OBSNML+"/obs_index_nml")
odb_nml=os.environ.get('odb_index_nml',OBSNML+"/odb_index_nml")
varno_nml=os.environ.get('odb_varno_nml',OBSNML+"/odb_varno_nml")
subtype_nml=os.environ.get('obs_subtype_nml',OBSNML+"/obs_subtype.nml")

def iri_load_cubes(infile,cnst=None,callback=None,stashcode=None,option=0):
	return(daview.iri_load_cubes(infile,cnst=cnst,callback=callback,stashcode=stashcode,option=option))

def irx_load_cubray(infile,varname,callback=None,stashcode=None,option=2,dims=None,coords=None):
	return(daview.irx_load_cubray(infile,varname,callback=callback,stashcode=stashcode,option=option,dims=dims,coords=coords))

def today(fmtstr="%Y%m%d"):
	return(obslib.today(fmtstr))

def globlist(string):
	return(obslib.globlist(string))

def DataFrame():
	return(obslib.DataFrame())

def reset_index(data):
	return(obslib.reset_index(data))

def ecbufr_decode_files(inpath,Tnode,slctstr,nmlfile,eleindxmaptbl=None,elemlist=None,subtype=None,keyfieldlst=[]):
	return(ecbufr.bufr_decode_files(inpath,Tnode,slctstr,nmlfile,eleindxmaptbl,elemlist=elemlist,subtype=subtype,keyfieldlst=keyfieldlst))

def ncbufr_decode_files(inpath,Tnode,slctstr,nmlfile):
	return(ncbufr.bufr_decode_files(inpath,Tnode,slctstr,nmlfile))

def ncbufr_test_read(slctstr,txtfile=None,field=None):
	return(ncbufr.test_read(slctstr,txtfile,field))

def get_nml(nltype):
    return{
    "obstore"	: obs_index_nml,
    "varobs"	: varobs_nml,
    "cx_surf"	: varcx_surf_nml,
    "cx_uair"   : varcx_uair_nml,
    }.get(nltype,None) 

diaglev=int(os.environ.get('GEN_MODE',0))
def errprint(*args, **kwargs):
    if diaglev > 0: print(*args, file=sys.stderr, **kwargs)

def pydate(datestring,hour=0,minute=0,second=0):
    return(obslib.pydate(datestring,hour,minute,second))

def cylcdate(datetime):
    return(obslib.cylcdate(datetime))

def odbdate(datetime):
    return(obslib.odbdate(datetime))

def binary_file_read(infilename,readpos=1,readwidth=1,reclen=1,readformat="q"):
    with open(infilename, "rb") as infileptr:
        data=obslib.binary_read_data(infileptr,readpos,readwidth,reclen,readformat)
    return(data)

def binary_file_record(infilename,readpos=1,numrec=1,reclen=1,readformat="q"):
    with open(infilename, "rb") as infileptr:
        data=obslib.binary_read_record(infileptr,readpos,numrec,reclen,readformat)
    return(data)

def ascii_file_write(data,outfile=None,option=0):
    outfile=obslib.obs_frame_ascii(data,outfile,option)
    return(outfile)

def obs_frame_ascii(data,outfile=None,option=0):
    return(obslib.obs_frame_ascii(data,outfile,option))

def binary_file_segment_read(infilename,sec_nam="lut"):
    with open(infilename, "rb") as infileptr:
       data = obslib.binary_file_segment_read(infileptr,sec_nam)
    return(data)

def varobs_get_varlist(infilename):
    with open(infilename, "rb") as infileptr:
        data = varobs.varobs_get_varlist(infileptr)
    return(data)

def binary_read_data_header(infilename):
    with open(infilename, "rb") as infileptr:
       (halp,hbet,hgam,ldc,rdc,cdc,lut,necklace)=obslib.binary_read_data_header(infileptr)
    return(halp,hbet,hgam,ldc,rdc,cdc,lut,necklace)

def varobs_batch_info(infilename,batchnum=1):
    with open(infilename, "rb") as infileptr:
       batch_header=varobs.varobs_batch_info(infileptr,batchnum)
    return(batch_header)

def varobs_read_record(infilename,recnum=1):
    with open(infilename, "rb") as infileptr:
       record=varobs.varobs_read_record(infileptr,recnum)
    return(record)

def varobs_rec_meta(infilename,recnum=1):
    with open(infilename, "rb") as infileptr:
       data=varobs.varobs_rec_meta(infileptr,recnum)
    return(data)

def varobs_get_data(infilename,batchnum=1,obsnum=1,varnum=1,colnum=1,levnum=1):
    with open(infilename, "rb") as infileptr:
       data=varobs.varobs_get_data(infileptr,batchnum,obsnum,varnum,colnum,levnum)
    return(data)

def varobs_read_batch(infilename,bthnum=1):
    with open(infilename, "rb") as infileptr:
       data=varobs.read_batch(infileptr,bthnum)
    return(data)

def varcx_read_batch(infilename,bthnum=1):
    with open(infilename, "rb") as infileptr:
       data=varcx.read_batch(infileptr,bthnum)
    return(data)

def varobs_read_batch_tail(infilename,bthnum=1):
    with open(infilename, "rb") as infileptr:
       data=varobs.varobs_read_batch_tail(infileptr,bthnum)
    return(data)

def varobs_read_cdc(infilename,bthnum=1):
    with open(infilename, "rb") as infileptr:
       data=varobs.varobs_read_cdc(infileptr,bthnum)
    return(data)

def varobs_get_latlon(infilename,bthnum=1):
    with open(infilename, "rb") as infileptr:
       data=varobs.get_latlon(infileptr,bthnum)
    return(data)

def varobs_get_time(infilename,bthnum=1):
    with open(infilename, "rb") as infileptr:
       data=varobs.get_time(infileptr,bthnum)
    return(data)

def proactive_stnlist(fsoi_info_path,outpath,freq_cutoff=None):
    return(fsoi.proactive_stnlist(fsoi_info_path,outpath,freq_cutoff))

def frame_data(values,columns,index=None,data=None):
    return(obslib.frame_data(values,columns,index,data))

def get_numerics(string_data):
    return(obslib.get_numerics(string_data))

def obsguiobj(mainbox=None, ROSE_SUITE_DIR=None, OUTFILE=None, obs_index_nml=None, odb_index_nml=None):
    if obs_index_nml is None: obs_index_nml=obs_nml
    if odb_index_nml is None: odb_index_nml=odb_nml
    return(obsgui.obsdata(mainbox, ROSE_SUITE_DIR, OUTFILE, obs_index_nml, odb_index_nml))

def symulate_obstore(outfile,infile,DT,nature_filevar,nmlfile=obs_nml):
    symobs.symulate_obstore(outfile,infile,DT,nature_filevar,nmlfile)

def obstore_read_batch_elements(obsfile_name,nmlfile=obs_nml):
    with open(obsfile_name, "rb") as obsfile:
       elist = obstore.obstore_read_batch_elements(obsfile,nmlfile)
    return(elist)

def obstore_batchinfo(obsfile_name,nmlfile=obs_nml):
    with open(obsfile_name, "rb") as obsfile:
        batchinfo=obstore.obstore_headerinfo(obsfile,nmlfile)
        obstore.print_batchinfo(batchinfo)
        for indx in batchinfo.Batch_Index.values:
            print(indx)
            if indx not in [0]: print(obstore.obstore_read_subhead_segment(obsfile,"lut",indx,1,LUTSIZE))
        #print(batchinfo[["Batch_Index","Batch_Subtype","Batch_Data_Position","Batch_Obs_Count","Batch_Col_Count","Batch_Data_Length","Batch_Data_End"]].loc[[0,1]])
        #print(batchinfo['elist_map'].loc[1])
        #print(obstore.binpos(obsfile,"ldc"))
        #print(obstore.binpos(obsfile,"rdc"))
        #print(obstore.binpos(obsfile,"cdc"))
        #print(obstore.binpos(obsfile,"lut"))
        #print(obstore.binpos(obsfile,"data")[0])
    return(batchinfo)

def obstore_read_index_subtype(obsfile,indx):
    return(obstore.obstore_read_index_subtype(obsfile,indx))
    
def obstore_read_subtype_index(obsfile,subtype):
    return(obstore.obstore_read_subtype_index(obsfile,subtype))

def query_obstore(obsfile,nmlfile=obs_nml,subtype=None,indx=None,selectlist=[],userquery=[]):
    return(sqlobs.query(obsfile,nmlfile,subtype ,indx,selectlist,userquery))

def odb_list_varno(odbfile):
    return(sqlodb.odb_list_varno(odbfile))

def odb_list_varname(odbfile,varno_nmlfile=varno_nml):
    varnolist=odb_list_varno(odbfile)
    return(obslib.listvarname(varno_nmlfile,varnolist))
    
def odb_get_varnolist(varnamlist,varno_nmlfile=varno_nml):
    return(obslib.listvarno(varno_nmlfile,varnamlist))

def odb_filter_varno(odbfile,varno):
    return(sqlodb.odb_filter_varno(odbfile,varno))

def getodbname(opsname,nmlfile=odb_nml):
    return(obslib.getodbname(nmlfile,opsname))

def getopsname(odbname,nmlfile=odb_nml):
    return(getopsname(nmlfile,odbname))

def odb_renamefield(data,nmlfile=odb_nml):
    return(obslib.odb_renamefield(data,nmlfile))
    
def query_odb(odbfile,nmlfile=odb_nml,subtype=None,selectlist=[],varnolist=[],userquery=[]):
    print(varnolist)
    data=sqlodb.query(odbfile,nmlfile,subtype,selectlist,varnolist,userquery)
    return(odb_renamefield(data))

def odb_list_subtype(odbfile):
    return(sqlodb.odb_list_subtype(odbfile))
    
def obstore_list_subtype(obsfile):
    return(obstore.obstore_read_subtype(obsfile))

def getelenams(obsfile,nmlfile=obs_nml,subtype=None,indx=None):
    return(obstore.getelenams(obsfile,nmlfile,subtype,indx))
    
def obstore_print_element_table(obsfile,nmlfile=obs_nml):
    with open(obsfile, "rb") as infile:
        subtype_list=obstore.obstore_read_subtype(infile)
        for indx,subtype in enumerate(subtype_list,start=1):
            print(indx,subtype)
            print(obstore.obstore_read_batch_elements(infile,indx,nmlfile))

def plot_gridmean(plotpath,data,cylcdatestr,prefix,varname,element,fill=False,extend="both",nmlfile=obs_nml):
    return(daview.mpl_plot_gridmean(plotpath,nmlfile,data,cylcdatestr,prefix,varname,element,fill=fill,extend=extend))

def print_obstore_lite(obs_file,nmlfile=obs_nml,option=diaglev):
    headerdiff=obstore.header_diffcheck(obs_file)
    with open(obs_file, "rb") as obsfile:
        #print(obstore.obstore_read_header(obsfile,1,305))
        #print(obstore.obstore_read_real(obsfile,306,34))
        if option > 1 : print(obstore.obstore_read_subtype(obsfile))
        for idx,subtype in enumerate(obstore.obstore_read_subtype(obsfile)[0:],start=1):
            print(idx, subtype)
            print(obstore.obstore_read_batch_elements(obsfile,idx,nmlfile))
            if option > 1 : obstore.print_data_batch(obsfile,nmlfile,int(subtype),idx)
        return(headerdiff)
    
def print_obstore(obs_file,nmlfile=obs_nml,selectlist=[],querystring="",option=diaglev):
    obstore.header_diffcheck(obs_file)
    with open(obs_file, "rb") as obsfile:
        #print(obstore.obstore_read_header(obsfile,1,305))
        #print(obstore.obstore_read_real(obsfile,306,34))
        if diaglev > 5 : print(obstore.obstore_read_subtype(obsfile))
        for idx,subtype in enumerate(obstore.obstore_read_subtype(obsfile)[0:],start=1):
            print(idx, subtype)
            print(obstore.obstore_read_batch_elements(obsfile,idx,nmlfile))
            data=sqlobs.query(obsfile,nmlfile,int(subtype),idx,selectlist,querystring)
            obslib.print_frame(data,option)
        return(data)

def obs_latlon_plot(datapath,plotpath,nmlpath=OBSNML,maxindx=MAXINDX,obstypelist=[],fltrkey="subtype",text="",title="",pltnam=None):
    for obstype in obstypelist:
        print(obstype)
        obstypedic=obsdic.obstype[obstype]
        filename=obstypedic["filename"]
        data_file=datapath+"/"+filename
	if pltnam is None:
        	plotfile=plotpath+"/"+filename+title+".png"
        	textfile=plotpath+"/"+filename+title+".txt"
	else:
		plotfile=plotpath+"/"+pltnam+".png"
        	textfile=plotpath+"/"+pltnam+".txt"
        latlon_data=obstore_read_latlon(data_file,fltrkey=fltrkey,maxindx=maxindx)
        obslib.obs_frame_ascii(latlon_data,textfile,option=diaglev)
        figure1=daview.mpl_plot_latlon(latlon_data,plotfile,fltrkey=fltrkey,text=text,title=title)

def elenam_list(infile,indx=1,nmlfile=obs_nml):
	elist=obstore.obstore_read_batch_elements(infile,indx,nmlfile)
	element_list=elist.Element.values
	return(element_list)

def obstore_read_latlon(obs_file,fltrkey="subtype",nmlfile=obs_nml,option=diaglev,maxindx=None):
    print("inside obstore_read_latlon function inside obsmod")
    #obstore.header_diffcheck(obs_file)
    data_list=[]
    with open(obs_file, "rb") as obsfile:
        maxindx=obstore.obstore_read_header(obsfile,111,1)
        print(maxindx)
        maxindx=obstore.obstore_read_header(obsfile,116,1)
        print(maxindx)
        maxindx=obstore.obstore_read_header(obsfile,121,1)
        print(maxindx)
        if fltrkey in elenam_list(obsfile):
           selectlist=["Latitude","Longitude",fltrkey]
        else:
           selectlist=["Latitude","Longitude"]
        print(selectlist)
        subtypelist=obstore.obstore_read_subtype(obsfile).flatten()
        print(subtypelist)
        for idx,subtype in enumerate(subtypelist[0:],start=1):
            print(idx,subtype)
            #print(obstore.obstore_read_batch_elements(obsfile,idx,nmlfile))
            data=sqlobs.query(obsfile,nmlfile,int(subtype),idx,selectlist=selectlist,maxindx=maxindx)
            data["subtype"]=pandas.Series([subtype for x in range(len(data.index)+1)]) 
            data_list=data_list+[data]
        latlon_data=obslib.obs_merge_batch(data_list)
    return(latlon_data)

def plotallvar_obs(plotpath,infile,cylcdatestr,obstype,nmlfile=obs_nml,subtype_nmlfile=subtype_nml,fill=False):
    with open(infile, "rb") as infile:
        subtype_list=obstore.obstore_read_subtype(infile)
        print(subtype_list)
        for indx,subtype in enumerate(subtype_list,start=1):
            subtype_name=obslib.get_subtype_name(subtype_nmlfile,subtype)
            prefix="obstore_"+obstype+"_"+subtype_name
            elist=obstore.obstore_read_batch_elements(infile,indx,nmlfile)
            element_list=elist.Element.values
            ldc=elist.LDC.values
            for i,element in enumerate(element_list):
                if element not in ["CharData","Year","Month", "Day", "Hour","Minute", "Second","Latitude","Longitude","WMOBlockNo","WMOStnNo", "WMORegNo", "StationReportType","RPRT_IDNY","CallSign","TailNumber"]:
                    data=sqlobs.query(infile,nmlfile,subtype,indx,["Latitude","Longitude",element],nanqlist=[element],nanvalue=-1073741824.0,minlev=10)
                    print(obslib.dfheader(data))
                    ncols=obstore.getldc(element,obsfile=infile,nmlfile=nmlfile,indx=indx)
                    if ncols > 1: varname=element+str(1)
                    else: varname=element
                    daview.mpl_scatterplot(plotpath,nmlfile,data,cylcdatestr,prefix,varname,element,fill=fill)
                    daview.mpl_plot_density(plotpath,nmlfile,data,cylcdatestr,prefix,varname,element,fill=fill)
                    daview.mpl_plot_gridmean(plotpath,nmlfile,data,cylcdatestr,prefix,varname,element,fill=fill)
        
def plotallvar_odb(plotpath,odbfile,cylcdatestr,obstype,odbnmlfile=odb_nml,varno_nmlfile=varno_nml,subtype_nmlfile=subtype_nml,fill=False):
    subtype_list=sqlodb.odb_list_subtype(odbfile)
    print(subtype_list)
    for subtype in subtype_list:
        subtype_name=obslib.get_subtype_name(subtype_nmlfile,subtype)
        prefix="odb_"+obstype+"_"+subtype_name
        varno_list=sqlodb.odb_list_varno(odbfile)
        print(varno_list)
        for varno in varno_list:
            nanvalue=-1073741824.0
            data = query_odb(odbfile,odbnmlfile,subtype,elenams=["lat","lon","varno","obsvalue","obs_error","fg_depar","an_depar"],querystring="varno="+ str(varno)+" and obsvalue != "+str(nanvalue))
            #sqlodb.sqlodb(odbfile,'select lat,lon,varno,obsvalue,obs_error,fg_depar,an_depar from "' + odbfile + '" where varno='+ str(varno) +' ;')
            varname=obslib.getvarname(varno_nmlfile,varno)
            long_name=obslib.getlongname(varno_nmlfile,varno)
            daview.mpl_scatterplot(plotpath,odbnmlfile,data,cylcdatestr,prefix,varname,long_name,fill=fill,extend="max")
            daview.mpl_plot_density(plotpath,odbnmlfile,data,cylcdatestr,prefix,varname,long_name,fill=fill,extend="max")
            daview.mpl_plot_gridmean(plotpath,odbnmlfile,data,cylcdatestr,prefix,varname,long_name,fill=fill,extend="both")
            daview.mpl_plot_depart_firstguess(plotpath,odbnmlfile,data,cylcdatestr,prefix,varname,long_name,fill=fill,extend="both")
            daview.mpl_plot_depart_anal(plotpath,odbnmlfile,data,cylcdatestr,prefix,varname,long_name,fill=fill,extend="both")

def obs_frame(datagroup=None,subtypegroup=None,outpath=None,filename="output",option=0,tagmark="",text="",maxindx=MAXINDX,obstore_info=None):
	if obstore_info is not None:
		if "datagroup" in obstore_info: datagroup=obstore_info["datagroup"]
		if "subtypegroup" in obstore_info: subtypegroup=obstore_info["subtypegroup"]
		if "outpath" in obstore_info: outpath=obstore_info["outpath"]
		if "filename" in obstore_info: filename=obstore_info["filename"].replace(".","_")
	plotfile=outpath+"/"+filename+".png"
	print(plotfile)
	if datagroup is not None: figure1=daview.mpl_plot_location(datagroup,plotfile,tagmark=tagmark,lblst=subtypegroup,text=text)
	#figure1.show()
	for idx in range(0,len(datagroup),1):
           data=datagroup[idx]
           textfile=outpath+"/"+filename+"_batch_"+str(idx+1)+".txt"
           print(textfile)
           if data is not None: obslib.obs_frame_ascii(data,textfile,option)
    
def symobs_main(Tnode,outpath,inpath,nmlpath,obstypelist=[],maxindx=MAXINDX,subtypelist=None):
    if len(obstypelist) == 0:
        obstypelist=obsdic.obstypelist
    for obstype in obstypelist:
        datagroup=symobs.symulate_obstore(outpath,inpath,nmlpath,Tnode,obstype,maxindx=maxindx,subtypelist=subtypelist)
        ######################
        
def obstore_write(data,keynmlfile,outpath,btchcnt=None,cntmax=None,DT=None,diagflag=0,missing_value=-1073741824.00000):
	outfile=obstore.obstore_write(data,keynmlfile,outpath,btchcnt=btchcnt,cntmax=cntmax,DT=DT,diagflag=diagflag,missing_value=missing_value)
	return(outfile)

    
def obstore_create_file(obstore_info,diagflg=0,callsignflag=False,filedata=None):
	obstore_info=obstore.obstore_create_file(obstore_info,diagflg=diagflg,callsignflag=callsignflag,filedata=filedata)
	obs_frame(obstore_info=obstore_info,option=1)
	return(obstore_info)


def create_obstore(obstore_info,diagflg=0,callsignflag=False,filedata=None):
	obstore_info=obstore_create_file(obstore_info,diagflg=diagflg,callsignflag=callsignflag,filedata=filedata)
	return(obstore_info)

def datset_extract(infile,varlst=None,dimlst=None,outpath=None,outfile=None,callback=None,stashcode=None,option=2,diagflg=0):
	datset=essio.datset_extract(infile,varlst=varlst,dimlst=dimlst,outpath=outpath,outfile=outfile,callback=callback,stashcode=stashcode,option=option,diagflg=diagflg)
	return(datset)

def datset_save(datset,outpath=None,outfile=None,infile=None):
	outfile=essio.datset_save(datset=datset,outpath=outpath,outfile=outfile,infile=infile)
	return(outfile)

def datset_print(datset):
	essio.datset_print(datset)

def mkdir(path):
	obslib.mkdir(path)

