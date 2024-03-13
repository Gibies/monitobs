#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 17:22:59 2018

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
obs_index_nml="obs_index_nml"
nmlfile="%s/%s" % (OBSNML,obs_index_nml)
import obslib
import obsdic
import obsheader
import obstore
import obsmod
import nature
import ncepradic
import netCDF4
import pandas
import numpy
import math
import itertools
 

MAXINDX=int(os.environ.get('MAXINDX',608))
diaglev=int(os.environ.get('GEN_MODE',0))

def nadir_latlon(Tnow,Tnode,Torbit,OrbInc,NodeShift=0):
    Tearth=1440
    nodlon=obslib.dawnlon(Tnode)
    OrP=obslib.to_minutes(Torbit)
    Tdelta=obslib.to_minutes(Tnow-Tnode)
    lat=(math.sin(Tdelta/OrP*2*math.pi)*math.degrees(math.asin(math.sin(math.radians(OrbInc)))))
    lon=nodlon-NodeShift-(Tdelta/Tearth)*360-(Tdelta/OrP)*360 
                    ###-((Tdelta/Tearth)*360*math.sin(math.radians(OrbInc)))  #+math.degrees(math.cos(2*math.pi*Tdelta/OrP)*math.cos(math.radians(OrbInc)))  #
    lon=degreelon(lon)
    data1={'Time':[Tnow],'Longitude':[lon],'Latitude':[lat]}
    return(pandas.DataFrame(data1,index=[Tdelta],columns=['Time','Latitude','Longitude']))

def los_latlon(SatID, idx, ObsAlt,Tnow,Tnode,Torbit,SMAxis,SatAlt,nadir,OrbInc,HVA,VVA):
    OrP=obslib.to_minutes(Torbit)
    Tdelta=obslib.to_minutes(Tnow-Tnode)
    nadirlon=float(nadir.Longitude.values[0])
    nadirlat=float(nadir.Latitude.values[0])
    RoE=6371  ###in km
    VLOS=(SatAlt-ObsAlt) ###in km
    HLOS=VLOS*math.tan(math.pi*VVA/180)
    OpticDist=math.hypot(HLOS,VLOS)
    losdelon=math.degrees(HLOS/RoE)*math.cos(math.radians(OrbInc+HVA))*math.cos(2*math.pi*Tdelta/OrP)
    losdelat=math.degrees(HLOS/RoE)*math.sin(math.radians(OrbInc+HVA))+math.degrees(HLOS/RoE)*math.cos(math.radians(OrbInc+HVA))*math.sin(2*math.pi*Tdelta/OrP)
    lon=nadirlon+losdelon
    lat=nadirlat+losdelat    
    azmh=get_ddd(losdelon,losdelat)
    ###Caution: Any change in azmh will affect HLOSWind calculation
    data1={'Index':[idx], 'Year':[Tnow.year], 'Month':[Tnow.month], 'Day':[Tnow.day], 'Hour':[Tnow.hour], 'Minute':[Tnow.minute], 'Second':[Tnow.second], 'Longitude':[lon],'Latitude':[lat], 'Altitude':[ObsAlt], 'Azimuth':[azmh]}
    data=pandas.DataFrame(data1,index=[idx],columns=['Index','Year','Month','Day','Hour','Minute','Second','Latitude','Longitude','Altitude','Azimuth'])
    return(data)

def fileindex(data,index,filedim,gph,altnam="Altitude",DT=None):
    data1=data[data.index == index]
    #print(index,data1)
    LTTD=data1.Latitude.values[0]
    LNGD=data1.Longitude.values[0]
    if altnam in data1.columns.values:
        ALTTD=data1.loc[:,altnam].values[0]*1000
    else:
        ALTTD=0
    if "Azimuth" in data1.columns.values: 
        AZMTH=data1.Azimuth.values[0]
    else:
        AZMTH=0
    #Year=data1.Year.values[0]
    time=filedim["time"]
    time_units=filedim["time_units"]
    lat=filedim["lat"]
    lon=filedim["lon"]
    DT=obslib.getdatetime(DT,data1)
    itime=find_nearest_time(time,DT,time_units)
    ilat=find_nearest_index(lat,LTTD)
    ilon=find_nearest_index(lon,LNGD)
    ilev=find_nearest_index(gph[itime,:,ilat,ilon],ALTTD)
    lev=filedim["lev"][ilev]
    findx1={'Index':[index], 'itime':[itime], 'ilat':[ilat], 'ilon':[ilon], 'ilev':[ilev], 'lev':[lev], 'azmh':[AZMTH]}
    findex=pandas.DataFrame(findx1,index=[index],columns=['Index','itime','ilat','ilon','ilev','lev','azmh'])
    return(findex)
    

def degreelon(lon):
    while lon < -180 :
        lon = lon + 360
    while lon > 180 :
        lon = lon - 360
    return(lon)

def ChanFreq(lev,chnldic):
    chlist=chnldic["Chanlist"]
    chtopdic=chnldic["Chantoplev"]
    chbotdic=chnldic["Chanbotlev"]
    chwavldic=chnldic["Chanwavlen"]
    chwluntfctr=chnldic["Chwlunitfctr"]
    spdlght=300000000
    for chaname in chlist:
        if (lev > chtopdic[chaname] and lev < chbotdic[chaname]) or lev == chbotdic[chaname]:
            wavlen=chwavldic[chaname]
            freq=spdlght*chwluntfctr/wavlen
    return(freq)

def dateindex(datetime,nctime,units):                                            
    num = int(netCDF4.date2num(datetime, units = units)   )
    num_list = nctime[:].tolist()
    idx = num_list.index(float(num))
    return idx

def indexdate(idx,nctime, units):
    return netCDF4.num2date(nctime[idx], units)

def find_nearest_time(nctime, DT, units):
    timevalue = netCDF4.date2num(DT, units = units)
    timelist = nctime[:].tolist()
    idx = find_nearest_index(timelist, timevalue)
    return(idx)

def find_nearest_index(a, value):
    array = numpy.asarray(a)
    idx = (numpy.abs(array - value)).argmin()
    return idx

def convert_mb_to_pa(a):
    if isinstance(a, (list,)):
        return([int(i)*100 for i in a])
    else:
        return(a*100)

def get_hloswnd(uwnd,vwnd,azmh):
    fff=get_fff(uwnd,vwnd)
    ddd=get_ddd(uwnd,vwnd)
    hloswnd=fff*math.cos(math.radians(ddd-azmh))*(-1)
    return(hloswnd)

def get_fff(uwnd,vwnd):
    return math.sqrt((uwnd*uwnd)+(vwnd*vwnd))

def get_ddd(uwnd,vwnd):
    ddd=math.degrees(math.atan(uwnd/vwnd))
    if vwnd > 0 : ddd += 180
    if ddd < 0 : ddd += 360
    return ddd

def get_xangle(uwnd,vwnd):
    xangle=math.degrees(math.atan(vwnd/uwnd))
    if uwnd < 0 : xangle += 180
    if xangle < 0 : xangle += 360
    return xangle

def get_uwnd(fff,ddd):
    return (-1.0*fff*math.sin(math.radians(ddd)))

def get_vwnd(fff,ddd):
    return (-1.0*fff*math.cos(math.radians(ddd)))

def get_dpt_onrhum(temp,rhum):
    LperRv=5423
    if rhum > 0:
        return (1.0/((1.0/temp)-(1.0/LperRv)*math.log(rhum/100)))
    else:
        return None

def get_dpt_onshum(shum):
    LperRv=5423
    T0=273
    if shum > 0:
        return (1.0/((1.0/T0)-(1.0/LperRv)*math.log(shum)))
    else:
        return None

def get_rhum(temp,dpt):
    LperRv=5423
    return ((math.exp(LperRv*((1.0/temp)-(1.0/dpt))))*100)

def conform_dims(data,axis,dims):
    if axis is "cell":
        return pandas.DataFrame([[data]*dims[1]]*dims[0], index=range(1,dims[0]+1,1), columns=range(1,dims[1]+1,1))
    else:
        return{
            'cell'      : pandas.DataFrame([[data]*dims[1]]*dims[0]),
            'index'     : pandas.DataFrame(numpy.array([[data[i]]*dims[1] for i in range(0,len(data),1)])),
            'columns'   : pandas.DataFrame([data]*dims[0]),
            'transpos'  : pandas.DataFrame(data[i].values for i in range(0,len(data.xs(0)),1))
            }.get(axis, pandas.DataFrame([data]*dims[0]*dims[1]))

def extend_dims(data,axis,length):
    if type(data) is int:
        return {
            'index'     : pandas.DataFrame([data]*length,index=range(1,length+1,1)),
            'columns'   : pandas.DataFrame([[data]*length],index=[1],columns=range(1,length+1,1))
            }.get(axis, pandas.DataFrame([data]))
    else:
        return {
            'index'     : pandas.DataFrame([data]*length,index=range(1,length+1,1)),
            'columns'   : pandas.DataFrame(numpy.array([[data[i]]*length for i in range(0,len(data),1)]))
            }.get(axis, pandas.DataFrame([data]))

def generate2Dindex(obs_cnt,LTTD,LNGD,datetime,time,time_units,lat,lon):
    itime=dateindex(datetime,time,time_units)
    for i in range(0,obs_cnt,1):
        ilat=find_nearest_index(lat,LTTD.values[i])
        ilon=find_nearest_index(lon,LNGD.values[i])
        if i in [0]: 
            idx=pandas.DataFrame([[itime,ilat,ilon]],index=[i+1],columns=["itime","ilat","ilon"])
        else:
            idx=idx.append(pandas.DataFrame([[itime,ilat,ilon]],index=[i+1],columns=["itime","ilat","ilon"]))
    return(idx)

def generate3Dindex(obs_cnt,LTTD,LNGD,ALTTD,datetime,time,time_units,lat,lon,lev,gph):
    itime=dateindex(datetime,time,time_units)
    for i in range(0,obs_cnt,1):
        ilat=find_nearest_index(lat,LTTD.values[i])
        ilon=find_nearest_index(lon,LNGD.values[i])
        ilev=find_nearest_index(gph[itime,:,ilat,ilon],ALTTD.values[i])
        if i in [0]: 
            idx=pandas.DataFrame([[itime,ilat,ilon,ilev]],index=[i+1],columns=["itime","ilat","ilon","ilev"])
        else:
            idx=idx.append(pandas.DataFrame([[itime,ilat,ilon,ilev]],index=[i+1],columns=["itime","ilat","ilon","ilev"]))
    return(idx)

def HLOSWind(location,findex,Tstart,outpath,nmlfile,obscount,maxindx=MAXINDX):   #,obsgroup,subtype,elist,data,batchcount=1,header_offset=339,maxindx=600,lut_ncols=LUTSIZE):
    data=location.join(findex.lev,on='ProfileNo')
    Year=location.Year.values[0]
    uwnd=ncepradic.getdata(Year,"uwnd",element="uwnd")
    vwnd=ncepradic.getdata(Year,"vwnd",element="vwnd")
    temp=ncepradic.getdata(Year,"tmp",element="tmp")
    ########################################################
    data=prepdata(data,"ProfileNo","SatID",obscount,findex,option="const",constdic=obsdic.HLOSWind)
    data=prepdata(data,"ProfileNo","SatZenithAngle",obscount,findex,"const",constdic=obsdic.HLOSWind)
    data=prepdata(data,"ProfileNo","PlevelsA",obscount,findex,option="PlevelsA")
    data=prepdata(data,"ProfileNo","AltPressLevel",obscount,findex,option="PlevelsA")
    data=prepdata(data,"ProfileNo","fff",obscount,findex,option="fff",uwnd=uwnd,vwnd=vwnd)
    data=prepdata(data,"ProfileNo","CmptVecWindSpeed",obscount,findex,option="fff",uwnd=uwnd,vwnd=vwnd)
    data=prepdata(data,"ProfileNo","ddd",obscount,findex,option="ddd",uwnd=uwnd,vwnd=vwnd)
    data=prepdata(data,"ProfileNo","CmptVecWindDirn",obscount,findex,option="ddd",uwnd=uwnd,vwnd=vwnd)
    data=prepdata(data,"ProfileNo","HLOSWIND",obscount,findex,option="hloswnd",uwnd=uwnd,vwnd=vwnd)
    data=prepdata(data,"ProfileNo","AltAirTempLevel",obscount,findex,option="temp",temp=temp)
    ########################################################
    filename=obsdic.HLOSWind["filename"]
    output_file="%s/%s" % (outpath,filename)
    obs_index=obsdic.HLOSWind["obs_index"]
    elist=obstore.obstore_create_element_table(nmlfile,obs_index)
    obsgroup=obsdic.HLOSWind["obsgroup"]
    subtype=obsdic.HLOSWind["subtype"]
    with open(output_file, "wb+") as outfile:
        obstore.create_obstore(Tstart,outfile,nmlfile,obsgroup,subtype,elist,data,batchcount=1,header_offset=339,maxindx=600,lut_ncols=LUTSIZE)
    #difflist=obstore.header_diffcheck(output_file)
    #print(difflist)
    batchinfo=obstore.batch_header_read(output_file,nmlfile,maxindx=maxindx,batchid=1)
    #print(batchinfo)
    return(data)

def LEOGEOAMV(location,findex,Tstart,nmlfile,maxindx=MAXINDX):   #,obsgroup,subtype,elist,data,batchcount=1,header_offset=339,maxindx=600,lut_ncols=LUTSIZE):
    data=location.join(findex.lev,on='ProfileNo')
    Year=location.Year.values[0]
    Month=location.Month.values[0]
    Day=location.Day.values[0]
    Hour=location.Hour.values[0]
    Minute=location.Minute.values[0]
    uwnd=ncepradic.getdata(Year,"uwnd",element="uwnd")
    vwnd=ncepradic.getdata(Year,"vwnd",element="vwnd")
    temp=ncepradic.getdata(Year,"tmp",element="tmp")
    subdic=obsdic.LEOGEOAMV
#    obs_index=subdic["obs_index"]
#    elist=obstore.obstore_create_element_table(nmlfile,obs_index)
#    elenams=elist.Element.values
    ########################################################
    data=prepdata(subdic,data,findex,"ProfileNo","ReceiptYear",option="const",constdic={"ReceiptYear" : Year})
    data=prepdata(subdic,data,findex,"ProfileNo","ReceiptMonth",option="const",constdic={"ReceiptMonth" : Month})
    data=prepdata(subdic,data,findex,"ProfileNo","ReceiptDay",option="const",constdic={"ReceiptDay" : Day})
    data=prepdata(subdic,data,findex,"ProfileNo","ReceiptHour",option="const",constdic={"ReceiptHour" : Hour})
    data=prepdata(subdic,data,findex,"ProfileNo","ReceiptMinute",option="const",constdic={"ReceiptMinute" : Minute})
    data=prepdata(subdic,data,findex,"ProfileNo","SatID",option="const",constdic=obsdic.LEOGEOAMV)
    data=prepdata(subdic,data,findex,"ProfileNo","OrigCtr",option="const",constdic=obsdic.LEOGEOAMV)
    data=prepdata(subdic,data,findex,"ProfileNo","CloudMotionMethod",option="const",constdic=obsdic.LEOGEOAMV)
    data=prepdata(subdic,data,findex,"ProfileNo","HeightAssMethod",option="const",constdic=obsdic.LEOGEOAMV)
    data=prepdata(subdic,data,findex,"ProfileNo","AltHeightAss",option="const",constdic=obsdic.LEOGEOAMV)
    data=prepdata(subdic,data,findex,"ProfileNo","SatZenithAngle",option="const",constdic=obsdic.LEOGEOAMV)
    data=prepdata(subdic,data,findex,"ProfileNo","PlevelsA",option="PlevelsA")
    data=prepdata(subdic,data,findex,"ProfileNo","AltPressLevel",option="PlevelsA")
    data=prepdata(subdic,data,findex,"ProfileNo","fff",option="fff",uwnd=uwnd,vwnd=vwnd)
    data=prepdata(subdic,data,findex,"ProfileNo","CmptVecWindSpeed",option="fff",uwnd=uwnd,vwnd=vwnd)
    data=prepdata(subdic,data,findex,"ProfileNo","ddd",option="ddd",uwnd=uwnd,vwnd=vwnd)
    data=prepdata(subdic,data,findex,"ProfileNo","CmptVecWindDirn",option="ddd",uwnd=uwnd,vwnd=vwnd)
    data=prepdata(subdic,data,findex,"ProfileNo","HLOSWIND",option="hloswnd",uwnd=uwnd,vwnd=vwnd)
    data=prepdata(subdic,data,findex,"ProfileNo","ChanCtralFreq","chnlfreq",constdic=obsdic.LEOGEOAMV)
    data=prepdata(subdic,data,findex,"ProfileNo","Temp",option="temp",temp=temp)
    data=prepdata(subdic,data,findex,"ProfileNo","AltAirTempLevel",option="temp",temp=temp)
    ########################################################
#    filename=obsdic.LEOGEOAMV["filename"]
#    output_file="%s/%s" % (outpath,filename)
#    obsgroup=obsdic.LEOGEOAMV["obsgroup"]
#    subtype=obsdic.LEOGEOAMV["subtype"]
#    obs_index=obsdic.LEOGEOAMV["obs_index"]
#    elist=obstore.obstore_create_element_table(nmlfile,obs_index)
#    with open(output_file, "wb+") as outfile:
#        (datapos,datalen,dataend)=obstore.create_obstore(Tstart,outfile,nmlfile,obsgroup,[subtype],[elist],[data],batchcount=1,header_offset=339,maxindx=maxindx,lut_ncols=LUTSIZE)
#    print(datapos,datalen,dataend)
    #difflist=obstore.header_diffcheck(output_file)
    #print(difflist)
    #batchinfo=obstore.batch_header_read(output_file,nmlfile,maxindx=600,batchid=1)
    #print(batchinfo)
    return(data)

def aladin_profile(Tnode,Tstart,Tstop,output_file,nmlpath,maxindx=MAXINDX):
    obs_index_nml="obs_index_nml"
    nmlfile="%s/%s" % (nmlpath,obs_index_nml)
    SatID=obsdic.Aeolus["SatID"]
    SMAxis=obsdic.Aeolus["SMAxis"]
    SatAlt=obsdic.Aeolus["SatAlt"]
    OrbInc=obsdic.Aeolus["OrbInc"]
    HVA=obsdic.Aladin["HVA"]
    VVA=obsdic.Aladin["VVA"]
    ObsAltList=obsdic.Aladin["ObsAltList"]
    profsec=obsdic.Aladin["profsec"]
    rptcyc=obsdic.Aeolus["repeat"]
    Torbit=obsdic.Aeolus["Torbit"]
    Tref=obsdic.Aeolus["Tref"]
    NodeShift=obslib.orbitalshift(Tnode,Tref,Torb=Torbit,Rcyc=rptcyc)
    Tprofile=obslib.timeperiod(ss=obsdic.Aladin["profsec"])
    Tobs=obslib.timeperiod(ss=obsdic.Aladin["obsec"])
    pid_start=obslib.profileid(Tstart,Tnode,profsec)
    pid_stop=obslib.profileid(Tstop,Tnode,profsec)
    ProfileIDs=numpy.arange(pid_start,pid_stop,1)
    for prfcount,pid in enumerate(ProfileIDs):                
        Tnow=Tnode+pid*Tprofile
        
        nadir=nadir_latlon(Tnow,Tnode,Torbit,OrbInc,NodeShift)
        nadirlon=float(nadir.Longitude.values[0])
        nadirlat=float(nadir.Latitude.values[0])
        LocalTime=obslib.localtime(Tnow,nadirlon)
        prf_start=prfcount*len(ObsAltList)
        for lev,ObsAltSI in enumerate(ObsAltList,1):
            ObsAlt=ObsAltSI/1000
            obscount=prf_start+lev
            if obscount in [1]:
                location=los_latlon(SatID, obscount, ObsAlt,Tnow,Tnode,Torbit,SMAxis,SatAlt,nadir,OrbInc,HVA,VVA)
                #findex=fileindex(location,obscount,filedim,gph,altnam="HeightCOG")
            else:
                location=location.append(los_latlon(SatID, obscount, ObsAlt,Tnow,Tnode,Torbit,SMAxis,SatAlt,nadir,OrbInc,HVA,VVA))
                #findex=findex.append(fileindex(location,obscount,filedim,gph,altnam="HeightCOG"))
    return(location)

def symulate_data(datetime,field,lat,lon,lev,time,time_units,obs_cnt,LTTD,LNGD):
    itime=dateindex(datetime,time,time_units)
    data=pandas.DataFrame(index=range(1,obs_cnt+1,1),columns=convert_mb_to_pa(lev))
    for i in range(0,obs_cnt,1):
        ilat=find_nearest_index(lat,LTTD.values[i])
        ilon=find_nearest_index(lon,LNGD.values[i])
        data.xs(i+1)[:]=field[itime,:,ilat,ilon]
    return(data)
    
def symulate_fff(datetime,uwnd,vwnd,lat,lon,lev,time,time_units,obs_cnt,LTTD,LNGD):
    idx=generate2Dindex(obs_cnt,LTTD,LNGD,datetime,time,time_units,lat,lon,lev)
    fff=pandas.DataFrame(index=range(1,obs_cnt+1,1),columns=convert_mb_to_pa(lev))
    for i in range(0,obs_cnt,1):
        for j,lvl in enumerate(lev):
            fff.xs(i+1)[j]=get_fff(uwnd[idx.itime[i],[j],idx.ilat[i],idx.ilon[i]],vwnd[idx.itime[i],[j],idx.ilat[i],idx.ilon[i]])
    return(fff)

def symulate_ddd(datetime,uwnd,vwnd,lat,lon,lev,time,time_units,obs_cnt,LTTD,LNGD):
    idx=generate2Dindex(obs_cnt,LTTD,LNGD,datetime,time,time_units,lat,lon)
    ddd=pandas.DataFrame(index=range(1,obs_cnt+1,1),columns=convert_mb_to_pa(lev))
    for i in range(0,obs_cnt,1):
        for j,lvl in enumerate(lev):
            ddd.xs(i+1)[j]=get_ddd(uwnd[idx.itime[i],[j],idx.ilat[i],idx.ilon[i]],vwnd[idx.itime[i],[j],idx.ilat[i],idx.ilon[i]])
    return(ddd)

def symulate_dpt_onrhum(datetime,temp,rhum,lat,lon,lev,time,time_units,obs_cnt,LTTD,LNGD):
    idx=generate2Dindex(obs_cnt,LTTD,LNGD,datetime,time,time_units,lat,lon)
    dpt=pandas.DataFrame(index=range(1,obs_cnt+1,1),columns=convert_mb_to_pa(lev))
    for i in range(0,obs_cnt,1):
        for j,lvl in enumerate(lev):
            dpt.xs(i+1)[j]=get_dpt_onrhum(temp[idx.itime[i],[j],idx.ilat[i],idx.ilon[i]],rhum[idx.itime[i],[j],idx.ilat[i],idx.ilon[i]])
    return(dpt)
    
def symulate_dpt_onshum(datetime,shum,lat,lon,lev,time,time_units,obs_cnt,LTTD,LNGD):
    idx=generate2Dindex(obs_cnt,LTTD,LNGD,datetime,time,time_units,lat,lon)
    dpt=pandas.DataFrame(index=range(1,obs_cnt+1,1),columns=convert_mb_to_pa(lev))
    for i in range(0,obs_cnt,1):
        for j,lvl in enumerate(lev):
            dpt.xs(i+1)[j]=get_dpt_onshum(shum[idx.itime[i],[j],idx.ilat[i],idx.ilon[i]])
    return(dpt)


def symulate_data(outfile,nmlfile,indx,element,fvidx,filevar,datetime,obs_cnt,LTTD,LNGD,option="data"):
    subtype=obstore.obstore_read_index_subtype(outfile,indx)
    elist = obstore.obstore_read_batch_elements(outfile,indx,nmlfile)
    if element in elist.Element.values:
        if diaglev > 0 : print(element)
        (rdc,cdc,ldc) = elist.query("Element == @element")[["RDC","CDC","LDC"]].reset_index(drop=True).values[0]
        if diaglev > 1 : print(element,rdc,cdc,ldc)
        if int(ldc) is 1: symulate_data_2dfield(outfile,nmlfile,indx,element,fvidx,filevar,datetime,obs_cnt,LTTD,LNGD)
        else: symulate_data_3dfield(outfile,nmlfile,indx,element,fvidx,filevar,datetime,obs_cnt,LTTD,LNGD)

def symulate_data_flightlev(outfile,nmlfile,batchindx,element,fvidx,filevar,datetime,obs_cnt,LTTD,LNGD,ALTTD,option="data"):
    subtype=obstore.obstore_read_index_subtype(outfile,batchindx)
    elist = obstore.obstore_read_batch_elements(outfile,batchindx,nmlfile)
    if element in elist.Element.values:
        with netCDF4.Dataset(filevar.file[fvidx], 'r',  format='NETCDF4_CLASSIC') as ncf, \
        netCDF4.Dataset(filevar.file["gph"], 'r',  format='NETCDF4_CLASSIC') as gphfile :
            gph = gphfile.variables[filevar.fvname["gph"]]
            filedata = ncf.variables[filevar.fvname[fvidx]]
            if diaglev > 20 : print(filedata)
            time=ncf.variables[filevar.fvname["time"]][:]
            time_units=ncf.variables[filevar.fvname["time"]].units
            lev=ncf.variables[filevar.fvname["lev"]][:]
            lat=ncf.variables[filevar.fvname["lat"]][:]
            lon=ncf.variables[filevar.fvname["lon"]][:]
            idx=generate3Dindex(obs_cnt,LTTD,LNGD,ALTTD,datetime,time,time_units,lat,lon,lev,gph)
            data=get_data(element,obs_cnt,idx,option=option,lev=lev,filedata=filedata,temp=[],rhum=[],uwnd=[],vwnd=[])
#            for i in range(0,obs_cnt,1):
#                idata=filedata[indx.itime[i+1],indx.ilev[i+1],indx.ilat[i+1],indx.ilon[i+1]]
#                if i in [0]: 
#                    data=pandas.DataFrame([[idata]],index=[i+1],columns=[element])
#                else:
#                    data=data.append(pandas.DataFrame([[idata]],index=[i+1],columns=[element]))
        obstore.obstore_write_data_element(outfile,nmlfile,batchindx,element,data)
        return(data)

def symulate_data_2dfield(outfile,nmlfile,batchindx,element,fvidx,filevar,datetime,obs_cnt,LTTD,LNGD,option="data"):
    with netCDF4.Dataset(filevar.file[fvidx], 'r',  format='NETCDF4_CLASSIC') as ncf:
        filedata = ncf.variables[filevar.fvname[fvidx]]
        if diaglev > 20 : print(filedata)
        time=ncf.variables[filevar.fvname["time"]][:]
        time_units=ncf.variables[filevar.fvname["time"]].units
        lat=ncf.variables[filevar.fvname["lat"]][:]
        lon=ncf.variables[filevar.fvname["lon"]][:]
        idx=generate2Dindex(obs_cnt,LTTD,LNGD,datetime,time,time_units,lat,lon)
        data=get_data(element,obs_cnt,idx,option=option,filedata=filedata,temp=[],rhum=[],uwnd=[],vwnd=[])
#        for i in range(0,obs_cnt,1):
#                idata=filedata[indx.itime[i+1],indx.ilat[i+1],indx.ilon[i+1]]
#                if i in [0]: 
#                    data=pandas.DataFrame([[idata]],index=[i+1],columns=[element])
#                else:
#                    data=data.append(pandas.DataFrame([[idata]],index=[i+1],columns=[element]))
#        data=pandas.DataFrame(index=range(1,obs_cnt+1,1),columns=[element])
#        for i in range(1,obs_cnt+1,1):
#            data.xs(i)[element]=filedata[idx.itime[i],idx.ilat[i],idx.ilon[i]]
    obstore.obstore_write_data_element(outfile,nmlfile,batchindx,element,data)
    return(data)
        
def symulate_data_3dfield(outfile,nmlfile,batchindx,element,fvidx,filevar,datetime,obs_cnt,LTTD,LNGD,option="data"):
    with netCDF4.Dataset(filevar.file[fvidx], 'r',  format='NETCDF4_CLASSIC') as ncf:
        filedata = ncf.variables[filevar.fvname[fvidx]]
        if diaglev > 20 : print(filedata)
        time=ncf.variables[filevar.fvname["time"]][:]
        time_units=ncf.variables[filevar.fvname["time"]].units
        lev=ncf.variables[filevar.fvname["lev"]][:]
        lat=ncf.variables[filevar.fvname["lat"]][:]
        lon=ncf.variables[filevar.fvname["lon"]][:]
        idx=generate2Dindex(obs_cnt,LTTD,LNGD,datetime,time,time_units,lat,lon)
        data=get_data(element,obs_cnt,idx,option=option,lev=lev,filedata=filedata,temp=[],rhum=[],uwnd=[],vwnd=[])
#        for i in range(0,obs_cnt,1):
#                idata=filedata[indx.itime[i+1],:,indx.ilat[i+1],indx.ilon[i+1]]
#                if i in [0]: 
#                    data=pandas.DataFrame([idata],index=[i+1],columns=convert_mb_to_pa(lev))
#                else:
#                    data=data.append(pandas.DataFrame([idata],index=[i+1],columns=convert_mb_to_pa(lev)))
#        data=pandas.DataFrame(index=range(1,obs_cnt+1,1),columns=convert_mb_to_pa(lev))
#        for i in range(0,obs_cnt,1):
#            data.xs(i+1)[:]=filedata[idx.itime[i],:,idx.ilat[i],idx.ilon[i]]
    obstore.obstore_write_data_element(outfile,nmlfile,batchindx,element,data)
    return(data)

def symulate_wind_2dfield(outfile,nmlfile,batchindx,element,fvidx,filevar,datetime,obs_cnt,LTTD,LNGD,option="data"):
    with netCDF4.Dataset(filevar.file[fvidx[0]], 'r',  format='NETCDF4_CLASSIC') as ufile, \
        netCDF4.Dataset(filevar.file[fvidx[1]], 'r',  format='NETCDF4_CLASSIC') as vfile:
        uwnd = ufile.variables[filevar.fvname[fvidx[0]]]
        vwnd = vfile.variables[filevar.fvname[fvidx[1]]]
        if diaglev > 20 : print(uwnd,vwnd)
        time=ufile.variables[filevar.fvname["time"]][:]
        time_units=ufile.variables[filevar.fvname["time"]].units
        lat=ufile.variables[filevar.fvname["lat"]][:]
        lon=ufile.variables[filevar.fvname["lon"]][:]
        idx=generate2Dindex(obs_cnt,LTTD,LNGD,datetime,time,time_units,lat,lon)
        data=get_data(element,obs_cnt,idx,option=option,filedata=[],temp=[],rhum=[],uwnd=uwnd,vwnd=vwnd)
#        for i in range(0,obs_cnt,1):
#                if option in ["uwnd"] : idata=uwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]]
#                if option in ["vwnd"] : idata=vwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]]
#                if option in ["ddd"] : idata=get_ddd(uwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]])
#                if option in ["fff"] : idata=get_fff(uwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]])
#                if i in [0]: 
#                    data=pandas.DataFrame([[idata]],index=[i+1],columns=[element])
#                else:
#                    data=data.append(pandas.DataFrame([[idata]],index=[i+1],columns=[element]))
#        data=pandas.DataFrame(index=range(1,obs_cnt+1,1),columns=[element])
#        for i in range(0,obs_cnt,1):
#            data.xs(i+1)[element]=get_ddd(uwnd[idx.itime[i],idx.ilat[i],idx.ilon[i]],vwnd[idx.itime[i],idx.ilat[i],idx.ilon[i]])
    obstore.obstore_write_data_element(outfile,nmlfile,batchindx,element,data)
    return(data)

def symulate_wind_3dfield(outfile,nmlfile,indx,element,fvidx,filevar,datetime,obs_cnt,LTTD,LNGD,option="data"):
    with netCDF4.Dataset(filevar.file[fvidx[0]], 'r',  format='NETCDF4_CLASSIC') as ufile, \
        netCDF4.Dataset(filevar.file[fvidx[1]], 'r',  format='NETCDF4_CLASSIC') as vfile:
        uwnd = ufile.variables[filevar.fvname[fvidx[0]]]
        vwnd = vfile.variables[filevar.fvname[fvidx[1]]]
        if diaglev > 20 : print(uwnd,vwnd)
        time=ufile.variables[filevar.fvname["time"]][:]
        time_units=ufile.variables[filevar.fvname["time"]].units
        lat=ufile.variables[filevar.fvname["lat"]][:]
        lon=ufile.variables[filevar.fvname["lon"]][:]
        lev=ufile.variables[filevar.fvname["lev"]][:]
        idx=generate2Dindex(obs_cnt,LTTD,LNGD,datetime,time,time_units,lat,lon)
        data=get_data(element,obs_cnt,idx,option=option,lev=lev,filedata=[],temp=[],rhum=[],uwnd=uwnd,vwnd=vwnd)
#        for i in range(0,obs_cnt,1):
#            idata=None*len(lev)
#            for j,lvl in enumerate(lev):
#                if element in ["uwnd"] : idata[j]=uwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]]
#                if element in ["vwnd"] : idata[j]=vwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]]
#                if element in ["ddd"] : idata[j]=get_ddd(uwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]])
#                if element in ["fff"] : idata[j]=get_fff(uwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]])
#            if i in [0]: 
#                data=pandas.DataFrame([idata],index=[i+1],columns=convert_mb_to_pa(lev))
#            else:
#                data=data.append(pandas.DataFrame([idata],index=[i+1],columns=convert_mb_to_pa(lev)))
        ########################################
#        data=pandas.DataFrame(index=range(1,obs_cnt+1,1),columns=convert_mb_to_pa(lev))
#        for i in range(0,obs_cnt,1):
#            for j,lvl in enumerate(lev):
#                data.xs(i+1)[j]=get_ddd(uwnd[idx.itime[i],[j],idx.ilat[i],idx.ilon[i]],vwnd[idx.itime[i],[j],idx.ilat[i],idx.ilon[i]])
    obstore.obstore_write_data_element(outfile,nmlfile,indx,element,data)
    return(data)
    
def symulate_wind_flightlev(outfile,nmlfile,batchindx,element,fvidx,filevar,datetime,obs_cnt,LTTD,LNGD,ALTTD,option="data"):
    subtype=obstore.obstore_read_index_subtype(outfile,batchindx)
    elist = obstore.obstore_read_batch_elements(outfile,batchindx,nmlfile)
    if element in elist.Element.values:
        with netCDF4.Dataset(filevar.file[fvidx[0]], 'r',  format='NETCDF4_CLASSIC') as ufile, \
        netCDF4.Dataset(filevar.file[fvidx[1]], 'r',  format='NETCDF4_CLASSIC') as vfile, \
        netCDF4.Dataset(filevar.file["gph"], 'r',  format='NETCDF4_CLASSIC') as gphfile :
            gph = gphfile.variables[filevar.fvname["gph"]]
            uwnd = ufile.variables[filevar.fvname[fvidx[0]]]
            vwnd = vfile.variables[filevar.fvname[fvidx[1]]]
            if diaglev > 20 : print(uwnd,vwnd)
            time=ufile.variables[filevar.fvname["time"]][:]
            time_units=ufile.variables[filevar.fvname["time"]].units
            lat=ufile.variables[filevar.fvname["lat"]][:]
            lon=ufile.variables[filevar.fvname["lon"]][:]
            lev=ufile.variables[filevar.fvname["lev"]][:]
            idx=generate3Dindex(obs_cnt,LTTD,LNGD,ALTTD,datetime,time,time_units,lat,lon,lev,gph)
            data=get_data(element,obs_cnt,idx,option=option,lev=lev,filedata=[],temp=[],rhum=[],uwnd=uwnd,vwnd=vwnd)
#            for i in range(0,obs_cnt,1):
#                if element in ["uwnd"] : idata=uwnd[indx.itime[i+1],indx.ilev[i+1],indx.ilat[i+1],indx.ilon[i+1]]
#                if element in ["vwnd"] : idata=vwnd[indx.itime[i+1],indx.ilev[i+1],indx.ilat[i+1],indx.ilon[i+1]]
#                if element in ["ddd"] : idata=get_ddd(uwnd[indx.itime[i+1],indx.ilev[i+1],indx.ilat[i+1],indx.ilon[i+1]],vwnd[indx.itime[i+1],indx.ilev[i+1],indx.ilat[i+1],indx.ilon[i+1]])
#                if element in ["fff"] : idata=get_fff(uwnd[indx.itime[i+1],indx.ilev[i+1],indx.ilat[i+1],indx.ilon[i+1]],vwnd[indx.itime[i+1],indx.ilev[i+1],indx.ilat[i+1],indx.ilon[i+1]])
#                if i in [0]: 
#                    data=pandas.DataFrame([[idata]],index=[i+1],columns=[element])
#                else:
#                    data=data.append(pandas.DataFrame([[idata]],index=[i+1],columns=[element]))
        obstore.obstore_write_data_element(outfile,nmlfile,batchindx,element,data)
        return(data)

#def symulate_fff_2dfield(outfile,nmlfile,indx,element,fvidx,filevar,datetime,obs_cnt,LTTD,LNGD):
#    with netCDF4.Dataset(filevar.file[fvidx[0]], 'r',  format='NETCDF4_CLASSIC') as ufile, \
#        netCDF4.Dataset(filevar.file[fvidx[1]], 'r',  format='NETCDF4_CLASSIC') as vfile:
#        uwnd = ufile.variables[filevar.fvname[fvidx[0]]]
#        vwnd = vfile.variables[filevar.fvname[fvidx[1]]]
#        if diaglev > 20 : print(uwnd,vwnd)
#        time=ufile.variables[filevar.fvname["time"]][:]
#        time_units=ufile.variables[filevar.fvname["time"]].units
#        lat=ufile.variables[filevar.fvname["lat"]][:]
#        lon=ufile.variables[filevar.fvname["lon"]][:]
#        idx=generate2Dindex(obs_cnt,LTTD,LNGD,datetime,time,time_units,lat,lon)
#        data=pandas.DataFrame(index=range(1,obs_cnt+1,1),columns=[element])
#        for i in range(0,obs_cnt,1):
#            data.xs(i+1)[element]=get_fff(uwnd[idx.itime[i],idx.ilat[i],idx.ilon[i]],vwnd[idx.itime[i],idx.ilat[i],idx.ilon[i]])
#    obstore.obstore_write_data_element(outfile,nmlfile,indx,element,data)
#    return(data)
#
#def symulate_fff_3dfield(outfile,nmlfile,indx,element,fvidx,filevar,datetime,obs_cnt,LTTD,LNGD):
#    with netCDF4.Dataset(filevar.file[fvidx[0]], 'r',  format='NETCDF4_CLASSIC') as ufile, \
#        netCDF4.Dataset(filevar.file[fvidx[1]], 'r',  format='NETCDF4_CLASSIC') as vfile:
#        uwnd = ufile.variables[filevar.fvname[fvidx[0]]]
#        vwnd = vfile.variables[filevar.fvname[fvidx[1]]]
#        if diaglev > 20 : print(uwnd,vwnd)
#        time=ufile.variables[filevar.fvname["time"]][:]
#        time_units=ufile.variables[filevar.fvname["time"]].units
#        lat=ufile.variables[filevar.fvname["lat"]][:]
#        lon=ufile.variables[filevar.fvname["lon"]][:]
#        lev=ufile.variables[filevar.fvname["lev"]][:]
#        idx=generate2Dindex(obs_cnt,LTTD,LNGD,datetime,time,time_units,lat,lon)
#        data=pandas.DataFrame(index=range(1,obs_cnt+1,1),columns=convert_mb_to_pa(lev))
#        for i in range(0,obs_cnt,1):
#            for j,lvl in enumerate(lev):
#                data.xs(i+1)[j]=get_fff(uwnd[idx.itime[i],[j],idx.ilat[i],idx.ilon[i]],vwnd[idx.itime[i],[j],idx.ilat[i],idx.ilon[i]])
#    obstore.obstore_write_data_element(outfile,nmlfile,indx,element,data)
#    return(data)

def symulate_dpt_2dfield(outfile,nmlfile,batchindx,element,fvidx,filevar,datetime,obs_cnt,LTTD,LNGD):
    if "shum" in fvidx:
        option = "dpt_onshum"
        with netCDF4.Dataset(filevar.file["shum"], 'r',  format='NETCDF4_CLASSIC') as ncf:
            filedata = ncf.variables[filevar.fvname["shum"]]
            if diaglev > 20 : print(filedata)
            time=ncf.variables[filevar.fvname["time"]][:]
            time_units=ncf.variables[filevar.fvname["time"]].units
            lat=ncf.variables[filevar.fvname["lat"]][:]
            lon=ncf.variables[filevar.fvname["lon"]][:]
            idx=generate2Dindex(obs_cnt,LTTD,LNGD,datetime,time,time_units,lat,lon)
            data=get_data(element,obs_cnt,idx,option=option,filedata=filedata,temp=[],rhum=[],uwnd=[],vwnd=[])
#            for i in range(0,obs_cnt,1):
#                if option in ["uwnd"] : idata=uwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]]
#                if option in ["vwnd"] : idata=vwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]]
#                if option in ["ddd"] : idata=get_ddd(uwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]])
#                if option in ["fff"] : idata=get_fff(uwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]])
#                if option in ["dpt_onshum"] : idata=get_dpt_onshum(filedata[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]])
#                if option in ["dpt_onrhum"] : idata=get_dpt_onrhum(temp[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]],rhum[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]])
#                if i in [0]: 
#                    data=pandas.DataFrame([[idata]],index=[i+1],columns=[element])
#                else:
#                    data=data.append(pandas.DataFrame([[idata]],index=[i+1],columns=[element]))
#            data=pandas.DataFrame(index=range(1,obs_cnt+1,1),columns=[element])
#            for i in range(0,obs_cnt,1):
#                    data.xs(i+1)[element]=get_dpt_onshum(filedata[idx.itime[i],idx.ilat[i],idx.ilon[i]])
    else:
        option = "dpt_onrhum"
        with netCDF4.Dataset(filevar.file["rhum"], 'r',  format='NETCDF4_CLASSIC') as rhfile, \
            netCDF4.Dataset(filevar.file["tmp"], 'r',  format='NETCDF4_CLASSIC') as tfile:
            temp = tfile.variables[filevar.fvname["tmp"]]
            rhum=rhfile.variables[filevar.fvname["rhum"]]
            if diaglev > 20 : print(temp,rhum)
            time=rhfile.variables[filevar.fvname["time"]][:]
            time_units=rhfile.variables[filevar.fvname["time"]].units
            lat=rhfile.variables[filevar.fvname["lat"]][:]
            lon=rhfile.variables[filevar.fvname["lon"]][:]
            idx=generate2Dindex(obs_cnt,LTTD,LNGD,datetime,time,time_units,lat,lon)
            data=get_data(element,obs_cnt,idx,option=option,filedata=filedata,temp=[],rhum=[],uwnd=[],vwnd=[])
#            for i in range(0,obs_cnt,1):
#                if option in ["uwnd"] : idata=uwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]]
#                if option in ["vwnd"] : idata=vwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]]
#                if option in ["ddd"] : idata=get_ddd(uwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]])
#                if option in ["fff"] : idata=get_fff(uwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]])
#                if option in ["dpt_onshum"] : idata=get_dpt_onshum(filedata[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]])
#                if option in ["dpt_onrhum"] : idata=get_dpt_onrhum(temp[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]],rhum[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]])
#                if i in [0]: 
#                    data=pandas.DataFrame([[idata]],index=[i+1],columns=[element])
#                else:
#                    data=data.append(pandas.DataFrame([[idata]],index=[i+1],columns=[element]))
#            data=pandas.DataFrame(index=range(1,obs_cnt+1,1),columns=convert_mb_to_pa(lev))
#            for i in range(0,obs_cnt,1):
#                data.xs(i+1)[element]=get_dpt_onrhum(temp[idx.itime[i],idx.ilat[i],idx.ilon[i]],rhum[idx.itime[i],idx.ilat[i],idx.ilon[i]])
    obstore.obstore_write_data_element(outfile,nmlfile,batchindx,element,data)
    return(data)


def symulate_dpt_3dfield(outfile,nmlfile,indx,element,fvidx,filevar,datetime,obs_cnt,LTTD,LNGD):
    if "shum" in fvidx:
        option = "dpt_onshum"
        with netCDF4.Dataset(filevar.file["shum"], 'r',  format='NETCDF4_CLASSIC') as ncf:
            filedata = ncf.variables[filevar.fvname["shum"]]
            if diaglev > 20 : print(filedata)
            time=ncf.variables[filevar.fvname["time"]][:]
            time_units=ncf.variables[filevar.fvname["time"]].units
            lat=ncf.variables[filevar.fvname["lat"]][:]
            lon=ncf.variables[filevar.fvname["lon"]][:]
            lev=ncf.variables[filevar.fvname["lev"]][:]
            idx=generate2Dindex(obs_cnt,LTTD,LNGD,datetime,time,time_units,lat,lon)
            data=get_data(element,obs_cnt,idx,option=option,lev=lev,filedata=filedata,temp=[],rhum=[],uwnd=[],vwnd=[])
#            for i in range(0,obs_cnt,1):
#                idata=[]*len(lev)
#                for j,lvl in enumerate(lev):
#                    if option in ["uwnd"] : idata[j]=uwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]]
#                    if option in ["vwnd"] : idata[j]=vwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]]
#                    if option in ["ddd"] : idata[j]=get_ddd(uwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]])
#                    if option in ["fff"] : idata[j]=get_fff(uwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]])
#                    if option in ["dpt_onshum"] : idata[j]=get_dpt_onshum(filedata[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]])
#                    if option in ["dpt_onrhum"] : idata[j]=get_dpt_onrhum(temp[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]],rhum[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]])
#                if i in [0]: 
#                    data=pandas.DataFrame([idata],index=[i+1],columns=convert_mb_to_pa(lev))
#                else:
#                    data=data.append(pandas.DataFrame([idata],index=[i+1],columns=convert_mb_to_pa(lev)))
#            data=pandas.DataFrame(index=range(1,obs_cnt+1,1),columns=convert_mb_to_pa(lev))
#            for i in range(0,obs_cnt,1):
#                for j,lvl in enumerate(lev):
#                    data.xs(i+1)[j]=get_dpt_onshum(filedata[idx.itime[i],[j],idx.ilat[i],idx.ilon[i]])
    else:
        option = "dpt_onrhum"
        with netCDF4.Dataset(filevar.file["rhum"], 'r',  format='NETCDF4_CLASSIC') as rhfile, \
            netCDF4.Dataset(filevar.file["tmp"], 'r',  format='NETCDF4_CLASSIC') as tfile:
            temp = tfile.variables[filevar.fvname["tmp"]]
            rhum=rhfile.variables[filevar.fvname["rhum"]]
            if diaglev > 20 : print(temp,rhum)
            time=rhfile.variables[filevar.fvname["time"]][:]
            time_units=rhfile.variables[filevar.fvname["time"]].units
            lat=rhfile.variables[filevar.fvname["lat"]][:]
            lon=rhfile.variables[filevar.fvname["lon"]][:]
            lev=rhfile.variables[filevar.fvname["lev"]][:]
            idx=generate2Dindex(obs_cnt,LTTD,LNGD,datetime,time,time_units,lat,lon)
            data=get_data(element,obs_cnt,idx,option=option,lev=lev,filedata=[],temp=temp,rhum=rhum,uwnd=[],vwnd=[])
#            for i in range(0,obs_cnt,1):
#                idata=[]*len(lev)
#                for j,lvl in enumerate(lev):
#                    if option in ["uwnd"] : idata[j]=uwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]]
#                    if option in ["vwnd"] : idata[j]=vwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]]
#                    if option in ["ddd"] : idata[j]=get_ddd(uwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]])
#                    if option in ["fff"] : idata[j]=get_fff(uwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]])
#                    if option in ["dpt_onshum"] : idata[j]=get_dpt_onshum(filedata[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]])
#                    if option in ["dpt_onrhum"] : idata[j]=get_dpt_onrhum(temp[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]],rhum[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]])
#                if i in [0]: 
#                    data=pandas.DataFrame([idata],index=[i+1],columns=convert_mb_to_pa(lev))
#                else:
#                    data=data.append(pandas.DataFrame([idata],index=[i+1],columns=convert_mb_to_pa(lev)))
#            data=pandas.DataFrame(index=range(1,obs_cnt+1,1),columns=convert_mb_to_pa(lev))
#            for i in range(0,obs_cnt,1):
#                for j,lvl in enumerate(lev):
#                    data.xs(i+1)[j]=get_dpt_onrhum(temp[idx.itime[i],[j],idx.ilat[i],idx.ilon[i]],rhum[idx.itime[i],[j],idx.ilat[i],idx.ilon[i]])
    obstore.obstore_write_data_element(outfile,nmlfile,indx,element,data)
    return(data)

def getaltlev(outfile,nmlfile,indx,filevar,datetime,obs_cnt,LTTD,LNGD,ALTTD):
    with netCDF4.Dataset(filevar.file["gph"], 'r',  format='NETCDF4_CLASSIC') as ncf:
        gph = ncf.variables[filevar.fvname["gph"]]
        time=ncf.variables[filevar.fvname["time"]][:]
        time_units=ncf.variables[filevar.fvname["time"]].units
        lat=ncf.variables[filevar.fvname["lat"]][:]
        lon=ncf.variables[filevar.fvname["lon"]][:]
        lev=ncf.variables[filevar.fvname["lev"]][:]
        indx=generate3Dindex(obs_cnt,LTTD,LNGD,ALTTD,datetime,time,time_units,lat,lon,lev,gph)
    return(lev[indx.ilev])

def obstore_copy_data_element(outfile,nmlfile,indx,element,infile):
    if element in obstore.obstore_read_batch_elements(infile,indx,nmlfile).Element.values:
        if diaglev > 0 : print(element)
        if diaglev > 1 : print(obstore.obstore_read_batch_elements(outfile,indx,nmlfile).query("Element == @element"))
        return(obstore.obstore_copy_data_element(outfile,nmlfile,indx,element,infile))    

def get_data(element,obs_cnt,idx,option="data",const=None,chnldic={},lev=numpy.empty(shape=[0]),lev_src_unit="mb",lev_req_unit="Pa",filedata=[],temp=[],rhum=[],uwnd=[],vwnd=[]):
    if lev.size == 0:
        for i in range(0,obs_cnt,1):
            if option in ["const"] : idata=const
            if option in ["data"] : idata=filedata[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]]
            if option in ["lev"] : idata=None
            if option in ["chnlfreq"] : idata=None
            if option in ["temp"] : idata=temp[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]]
            if option in ["uwnd"] : idata=uwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]]
            if option in ["vwnd"] : idata=vwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]]
            if option in ["ddd"] : idata=get_ddd(uwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]])
            if option in ["fff"] : idata=get_fff(uwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]])
            if option in ["hloswnd"] : idata=get_hloswnd(uwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]],idx.azmh[i+1])
            if option in ["dpt_onshum"] : idata=get_dpt_onshum(filedata[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]])
            if option in ["dpt_onrhum"] : idata=get_dpt_onrhum(temp[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]],rhum[idx.itime[i+1],idx.ilat[i+1],idx.ilon[i+1]])
            if i in [0]: 
                data=pandas.DataFrame([idata],index=[i+1],columns=[element])
            else:
                data=data.append(pandas.DataFrame([idata],index=[i+1],columns=[element]))
    elif "lev" in obslib.dfheader(idx):
        for i in range(0,obs_cnt,1):
            if option in ["const"] : idata=const
            if option in ["data"] : idata=filedata[idx.itime[i+1],idx.ilev[i+1],idx.ilat[i+1],idx.ilon[i+1]]
            if option in ["lev"] : idata=obslib.unit_convert(idx.lev[i+1],lev_src_unit,lev_req_unit)
            if option in ["chnlfreq"] : idata=ChanFreq(idx.lev[i+1],chnldic)
            if option in ["temp"] : idata=temp[idx.itime[i+1],idx.ilev[i+1],idx.ilat[i+1],idx.ilon[i+1]]
            if option in ["uwnd"] : idata=uwnd[idx.itime[i+1],idx.ilev[i+1],idx.ilat[i+1],idx.ilon[i+1]]
            if option in ["vwnd"] : idata=vwnd[idx.itime[i+1],idx.ilev[i+1],idx.ilat[i+1],idx.ilon[i+1]]
            if option in ["ddd"] : idata=get_ddd(uwnd[idx.itime[i+1],idx.ilev[i+1],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],idx.ilev[i+1],idx.ilat[i+1],idx.ilon[i+1]])
            if option in ["fff"] : idata=get_fff(uwnd[idx.itime[i+1],idx.ilev[i+1],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],idx.ilev[i+1],idx.ilat[i+1],idx.ilon[i+1]])
            if option in ["hloswnd"] : idata=get_hloswnd(uwnd[idx.itime[i+1],idx.ilev[i+1],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],idx.ilev[i+1],idx.ilat[i+1],idx.ilon[i+1]],idx.azmh[i+1])
            if option in ["dpt_onshum"] : idata=get_dpt_onshum(filedata[idx.itime[i+1],idx.ilev[i+1],idx.ilat[i+1],idx.ilon[i+1]])
            if option in ["dpt_onrhum"] : idata=get_dpt_onrhum(temp[idx.itime[i+1],idx.ilev[i+1],idx.ilat[i+1],idx.ilon[i+1]],rhum[idx.itime[i+1],idx.ilev[i+1],idx.ilat[i+1],idx.ilon[i+1]])
            if i in [0]: 
                data=pandas.DataFrame([idata],index=[i+1],columns=[element])
            else:
                data=data.append(pandas.DataFrame([idata],index=[i+1],columns=[element]))
    else:
        for i in range(0,obs_cnt,1):
                idata=[None]*len(lev)
                for j,lvl in enumerate(lev):
                    if option in ["const"] : idata=const
                    if option in ["data"] : idata[j]=filedata[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]]
                    if option in ["lev"] : idata=obslib.unit_convert(lvl,lev_src_unit,lev_req_unit)
                    if option in ["chnlfreq"] : idata=const
                    if option in ["temp"] : idata=temp[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]]
                    if option in ["uwnd"] : idata[j]=uwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]]
                    if option in ["vwnd"] : idata[j]=vwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]]
                    if option in ["ddd"] : idata[j]=get_ddd(uwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]])
                    if option in ["fff"] : idata[j]=get_fff(uwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]])
                    if option in ["hloswnd"] : idata=get_hloswnd(uwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]],vwnd[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]],idx.azmh[i+1])
                    if option in ["dpt_onshum"] : idata[j]=get_dpt_onshum(filedata[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]])
                    if option in ["dpt_onrhum"] : idata[j]=get_dpt_onrhum(temp[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]],rhum[idx.itime[i+1],[j],idx.ilat[i+1],idx.ilon[i+1]])
                if i in [0]: 
                    data=pandas.DataFrame([idata],index=[i+1],columns=convert_mb_to_pa(lev))
                else:
                    data=data.append(pandas.DataFrame([idata],index=[i+1],columns=convert_mb_to_pa(lev)))
    return(data)

def prepdata(subdic,data,findex,idxnam,elenam,option="const",DT=None,lev=None,lev_src_unit="mb",lev_req_unit="Pa",constdic={},filedata=[],temp=[],rhum=[],uwnd=[],vwnd=[]):
    obs_index=subdic["obs_index"]
    elist=obstore.obstore_create_element_table(nmlfile,obs_index)
    elenamlist=elist.Element.values
    obscount=len(data.index)

#    
#    filedim=getfiledim(DT)

    #gph=ncepradic.getdata(Year,"gph",element="gph") 
#    for
#    data1=data[data.index == index]
#    if obscount in [1]:
#        findex=fileindex(data,obscount,filedim,gph)
#    else:
#        findex=findex.append(fileindex(data,obscount,filedim,gph))
#    data=data.join(findex.lev,on='ProfileNo')
    if option in ["index"]:
        data=data.join(pandas.DataFrame(data.index,index=data.index,columns=["Index"]))
    if elenam in elenamlist:
        print(elenam)
        if lev is None:
            if "lev" in data.columns.values:
                lev=data.lev
            elif "lev" in findex.columns.values:
                lev=findex.lev
            else:
                lev=numpy.empty(shape=[0])
        if elenam in ["t2","td2","ddd10","fff10"]:
            lev=numpy.empty(shape=[0])
                
        #lev=obslib.unit_convert(lev,lev_src_unit,lev_req_unit)
        #print(lev)

        if option in ["date"]:
            Year=DT.year
            Month=DT.month
            Day=DT.day
            Hour=DT.hour
            Minute=DT.minute
            constdic={"Year" : Year, "Month" : Month, "Day" : Day, "Hour" : Hour, "Minute" : Minute,
                      "ReceiptYear" : Year, "ReceiptMonth" : Month, "ReceiptDay" : Day, "ReceiptHour" : Hour, "ReceiptMinute" : Minute}
            option="const"
        
        if option in ["const","chnlfreq"]:
            const=constdic[elenam]
        else:
            const=None
        if option in ["chnlfreq"]:
            chnldic=constdic["Chnldic"]
        else:
            chnldic={}
        if option in data.columns:
            data1=data[[option]].rename({option:elenam}, axis='columns')
            data=data.join(data1 ,on='Index')
        else:
            element=get_data(elenam,obscount,findex,option=option,const=const,chnldic=chnldic,lev=lev,lev_src_unit=lev_src_unit,lev_req_unit=lev_req_unit,filedata=filedata,temp=temp,rhum=rhum,uwnd=uwnd,vwnd=vwnd)
            data=data.join(element,on=idxnam)
    return(data)

def symulate_subtype(obs_info,obstypedic=None,maxindx=MAXINDX,temp=[],rhum=[],uwnd=[],vwnd=[]):
    DT = obs_info["timeinfo"]
    data = obs_info["data"]
    elist = obs_info["elist"]
    subtype = obs_info["subtype"]
    obstype = obs_info["obstype"]
    if obstypedic is None: obstypedic=obsdic.obstype[obstype]
    if subtype in obstypedic["subtype"]:
        indxnam="Index"
        altnam="Altitude"
        print("simulating subtype "+str(subtype))
        if diaglev > 0 : print(elist)
        elenams=elist.Element.values
        if diaglev > 5 : print(elenams)
        #########################################
        Year=DT.year
        Month=DT.month
        Day=DT.day
        Hour=DT.hour
        Minute=DT.minute
        DT=obslib.pydatetime(Year,Month,Day,Hour,Minute,0)
	print(DT)
        Year,filedim=nature.getfiledim(DT)
        #lev=nature.getdata(Year,"lev")
        gph=nature.getdata("gph",year=Year,month=Month,day=Day)
        temp=nature.getdata("tmp",year=Year,month=Month,day=Day)
        rhum=nature.getdata("rhum",year=Year,month=Month,day=Day)  
        u10m=nature.getdata("u10m",year=Year,month=Month,day=Day)
        v10m=nature.getdata("v10m",year=Year,month=Month,day=Day)
        slp=nature.getdata("slp",year=Year,month=Month,day=Day)
        psfc=nature.getdata("psfc",year=Year,month=Month,day=Day)
        t2m=nature.getdata("t2m",year=Year,month=Month,day=Day)
        rh2m=nature.getdata("rh2m",year=Year,month=Month,day=Day)  
        sh2m=nature.getdata("sh2m",year=Year,month=Month,day=Day) 
        #uwnd=nature.getdata("uwnd",year=Year,month=Month,day=Day)
        vwnd=nature.getdata("vwnd",year=Year,month=Month,day=Day)
    
#    if subtype is 22501 : 
#        altnam="HeightCOG"
#        indxnam="ProfileNo"
#    else:
        
    
        obsgroup=obstypedic["obsgroup"]
#        print(obsgroup)
#        if obsgroup in [1, 2, 3, 5]:
	print(len(data.index))
	print(max(data.index))
	print(data.columns)
        for indx in range(1,len(data.index)+1,1):
            DT=obslib.getdatetime(DT,data,indx)
            if Year != DT.year:
                Year,filedim=getfiledim(DT)
                gph=nature.getdata("gph",year=Year,month=Month,day=Day)
            if indx in [1]:
                findex=fileindex(data,indx,filedim,gph,altnam,DT)
            else:
                findex=findex.append(fileindex(data,indx,filedim,gph,altnam,DT))
#        if "Index" in data.columns.values: data=data.join(pandas.DataFrame(data.index,index=data.index),columns=["Index"])
#        print(data.columns.values, findex.columns.values)
#        if obsgroup not in [1,3] :
#            if "lev" in findex.columns.values: data=data.join(findex.lev,on='Index')
        subdic={"obs_index" : obstypedic[str(subtype)]}      #elist.obs_index.values with repeatition
        #data=LEOGEOAMV(data,findex,DT,nmlfile,maxindx=maxindx)
        if obsgroup not in [2]:
            data=prepdata(subdic,data,findex,indxnam,"Index",option="index")
            data=prepdata(subdic,data,findex,indxnam,"Year",option="date",DT=DT)
            data=prepdata(subdic,data,findex,indxnam,"Month",option="date",DT=DT)
            data=prepdata(subdic,data,findex,indxnam,"Day",option="date",DT=DT)
            data=prepdata(subdic,data,findex,indxnam,"Hour",option="date",DT=DT)
            data=prepdata(subdic,data,findex,indxnam,"Minute",option="date",DT=DT)
            
        data=prepdata(subdic,data,findex,indxnam,"ReceiptYear",option="date",DT=DT)
        data=prepdata(subdic,data,findex,indxnam,"ReceiptMonth",option="date",DT=DT)
        data=prepdata(subdic,data,findex,indxnam,"ReceiptDay",option="date",DT=DT)
        data=prepdata(subdic,data,findex,indxnam,"ReceiptHour",option="date",DT=DT)
        data=prepdata(subdic,data,findex,indxnam,"ReceiptMinute",option="date",DT=DT)
        data=prepdata(subdic,data,findex,indxnam,"SatID",option="const",constdic=obstypedic)
        data=prepdata(subdic,data,findex,indxnam,"ProfileNo",option="Index",constdic=obstypedic)
        data=prepdata(subdic,data,findex,indxnam,"OrigCtr",option="const",constdic=obstypedic)
        data=prepdata(subdic,data,findex,indxnam,"CloudMotionMethod",option="const",constdic=obstypedic)
        data=prepdata(subdic,data,findex,indxnam,"HeightAssMethod",option="const",constdic=obstypedic)
        data=prepdata(subdic,data,findex,indxnam,"HeightCOG",option="Altitude",constdic=obstypedic)
        data=prepdata(subdic,data,findex,indxnam,"AzimuthCOG",option="Azimuth",constdic=obstypedic)
        data=prepdata(subdic,data,findex,indxnam,"AltHeightAss",option="const",constdic=obstypedic)
        data=prepdata(subdic,data,findex,indxnam,"SatZenithAngle",option="const",constdic=obstypedic)
        data=prepdata(subdic,data,findex,indxnam,"AltPressLevel",option="lev")
        data=prepdata(subdic,data,findex,indxnam,"CmptVecWindSpeed",option="fff",uwnd=uwnd,vwnd=vwnd)
        data=prepdata(subdic,data,findex,indxnam,"CmptVecWindDirn",option="ddd",uwnd=uwnd,vwnd=vwnd)
        data=prepdata(subdic,data,findex,indxnam,"HLOSWIND",option="hloswnd",uwnd=uwnd,vwnd=vwnd)
        data=prepdata(subdic,data,findex,indxnam,"ChannelNumber","const",constdic=obstypedic)
        data=prepdata(subdic,data,findex,indxnam,"LidarClass","const",constdic=obstypedic)
        data=prepdata(subdic,data,findex,indxnam,"WindErr","const",constdic=obstypedic)
        data=prepdata(subdic,data,findex,indxnam,"AladinConfFlag","const",constdic=obstypedic)
        data=prepdata(subdic,data,findex,indxnam,"Temp",option="temp",temp=temp)
        data=prepdata(subdic,data,findex,indxnam,"AltAirTempLevel",option="temp",temp=temp)
        data=prepdata(subdic,data,findex,indxnam,"ShipDirec",option="const",constdic={"ShipDirec":0})
        data=prepdata(subdic,data,findex,indxnam,"ShipSpeed",option="const",constdic={"ShipSpeed":0})
        data=prepdata(subdic,data,findex,indxnam,"BuoyDirection",option="const",constdic={"BuoyDirection":0})
        data=prepdata(subdic,data,findex,indxnam,"BuoySpeed",option="const",constdic={"BuoySpeed":0})
        data=prepdata(subdic,data,findex,indxnam,"Pmsl",option="data",filedata=slp,lev=numpy.empty(shape=[0]))
        #data=prepdata(subdic,data,findex,indxnam,"Pmsl",option="data",filedata=slp,lev=numpy.empty(shape=[0]))
        data=prepdata(subdic,data,findex,indxnam,"NesdisQI1",option="const",constdic={"NesdisQI1":95})
        data=prepdata(subdic,data,findex,indxnam,"NesdisQI2",option="const",constdic={"NesdisQI2":95})
        data=prepdata(subdic,data,findex,indxnam,"NesdisQI3",option="const",constdic={"NesdisQI3":95})
       
        if subtype in range(10100,11100,100):
            data=prepdata(subdic,data,findex,indxnam,"t2",option="temp",temp=t2m,lev=numpy.empty(shape=[0]))
            data=prepdata(subdic,data,findex,indxnam,"td2",option="dpt_onshum",temp=t2m,filedata=sh2m,lev=numpy.empty(shape=[0]))
            data=prepdata(subdic,data,findex,indxnam,"ddd10",option="ddd",uwnd=u10m,vwnd=v10m,lev=numpy.empty(shape=[0]))
            data=prepdata(subdic,data,findex,indxnam,"fff10",option="fff",uwnd=u10m,vwnd=v10m,lev=numpy.empty(shape=[0]))
            data=prepdata(subdic,data,findex,indxnam,"DewPoint",option="dpt_onshum",temp=t2m,filedata=sh2m,lev=numpy.empty(shape=[0]))
        elif subtype in range(50100,50500,100):
            data=prepdata(subdic,data,findex,indxnam,"PlevelsA",option="lev",lev=lev)
            data=prepdata(subdic,data,findex,indxnam,"fff",option="fff",uwnd=uwnd,vwnd=vwnd,lev=lev)
            data=prepdata(subdic,data,findex,indxnam,"ddd",option="ddd",uwnd=uwnd,vwnd=vwnd,lev=lev)
            data=prepdata(subdic,data,findex,indxnam,"DewPoint",option="dpt_onrhum",temp=temp,rhum=rhum,lev=lev)
        else:
            data=prepdata(subdic,data,findex,indxnam,"PlevelsA",option="lev")
            data=prepdata(subdic,data,findex,indxnam,"fff",option="fff",uwnd=uwnd,vwnd=vwnd)
            data=prepdata(subdic,data,findex,indxnam,"ddd",option="ddd",uwnd=uwnd,vwnd=vwnd)
            data=prepdata(subdic,data,findex,indxnam,"DewPoint",option="dpt_onrhum",temp=temp,rhum=rhum)
            
    print(data.columns.values)
    return(data)
    
#    #########################################
#    (indx,pos_data,obs_cnt,tcols,data_len,data_end)=obstore.obstore_read_batchinfo(infile,indx)
#    if "ShipDirec" in elenams: obstore.obstore_write_data_element(outfile,nmlfile,indx,"ShipDirec",conform_dims(0,'cell',[obs_cnt,1]))
#    if "ShipSpeed" in elenams: obstore.obstore_write_data_element(outfile,nmlfile,indx,"ShipSpeed",conform_dims(0,'cell',[obs_cnt,1]))
#    if "BuoyDirection" in elenams: obstore.obstore_write_data_element(outfile,nmlfile,indx,"BuoyDirection",conform_dims(0,'cell',[obs_cnt,1]))
#    if "BuoySpeed" in elenams: obstore.obstore_write_data_element(outfile,nmlfile,indx,"BuoySpeed",conform_dims(0,'cell',[obs_cnt,1]))
#    #######################################
#    if "t2" in elenams: symulate_data(outfile,nmlfile,indx,"t2","tmp",filevar,DT,obs_cnt,Latitude,Longitude)
#    if "ddd10" in elenams: symulate_wind_2dfield(outfile,nmlfile,indx,"ddd10",["uwnd","vwnd"],filevar,DT,obs_cnt,Latitude,Longitude,option="ddd")
#    if "ddd10AmbWind" in elenams: symulate_wind_2dfield(outfile,nmlfile,indx,"ddd10AmbWind",["uwnd","vwnd"],filevar,DT,obs_cnt,Latitude,Longitude,option="ddd")
#    if "fff10" in elenams: symulate_wind_2dfield(outfile,nmlfile,indx,"fff10",["uwnd","vwnd"],filevar,DT,obs_cnt,Latitude,Longitude,option="fff")
#    if "WIND_SPED" in elenams: symulate_wind_2dfield(outfile,nmlfile,indx,"WIND_SPED",["uwnd","vwnd"],filevar,DT,obs_cnt,Latitude,Longitude,option="fff")
#    if "td2" in elenams: symulate_dpt_2dfield(outfile,nmlfile,indx,"td2",["shum"],filevar,DT,obs_cnt,Latitude,Longitude)
#    if "Pmsl" in elenams: symulate_data(outfile,nmlfile,indx,"Pmsl","slp",filevar,DT,obs_cnt,Latitude,Longitude)
#    ###############################
#
#            
#    if subtype in range(30100,30500,100): 
##        with netCDF4.Dataset(filevar.file["gph"], 'r',  format='NETCDF4_CLASSIC') as hfile:
##            filelev=hfile.variables[filevar.fvname["lev"]][:]
#        if "PlevelsA" in elenams: obstore_copy_data_element(outfile,nmlfile,indx,"PlevelsA",infile)
#        if "LEVL_HGHT" in elenams: obstore_copy_data_element(outfile,nmlfile,indx,"LEVL_HGHT",infile)
#        if "Altitude" in elenams: Altitude = obstore_copy_data_element(outfile,nmlfile,indx,"Altitude",infile)
##           print(obstore.getelenams(infile,nmlfile,subtype,indx))
# #          print(obstore.getelenams(outfile,nmlfile,subtype,indx))
#         #print(getaltlev(outfile,nmlfile,indx,filevar,DT,obs_cnt,Latitude,Longitude,Altitude))
####'TailNumber' 'CallSign' 'AcarsStatus' 'Year' 'Month' 'Day' 'Hour' 'Minute' 'Latitude' 'Longitude' 'FlightPhase' 'Altitude' 'PlevelsA' 'ddd' 'fff' 'TurbulenceDegree' 'Temp' 'ICE_DEGR' 'CharData'
#        if "Temp" in elenams: symulate_data_flightlev(outfile,nmlfile,indx,"Temp","tmp",filevar,DT,obs_cnt,Latitude,Longitude,Altitude)
#        if "ddd" in elenams: symulate_wind_flightlev(outfile,nmlfile,indx,"ddd",["uwnd","vwnd"],filevar,DT,obs_cnt,Latitude,Longitude,Altitude,option="ddd")
#        if "fff" in elenams: symulate_wind_flightlev(outfile,nmlfile,indx,"fff",["uwnd","vwnd"],filevar,DT,obs_cnt,Latitude,Longitude,Altitude,option="fff")
#        #if "DewPoint" in elenams: symulate_dpt_flightlev(outfile,nmlfile,indx,"DewPoint",["rhum","tmp"],filevar,DT,obs_cnt,Latitude,Longitude,Altitude)
#        if "Ucomponent" in elenams: symulate_data_flightlev(outfile,nmlfile,indx,"Ucomponent","uwnd",filevar,DT,obs_cnt,Latitude,Longitude,Altitude,option="uwnd")
#        if "Vcomponent" in elenams: symulate_data_flightlev(outfile,nmlfile,indx,"Vcomponent","vwnd",filevar,DT,obs_cnt,Latitude,Longitude,Altitude,option="vwnd")
#    elif subtype in range(50100,50500,100): 
#        with netCDF4.Dataset(filevar.file["gph"], 'r',  format='NETCDF4_CLASSIC') as hfile:
#            lev=hfile.variables[filevar.fvname["lev"]][:]
#            obstore.obstore_write_data_element(outfile,nmlfile,indx,"LEVL_RPLTN_CONT",extend_dims(len(lev),'index',obs_cnt))
#            obstore.obstore_write_data_element(outfile,nmlfile,indx,"PlevelsA",extend_dims(convert_mb_to_pa(lev),'index',obs_cnt))
#            obstore.obstore_write_data_element(outfile,nmlfile,indx,"LEVL_IDNY",conform_dims(32,'cell',[obs_cnt,17]))
#        if "LEVL_HGHT" in elenams: symulate_data(outfile,nmlfile,indx,"LEVL_HGHT","gph",filevar,DT,obs_cnt,Latitude,Longitude)
#        if "Temp" in elenams: symulate_data(outfile,nmlfile,indx,"Temp","tmp",filevar,DT,obs_cnt,Latitude,Longitude)
#        if "ddd" in elenams: symulate_wind_3dfield(outfile,nmlfile,indx,"ddd",["uwnd","vwnd"],filevar,DT,obs_cnt,Latitude,Longitude,option="ddd")
#        if "fff" in elenams: symulate_wind_3dfield(outfile,nmlfile,indx,"fff",["uwnd","vwnd"],filevar,DT,obs_cnt,Latitude,Longitude,option="fff")
#        #if "DewPoint" in elenams: symulate_dpt_3dfield(outfile,nmlfile,indx,"DewPoint",["rhum","tmp"],filevar,DT,obs_cnt,Latitude,Longitude)
#        if "Ucomponent" in elenams: symulate_data(outfile,nmlfile,indx,"Ucomponent","uwnd",filevar,DT,obs_cnt,Latitude,Longitude,option="uwnd")
#        if "Vcomponent" in elenams: symulate_data(outfile,nmlfile,indx,"Vcomponent","vwnd",filevar,DT,obs_cnt,Latitude,Longitude,option="vwnd")
#    else:
#        if "PlevelsA" in elenams: symulate_data(outfile,nmlfile,indx,"PlevelsA","pres",filevar,DT,obs_cnt,Latitude,Longitude)
#        if "LEVL_HGHT" in elenams: symulate_data(outfile,nmlfile,indx,"LEVL_HGHT","gph",filevar,DT,obs_cnt,Latitude,Longitude)
#        if "Temp" in elenams: symulate_data(outfile,nmlfile,indx,"Temp","tmp",filevar,DT,obs_cnt,Latitude,Longitude)
#        if "ddd" in elenams: symulate_wind_3dfield(outfile,nmlfile,indx,"ddd",["uwnd","vwnd"],filevar,DT,obs_cnt,Latitude,Longitude,option="ddd")
#        if "fff" in elenams: symulate_wind_3dfield(outfile,nmlfile,indx,"fff",["uwnd","vwnd"],filevar,DT,obs_cnt,Latitude,Longitude,option="fff")
#        if "DewPoint" in elenams: symulate_dpt_3dfield(outfile,nmlfile,indx,"DewPoint",["rhum","tmp"],filevar,DT,obs_cnt,Latitude,Longitude)
#        if "Ucomponent" in elenams: symulate_data(outfile,nmlfile,indx,"Ucomponent","uwnd",filevar,DT,obs_cnt,Latitude,Longitude,option="uwnd")
#        if "Vcomponent" in elenams: symulate_data(outfile,nmlfile,indx,"Vcomponent","vwnd",filevar,DT,obs_cnt,Latitude,Longitude,option="vwnd")
def subtype_filter(subtypegroup, subtypelist):
	indxlist=[]
	if subtypelist is None:
		indxlist = range(0,len(subtypegroup),1)
	else:
	    for indx,subtype in enumerate(subtypegroup):
		if subtype in subtypelist: 
			indxlist=indxlist+[indx]
	return(indxlist)
    

def symulate_obstore(outpath,inpath,nmlpath,DT,obstype,filevar=None,maxindx=MAXINDX,subtypelist=None):
    obstypedic=obsdic.obstype[obstype]
    filename=obstypedic["filename"]
    output_file="%s/%s" % (outpath,filename)
    input_file="%s/%s" % (inpath,filename)
    obs_index_nml="obs_index_nml"
    nmlfile="%s/%s" % (nmlpath,obs_index_nml)
    obs_subtype_nml="obs_subtype.nml"
    subtype_nmlfile="%s/%s" % (nmlpath,obs_subtype_nml)
    obsgroup=obstypedic["obsgroup"]
    if obstype in ["aladin","leogeo","hlosw"] : 
        batchcount=3
        subtype=obstypedic["subtype"][0]
        subtypegroup=[subtype]*batchcount
        obs_index=obstypedic[str(subtype)]
        elistgroup=[obstore.obstore_create_element_table(nmlfile,obs_index)]*batchcount
        print(elistgroup[0])
        Twindow=obslib.timeperiod(year=0,mn=0,dy=0,hh=6,mm=0,ss=0)
        Tstart=DT-Twindow/2
        Tdelta=Twindow/batchcount
        datagroup=[None]*batchcount
        for indx in range(1,batchcount+1,1):
            Tstop=Tstart+Tdelta
            location=aladin_profile(DT,Tstart,Tstop,outpath,nmlpath,maxindx=maxindx)
            datagroup[indx-1]=location
            #datagroup[indx-1]=LEOGEOAMV(location,findex,Tstart,output_file,nmlfile,obscount,maxindx=maxindx)
            #datagroup[indx-1]=symulate_subtype(location,elistgroup[indx],obstypedic,maxindx=maxindx)
            print(DT,Tstart,Tstop)
            Tstart=Tstop
    else:
        exclude_subtype_list=[]		#[11100,40100,50400]
        with open(input_file, "rb") as infile:
            obsgroup=obstore.obstore_read_header(infile,287,1)
            print(obsgroup)
#        head_305 = obstore.obstore_read_header(infile,1,305)
#        head_306_339 = obstore.obstore_read_real(infile,306,34)
#        #obstore.obstore_write_formatted(itertools.chain(head_305,head_306_339),">305q34d",outfile)
#        batchinfo=obstore.obstore_headerinfo(infile,nmlfile)
#        #obstore.obstore_write_header(outfile,nmlfile,batchinfo)
#        #obstore.obstore_erase_data(outfile)
#        ###### Erase function allocate the space for the Data to be written #########
#        #obstore.obstore_write_datetime(DT,outfile,nmlfile)
#        ##### Keep the header for all batches before writing the date time information
#        for indx in batchinfo.Batch_Index.values:
#            print(indx)
#            if indx not in [0]: print(obstore.obstore_read_subhead_segment(infile,"lut",indx,1,LUTSIZE))
#        if diaglev > 5 : print(obstore.obstore_read_subtype(infile))
            subtypegroup=obstore.obstore_read_subtype(infile)
	    indxlist=subtype_filter(subtypegroup,subtypelist)
	    print(indxlist)
	    subtypegroup=subtypegroup[indxlist[:]]
	    print(subtypegroup)
            batchcount=len(subtypegroup)
            elistgroup=[None]*batchcount
            datagroup=[None]*batchcount
            for indx,preindx in enumerate(indxlist,start=1):
		subtype=subtypegroup[indx-1]
		subtype_name=obslib.get_subtype_name(subtype_nmlfile,subtype)
                print(indx,subtype_name)
                elist=obstore.obstore_read_batch_elements(infile,preindx+1,nmlfile,maxindx=maxindx)
                elistgroup[indx-1]=elist
                print(elist.obs_index.values)
                #print(elist)
                count=0
                for element in obsdic.stationelist:
                    elemdata = obstore.obstore_read_data_element(infile,nmlfile,preindx+1,element,maxindx=maxindx)
                    if elemdata is not None:
                        count+=1
                        if count == 1:
                            location=elemdata
                        else:
                            location=location.join(elemdata)
#                data=prepdata(subdic,data,findex=None,"ProfileNo","Year",option="date",DT=DT)
#                data=prepdata(subdic,data,findex=None,"ProfileNo","Month",option="date",DT=DT)
#                data=prepdata(subdic,data,findex=None,"ProfileNo","Day",option="date",DT=DT)
#                data=prepdata(subdic,data,findex=None,"ProfileNo","Hour",option="date",DT=DT)
#                data=prepdata(subdic,data,findex=None,"ProfileNo","Minute",option="date",DT=DT)
                datagroup[indx-1]=location
                
#            if subtype not in exclude_subtype_list:
#                print(indx, subtype)
#                if diaglev > 5 : print(obstore.obstore_read_batch_elements(infile,indx,nmlfile))
#                symulate_subtype(outfile,DT,nmlfile,infile,filevar,batchinfo,int(subtype),indx)         
#                if diaglev > 5 : print(obstore.obstore_read_batch_elements(outfile,indx,nmlfile))
#            else:
#                #obstore.obstore_copy_batchinfo(outfile,nmlfile,indx,int(subtype),infile)
#                print("Skipping observation subtype %s"%(subtype))
#                elist=obstore.obstore_read_batch_elements(infile,indx,nmlfile)
#                elistgroup=elistgroup.append(elist)
#                #obstore.obstore_write_subtype_index(outfile,nmlfile,indx,subtype,0,elist)
#                obstore.erase_data_batch(outfile,nmlfile,subtype,indx)
    #print(subtypegroup)
    #print(elistgroup)
    for indx,subtype in enumerate(subtypegroup[0:batchcount],start=1):
	obs_info={
		"timeinfo" : DT,
		"data" : datagroup[indx-1],
		"elist" : elistgroup[indx-1],
		"subtype" : subtypegroup[indx-1],
		"obstype" : obstype,
		}
        datagroup[indx-1]=symulate_subtype(obs_info,obstypedic,maxindx=maxindx)
        print(datagroup[indx-1])
    
    print("Writting "+str(batchcount)+" batches of data to "+ output_file)
    with open(output_file, "wb+") as outfile:
        (datapos,datalen,dataend)=obstore.create_obstore(DT,outfile,nmlfile,obsgroup,subtypegroup,elistgroup,datagroup,batchcount=batchcount,header_offset=339,maxindx=maxindx,lut_ncols=LUTSIZE)
        print(datapos,datalen,dataend)
    print("Writting to "+output_file+ " is completed")
    obsmod.obs_frame(datagroup,subtypegroup,outpath,filename=obstype,option=1)
    return(datagroup)

def sose_merge_data(infodic):
    obstore_info={}
    obstore_info["inpath"] =  infodic["inpath"]
    obstore_info["outpath"] =  infodic["outpath"]
    obstore_info["nmlpath"] = infodic["nmlpath"]
    obstore_info["obstype"] = infodic["obstype"]
    obstore_info["timeinfo"] = infodic["timeinfo"]
    obstore_info["header_offset"] = infodic["header_offset"]
    obstore_info["lut_ncols"] = infodic["lut_ncols"]
    obstore_info["maxindx"] = infodic["maxindx"]
    obstore_info["subtypelist"] = infodic["subtypelist"]
    datafilter={}
    datafilter["synbuoyloc"] = infodic["synbuoyloc"]
    datafilter["array_weight"] = infodic["array_weight"]
    datafilter["latmin"] = infodic["latmin"]
    datafilter["latmax"] = infodic["latmax"]
    datafilter["lonmin"] = infodic["lonmin"]
    datafilter["lonmax"] = infodic["lonmax"] 
    datafilter["subtypelist"] = infodic["subtypelist"]
    obstore_info=filtered_obstore_read(obstore_info,datafilter)
    subtypegroup=obstore_info["subtypegroup"]
    batchcount=len(subtypegroup)
    for indx,subtype in enumerate(subtypegroup[0:batchcount],start=1):
	if subtype in infodic["syntype"] :
    		obs_info={}
    		obs_info["obstype"] = infodic["obstype"]
    		obs_info["syntype"] = infodic["syntype"]
    		obs_info["timeinfo"] = infodic["timeinfo"]
		obs_info["subtype"] = subtype
		fltrdata=obstore_info["datagroup"][indx-1]
    		obs_info["elist"] = obstore_info["elistgroup"][indx-1]
		fltrdata=fltrdata.append(synbuoy(datafilter,obs_info),ignore_index=True)
    		fltrdata=indexreset(fltrdata)
		obstore_info["datagroup"][indx-1]=fltrdata
	else:
		fltrdata=obstore_info["datagroup"][indx-1]
		print(subtype)
    #obstore_info["datagroup"] = datagroup
    obstore_info=obsmod.create_obstore(obstore_info)
    return(obstore_info)

def filtered_obstore_read(obstore_info,datafilter):
    inpath = obstore_info["inpath"]
    outpath = obstore_info["outpath"]
    nmlpath = obstore_info["nmlpath"]
    obstype = obstore_info["obstype"]
    DT = obstore_info["timeinfo"]
    header_offset = obstore_info["header_offset"]
    lut_ncols = obstore_info["lut_ncols"]
    maxindx = obstore_info["maxindx"]
    obstypedic=obsdic.obstype[obstype]
    filename=obstypedic["filename"]
    output_file="%s/%s" % (outpath,filename)
    input_file="%s/%s" % (inpath,filename)
    obs_index_nml="obs_index_nml"
    nmlfile="%s/%s" % (nmlpath,obs_index_nml)
    obs_subtype_nml="obs_subtype.nml"
    subtype_nmlfile="%s/%s" % (nmlpath,obs_subtype_nml)
    obsgroup=obstypedic["obsgroup"]
    if obstype in ["aladin","leogeo","hlosw"] : 
        batchcount=3
        subtype=obstypedic["subtype"][0]
        subtypegroup=[subtype]*batchcount
        obs_index=obstypedic[str(subtype)]
        elistgroup=[obstore.obstore_create_element_table(nmlfile,obs_index)]*batchcount
        print(elistgroup[0])
        Twindow=obslib.timeperiod(year=0,mn=0,dy=0,hh=6,mm=0,ss=0)
        Tstart=DT-Twindow/2
        Tdelta=Twindow/batchcount
        datagroup=[None]*batchcount
        for indx in range(1,batchcount+1,1):
            Tstop=Tstart+Tdelta
            location=aladin_profile(DT,Tstart,Tstop,outpath,nmlpath,maxindx=maxindx)
            datagroup[indx-1]=location
            print(DT,Tstart,Tstop)
            Tstart=Tstop
    else:
        exclude_subtype_list=[]		#[11100,40100,50400]
        with open(input_file, "rb") as infile:
            obsgroup=obstore.obstore_read_header(infile,287,1)
            print(obsgroup)
            subtypegroup=obstore.obstore_read_subtype(infile)
            subtypelist = datafilter["subtypelist"]
	    indxlist=subtype_filter(subtypegroup,subtypelist)
	    print(indxlist)
	    subtypegroup=subtypegroup[indxlist[:]]
	    print(subtypegroup)
            batchcount=len(subtypegroup)
            elistgroup=[None]*batchcount
            datagroup=[None]*batchcount
            for indx,preindx in enumerate(indxlist,start=1):
		subtype=subtypegroup[indx-1]
		subtype_name=obslib.get_subtype_name(subtype_nmlfile,subtype)
                print(indx,subtype_name)
                elist=obstore.obstore_read_batch_elements(infile,preindx+1,nmlfile,maxindx=maxindx)
                elistgroup[indx-1]=elist
                print(elist.obs_index.values)
                #print(elist)
                count=0
                for element in obsdic.station_call_list:
                    elemdata = obstore.obstore_read_data_element(infile,nmlfile,preindx+1,element,maxindx=maxindx)
                    if elemdata is not None:
                        count+=1
                        if count == 1:
                            location=elemdata
                        else:
                            location=location.join(elemdata)
                datagroup[indx-1]=location
    #### Synthetic Obstore Data Selection for all Batches
    for indx,subtype in enumerate(subtypegroup[0:batchcount],start=1):
        #print(datagroup[indx-1])
	obs_info={
		"timeinfo" : DT,
		"data" : datagroup[indx-1],
		"elist" : elistgroup[indx-1],
		"subtype" : subtypegroup[indx-1],
		"obstype" : obstype,
		}
        datagroup[indx-1]=symulate_subtype(obs_info,obstypedic,maxindx=maxindx)
        #print(datagroup[indx-1])
    obstore_info={
	"obstype" : obstype,
	"obsgroup" : obsgroup,
	"subtypegroup" : subtypegroup,
	"elistgroup" : elistgroup,
	"datagroup" : datagroup,
	"timeinfo" : DT,
	"outpath" : outpath,
	"filename" : filename,
	"nmlpath" : nmlpath,
	"maxindx" : maxindx,
	"batchcount" : batchcount,
	"header_offset" : header_offset,
	"lut_ncols" : lut_ncols,
		}
    if datafilter is not None: obstore_info=filter_data(obstore_info,datafilter)
    datagroup = obstore_info["datagroup"]
    return(obstore_info)
    
def filter_data(obstore_info,datafilter):
    obstype = obstore_info["obstype"]
    obsgroup = obstore_info["obsgroup"]
    subtypegroup = obstore_info["subtypegroup"]
    elistgroup = obstore_info["elistgroup"]
    datagroup = obstore_info["datagroup"]
    DT = obstore_info["timeinfo"]
    outpath = obstore_info["outpath"]
    filename = obstore_info["filename"]
    nmlpath = obstore_info["nmlpath"]
    maxindx = obstore_info["maxindx"]
    batchcount = obstore_info["batchcount"]
    header_offset = obstore_info["header_offset"]
    lut_ncols = obstore_info["lut_ncols"]
    subtype_select = datafilter["subtypelist"]
    latmin = datafilter["latmin"]
    latmax = datafilter["latmax"]
    lonmin = datafilter["lonmin"]
    lonmax = datafilter["lonmax"]
    for indx,subtype in enumerate(subtypegroup[0:batchcount],start=1):
	if subtype in subtype_select :
		fltrdata=datagroup[indx-1]
		fltrdata=datamask(fltrdata,latmin,latmax,lonmin,lonmax)
		print(subtype)
		obs_info={
			"obstype" : obstype,
			"timeinfo" : DT,
			"subtype" : subtype,
			"elist"	: elistgroup[indx-1],
			}
		datagroup[indx-1]=fltrdata
    obstore_info["datagroup"] = datagroup
    return(obstore_info)

def datamask(data,latmin,latmax,lonmin,lonmax,maskvalue=-1.07374e+09):
	#print(data.index)
	size=len(data.Latitude)
	maskindx=data[(data.Latitude > numpy.repeat(latmin,size)) & (data.Latitude < numpy.repeat(latmax,size)) & (data.Longitude > numpy.repeat(lonmin,size)) & (data.Longitude < numpy.repeat(lonmax,size)) ].index
	data = data.drop(maskindx)
	data=indexreset(data)
	#print(data.index)
	return(data)

def synbuoy(datafilter,obs_info):
    synbuoyloc = datafilter["synbuoyloc"]
    arrywght = datafilter["array_weight"] 
    locinfo=pandas.read_table(synbuoyloc).rename(columns={"buoyid":"BuoyID","lat":"Latitude","lon":"Longitude",})
    locinfo=locinfo[locinfo.priority >= arrywght]
    data=locinfo[["BuoyID","Latitude","Longitude"]]    
    data=indexreset(data)
    obs_info["data"]=data	
    data=symulate_subtype(obs_info)
    #print(data)
    return(data)

def indexreset(data):
	data = data.reset_index(drop=True)	# Fresh index start with 0
	data = data.shift()[1:len(data)+1]	# To start index with 1
	return(data)
