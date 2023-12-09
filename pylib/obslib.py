#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 19:32:07 2019

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

#import obsdic
import pandas
import numpy
import struct
import datetime
import glob
from itertools import chain
NAN_VAL_INT=-32768
NAN_VAL=-1.07374182e+09

diaglev=int(os.environ.get('GEN_MODE',0))
def errprint(*args, **kwargs):
    if diaglev > 0: print(*args, file=sys.stderr, **kwargs)

def DataFrame(data=None, index=None, columns=None, dtype=None, copy=None):
	return(pandas.DataFrame(data, index, columns, dtype, copy))

def globlist(string):
	return(glob.glob(string)[::-1])

def dic_print(dicid):
    for (key,val) in dicid.items():
       print(key+" : "+str(val))

def printall(array1):
	with numpy.printoptions(threshold=numpy.inf):
		print(array1)

########Golden key for binary write 20221027########################
def binary_write(data,pos,obsfile,binwidth=8):
    obsfile.seek((pos-1)*binwidth,0)
    if type(data) is int or type(data) is numpy.int64 :
    	dlen=1
    	fmtstr=struct.Struct(">"+str(dlen)+"q")
        status=obsfile.write(fmtstr.pack(data))
    elif type(data) is float or type(data) is numpy.float64 :
    	dlen=1
    	fmtstr=struct.Struct(">"+str(dlen)+"d")
        status=obsfile.write(fmtstr.pack(data))
    else:
    	dlen=len(data)
    	if type(data[0]) is int or isinstance(data[0], numpy.int64) : 
		fmtstr=struct.Struct(">"+str(dlen)+"q")
    	elif type(data[0]) is float or isinstance(data[0], numpy.float64) : 
		fmtstr=struct.Struct(">"+str(dlen)+"d")
	else :
		print("Datatype "+str(type(data[0]))+" is not yet supported")
        status=obsfile.write(fmtstr.pack(*data))
    return(status)
####################################################################

#def binary_write_formatted(obsfile,pos,word_len,data_fmt,data):
#    obsfile.seek((pos-1)*word_len,0)
#    frmtstr=struct.Struct(data_fmt)
#    obsfile.write(frmtstr.pack(*data))

#def binary_write(obsfile,pos,size,fmtkey,data,binwidth=8):
#    obsfile.seek((pos-1)*binwidth,0)
#    data_fmt=">"+str(size)+fmtkey
#    frmtstr=struct.Struct(data_fmt)
#    if type(data) is int:
#        obsfile.write(frmtstr.pack(data))
#    else:
#        obsfile.write(frmtstr.pack(*data))

def binary_read(obsfile,pos,size=1,fmtkey="q",binwidth=8):
    obsfile.seek((pos-1)*binwidth,0)
    data = obsfile.read(size*binwidth)
    data_fmt=">"+str(size)+fmtkey
    val = struct.unpack(data_fmt, data)
    if size is 1:
        return(val)[0]
    else:
        return(numpy.asarray(val))

def binary_read_header(obsfile,pos,size=1,fmtkey="q",binwidth=8):
    data=binary_read(obsfile,pos,size,fmtkey,binwidth)
    return(data)

def binary_read_data(obsfile,pos,size=1,reclen=1,fmtkey="d",binwidth=8):
    data=binary_read(obsfile,pos,size,fmtkey,binwidth)
    numrec=size/reclen
    data=numpy.reshape(data,[numrec,reclen])
    return(data)

def binary_read_record(obsfile,pos,numrec=1,reclen=1,fmtkey="d",binwidth=8):
    size=numrec*reclen
    data=binary_read(obsfile,pos,size,fmtkey,binwidth)
    data=numpy.reshape(data,[numrec,reclen])
    return(data)

def get_datapos(obsfile,fmtkey="q",binwidth=8):
    (pos,rowlen,rcnt)=binpos(obsfile,"data",fmtkey,binwidth)
    datapos=binary_read(obsfile,pos)
    #print("Data position : "+str(datapos))
    return(datapos)

def binary_segment_locate(obsfile,p1,p2,p3=None,fmtkey="q",binwidth=8):
    pos = binary_read_header(obsfile,p1,1,fmtkey,binwidth)
    tcol = binary_read_header(obsfile,p2,1,fmtkey,binwidth)
    if p3 is not None :
       nrow = binary_read_header(obsfile,p3,1,fmtkey,binwidth)
    else :
       nrow = 1
    return(pos,tcol,nrow)

def binfmtkey(obsfile,sec_nam):
    return {
       'halp': "q",
       'hbet': "q",
       'hgam': "d",
       'ldc' : "d",
       'rdc' : "d",
       'cdc' : "d",
       'lut' : "q",
       'data': "d",
    }.get(sec_nam,"d") 

def binpos(obsfile,sec_nam,fmtkey="q",binwidth=8):
    return {
        'halp': (1,256,1),
        'hbet': binary_segment_locate(obsfile,100,101,None,fmtkey,binwidth),
        'hgam': binary_segment_locate(obsfile,105,106,None,fmtkey,binwidth), 
        'ldc' : binary_segment_locate(obsfile,110,111,112,fmtkey,binwidth),
        'rdc' : binary_segment_locate(obsfile,115,116,117,fmtkey,binwidth),
        'cdc' : binary_segment_locate(obsfile,120,121,122,fmtkey,binwidth),
        'lut' : binary_segment_locate(obsfile,150,151,152,fmtkey,binwidth),
        'data': binary_segment_locate(obsfile,160,161,None,fmtkey,binwidth)
    }.get(sec_nam, (0,0,0))
    
def binary_get_segment_endpos(obsfile,sec_nam="halp",fmtkey="q",binwidth=8):
    (pos,rowlen,rcnt)=binpos(obsfile,sec_nam,fmtkey,binwidth)
    if pos > 0 :
       seglen=rowlen*rcnt
       segendpos=pos+seglen-1
    else :
       segendpos=None
    return(segendpos)

def binary_file_segment_read(obsfile,sec_nam="halp",recno=None,fmtkey="q",binwidth=8):
    (pos,rowlen,rcnt) = binpos(obsfile,sec_nam)
    fmtkey = binfmtkey(obsfile,sec_nam)
    seglen=rowlen*rcnt
    if pos > 0 :
      # print(pos,seglen)
       segment=binary_read(obsfile,pos,seglen,fmtkey,binwidth)
       if seglen > rowlen : segment=numpy.reshape(segment,[rcnt,rowlen])
    else :
       pos = None
       segment=[]
    if recno is not None:
       if recno <= rcnt : segment= segment[recno-1]
    segendpos=binary_get_segment_endpos(obsfile,sec_nam)
    print("Read "+sec_nam+" "+str(pos)+":"+str(segendpos))
    #print(segment)
    return(segment)

def read_lut(obsfile,recno=None):
    lut=binary_file_segment_read(obsfile,"lut",recno)
    return(lut) 

def read_cdc(obsfile,recno=None):
    cdc=binary_file_segment_read(obsfile,"cdc",recno)
    return(cdc) 

def read_rdc(obsfile,recno=None):
    rdc=binary_file_segment_read(obsfile,"rdc",recno)
    return(rdc) 

def read_ldc(obsfile,recno=None):
    ldc=binary_file_segment_read(obsfile,"ldc",recno)
    return(ldc) 

def binary_read_data_header(obsfile):
    halp = binary_file_segment_read(obsfile,"halp")
    hbet = binary_file_segment_read(obsfile,"hbet")
    hgam = binary_file_segment_read(obsfile,"hgam",fmtkey="d")
    ldc = read_ldc(obsfile)
    rdc = read_rdc(obsfile)
    cdc = read_cdc(obsfile)
    lut = read_lut(obsfile)
    neckpoint = binary_get_segment_endpos(obsfile,"lut")
    datapoint=get_datapos(obsfile)
    necklen=datapoint-neckpoint
    necklace = binary_read_data(obsfile,neckpoint,necklen)
    return(halp,hbet,hgam,ldc,rdc,cdc,lut,necklace)

#####################################################################################################
### Syntax:
###        nmldata=obslib.frame_data([[subtype,obstype]],columns=["subtype","obstype"],data=nmldata) 
#####################################################################################################
def frame_data(values,columns=None,index=None,data=None) :
    if isinstance(data, pandas.DataFrame):
       data = data.append(pandas.DataFrame(values,index=index,columns=columns))
    else :
       data = pandas.DataFrame(values,index=index,columns=columns)
    return(data)

def getascii(string):
    return([ord(c) for c in string])

def getstring(ascii):
    ascii=list(ascii)
    for i,a in enumerate(ascii):
        if a in range(1,256,1):
            ascii[i]=int(a)
        else:
            ascii[i]=32
    string=''.join(chr(int(a)) for a in ascii)
    return(string)

def str_list_index(val,listin):
	indx = numpy.where(numpy.array(listin)==str(val))[0][0]
	return(indx)

def nan_sort(a):
    temp = a.copy()
    temp= numpy.ma.masked_where(temp ==-3.2768e+04, temp)
    return temp.sort()

def mask_array(a,missing=numpy.nan):
    temp = a.copy()
    if str(temp.dtype) == "int16": temp = numpy.ma.masked_where(temp==32767, temp)
    if str(temp.dtype) == "uint16": temp = numpy.ma.masked_where(temp==65535, temp)
    temp = numpy.ma.masked_where(temp==missing, temp)
    temp = temp.astype('float64')
    arry = numpy.ma.filled(temp, numpy.nan)
    #temp = numpy.ma.masked_object(temp,missing)
    return(arry)

def binsort(a,binmin=numpy.NINF,binmax=numpy.inf,missing=numpy.nan):
    temp = a.copy()
    temp= numpy.ma.masked_where(temp =="nan", temp)
    temp= numpy.ma.masked_where(temp =='nan', temp)
    temp= numpy.ma.masked_where(temp ==missing, temp)
    #temp= numpy.ma.masked_where(temp ==str(missing), temp)
    temp= numpy.ma.masked_where(temp < binmin, temp)
    temp= numpy.ma.masked_where(temp > binmax, temp)
    temp=temp[~temp.mask] 
    temp.sort()
    return(temp)

def unique_int(array,missing=numpy.nan):
    intlst=[]
    for elem in array:
        if elem != missing : 
           intlst=intlst+[int(elem)]
    int_frm=pandas.DataFrame(intlst,columns=["Data"])
    int_arr=int_frm.Data.unique()
    return(int_arr)

def str_zfill(intval,padlen):
    return(str(intval).zfill(padlen))
    
def compare(la,ra):
    for i in range(0,len(ra)):
        print(i,la[i],ra[i])
    return(arraydiff(la,ra))

def arraydiff(la,ra):
    #print(la,ra)
    difflist=[]
    for i in range(0,len(la),1):
        if la[i] != ra[i]:
            difflist.append(i)
    return(difflist)

def arraymatch(la,ra):
    #print(la,ra)
    matchlist=[]
    for i in range(0,len(la),1):
        if la[i] == ra[i]:
            matchlist.append(i)
    return(matchlist)

def mean_of(a,absol=False):
	if absol: a=numpy.absolute(a)
	mean=float("{:.4f}".format(numpy.nanmean(a)))
	return(mean)

def append(df,ary):
    df.loc[-1] = numpy.array(ary)  # adding a row
    df.index = df.index + 1  # shifting index
    df = df.sort_index()  # sorting by index
    return(df)
    
def unit_convert(values,src_unit=None,req_unit=None,unit_fctr=None):
    if src_unit is "mb" and req_unit is "Pa" : unit_fctr=100
    if unit_fctr is not None:
        if isinstance(values, (list,)):
            return([int(i)*unit_fctr for i in values])
        else:
            return(values*unit_fctr)

def frame_window_filter(data,item,minval,maxval):
	newdf=data.query(item+" > @minval and "+item+" < @maxval")
	return(newdf)

def frame_select_filter(data,item,value):
	newdf=data.query(item+" == @value")
	return(newdf)

###############################################################
#### Date-time related functions
###############################################################

def today(fmtstr="%Y%m%d"):
	DT=datetime.date.today()
	today_fmted=DT.strftime(fmtstr)
	return(today_fmted)

def mins_sinceTref(DT,Tref=None):
	if Tref is not None : 
		Tref = pydate(date=str(Tref))
	else :
		Tref = pydate(date=19700101)
	delta = DT - Tref
	delta_mins = to_minutes(delta)
	return(delta_mins)

def localtime(DT,lon):
    time=lon/0.25
    seconds=time*60
    delta=datetime.timedelta(seconds=seconds)
    if delta.days < 0:
        delta = datetime.timedelta(seconds=delta.total_seconds() + 3600*24)
    if delta.days > 0:
        delta = datetime.timedelta(seconds=delta.total_seconds() - 3600*24)
    LT=DT+delta
    return(LT)

def timezoneshift(DT,zone):
    return(localtime(DT,obsdic.timezone[zone]))

def midnightlon(DT):
    hour=DT.hour
    minute=DT.minute
    second=DT.second
    time=(hour*60+minute+second/60.0)
    lon=-0.25*time
    while lon < -180 :
        lon = lon + 360
    return(lon)

def dawnlon(DT):
    hour=DT.hour
    minute=DT.minute
    second=DT.second
    time=(hour*60+minute+second/60.0)
    lon=-0.25*time+90
    while lon < -180 :
        lon = lon + 360
    return(lon)
    
def timeperiod(year=0,mn=0,dy=0,hh=0,mm=0,ss=0):
    days=int(year)*365+int(mn)*30+int(dy)
    seconds=(int(hh)*3600+int(mm)*60+int(ss))
    return(datetime.timedelta(days=days,seconds=seconds))

def to_minutes(TimeDelta):
    mins=(TimeDelta.days*1440+TimeDelta.seconds/60.0)
    return(mins)
    
def to_seconds(TimeDelta):
    return(TimeDelta.days*1440*60.0+TimeDelta.seconds)

def profileid(Tnow,Tnode,Profsec=None,Tp=None):
    if Tp is None: Tp=timeperiod(ss=Profsec)
    pid=int(to_minutes(Tnow-Tnode)/to_minutes(Tp))
    return(pid)

def orbitalshift(Tnow,Tnode,Profsec=None,Torb=None,Rcyc=1):
    if Torb is None: Torb=timeperiod(ss=Profsec)
    orbid=int(to_minutes(Tnow-Tnode)/to_minutes(Torb))
    mapcount=orbid/Rcyc
    shift=orbid/Rcyc
    Tearth=timeperiod(hh=24)
    deltalon=360/to_seconds(Tearth)*to_seconds(Torb)
    shiftlon=deltalon/shift
    return(shiftlon)
        
def cylcdate(DT):
    year=DT.year
    month=DT.month
    day=DT.day
    hour=DT.hour
    minute=DT.minute
    second=DT.second
    return(str_zfill(year,4)+str_zfill(month,2)+str_zfill(day,2)+"T"+str_zfill(hour,2)+str_zfill(minute,2)+"Z")
    
def odbdate(DT):
    year=DT.year
    month=DT.month
    day=DT.day
    hour=DT.hour
    minute=DT.minute
    second=DT.second
    return(str_zfill(year,4)+str_zfill(month,2)+str_zfill(day,2))

def odbtime(DT):
    year=DT.year
    month=DT.month
    day=DT.day
    hour=DT.hour
    minute=DT.minute
    second=DT.second
    return(str_zfill(hour,2)+str_zfill(minute,2)+str_zfill(second,2))

def cylcdate_get_hour(cylc):
	hour=str_zfill(int(cylc[9:11]),2)
	return(hour)

def pydatetime(year,month,day,hour=0,minute=0,second=0):
    return(datetime.datetime(year,month,day,hour,minute,second))

def pydate(date=None, year=None, month=None, day=None, hour=00, minute=00, second=00):
    if date is not None :
	datestr=str(date)
	year=int(datestr[0:4])
	month=int(datestr[4:6])
	day=int(datestr[6:8])
    hour=int(hour)
    minute=int(minute)
    second=int(second)
    return(pydatetime(year,month,day,hour,minute,second))
  
def fmtdatetime(fmtstr,date=None, year=None, month=None, day=None, hour=00, minute=00, second=00):
	DT=pydate(date,year,month,day,hour,minute,second)
	dtfmted=DT.strftime(fmtstr)
	return(dtfmted)
 
def cylcdate_to_pydate(cylc):
	year=int(cylc[0:4])
	month=int(cylc[4:6])
	day=int(cylc[6:8])
	hour=int(cylc[9:11])
	minute=int(cylc[11:13])
	second=0
	return(datetime.datetime(year,month,day,hour,minute,second))
 
def odb_dt_to_pydate(date, time):
  year = int(date[0:4])
  month = int(date[4:6])
  day = int(date[6:8])
  hour = int(time[0:2])
  minute = int(time[2:4])
  second = int(time[4:6])
  return(datetime.datetime(year,month,day,hour,minute,second))

def getdatetime(DT,data,idx=1):
    if "Year" in data.columns.values:
        Year=data.Year.values[idx-1]
    else:
        Year=DT.year
    if "Month" in data.columns.values:    
        Month=data.Month.values[idx-1]
    else:
        Month=DT.month
    if "Day" in data.columns.values:
        Day=data.Day.values[idx-1]
    else:
        Day=DT.day
    if "Hour" in data.columns.values:
        Hour=data.Hour.values[idx-1]
    else:
        Hour=DT.hour
    if "Minute" in data.columns.values:
        Minute=data.Minute.values[idx-1]
    else:
        Minute=DT.minute
    if "Second" in data.columns.values:
        Second=data.Second.values[idx-1]
    else:
        Second=DT.second
    return(datetime.datetime(Year,Month,Day,Hour,Minute,Second))

def get_date_info(data):
	year=int(data.iloc[0].Year)
	month=int(data.iloc[0].Month)
	day=int(data.iloc[0].Day)
	hour=int(data.iloc[0].Hour)
	DT=pydate(None,year,month,day,hour)	
	return(DT)

def pandas_dtfmt(data,dtfmt,fldnamlst=None,fldtype=None):
    if fldtype is not None : 
    	fldnamlst=list(chain(fldtype["year"],fldtype["month"],fldtype["day"],fldtype["hour"],fldtype["minute"],fldtype["second"]))
    else :
	fldtype = {"year":[],"month":[],"day":[],"hour":[],"minute":[],"second":[]}
    	if fldnamlst is None: fldnamlst = data.columns
	if "Year" in fldnamlst : 
		fldnam = "Year"
		fldtype["year"].append(fldnam)
	if "Month" in fldnamlst : 
		fldnam = "Month"
		fldtype["month"].append(fldnam)
	if "Day" in fldnamlst : 
		fldnam = "Day"
		fldtype["day"].append(fldnam)
	if "Hour" in fldnamlst : 
		fldnam = "Hour"
		fldtype["hour"].append(fldnam)
	if "Minute" in fldnamlst : 
		fldnam = "Minute"
		fldtype["minute"].append(fldnam)
	if "Second" in fldnamlst : 
		fldnam = "Second"
		fldtype["second"].append(fldnam)
	fldnamlst=list(chain(fldtype["year"],fldtype["month"],fldtype["day"],fldtype["hour"],fldtype["minute"],fldtype["second"]))
    for fldnam in fldnamlst:
	if fldnam in fldtype["year"]:
		data[fldnam]=pandas.to_datetime(data[fldnam], format=dtfmt).dt.year
	if fldnam in fldtype["month"]:
		data[fldnam]=pandas.to_datetime(data[fldnam], format=dtfmt).dt.month
	if fldnam in fldtype["day"]:
		data[fldnam]=pandas.to_datetime(data[fldnam], format=dtfmt).dt.day
	if fldnam in fldtype["hour"]:
		data[fldnam]=pandas.to_datetime(data[fldnam], format=dtfmt).dt.hour
	if fldnam in fldtype["minute"]:
		data[fldnam]=pandas.to_datetime(data[fldnam], format=dtfmt).dt.minute
	if fldnam in fldtype["second"]:
		data[fldnam]=pandas.to_datetime(data[fldnam], format=dtfmt).dt.second
    return(data)

def datetimeframe(datain,dtfmt):
        data1 = numpy.array(datain).flatten()
	fldnam="date"
	data=pandas.DataFrame(data1,columns=[fldnam])
	dfnew=DataFrame()
	dfnew["year"]=pandas.to_datetime(data[fldnam], format=dtfmt).dt.year
	dfnew["month"]=pandas.to_datetime(data[fldnam], format=dtfmt).dt.month
	dfnew["day"]=pandas.to_datetime(data[fldnam], format=dtfmt).dt.day
	dfnew["hour"]=pandas.to_datetime(data[fldnam], format=dtfmt).dt.hour
	dfnew["minute"]=pandas.to_datetime(data[fldnam], format=dtfmt).dt.minute
	dfnew["second"]=pandas.to_datetime(data[fldnam], format=dtfmt).dt.second
	return(dfnew)

def dtextract(dtframe,indx=0):
	year=int(dtframe.year.values[indx])
	month=int(dtframe.month.values[indx])
	day=int(dtframe.day.values[indx])
	hour=int(dtframe.hour.values[indx])
	minute=int(dtframe.minute.values[indx])
	second=int(dtframe.second.values[indx])
	return(datetime.datetime(year,month,day,hour,minute,second))

def pandas_strcrop(data,field,chaaa,chzzz):
	data[field]=data[field].str[int(chaaa):int(chzzz)].values

def get_numerics(string_data):
    numeric_data=[]
    num_string = ""
    for c in string_data:
        if c.isdigit():
           num_string = num_string + c
        else :
           if len(num_string) > 0:
              numeric_data = numeric_data + [int(num_string)]
              num_string = ""
    return(numeric_data)

def dtlist(dtstart,tdelta,count):
	dtli=[]
	for i in range(0,count,1):
		dtli.append(dtstart+(i*tdelta))
	return(dtli) 
        
def dtlist_backward(dtend,tdelta,count):
	count=int(count)
	tspan=tdelta*(count-1)
	dtstart=dtend-tspan
	print(dtstart,tspan, dtend)
	dtli=dtlist(dtstart,tdelta,count)
	return(dtli)
        
def obs_merge_batch(dataframelist):
    print(dataframelist)
    for df in dataframelist:
	print(len(df))
	print(df.columns)
    return(pandas.concat(dataframelist, ignore_index=True))

def obsdfcat(dataframelist):
    return(pandas.concat(dataframelist, axis='columns'))

def to_string(data,dtype=None):
    if isinstance(data,numpy.ndarray):
       if dtype is 'float_kind': data=numpy.array2string(data, precision=2, separator=',', threshold=100000, formatter={'float_kind':lambda x: "%.2f" % x})
       if dtype is 'int' : data=numpy.array2string(data, precision=2, separator=',', threshold=100000, formatter={'int':lambda x: "%0>5d" % x})
    else :
       if isinstance(data, tuple) :
          data=str(data)
       elif isinstance(data, pandas.Series):
	  data=str(data.values)
       else:
          if isinstance(data, (int,float)) :
             data=str(data)
          else :
             print("Data type not supported")
    return(data)

def print_frame(data,option=0):
    print("Option : "+str(option))
    if option in [0] : 
        print(data)
    if option in [1] : 
        if isinstance(data, pandas.DataFrame):
           for i in data.index.values: print(data.xs(i))
        else :
           if isinstance(data,numpy.ndarray):
              for i in range(0,len(data),1):
                print(data[i])
           else :
	      print(to_string(data))
    if option in [2] :
        if isinstance(data, pandas.DataFrame):
           for i in data.index.values:
              for j in data.columns.values:
                 print(i,j,data.xs(i)[j])
        else :
           if isinstance(data,numpy.ndarray):
              for i in range(0,len(data),1):
                 for j in range(0,len(data[i]),1):
                    print(data[i,j])
           else :
	      print(to_string(data))

def mkdir(path):
    paretdir=path.rsplit("/",1)[0]
    if not os.path.isdir(paretdir):
       mkdir(paretdir)
    if not os.path.isdir(path) :
       os.mkdir(path)

def obs_frame_ascii(data,outfile=None,option=0):
    if outfile is None:
        print_frame(data,option=option)
    else:
        filedir=outfile.rsplit("/",1)[0]
        mkdir(filedir)
        print(filedir) 
        with open(outfile,"w") as op:
          if isinstance(data, pandas.DataFrame):
             data.to_string(op)
          else :
             data=to_string(data)
             op.write(data)
    return(outfile)

def getopsname(nmlfile,odbname):
    with open(nmlfile, "r") as nml:
        try:opsname = str(pandas.read_table(nml, skiprows=None, header=0).query("odbname == @odbname").opsname.values[0])
        except:opsname=str(odbname)
    return(opsname)

def getodbname(nmlfile,opsname):
    with open(nmlfile, "r") as nml:
        try:odbname = str(pandas.read_table(nml, skiprows=None, header=0).query("opsname == @opsname").odbname.values[0])
        except:odbname=str(opsname)
    return(odbname)

def get_subtype_name(nmlfile,subtype):
    with open(nmlfile, "r") as nml:
        try:name = str(pandas.read_table(nml, skiprows=None, header=0).query("subtype == @subtype").obstypnam.values[0])
        except:name=str(subtype)
    return(name)

def get_subtype_code(nmlfile,obstypnam):
    with open(nmlfile, "r") as nml:
        try:stypcode = int(pandas.read_table(nml, skiprows=None, header=0).query("obstypnam == @obstypnam").subtype.values[0])
        except:stypcode=int(0)
    return(stypcode)

def dfheader(data):
    cnt=len(data.columns.values)
    data_header=pandas.DataFrame(data.columns.values,index=range(1,cnt+1,1))
    return(data_header.values[:,0])
    
def odb_renamefield(data,nmlfile):
    elenams=data.columns.values
    for odbname in elenams:
            opsname=getopsname(nmlfile,odbname)
            data=data.rename(index=str,columns={str(odbname):str(opsname)})
    return(data)

def clevgen(long_name,field="Obsvalue"):
    if field in ["Obsvalue"]:
        if long_name.split("_")[0] in ["Pressure"]: 
            if long_name.split("_")[1] in ["MSL"]: clevs=numpy.arange(100000, 105000, 100)
            else: clevs=numpy.arange(70000, 100000, 1000)
        elif long_name.split("_")[0] in ["WindDirection"]: clevs=numpy.arange(0, 360, 10)
        elif long_name.split("_")[0] in ["WindSpeed"]: clevs=numpy.arange(0, 20, 2)
        elif long_name.split("_")[0] in ["Uwind","Vwind"]: clevs=numpy.arange(-5, 5, 0.5)
        elif long_name.split("_")[0] in ["Visibility"]: clevs=numpy.arange(0, 1000, 50)
        elif long_name.split("_")[0] in ["RelativeHumidity"]: clevs=numpy.arange(0, 100, 5)
        elif long_name.split("_")[0] in ["Temperature"]: clevs=numpy.arange(260, 300, 1)
        else: clevs=numpy.arange(260, 300, 1) 
    elif field in ["Depar","FGDep","AnalDep"]:
        if long_name.split("_")[0] in ["Pressure","WindDirection"]: clevs=numpy.arange(-100, 100, 10)
        elif long_name.split("_")[0] in ["Visibility"]: clevs=numpy.arange(-1000, 1000, 10)
        elif long_name.split("_")[0] in ["WindSpeed","Temperature","Uwind","Vwind"]: clevs=numpy.arange(-2, 2, 0.2)
        elif long_name.split("_")[0] in ["RelativeHumidity","CloudAmount"]: clevs=numpy.arange(-10, 10, 1)
        else: clevs=numpy.arange(-10, 10, 1) 
    elif field in ["ObsDensity"]: clevs=numpy.arange(0, 10, 1)      #
    else: errprint("No rules to generate contour levels")
    return(clevs)

def dataunit(long_name):
    if long_name.split("_")[0] in ["Pressure"]: dataunit="Pa"
    else: dataunit=""
    return(dataunit)

def getvarno(varno_nmlfile,varnam):
    with open(varno_nmlfile, "r") as nml:
        varno = str(pandas.read_table(nml, skiprows=None, header=0).query("varnam == @varnam").varno.values[0])
        #if varname in ["--","-",""]: varname="varno_"+str(varno)
    return(varno)

def getvarname(varno_nmlfile,varno):
    with open(varno_nmlfile, "r") as nml:
        varname = str(pandas.read_table(nml, skiprows=None, header=0).query("varno == @varno").varnam.values[0])
        if varname in ["--","-",""]: varname="varno_"+str(varno)
    return(varname)

def listvarname(varno_nmlfile,varnolist):
    return([getvarname(varno_nmlfile,varno) for varno in varnolist])

def listvarno(varno_nmlfile,varnamlist):
    return([getvarno(varno_nmlfile,varnam) for varnam in varnamlist])

def getlongname(varno_nmlfile,varno):
    with open(varno_nmlfile, "r") as nml:
        long_name = str(pandas.read_table(nml, skiprows=None, header=0).query("varno == @varno").longnam.values[0])
        if long_name in ["--","-",""]: long_name="varno_"+str(varno)
    return(long_name)

def get_subtype_name(subtype_nmlfile,subtype):
    with open(subtype_nmlfile, "r") as nml:
        subtype_name = str(pandas.read_table(nml, skiprows=None, header=0).query("subtype == @subtype").stname.values[0])
        if subtype_name in ["--","-",""]: subtype_name="subtype_"+str(subtype)
    return(subtype_name)

def gridded_count_1x1deg(obsframe,varname=None):
    nrows=len(obsframe)
    latmin=int(-90)
    latmax=int(90)
    lonmin=int(-180)
    lonmax=int(180)
    glat=range(latmin,latmax,1)
    glon=range(lonmin,lonmax,1)
    gridded_count = numpy.zeros((len(glat),len(glon)), dtype=numpy.int)
    for row in obsframe.index:
        gridded_count[int(obsframe.Latitude[row])-latmin-1,int(obsframe.Longitude[row])-lonmin-1] += 1
    gridded_count[gridded_count==0]=-99999
    datacount=pandas.DataFrame(gridded_count,columns=glon,index=glat)
    return(datacount)
    
def gridded_sum_1x1deg(obsframe,varname,fieldname="Obsvalue"):
    nrows=len(obsframe)
    latmin=int(-90)
    latmax=int(90)
    lonmin=int(-180)
    lonmax=int(180)
    glat=range(latmin,latmax,1)
    glon=range(lonmin,lonmax,1)
    gridded_sum = numpy.zeros((len(glat),len(glon)), dtype=numpy.int)
    try: data=obsframe[str(fieldname)]
    except: data=obsframe[str(varname)]
    for row in obsframe.index:
        gridded_sum[int(obsframe.Latitude[row])-latmin-1,int(obsframe.Longitude[row])-lonmin-1] += data[row]
    gridded_sum[gridded_sum==0]=numpy.nan
    data=pandas.DataFrame(gridded_sum,columns=glon,index=glat)
    return(data)

def gridded_mean_1x1deg(obsframe,varname,fieldname="Obsvalue"):
    nrows=len(obsframe)
    latmin=int(-90)
    latmax=int(90)
    lonmin=int(-180)
    lonmax=int(180)
    glat=range(latmin,latmax,1)
    glon=range(lonmin,lonmax,1)
    gridded_count = numpy.zeros((len(glat),len(glon)), dtype=numpy.float64)
    gridded_sum = numpy.zeros((len(glat),len(glon)), dtype=numpy.float64)
    try: data=obsframe[str(fieldname)]
    except: data=obsframe[str(varname)]
    for row in obsframe.index:
        try:
            gridded_count[int(obsframe.loc[row,"Latitude"])-latmin-1,int(obsframe.loc[row,"Longitude"])-lonmin-1] += 1
            gridded_sum[int(obsframe.loc[row,"Latitude"])-latmin-1,int(obsframe.loc[row,"Longitude"])-lonmin-1] += data[row]
        except:
           message="Skipped row :"+str(row) +" ("+str(obsframe.loc[row,"Latitude"])+","+str(obsframe.loc[row,"Longitude"])+")"
           print(message)
    gridded_sum[gridded_count==0]=numpy.nan
    gridded_count[gridded_count==0]=numpy.nan
    gridded_mean=numpy.divide(gridded_sum,gridded_count)
    data=pandas.DataFrame(gridded_mean,columns=glon,index=glat)
    return(data)

def gridded_rms_1x1deg(obsframe,varname,fieldname="Obsvalue"):
    nrows=len(obsframe)
    latmin=int(-90)
    latmax=int(90)
    lonmin=int(-180)
    lonmax=int(180)
    glat=range(latmin,latmax,1)
    glon=range(lonmin,lonmax,1)
    gridded_count = numpy.zeros((len(glat),len(glon)), dtype=numpy.float64)
    gridded_sqrsum = numpy.zeros((len(glat),len(glon)), dtype=numpy.float64)
    try: data=obsframe[str(fieldname)]
    except: data=obsframe[str(varname)]
    for row in obsframe.index:
        try:
            gridded_count[int(obsframe.loc[row,"Latitude"])-latmin-1,int(obsframe.loc[row,"Longitude"])-lonmin-1] += 1
            gridded_sqrsum[int(obsframe.loc[row,"Latitude"])-latmin-1,int(obsframe.loc[row,"Longitude"])-lonmin-1] += (data[row]*data[row])
        except:
           message="Skipped row :"+str(row) +" ("+str(obsframe.loc[row,"Latitude"])+","+str(obsframe.loc[row,"Longitude"])+")"
           print(message)
    gridded_sqrsum[gridded_count==0]=numpy.nan
    gridded_count[gridded_count==0]=numpy.nan
    gridded_rms=numpy.sqrt(numpy.divide(gridded_sqrsum,gridded_count))
    data=pandas.DataFrame(gridded_rms,columns=glon,index=glat)
    return(data)

#def gridded_mean_1x1deg(obsframe):
#    nrows=len(obsframe)
#    latmin=int(-90)
#    latmax=int(90)
#    lonmin=int(-180)
#    lonmax=int(180)
#    glat=range(latmin,latmax,1)
#    glon=range(lonmin,lonmax,1)
#    gridded_count = numpy.zeros((len(glat),len(glon)), dtype=numpy.int)
#    #x,y=numpy.meshgrid(glon,glat)
#    for row in range(0,nrows-1,1):
#        gridded_count[int(obsframe.lat@hdr[row])-latmin-1,int(obsframe.lon@hdr[row])-lonmin-1] += 1
#        gridded_sum[int(obsframe.lat@hdr[row])-latmin-1,int(obsframe.lon@hdr[row])-lonmin-1] += obsframe.obsvalue[row]

def get_numerics(string_data):
    numeric_data=[]
    num_string = ""
    for c in string_data:
        if c.isdigit():
           num_string = num_string + c
        else :
           if len(num_string) > 0:
              numeric_data = numeric_data + [int(num_string)]
              num_string = ""
    return(numeric_data)

def unique_int(array,mfactor=1,missing=numpy.nan):
    intlst=[]
    for elem in array:
        if elem != missing : 
           if mfactor > 1 : intlst=intlst+[int(elem/mfactor+0.5)*mfactor]
           else : intlst=intlst+[int(elem)]
    int_frm=pandas.DataFrame(intlst,columns=["Data"])
    int_arr=int_frm.Data.unique()
    return(int_arr)

def obs_clock_hour(obstime,cylchour):
    sec_of_day=int(cylchour)*3600+obstime
    hour_of_day=sec_of_day/3600
    return(clock_24_hour(hour_of_day))

def clock_24_hour(hour):
    if hour > 24 : hour=hour-25
    if hour < 0 : hour=hour+24
    if hour < 0 or hour > 24 : clock_24_hour(hour)
    return(hour)

def get_key_info(nmlfile,key="obsgroup"):
    print(nmlfile)
    nmlinfo=pandas.read_csv(nmlfile, delimiter=': ',engine='python')
    if key in nmlinfo["keys"].values:
    	keyinfo=nmlinfo.query("keys == @key").information.values[0]
    else:
	print("Key '"+str(key)+"' not found")
	print(nmlinfo["keys"].values)
    return(keyinfo)

def get_key_list_info(nmlfile,key):
	strinfo=get_key_info(nmlfile,key)
	lstinfo=list(numpy.fromstring(strinfo[1:-1],sep=',',dtype=int))
	return(lstinfo)
	
def get_key_dic(keyinfofile,keylist,infodic=None):
   #print(keyinfofile,keylist,infodic)
   if infodic is None : infodic={}
   if os.path.exists(keyinfofile): 
	#print("Readimng "+keyinfofile)
	for key in keylist:
		infodic[key]=get_key_info(keyinfofile,key)
   else:
	#print("File not found: "+keyinfofile)
	srcdic=obsdic.obstype[obstype]
	infodic=srcdic.copy()
   return(infodic)

def get_elist(obstypnam,obstypnml,keynml):
	subtype=get_subtype_code(obstypnml,obstypnam)
	elist=get_key_list_info(keynml,["elemlist_"+str(subtype)])
	elist=list(elist)
	return(elist)

def datframe_add_subtype(data,obstypnam,obstypnml):
	subtype=get_subtype_code(obstypnml,obstypnam)
	data=data.assign(subtype=[int(subtype)]*len(data))
	return(data)

def reset_index(data,index=None):
    if index is None :
	count=len(data)
	index=pandas.Series(range(1,(count+1),1))
    data=data.set_index(index)
    return(data)

def data_thinning(data,cntmax=500000,callsign=None,fillval=NAN_VAL):
	if callsign is not None :
		validindx=data[data[callsign] != fillval].index
		data=data.loc[validindx,:]
	print(data)
	obscnt=len(data)
	if obscnt > cntmax:
		skpintvl=int(obscnt/cntmax)+1
		fltrindx=pandas.Series(range(1,obscnt,skpintvl))
		data_new=data.loc[fltrindx,:].copy(deep=True)
		data=reset_index(data_new,index=None)
	return(data)

