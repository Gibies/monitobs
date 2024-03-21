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

import obslib
import essio


def nextday(cylcdate):
	currdate=obslib.pydate(cylcdate=cylcdate)
	nextdate=currdate+obslib.timeperiod(hh=24)
	return(obslib.cylcdate(nextdate))

def get_cylc_file_list(filepath,filenam,cylcstrt,cylcfinl,cylcinvl):
	dtstrt=obslib.pydate(cylcdate=cylcstrt)
	dtfinl=obslib.pydate(cylcdate=cylcfinl)
	dtinvl=obslib.timeperiod(hh=cylcinvl)
	filelst=[]
	dtrecrd=[]
	dtarray=[]
	DT=dtstrt
	while (DT < dtfinl):
		cylcdt=obslib.cylcdate(DT)
		timval=obslib.to_units(DT)
		dtarray=dtarray+[DT]
		filelst=filelst+[filepath+"/"+str(cylcdt)+"/"+filenam]
		dtrecrd=dtrecrd+[timval]
		DT = DT + dtinvl
	data_info={}
	data_info.update({"filelst":filelst})
	data_info.update({"dtrecrd":dtrecrd})
	data_info.update({"dtarray":dtarray})
	return(data_info)


def datset_read_daily(filepath,filenam,date,cylcinvl,recdim,varlst,dimlst):
	cylcstrt=str(date)+"T0000Z"
	cylcfinl=nextday(cylcdate=cylcstrt)
	data_info=get_cylc_file_list(filepath,filenam,cylcstrt,cylcfinl,cylcinvl)
	infile=data_info["filelst"]
	recrds=data_info["dtrecrd"]
	reclen=len(recrds)
	datset=essio.datset_append(infile,recdim=recdim,recrds=recrds,reclen=reclen,varlst=varlst,dimlst=dimlst)
	return(datset)

