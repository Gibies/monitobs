#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 01 11:02:15 2022

@author: gibies
"""

import os,sys

USER=os.environ.get('USER',"myhome")
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
PKGHOME=os.path.dirname(CURR_PATH)
OBSLIB=os.environ.get('OBSLIB',PKGHOME+"/pylib")
sys.path.append(OBSLIB)
OBSDIC=os.environ.get('OBSDIC',PKGHOME+"/pydic")
sys.path.append(OBSDIC)
OBSNML=os.environ.get('OBSNML',PKGHOME+"/nml")
sys.path.append(OBSNML)

#print(sys.path)

import obsmod
import obstore
import obsdic
import datetime

today = datetime.date.today()
TODATE=today.strftime("%Y%m%d")
print(TODATE)
ystrdy=today - datetime.timedelta(days = 1)
YSTRDY=ystrdy.strftime("%Y%m%d")
TDATE=os.environ.get('PDY',YSTRDY)
print(TDATE)
REFDATE=os.environ.get('PDY',"20220920")
Tnode=obsmod.pydate(TDATE)
datacfldrname=obsmod.cylcdate(Tnode)
#inpath="/scratch/gibies/tmpFldr/obstore/output/satwind_umprod_gdas_dump_20220923_00"
#inpath="/scratch/gibies/tmpFldr/obstore/output/grndgps_umprod_gdas_dump_20220828_00"
#inpath="/home/gibies/init_cdas/obs_atmos/"+obsmod.cylcdate(obsmod.pydate(REFDATE))
#inpath="/home/umprod/cylc-run/PS43_Hybrid/share/cycle/"+datacfldrname+"/glu_obstore"
if len(sys.argv) > 1:
    inpath = sys.argv[1]
else :
    print("Input file path is not provided")

maxindx=630	#438	#512
#obstypelist=["surface","sonde","aircraft","goesabi","fy3mws","satwind","grndgps","seviri","mwri"]
obstypelist=["scatwind"]
subtypelist=None
plotpath="/home/"+USER+"/plots/obstore"

obsmod.obs_latlon_plot(inpath,plotpath,nmlpath=OBSNML,maxindx=maxindx,obstypelist=obstypelist)	#,fltrkey="SITE_NAME")

print(plotpath)
