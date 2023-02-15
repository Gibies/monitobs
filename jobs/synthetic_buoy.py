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

syn_buoy_loc=OBSNML+"/osse_synbuoy.nml"
REFDATE=os.environ.get('PDY',"20200609")

TDATE=os.environ.get('PDY',"20210101")
Tnode=obsmod.pydate(TDATE)
datacfldrname=obsmod.cylcdate(Tnode)
inpath="/home/gibies/init_cdas/obs_atmos/"+obsmod.cylcdate(obsmod.pydate(REFDATE))
#inpath="/home/umprod/cylc-run/PS43_Hybrid/share/cycle/"+datacfldrname+"/glm_obstore"

obs_index_max=512
obstypelist=["surface"]
subtypelist=[10300,]	 #, 10800, 11100]
#obstypelist=["sonde"]
outpath="/home/"+USER+"/data/research/osse_buoy"

obsmod.symobs_buoy_nio(Tnode,outpath,inpath,nmlpath=OBSNML,obs_index_max=obs_index_max,obstypelist=obstypelist,subtypelist=subtypelist,synbuoyloc=syn_buoy_loc)

plotpath="/home/"+USER+"/plots/research/osse_buoy/infile"
#obsmod.obs_latlon_plot(inpath,plotpath,nmlpath=OBSNML,obs_index_max=obs_index_max,obstypelist=obstypelist)
plotpath="/home/"+USER+"/plots/research/osse_buoy/outfile"
obsmod.obs_latlon_plot(outpath,plotpath,nmlpath=OBSNML,obs_index_max=obs_index_max,obstypelist=obstypelist)

#obstype=obstypelist[0]
#obstypedic=obsdic.obstype[obstype]
#filename=obstypedic["filename"]
#obstore_out_file=outpath+"/"+filename
#textfile=outpath+"/"+filename+".txt"
#obsmod.obstore_read_latlon(obstore_out_file,textfile)
#obstype=obstypelist[0]
#obstypedic=obsdic.obstype[obstype]
#filename=obstypedic["filename"]
#with open(inpath+"/"+filename, "rb") as infile:
#	elist=obstore.obstore_read_batch_elements(infile,3,obs_index_max=obs_index_max)
#	#data=obstore.obstore_read_header(infile,10324,100,"d")
#	data=obstore.obstore_read_data_record(infile,3,[1,2,3])
#	print(data)
#	print(elist)
#
#with open(outpath+"/"+filename, "rb") as infile:
#	elist=obstore.obstore_read_batch_elements(infile,3,obs_index_max=obs_index_max)
#	#data=obstore.obstore_read_header(infile,3666,100,"d")
#	data=obstore.obstore_read_data_record(infile,3,[1,2,3])
#	print(data)
#	print(elist)
#
#	data=obstore.obstore_read_header(infile,10324,100,"d")
#	print(data)
#	data=obstore.obstore_read_header(infile,1499475,100,"d")
#	print(data)
#	data=obstore.obstore_read_header(infile,1990046,100,"d")
#	print(data)
#	data=obstore.obstore_read_header(infile,2410488,100,"d")
#	print(data)
##	data=obstore.obstore_read_header(infile,2410499,100,"d")
###	print(data)
