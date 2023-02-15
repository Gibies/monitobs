#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 17:22:59 2021

@author: gibies
"""
from __future__ import print_function
import sys
import os
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
CYLCROOT=os.path.dirname(os.path.dirname(os.path.dirname(CURR_PATH)))
CYLCPATH=os.environ.get('CYLCPATH',CYLCROOT)
MONITOBS=os.environ.get('MONITOBS',CYLCPATH+"/modules/monitobs")
OBSLIB=os.environ.get('OBSLIB',MONITOBS+"/pylib")
sys.path.append(OBSLIB)
OBSDIC=os.environ.get('OBSDIC',MONITOBS+"/pydic")
sys.path.append(OBSDIC)
OBSNML=os.environ.get('OBSDIC',MONITOBS+"/nml")
sys.path.append(OBSNML)
obs_index_nml="obs_index_nml"
nmlfile="%s/%s" % (OBSNML,obs_index_nml)
import obslib
import obsdic
import ncepradic
import imdaadic



def getdata(Year,var="time",element="gph"):
	data=ncepradic.getdata(Year,var,element)
#	data=imdaadic.getdata(Year,var,element)
	return(data)
