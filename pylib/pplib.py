#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 12:29:36 2020

@author: gibies
"""

import os, sys
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
CYLCROOT=os.path.dirname(os.path.dirname(os.path.dirname(CURR_PATH)))
CYLCPATH=os.environ.get('ROSE_SUITE_DIR',CYLCROOT)
MONITOBS=os.environ.get('MONITOBS',CYLCPATH+"/modules/monitobs")
LIB=os.environ.get('LIB',MONITOBS+"/pylib")
sys.path.append(LIB)
DIC=os.environ.get('DIC',MONITOBS+"/pydic")
sys.path.append(DIC)
NML=os.environ.get('NML',MONITOBS+"/nml")
sys.path.append(NML)


import numpy, pandas

def get_2dlat(lat,shape):
	lat=lat.reshape(len(lat),1)
	lat2d=lat*numpy.ones(shape)
	return(lat2d)

def get_2dlon(lon,shape):
	lon=lon.reshape(1,len(lon))
	lon2d=lon*numpy.ones(shape)
	return(lon2d)
