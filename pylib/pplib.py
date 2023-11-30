#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 12:29:36 2020

@author: gibies
"""

import os, sys
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
PKGHOME=os.path.dirname(CURR_PATH)
OBSLIB=os.environ.get('OBSLIB',PKGHOME+"/pylib")
sys.path.append(OBSLIB)
OBSDIC=os.environ.get('OBSDIC',PKGHOME+"/pydic")
sys.path.append(OBSDIC)
OBSNML=os.environ.get('OBSNML',PKGHOME+"/nml")
sys.path.append(OBSNML)


import numpy, pandas

def get_2dlat(lat,shape):
	lat=lat.reshape(len(lat),1)
	lat2d=lat*numpy.ones(shape)
	return(lat2d)

def get_2dlon(lon,shape):
	lon=lon.reshape(1,len(lon))
	lon2d=lon*numpy.ones(shape)
	return(lon2d)
