#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu May 21 23:08:37 2020

@author: gibies
"""
from __future__ import print_function
import sys
import os
import re 
import string
import pandas
CYLCPATH=os.environ.get('CYLCPATH',"../../..")
MONITOBS=os.environ.get('MONITOBS',CYLCPATH+"/modules/monitobs")
OBSLIB=os.environ.get('OBSLIB',MONITOBS+"/pylib")
sys.path.append(OBSLIB)
OBSDIC=os.environ.get('OBSDIC',MONITOBS+"/pydic")
sys.path.append(OBSDIC)
OBSNML=os.environ.get('OBSDIC',MONITOBS+"/nml")
sys.path.append(OBSNML)
import bufrdic
import obsdic

obstype={
        "surface" : obsdic.Surface,
        "sonde" : obsdic.Sonde,
        "aircraft" : obsdic.Aircraft,
        "leogeo" : obsdic.LEOGEOAMV,           
        "hlosw" : obsdic.HLOSWind,
        }