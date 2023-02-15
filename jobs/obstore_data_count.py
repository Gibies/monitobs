#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 02 11:02:15 2022

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

obs_file="/home/umprod/cylc-run/PS43_Hybrid/share/cycle/20220901T0000Z/glu_obstore/ABIClr.obstore"
outfile="/scratch/"+USER+"/test/goes_abi_latlon.dat"
obsmod.obstore_read_latlon(obs_file,outfile)


