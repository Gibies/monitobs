#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 15:17:25 2019

@author: gibies
"""
import sys
import os
USER=os.environ.get('USER',"myhome")
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
PKGHOME=os.path.dirname(CURR_PATH)
OBSLIB=os.environ.get('OBSLIB',PKGHOME+"/pylib")
sys.path.append(OBSLIB)
OBSDIC=os.environ.get('OBSDIC',PKGHOME+"/pydic")
sys.path.append(OBSDIC)
OBSNML=os.environ.get('OBSNML',PKGHOME+"/nml")
sys.path.append(OBSNML)

import obsmod


obs_index_nml="${OBSNML}/obs_index_nml"
odb_index_nml="${OBSNML}/odb_index_nml"

ROSE_SUITE_DIR=""
#obstore_file=ROSE_SUITE_DIR+"/share/cycle/20190110T0000Z/gl_obstore/Surface.obstore"
#odb_file=ROSE_SUITE_DIR+"/share/cycle/20190110T0000Z/glu_odb2/surface.odb"

ascii_out="/scratch/"+USER+"/log/monitobs/test_out.txt"

obj=obsmod.obsguiobj(ROSE_SUITE_DIR=ROSE_SUITE_DIR,OUTFILE=ascii_out)

#print(obj.data)
### Function to write dataframe on asccii file
#obsmod.obs_frame_ascii(asccii_out,obj.data)
