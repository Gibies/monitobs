#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 15:17:25 2019

@author: gibies
"""
import sys
import os
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
CYLCROOT=os.path.dirname(os.path.dirname(CURR_PATH))
CYLCPATH=os.environ.get('CYLCPATH',CYLCROOT)
MONITOBS=os.environ.get('MONITOBS',CYLCPATH+"/modules/monitobs")
OBSLIB=os.environ.get('OBSLIB',MONITOBS+"/pylib")
sys.path.append(OBSLIB)
import obsmod


obs_index_nml="${MONITOBS}/nml/obs_index_nml"
odb_index_nml="${MONITOBS}/nml/odb_index_nml"

ROSE_SUITE_DIR=CYLCPATH
#obstore_file=ROSE_SUITE_DIR+"/share/cycle/20190110T0000Z/gl_obstore/Surface.obstore"
#odb_file=ROSE_SUITE_DIR+"/share/cycle/20190110T0000Z/glu_odb2/surface.odb"

ascii_out=CYLCPATH+"/share/data/test_out.txt"

obj=obsmod.obsguiobj(ROSE_SUITE_DIR=ROSE_SUITE_DIR,OUTFILE=ascii_out)

#print(obj.data)
### Function to write dataframe on asccii file
#obsmod.obs_frame_ascii(asccii_out,obj.data)
