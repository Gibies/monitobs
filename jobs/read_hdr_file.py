

from __future__ import print_function
import traceback
from eccodes import *

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
print(OBSLIB)
import obsmod
import obstore

if len(sys.argv) > 1:
	pathstring=sys.argv[1]
else:
	pathstring="/home/umprod/cylc-run/PS43_Hybrid/share/cycle/20221104T0000Z/glu_obstore/MWRI.obstore"
inputfile=obsmod.globlist(pathstring)[0]
#maxindx=512
maxindx=621
hdr_info=obstore.obs_hdr_read(inputfile,maxindx=maxindx)
for itm in ["alpha","beeta","gamma","elist"]:
	if itm in hdr_info: print(hdr_info[itm])

#print(alpha)
#print(beeta)
#print(gamma)
#print(LDC)
#print(RDC)
#print(CDC)
#print(LUT)

#with open(inputfile, "rb") as obsfile:
#	data=obstore.obstore_read_data_record(obsfile,1,[1,2,3])
#print(data)

#nmlfile=OBSNML+"/obs_index_nml"
#indx=1
#elenams=["Latitude","Longitude"]
#with open(inputfile, "rb") as obsfile:
#	data=obstore.frame_data_batch(obsfile,nmlfile,indx,elenams,maxindx=512)
#print(data)
