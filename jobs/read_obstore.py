

from __future__ import print_function
import traceback
#from eccodes import *

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

print(inputfile)

if len(sys.argv) > 2:
	pos=int(sys.argv[2])
else:
	pos=2003
if len(sys.argv) > 3:
	size=int(sys.argv[3])
else:
	size=19

if len(sys.argv) > 4:
	fmtkey=str(sys.argv[4])
else:
	fmtkey="d"

with open(inputfile, "rb") as obsfile:
	hdr_info=obstore.obstore_read_batch_header(obsfile)
alpha=hdr_info["alpha"]
beeta=hdr_info["beeta"]
gamma=hdr_info["gamma"]
ldc=hdr_info["ldc"]
rdc=hdr_info["rdc"]
cdc=hdr_info["cdc"]
lut=hdr_info["lut"]
print(alpha)
print(beeta)
print(gamma)
print(ldc)
print(rdc)
print(cdc)
print(lut)

with open(inputfile, "rb") as obsfile:
	data=obstore.obstore_read_real(obsfile,pos,size,fmtkey=fmtkey)
	#data=obstore.obstore_read_data_record(obsfile,1,[1,2,3])
print(data)

nmlfile=OBSNML+"/obs_index_nml"
indx=1
elenams=["Latitude","Longitude"]
with open(inputfile, "rb") as obsfile:
	datptr=2003
	datsze=19
	data=obstore.obstore_read_bin8(obsfile,datptr,datsze,sec_nam="data")
	#data=obstore.frame_data_batch(obsfile,nmlfile,indx,elenams,maxindx=512)
#print(data)
