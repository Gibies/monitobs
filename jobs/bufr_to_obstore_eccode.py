

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

if len(sys.argv) > 1:
    PDY = sys.argv[1]
else:
    PDY = obsmod.today()

if len(sys.argv) > 2:
    CYC = sys.argv[2]
else:
    CYC = "00"

obsname=os.environ.get('obsname', "")
cntmax=int(os.environ.get('OBS_CNT_MAX', 500000))
btchcnt=int(os.environ.get('BTCH_CNT_MAX',5))
TDATE=os.environ.get('PDY',PDY)
print(TDATE)
Tnode=obsmod.pydate(TDATE)
datacfldrname=obsmod.cylcdate(Tnode)
outpath=os.environ.get('WRK_OBSTORE',"/scratch/"+USER+"/obstore/work/"+obsname+"/"+PDY+"/"+CYC+"/workdir_obstore")
inpath=os.environ.get('BUFRDIR',"")
slctstr=os.environ.get('BUFRFILESTR',"")
nmlfile=os.environ.get('KEYNMLFILE',OBSNML+"/keys_"+obsname+".nml")

data=obsmod.ecbufr_decode_files(inpath,Tnode,slctstr,nmlfile)
obsmod.ascii_file_write(data,outfile=outpath+"/data_"+obsname+".txt",option=1)
outfile=obsmod.obstore_write(data,nmlfile,outpath,btchcnt=btchcnt,cntmax=cntmax,DT=Tnode,diagflag=0)


#print(outfile)

#out=obstore.obstore_read_file(outpath,"seviri")

#print(out.data)
