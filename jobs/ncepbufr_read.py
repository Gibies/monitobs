
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

import ncepbufr

obsname="grndgps"
PDY="20180530"
TDATE=os.environ.get('PDY',PDY)
print(TDATE)
Tnode=obsmod.pydate(TDATE)

if len(sys.argv) > 1:
	pathstring=sys.argv[1]
else:
	print("Input bufr file is not provided")

nmlfile=os.environ.get('KEYNMLFILE',OBSNML+"/keys_"+obsname+".nml")
txtfile="/scratch/"+USER+"/print/bufr/dump/"+obsname+".txt"

data=obsmod.ncbufr_test_read(slctstr=pathstring,txtfile=txtfile)

#data=obsmod.ncbufr_decode_files("",Tnode,slctstr=pathstring,nmlfile=nmlfile)
