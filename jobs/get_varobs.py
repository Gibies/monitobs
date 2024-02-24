import sys
import os
diaglev=int(os.environ.get('GEN_MODE',0))
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
PKGHOME=os.path.dirname(CURR_PATH)
OBSLIB=os.environ.get('OBSLIB',PKGHOME+"/pylib")
sys.path.append(OBSLIB)
OBSDIC=os.environ.get('OBSDIC',PKGHOME+"/pydic")
sys.path.append(OBSDIC)
OBSNML=os.environ.get('OBSNML',PKGHOME+"/nml")
sys.path.append(OBSNML)
varobs_nml=OBSNML+"/varobs_nml"
import varobs


infilename="/home/gibies/data/ose/PS45_withSatTCWV/share/cycle/20230601T0000Z/glm_varobs/SatTCWV.varobs"
latlon=varobs.read_latlon(infilename,btchid=1)
print(latlon)
