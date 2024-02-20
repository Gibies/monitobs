import os,sys
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
PKGHOME=os.environ.get('PKGHOME',os.path.dirname(CURR_PATH))
PKGNAME=os.path.basename(PKGHOME)
LIB=os.environ.get('LIB',PKGHOME+"/pylib")
sys.path.append(LIB)
DIC=os.environ.get('DIC',PKGHOME+"/pydic")
sys.path.append(DIC)
NML=os.environ.get('NML',PKGHOME+"/nml")
sys.path.append(NML)
import daview

usernam=os.environ.get('USER',"user")

data1="/home/meena/q_20230610_00_ctl"
data2="/home/meena/q_20230610_00_exp"
plotfile="/scratch/"+usernam+"/test.png"

dataset={
	"data1" : data1,
	"data2" : data2,
	}


daview.mpl_plot_indian(dataset,plotfile)
