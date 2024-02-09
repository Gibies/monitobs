import sys
sys.path.append("/home/meena/packages/monitobs/pylib")

import obsmod
#import obsplot

datapath="/scratch/meena/cylc-run/PS45_withSatTCWV/share/cycle/20230601T0000Z/glm_obstore"
plotpath="/scratch/meena/plots"
nmlpath="/home/meena/packages/obspreproc/opps_sattcwv/nml"
obstypelist = ["sattcwv"]
text="20230601T0000Z SatTCWV"
title="20230601T0000Z SatTCWV"

obsmod.obs_latlon_plot(datapath,plotpath,nmlpath,obstypelist=obstypelist,text=text,title=title)
