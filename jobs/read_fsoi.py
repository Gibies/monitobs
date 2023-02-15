

import sys,os
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
CYLCROOT=os.path.dirname(os.path.dirname(os.path.dirname(CURR_PATH)))
CYLCPATH=os.environ.get('CYLCPATH',CYLCROOT)
HOMEDIR=os.environ.get('HOMEDIR',CYLCPATH+"/modules/monitobs")
OBSLIB=os.environ.get('OBSLIB',HOMEDIR+"/pylib")
OBSDIC=os.environ.get('OBSLIB',HOMEDIR+"/pydic")
OBSNML=os.environ.get('OBSNML',HOMEDIR+"/nml")
sys.path.append(OBSLIB)
sys.path.append(OBSDIC)
sys.path.append(OBSNML)
import fsoi
import obslib
import obsmod

USER=os.environ.get('USER',"")
SCRATCH='/scratch/'+USER

taskname=os.environ.get('CYLC_TASK_NAME',"glu_fsoi_plot")
varnml=os.environ.get('VarNml',OBSNML+"/varobs_nml")
#stnlst_path=os.environ.get('StnLstDir',CYLCPATH+"/share/data/etc/stationlists/atmos")
cylctime=os.environ.get('CYLC_TASK_CYCLE_POINT',"20210801T0000Z")
inputdir=os.environ.get('INPUTDIR',CYLCPATH+"/share/data/fsoimpact")
outputdir=os.environ.get('OUTPUTDIR',CYLCPATH+"/share/cycle/"+cylctime+"/fsoi_plots")
wrkdir=os.environ.get('WRKDIR',SCRATCH+"/work/"+taskname)


yyyymmdd=cylctime[0:8]
hh=cylctime[9:11]

search_string=os.environ.get('FILE_SEARCH_STRING',inputdir+"/"+cylctime+"/"+yyyymmdd+hh+".FSO")
print(search_string)
subtype=os.environ.get('subtype',"00000")
infile=inputdir+"/"+cylctime+"/"+yyyymmdd+hh+".FSO"
outfile=outputdir+"/data_"+subtype+"_"+yyyymmdd+"T"+hh+"00Z.fso"
print(infile,outfile)
print(search_string,outputdir,wrkdir,yyyymmdd,hh)

data=fsoi.read_fsoi(infile,outputdir,wrkdir,yyyymmdd,hh)
nmldata=fsoi.fsoi_station_nml(data)
print(nmldata)

#subtype_list=[20700]
#for subtype in subtype_list:
#    subtype_data=fsoi.fsoi_subtype_read(subtype,infile,outputdir,wrkdir,yyyymmdd,hh)
#    outfile=outputdir+"/data_"+str(subtype)+"_"+str(yyyymmdd)+"T"+str(hh)+"00Z.fso"
#    obsmod.ascii_file_write(subtype_data,outfile)
#    print(subtype)

#data=fsoi.read_fsoi(infile,outputdir,wrkdir,yyyymmdd,hh)

#print(data)
#subtype_list=data.subtype.unique()
#for subtype in subtype_list:
#    subtype_data=obslib.frame_select_filter(data,"subtype",int(subtype))
#print(fsoi.fsoi_subtype_nml(["SYNOP","SURF","TEMP","BUOY","HLOS"],infile,outputdir,wrkdir,yyyymmdd,hh))
#subtype_list=fsoi.fsoi_subtype_list(infile,outputdir,wrkdir,yyyymmdd,hh)
    

