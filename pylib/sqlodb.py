#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 06:45:11 2019

@author: gibies
"""
from __future__ import print_function
import os
import sys
OBSLIB=os.environ.get('OBSLIB',"${MONITOBS}/pylib")
sys.path.append(OBSLIB)
import obslib
import subprocess
subprocess.call("module unload PrgEnv-cray", shell=True)
subprocess.call("module load PrgEnv-intel/6.0.4", shell=True)
subprocess.call("module load intel/odbserver/0.16.2.omp.1", shell=True)
import odb
import pandas
import numpy
import datetime

diaglev=int(os.environ.get('GEN_MODE',0))
def errprint(*args, **kwargs):
    if diaglev > 0: print(*args, file=sys.stderr, **kwargs)
    
def sqlodb(odbfile,sqlquerystring):
    if os.path.getsize(odbfile) > 0:
        data=pandas.read_sql_query(str(sqlquerystring),odb.connect(odbfile))
        return(data)

def query(odbfile,nmlfile,subtype=None,elenams=[],varnolist=[],userquery=[]):
    print(varnolist)
    odbname=[None]*len(elenams)
    for i,opsname in enumerate(elenams):
            odbname[i]=obslib.getodbname(nmlfile,opsname)
    if not odbname: selectstring = "*"
    else: selectstring = ','.join(odbname+["ops_subtype"])
    if not subtype: subtypequery = ""
    else: subtypequery = "where ops_subtype = "+str(subtype)
    if not userquery: querystring= subtypequery
    else: querystring=' and '.join( [subtypequery]+ userquery )
    varnoquery=queryvarno(varnolist)
    if not varnolist: print(querystring)
    else: querystring= ' AND '.join( [querystring] + varnoquery  )
    data=sqlodb(odbfile,'select ' + selectstring + ' from "' + odbfile + '" ' + querystring + ';')
    return(data)

def queryvarno(varnolist):
    varnoquery=' OR varno = '.join(["varno = "+str(varnolist[0])] + [str(i) for i in varnolist[1:]] )
    return(["( " +varnoquery +")"])
    
def odb_readdata(odbfile):
    return(sqlodb(odbfile,'select * from "' + odbfile + '" ;'))
#    if os.path.getsize(odbfile) > 0:
#        data=pandas.read_sql_query('select * from "' + odbfile + '" ;',odb.connect(odbfile))
#        return(data)
        
def odb_sqlwhere(odbfile,querystring):
    return(sqlodb(odbfile,'select * from "' + odbfile + '" where '+querystring+' ;'))
        
def odb_sqlselect(odbfile,elenams):
    elestr=','.join(elenams)
    return(sqlodb(odbfile,'select '+ elestr +' from "' + odbfile + '" ;'))

def odb_list_varno(odbfile):
    return(sqlodb(odbfile,'select distinct varno from "' + odbfile + '";').values[:,0])

def odb_list_subtype(odbfile):
    return(pandas.read_sql_query('select distinct ops_subtype from "' +odbfile + '" ;',odb.connect(odbfile)).values[:,0])

def odb_filter_varno(odbfile,varno):
    return(sqlodb(odbfile,'select * from "' + odbfile + '" where varno ='+ str(varno) +' ;'))

def odb_readfield(odbfile,nmlfile,element):
    odbname = obslib.getodbname(nmlfile,element)
    data=sqlodb(odbfile,'select '+odbname+' from "' + odbfile + '" ;')
    return(data.rename(index=str,columns={odbname:element}))
    
def obs_frametable(odbfile,odbnmlfile,elenams):
    dataframelist=[None]*len(elenams)
    for i,element in enumerate(elenams):
        print(element)
        dataframelist[i]=odb_readfield(odbfile,odbnmlfile,element)
    return(obslib.obsdfcat(dataframelist))

