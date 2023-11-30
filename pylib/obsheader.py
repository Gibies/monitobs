#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 19:32:38 2020

@author: gibies
"""
from __future__ import print_function
import sys
import os
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
PKGHOME=os.path.dirname(CURR_PATH)
OBSLIB=os.environ.get('OBSLIB',PKGHOME+"/pylib")
sys.path.append(OBSLIB)
OBSDIC=os.environ.get('OBSDIC',PKGHOME+"/pydic")
sys.path.append(OBSDIC)
OBSNML=os.environ.get('OBSNML',PKGHOME+"/nml")
sys.path.append(OBSNML)
import obslib
import fixheader
import obstore
import numpy
import itertools
MAXINDX=int(os.environ.get('MAXINDX',fixheader.MAXINDX))
HDRSIZE=int(os.environ.get('HDRSIZE',fixheader.HDRSIZE))
LUTSIZE=int(os.environ.get('LUTSIZE',fixheader.LUTSIZE))
HBpos=int(os.environ.get('HBpos',fixheader.HBpos))
HBlen=int(os.environ.get('HBlen',fixheader.HBlen))
HCpos=int(os.environ.get('HCpos',fixheader.HCpos))
HClen=int(os.environ.get('HClen',fixheader.HClen))
HDR20=fixheader.HDR20
HDRgam=fixheader.HDRgam
FIXHDR=fixheader.FIXHDR

#def getalphaheader(alphalen):
#    alpha=numpy.empty(shape=[alphalen])
#    for i,hdr in enumerate(fixheader.headerlist[0:alphalen]):
#        alpha[i]=fixheader.headerdic[hdr]
#    alpha=alpha.tolist()
#    return(alpha)

def getalphaheader(alphalen,callsignflag=False):
    alpha=numpy.empty(shape=[alphalen])
    for i,hdr in enumerate(fixheader.headerlist[0:alphalen]):
        alpha[i]=fixheader.headerdic[hdr]
    if callsignflag :
	print(callsignflag)
	alpha[0:12]=fixheader.headerdic_with_callsign[fixheader.headerlist[0:12]]
    alpha=alpha.tolist()
    return(alpha)

def readalphaheader(input_file,alphalen):
    alphapos=1
    with open(input_file, "rb") as infile:
        alpha=obstore.obstore_read_header(infile,alphapos,alphalen)
    return(alpha)
#def readalphaheader(input_file,alphalen):
#    alphapos=1
#    with open(input_file, "rb") as infile:
#        alpha=obstore_read_header(infile,alphapos,alphalen)
#    return(alpha)
#    
    
def checkalpha(input_file,alphalen,deepcheck):
    alphalocal=readalphaheader(input_file,alphalen)
    alphafix=getalphaheader(alphalen)
    #print(alphalocal,alphafix)
    if deepcheck:
        obslib.compare(alphalocal,alphafix)
    checklist=obslib.arraydiff(alphalocal,alphafix)
    print(checklist)
    details=list()
    values=list()
    for i,pos in enumerate(checklist):
        details=details + [fixheader.headerlist[pos]]
        values=values + [str(alphalocal[pos])]
    return(details)
#def checkalpha(input_file,alphalen,deepcheck):
#    alphalocal=readalphaheader(input_file,alphalen)
#    alphafix=getalphaheader(alphalen)
#    #print(alphalocal,alphafix)
#    if deepcheck:
#        obslib.compare(alphalocal,alphafix)
#    checklist=obslib.arraydiff(alphalocal,alphafix)
#    print(checklist)
#    details=list()
#    values=list()
#    for i,pos in enumerate(checklist):
#        print(pos,obsheader.headerlist[pos])
#        details=details + [obsheader.headerlist[pos]]
#        values=values + [str(alphalocal[pos])]
#    return(details)

def getbeetaheader(alphalen):
    alpha=getalphaheader(alphalen)
    beetapos=int(alpha[100-1])
    beetalen=int(alpha[101-1])
    beeta=numpy.empty(shape=[beetalen])
    for i,hdr in enumerate(fixheader.headerlist[beetapos:(beetapos+beetalen)]):
        beeta[i]=int(fixheader.headerdic[hdr])
    beeta=beeta.tolist()
    return(beeta)
#def getbeetaheader(alphalen):
#    alpha=getalphaheader(alphalen)
#    beetapos=int(alpha[100-1])
#    beetalen=int(alpha[101-1])
#    beeta=numpy.empty(shape=[beetalen])
#    for i,hdr in enumerate(obsheader.headerlist[alphalen:(alphalen+beetalen)]):
#        beeta[i]=int(obsheader.headerdic[hdr])
#    beeta=beeta.tolist()
#    return(beeta)
#    
    
def readbeetaheader(input_file,alphalen):
    alpha=readalphaheader(input_file,alphalen)
    beetapos=int(alpha[100-1])
    beetalen=int(alpha[101-1])
    with open(input_file, "rb") as infile:
        beeta=obstore.obstore_read_header(infile,beetapos,beetalen)
    return(beeta)
#def readbeetaheader(input_file,alphalen):
#    alpha=readalphaheader(input_file,alphalen)
#    beetapos=int(alpha[100-1])
#    beetalen=int(alpha[101-1])
#    with open(input_file, "rb") as infile:
#        beeta=obstore_read_header(infile,beetapos,beetalen)
#    return(beeta)
#

def checkbeeta(input_file,alphalen,deepcheck):
    alpha=readalphaheader(input_file,alphalen)
    beetapos=int(alpha[100-1])
    #print(alpha)
    beetalocal=readbeetaheader(input_file,alphalen)
    beetafix=getbeetaheader(alphalen)
    if deepcheck:
        obslib.compare(beetalocal,beetafix)
    checklist=obslib.arraydiff(beetalocal,beetafix)
    print(checklist)
    details=list()
    values=list()
    for i,pos in enumerate(checklist):
        abspos=beetapos+pos
        details=details + [fixheader.headerlist[abspos]]
    return(details)
#def checkbeeta(input_file,alphalen,deepcheck):
#    alpha=readalphaheader(input_file,alphalen)
#    #print(alpha)
#    beetalocal=readbeetaheader(input_file,alphalen)
#    beetafix=getbeetaheader(alphalen)
#    if deepcheck:
#        obslib.compare(beetalocal,beetafix)
#    checklist=obslib.arraydiff(beetalocal,beetafix)
#    print(checklist)
#    details=list()
#    values=list()
#    for i,pos in enumerate(checklist):
#        abspos=alphalen+pos
#        print(abspos,obsheader.headerlist[abspos],beetalocal[pos],beetafix[pos])
#        details=details + [obsheader.headerlist[abspos]]
#    return(details)
#    
    
def getgamaheader(alphalen):
    alpha=getalphaheader(alphalen)
    gamapos=int(alpha[105-1])
    gamalen=int(alpha[106-1])
    gama=numpy.empty(shape=[gamalen])
    for i,hdr in enumerate(fixheader.headerlist[gamapos:(gamapos+gamalen)]):
        gama[i]=int(fixheader.headerdic[hdr])
    gama=gama.tolist()
    return(gama)
#def getgamaheader(alphalen):
#    alpha=getalphaheader(alphalen)
#    gamapos=int(alpha[105-1])
#    gamalen=int(alpha[106-1])
#    gama=numpy.empty(shape=[gamalen])
#    for i,hdr in enumerate(obsheader.headerlist[(gamapos-1):(gamapos+gamalen-1)]):
#        gama[i]=float(obsheader.headerdic[hdr])
#    gama=gama.tolist()
#    return(gama)
#

def readgamaheader(input_file,alphalen):
    alpha=readalphaheader(input_file,alphalen)
    gamapos=int(alpha[105-1])
    gamalen=int(alpha[106-1])
    with open(input_file, "rb") as infile:
        gama=obstore.obstore_read_real(infile,gamapos,gamalen)
    return(gama)
#def readgamaheader(input_file,alphalen):
#    alpha=readalphaheader(input_file,alphalen)
#    gamapos=int(alpha[105-1])
#    gamalen=int(alpha[106-1])
#    with open(input_file, "rb") as infile:
#        gama=obstore_read_real(infile,gamapos,gamalen)
#    return(gama)
#

def checkgama(input_file,alphalen,deepcheck):
    alpha=getalphaheader(alphalen)
    gamapos=int(alpha[105-1])
    gamalocal=readgamaheader(input_file,alphalen)
    gamafix=getgamaheader(alphalen)
    if deepcheck:
        obslib.compare(gamalocal,gamafix)
    checklist=obslib.arraydiff(gamalocal,gamafix)
    print(checklist)
    details=list()
    values=list()
    for i,pos in enumerate(checklist):
        abspos=gamapos+pos
        details=details + [fixheader.headerlist[abspos]]
    return(details)
#def checkgama(input_file,alphalen,deepcheck):
#    alpha=getalphaheader(alphalen)
#    beetalen=int(alpha[101-1])
#    gamalocal=readgamaheader(input_file,alphalen)
#    gamafix=getgamaheader(alphalen)
#    if deepcheck:
#        obslib.compare(gamalocal,gamafix)
#    checklist=obslib.arraydiff(gamalocal,gamafix)
#    print(checklist)
#    details=list()
#    values=list()
#    for i,pos in enumerate(checklist):
#        abspos=alphalen+beetalen+pos
#        print(abspos,obsheader.headerlist[abspos],gamalocal[pos],gamafix[pos])
#        details=details + [obsheader.headerlist[abspos]]
#    return(details)
#

def read_obsheader(input_file):
    alphalen=256
    alpha=readalphaheader(input_file,alphalen)
    beeta=readbeetaheader(input_file,alphalen)
    gama=readgamaheader(input_file,alphalen)
    print(len(alpha),len(beeta),len(gama))
#def read_obsheader(input_file):
#    alphalen=256
#    alpha=readalphaheader(input_file,alphalen)
#    beeta=readbeetaheader(input_file,alphalen)
#    gama=readgamaheader(input_file,alphalen)
#    print(len(alpha),len(beeta),len(gama))
#

def header_diffcheck(input_file,deepcheck=False):
    alphalen=256
    difflist=checkalpha(input_file,alphalen,deepcheck)
    difflist=difflist+checkbeeta(input_file,alphalen,deepcheck)
    difflist=difflist+checkgama(input_file,alphalen,deepcheck)
    print(difflist)
    return(difflist)
#def header_diffcheck(input_file,deepcheck=False):
#    alphalen=256
#    difflist=checkalpha(input_file,alphalen,deepcheck)
#    difflist=difflist+checkbeeta(input_file,alphalen,deepcheck)
#    difflist=difflist+checkgama(input_file,alphalen,deepcheck)
#    return(difflist)
#


#def write_obsheader(outfile):
#    alphalen=256
#    alpha=getalphaheader(alphalen)
#    beeta=getbeetaheader(alphalen)
#    gama=getgamaheader(alphalen)
#    obstore.obstore_write_formatted(itertools.chain(alpha,beeta,gama),">"+str(len(alpha))+"q"+str(len(beeta))+"q"+str(len(gama))+"d",outfile)
#    return(0)

def write_obsheader(obsfile,nmlfile,obsgroup,maxindx=MAXINDX,callsignflag=False):
    alphalen=256
    alpha=getalphaheader(alphalen,callsignflag)
    beeta=getbeetaheader(alphalen)
    gama=getgamaheader(alphalen)
    datafmt=">"+str(len(alpha))+"q"+str(len(beeta))+"q"+str(len(gama))+"d"
    obstore.obstore_write_formatted(itertools.chain(alpha,beeta,gama),datafmt,obsfile)
    obstore.obstore_write_obsgroup(obsgroup,obsfile)
    hdrlen=len(alpha)+len(beeta)+len(gama)
    obslib.binary_write(FIXHDR,1,obsfile)
    obslib.binary_write(HDR20,1,obsfile)
    obslib.binary_write(HDRgam,306,obsfile)
    return(hdrlen)

def write_batchheader(outfile,batchcount=1,header_offset=339,maxindx=MAXINDX,lut_ncols=LUTSIZE):
    obstore.obstore_set_batchpos(outfile,batchcount,header_offset,maxindx,lut_ncols)
    return(0)
