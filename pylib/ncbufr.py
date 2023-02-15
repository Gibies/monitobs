#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 28 13:45:39 2022

@author: gibies
"""
from __future__ import print_function
import traceback
import sys,os
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
PKGHOME=os.path.dirname(CURR_PATH)
OBSLIB=os.environ.get('OBSLIB',PKGHOME+"/pylib")
sys.path.append(OBSLIB)
OBSDIC=os.environ.get('OBSDIC',PKGHOME+"/pydic")
sys.path.append(OBSDIC)
OBSNML=os.environ.get('OBSNML',PKGHOME+"/nml")
sys.path.append(OBSNML)
NCBUFRNML=OBSNML+"/ncbufr_fieldname.nml"
import obslib
import ncepbufr
import pandas
import numpy
import collections

#################################################################################
### Common Methords used are as follows.
#################################################################################
#bufr.advance()
#bufr.load_subset()
#bufr.read_subset()
#bufr.rewind()
#bufr.close()
#################################################################################

#################################################################################
def test_read(slctstr=None,txtfile=None,field=None):
   inputfile=obslib.globlist(slctstr)[0]
   if txtfile is not None: obslib.mkdir(txtfile.rsplit("/",1)[0])
   bufr = ncepbufr.open(inputfile)
   while bufr.advance() == 0:
   	print(bufr.msg_counter, bufr.msg_type, bufr.msg_date)
   	while bufr.load_subset() == 0:
		if txtfile is not None: 
			bufr.dump_subset(txtfile)
		else:
            		bufr.print_subset(verbose=verbose)
   		#hdrstr ='YEAR MNTH DAYS HOUR MINU PCCF ELRC SAID PTID GEODU'
   		#hdr = bufr.read_subset(hdrstr).squeeze()
		if field is not None :
   			data = bufr.read_subset(field)[:]
   			print(data)
   bufr.close()
   print(txtfile)
#################################################################################

def get_ncbufr_fieldlist(elemlist=None,nmlfile=None,ncbufrnml=NCBUFRNML):
	if nmlfile is None : nmlfile=OBSNML+"/keys_monitor.nml"
	if elemlist is None : elemlist=obslib.get_key_info(nmlfile,key="elemlist")
	elemlist=numpy.fromstring(elemlist[1:-1],sep=',',dtype=int)
	elist=pandas.DataFrame(collections.Counter(elemlist).items(),columns=["elem","elcnt"])
	obs_fieldlist=[]
	ncb_fieldlist=[]
	obsindxlist=[]
	levlist=[]
	for elem in elist.elem:
		#print(elem)
		fieldname=pandas.read_table(ncbufrnml).query("indx == @elem").fieldname.values[0]
		elename=pandas.read_table(ncbufrnml).query("indx == @elem").elename.values[0]
		count=elist.query("elem == @elem").elcnt.values[0]
		#print(count)
		if count == 1 :
		    if fieldname is not numpy.nan:
			#print(fieldname)
			if elem > 261 :
				ncb_fieldlist=ncb_fieldlist+["#"+str(count)+"#"+fieldname]
			else :
				ncb_fieldlist=ncb_fieldlist+[fieldname]
			obs_fieldlist=obs_fieldlist+[elename]
			obsindxlist=obsindxlist+[elem]
			levlist=levlist+[1]
		else :
			chnlmap=obslib.get_key_info(nmlfile,key=fieldname)
			chnlmap=numpy.fromstring(chnlmap[1:-1],sep=',',dtype=int)
			#print(chnlmap)
			for i,indx in enumerate(chnlmap,start=1):
			    if indx > 0:
				#print(i)
				if fieldname is not numpy.nan: 
					ncb_fieldlist=ncb_fieldlist+["#"+str(indx)+"#"+fieldname]
					obs_fieldlist=obs_fieldlist+[elename+"_"+str(i)]
					obsindxlist=obsindxlist+[elem]
					levlist=levlist+[i]
	fieldlist=pandas.DataFrame({"obsindx":obsindxlist,"chnlev":levlist,"elename":obs_fieldlist,"fieldname":ncb_fieldlist})
	fieldlist=fieldlist.set_index(["obsindx","chnlev"]).sort_index(ascending=True)
	#print(fieldlist)
	return(fieldlist)

def element_read(bufr,field,data,ncbufrnml=NCBUFRNML):
	print(field)
	data[field] = bufr.read_subset(field)
	return(data)


def message_read(bufr,nmlfile,fieldlist=None):
	data=pandas.DataFrame()
	if fieldlist is None:
		fieldlist = get_ncbufr_fieldlist(nmlfile=nmlfile)
	print(fieldlist)
	while bufr.load_subset() == 0:
		for indx,field in fieldlist.iterrows() :
			data=element_read(bufr,field,data)
	print(data)
	return(data)
		

def bufr_decode(infile,nmlfile,fieldlist=None):
	data=pandas.DataFrame()
	bufr = ncepbufr.open(infile)
	while bufr.advance() == 0:
		data1=message_read(bufr,nmlfile,fieldlist=fieldlist)
    		data=data.append(data1, ignore_index=True)
		data=obslib.reset_index(data)
	bufr.close()
	subtype=obslib.get_key_info(nmlfile,"obsubtyp")
	data=data.assign(subtype=[int(subtype)]*len(data))
	return(data)

def bufr_decode_files(inpath,Tnode,slctstr,nmlfile,):
	searchstring=inpath+"/"+slctstr
	infiles=obslib.globlist(searchstring)
	if len(infiles) == 0: print("File not found: "+searchstring)
	print("Received "+str(len(infiles))+" bufr files")
	data=obslib.DataFrame()
	for infile in infiles[:] :
		print(infile)
		data1=bufr_decode(infile,nmlfile)
		data=data.append(data1)
	data=obslib.reset_index(data)
	#print(data)
	return(data)
