#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 15:45:39 2020

@author: gibies
"""

import os,sys
CURR_PATH=os.path.dirname(os.path.abspath(__file__))
PKGHOME=os.path.dirname(CURR_PATH)
OBSLIB=os.environ.get('OBSLIB',PKGHOME+"/pylib")
sys.path.append(OBSLIB)
OBSDIC=os.environ.get('OBSDIC',PKGHOME+"/pydic")
sys.path.append(OBSDIC)
OBSNML=os.environ.get('OBSNML',PKGHOME+"/nml")
sys.path.append(OBSNML)

from itertools import islice

def display_text(infile):
	with open(infile, 'r') as content_file:
	    content = content_file.read()
	print(content)

def numeric(strng):
	newstr = [character for character in strng if character.isalnum()]
	word="".join(newstr)
	word=int(word) if int(word) == float(word) else float(word)
	return(word)

def alnum(strng):
	alphanumeric = [character for character in strng if character.isalnum()]
	word="".join(alphanumeric)
	return(word)

def get_file_data(infile):
	with open(infile) as f: data = f.readlines()
	return(data)

def get_file_length(infile):
	data=get_file_data(infile)
	return(len(data))

def get_line_number(infile,strng):
	data=get_file_data(infile)
	for line in data: 
		if strng in line: stringline=line
	line_no = data.index(stringline) + 1
	return(line_no)

def get_lines(infile,line_start,num_lines=1):
	start_indx=(line_start-1)
	end_indx=(line_start+num_lines-1)
	data=[]
	with open(infile) as lines: 
	    for line in islice(lines, (start_indx), (end_indx)):
		data.append(line.rstrip("\n"))
	return(data)

def get_visual_block(infile,start_string,end_string=None):
	start_line_no=get_line_number(infile,start_string)
	start_indx=(start_line_no-1)
	data=get_file_data(infile)
	if end_string != None:
		for line in data[start_line_no:]:
			if end_string in line: end_string=line
	else:
		end_string = "\n"
	end_line_no=data[(start_indx):].index(end_string) + 1
	end_indx=(start_line_no+end_line_no-1)
	return(data[(start_indx):(end_indx)])

def crop_header(infile,header_line_count):
	data=get_file_data(infile)
	return(data[header_line_count:])

def column(block,col):
	colval=[]
	for line in block:
		row=line.split()
		colval.append(alnum(row[col]))
	return(colval)

def column_numeric(block,col):
	colval=[]
	for line in block:
		row=line.split()
		colval.append(numeric(row[col]))
	return(colval)

def visblk_print(block):
	for line in block:
		print(line)

