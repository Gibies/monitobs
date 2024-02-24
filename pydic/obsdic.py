#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 20:47:30 2020

@author: gibies
"""
import sys
import os
MONITOBS=os.environ.get('MONITOBS',"../")
OBSLIB=os.environ.get('OBSLIB',MONITOBS+"/pylib")
sys.path.append(OBSLIB)
OBSDIC=os.environ.get('OBSDIC',MONITOBS+"/pydic")
sys.path.append(OBSDIC)
import obslib
import spectradic
import numpy
import datetime

def pydate(year,month,day,hour,minute,second):
    return(datetime.datetime(year,month,day,hour,minute,second))
    
def timeperiod(year=0,mn=0,dy=0,hh=0,mm=0,ss=0):
    days=year*365+mn*30+dy
    seconds=(hh*3600+mm*60+ss)
    return(datetime.timedelta(days=days,seconds=seconds))

timezone = {
        "UTC" : 0,
        "IST" : 82.5
        }

Aeolus = {
        "SatID" : 48,
        "SMAxis" : 6689,
        "SatAlt" : 320,
        "OrbInc" : -96.7,
        "horbit" : 1,
        "morbit" : 30,
        "sorbit" : 42,
        "nodlon" : 90,
        "nodlat" : 0,
        "repeat" : 111,
        "Torbit" : timeperiod(hh=1,mm=30,ss=42),
        "Tref" : pydate(2020,1,20,00,00,00),
        }
Aladin = {
        "SatZenithAngle" : 35,
        "HVA" : 90,
        "VVA" : 35,
        "profsec" : 28,
        "obsec" : 7,
        "ObsAltList" : numpy.concatenate((numpy.arange(0.0,2000.0,500.0),numpy.arange(2000.0,16000.0,1000.0),numpy.arange(16000.0,26000.0,2000.0)),axis=None), ###in meters
        "20700" : [27,3,4,5,6,7,8,1,2,580,584,220,587,579,578,583,588,586,585],
        "22501" : [27 ,3, 4, 5, 6, 7, 257, 258, 259, 260, 261, 1, 2, 187, 93, 32, 52, 53, 192, 55, 194, 197] + [199,200]*4 + [201,202,203]*10,
        }

HLOSWind= {
        "filename" : "HLOSwind.obstore",
        "obsgroup" : 2,
        "subtype" : [20700,],
        "SatID" : 48,
        "SatZenithAngle" : 35,
        "ChannelNumber" : 1,
        "LidarClass" : 0,
        "WindErr" : 1.0,
        "AladinConfFlag" : 0,
        "20700" : [27,3,4,5,6,7,8,1,2,580,584,220,587,579,578,583,588,586,585],
        }

LEOGEOAMV = {
        "filename" : "Satwind.obstore",
        "obsgroup" : 2,
        "subtype" : [22500,],
        "SatID" : 854,
        "OrigCtr" : 160,
        "SatZenithAngle" : 35,
        "CloudMotionMethod" : 1,
        "HeightAssMethod" : 1,
        "AltHeightAss" : 1,
        "ChanCtralFreq" : 28037400000000,
        "Chnldic" : spectradic.ABIdic,
        "22500" : [27 ,3, 4, 5, 6, 7, 257, 258, 259, 260, 261, 1, 2, 187, 93, 32, 52, 53, 192, 55, 194, 197] + [199,200]*4 + [201,202,203]*10 + [281,279,280]*10,
        }

METOPAMV = {
        "filename" : "Satwind.obstore",
        "obsgroup" : 2,
        "subtype" : 22500,
        "SatID" : 852,
        "CloudMotionMethod" : 3,
        "HeightAssMethod" : 2,
        "AltHeightAss" : 2,
        "obs_index" : [27,3, 4, 5, 6, 7, 257, 258, 259, 260, 261, 1, 2, 187, 93, 32, 52, 53, 192, 55, 194, 197] + [199,200]*4 + [201,202,203]*10
        }

INSAT_AMV = {
        "obsgroup" : 2,
        "filename" : "Satwind.obstore",
        "obs_index" : [27,3, 4, 5, 6, 7, 257, 258, 259, 260, 261, 1, 2, 187, 93, 32, 52, 53, 192, 55, 194, 197] + [199,200]*4 + [201,202,203]*10
        }

MODIS_AMV = {
        "obsgroup" : 2,
        "filename" : "Satwind.obstore",
        "obs_index" : [27,3, 4, 5, 6, 7, 257, 258, 259, 260, 261, 1, 2, 187, 93, 32, 52, 53, 192, 55, 194, 197] + [199,200]*4 + [201,202,203]*10 + [281, 279, 280]
        }

ESAHRWVW_AMV = {
        "obsgroup" : 2,
        "filename" : "Satwind.obstore",
        "obs_index" : [27,3, 4, 5, 6, 7, 257, 258, 259, 260, 261, 1, 2, 187, 93, 32, 52, 53, 192, 55, 194, 197] + [199,200]*4 + [201,202,203]*10 + [204]*4 + [215]*4
        }

subtype10100={
        "obs_index" : [3, 4, 5, 6, 7, 1, 2,]
        }

subtype10200={
        "obs_index" : [3, 4, 5, 6, 7, 1, 2,]
        }

subtype10300={
        "obs_index" : [3, 4, 5, 6, 7, 1, 2,]
        }

subtype10400={
        "obs_index" : [3, 4, 5, 6, 7, 1, 2,]
        }

subtype10500={
        "obs_index" : [3, 4, 5, 6, 7, 1, 2,]
        }

Surface={
        "obsgroup" : 1,
        "filename" : "Surface.obstore",
        "subtype" : [10100,10200,10300,10800,11100,40100,],
        "10100" : [59,60,61,105,1,2,36,80,3, 4, 5, 6, 7,46,47,45,103,62,63,64,242,65,44,43,37,],
        "10200" : [30,105,1,2,36,80,3,4,5, 6, 7,46,47,45,103,62,63,64,242,65,34,67,68,57,],
        "10300" : [104,3, 4, 5, 6, 7, 1, 2,120,121,46,47,45,103,62,44,65,100,101,],
        "10400" : [3, 4, 5, 6, 7, 1, 2,],
        "10500" : [3, 4, 5, 6, 7, 1, 2,],
        "10600" : [3, 4, 5, 6, 7, 1, 2,],
        "10700" : [3, 4, 5, 6, 7, 1, 2,],
        "10800" : [30,105,1,2,36,80,3, 4, 5, 6, 7,46,47,45,103,62,63,64,65,44,43,37,57,],
        "10900" : [3, 4, 5, 6, 7, 1, 2,],
        "11000" : [3, 4, 5, 6, 7, 1, 2,],
        "11100" : [30,105,1,2,36,3, 4, 5, 6, 7, 46,47,45,103,63,64,65,57,],
        "40100" : [3, 4, 5, 6, 7, 1, 2,65,46,47,437,438,],
        }

Sonde={
       "obsgroup" : 5,
       "filename" : "Sonde.obstore",
       "subtype" : [50100,50200,50300,50400,],
       "50100" : [92,59,60,30,1,2,36,80,3, 4, 5, 6, 7, 81, 82,83,84,85,88,89,90,32,91,55,56,52,53,57,],
       "50200" : [92,59,60,30,1,2,36,80,3, 4, 5, 6, 7, 81, 82,83,84,85,88,89,90,32,91,55,56,52,53,57,],
       "50300" : [92,59,60,30,1,2,36,80,3, 4, 5, 6, 7, 81, 82,83,84,85,88,89,90,32,91,55,56,52,53,57,],
       "50400" : [59,60,105,3, 4, 5, 6, 7, 1, 2,36,88,91,159,160,208,207,],
       "50500" : [3, 4, 5, 6, 7, 1, 2,],
       "50600" : [3, 4, 5, 6, 7, 1, 2,],
       "50700" : [3, 4, 5, 6, 7, 1, 2,],
       "50800" : [3, 4, 5, 6, 7, 1, 2,],
       "50900" : [3, 4, 5, 6, 7, 1, 2,],
       "51000" : [3, 4, 5, 6, 7, 1, 2,],
       "51100" : [3, 4, 5, 6, 7, 1, 2,],
       "51200" : [3, 4, 5, 6, 7, 1, 2,],
       }

Aircraft={
        "obsgroup" : 3,
        "filename" : "Aircraft.obstore",
        "subtype" : [30100,30200,],
        "30100" : [40,30,42,3, 4, 5, 6, 7, 1, 2,41,58,32,52,53,123,55,262,57,],
        "30200" : [30,96,3, 4, 5, 6, 7, 1, 2,58,52,53,123,55,262,97,29,57,],
        "30300" : [3, 4, 5, 6, 7, 1, 2,],
        "30400" : [3, 4, 5, 6, 7, 1, 2,],
        "30500" : [3, 4, 5, 6, 7, 1, 2,],
        "30600" : [3, 4, 5, 6, 7, 1, 2,],
        "30700" : [3, 4, 5, 6, 7, 1, 2,],
        "30800" : [3, 4, 5, 6, 7, 1, 2,],
        "30900" : [3, 4, 5, 6, 7, 1, 2,],
        "31000" : [3, 4, 5, 6, 7, 1, 2,],
        "31100" : [3, 4, 5, 6, 7, 1, 2,],
        }

ahiclr={
	"filename" : "AHIClr.obstore",
	}

aod =	{
	"filename" : "AOD.obstore",
	}

abiclr = {
	"filename" : "ABIClr.obstore",
	}

saphir = {
	"filename" : "SAPHIR.obstore",
	}

fy3mws = {
	"filename" : "MWSFY3.obstore",
	}

mwsfy3 = {
	"filename" : "MWSFY3.obstore",
	}

mwsfy3b = {
	"filename" : "MWSFY3B.obstore",
	}

mwsfy3c = {
	"filename" : "MWSFY3C.obstore",
	}

mwri = {
	"filename" : "MWRI.obstore",
	}

obstype={
        "surface" 	: Surface,
        "sonde"		: Sonde,
        "aircraft"	: Aircraft,
        "leogeo" 	: LEOGEOAMV,           
        "hlosw" 	: HLOSWind,
	"goesabi" 	: { "filename" : "ABIClr.obstore",},
	"fy3mws" 	: { "filename" : "MWSFY3.obstore",},
	"satwind" 	: { "filename" : "Satwind.obstore",},
	"scatwind" 	: { "filename" : "Scatwind.obstore",},
	"grndgps" 	: { "filename" : "GroundGPS.obstore",},
	"seviri"  	: { "filename" : "SEVIRIClr.obstore", },
	"mwri"		: {"filename" : "MWRI.obstore", },
	"sattcwv" 	: {"filename" : "SatTCWV.obstore", },
        }

obstypelist=["aladin",]  #"aircraft","sonde","surface",

station_call_list=[
        "AcarsStatus",
        "Altitude",
        "BuoyID",
        "CallSign",
        "CharData",
        "CrossTrackCell",
        "FlightPhase",
        "Latitude",
        "Longitude",
        "ObPractice",
        "PESR_SNSR_HGHT",
        "RADI_SNDE_TYPE",
        "RADTN_CORTN",
        "RPRT_IDNY",
        "SatDirectionOfMotion",
        "SatID",
        "SatInst",
        "SatProc",
        "StationReportType",
        "TailNumber",
        "TRCKG_SYTM",
        "WMOBlockNo", #59
        "WMORegNo", #60
        "WMOStnNo", #61
        "Zstation", #37
        ]
