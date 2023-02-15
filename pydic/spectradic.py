#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 20:40:01 2020

@author: gibies
"""



ABIChanlist=["Ch15","Ch16","Ch10","Ch9","Ch8","Ch12"]

ABItoplev={
        "Ch8" : 150,
        "Ch9" : 400,
        "Ch10" : 550,
        "Ch12" : 0,
        "Ch15" : 900,
        "Ch16" : 750,
        }

ABIbotlev={
        "Ch8" : 400,
        "Ch9" : 550,
        "Ch10" : 750,
        "Ch12" : 150,
        "Ch15" : 1200,
        "Ch16" : 900,
        }

ABIspecband={
        "Ch1" : 0.47,
        "Ch2" : 0.64,
        "Ch3" : 0.86,
        "Ch4" : 1.38,
        "Ch5" : 1.61,
        "Ch6" : 2.26,
        "Ch7" : 3.90,
        "Ch8" : 6.19,
        "Ch9" : 6.95,
        "Ch10" : 7.34,
        "Ch11" : 8.5,
        "Ch12" : 9.61,
        "Ch13" : 10.35,
        "Ch14" : 11.2,
        "Ch15" : 12.3,
        "Ch16" : 13.3,
        }

ABIdic={
        "Chanlist" : ABIChanlist,
        "Chantoplev" : ABItoplev,
        "Chanbotlev" : ABIbotlev,
        "Chanwavlen" : ABIspecband,
        "Chwlunitfctr" : 1000000,
        }