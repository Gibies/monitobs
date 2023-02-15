#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 09:49:00 2020

@author: bushair
Adopted to Monitobs on Tue Nov 10 10:05:00 2020
"""
import numpy as np
import matplotlib.pyplot as plt
import math

def rmse(predictions, targets):
    differences = predictions - targets                       
    differences_squared = differences ** 2                    
    mean_of_differences_squared = differences_squared.mean()  
    rmse_val = np.sqrt(mean_of_differences_squared)           
    return rmse_val 
def estimate_coef(x, y): 
    n = np.size(x) 
    m_x, m_y = np.mean(x), np.mean(y) 
    SS_xy = np.sum(y*x) - n*m_y*m_x 
    SS_xx = np.sum(x*x) - n*m_x*m_x 
    b_1 = SS_xy / SS_xx 
    b_0 = m_y - b_1*m_x 
    return(b_0, b_1) 
def plot_regression_line(x, y, b):
    y_pred = b[0] + b[1]*x
    bn=20 
    bins = [bn, bn]
    hh, locx, locy = np.histogram2d(x, y, bins=bins)
    z = np.array([hh[np.argmax(a<=locx[1:]),np.argmax(b<=locy[1:])] for a,b in zip(x,y)])
    idx = z.argsort()
    x2, y2, z2 = x[idx], y[idx], z[idx]
    plt.scatter(x2, y2, c=z2, cmap='jet', marker='.') 
    plt.colorbar(extend='max') 
    plt.plot(x, y_pred, color = "g") 
