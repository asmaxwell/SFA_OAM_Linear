#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 11 11:37:03 2021

@author: Andrew Maxwell

Script to Plot SFA Linear Data
"""

#libraries to use
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size

mpl.rcParams['savefig.pad_inches'] = 0


#Code for file input output (using pickle files)
import pickle
import os
def loadDataDump(DataKey, DumpFileName):
    '''Function to retrieve data from pickle dump'''
    DataDumpFile = open(DumpFileName,'rb')
    DataOut = pickle.load(DataDumpFile)[DataKey]
    DataDumpFile.close()
    return DataOut

#Filenames and paths
DumpFileName = 'SF_ML_DataDump' 
ExportPath = "/home/amaxwell/Dropbox/Documents/2020_PostDoc_ICFO/SF_ML/Documents/Figures/" #chagne this!


[NPulses, pzList, pxList, MGrids, MMaxs] = loadDataDump('MomentumPlot', DumpFileName)

plt.style.use('default')

print("Making Momentum figures")
plt.figure(num=None, figsize=(14, 5*NPulses), dpi=80, facecolor='w', edgecolor='k')
plt.autoscale(tight=True)
for i in range(0, NPulses):
    ax = plt.subplot(NPulses, 1, i+1)
#     ax.text(0.025, 0.95, panelLabels[2][i], transform=ax.transAxes,
#       fontsize=40, fontweight='bold', va='top', color ='white')
    plt.ylabel('$p_{x}$ a.u.', fontsize=(24))
    plt.xlabel('$p_{z}$ a.u.', fontsize=(24))
    plt.xticks(fontsize=(19))
    plt.yticks(fontsize=(19))
    
    #ax = plt.gca()
    im = ax.imshow(np.flip(MGrids[i],0), extent = (np.amin(pzList), np.amax(pzList), np.amin(pxList), np.amax(pxList)),
              cmap=cm.inferno, norm=LogNorm(vmin=MMaxs[i]*1e-6, vmax=MMaxs[i]), 
                interpolation = 'bicubic', aspect = 1.)

    
    aspect = 15
    pad_fraction = 0.5
    divider = make_axes_locatable(ax)
    width = axes_size.AxesY(ax, aspect=1./aspect)
    pad = axes_size.Fraction(pad_fraction, width)
    cax = divider.append_axes("right", size=width, pad=pad)
    cbar = plt.colorbar(im, cax=cax)
    cbar.ax.tick_params(labelsize=20)

print("Exporting Momentum figures")
plt.savefig(ExportPath + "MomPlots.pdf")


print("Making spectra figure")
#spectra plots
[NPulses, EList, SpectraList] = loadDataDump('SpectraPlot', DumpFileName)
plt.figure(num=None, figsize=(20, 12), dpi=80, facecolor='w', edgecolor='k')
plt.ylabel('Signal arb.', fontsize=20)
plt.xlabel('E a.u.', fontsize=20)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.yscale("log")
plt.xlim(0.01, 1.75)
plt.ylim(1e-6, 1e2)
for Spec in SpectraList:
    plt.plot(EList,Spec,'-')
    
print("Exporintg spectra figure")
plt.savefig(ExportPath + "SpectraPlots.pdf")
