#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 29 13:56:23 2022
Script adapted from code written by Xavier Barcons Planas

Generation script for data from:
    Ultrafast imaging of molecular chirality with photoelectron vortices
Xavier Barcons Planas, Andrés Ordóñez, Maciej Lewenstein, Andrew Stephen Maxwell

Preprint available at: https://arxiv.org/abs/2202.07289
--------------------------------------------------------------------------------

You must install SFA_OAM_Linear before running this script
Find instruction available on the github repository https://github.com/asmaxwell/SFA_OAM_Linear

Please note this script will loop twice to generate all the data, once with 'Gauge=0' for the length gauge and once with 
change dp to reduce computation by having a smaller grid
"""

import numpy as np
from itertools import repeat
import multiprocessing

try:
    cpus = multiprocessing.cpu_count()
except NotImplementedError:
    cpus = 8   # arbitrary default
    
pool = multiprocessing.Pool(processes=cpus)

Pi = np.pi
I = complex(0.,1.)

### Dump File and Pickling
#currently this is how Data may be outputted, as a pickle file, we may want to change this to JSON? However, pickle files are very compact.
import pickle
import os
DumpFileName = 'OAM_Lin_DataDump'

# Make momentum grid
dp=0.01
pzList = np.arange(-1.7, 1.7, dp)
pxList = np.arange(0, 1., dp)
pzGrid, pxGrid = np.meshgrid(pzList, pxList)

def saveDataDump(DataKey, Data, DumpFileName):
    '''Function saves data to file via pickling for easy access and loading'''
    Dics={"Null" : 0}
    #load full database
    if(os.path.exists(DumpFileName)):
        with open(DumpFileName,'rb') as rfp:
            Dics.update(pickle.load(rfp))
            
    #Make new dictionary of data
    OutputDataDict={DataKey : Data}
    #Append new dicitonary to all dictionarys
    Dics.update(OutputDataDict)
    #open file and dump all dicitionarys
    DataDumpFile = open(DumpFileName,'wb')
    pickle.dump(Dics, DataDumpFile)
    DataDumpFile.close()
    
def loadDataDump(DataKey, DumpFileName):
    '''Function to retrieve data from pickle dump'''
    DataDumpFile = open(DumpFileName,'rb')
    DataOut = pickle.load(DataDumpFile)[DataKey]
    DataDumpFile.close()
    return DataOut

## Import SFA Linear Pulse
import SFALinearPulse as SFA_Lin
# Generate Figure paper data


#loop to do both gauges
for Gauge in [0,1]: #Velocity gauge: Gauge=0 \\ Length gauge: Gauge=1
    #Make class instance
    omegaIn = 0.057
    IpIn = 0.579
    UpIn = 0.44
    NIn = 2
    
    CEPIn = [0, Pi/2, Pi, 3*Pi/2]
    ts = Pi/omegaIn + I* Pi/(2*omegaIn)
    
    #Velocity gauge: Gauge=0 \\ Length gauge: Gauge=1
    #Gauge=0 #Gauge set in loop now
    
    enant=0
    nIn=[4,4,4,4]
    lIn=[2,3,2,3]
    mIn=[1,1,-1,-1]
    constIn=[1,1j,-1,1j]
    
    #make SFA pulse instance
    SP = np.array([SFA_Lin.SFALinearPulse(IpIn, UpIn, omegaIn, NIn, nIn, lIn, mIn, constIn, CEPi, Gauge) 
                   for CEPi in CEPIn])
    NPulses = len(SP)
    
    # OAM + OAM residual difference distributions data
    OAMList=[-1,1]
    MGrids_align_CEP = np.array([[pool.starmap(SP[i].aligned_av_list, zip(pxGrid, pzGrid, repeat(OAM), repeat(enant)) ) 
                         for OAM in OAMList] for i in range(0, NPulses)])
    OAMList=[-2,2,1]
    MGrids_oriAv_CEP = np.array([[pool.starmap(SP[i].orientation_av_list, zip(pxGrid, pzGrid, repeat(OAM), repeat(enant)) ) 
                         for OAM in OAMList] for i in range(0, NPulses)])
    
    #CEP averaging
    MGrids_align = np.sum(MGrids_align_CEP, 0)/NPulses
    MGrids_oriAv = np.sum(MGrids_oriAv_CEP, 0)/NPulses
    
    #insets data
    M_align_diff = 2*(-MGrids_align[1]+MGrids_align[0])/(MGrids_align[1]+MGrids_align[0]+1e-32)
    M_oriAv_diff = 2*(-MGrids_oriAv[1]+MGrids_oriAv[0])/(MGrids_oriAv[1]+MGrids_oriAv[0]+1e-32)
    
    Fig2InsetData = [[MGrids_align[1], MGrids_oriAv[2]], [M_align_diff, M_oriAv_diff]]
    
    #May be used to compute OAM spectra but this part has been switched off as it was not used in the final publication
    # upper hemisphere integrals
    # ELim = 0.5 if Gauge==0 else 1.
    # dE = 0.002
    # EList = np.arange(0, ELim, dE)
    
    # upper hemisphere integrals
    
    # OAMList = [-2,-1,0,1,2]
    # Mint_align_CEP = np.array([[pool.starmap(SP[i].M_integration_align, zip(EList, repeat(0.), repeat(Pi/6), repeat(OAM), repeat(1.0e-02), repeat(500), repeat(enant) )) 
    #                        for OAM in OAMList] for i in range(0, NPulses)])
    # OAMList = [-2,-1,0,1,2]
    # Mint_oriAv_CEP = np.array([[pool.starmap(SP[i].M_integration_av, zip(EList, repeat(0.), repeat(Pi/6), repeat(OAM), repeat(1.0e-02), repeat(500), repeat(enant) )) 
    #                        for OAM in OAMList] for i in range(0, NPulses)])
    
    #CEP averaging
    # Mint_align = np.sum(Mint_align_CEP, 0)/NPulses
    # Mint_oriAv = np.sum(Mint_oriAv_CEP, 0)/NPulses
    
    
    #subplots data
    # Mint_align_diff1 = 2*(-Mint_align[4]+Mint_align[0])/(Mint_align[4]+Mint_align[0]+1e-32)
    # Mint_align_diff2 = 2*(-Mint_align[3]+Mint_align[1])/(Mint_align[3]+Mint_align[1]+1e-32)
    # Mint_align_diff3 = -Mint_align_diff1
    # Mint_align_diff4 = -Mint_align_diff2
    # Mint_align_diff =[Mint_align_diff1, Mint_align_diff2, Mint_align_diff3, Mint_align_diff4]
    
    # Mint_oriAv_diff1 = 2*(-Mint_oriAv[4]+Mint_oriAv[0])/(Mint_oriAv[4]+Mint_oriAv[0]+1e-32)
    # Mint_oriAv_diff2 = 2*(-Mint_oriAv[3]+Mint_oriAv[1])/(Mint_oriAv[3]+Mint_oriAv[1]+1e-32)
    # Mint_oriAv_diff3 = -Mint_oriAv_diff1
    # Mint_oriAv_diff4 = -Mint_oriAv_diff2
    # Mint_oriAv_diff =[Mint_oriAv_diff1, Mint_oriAv_diff2, Mint_oriAv_diff3, Mint_oriAv_diff4]
    
    #Mint_oriAv[3]=np.zeros(len(Mint_oriAv[0]))
    
    # Fig2Data = [[Mint_align, Mint_oriAv], [Mint_align_diff, Mint_oriAv_diff]]
    Fig2Data = []
    
    #save data:
    DataName = 'Fig_Paper_Velocity' if Gauge==0 else 'Fig_Paper_Length'
    saveDataDump(DataName, [Fig2InsetData, Fig2Data], DumpFileName)
    #saveDataDump('Fig_Paper_Length', [Fig2InsetData, Fig2Data], DumpFileName)
