# -*- coding: utf-8 -*-
"""
Spyder Editor

Script to gnerate Linear SFA spectra and mometum distribusions

Andrew Maxwell 11/05/2021
"""

#libraries to use
import numpy as np
import functools
import time
from itertools import repeat
from itertools import product

import multiprocessing
#import ipyparallel as ipp

try:
    cpus = multiprocessing.cpu_count()
except NotImplementedError:
    cpus = 8   # arbitrary default
    
pool = multiprocessing.Pool(processes=cpus)

Pi = np.pi
I = complex(0.,1.)

#Code for file input output (using pickle files)
import pickle
import os
#Name of your DataDump, chose wisely
#Can be the same fie if different DataKey is used
DumpFileName = 'SF_ML_DataDump' 

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


#Import SFA Linear Library
print("Loading SFA Library")
import SFALinearPulse as SFA_Lin
omegaIn = 0.057

print("Setting up pulses")
# Create list of pulse and target parameters
Ns = np.array([6, 6, 6, 6, 6, 6])
CEPs = np.array([0., Pi/2, 0., Pi/2, 0., Pi/2])
Ups = np.array([0.22, 0.22, 0.44, 0.44, 0.66, 0.66])
Ips = np.array([0.579, 0.579, 0.579, 0.579, 0.904, 0.904])
TargetList =np.array([3, 3, 3, 3, 0, 0]) # currently 3=Ar and 0=He, can make more user friendly selction in future
SPs = np.array([SFA_Lin.SFALinearPulse(Ipi, Upi, omegaIn, Ni, CEPi, Targi) for Ipi, Upi, Ni, CEPi, Targi in zip(Ips, Ups, Ns, CEPs, TargetList)])
NPulses = len(SPs)

# Compute Momentum Distribusions¶
dp = 4./250
pzList = np.arange(-1.7, 1.7, dp)
pxList = np.arange(0, 1.7, dp)
pzGrid, pxGrid = np.meshgrid(pzList, pxList)
py = 0.1

print("computing momentum distrbusions for ", NPulses, " pulses, the time of computation will be printed")
#Here we compute the transition amplitudes from scratch and print the computational time required for all pulses in total
t1 = time.time()
MGrids = [ np.array(pool.starmap(SPs[i].Mxz_List, zip(pxGrid, repeat(py),  pzGrid) ))
          for i in range(0, NPulses)]
t2 = time.time()
print(t2 - t1)

MGrids = [np.abs(MGrids[i])**2 for i in range(0, NPulses)]
MMaxs = [np.max(MGrids[i]) for i in range(0, NPulses)]

#Change Key in Data dump to save additional data
DataKeyMomentum = 'MomentumPlot'
print("saving momentum data in file ",DumpFileName,", with key ", DataKeyMomentum)
saveDataDump(DataKeyMomentum, [NPulses, pzList, pxList, MGrids, MMaxs], DumpFileName)

# Spectra Computation
ELim = 2.
dE = ELim/250
EList = np.arange(0, ELim, dE)

print("computing spectra for ", NPulses, " pulses, the time of computation will be printed")
#Here we compute the spectra from scratch  and print the computational time over all pulses in total
t1 = time.time()
SpectraList = [pool.starmap(SPs[i].Spectra, zip(EList, repeat(0.), repeat(np.inf) , repeat(1.0e-02) ))
               for i in range(0,NPulses)]
t2 = time.time()
print(t2 - t1)
DataKeySpectra = 'SpectraPlot'
print("saving spectra data in file ",DumpFileName,", with key ", DataKeySpectra)
saveDataDump(DataKeySpectra, [NPulses, EList, SpectraList], DumpFileName)










