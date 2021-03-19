#!/bin/env python
#
# PyROOT study of standard selector performance using sim-hit matching 
# to identify fake and signal muons
#
import os, re, ROOT, sys, pickle, time
from pprint import pprint
from math import *
from array import array
from DataFormats.FWLite import Events, Handle
import numpy as np


##
## User Input
##
def getPFNs(lfns):
    files = []
    for file in lfns:
        fullpath = "/eos/cms/" + file
        if os.path.exists(fullpath):
            files.append(fullpath)
        else:
            raise Exception("File not found: %s" % fullpath)
    return files


def IsMatched(muon1,muon2,sharedFrac=0.5):
    # first check if muon share Wh/Se/St 
    if (muon1.whNum()!=muon2.whNum()): return False 
    if (muon1.scNum()!=muon2.scNum()): return False 
    if (muon1.stNum()!=muon2.stNum()): return False     
    
    # now count the number of shared hits: 
    numShared=0.
    totMuon1=0.
    for ly in range(0,7):
        if (muon1.pathWireId(ly)>=0): 
            totMuon1=totMuon1+1. 
        else:  
            continue
        if (muon1.pathWireId(ly)!=muon2.pathWireId(ly)): continue
        if (muon1.pathTDC(ly)!=muon2.pathTDC(ly)): continue
        
        numShared = numShared+1.

    if (numShared/totMuon1 >= sharedFrac): return True

    return False
# files = getPFNs(lfns)

muoBayesHandle, muoBayesLabel = Handle("L1Phase2MuDTExtPhContainer"), ("dtTriggerPhase2BayesPrimitiveDigis","","L1DTTrigPhase2Prod")
muoStdHandle, muoStdLabel = Handle("L1Phase2MuDTExtPhContainer"), ("dtTriggerPhase2StdPrimitiveDigis","","L1DTTrigPhase2Prod")
genHandle, genLabel = Handle("vector<reco::GenParticle>"), "genParticles"


ROOT.gROOT.SetBatch(True)

##
## Main part
##
files = ['/afs/cern.ch/user/n/ntrevisa/work/DT/CMSSW_11_2_0_pre2/src/DTTriggerPhase2Primitives.root']

print "Number of files: %d" % len(files)
events = Events(files)    
print "we got the events" 
## load some histograms (and efficiencies): 
outputDict = {} 
dumpToFile = True

for frac in [0.25,0.5,0.75,1.00]:
    print "shared fraction is %0.2f" %(frac)
    fracname="shared%i" %(frac*100)

    hPhiRes_q1 = []
    hPhiRes_q3 = []
    hPhiRes_q5 = []
    hPhiRes_q8 = []
    
    hPhiBRes_q1 = []
    hPhiBRes_q3 = []
    hPhiBRes_q5 = []
    hPhiBRes_q8 = []
    
    hChi2Res_q1 = []
    hChi2Res_q3 = []
    hChi2Res_q5 = []
    hChi2Res_q8 = []
    
    hBxRes_q1 = []
    hBxRes_q3 = []
    hBxRes_q5 = []
    hBxRes_q8 = []
    
    hTimeRes_q1 = []
    hTimeRes_q3 = []
    hTimeRes_q5 = []
    hTimeRes_q8 = []

    for st in range(1,5):    
        hPhiRes_q1.append(ROOT.TH1F("hPhiRes_MB%i_q1_%s" %(st,fracname),"",20,-5000,5000.))
        hPhiRes_q3.append(ROOT.TH1F("hPhiRes_MB%i_q3_%s" %(st,fracname),"",20,-5000,5000.))
        hPhiRes_q5.append(ROOT.TH1F("hPhiRes_MB%i_q5_%s" %(st,fracname),"",20,-5000,5000.))
        hPhiRes_q8.append(ROOT.TH1F("hPhiRes_MB%i_q8_%s" %(st,fracname),"",20,-5000,5000.))
        
        hPhiBRes_q1.append(ROOT.TH1F("hPhiBRes_MB%i_q1_%s" %(st,fracname),"",20,-100,100.))
        hPhiBRes_q3.append(ROOT.TH1F("hPhiBRes_MB%i_q3_%s" %(st,fracname),"",20,-100,100.))
        hPhiBRes_q5.append(ROOT.TH1F("hPhiBRes_MB%i_q5_%s" %(st,fracname),"",20,-100,100.))
        hPhiBRes_q8.append(ROOT.TH1F("hPhiBRes_MB%i_q8_%s" %(st,fracname),"",20,-100,100.))
        
        hChi2Res_q1.append(ROOT.TH1F("hChi2Res_MB%i_q1_%s" %(st,fracname),"",20,-5000,5000.))
        hChi2Res_q3.append(ROOT.TH1F("hChi2Res_MB%i_q3_%s" %(st,fracname),"",20,-5000,5000.))
        hChi2Res_q5.append(ROOT.TH1F("hChi2Res_MB%i_q5_%s" %(st,fracname),"",20,-5000,5000.))
        hChi2Res_q8.append(ROOT.TH1F("hChi2Res_MB%i_q8_%s" %(st,fracname),"",20,-5000,5000.))
        
        hBxRes_q1.append(ROOT.TH1F("hBxRes_MB%i_q1_%s" %(st,fracname),"",20,-10,10.))
        hBxRes_q3.append(ROOT.TH1F("hBxRes_MB%i_q3_%s" %(st,fracname),"",20,-10,10.))
        hBxRes_q5.append(ROOT.TH1F("hBxRes_MB%i_q5_%s" %(st,fracname),"",20,-10,10.))
        hBxRes_q8.append(ROOT.TH1F("hBxRes_MB%i_q8_%s" %(st,fracname),"",20,-10,10.))
        
        hTimeRes_q1.append(ROOT.TH1F("hTimeRes_MB%i_q1_%s" %(st,fracname),"",20,-100,100.))
        hTimeRes_q3.append(ROOT.TH1F("hTimeRes_MB%i_q3_%s" %(st,fracname),"",20,-100,100.))
        hTimeRes_q5.append(ROOT.TH1F("hTimeRes_MB%i_q5_%s" %(st,fracname),"",20,-100,100.))
        hTimeRes_q8.append(ROOT.TH1F("hTimeRes_MB%i_q8_%s" %(st,fracname),"",20,-100,100.))
    
    outputDict[fracname] = {}
    for st in range(1,5):
        outputDict[fracname]['hPhiRes_MB%i_q1'%st] = hPhiRes_q1[st-1]
        outputDict[fracname]['hPhiRes_MB%i_q3'%st] = hPhiRes_q3[st-1]
        outputDict[fracname]['hPhiRes_MB%i_q5'%st] = hPhiRes_q5[st-1]
        outputDict[fracname]['hPhiRes_MB%i_q8'%st] = hPhiRes_q8[st-1]       

        outputDict[fracname]['hPhiBRes_MB%i_q1'%st] = hPhiBRes_q1[st-1]
        outputDict[fracname]['hPhiBRes_MB%i_q3'%st] = hPhiBRes_q3[st-1]
        outputDict[fracname]['hPhiBRes_MB%i_q5'%st] = hPhiBRes_q5[st-1]
        outputDict[fracname]['hPhiBRes_MB%i_q8'%st] = hPhiBRes_q8[st-1]
    
        outputDict[fracname]['hChi2Res_MB%i_q1'%st] = hChi2Res_q1[st-1]
        outputDict[fracname]['hChi2Res_MB%i_q3'%st] = hChi2Res_q3[st-1]
        outputDict[fracname]['hChi2Res_MB%i_q5'%st] = hChi2Res_q5[st-1]
        outputDict[fracname]['hChi2Res_MB%i_q8'%st] = hChi2Res_q8[st-1]       

        outputDict[fracname]['hBxRes_MB%i_q1'%st] = hBxRes_q1[st-1]
        outputDict[fracname]['hBxRes_MB%i_q3'%st] = hBxRes_q3[st-1]
        outputDict[fracname]['hBxRes_MB%i_q5'%st] = hBxRes_q5[st-1]
        outputDict[fracname]['hBxRes_MB%i_q8'%st] = hBxRes_q8[st-1]       

        outputDict[fracname]['hTimeRes_MB%i_q1'%st] = hTimeRes_q1[st-1]
        outputDict[fracname]['hTimeRes_MB%i_q3'%st] = hTimeRes_q3[st-1]
        outputDict[fracname]['hTimeRes_MB%i_q5'%st] = hTimeRes_q5[st-1]
        outputDict[fracname]['hTimeRes_MB%i_q8'%st] = hTimeRes_q8[st-1]       

    # loop over events
    print 'here :)'
    count = 0
    
    if (dumpToFile): 
        f= open("EventDumpList_BayesToStd.log","w+")
    
    for ev in events:
        if not count%10:  print count, events.size()
        count = count+1
        ev.getByLabel(muoBayesLabel, muoBayesHandle)
        ev.getByLabel(muoStdLabel, muoStdHandle)
        ev.getByLabel(genLabel, genHandle)

        muon_bayes = muoBayesHandle.product().getContainer()
        muon_std = muoStdHandle.product().getContainer()
    
        if (dumpToFile):
            f.write( "\nInspecting Event Number %i \n" %(ev.eventAuxiliary().id().event())  )
            f.write( "         Wh   Se   St  | w1 w2 w3 w4 w5 w6 w7 w8 |  tdc1  tdc2  tdc3  tdc4  tdc5  tdc6  tdc7  tdc8 | Q     phi  phib   chi2  bX\n" )
            f.write( "--------------------------------------------------------------------------------------------------------------------------------\n" )
            for muon in muon_bayes:
                f.write( "[Bayes]: Wh%2d Se%2d St%1d | %2d %2d %2d %2d %2d %2d %2d %2d | %5d %5d %5d %5d %5d %5d %5d %5d | Q%1d %6d %5d %7d %2d\n" %(muon.whNum(),muon.scNum(),muon.stNum(),muon.pathWireId(0),muon.pathWireId(1),muon.pathWireId(2),muon.pathWireId(3),muon.pathWireId(4),muon.pathWireId(5),muon.pathWireId(6),muon.pathWireId(7),muon.pathTDC(0),muon.pathTDC(1),muon.pathTDC(2),muon.pathTDC(3),muon.pathTDC(4),muon.pathTDC(5),muon.pathTDC(6),muon.pathTDC(7),muon.quality(),muon.phi(),muon.phiBend(),muon.chi2(),muon.bxNum()-20) )
        
        
        for muon in muon_std:
            if (dumpToFile): 
                f.write( "[Std  ]: Wh%2d Se%2d St%1d | %2d %2d %2d %2d %2d %2d %2d %2d | %5d %5d %5d %5d %5d %5d %5d %5d | Q%1d %6d %5d %7d %2d\n" %(muon.whNum(),muon.scNum(),muon.stNum(),muon.pathWireId(0),muon.pathWireId(1),muon.pathWireId(2),muon.pathWireId(3),muon.pathWireId(4),muon.pathWireId(5),muon.pathWireId(6),muon.pathWireId(7),muon.pathTDC(0),muon.pathTDC(1),muon.pathTDC(2),muon.pathTDC(3),muon.pathTDC(4),muon.pathTDC(5),muon.pathTDC(6),muon.pathTDC(7),muon.quality(),muon.phi(),muon.phiBend(),muon.chi2(), muon.bxNum()-20)  )
        
        ## now match with the previous 
        for muon1 in muon_bayes:
            st = muon1.stNum()-1
            for muon2 in muon_std: 
                if not IsMatched(muon1,muon2,frac): continue
                if (muon1.quality()>=1) :  
                    hPhiRes_q1[st]  .Fill( (muon1.phi()-muon2.phi()) )
                    hPhiBRes_q1[st] .Fill( (muon1.phiBend()-muon2.phiBend()) )
                    hChi2Res_q1[st] .Fill( (muon1.chi2()-muon2.chi2()) )
                    hBxRes_q1[st]   .Fill( (muon1.bxNum()-muon2.bxNum()) )
                    hTimeRes_q1[st] .Fill( (muon1.t0()-muon2.t0()) )
                    
                if (muon1.quality()>=3) :
                    hPhiRes_q3[st]  .Fill( (muon1.phi()-muon2.phi()) )
                    hPhiBRes_q3[st] .Fill( (muon1.phiBend()-muon2.phiBend()) )
                    hChi2Res_q3[st] .Fill( (muon1.chi2()-muon2.chi2()) )
                    hBxRes_q3[st]   .Fill( (muon1.bxNum()-muon2.bxNum()) )
                    hTimeRes_q3[st] .Fill( (muon1.t0()-muon2.t0()) )

                if (muon1.quality()>=5) :
                    hPhiRes_q5[st]  .Fill( (muon1.phi()-muon2.phi()) ) 
                    hPhiBRes_q5[st] .Fill( (muon1.phiBend()-muon2.phiBend()) )
                    hChi2Res_q5[st] .Fill( (muon1.chi2()-muon2.chi2()) )
                    hBxRes_q5[st]   .Fill( (muon1.bxNum()-muon2.bxNum()) )
                    hTimeRes_q5[st] .Fill( (muon1.t0()-muon2.t0()) )
               
                if (muon1.quality()>=8) :
                    hPhiRes_q8[st]  .Fill( (muon1.phi()-muon2.phi()) )
                    hPhiBRes_q8[st] .Fill( (muon1.phiBend()-muon2.phiBend()) )
                    hChi2Res_q8[st] .Fill( (muon1.chi2()-muon2.chi2()) )
                    hBxRes_q8[st]   .Fill( (muon1.bxNum()-muon2.bxNum()) )
                    hTimeRes_q8[st] .Fill( (muon1.t0()-muon2.t0()) )                

    if (dumpToFile): f.close()
    dumpToFile=False
    ev.toBegin()


import pickle 
with open('GroupingComparison_BayesToStd.pickle', 'wb') as handle:
    pickle.dump(outputDict, handle, protocol=pickle.HIGHEST_PROTOCOL)
