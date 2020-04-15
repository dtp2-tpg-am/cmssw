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

# files = getPFNs(lfns)


muonHandle, muonLabel = Handle("std::vector<pat::Muon>"), "slimmedMuons"
muoBayesHandle, muoBayesLabel = Handle(),""
muoStdHandle, muonStdLabel = Handle(),""
genHandle, genLabel = Handle("vector<reco::GenParticle>"), "genParticles"

minPt = 20
maxPt = 1e9

n_events_limit = None
# n_events_limit = 10000

ROOT.gROOT.SetBatch(True)

##
## Main part
##


files = ['/afs/cern.ch/work/s/sesanche/private/Muonico/CMSSW_10_1_4/src/MuonAnalysis/MuonFastFeedbackTools/data/08540076-DB76-E811-A01E-FA163EF8C86B.root']

print "Number of files: %d" % len(files)

events = Events(files)



# loop over events
for event in events:
    event.getByLabel(muonLabel, muonHandle)
    muons = muonHandle.product()
    for muon in muons:
        print muon.miniPFIsolation().chargedHadronIso()
