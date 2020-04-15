#!/bin/env python
#
#

import os, re, ROOT, sys, pickle, time
from pprint import pprint
from math import *
from array import array
from DataFormats.FWLite import Events, Handle
import numpy as np

datasets = { 
    #'SingleMu_NOPU'  : { 'path' : '/eos/cms/store/group/phys_muon/sesanche/overlap_oct17_2/Mu_FlatPt2to100-pythia8-gun/SingleMu_NOPU/191017_094114/0000/',},
    'SingleMu_PU200_P2' : { 'path' : '/eos/cms/store/user/folguera/OMTF/Feb20_P2TP_Mar8/Mu_FlatPt2to100-pythia8-gun/crab_L1MuPhase2Ntuples_Mu_FlatPt2to100_PU200/200308_212910/0000/',},
    'SingleMu_PU200_P1' : { 'path' : '/eos/cms/store/user/folguera/OMTF/Feb20_P1TP/Mu_FlatPt2to100-pythia8-gun/crab_L1MuPhase2Ntuples_Mu_FlatPt2to100_PU200/200204_162831/0000/',},
    #'SingleMu_NOPU_aging1000'  : { 'path' : '/eos/cms/store/group/phys_muon/cericeci/overlap_oct22/Mu_FlatPt2to100-pythia8-gun/SingleMu_NOPU/191031_163317/0000/',},
    #'SingleMu_PU200_aging1000' : { 'path' : '/eos/cms/store/group/phys_muon/cericeci/overlap_oct22/Mu_FlatPt2to100-pythia8-gun/SingleMu_PU200/191031_155204/0000/',},
    #'SingleMu_NOPU_aging3000'  : { 'path' : '/eos/cms/store/group/phys_muon/cericeci/overlap_oct22/Mu_FlatPt2to100-pythia8-gun/SingleMu_NOPU/191031_155658/0000/',},
    #'SingleMu_PU200_aging3000' : { 'path' : '/eos/cms/store/group/phys_muon/cericeci/overlap_oct22/Mu_FlatPt2to100-pythia8-gun/SingleMu_PU200/191031_155834/0000/',},


}

def getEta(etaHW):
    return  etaHW/240.*2.61

def getPt(ptHW):
    return (ptHW-1.)/2.

def modulo(val):
    while val > 2*pi: val -= 2*pi
    while val <0: val += 2*pi
    return val

def getPhi(obj):
    return modulo(((15.+obj.processor()*60.)/360. + obj.hwPhi()/576.)*2*pi)

def getP4(obj):
    pt = obj.pt()
    eta= obj.eta()
    phi= obj.phi()
    v = ROOT.TLorentzVector()
    v.SetPtEtaPhiM( pt, eta, phi, 0)
    return v

def getP4FromHW(obj):
    pt = obj.hwPt()
    eta= obj.hwEta()
    phi= obj.hwPhi()
    v = ROOT.TLorentzVector()
    v.SetPtEtaPhiM( getPt(pt), getEta(eta), getPhi(obj), 0)
    return v






for dataset in datasets:
    thefiles = [] 
    for files in os.listdir(datasets[dataset]['path']):
        #fil = filter( lambda x : '.root' in x and "super" in x, files)
        if not(".root" in files and "l1tomtf" in files): continue
        thefiles.append( datasets[dataset]['path'] + '/' + files)
    datasets[dataset]['files'] = thefiles
        
muonHandle, muonLabel = Handle("BXVector<l1t::RegionalMuonCand>"), ("simBayesOmtfDigis", "OMTF", "L1TMuonEmulation" )
genHandle, genLabel = Handle("vector<reco::GenParticle>"), "genParticles"

outputDict = {} 
maxEvents = -1

for dataset in datasets:
    print 'starting to process', dataset
    events = Events(datasets[dataset]['files'])
#    events = Events("/afs/cern.ch/user/f/folguera/workdir/Upgrade/L1T/CMSSW_10_6_1_patch2_omtfFeb20_P1TP/src/L1Trigger/L1TMuonBayes/test/crab/l1tomtf_superprimitives1.root")
    print 'we got the events'
    hEta20_q12_cut20   = ROOT.TEfficiency("hEta20_%s_q12"%dataset,"",240,-2.4,2.4)

    hEta5_20_q12_cut20  = ROOT.TEfficiency("hEta20_%s_q12"%dataset,"",240,-2.4,2.4)
    hEta5_20_q12_cut0   = ROOT.TEfficiency("hEta20_%s_q12"%dataset,"",240,-2.4,2.4)
    hEta10_q12_cut5     = ROOT.TEfficiency("hEta20_%s_q12"%dataset,"",240,-2.4,2.4)
    hEta5_q12_cut5      = ROOT.TEfficiency("hEta20_%s_q12"%dataset,"",240,-2.4,2.4)
    hEta7_10_q12_cut5   = ROOT.TEfficiency("hEta20_%s_q12"%dataset,"",240,-2.4,2.4)
    hEta7_15_q12_cut5   = ROOT.TEfficiency("hEta20_%s_q12"%dataset,"",240,-2.4,2.4)
    hEta7_20_q12_cut5   = ROOT.TEfficiency("hEta20_%s_q12"%dataset,"",240,-2.4,2.4)


    hPt3_q12_5GeV     = ROOT.TEfficiency("hPt3_%s_q12_5"%dataset,"",20,0,100)
    hPt10_q12_5GeV    = ROOT.TEfficiency("hPt10_%s_q12_5"%dataset,"",20,0,100)
    hPt20_q12_5GeV    = ROOT.TEfficiency("hPt20_%s_q12_5"%dataset,"",20,0,100)
    hPt3_q12_3GeV     = ROOT.TEfficiency("hPt3_%s_q12_3"%dataset,"",34,0,102)
    hPt10_q12_3GeV    = ROOT.TEfficiency("hPt10_%s_q12_3"%dataset,"",34,0,102)
    hPt20_q12_3GeV    = ROOT.TEfficiency("hPt20_%s_q12_3"%dataset,"",34,0,102)
    hPt3_q12_2GeV     = ROOT.TEfficiency("hPt3_%s_q12_2"%dataset,"",50,0,100)
    hPt10_q12_2GeV    = ROOT.TEfficiency("hPt10_%s_q12_2"%dataset,"",50,0,100)
    hPt20_q12_2GeV    = ROOT.TEfficiency("hPt20_%s_q12_2"%dataset,"",50,0,100)
    hPt3_q12_1GeV     = ROOT.TEfficiency("hPt3_%s_q12_1"%dataset,"",100,0,100)
    hPt10_q12_1GeV    = ROOT.TEfficiency("hPt10_%s_q12_1"%dataset,"",100,0,100)
    hPt20_q12_1GeV    = ROOT.TEfficiency("hPt20_%s_q12_1"%dataset,"",100,0,100)

    outputDict[dataset] = { }
    outputDict[dataset]['hEta20_q12_cut20']   = hEta20_q12_cut20
    outputDict[dataset]['hEta5_20_q12_cut20'] = hEta5_20_q12_cut20
    outputDict[dataset]['hEta5_20_q12_cut0']  = hEta5_20_q12_cut0
    outputDict[dataset]['hEta10_q12_cut5']    = hEta10_q12_cut5
    outputDict[dataset]['hEta5_q12_cut5']     = hEta5_q12_cut5
    outputDict[dataset]['hEta7_10_q12_cut5'] = hEta7_10_q12_cut5
    outputDict[dataset]['hEta7_15_q12_cut5'] = hEta7_15_q12_cut5
    outputDict[dataset]['hEta7_20_q12_cut5'] = hEta7_20_q12_cut5

    outputDict[dataset]['eff_pt3_q12_5GeV']     = hPt3_q12_5GeV
    outputDict[dataset]['eff_pt10_q12_5GeV']    = hPt10_q12_5GeV
    outputDict[dataset]['eff_pt20_q12_5GeV']    = hPt20_q12_5GeV

    outputDict[dataset]['eff_pt3_q12_3GeV']     = hPt3_q12_3GeV
    outputDict[dataset]['eff_pt10_q12_3GeV']    = hPt10_q12_3GeV
    outputDict[dataset]['eff_pt20_q12_3GeV']    = hPt20_q12_3GeV

    outputDict[dataset]['eff_pt3_q12_2GeV']     = hPt3_q12_2GeV
    outputDict[dataset]['eff_pt10_q12_2GeV']    = hPt10_q12_2GeV
    outputDict[dataset]['eff_pt20_q12_2GeV']    = hPt20_q12_2GeV

    outputDict[dataset]['eff_pt3_q12_1GeV']     = hPt3_q12_1GeV
    outputDict[dataset]['eff_pt10_q12_1GeV']    = hPt10_q12_1GeV
    outputDict[dataset]['eff_pt20_q12_1GeV']    = hPt20_q12_1GeV


    print 'here :)'
    count = 0
    for ev in events:
        if not count%1000:  print count, events.size()
        count = count + 1 
        if count > maxEvents and maxEvents > 0: break
        ev.getByLabel(muonLabel, muonHandle)
        ev.getByLabel(genLabel, genHandle)
        muons = muonHandle.product()
        gens  = genHandle.product()
        goodmuons = [] 
        goodIndex = 0
        for bxNumber in range(muons.getFirstBX(), muons.getLastBX()+1):
           
            size = muons.size(bxNumber)
            for i in range(size):
                muon = muons[i+goodIndex]
                if muon.trackFinderType() not in [1,2]: continue
                #if muon.hwQual() < 12: continue
    
                goodmuons.append( muon) 

            goodIndex = goodIndex + size
        #print "Next event!"    
        for gen in gens:
            if abs(gen.pdgId()) != 13: continue
            passes_0_12   = False
            passes_3_12  = False
            passes_5_12  = False
            passes_10_12 = False
            passes_20_12 = False
            v_gen = getP4(gen)
            toprint = v_gen.Pt() > 60 and v_gen.Pt() < 80 and v_gen.Eta() >= 0.8 and v_gen.Eta() <= 1.3
            #if toprint: print "GEN: ", v_gen.Pt(), v_gen.Eta(), v_gen.Phi() 
            rematchedPt = 0
            for mu in goodmuons:
                v_mu = getP4FromHW(mu)
                #if toprint: print "L1 Cand: ", v_mu.Pt(), v_mu.Eta(), v_mu.Phi(), muon.hwQual()

                if v_mu.DeltaR( v_gen ) < 0.5 and mu.hwQual() >= 12: 
                    if v_mu.Pt() >= 20: passes_20_12 = True
                    if v_mu.Pt() >= 10: passes_10_12 = True
                    if v_mu.Pt() >= 5 : passes_5_12  = True
                    if v_mu.Pt() >= 3 : passes_3_12  = True
                    rematchedPt = v_mu.Pt() 
                    passes_0_12 = True
                    if passes_20_12: break

            if abs(gen.eta()) < 1.24 and abs(gen.eta()) > 0.82: hPt3_q12_5GeV.Fill( passes_3_12, v_gen.Pt())
            if abs(gen.eta()) < 1.24 and abs(gen.eta()) > 0.82: hPt10_q12_5GeV.Fill( passes_10_12, v_gen.Pt())
            if abs(gen.eta()) < 1.24 and abs(gen.eta()) > 0.82: hPt20_q12_5GeV.Fill( passes_20_12, v_gen.Pt())
            if abs(gen.eta()) < 1.24 and abs(gen.eta()) > 0.82: hPt3_q12_3GeV.Fill( passes_3_12, v_gen.Pt())
            if abs(gen.eta()) < 1.24 and abs(gen.eta()) > 0.82: hPt10_q12_3GeV.Fill( passes_10_12, v_gen.Pt())
            if abs(gen.eta()) < 1.24 and abs(gen.eta()) > 0.82: hPt20_q12_3GeV.Fill( passes_20_12, v_gen.Pt())
            if abs(gen.eta()) < 1.24 and abs(gen.eta()) > 0.82: hPt3_q12_2GeV.Fill( passes_3_12, v_gen.Pt())
            if abs(gen.eta()) < 1.24 and abs(gen.eta()) > 0.82: hPt10_q12_2GeV.Fill( passes_10_12, v_gen.Pt())
            if abs(gen.eta()) < 1.24 and abs(gen.eta()) > 0.82: hPt20_q12_2GeV.Fill( passes_20_12, v_gen.Pt())
            if abs(gen.eta()) < 1.24 and abs(gen.eta()) > 0.82: hPt3_q12_1GeV.Fill( passes_3_12, v_gen.Pt())
            if abs(gen.eta()) < 1.24 and abs(gen.eta()) > 0.82: hPt10_q12_1GeV.Fill( passes_10_12, v_gen.Pt())
            if abs(gen.eta()) < 1.24 and abs(gen.eta()) > 0.82: hPt20_q12_1GeV.Fill( passes_20_12, v_gen.Pt())
            #print passes_5_0, v_gen.Pt()
            #Now to the fun part with the eta
            if (v_gen.Pt() > 25):
                hEta20_q12_cut20.Fill(passes_20_12, abs(v_gen.Eta()))

            if (v_gen.Pt() > 5 and v_gen.Pt() < 20):
                hEta5_20_q12_cut20.Fill(passes_20_12, abs(v_gen.Eta()))

            if (v_gen.Pt() > 5 and v_gen.Pt() < 20):
                hEta5_20_q12_cut0.Fill(passes_0_12, abs(v_gen.Eta()))

            if (v_gen.Pt() > 10):
                hEta10_q12_cut5.Fill(passes_5_12, abs(v_gen.Eta()))

            if (v_gen.Pt() > 5):
                hEta5_q12_cut5.Fill(passes_5_12, abs(v_gen.Eta()))

            if (v_gen.Pt() > 7 and v_gen.Pt() < 10):
                hEta7_10_q12_cut5.Fill(passes_5_12, abs(v_gen.Eta()))

            if (v_gen.Pt() > 7 and v_gen.Pt() < 15):
                hEta7_15_q12_cut5.Fill(passes_5_12, abs(v_gen.Eta()))

            if (v_gen.Pt() > 7 and v_gen.Pt() < 20):
                hEta7_20_q12_cut5.Fill(passes_5_12, abs(v_gen.Eta()))

import pickle 
with open('tdr_eff_moreetabins_PU200_comparison_Mar09.pickle', 'wb') as handle:
    pickle.dump(outputDict, handle, protocol=pickle.HIGHEST_PROTOCOL)

