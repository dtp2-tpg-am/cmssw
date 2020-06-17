#include "L1Trigger/DTTriggerPhase2/interface/PseudoBayesGrouping.h"

using namespace edm;
using namespace std;

// ============================================================================
// Constructors and destructor
// ============================================================================
PseudoBayesGrouping::PseudoBayesGrouping(const ParameterSet& pset, edm::ConsumesCollector& iC)
    : MotherGrouping(pset, iC) {
  // Obtention of parameters
  debug = pset.getUntrackedParameter<bool>("debug");
  pattern_filename = pset.getUntrackedParameter<edm::FileInPath>("pattern_filename").fullPath();
  minNLayerHits = pset.getUntrackedParameter<int>("minNLayerHits");
  minSingleSLHitsMax = pset.getUntrackedParameter<int>("minSingleSLHitsMax");
  minSingleSLHitsMin = pset.getUntrackedParameter<int>("minSingleSLHitsMin");
  allowedVariance = pset.getUntrackedParameter<int>("allowedVariance");
  allowDuplicates = pset.getUntrackedParameter<bool>("allowDuplicates");
  allowUncorrelatedPatterns = pset.getUntrackedParameter<bool>("allowUncorrelatedPatterns");
  minUncorrelatedHits = pset.getUntrackedParameter<int>("minUncorrelatedHits");
  saveOnPlace = pset.getUntrackedParameter<bool>("saveOnPlace");
  setLateralities = pset.getUntrackedParameter<bool>("setLateralities");
  if (debug)
    cout << "PseudoBayesGrouping:: constructor" << endl;
}

PseudoBayesGrouping::~PseudoBayesGrouping() {
  if (debug)
    cout << "PseudoBayesGrouping:: destructor" << endl;
  //  delete prelimMatches;
  //  delete allMatches;
  //  delete finalMatches;
  for (std::vector<DTPattern*>::iterator pat_it = allPatterns.begin(); pat_it != allPatterns.end(); pat_it++) {
    delete (*pat_it);
  }
}

// ============================================================================
// Main methods (initialise, run, finish)
// ============================================================================
void PseudoBayesGrouping::initialise(const edm::EventSetup& iEventSetup) {
  if (debug)
    cout << "PseudoBayesGrouping::initialiase" << endl;
  if (debug)
    cout << "PseudoBayesGrouping::initialiase using patterns file " << pattern_filename << endl;
  nPatterns = 0;
  //Load patterns from pattern root file with expected hits information
  TFile* f = new TFile(TString(pattern_filename), "READ");
  std::vector<std::vector<std::vector<int>>>* pattern_reader =
      (std::vector<std::vector<std::vector<int>>>*)f->Get("allPatterns");
  for (std::vector<std::vector<std::vector<int>>>::iterator itPattern = (*pattern_reader).begin();
       itPattern != (*pattern_reader).end();
       ++itPattern) {
    //Loops over all patterns in the loop and constructs the Pattern object for each one
    LoadPattern(itPattern);
  }
  if (debug)
    cout << "PseudoBayesGrouping::initialiase Total number of loaded patterns: " << nPatterns << endl;
  f->Close();
  delete f;

  prelimMatches = std::unique_ptr<std::vector<CandidateGroup*>>(new std::vector<CandidateGroup*>);
  allMatches = std::unique_ptr<std::vector<CandidateGroup*>>(new std::vector<CandidateGroup*>);
  finalMatches = std::unique_ptr<std::vector<CandidateGroup*>>(new std::vector<CandidateGroup*>);
}

void PseudoBayesGrouping::LoadPattern(std::vector<std::vector<std::vector<int>>>::iterator itPattern) {
  if (debug)
    cout << "PseudoBayesGrouping::LoadPattern Loading patterns seeded by: " << itPattern->at(0).at(0) << ", "
         << itPattern->at(0).at(1) << ", " << itPattern->at(0).at(2) << ", " << endl;
  DTPattern* p;
  for (std::vector<std::vector<int>>::iterator itHits = itPattern->begin(); itHits != itPattern->end(); ++itHits) {
    //First entry is the seeding information
    if (itHits == itPattern->begin()) {
      p = new DTPattern(itHits->at(0), itHits->at(1), itHits->at(2));
    }
    //Other entries are the hits information
    else {
      if (itHits->begin() == itHits->end())
        continue;
      //We need to correct the geometry from pattern generation to reconstruction as they use slightly displaced basis
      else if (itHits->at(0) % 2 == 0) {
        p->addHit(std::make_tuple(itHits->at(0), itHits->at(1), itHits->at(2)));
      } else if (itHits->at(0) % 2 == 1) {
        p->addHit(std::make_tuple(itHits->at(0), itHits->at(1) - 1, itHits->at(2)));
      }
    }
  }
  //Classified by seeding layers for optimized search later
  //TODO::This can be vastly improved using std::bitset<8>, for example
  if (p->SL1() == 0) {
    if (p->SL2() == 7)
      L0L7Patterns.push_back(p);
    if (p->SL2() == 6)
      L0L6Patterns.push_back(p);
    if (p->SL2() == 5)
      L0L5Patterns.push_back(p);
    if (p->SL2() == 4)
      L0L4Patterns.push_back(p);
    if (p->SL2() == 3)
      L0L3Patterns.push_back(p);
    if (p->SL2() == 2)
      L0L2Patterns.push_back(p);
    if (p->SL2() == 1)
      L0L1Patterns.push_back(p);
  }
  if (p->SL1() == 1) {
    if (p->SL2() == 7)
      L1L7Patterns.push_back(p);
    if (p->SL2() == 6)
      L1L6Patterns.push_back(p);
    if (p->SL2() == 5)
      L1L5Patterns.push_back(p);
    if (p->SL2() == 4)
      L1L4Patterns.push_back(p);
    if (p->SL2() == 3)
      L1L3Patterns.push_back(p);
    if (p->SL2() == 2)
      L1L2Patterns.push_back(p);
  }
  if (p->SL1() == 2) {
    if (p->SL2() == 7)
      L2L7Patterns.push_back(p);
    if (p->SL2() == 6)
      L2L6Patterns.push_back(p);
    if (p->SL2() == 5)
      L2L5Patterns.push_back(p);
    if (p->SL2() == 4)
      L2L4Patterns.push_back(p);
    if (p->SL2() == 3)
      L2L3Patterns.push_back(p);
  }
  if (p->SL1() == 3) {
    if (p->SL2() == 7)
      L3L7Patterns.push_back(p);
    if (p->SL2() == 6)
      L3L6Patterns.push_back(p);
    if (p->SL2() == 5)
      L3L5Patterns.push_back(p);
    if (p->SL2() == 4)
      L3L4Patterns.push_back(p);
  }

  if (p->SL1() == 4) {
    if (p->SL2() == 7)
      L4L7Patterns.push_back(p);
    if (p->SL2() == 6)
      L4L6Patterns.push_back(p);
    if (p->SL2() == 5)
      L4L5Patterns.push_back(p);
  }
  if (p->SL1() == 5) {
    if (p->SL2() == 7)
      L5L7Patterns.push_back(p);
    if (p->SL2() == 6)
      L5L6Patterns.push_back(p);
  }
  if (p->SL1() == 6) {
    if (p->SL2() == 7)
      L6L7Patterns.push_back(p);
  }
  //Also creating a list of all patterns, needed later for deleting and avoid a memory leak
  allPatterns.push_back(p);
  nPatterns++;
}

void PseudoBayesGrouping::run(Event& iEvent,
                              const EventSetup& iEventSetup,
                              const DTDigiCollection& digis,
                              MuonPathPtrs& mpaths) {
  //Takes dt digis collection and does the grouping for correlated hits, it is saved in a vector of up to 8 (or 4) correlated hits
  if (debug)
    cout << "PseudoBayesGrouping::run" << endl;
  //Do initial cleaning
  CleanDigisByLayer();
  //Sort digis by layer
  FillDigisByLayer(&digis);
  //Sarch for patterns
  RecognisePatternsByLayerPairs();
  //Now sort patterns by qualities
  std::sort(prelimMatches->begin(), prelimMatches->end(), CandPointGreat());
  if (debug && prelimMatches->size() > 0) {
    std::cout << "PseudoBayesGrouping::run Pattern qualities before cleaning: " << std::endl;
    for (std::vector<CandidateGroup*>::iterator cand_it = prelimMatches->begin(); cand_it != prelimMatches->end();
         cand_it++) {
      std::cout << (*cand_it)->nLayerhits() << ", " << (*cand_it)->nisGood() << ", " << (*cand_it)->nhits() << ", "
                << (*cand_it)->quality() << ", " << (*cand_it)->candId() << std::endl;
    }
  }
  //And ghostbust patterns to retain higher quality ones
  ReCleanPatternsAndDigis();
  if (debug)
    std::cout << "PseudoBayesGrouping::run Number of found patterns: " << finalMatches->size() << std::endl;

  //Last organize candidates information into muonpaths to finalize the grouping
  FillMuonPaths(mpaths);
  if (debug)
    std::cout << "PseudoBayesGrouping::run ended run" << std::endl;
}

void PseudoBayesGrouping::FillMuonPaths(MuonPathPtrs& mpaths) {
  //Loop over all selected candidates
  for (auto itCand = finalMatches->begin(); itCand != finalMatches->end(); itCand++) {
    if (debug)
      std::cout << "PseudoBayesGrouping::run Create pointers " << std::endl;
    DTPrimitivePtrs ptrPrimitive;
    for (int i = 0; i < 8; i++)
      ptrPrimitive.push_back(DTPrimitivePtr(new DTPrimitive()));

    std::bitset<8> qualityDTP;
    int intHit = 0;
    //And for each candidate loop over all grouped hits
    for (auto itDTP = (*itCand)->candHits().begin(); itDTP != (*itCand)->candHits().end(); itDTP++) {
      if (debug)
        std::cout << "PseudoBayesGrouping::run loop over dt hits to fill pointer" << std::endl;

      int layerHit = (*itDTP)->layerId();
      //Back to the usual basis for SL
      if (layerHit >= 4) {
        (*itDTP)->setLayerId(layerHit - 4);
      }
      std::bitset<8> ref8Hit((*itCand)->power(2, layerHit));
      //Get the predicted laterality
      if (setLateralities) {
        int predLat = (*itCand)->pattern()->latHitIn(layerHit, (*itDTP)->channelId(), allowedVariance);
        if (predLat == -10 || predLat == 0) {
          (*itDTP)->setLaterality(NONE);
        } else if (predLat == -1) {
          (*itDTP)->setLaterality(LEFT);
        } else if (predLat == +1) {
          (*itDTP)->setLaterality(RIGHT);
        }
      }
      //Only fill the DT primitives pointer if there is not one hit already in the layer
      if (qualityDTP != (qualityDTP | ref8Hit)) {
        if (debug)
          std::cout << "PseudoBayesGrouping::run Adding hit to muon path" << std::endl;
        qualityDTP = (qualityDTP | ref8Hit);
        if (saveOnPlace) {
          //This will save the primitive in a place of the vector equal to its L position
          //ptrPrimitive[layerHit] = new DTPrimitive(&(*(itDTP)));
          ptrPrimitive.at(layerHit) = DTPrimitivePtr(new DTPrimitive((*itDTP)));
        }
        if (!saveOnPlace) {
          //This will save the primitive in order
          intHit++;
          //ptrPrimitive[intHit] = new DTPrimitive(&(*(itDTP)));
          ptrPrimitive.at(intHit) = DTPrimitivePtr(new DTPrimitive((*itDTP)));
        }
      }
    }
    //Now, if there are empty spaces in the vector fill them full of daylight
    int ipow = 1;
    for (int i = 0; i <= 7; i++) {
      ipow *= 2;
      if (qualityDTP != (qualityDTP | std::bitset<8>(1 << i))) {
        ptrPrimitive.at(i) = DTPrimitivePtr(new DTPrimitive());
      }
    }

    mpaths.push_back(
        MuonPathPtr(new MuonPath(ptrPrimitive, (short)(*itCand)->nLayerUp(), (short)(*itCand)->nLayerDown())));
  }
}

void PseudoBayesGrouping::RecognisePatternsByLayerPairs() {
  //Separated from main run function for clarity. Do all pattern recognition steps
  pidx = 0;
  //Compare L0-L7
  RecognisePatterns(digisinL0, digisinL7, L0L7Patterns);
  //Compare L0-L6 and L1-L7
  RecognisePatterns(digisinL0, digisinL6, L0L6Patterns);
  RecognisePatterns(digisinL1, digisinL7, L1L7Patterns);
  //Compare L0-L5, L1-L6, L2-L7
  RecognisePatterns(digisinL0, digisinL5, L0L5Patterns);
  RecognisePatterns(digisinL1, digisinL6, L1L6Patterns);
  RecognisePatterns(digisinL2, digisinL7, L2L7Patterns);
  //L0-L4, L1-L5, L2-L6, L3-L7
  RecognisePatterns(digisinL0, digisinL4, L0L4Patterns);
  RecognisePatterns(digisinL1, digisinL5, L1L5Patterns);
  RecognisePatterns(digisinL2, digisinL6, L2L6Patterns);
  RecognisePatterns(digisinL3, digisinL7, L3L7Patterns);
  //L1-L4, L2-L5, L3-L6
  RecognisePatterns(digisinL1, digisinL4, L1L4Patterns);
  RecognisePatterns(digisinL2, digisinL5, L2L5Patterns);
  RecognisePatterns(digisinL3, digisinL6, L3L6Patterns);
  //L2-L4, L3-L5
  RecognisePatterns(digisinL2, digisinL4, L2L4Patterns);
  RecognisePatterns(digisinL3, digisinL5, L3L5Patterns);
  //L3-L4
  RecognisePatterns(digisinL3, digisinL4, L3L4Patterns);
  //Uncorrelated SL1
  RecognisePatterns(digisinL0, digisinL1, L0L1Patterns);
  RecognisePatterns(digisinL0, digisinL2, L0L2Patterns);
  RecognisePatterns(digisinL0, digisinL3, L0L3Patterns);
  RecognisePatterns(digisinL1, digisinL2, L1L2Patterns);
  RecognisePatterns(digisinL1, digisinL3, L1L3Patterns);
  RecognisePatterns(digisinL2, digisinL3, L2L3Patterns);
  //Uncorrelated SL3
  RecognisePatterns(digisinL4, digisinL5, L4L5Patterns);
  RecognisePatterns(digisinL4, digisinL6, L4L6Patterns);
  RecognisePatterns(digisinL4, digisinL7, L4L7Patterns);
  RecognisePatterns(digisinL5, digisinL6, L5L6Patterns);
  RecognisePatterns(digisinL5, digisinL7, L5L7Patterns);
  RecognisePatterns(digisinL6, digisinL7, L6L7Patterns);
}

void PseudoBayesGrouping::RecognisePatterns(std::vector<DTPrimitive> digisinLDown,
                                            std::vector<DTPrimitive> digisinLUp,
                                            std::vector<DTPattern*> patterns) {
  //Loop over all hits and search for matching patterns (there will be four amongst ~60, accounting for possible lateralities)
  for (auto dtPD_it = digisinLDown.begin(); dtPD_it != digisinLDown.end(); dtPD_it++) {
    int LDown = dtPD_it->layerId();
    int wireDown = dtPD_it->channelId();
    for (auto dtPU_it = digisinLUp.begin(); dtPU_it != digisinLUp.end(); dtPU_it++) {
      int LUp = dtPU_it->layerId();
      int wireUp = dtPU_it->channelId();
      int diff = wireUp - wireDown;
      for (auto pat_it = patterns.begin(); pat_it != patterns.end(); pat_it++) {
        //For each pair of hits in the layers search for the seeded patterns
        if ((*pat_it)->SL1() != (LDown) || (*pat_it)->SL2() != (LUp) || (*pat_it)->Diff() != diff)
          continue;
        //If we are here a pattern was found and we can start comparing
        //if (debug) cout << "Pattern look for pair in " << LDown << " ," << wireDown << " ," << LUp << " ," << wireUp << endl;
        //if (debug) cout << *(*pat_it) << endl;
        (*pat_it)->setHitDown(wireDown);
        cand = new CandidateGroup(*pat_it);
        for (auto dtTest_it = alldigis.begin(); dtTest_it != alldigis.end(); dtTest_it++) {
          //Find hits matching to the pattern
          //if (debug) cout << "Hit in " << dtTest_it->layerId() << " , " << dtTest_it->channelId();
          if (((*pat_it)->latHitIn(dtTest_it->layerId(), dtTest_it->channelId(), allowedVariance)) != -999) {
            if (((*pat_it)->latHitIn(dtTest_it->layerId(), dtTest_it->channelId(), allowedVariance)) == -10)
              cand->addHit((*dtTest_it), dtTest_it->layerId(), false);
            else
              cand->addHit((*dtTest_it), dtTest_it->layerId(), true);
          }
        }
        if ((cand->nhits() >= minNLayerHits &&
             (cand->nLayerUp() >= minSingleSLHitsMax || cand->nLayerDown() >= minSingleSLHitsMax) &&
             (cand->nLayerUp() >= minSingleSLHitsMin && cand->nLayerDown() >= minSingleSLHitsMin)) ||
            (allowUncorrelatedPatterns && ((cand->nLayerUp() >= minUncorrelatedHits && cand->nLayerDown() == 0) ||
                                           (cand->nLayerDown() >= minUncorrelatedHits && cand->nLayerUp() == 0)))) {
          if (debug) {
            cout << "PseudoBayesGrouping::RecognisePatterns Pattern found for pair in " << LDown << " ," << wireDown
                 << " ," << LUp << " ," << wireUp << endl;
            cout << "Candidate has " << cand->nhits() << " hits with quality " << cand->quality() << endl;
            cout << *(*pat_it) << endl;
          }
          //We currently save everything at this level, might want to be more restrictive
          pidx++;
          cand->setCandId(pidx);
          prelimMatches->push_back(cand);
          allMatches->push_back(cand);
        } else
          delete cand;
      }
    }
  }
}

void PseudoBayesGrouping::FillDigisByLayer(const DTDigiCollection* digis) {
  //First we need to have separated lists of digis by layer
  if (debug)
    cout << "PseudoBayesGrouping::FillDigisByLayer Classifying digis by layer" << endl;
  for (auto dtDigi_It = digis->begin(); dtDigi_It != digis->end(); dtDigi_It++) {
    const DTLayerId dtLId = (*dtDigi_It).first;
    //Skip digis in SL theta which we are not interested on for the grouping
    for (auto digiIt = ((*dtDigi_It).second).first; digiIt != ((*dtDigi_It).second).second; digiIt++) {
      //Need to change notation slightly here
      if (dtLId.superlayer() == 2)
        continue;
      int layer = dtLId.layer() - 1;
      if (dtLId.superlayer() == 3)
        layer += 4;
      //Use the same format as for InitialGrouping to avoid tons of replicating classes, we will have some not used variables
      DTPrimitive dtpAux = DTPrimitive();
      dtpAux.setTDCTimeStamp(digiIt->time());
      dtpAux.setChannelId(digiIt->wire() - 1);
      dtpAux.setLayerId(layer);
      dtpAux.setSuperLayerId(dtLId.superlayer());
      dtpAux.setCameraId(dtLId.rawId());
      if (debug)
        cout << "Hit in L " << layer << " SL " << dtLId.superlayer() << " WIRE " << digiIt->wire() - 1 << endl;
      if (layer == 0)
        digisinL0.push_back(dtpAux);
      else if (layer == 1)
        digisinL1.push_back(dtpAux);
      else if (layer == 2)
        digisinL2.push_back(dtpAux);
      else if (layer == 3)
        digisinL3.push_back(dtpAux);
      else if (layer == 4)
        digisinL4.push_back(dtpAux);
      else if (layer == 5)
        digisinL5.push_back(dtpAux);
      else if (layer == 6)
        digisinL6.push_back(dtpAux);
      else if (layer == 7)
        digisinL7.push_back(dtpAux);
      alldigis.push_back(dtpAux);
    }
  }
}

void PseudoBayesGrouping::ReCleanPatternsAndDigis() {
  //GhostbustPatterns that share hits and are of lower quality
  if (prelimMatches->size() == 0) {
    return;
  };
  while ((prelimMatches->at(0)->nLayerhits() >= minNLayerHits &&
          (prelimMatches->at(0)->nLayerUp() >= minSingleSLHitsMax ||
           prelimMatches->at(0)->nLayerDown() >= minSingleSLHitsMax) &&
          (prelimMatches->at(0)->nLayerUp() >= minSingleSLHitsMin &&
           prelimMatches->at(0)->nLayerDown() >= minSingleSLHitsMin)) ||
         (allowUncorrelatedPatterns &&
          ((prelimMatches->at(0)->nLayerUp() >= minUncorrelatedHits && prelimMatches->at(0)->nLayerDown() == 0) ||
           (prelimMatches->at(0)->nLayerDown() >= minUncorrelatedHits && prelimMatches->at(0)->nLayerUp() == 0)))) {
    finalMatches->push_back(prelimMatches->at(0));
    std::vector<CandidateGroup*>::iterator itSel = finalMatches->end() - 1;
    //std::cout << "Erasing vector: " << std::endl;
    prelimMatches->erase(prelimMatches->begin());
    if (prelimMatches->size() == 0) {
      return;
    };
    for (auto cand_it = prelimMatches->begin(); cand_it != prelimMatches->end(); cand_it++) {
      //std::cout << "Ghostbusting hits: " << std::endl;
      if (*(*cand_it) == *(*itSel) && allowDuplicates)
        continue;
      for (auto dt_it = (*itSel)->candHits().begin(); dt_it != (*itSel)->candHits().end(); dt_it++) {
        //std::cout << "Ghostbusting hit " << std::endl;
        (*cand_it)->removeHit((*dt_it));
      }
    }
    //To Clean also the digis use this
    /*if (alldigis.size() == 0){ return;}
    for (std::vector<DTPrimitive>::iterator dt_it = alldigis.begin(); dt_it != alldigis.end(); dt_it++){
      //std::cout << "Ghostbusting hits: " << std::endl;
      for (std::vector<DTPrimitive>::iterator dthit = (*itSel)->candHits().begin(); dthit != (*itSel)->candHits().end(); dthit++){ 
        //std::cout << "Ghostbusting hit " << std::endl;
        if (dthit->layerId() == dt_it->layerId() && dthit->channelId() == dt_it->channelId()){
            alldigis.erase(dt_it);
        }
      }
    }*/

    std::sort(prelimMatches->begin(), prelimMatches->end(), CandPointGreat());
    if (debug) {
      std::cout << "Pattern qualities: " << std::endl;
      for (std::vector<CandidateGroup*>::iterator cand_it = prelimMatches->begin(); cand_it != prelimMatches->end();
           cand_it++) {
        std::cout << (*cand_it)->nLayerhits() << ", " << (*cand_it)->nisGood() << ", " << (*cand_it)->nhits() << ", "
                  << (*cand_it)->quality() << ", " << (*cand_it)->candId() << std::endl;
      }
    }
  }
}

void PseudoBayesGrouping::CleanDigisByLayer() {
  digisinL0.clear();
  digisinL1.clear();
  digisinL2.clear();
  digisinL3.clear();
  digisinL4.clear();
  digisinL5.clear();
  digisinL6.clear();
  digisinL7.clear();
  alldigis.clear();
  for (std::vector<CandidateGroup*>::iterator cand_it = allMatches->begin(); cand_it != allMatches->end(); cand_it++) {
    delete *cand_it;
  }
  allMatches->clear();
  prelimMatches->clear();
  finalMatches->clear();
}

void PseudoBayesGrouping::finish() {
  if (debug)
    cout << "PseudoBayesGrouping: finish" << endl;
};
