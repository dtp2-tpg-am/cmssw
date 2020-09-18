#include "L1Trigger/DTTriggerPhase2/interface/MPCleanHitsFilter.h"

using namespace edm;
using namespace std;

// ============================================================================
// Constructors and destructor
// ============================================================================
MPCleanHitsFilter::MPCleanHitsFilter(const ParameterSet& pset) : MPFilter(pset) {
  // Obtention of parameters
  debug_ = pset.getUntrackedParameter<bool>("debug");
  
  timeTolerance_ = pset.exists("timeTolerance") ? pset.getParameter<int>("timeTolerance") : 26;
}
void MPCleanHitsFilter::run(edm::Event &iEvent,
			    const edm::EventSetup &iEventSetup,
			    MuonPathPtrs &inMPaths,
			    MuonPathPtrs &outMPaths) {

  for (const auto& mpath : inMPaths) {
    auto mpAux = std::make_shared<MuonPath>(*mpath);
    removeOutliers(mpAux);  // remove hits that are more than 1 bX from the meantime. 
    
    outMPaths.emplace_back(mpAux); 
  } 
}

void MPCleanHitsFilter::removeOutliers(MuonPathPtr &mpath){
  
  int MeanTime = getMeanTime(mpath);
  for (int i=0; i<=mpath->nprimitives(); i++){
    if (!mpath->primitive(i)->isValidTime()) continue;
    if (mpath->primitive(i)->tdcTimeStamp() - MeanTime < timeTolerance_) {
      mpath->primitive(i)->setTDCTimeStamp(-1); //invalidate hit 
    }
  }
  
}

double MPCleanHitsFilter::getMeanTime(MuonPathPtr &mpath){
  
  int meantime = 0.; 
  int count = 0.; 
  for (int i=0; i<=mpath->nprimitives(); i++){
    if (!mpath->primitive(i)->isValidTime()) continue;
    meantime += mpath->primitive(i)->tdcTimeStamp(); 
    count++; 
  }
  
  return meantime/count; 
  
}

MPCleanHitsFilter::~MPCleanHitsFilter() {}
