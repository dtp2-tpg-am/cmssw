#include "L1Trigger/DTTriggerPhase2/interface/MuonPathAnalyzerInChamber.h"
#include <cmath>
#include <memory>

using namespace edm;
using namespace std;
using namespace cmsdt;
// ============================================================================
// Constructors and destructor
// ============================================================================
MuonPathAnalyzerInChamber::MuonPathAnalyzerInChamber(const ParameterSet &pset, edm::ConsumesCollector &iC)
    : MuonPathAnalyzer(pset, iC),
      debug_(pset.getUntrackedParameter<bool>("debug")),
      chi2Th_(pset.getUntrackedParameter<double>("chi2Th")),
      shift_filename_(pset.getParameter<edm::FileInPath>("shift_filename")),
      bxTolerance_(30),
      minQuality_(H3PLUS0),
      chiSquareThreshold_(50),
      minHits4Fit_(pset.getUntrackedParameter<int>("minHits4Fit")),
      splitPathPerSL_(pset.getUntrackedParameter<bool>("splitPathPerSL")){
  // Obtention of parameters

  if (debug_)
    LogDebug("MuonPathAnalyzerInChamber") << "MuonPathAnalyzer: constructor";

  setChiSquareThreshold(chi2Th_ * 10000.);

  //shift
  int rawId;
  std::ifstream ifin3(shift_filename_.fullPath());
  double shift;
  if (ifin3.fail()) {
    throw cms::Exception("Missing Input File")
        << "MuonPathAnalyzerInChamber::MuonPathAnalyzerInChamber() -  Cannot find " << shift_filename_.fullPath();
  }
  while (ifin3.good()) {
    ifin3 >> rawId >> shift;
    shiftinfo_[rawId] = shift;
  }

  dtGeomH = iC.esConsumes<DTGeometry, MuonGeometryRecord, edm::Transition::BeginRun>();
}

MuonPathAnalyzerInChamber::~MuonPathAnalyzerInChamber() {
  if (debug_)
    LogDebug("MuonPathAnalyzerInChamber") << "MuonPathAnalyzer: destructor";
}

// ============================================================================
// Main methods (initialise, run, finish)
// ============================================================================
void MuonPathAnalyzerInChamber::initialise(const edm::EventSetup &iEventSetup) {
  if (debug_)
    LogDebug("MuonPathAnalyzerInChamber") << "MuonPathAnalyzerInChamber::initialiase";

  const MuonGeometryRecord &geom = iEventSetup.get<MuonGeometryRecord>();
  dtGeo_ = &geom.get(dtGeomH);
}

void MuonPathAnalyzerInChamber::run(edm::Event &iEvent,
                                    const edm::EventSetup &iEventSetup,
                                    MuonPathPtrs &muonpaths,
                                    MuonPathPtrs &outmuonpaths) {
  if (debug_)
    LogDebug("MuonPathAnalyzerInChamber") << "MuonPathAnalyzerInChamber: run";

  // fit per SL (need to allow for multiple outputs for a single mpath)
  int nMuonPath_counter = 0;
  for (auto muonpath = muonpaths.begin(); muonpath != muonpaths.end(); ++muonpath) {

    if (debug_) {
      LogDebug("MuonPathAnalyzerInChamber")
	<< "Full path: "
	<< nMuonPath_counter << " , " << muonpath->get()->nprimitives() << " , " 
	<< muonpath->get()->nprimitivesUp() << " , " << muonpath->get()->nprimitivesDown();

      cout << "Full path: "
	<< nMuonPath_counter << " , " << muonpath->get()->nprimitives() << " , " 
	<< muonpath->get()->nprimitivesUp() << " , " << muonpath->get()->nprimitivesDown() 
	   << endl;
    }
    ++nMuonPath_counter;

    // Define muonpaths for up/down SL only
    MuonPathPtr muonpathUp_ptr = std::make_shared<MuonPath>();
    muonpathUp_ptr->setNPrimitives(8);
    muonpathUp_ptr->setNPrimitivesUp(muonpath->get()->nprimitivesUp());
    muonpathUp_ptr->setNPrimitivesDown(0);

    MuonPathPtr muonpathDown_ptr = std::make_shared<MuonPath>();
    muonpathDown_ptr->setNPrimitives(8);
    muonpathDown_ptr->setNPrimitivesUp(0);
    muonpathDown_ptr->setNPrimitivesDown(muonpath->get()->nprimitivesDown());

    for (int n = 0; n < muonpath->get()->nprimitives(); ++n){
      DTPrimitivePtr prim = muonpath->get()->primitive(n);
      // UP
      if (prim->superLayerId() == 3){
    	muonpathUp_ptr->setPrimitive(prim, n);
      }
      // DOWN
      else if (prim->superLayerId() == 1){
    	muonpathDown_ptr->setPrimitive(prim, n);
      }
      // NOT UP NOR DOWN
      else continue;
    }

    // cout << "Fitting " << endl;

    analyze(*muonpath, outmuonpaths);

    // cout << "Full path fitted" << endl; 
    
    if (splitPathPerSL_){
      // Analyze only if original muonpath was correlated or confirmed.
      // Otherwise, it's just repeating the original muonpath fit
      
      if (muonpathUp_ptr->nprimitivesUp() > 1 && muonpath->get()->nprimitivesDown() > 0)
	analyze(muonpathUp_ptr, outmuonpaths);

      if (muonpathDown_ptr->nprimitivesDown() > 1 && muonpath->get()->nprimitivesUp() > 0)
	analyze(muonpathDown_ptr, outmuonpaths);
    }
  }
}

void MuonPathAnalyzerInChamber::finish() {
  if (debug_)
    LogDebug("MuonPathAnalyzerInChamber") << "MuonPathAnalyzer: finish";
};

//------------------------------------------------------------------
//--- Private method
//------------------------------------------------------------------
void MuonPathAnalyzerInChamber::analyze(MuonPathPtr &inMPath, MuonPathPtrs &outMPath) {
  if (debug_)
    LogDebug("MuonPathAnalyzerInChamber") << "DTp2:analyze \t\t\t\t starts";

  // Clone the analyzed object
  if (debug_)
    LogDebug("MuonPathAnalyzerInChamber") << inMPath->nprimitives();
  auto mPath = std::make_shared<MuonPath>(inMPath);

  if (debug_) {
    LogDebug("MuonPathAnalyzerInChamber") << "DTp2::analyze, looking at mPath: ";
    for (int i = 0; i < mPath->nprimitives(); i++){
      LogDebug("MuonPathAnalyzerInChamber")
          << mPath->primitive(i)->layerId() << " , " << mPath->primitive(i)->superLayerId() << " , "
          << mPath->primitive(i)->channelId() << " , " << mPath->primitive(i)->laterality();

      cout << "Layer = " << mPath->primitive(i)->layerId() << " , " 
	   << "SL = "    << mPath->primitive(i)->superLayerId() << " , "
	   << "channel ID = " << mPath->primitive(i)->channelId() << " , " 
	   << "laterality = " << mPath->primitive(i)->laterality() << " , "
	   << "orbit = " << mPath->primitive(i)->orbit() << " , "
	   << "time correction = " << mPath->primitive(i)->timeCorrection() << " , "
	   << "cameraId = " << mPath->primitive(i)->cameraId() << " , "
	   << endl;
    }
  }

  if (debug_)
    LogDebug("MuonPathAnalyzerInChamber") << "DTp2:analyze \t\t\t\t\t is Analyzable? ";
  if (!mPath->isAnalyzable())
    return;
  if (debug_)
    LogDebug("MuonPathAnalyzerInChamber") << "DTp2:analyze \t\t\t\t\t yes it is analyzable " << mPath->isAnalyzable();


  // cout << "Building lateralities and setting wires positions" << endl;

  // first of all, get info from primitives, so we can reduce the number of latereralities:
  buildLateralities(mPath);
  setWirePosAndTimeInMP(mPath);

  std::shared_ptr<MuonPath> mpAux;
  int bestI = -1;
  float best_chi2 = 99999.;
  for (int i = 0; i < totalNumValLateralities_; i++) {  // LOOP for all lateralities:
    if (debug_)
      LogDebug("MuonPathAnalyzerInChamber") << "DTp2:analyze \t\t\t\t\t Start with combination " << i;
    int NTotalHits = NUM_LAYERS_2SL;
    float xwire[NUM_LAYERS_2SL];
    int present_layer[NUM_LAYERS_2SL];
    for (int ii = 0; ii < 8; ii++) {
      xwire[ii] = mPath->xWirePos(ii);
      if (xwire[ii] == 0) {
        present_layer[ii] = 0;
        NTotalHits--;
      } 
      else {
        present_layer[ii] = 1;
      }
    }
    
    // why do we repeat the fit reducing NTotalHits? What changes?
    while (NTotalHits >= minHits4Fit_) {
      cout << "Total hits: " << NTotalHits << endl;
      mPath->setChiSquare(0);
      calculateFitParameters(mPath, lateralities_[i], present_layer);
      if (mPath->chiSquare() != 0)
        break;
      NTotalHits--;
    }
    
    cout << "" << endl;
    cout << " Chi2 = " << mPath->chiSquare() << ", to be compared with a best Chi2 of " << best_chi2 
	 <<  " and a Chi2 threshold of " << chiSquareThreshold_ << endl;
    if (mPath->chiSquare() > chiSquareThreshold_)
      continue;
    
    cout << "Chi2 cut passed" << endl;

    evaluateQuality(mPath);

    cout << "Quality = " << mPath->quality() << endl;

    if (mPath->quality() < minQuality_)
      continue;

    cout << "Quality cut passed" << endl;

    double z = 0.;
    
    // mPath->horizPos() is the position wrt the center of the chamber
    double jm_x = (mPath->horizPos());
    int selected_Id = 0;
    // Lower layer with a hit
    for (int i = 0; i < mPath->nprimitives(); i++) {
      if (mPath->primitive(i)->isValidTime()) {
        selected_Id = mPath->primitive(i)->cameraId();
        mPath->setRawId(selected_Id);
        break;
      }
    }

    DTLayerId thisLId(selected_Id);

    DTChamberId ChId(thisLId.wheel(), thisLId.station(), thisLId.sector());

    // cout << "MB = " << thisLId.station() << ", SL = " << thisLId.superLayer() << ",  layer = " << thisLId.layer() << endl;    
    
    if (thisLId.station() >= 3)
      // z = Z_SHIFT_MB4;
      z += Z_SHIFT_MB4;

    DTSuperLayerId MuonPathSLId(thisLId.wheel(), thisLId.station(), thisLId.sector(), thisLId.superLayer());

    // cout << "Jm_x = " << jm_x << endl;
    // cout << "Propagating x to SL1: " << jm_x + mPath->tanPhi() * (11.1 + 0.65) << endl;
    // cout << "Propagating x to SL3: " << jm_x - mPath->tanPhi() * (11.1 + 0.65) << endl;

    // Count hits in each SL
    int hits_in_SL1 = 0;
    int hits_in_SL3 = 0;
    for (int i = 0; i < mPath->nprimitives(); i++) {
      if (mPath->primitive(i)->isValidTime()) {
	if (i <= 3) ++hits_in_SL1;
	else if (i > 3) ++hits_in_SL3;
      }
    }

    cout << "Hits in SL1 = " << hits_in_SL1 << ", hits in SL3 = " << hits_in_SL3 << endl;

    GlobalPoint jm_x_cmssw_global;
    // Depending on which SL has hits, propagate jm_x to SL1, SL3, or to the center of the chamber
    if (hits_in_SL1 > 2 && hits_in_SL3 < 2){
      // Uncorrelated or confirmed with 3 or 4 hits in SL1: propagate to SL1
      cout << "Uncorrelated or confirmed with 3 or 4 hits in SL1: propagate to SL1" << endl;
      jm_x += mPath->tanPhi() * (11.1 + 0.65); 
      jm_x_cmssw_global = dtGeo_->chamber(ChId)->toGlobal(LocalPoint(jm_x, 0., z + 11.75));
      // cout << "Local Point = " << LocalPoint(jm_x, 0., z + 11.75) << endl;
    }
    else if (hits_in_SL1 < 2 && hits_in_SL3 > 2){
      // Uncorrelated or confirmed with 3 or 4 hits in SL3: propagate to SL3
      cout << "Uncorrelated or confirmed with 3 or 4 hits in SL3: propagate to SL3" << endl;
      jm_x -= mPath->tanPhi() * (11.1 + 0.65); 
      jm_x_cmssw_global = dtGeo_->chamber(ChId)->toGlobal(LocalPoint(jm_x, 0., z - 11.75));
      // cout << "Local Point = " << LocalPoint(jm_x, 0., z - 11.75) << endl;
    }
    else if (hits_in_SL1 > 2 && hits_in_SL3 > 2){ // || hits_in_SL1 >= 2 && hits_in_SL3 > 2){
      // Correlated: stay at chamber center
      cout << "Correlated: stay at chamber center" << endl;
      jm_x_cmssw_global = dtGeo_->chamber(ChId)->toGlobal(LocalPoint(jm_x, 0., z));
      // cout << "Local Point = " << LocalPoint(jm_x, 0., z) << endl;
    }
    else {
      // Not interesting
      cout << "Not interesting" << endl;
      continue;
    }

    // Updating muon-path horizontal position
    mPath->setHorizPos(jm_x);

    
    // cout << "Chamber ID = " << ChId << endl;
    // cout << "SuperLy ID = " << MuonPathSLId << endl;

    int thisec = MuonPathSLId.sector();
    if (thisec == 13)
      thisec = 4;
    if (thisec == 14)
      thisec = 10;

    // cout << "Global phi = " << jm_x_cmssw_global.phi() << endl;
    double phi = jm_x_cmssw_global.phi() - PHI_CONV * (thisec - 1);
    // cout << "phi = " << phi * 65536. / 0.8 << endl;

    double psi = atan(mPath->tanPhi());
    // cout << "psi = " << psi << endl;

    mPath->setPhi(jm_x_cmssw_global.phi() - PHI_CONV * (thisec - 1));
    mPath->setPhiB(hasPosRF(MuonPathSLId.wheel(), MuonPathSLId.sector()) ? psi - phi : -psi - phi);

    //    cout << "PhiB = " << mPath->phiB() * 2048. / 1.4 << endl;

    //     cout << "Evaluating Chi2" << endl;

    
    cout << " Chi2 = " << mPath->chiSquare() << ", to be compared with a best Chi2 of " << best_chi2 << endl;
    
    if (mPath->chiSquare() < best_chi2)
      cout << "Chi2 is smaller than best_chi2" << endl;
    else 
      cout << "Chi2 is larger than best_chi2" << endl;

    if (mPath->chiSquare() > 0)
      cout << "Chi2 is positive" << endl;
    else
      cout << "Chi2 is 0 or negative" << endl;


    if (mPath->chiSquare() < best_chi2 && mPath->chiSquare() > 0) {
      cout << "New mpAux!" << endl;
      mpAux = std::make_shared<MuonPath>(mPath);
      bestI = i;
      best_chi2 = mPath->chiSquare();
    }
  }

  cout << "Best laterality: ";
  for (int p = 0; p < 8; ++p){
    cout << mPath->primitive(p)->laterality() << " ";
  }
  cout << "" << endl;
 
  // cout << "Pushing back to output" << endl;

  if (mpAux != nullptr) {
    outMPath.push_back(std::move(mpAux));
    if (debug_){
      LogDebug("MuonPathAnalyzerInChamber")
          << "DTp2:analize \t\t\t\t\t Laterality " << bestI << " is the one with smaller chi2";
      cout << "DTp2:analize \t\t\t\t\t Laterality " << bestI << " is the one with smaller chi2" << endl;
    }
  } 
  else {
    if (debug_){
      LogDebug("MuonPathAnalyzerInChamber")
          << "DTp2:analize \t\t\t\t\t No Laterality found with chi2 smaller than threshold";
      cout << "DTp2:analize \t\t\t\t\t No Laterality found with chi2 smaller than threshold" << endl;
    }
  }
  if (debug_)
    LogDebug("MuonPathAnalyzerInChamber") << "DTp2:analize \t\t\t\t\t Ended working with this set of lateralities";

  // cout << "Finishing fitting" << endl;

}

void MuonPathAnalyzerInChamber::setCellLayout(MuonPathPtr &mpath) {
  for (int i = 0; i <= mpath->nprimitives(); i++) {
    if (mpath->primitive(i)->isValidTime())
      cellLayout_[i] = mpath->primitive(i)->channelId();
    else
      cellLayout_[i] = -99;
  }

  // putting info back into the mpath:
  mpath->setCellHorizontalLayout(cellLayout_);
  for (int i = 0; i <= mpath->nprimitives(); i++) {
    if (cellLayout_[i] >= 0) {
      mpath->setBaseChannelId(cellLayout_[i]);
      break;
    }
  }
}

/**
 * For a combination of up to 8 cells, build all the lateralities to be tested,
 */
void MuonPathAnalyzerInChamber::buildLateralities(MuonPathPtr &mpath) {
  if (debug_)
    LogDebug("MuonPathAnalyzerInChamber") << "MPAnalyzer::buildLateralities << setLateralitiesFromPrims ";

  // cout << "Inside buildLateralities" << endl;

  mpath->setLateralCombFromPrimitives();

  totalNumValLateralities_ = 0;
  lateralities_.clear();
  latQuality_.clear();

  /* We generate all the possible laterality combinations compatible with the built 
     group in the previous step */
  lateralities_.push_back(TLateralities());
  for (int ilat = 0; ilat < NLayers; ilat++) {
    // Get value from input
    LATERAL_CASES lr = (mpath->lateralComb())[ilat];
    if (debug_)
      LogDebug("MuonPathAnalyzerInChamber") << "[DEBUG_] Input[" << ilat << "]: " << lr;

    // If left/right fill number
    if (lr != NONE) {
      if (debug_)
        LogDebug("MuonPathAnalyzerInChamber") << "[DEBUG_]   - Adding it to " << lateralities_.size() << " lists...";
      for (unsigned int iall = 0; iall < lateralities_.size(); iall++) {
        lateralities_[iall][ilat] = lr;
      }
    }
    // both possibilites
    else {
      // Get the number of possible options now
      auto ncurrentoptions = lateralities_.size();

      // Duplicate them
      if (debug_)
        LogDebug("MuonPathAnalyzerInChamber") << "[DEBUG_]   - Duplicating " << ncurrentoptions << " lists...";
      copy(lateralities_.begin(), lateralities_.end(), back_inserter(lateralities_));
      if (debug_)
        LogDebug("MuonPathAnalyzerInChamber") << "[DEBUG_]   - Now we have " << lateralities_.size() << " lists...";

      // Asign LEFT to first ncurrentoptions and RIGHT to the last
      for (unsigned int iall = 0; iall < ncurrentoptions; iall++) {
        lateralities_[iall][ilat] = LEFT;
        lateralities_[iall + ncurrentoptions][ilat] = RIGHT;
      }
    }  // else
  }    // Iterate over input array

  totalNumValLateralities_ = (int)lateralities_.size();

  if (totalNumValLateralities_ > 128) {
    // ADD PROTECTION!
    LogDebug("MuonPathAnalyzerInChamber") << "[WARNING]: TOO MANY LATERALITIES TO CHECK !!";
    LogDebug("MuonPathAnalyzerInChamber") << "[WARNING]: skipping this muon";
    lateralities_.clear();
    latQuality_.clear();
    totalNumValLateralities_ = 0;
  }

  // Dump values
  if (debug_) {
    for (unsigned int iall = 0; iall < lateralities_.size(); iall++) {
      LogDebug("MuonPathAnalyzerInChamber") << iall << " -> [";
      for (int ilat = 0; ilat < NLayers; ilat++) {
        if (ilat != 0)
          LogDebug("MuonPathAnalyzerInChamber") << ",";
        LogDebug("MuonPathAnalyzerInChamber") << lateralities_[iall][ilat];
      }
      LogDebug("MuonPathAnalyzerInChamber") << "]";
    }
  }
}

void MuonPathAnalyzerInChamber::setLateralitiesInMP(MuonPathPtr &mpath, TLateralities lat) {
  LATERAL_CASES tmp[NUM_LAYERS_2SL];
  cout << "Current laterality: ";
  for (int i = 0; i < 8; i++){
    //tmp[i] = lat[i];
    mpath->primitive(i)->setLaterality(lat[i]);
    cout << " " << lat[i] << " ";
  }
  cout << "" <<endl;
  // Is this working? setLateralComb wants as input a vector of 4 elements
  //  mpath->setLateralComb(tmp);
}

void MuonPathAnalyzerInChamber::setWirePosAndTimeInMP(MuonPathPtr &mpath) {
  // Get lower layer with a hit
  int selected_Id = 0;
  for (int i = 0; i < mpath->nprimitives(); i++) {
    if (mpath->primitive(i)->isValidTime()) {
      selected_Id = mpath->primitive(i)->cameraId();
      mpath->setRawId(selected_Id);
      break;
    }
  }
  // Layer ID of lower layer with a hit
  DTLayerId thisLId(selected_Id);
  DTChamberId chId(thisLId.wheel(), thisLId.station(), thisLId.sector());
  if (debug_)
    LogDebug("MuonPathAnalyzerInChamber")
        << "Id " << chId.rawId() << " Wh " << chId.wheel() << " St " << chId.station() << " Se " << chId.sector();
  mpath->setRawId(chId.rawId());

  DTSuperLayerId MuonPathSLId1(thisLId.wheel(), thisLId.station(), thisLId.sector(), 1);
  DTSuperLayerId MuonPathSLId3(thisLId.wheel(), thisLId.station(), thisLId.sector(), 3);
  DTWireId wireId1(MuonPathSLId1, 2, 1);
  DTWireId wireId3(MuonPathSLId3, 2, 1);

  if (debug_)
    LogDebug("MuonPathAnalyzerInChamber")
        << "shift1=" << shiftinfo_[wireId1.rawId()] << " shift3=" << shiftinfo_[wireId3.rawId()];

  float delta = 42000;                                                                      //um
  float zwire[NUM_LAYERS_2SL] = {-13.7, -12.4, -11.1, -9.8002, 9.79999, 11.1, 12.4, 13.7};  // mm
  for (int i = 0; i < mpath->nprimitives(); i++) {
    if (mpath->primitive(i)->isValidTime()) {
      if (i < 4)
        mpath->setXWirePos(10000 * shiftinfo_[wireId1.rawId()] +
                               (mpath->primitive(i)->channelId() + 0.5 * (double)((i + 1) % 2)) * delta,
                           i);
      if (i >= 4)
        mpath->setXWirePos(10000 * shiftinfo_[wireId3.rawId()] +
                               (mpath->primitive(i)->channelId() + 0.5 * (double)((i + 1) % 2)) * delta,
                           i);
      mpath->setZWirePos(zwire[i] * 1000, i);  // in um
      mpath->setTWireTDC(mpath->primitive(i)->tdcTimeStamp() * DRIFT_SPEED, i);
    } else {
      mpath->setXWirePos(0., i);
      mpath->setZWirePos(0., i);
      mpath->setTWireTDC(-1 * DRIFT_SPEED, i);
    }
    if (debug_)
      LogDebug("MuonPathAnalyzerInChamber") << mpath->primitive(i)->tdcTimeStamp() << " ";
  }
  if (debug_)
    LogDebug("MuonPathAnalyzerInChamber");
}

void MuonPathAnalyzerInChamber::calculateFitParameters(MuonPathPtr &mpath,
                                                       TLateralities laterality,
                                                       int present_layer[NUM_LAYERS_2SL]) {
  cout << "Calculating fit parameters" << endl;
  
  // First prepare mpath for fit:
  float xwire[NUM_LAYERS_2SL], zwire[NUM_LAYERS_2SL], tTDCvdrift[NUM_LAYERS_2SL];
  double b[NUM_LAYERS_2SL];
  for (int i = 0; i < 8; i++) {
    xwire[i] = mpath->xWirePos(i);
    zwire[i] = mpath->zWirePos(i);
    tTDCvdrift[i] = mpath->tWireTDC(i);
    b[i] = 1;
  }

  //// NOW Start FITTING:

  // fill hit position
  float xhit[NUM_LAYERS_2SL];
  for (int lay = 0; lay < 8; lay++) {
    if (debug_)
      LogDebug("MuonPathAnalyzerInChamber") << "In fitPerLat " << lay << " xwire " << xwire[lay] << " zwire "
                                            << zwire[lay] << " tTDCvdrift " << tTDCvdrift[lay];
    xhit[lay] = xwire[lay] + (-1 + 2 * laterality[lay]) * 1000 * tTDCvdrift[lay];
    cout << "xwire = " << xwire[lay]
    	 << ", zwire = " << zwire[lay]
         << ", tTDCvdrift = " << tTDCvdrift[lay]
	 << ", laterality = " << laterality[lay] << endl;

    if (debug_)
      LogDebug("MuonPathAnalyzerInChamber") << "In fitPerLat " << lay << " xhit " << xhit[lay];
  }

  //Proceed with calculation of fit parameters
  double cbscal = 0.0;
  double zbscal = 0.0;
  double czscal = 0.0;
  double bbscal = 0.0;
  double zzscal = 0.0;
  double ccscal = 0.0;

  for (int lay = 0; lay < 8; lay++) {
    if (present_layer[lay] == 0)
      continue;
    if (debug_)
      LogDebug("MuonPathAnalyzerInChamber")
	<< " For layer " << lay + 1 << " xwire[lay] " << xwire[lay] << " zwire " << zwire[lay] << " b " << b[lay];
    if (debug_)
      LogDebug("MuonPathAnalyzerInChamber") << " xhit[lat][lay] " << xhit[lay];
    cbscal = (-1 + 2 * laterality[lay]) * b[lay] + cbscal;
    zbscal = zwire[lay] * b[lay] + zbscal;      // it actually does not depend on laterality
    czscal = (-1 + 2 * laterality[lay]) * zwire[lay] + czscal;

    bbscal = b[lay] * b[lay] + bbscal;          // it actually does not depend on laterality
    zzscal = zwire[lay] * zwire[lay] + zzscal;  // it actually does not depend on laterality
    ccscal = (-1 + 2 * laterality[lay]) * (-1 + 2 * laterality[lay]) + ccscal;
  }

  double cz = 0.0;
  double cb = 0.0;
  double zb = 0.0;
  double zc = 0.0;
  double bc = 0.0;
  double bz = 0.0;

  cz = (cbscal * zbscal - czscal * bbscal) / (zzscal * bbscal - zbscal * zbscal);
  cb = (czscal * zbscal - cbscal * zzscal) / (zzscal * bbscal - zbscal * zbscal);

  zb = (czscal * cbscal - zbscal * ccscal) / (bbscal * ccscal - cbscal * cbscal);
  zc = (zbscal * cbscal - czscal * bbscal) / (bbscal * ccscal - cbscal * cbscal);

  bc = (zbscal * czscal - cbscal * zzscal) / (ccscal * zzscal - czscal * czscal);
  bz = (cbscal * czscal - zbscal * ccscal) / (ccscal * zzscal - czscal * czscal);

  double c_tilde[NUM_LAYERS_2SL];
  double z_tilde[NUM_LAYERS_2SL];
  double b_tilde[NUM_LAYERS_2SL];

  for (int lay = 0; lay < 8; lay++) {
    if (present_layer[lay] == 0)
      continue;
    if (debug_)
      LogDebug("MuonPathAnalyzerInChamber")
          << " For layer " << lay + 1 << " xwire[lay] " << xwire[lay] << " zwire " << zwire[lay] << " b " << b[lay];
    c_tilde[lay] = (-1 + 2 * laterality[lay]) + cz * zwire[lay] + cb * b[lay];
    z_tilde[lay] = zwire[lay] + zb * b[lay] + zc * (-1 + 2 * laterality[lay]);
    b_tilde[lay] = b[lay] + bc * (-1 + 2 * laterality[lay]) + bz * zwire[lay];
  }

  //Calculate results per lat
  double xctilde = 0.0;
  double xztilde = 0.0;
  double xbtilde = 0.0;
  double ctildectilde = 0.0;
  double ztildeztilde = 0.0;
  double btildebtilde = 0.0;

  double rect0vdrift = 0.0;
  double recslope = 0.0;
  double recpos = 0.0;

  for (int lay = 0; lay < 8; lay++) {
    if (present_layer[lay] == 0)
      continue;
    xctilde = xhit[lay] * c_tilde[lay] + xctilde;
    ctildectilde = c_tilde[lay] * c_tilde[lay] + ctildectilde;
    xztilde = xhit[lay] * z_tilde[lay] + xztilde;
    ztildeztilde = z_tilde[lay] * z_tilde[lay] + ztildeztilde;
    xbtilde = xhit[lay] * b_tilde[lay] + xbtilde;
    btildebtilde = b_tilde[lay] * b_tilde[lay] + btildebtilde;
  }

  //Results for t0vdrift (BX), slope and position per lat
  rect0vdrift = xctilde / ctildectilde;
  recslope = xztilde / ztildeztilde;
  recpos = xbtilde / btildebtilde;
  if (debug_) {
    LogDebug("MuonPathAnalyzerInChamber") << " In fitPerLat Reconstructed values per lat "
                                          << " rect0vdrift " << rect0vdrift;
    LogDebug("MuonPathAnalyzerInChamber")
        << "rect0 " << rect0vdrift / DRIFT_SPEED << " recBX " << rect0vdrift / DRIFT_SPEED / 25 << " recslope "
        << recslope << " recpos " << recpos;
  }

  //Get t*v and residuals per layer, and chi2 per laterality
  double rectdriftvdrift[NUM_LAYERS_2SL];
  double recres[NUM_LAYERS_2SL];
  double recchi2 = 0.0;
  int sign_tdriftvdrift = {0};
  int incell_tdriftvdrift = {0};
  int physical_slope = {0};

  // Select the worst hit in order to get rid of it
  double maxDif = -1;
  int maxInt = -1;

  for (int lay = 0; lay < 8; lay++) {
    if (present_layer[lay] == 0)
      continue;
    rectdriftvdrift[lay] = tTDCvdrift[lay] - rect0vdrift / 1000;
    if (debug_)
      LogDebug("MuonPathAnalyzerInChamber") << rectdriftvdrift[lay];
    recres[lay] = xhit[lay] - zwire[lay] * recslope - b[lay] * recpos - (-1 + 2 * laterality[lay]) * rect0vdrift;
    if ((present_layer[lay] == 1) && (rectdriftvdrift[lay] < -0.1)) {
      sign_tdriftvdrift = -1;
      if (-0.1 - rectdriftvdrift[lay] > maxDif) {
        maxDif = -0.1 - rectdriftvdrift[lay];
        maxInt = lay;
      }
    }
    if ((present_layer[lay] == 1) && (abs(rectdriftvdrift[lay]) > 21.1)) {
      incell_tdriftvdrift = -1;  //Changed to 21.1 to account for resolution effects
      if (rectdriftvdrift[lay] - 21.1 > maxDif) {
        maxDif = rectdriftvdrift[lay] - 21.1;
        maxInt = lay;
      }
    }
  }

  if (fabs(recslope / 10) > 1.3)
    physical_slope = -1;

  if (physical_slope == -1 && debug_)
    LogDebug("MuonPathAnalyzerInChamber") << "Combination with UNPHYSICAL slope ";
  if (sign_tdriftvdrift == -1 && debug_)
    LogDebug("MuonPathAnalyzerInChamber") << "Combination with negative tdrift-vdrift ";
  if (incell_tdriftvdrift == -1 && debug_)
    LogDebug("MuonPathAnalyzerInChamber") << "Combination with tdrift-vdrift larger than half cell ";

  for (int lay = 0; lay < 8; lay++) {
    if (present_layer[lay] == 0)
      continue;
    recchi2 = recres[lay] * recres[lay] + recchi2;
  }
  if (debug_)
    LogDebug("MuonPathAnalyzerInChamber")
        << "In fitPerLat Chi2 " << recchi2 << " with sign " << sign_tdriftvdrift << " within cell "
        << incell_tdriftvdrift << " physical_slope " << physical_slope;

  // LATERALITY IS NOT VALID - is this actually used?
  if (true && maxInt != -1) {
    present_layer[maxInt] = 0;
    if (debug_)
      LogDebug("MuonPathAnalyzerInChamber") << "We get rid of hit in layer " << maxInt;
  }

  // LATERALITY IS VALID...
  if (!(sign_tdriftvdrift == -1) && !(incell_tdriftvdrift == -1) && !(physical_slope == -1)) {
    mpath->setBxTimeValue((rect0vdrift / DRIFT_SPEED) / 1000);
    mpath->setTanPhi(-1 * recslope / 10);
    mpath->setHorizPos(recpos / 10000);
    mpath->setChiSquare(recchi2 / 100000000);
    setLateralitiesInMP(mpath, laterality);

    if (debug_)
      LogDebug("MuonPathAnalyzerInChamber")
          << "In fitPerLat "
          << "t0 " << mpath->bxTimeValue() << " slope " << mpath->tanPhi() << " pos " << mpath->horizPos() << " chi2 "
          << mpath->chiSquare() << " rawId " << mpath->rawId();

    cout << "In fitPerLat "
	 << "t0 " << mpath->bxTimeValue() << " slope " << mpath->tanPhi() 
	 << " pos " << mpath->horizPos() << " chi2 "
	 << mpath->chiSquare() << " rawId " << mpath->rawId();
  }
}

void MuonPathAnalyzerInChamber::evaluateQuality(MuonPathPtr &mPath) {
  mPath->setQuality(NOPATH);
  
  cout << " " << endl;  
  cout << "Lateralities after the fit: ";
  
  int validHits(0),nPrimsUp(0),nPrimsDown(0);
  for (int i=0; i<NUM_LAYERS_2SL; i++) {

    cout << " " << mPath->primitive(i)->laterality() << " ";
    if (mPath->primitive(i)->isValidTime()){
      validHits++;       
      if      ( i<4  ) nPrimsDown++; 
      else if ( i>=4 ) nPrimsUp++; 
    }
  }

  cout << " " << endl;

  mPath->setNPrimitivesUp(nPrimsUp);
  mPath->setNPrimitivesDown(nPrimsDown);
  

  // NEW QUALITY ASSIGNMENT
  // 4 + 4 --> H4PLUS4 --> Q = 8
  if (mPath->nprimitivesUp() >= 4 && mPath->nprimitivesDown() >= 4) {
    mPath->setQuality(H4PLUS4);
  }
  // 4 + 3 --> H4PLUS3 --> Q = 7
  else if ((mPath->nprimitivesUp() == 4 && mPath->nprimitivesDown() == 3) ||
	   (mPath->nprimitivesUp() == 3 && mPath->nprimitivesDown() == 4)) {
    mPath->setQuality(H4PLUS3);
  }
  // 3 + 3 --> H3PLUS3 --> Q = 6
  else if (mPath->nprimitivesUp() == 3 && mPath->nprimitivesDown() == 3) {
    mPath->setQuality(H3PLUS3);
  }
  // 4 + 2 --> H4PLUS2 --> Q = 4
  else if ((mPath->nprimitivesUp() == 4 && mPath->nprimitivesDown() == 2) ||
	   (mPath->nprimitivesUp() == 2 && mPath->nprimitivesDown() == 4)) {
    mPath->setQuality(H4PLUS2);
  }
  // (4 + 0) or (4 + 1) --> H4PLUS0 --> Q = 3
  else if ((mPath->nprimitivesUp() == 4 && 
	    (mPath->nprimitivesDown() == 1 || mPath->nprimitivesDown() == 0)) ||
	   ((mPath->nprimitivesUp() == 1 || mPath->nprimitivesUp() == 0) &&
	    mPath->nprimitivesDown() == 4)){
    mPath->setQuality(H4PLUS0);
  }
  // 3 + 2 --> H3PLUS2 --> Q = 2
  else if ((mPath->nprimitivesUp() == 3 && mPath->nprimitivesDown() == 2) ||
	   (mPath->nprimitivesUp() == 2 && mPath->nprimitivesDown() == 3)){
    mPath->setQuality(H3PLUS2);
  }
  // (3 + 0) or (3 + 1) --> H3PLUS0 --> Q = 1
  else if ((mPath->nprimitivesUp() == 3 && 
	    (mPath->nprimitivesDown() == 1 || mPath->nprimitivesDown() == 0)) ||
	   ((mPath->nprimitivesUp() == 1 || mPath->nprimitivesUp() == 0) &&
	    mPath->nprimitivesDown() == 3)){
    mPath->setQuality(H3PLUS0);
  }

}
