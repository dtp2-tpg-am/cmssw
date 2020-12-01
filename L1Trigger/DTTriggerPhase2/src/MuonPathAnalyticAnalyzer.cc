#include "L1Trigger/DTTriggerPhase2/interface/MuonPathAnalyticAnalyzer.h"
#include <cmath>
#include <memory>

using namespace edm;
using namespace std;
using namespace cmsdt;
// ============================================================================
// Constructors and destructor
// ============================================================================
MuonPathAnalyticAnalyzer::MuonPathAnalyticAnalyzer(const ParameterSet &pset, edm::ConsumesCollector &iC)
    : MuonPathAnalyzer(pset, iC),
      bxTolerance_(30),
      minQuality_(LOWQGHOST),
      chiSquareThreshold_(50),
      debug_(pset.getUntrackedParameter<bool>("debug")),
      chi2Th_(pset.getUntrackedParameter<double>("chi2Th")),
      tanPhiTh_(pset.getUntrackedParameter<double>("tanPhiTh")),
      use_LSB_(pset.getUntrackedParameter<bool>("use_LSB")),
      tanPsi_precision_(pset.getUntrackedParameter<double>("tanPsi_precision")),
      x_precision_(pset.getUntrackedParameter<double>("x_precision")) {
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "MuonPathAnalyzer: constructor";

  setChiSquareThreshold(chi2Th_ * 100.);
  fillLAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER();
  std::cout << "Filled" << std::endl;

  //shift
  int rawId;
  shift_filename_ = pset.getParameter<edm::FileInPath>("shift_filename");
  std::ifstream ifin3(shift_filename_.fullPath());
  double shift;
  if (ifin3.fail()) {
    throw cms::Exception("Missing Input File")
        << "MuonPathAnalyticAnalyzer::MuonPathAnalyticAnalyzer() -  Cannot find " << shift_filename_.fullPath();
  }
  while (ifin3.good()) {
    ifin3 >> rawId >> shift;
    shiftinfo_[rawId] = shift;
  }

  chosen_sl_ = pset.getUntrackedParameter<int>("trigger_with_sl");

  if (chosen_sl_ != 1 && chosen_sl_ != 3 && chosen_sl_ != 4) {
    LogDebug("MuonPathAnalyticAnalyzer") << "chosen sl must be 1,3 or 4(both superlayers)";
    assert(chosen_sl_ != 1 && chosen_sl_ != 3 && chosen_sl_ != 4);  //4 means run using the two superlayers
  }

  dtGeomH = iC.esConsumes<DTGeometry, MuonGeometryRecord, edm::Transition::BeginRun>();
}

MuonPathAnalyticAnalyzer::~MuonPathAnalyticAnalyzer() {
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "MuonPathAnalyzer: destructor";
}

// ============================================================================
// Main methods (initialise, run, finish)
// ============================================================================
void MuonPathAnalyticAnalyzer::initialise(const edm::EventSetup &iEventSetup) {
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "MuonPathAnalyticAnalyzer::initialiase";

  const MuonGeometryRecord &geom = iEventSetup.get<MuonGeometryRecord>();
  dtGeo_ = &geom.get(dtGeomH);
}

void MuonPathAnalyticAnalyzer::run(edm::Event &iEvent,
                                const edm::EventSetup &iEventSetup,
                                MuonPathPtrs &muonpaths,
                                std::vector<metaPrimitive> &metaPrimitives) {
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "MuonPathAnalyticAnalyzer: run";

  // fit per SL (need to allow for multiple outputs for a single mpath)
  for (auto &muonpath : muonpaths) {
    analyze(muonpath, metaPrimitives);
  }
}

void MuonPathAnalyticAnalyzer::finish() {
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "MuonPathAnalyzer: finish";
};

constexpr int MuonPathAnalyticAnalyzer::LAYER_ARRANGEMENTS_[NUM_LAYERS][NUM_CELL_COMB] = {
    {0, 1, 2},
    {1, 2, 3},  // Consecutive groups
    {0, 1, 3},
    {0, 2, 3}  // Non-consecutive groups
};

//------------------------------------------------------------------
//--- Métodos privados
//------------------------------------------------------------------

void MuonPathAnalyticAnalyzer::analyze(MuonPathPtr &inMPath, std::vector<metaPrimitive> &metaPrimitives) {
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:analyze \t\t\t\t starts";

  // LOCATE MPATH
  int selected_Id = 0;
  if (inMPath->primitive(0)->tdcTimeStamp() != -1)
    selected_Id = inMPath->primitive(0)->cameraId();
  else if (inMPath->primitive(1)->tdcTimeStamp() != -1)
    selected_Id = inMPath->primitive(1)->cameraId();
  else if (inMPath->primitive(2)->tdcTimeStamp() != -1)
    selected_Id = inMPath->primitive(2)->cameraId();
  else if (inMPath->primitive(3)->tdcTimeStamp() != -1)
    selected_Id = inMPath->primitive(3)->cameraId();

  DTLayerId thisLId(selected_Id);
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "Building up MuonPathSLId from rawId in the Primitive";
  DTSuperLayerId MuonPathSLId(thisLId.wheel(), thisLId.station(), thisLId.sector(), thisLId.superLayer());
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "The MuonPathSLId is" << MuonPathSLId;

  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer")
        << "DTp2:analyze \t\t\t\t In analyze function checking if inMPath->isAnalyzable() " << inMPath->isAnalyzable();

  if (chosen_sl_ < 4 && thisLId.superLayer() != chosen_sl_)
    return;  // avoid running when mpath not in chosen SL (for 1SL fitting)

  auto mPath = std::make_shared<MuonPath>(inMPath);
  mPath->setQuality(NOPATH);
  buildLateralities();
  superlayer_datapath(mPath);
  
  int wi[8], tdc[8], lat[8];
  DTPrimitivePtr Prim0(mPath->primitive(0));
  wi[0] = Prim0->channelId();
  tdc[0] = Prim0->tdcTimeStamp();
  DTPrimitivePtr Prim1(mPath->primitive(1));
  wi[1] = Prim1->channelId();
  tdc[1] = Prim1->tdcTimeStamp();
  DTPrimitivePtr Prim2(mPath->primitive(2));
  wi[2] = Prim2->channelId();
  tdc[2] = Prim2->tdcTimeStamp();
  DTPrimitivePtr Prim3(mPath->primitive(3));
  wi[3] = Prim3->channelId();
  tdc[3] = Prim3->tdcTimeStamp();
  for (int i = 4; i < 8; i++) {
    wi[i] = -1;
    tdc[i] = -1;
    lat[i] = -1;
  }

  DTWireId wireId(MuonPathSLId, 2, 1);

  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:analyze \t\t\t\t checking if it passes the min quality cut "
                                      << mPath->quality() << ">" << minQuality_;
  if (mPath->quality() >= minQuality_) {
    if (debug_)
      LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:analyze \t\t\t\t min quality achievedCalidad: " << mPath->quality();
    for (int i = 0; i <= 3; i++) {
      if (debug_)
        LogDebug("MuonPathAnalyticAnalyzer")
            << "DTp2:analyze \t\t\t\t  Capa: " << mPath->primitive(i)->layerId()
            << " Canal: " << mPath->primitive(i)->channelId() << " TDCTime: " << mPath->primitive(i)->tdcTimeStamp();
    }
    if (debug_)
      LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:analyze \t\t\t\t Starting lateralities loop, totalNumValLateralities: "
                                        << totalNumValLateralities_;

    double best_chi2 = 99999.;
    double chi2_jm_tanPhi = 999;
    double chi2_jm_x = -1;
    double chi2_jm_t0 = -1;
    double chi2_phi = -1;
    double chi2_phiB = -1;
    double chi2_chi2 = -1;
    int chi2_quality = -1;
    int bestLat[8];
    for (int i = 0; i < 8; i++) {
      bestLat[i] = -1;
    }

    for (int i = 0; i < totalNumValLateralities_; i++) {  //here
      if (debug_)
        LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:analyze \t\t\t\t\t laterality #- " << i;
      if (debug_)
        LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:analyze \t\t\t\t\t laterality #- " << i << " checking quality:";
      if (debug_)
        LogDebug("MuonPathAnalyticAnalyzer")
            << "DTp2:analyze \t\t\t\t\t laterality #- " << i << " checking mPath Quality=" << mPath->quality();
      if (debug_)
        LogDebug("MuonPathAnalyticAnalyzer")
            << "DTp2:analyze \t\t\t\t\t laterality #- " << i << " latQuality_[i].val=" << latQuality_[i].valid;
      if (debug_)
        LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:analyze \t\t\t\t\t laterality #- " << i << " before if:";

      if (latQuality_[i].valid and
          (((mPath->quality() == HIGHQ or mPath->quality() == HIGHQGHOST) and latQuality_[i].quality == HIGHQ) or
           ((mPath->quality() == LOWQ or mPath->quality() == LOWQGHOST) and latQuality_[i].quality == LOWQ))) {
        if (debug_)
          LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:analyze \t\t\t\t\t laterality #- " << i << " inside if";
        mPath->setBxTimeValue(latQuality_[i].bxValue);
        if (debug_)
          LogDebug("MuonPathAnalyticAnalyzer")
              << "DTp2:analyze \t\t\t\t\t laterality #- " << i << " settingLateralCombination";
        mPath->setLateralComb(lateralities_[i]);
        if (debug_)
          LogDebug("MuonPathAnalyticAnalyzer")
              << "DTp2:analyze \t\t\t\t\t laterality #- " << i << " done settingLateralCombination";

        // Clonamos el objeto analizado.
        auto mpAux = std::make_shared<MuonPath>(mPath);
        lat[0] = mpAux->lateralComb()[0];
        lat[1] = mpAux->lateralComb()[1];
        lat[2] = mpAux->lateralComb()[2];
        lat[3] = mpAux->lateralComb()[3];

        int wiOk[NUM_LAYERS], tdcOk[NUM_LAYERS], latOk[NUM_LAYERS];
        for (int lay = 0; lay < 4; lay++) {
          if (latQuality_[i].invalidateHitIdx == lay) {
            wiOk[lay] = -1;
            tdcOk[lay] = -1;
            latOk[lay] = -1;
          } else {
            wiOk[lay] = wi[lay];
            tdcOk[lay] = tdc[lay];
            latOk[lay] = lat[lay];
          }
        }

        int idxHitNotValid = latQuality_[i].invalidateHitIdx;
        if (idxHitNotValid >= 0) {
          auto dtpAux = std::make_shared<DTPrimitive>();
          mpAux->setPrimitive(dtpAux, idxHitNotValid);
        }

        if (debug_)
          LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:analyze \t\t\t\t\t  calculating parameters ";
        calculatePathParameters(mpAux);
        /* 
		 * After calculating the parameters, and if it is a 4-hit fit,
		 * if the resultant chi2 is higher than the programmed threshold, 
		 * the mpath is eliminated and we go to the next element
		 */
        if ((mpAux->quality() == HIGHQ or mpAux->quality() == HIGHQGHOST) &&
            mpAux->chiSquare() > chiSquareThreshold_) {  //check this if!!!
          if (debug_)
            LogDebug("MuonPathAnalyticAnalyzer")
                << "DTp2:analyze \t\t\t\t\t  HIGHQ or HIGHQGHOST but min chi2 or Q test not satisfied ";
        } else {
          if (debug_)
            LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:analyze \t\t\t\t\t  inside else, returning values: ";
          if (debug_)
            LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:analyze \t\t\t\t\t  BX Time = " << mpAux->bxTimeValue();
          if (debug_)
            LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:analyze \t\t\t\t\t  BX Id   = " << mpAux->bxNumId();
          if (debug_)
            LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:analyze \t\t\t\t\t  XCoor   = " << mpAux->horizPos();
          if (debug_)
            LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:analyze \t\t\t\t\t  tan(Phi)= " << mpAux->tanPhi();
          if (debug_)
            LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:analyze \t\t\t\t\t  chi2= " << mpAux->chiSquare();
          if (debug_)
            LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:analyze \t\t\t\t\t  lateralities = "
                                              << " " << mpAux->lateralComb()[0] << " " << mpAux->lateralComb()[1] << " "
                                              << mpAux->lateralComb()[2] << " " << mpAux->lateralComb()[3];

          DTChamberId ChId(MuonPathSLId.wheel(), MuonPathSLId.station(), MuonPathSLId.sector());

          double jm_tanPhi = -1. * mpAux->tanPhi();  //testing with this line
          if (use_LSB_)
            jm_tanPhi = floor(jm_tanPhi / tanPsi_precision_) * tanPsi_precision_;
          double jm_x =
              (((double)mpAux->horizPos()) / 10.) + x_precision_ * (round(shiftinfo_[wireId.rawId()] / x_precision_));
          if (use_LSB_)
            jm_x = ((double)round(((double)jm_x) / x_precision_)) * x_precision_;
          //changing to chamber frame or reference:
          double jm_t0 = mpAux->bxTimeValue();
          int quality = mpAux->quality();

          //computing phi and phiB
          double z = 0;
          double z1 = Z_POS_SL;
          double z3 = -1. * z1;
          if (ChId.station() == 3 or ChId.station() == 4) {
            z1 = z1 + Z_SHIFT_MB4;
            z3 = z3 + Z_SHIFT_MB4;
          } else if (MuonPathSLId.superLayer() == 1)
            z = z1;
          else if (MuonPathSLId.superLayer() == 3)
            z = z3;

          GlobalPoint jm_x_cmssw_global = dtGeo_->chamber(ChId)->toGlobal(LocalPoint(jm_x, 0., z));
          int thisec = MuonPathSLId.sector();
          if (thisec == 13)
            thisec = 4;
          if (thisec == 14)
            thisec = 10;
          double phi = jm_x_cmssw_global.phi() - PHI_CONV * (thisec - 1);
          double psi = atan(jm_tanPhi);
          double phiB = hasPosRF(MuonPathSLId.wheel(), MuonPathSLId.sector()) ? psi - phi : -psi - phi;
          double chi2 = mpAux->chiSquare() * 0.01;  //in cmssw we need cm, 1 cm^2 = 100 mm^2

          if (debug_)
            LogDebug("MuonPathAnalyticAnalyzer")
                << "DTp2:analyze \t\t\t\t\t\t\t\t  pushing back metaPrimitive at x=" << jm_x << " tanPhi:" << jm_tanPhi
                << " t0:" << jm_t0;

          if (mpAux->quality() == HIGHQ or
              mpAux->quality() == HIGHQGHOST) {  //keep only the values with the best chi2 among lateralities
            if ((chi2 < best_chi2) && (std::abs(jm_tanPhi) <= tanPhiTh_)) {
              chi2_jm_tanPhi = jm_tanPhi;
              chi2_jm_x = (mpAux->horizPos() / 10.) + shiftinfo_[wireId.rawId()];
              chi2_jm_t0 = mpAux->bxTimeValue();
              chi2_phi = phi;
              chi2_phiB = phiB;
              chi2_chi2 = chi2;
              best_chi2 = chi2;
              chi2_quality = mpAux->quality();
              for (int i = 0; i < 4; i++) {
                bestLat[i] = lat[i];
              }
            }
          } else if (std::abs(jm_tanPhi) <=
                     tanPhiTh_) {  //write the metaprimitive in case no HIGHQ or HIGHQGHOST and tanPhi range
            if (debug_)
              LogDebug("MuonPathAnalyticAnalyzer")
                  << "DTp2:analyze \t\t\t\t\t\t\t\t  pushing back metaprimitive no HIGHQ or HIGHQGHOST";
            metaPrimitives.emplace_back(metaPrimitive({MuonPathSLId.rawId(),
                                                       jm_t0,
                                                       jm_x,
                                                       jm_tanPhi,
                                                       phi,
                                                       phiB,
                                                       chi2,
                                                       quality,
                                                       wiOk[0],
                                                       tdcOk[0],
                                                       latOk[0],
                                                       wiOk[1],
                                                       tdcOk[1],
                                                       latOk[1],
                                                       wiOk[2],
                                                       tdcOk[2],
                                                       latOk[2],
                                                       wiOk[3],
                                                       tdcOk[3],
                                                       latOk[3],
                                                       wi[4],
                                                       tdc[4],
                                                       lat[4],
                                                       wi[5],
                                                       tdc[5],
                                                       lat[5],
                                                       wi[6],
                                                       tdc[6],
                                                       lat[6],
                                                       wi[7],
                                                       tdc[7],
                                                       lat[7],
                                                       -1}));
            if (debug_)
              LogDebug("MuonPathAnalyticAnalyzer")
                  << "DTp2:analyze \t\t\t\t\t\t\t\t  done pushing back metaprimitive no HIGHQ or HIGHQGHOST";
          }
        }
      } else {
        if (debug_)
          LogDebug("MuonPathAnalyticAnalyzer")
              << "DTp2:analyze \t\t\t\t\t\t\t\t  latQuality_[i].valid and (((mPath->quality()==HIGHQ or "
                 "mPath->quality()==HIGHQGHOST) and latQuality_[i].quality==HIGHQ) or  ((mPath->quality() "
                 "== LOWQ or mPath->quality()==LOWQGHOST) and latQuality_[i].quality==LOWQ)) not passed";
      }
    }
    if (chi2_jm_tanPhi != 999 and std::abs(chi2_jm_tanPhi) < tanPhiTh_) {  //
      if (debug_)
        LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:analyze \t\t\t\t\t\t\t\t  pushing back best chi2 metaPrimitive";
      metaPrimitives.emplace_back(metaPrimitive({MuonPathSLId.rawId(),
                                                 chi2_jm_t0,
                                                 chi2_jm_x,
                                                 chi2_jm_tanPhi,
                                                 chi2_phi,
                                                 chi2_phiB,
                                                 chi2_chi2,
                                                 chi2_quality,
                                                 wi[0],
                                                 tdc[0],
                                                 bestLat[0],
                                                 wi[1],
                                                 tdc[1],
                                                 bestLat[1],
                                                 wi[2],
                                                 tdc[2],
                                                 bestLat[2],
                                                 wi[3],
                                                 tdc[3],
                                                 bestLat[3],
                                                 wi[4],
                                                 tdc[4],
                                                 bestLat[4],
                                                 wi[5],
                                                 tdc[5],
                                                 bestLat[5],
                                                 wi[6],
                                                 tdc[6],
                                                 bestLat[6],
                                                 wi[7],
                                                 tdc[7],
                                                 bestLat[7],
                                                 -1}));
    }
  }
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:analyze \t\t\t\t finishes";
}

void MuonPathAnalyticAnalyzer::setCellLayout(const int layout[NUM_LAYERS]) {
  memcpy(cellLayout_, layout, 4 * sizeof(int));

  buildLateralities();
}

/**
 * For a given 4-cell combination (one per layer), all the possible lateralities 
 * combinations that are compatible with a straight line are generated. 
 */
void MuonPathAnalyticAnalyzer::buildLateralities(void) {
  LATERAL_CASES(*validCase)[NUM_LAYERS], sideComb[NUM_LAYERS];

  totalNumValLateralities_ = 0;
  /* We generate all the possible lateralities combination for a given group 
     of cells */
  for (int lowLay = LEFT; lowLay <= RIGHT; lowLay++)
    for (int midLowLay = LEFT; midLowLay <= RIGHT; midLowLay++)
      for (int midHigLay = LEFT; midHigLay <= RIGHT; midHigLay++)
        for (int higLay = LEFT; higLay <= RIGHT; higLay++) {
          sideComb[0] = static_cast<LATERAL_CASES>(lowLay);
          sideComb[1] = static_cast<LATERAL_CASES>(midLowLay);
          sideComb[2] = static_cast<LATERAL_CASES>(midHigLay);
          sideComb[3] = static_cast<LATERAL_CASES>(higLay);

          /* If a laterality combination is valid, we store it  */
          if (isStraightPath(sideComb)) {
            validCase = lateralities_ + totalNumValLateralities_;
            memcpy(validCase, sideComb, 4 * sizeof(LATERAL_CASES));

            latQuality_[totalNumValLateralities_].valid = false;
            latQuality_[totalNumValLateralities_].bxValue = 0;
            latQuality_[totalNumValLateralities_].quality = NOPATH;
            latQuality_[totalNumValLateralities_].invalidateHitIdx = -1;

            totalNumValLateralities_++;
          }
        }
}

/**
 * This method checks whether a given combination conform a straight line or not
 */
bool MuonPathAnalyticAnalyzer::isStraightPath(LATERAL_CASES sideComb[NUM_LAYERS]) {
  return true;  //trying with all lateralities to be confirmed

  int i, ajustedLayout[NUM_LAYERS], pairDiff[3], desfase[3];

  for (i = 0; i <= 3; i++)
    ajustedLayout[i] = cellLayout_[i] + sideComb[i];
  for (i = 0; i <= 2; i++)
    pairDiff[i] = ajustedLayout[i + 1] - ajustedLayout[i];
  for (i = 0; i <= 1; i++)
    desfase[i] = abs(pairDiff[i + 1] - pairDiff[i]);
  desfase[2] = abs(pairDiff[2] - pairDiff[0]);
  bool resultado = (desfase[0] > 1 or desfase[1] > 1 or desfase[2] > 1);

  return (!resultado);
}

int MuonPathAnalyticAnalyzer::eqMainTerm(LATERAL_CASES sideComb[2], int layerIdx[2], MuonPathPtr &mPath, int bxValue) {
  int eqTerm = 0, coefs[2];

  lateralCoeficients(sideComb, coefs);

  if (!use_LSB_)
    eqTerm = coefs[0] * (mPath->primitive(layerIdx[0])->tdcTimeStampNoOffset() - bxValue) +
             coefs[1] * (mPath->primitive(layerIdx[1])->tdcTimeStampNoOffset() - bxValue);
  else
    eqTerm = coefs[0] * floor((DRIFT_SPEED / (10 * x_precision_)) *
                              (mPath->primitive(layerIdx[0])->tdcTimeStampNoOffset() - bxValue)) +
             coefs[1] * floor((DRIFT_SPEED / (10 * x_precision_)) *
                              (mPath->primitive(layerIdx[1])->tdcTimeStampNoOffset() - bxValue));

  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:\t\t\t\t\t EQTerm(Main): " << eqTerm;

  return (eqTerm);
}

void MuonPathAnalyticAnalyzer::lateralCoeficients(LATERAL_CASES sideComb[2], int *coefs) {
  if ((sideComb[0] == LEFT) && (sideComb[1] == LEFT)) {
    *(coefs) = +1;
    *(coefs + 1) = -1;
  } else if ((sideComb[0] == LEFT) && (sideComb[1] == RIGHT)) {
    *(coefs) = +1;
    *(coefs + 1) = +1;
  } else if ((sideComb[0] == RIGHT) && (sideComb[1] == LEFT)) {
    *(coefs) = -1;
    *(coefs + 1) = -1;
  } else if ((sideComb[0] == RIGHT) && (sideComb[1] == RIGHT)) {
    *(coefs) = -1;
    *(coefs + 1) = +1;
  }
}

/** Calculate the parameters of the detected trayectories */
void MuonPathAnalyticAnalyzer::calculatePathParameters(MuonPathPtr &mPath) {
  // The order is important.
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer")
        << "DTp2:calculatePathParameters \t\t\t\t\t\t  calculating calcCellDriftAndXcoor(mPath) ";
  calcCellDriftAndXcoor(mPath);
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:calculatePathParameters \t\t\t\t\t\t  checking mPath->quality() "
                                      << mPath->quality();
  if (mPath->quality() == HIGHQ or mPath->quality() == HIGHQGHOST) {
    if (debug_)
      LogDebug("MuonPathAnalyticAnalyzer")
          << "DTp2:calculatePathParameters \t\t\t\t\t\t\t  Quality test passed, now calcTanPhiXPosChamber4Hits(mPath) ";
    calcTanPhiXPosChamber4Hits(mPath);
  } else {
    if (debug_)
      LogDebug("MuonPathAnalyticAnalyzer")
          << "DTp2:calculatePathParameters \t\t\t\t\t\t\t  Quality test NOT passed calcTanPhiXPosChamber3Hits(mPath) ";
    calcTanPhiXPosChamber3Hits(mPath);
  }

  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:calculatePathParameters \t\t\t\t\t\t calcChiSquare(mPath) ";
  calcChiSquare(mPath);
}

void MuonPathAnalyticAnalyzer::calcTanPhiXPosChamber(MuonPathPtr &mPath) {
  int layerIdx[2];
  /*
      To calculate path's angle are only necessary two valid primitives.
      This method should be called only when a 'MuonPath' is determined as valid,
      so, at least, three of its primitives must have a valid time.
      With this two comparitions (which can be implemented easily as multiplexors
      in the FPGA) this method ensures to catch two of those valid primitives to
      evaluate the angle.

      The first one is below the middle line of the superlayer, while the other
      one is above this line
    */
  if (mPath->primitive(0)->isValidTime())
    layerIdx[0] = 0;
  else
    layerIdx[0] = 1;

  if (mPath->primitive(3)->isValidTime())
    layerIdx[1] = 3;
  else
    layerIdx[1] = 2;

  /* We identify along which cells' sides the muon travels */
  LATERAL_CASES sideComb[2];
  sideComb[0] = (mPath->lateralComb())[layerIdx[0]];
  sideComb[1] = (mPath->lateralComb())[layerIdx[1]];

  /* Horizontal gap between cells in cell's semi-length units */
  int dHoriz = (mPath->cellLayout())[layerIdx[1]] - (mPath->cellLayout())[layerIdx[0]];

  /* Vertical gap between cells in cell's height units */
  int dVert = layerIdx[1] - layerIdx[0];

  /*-----------------------------------------------------------------*/
  /*--------------------- Phi angle calculation ---------------------*/
  /*-----------------------------------------------------------------*/
  float num = CELL_SEMILENGTH * dHoriz + DRIFT_SPEED * eqMainTerm(sideComb, layerIdx, mPath, mPath->bxTimeValue());

  float denom = CELL_HEIGHT * dVert;
  float tanPhi = num / denom;

  mPath->setTanPhi(tanPhi);

  /*-----------------------------------------------------------------*/
  /*----------------- Horizontal coord. calculation -----------------*/
  /*-----------------------------------------------------------------*/

  /*
      Using known coordinates, relative to superlayer axis reference, (left most
      superlayer side, and middle line between 2nd and 3rd layers), calculating
      horizontal coordinate implies using a basic line equation:
      (y - y0) = (x - x0) * cotg(Phi)
      This horizontal coordinate can be obtained setting y = 0 on last equation,
      and also setting y0 and x0 with the values of a known muon's path cell
      position hit.
      It's enough to use the lower cell (layerIdx[0]) coordinates. So:
      xC = x0 - y0 * tan(Phi)
    */
  float lowerXPHorizPos = mPath->xCoorCell(layerIdx[0]);

  float lowerXPVertPos = 0;  // This is only the absolute value distance.
  if (layerIdx[0] == 0)
    lowerXPVertPos = CELL_HEIGHT + CELL_SEMIHEIGHT;
  else
    lowerXPVertPos = CELL_SEMIHEIGHT;

  mPath->setHorizPos(lowerXPHorizPos + lowerXPVertPos * tanPhi);
}

/**
 * Coordinate and angle calculations for a 4 HITS cases
 */
void MuonPathAnalyticAnalyzer::calcTanPhiXPosChamber4Hits(MuonPathPtr &mPath) {
  int x_prec_inv = (int)(1. / (10. * x_precision_));
  int numberOfBits = (int)(round(std::log(x_prec_inv) / std::log(2.)));
  int numerator = 3 * (int)round(mPath->xCoorCell(3) / (10 * x_precision_)) +
                  (int)round(mPath->xCoorCell(2) / (10 * x_precision_)) -
                  (int)round(mPath->xCoorCell(1) / (10 * x_precision_)) -
                  3 * (int)round(mPath->xCoorCell(0) / (10 * x_precision_));
  int CELL_HEIGHT_JM = pow(2, 15) / ((int)(10 * CELL_HEIGHT));
  int tanPhi_x4096 = (numerator * CELL_HEIGHT_JM) >> (3 + numberOfBits);
  mPath->setTanPhi(tanPhi_x4096 * tanPsi_precision_);

  float XPos = (mPath->xCoorCell(0) + mPath->xCoorCell(1) + mPath->xCoorCell(2) + mPath->xCoorCell(3)) / 4;
  mPath->setHorizPos(floor(XPos / (10 * x_precision_)) * 10 * x_precision_);
}

/**
 *  3 HITS cases
 */
void MuonPathAnalyticAnalyzer::calcTanPhiXPosChamber3Hits(MuonPathPtr &mPath) {
  int layerIdx[2];
  int x_prec_inv = (int)(1. / (10. * x_precision_));
  int numberOfBits = (int)(round(std::log(x_prec_inv) / std::log(2.)));

  if (mPath->primitive(0)->isValidTime())
    layerIdx[0] = 0;
  else
    layerIdx[0] = 1;

  if (mPath->primitive(3)->isValidTime())
    layerIdx[1] = 3;
  else
    layerIdx[1] = 2;

  /*-----------------------------------------------------------------*/
  /*--------------------- Phi angle calculation ---------------------*/
  /*-----------------------------------------------------------------*/

  int tan_division_denominator_bits = 16;

  int num =
      ((int)((int)(x_prec_inv * mPath->xCoorCell(layerIdx[1])) - (int)(x_prec_inv * mPath->xCoorCell(layerIdx[0])))
       << (12 - numberOfBits));
  int denominator = (layerIdx[1] - layerIdx[0]) * CELL_HEIGHT;
  int denominator_inv = ((int)(0.5 + pow(2, tan_division_denominator_bits) / float(denominator)));

  float tanPhi = ((num * denominator_inv) >> tan_division_denominator_bits) / ((1. / tanPsi_precision_));

  mPath->setTanPhi(tanPhi);

  /*-----------------------------------------------------------------*/
  /*----------------- Horizontal coord. calculation -----------------*/
  /*-----------------------------------------------------------------*/
  float XPos = 0;
  if (mPath->primitive(0)->isValidTime() and mPath->primitive(3)->isValidTime())
    XPos = (mPath->xCoorCell(0) + mPath->xCoorCell(3)) / 2;
  else
    XPos = (mPath->xCoorCell(1) + mPath->xCoorCell(2)) / 2;

  mPath->setHorizPos(floor(XPos / (10 * x_precision_)) * 10 * x_precision_);
}

/**
 * Calculate the drift distances of each wire and the horizontal position 
 */
void MuonPathAnalyticAnalyzer::calcCellDriftAndXcoor(MuonPathPtr &mPath) {
  long int drift_speed_new = 889;
  long int drift_dist_um_x4;
  long int wireHorizPos_x4;
  long int pos_mm_x4;
  int x_prec_inv = (int)(1. / (10. * x_precision_));

  for (int i = 0; i <= 3; i++)
    if (mPath->primitive(i)->isValidTime()) {
      drift_dist_um_x4 =
          drift_speed_new * ((long int)mPath->primitive(i)->tdcTimeStampNoOffset() - (long int)mPath->bxTimeValue());
      wireHorizPos_x4 = (long)(mPath->primitive(i)->wireHorizPos() * x_prec_inv);

      if ((mPath->lateralComb())[i] == LEFT)
        pos_mm_x4 = wireHorizPos_x4 - (drift_dist_um_x4 >> 10);
      else
        pos_mm_x4 = wireHorizPos_x4 + (drift_dist_um_x4 >> 10);

      mPath->setXCoorCell(pos_mm_x4 * (10 * x_precision_), i);
      mPath->setDriftDistance(((float)(drift_dist_um_x4 >> 10)) * (10 * x_precision_), i);
    }
}

/**
 * Calculate the quality estimator of each trayectory.
 */
void MuonPathAnalyticAnalyzer::calcChiSquare(MuonPathPtr &mPath) {
  int x_prec_inv = (int)(1. / (10. * x_precision_));
  int numberOfBits = (int)(round(std::log(x_prec_inv) / std::log(2.)));
  long int Z_FACTOR[NUM_LAYERS] = {-6, -2, 2, 6};
  for (int i = 0; i < 4; i++) {
    Z_FACTOR[i] = Z_FACTOR[i] * (long int)CELL_HEIGHT;
  }
  long int sum_A = 0, sum_B = 0;
  long int chi2_mm2_x1024 = 0;
  for (int i = 0; i < 4; i++) {
    if (mPath->primitive(i)->isValidTime()) {
      sum_A = (((int)(mPath->xCoorCell(i) / (10 * x_precision_))) - ((int)(mPath->horizPos() / (10 * x_precision_))))
              << (14 - numberOfBits);
      sum_B = Z_FACTOR[i] * ((int)(mPath->tanPhi() / tanPsi_precision_));
      chi2_mm2_x1024 += (sum_A - sum_B) * (sum_A - sum_B);
    }
  }
  chi2_mm2_x1024 = chi2_mm2_x1024 >> 18;

  mPath->setChiSquare(((double)chi2_mm2_x1024 / 1024.));
}

int MuonPathAnalyticAnalyzer::omittedHit(int idx) {
  switch (idx) {
    case 0:
      return 3;
    case 1:
      return 0;
    case 2:
      return 2;
    case 3:
      return 1;
  }

  return -1;
}

void MuonPathAnalyticAnalyzer::superlayer_datapath(MuonPathPtr &mPath) {
  int validCells = 0;
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") <<  ">>>>DEBUG: SUPERLAYER_DATAPATH " << endl;

  for (int j = 0; j < 4; j++)
    if (mPath->primitive(j)->isValidTime())
      validCells++;

  if (validCells < 3) return;
  int wires[4], t0s[4], valids[4], cell_horiz_layout[4];
  for (int j = 0; j < 4; j++) {
    if (mPath->primitive(j)->isValidTime()){
      wires[j] = mPath->primitive(j)->channelId();
      t0s[j] = mPath->primitive(j)->tdcTimeStamp();
      valids[j] = 1; 
    } else {
      wires[j] = -1;
      t0s[j] = -1;
      valids[j] = 0; 
    }
  }

  if (wires[0] < 0) wires[0] = wires[1]; 
  if (wires[1] < 0) wires[1] = wires[0]; 
  if (wires[2] < 0) wires[2] = wires[1] - 1; 
  if (wires[3] < 0) wires[3] = wires[2]; 

  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << ">>>>DEBUG: SUPERLAYER_DATAPATH" << endl
      << "wires[0]:" << wires[0] << " wires[1]:" << wires[1] << " wires[2]:" << wires[2] <<" wires[3]:" << wires[3] << endl 
      << "t0[0]:" << t0s[0] << " t0[1]:" << t0s[1] << " t0[2]:" << t0s[2] <<" t0[3]:" << t0s[3] << endl 
      << "valid[0]:" << valids[0] << " valid[1]:" << valids[1] << " valid[2]:" << valids[2] <<" valid[3]:" << valids[3] << endl 
      << ">>>>>>>>>>>" << endl; 


  for (int j = 0; j < 4; j++) {
    cell_horiz_layout[j] = 2*(wires[j] - wires[0]);
    if (j == 1 || j == 3) cell_horiz_layout[j]--; 
  } 

  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "Starting segment composer with cell_horiz_layout " << cell_horiz_layout[0] << cell_horiz_layout[1] <<cell_horiz_layout[2] <<cell_horiz_layout[3] << endl;
  segment_composer(valids, wires, t0s, cell_horiz_layout, mPath); 

}

void MuonPathAnalyticAnalyzer::segment_composer(int valids[NUM_LAYERS], int wires[NUM_LAYERS],
    int t0s[NUM_LAYERS], int cell_horiz_layout[NUM_LAYERS], MuonPathPtr &mPath) {

  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << ">>>>DEBUG: Inside segment composer" << endl
      << "wires[0]:" << wires[0] << " wires[1]:" << wires[1] << " wires[2]:" << wires[2] <<" wires[3]:" << wires[3]<< endl 
      << "t0[0]:" << t0s[0] << " t0[1]:" << t0s[1] << " t0[2]:" << t0s[2] <<" t0[3]:" << t0s[3]<< endl 
      << "valid[0]:" << valids[0] << " valid[1]:" << valids[1] << " valid[2]:" << valids[2] <<" valid[3]:" << valids[3]<< endl 
      << ">>>>>>>>>>>" << endl; 
  
  bool notAllValid = false;
  int omitted_layer = -1;  
  for (int i = 0; i < 4; i++){
    if (valids[i] == 0) { 
      notAllValid = true;
      omitted_layer = i; 
    }
  }

  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "notAllValid=" << notAllValid << endl; 

  std::vector <latcomb_consts> my_latcomb_consts = {};
  for (auto & elem : LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER) 
    if (elem.first.valid[0] == valids[0] 
        && elem.first.valid[1] == valids[1]
        && elem.first.valid[2] == valids[2]
        && elem.first.valid[3] == valids[3]
        && elem.first.cell_horiz_layout[0] == cell_horiz_layout[0]
        && elem.first.cell_horiz_layout[1] == cell_horiz_layout[1]
        && elem.first.cell_horiz_layout[2] == cell_horiz_layout[2]
        && elem.first.cell_horiz_layout[3] == cell_horiz_layout[3]) {

      for (auto & ind_latcomb_consts : elem.second)
        my_latcomb_consts.push_back(ind_latcomb_consts);
    } 
  for (auto & ind_latcomb_consts : my_latcomb_consts ){
    if (debug_)
      LogDebug("MuonPathAnalyticAnalyzer") << "Computing t0 with latcomb_consts"  << ind_latcomb_consts.latcomb <<","<<  ind_latcomb_consts.numer_const << ",{" << ind_latcomb_consts.numer_coeff[0] << "," <<ind_latcomb_consts.numer_coeff[1] << "," <<ind_latcomb_consts.numer_coeff[2] << "," <<ind_latcomb_consts.numer_coeff[3] << "}," << ind_latcomb_consts.denom << endl; 
    int prim_t0 = validate_4hits_one_latcomb_reduced(valids, t0s, cell_horiz_layout, ind_latcomb_consts);
   
    int lat = changeLatToEmu (ind_latcomb_consts.latcomb); 
    if (prim_t0 != 0) {
      if (debug_)
        LogDebug("MuonPathAnalyticAnalyzer") << "Will fill Alvaro's laterality " <<  ind_latcomb_consts.latcomb << "= Emu lat " << lat << endl; 
      latQuality_[lat].valid = true; 
      latQuality_[lat].bxValue = prim_t0; 
      
      if (notAllValid) latQuality_[lat].quality = LOWQ; 
      else latQuality_[lat].quality = HIGHQ;
      
      if (debug_)
        LogDebug("MuonPathAnalyticAnalyzer") << "Obtained a combination with quality " << latQuality_[lat].quality << endl; 
      if (debug_)
        LogDebug("MuonPathAnalyticAnalyzer") << "My mpath quality before update is " << mPath->quality() << endl;  
      if (mPath->quality() < latQuality_[lat].quality) mPath->setQuality(latQuality_[lat].quality);
      if (debug_)
        LogDebug("MuonPathAnalyticAnalyzer") << "My mpath quality after update is " << mPath->quality() << endl;  

      latQuality_[lat].invalidateHitIdx = omitted_layer; 
    }
  }
 
}

int MuonPathAnalyticAnalyzer::validate_4hits_one_latcomb_reduced(int valids[cmsdt::NUM_LAYERS],
    int t0s[cmsdt::NUM_LAYERS], int cell_horiz_layout[cmsdt::NUM_LAYERS], latcomb_consts my_latcomb_consts) {

  int t_perfect, t_tucked;  
  bool notAllValid = false; 
  for (int i = 0; i < cmsdt::NUM_LAYERS; i++){
    if (valids[i] == 0) notAllValid = true;
  }

  if (my_latcomb_consts.denom == 0) {
    t_perfect = 0; 
    t_tucked = 0;
  } else {
    t_perfect = my_latcomb_consts.numer_const;
    for (int i = 0; i < cmsdt::NUM_LAYERS; i++){ 
      t_perfect = t_perfect + my_latcomb_consts.numer_coeff[i] * t0s[i];
    }
    t_perfect = t_perfect / my_latcomb_consts.denom; 
    t_tucked = tuck_t0(t_perfect, t0s, valids);
    if (debug_)
      LogDebug("MuonPathAnalyticAnalyzer") << "Perfect time: " << t_perfect << ", Tucked time: " << t_tucked << ", Equal = "<< (t_tucked == t_perfect) << endl; 
    if (t_tucked != t_perfect && notAllValid) t_tucked = 0;
    if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "Tucked time after comparison: " << t_tucked << endl; 
  }
  return int(t_tucked);
}

int MuonPathAnalyticAnalyzer::tuck_t0(int t, int t0s[4], int valids[4]){
    
  double MAX_DRIFT_TIME_F = 386.75; 
  int min_t = -1, max_t = 99999;  
  for (int i = 0; i < 4; i++){
    if (valids[i] == 0) continue; 
    if (t0s[i] < max_t) max_t = t0s[i];
    int tmp_min_t =  t0s[i] - round(MAX_DRIFT_TIME_F);
    if (tmp_min_t>min_t) min_t = tmp_min_t; 
  }
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "min_t:" << min_t << " t:" << t << " max_t:" << max_t << endl;
  if (min_t > max_t || t < min_t-bxTolerance_ || t>max_t+bxTolerance_) return 0;
  else {
    std::vector <int> times = {min_t, t, max_t}; 
    std::sort (times.begin(), times.end()); 
    return times.at(1);
  } 

}

int MuonPathAnalyticAnalyzer::changeLatToEmu (int lat) {
  std::vector <int> binaryNum = {};
  while (lat > 1) {
    binaryNum.push_back(lat % 2);
    lat = lat / 2;
  }
  binaryNum.push_back(lat);
  while (binaryNum.size() < 4) binaryNum.push_back(0);
 

  int emuLat = 0, factor = 1; 
  for (int i = binaryNum.size()-1; i>=0; i--){
     int digit = binaryNum.at(i);
     emuLat = emuLat + digit*factor; 
     factor = 2*factor; 
  }
  return emuLat; 
}

void MuonPathAnalyticAnalyzer::fillLAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER(){
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,0,1},{1,1,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({2,-3094.0,{4,7,2,1},14}),
  latcomb_consts({6,-773.5,{1,1,1,1},4}),
  latcomb_consts({14,-2320.5,{3,4,1,-2},6}),
  latcomb_consts({0,0,{0,0,0,0},0}),
  latcomb_consts({0,0,{0,0,0,0},0}),
  latcomb_consts({0,0,{0,0,0,0},0}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,-2,-1},{1,1,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({4,-3094.0,{1,2,7,4},14}),
  latcomb_consts({6,-773.5,{1,1,1,1},4}),
  latcomb_consts({7,-2320.5,{-2,1,4,3},6}),
  latcomb_consts({0,0,{0,0,0,0},0}),
  latcomb_consts({0,0,{0,0,0,0},0}),
  latcomb_consts({0,0,{0,0,0,0},0}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,0,1},{1,1,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({4,-2320.5,{1,2,7,4},14}),
  latcomb_consts({12,-1547.0,{-1,3,3,-1},4}),
  latcomb_consts({5,-1547.0,{1,3,3,1},8}),
  latcomb_consts({13,-2320.5,{4,7,2,1},14}),
  latcomb_consts({0,0,{0,0,0,0},0}),
  latcomb_consts({0,0,{0,0,0,0},0}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,0,-1},{1,1,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({1,-2320.5,{3,4,1,-2},6}),
  latcomb_consts({9,-773.5,{1,1,1,1},4}),
  latcomb_consts({13,-3094.0,{4,7,2,1},14}),
  latcomb_consts({0,0,{0,0,0,0},0}),
  latcomb_consts({0,0,{0,0,0,0},0}),
  latcomb_consts({0,0,{0,0,0,0},0}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,0,-1},{1,1,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({2,-2320.5,{4,7,2,1},14}),
  latcomb_consts({10,-1547.0,{1,3,3,1},8}),
  latcomb_consts({3,-1547.0,{-1,3,3,-1},4}),
  latcomb_consts({11,-2320.5,{1,2,7,4},14}),
  latcomb_consts({0,0,{0,0,0,0},0}),
  latcomb_consts({0,0,{0,0,0,0},0}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,2,1},{1,1,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({8,-2320.5,{-2,1,4,3},6}),
  latcomb_consts({9,-773.5,{1,1,1,1},4}),
  latcomb_consts({11,-3094.0,{1,2,7,4},14}),
  latcomb_consts({0,0,{0,0,0,0},0}),
  latcomb_consts({0,0,{0,0,0,0},0}),
  latcomb_consts({0,0,{0,0,0,0},0}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,-2,-3},{1,1,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({8,0.0,{-2,1,4,3},6}),
  latcomb_consts({12,0.0,{-1,3,3,-1},4}),
  latcomb_consts({14,0.0,{3,4,1,-2},6}),
  latcomb_consts({1,0.0,{3,4,1,-2},6}),
  latcomb_consts({3,0.0,{-1,3,3,-1},4}),
  latcomb_consts({7,0.0,{-2,1,4,3},6}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,2,3},{1,1,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({8,0.0,{-2,1,4,3},6}),
  latcomb_consts({12,0.0,{-1,3,3,-1},4}),
  latcomb_consts({14,0.0,{3,4,1,-2},6}),
  latcomb_consts({1,0.0,{3,4,1,-2},6}),
  latcomb_consts({3,0.0,{-1,3,3,-1},4}),
  latcomb_consts({7,0.0,{-2,1,4,3},6}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,0,1},{0,1,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({2,-0.0,{0,1,2,-1},2}),
  latcomb_consts({4,0.0,{0,1,2,1},4}),
  latcomb_consts({6,0.0,{0,-1,2,1},2}),
  latcomb_consts({8,-0.0,{0,-1,2,1},2}),
  latcomb_consts({10,-0.0,{0,1,2,1},4}),
  latcomb_consts({12,0.0,{0,1,2,-1},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,-2,-1},{0,1,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({2,773.5,{0,1,2,-1},2}),
  latcomb_consts({4,-773.5,{0,1,2,1},4}),
  latcomb_consts({6,-773.5,{0,-1,2,1},2}),
  latcomb_consts({8,773.5,{0,-1,2,1},2}),
  latcomb_consts({10,773.5,{0,1,2,1},4}),
  latcomb_consts({12,-773.5,{0,1,2,-1},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,0,1},{0,1,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({2,773.5,{0,1,2,-1},2}),
  latcomb_consts({4,-773.5,{0,1,2,1},4}),
  latcomb_consts({6,-773.5,{0,-1,2,1},2}),
  latcomb_consts({8,773.5,{0,-1,2,1},2}),
  latcomb_consts({10,773.5,{0,1,2,1},4}),
  latcomb_consts({12,-773.5,{0,1,2,-1},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,0,-1},{0,1,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({2,-0.0,{0,1,2,-1},2}),
  latcomb_consts({4,0.0,{0,1,2,1},4}),
  latcomb_consts({6,0.0,{0,-1,2,1},2}),
  latcomb_consts({8,-0.0,{0,-1,2,1},2}),
  latcomb_consts({10,-0.0,{0,1,2,1},4}),
  latcomb_consts({12,0.0,{0,1,2,-1},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,0,-1},{0,1,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({2,-773.5,{0,1,2,-1},2}),
  latcomb_consts({4,773.5,{0,1,2,1},4}),
  latcomb_consts({6,773.5,{0,-1,2,1},2}),
  latcomb_consts({8,-773.5,{0,-1,2,1},2}),
  latcomb_consts({10,-773.5,{0,1,2,1},4}),
  latcomb_consts({12,773.5,{0,1,2,-1},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,2,1},{0,1,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({2,-773.5,{0,1,2,-1},2}),
  latcomb_consts({4,773.5,{0,1,2,1},4}),
  latcomb_consts({6,773.5,{0,-1,2,1},2}),
  latcomb_consts({8,-773.5,{0,-1,2,1},2}),
  latcomb_consts({10,-773.5,{0,1,2,1},4}),
  latcomb_consts({12,773.5,{0,1,2,-1},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,-2,-3},{0,1,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({2,-0.0,{0,1,2,-1},2}),
  latcomb_consts({4,0.0,{0,1,2,1},4}),
  latcomb_consts({6,0.0,{0,-1,2,1},2}),
  latcomb_consts({8,-0.0,{0,-1,2,1},2}),
  latcomb_consts({10,-0.0,{0,1,2,1},4}),
  latcomb_consts({12,0.0,{0,1,2,-1},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,2,3},{0,1,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({2,-0.0,{0,1,2,-1},2}),
  latcomb_consts({4,0.0,{0,1,2,1},4}),
  latcomb_consts({6,0.0,{0,-1,2,1},2}),
  latcomb_consts({8,-0.0,{0,-1,2,1},2}),
  latcomb_consts({10,-0.0,{0,1,2,1},4}),
  latcomb_consts({12,0.0,{0,1,2,-1},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,0,1},{1,0,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({1,773.5,{1,0,3,-2},2}),
  latcomb_consts({4,-773.5,{1,0,3,2},6}),
  latcomb_consts({5,-773.5,{-1,0,3,2},4}),
  latcomb_consts({8,773.5,{-1,0,3,2},4}),
  latcomb_consts({9,773.5,{1,0,3,2},6}),
  latcomb_consts({12,-773.5,{1,0,3,-2},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,-2,-1},{1,0,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({1,1547.0,{1,0,3,-2},2}),
  latcomb_consts({4,-1547.0,{1,0,3,2},6}),
  latcomb_consts({5,-1547.0,{-1,0,3,2},4}),
  latcomb_consts({8,1547.0,{-1,0,3,2},4}),
  latcomb_consts({9,1547.0,{1,0,3,2},6}),
  latcomb_consts({12,-1547.0,{1,0,3,-2},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,0,1},{1,0,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({1,773.5,{1,0,3,-2},2}),
  latcomb_consts({4,-773.5,{1,0,3,2},6}),
  latcomb_consts({5,-773.5,{-1,0,3,2},4}),
  latcomb_consts({8,773.5,{-1,0,3,2},4}),
  latcomb_consts({9,773.5,{1,0,3,2},6}),
  latcomb_consts({12,-773.5,{1,0,3,-2},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,0,-1},{1,0,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({1,-773.5,{1,0,3,-2},2}),
  latcomb_consts({4,773.5,{1,0,3,2},6}),
  latcomb_consts({5,773.5,{-1,0,3,2},4}),
  latcomb_consts({8,-773.5,{-1,0,3,2},4}),
  latcomb_consts({9,-773.5,{1,0,3,2},6}),
  latcomb_consts({12,773.5,{1,0,3,-2},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,0,-1},{1,0,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({1,-773.5,{1,0,3,-2},2}),
  latcomb_consts({4,773.5,{1,0,3,2},6}),
  latcomb_consts({5,773.5,{-1,0,3,2},4}),
  latcomb_consts({8,-773.5,{-1,0,3,2},4}),
  latcomb_consts({9,-773.5,{1,0,3,2},6}),
  latcomb_consts({12,773.5,{1,0,3,-2},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,2,1},{1,0,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({1,-1547.0,{1,0,3,-2},2}),
  latcomb_consts({4,1547.0,{1,0,3,2},6}),
  latcomb_consts({5,1547.0,{-1,0,3,2},4}),
  latcomb_consts({8,-1547.0,{-1,0,3,2},4}),
  latcomb_consts({9,-1547.0,{1,0,3,2},6}),
  latcomb_consts({12,1547.0,{1,0,3,-2},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,-2,-3},{1,0,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({1,-0.0,{1,0,3,-2},2}),
  latcomb_consts({4,0.0,{1,0,3,2},6}),
  latcomb_consts({5,0.0,{-1,0,3,2},4}),
  latcomb_consts({8,-0.0,{-1,0,3,2},4}),
  latcomb_consts({9,-0.0,{1,0,3,2},6}),
  latcomb_consts({12,0.0,{1,0,3,-2},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,2,3},{1,0,1,1}}),std::vector<latcomb_consts>({
  latcomb_consts({1,-0.0,{1,0,3,-2},2}),
  latcomb_consts({4,0.0,{1,0,3,2},6}),
  latcomb_consts({5,0.0,{-1,0,3,2},4}),
  latcomb_consts({8,-0.0,{-1,0,3,2},4}),
  latcomb_consts({9,-0.0,{1,0,3,2},6}),
  latcomb_consts({12,0.0,{1,0,3,-2},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,0,1},{1,1,0,1}}),std::vector<latcomb_consts>({
  latcomb_consts({1,1547.0,{2,3,0,-1},4}),
  latcomb_consts({2,-1547.0,{2,3,0,1},6}),
  latcomb_consts({3,-1547.0,{-2,3,0,1},2}),
  latcomb_consts({8,1547.0,{-2,3,0,1},2}),
  latcomb_consts({9,1547.0,{2,3,0,1},6}),
  latcomb_consts({10,-1547.0,{2,3,0,-1},4}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,-2,-1},{1,1,0,1}}),std::vector<latcomb_consts>({
  latcomb_consts({1,773.5,{2,3,0,-1},4}),
  latcomb_consts({2,-773.5,{2,3,0,1},6}),
  latcomb_consts({3,-773.5,{-2,3,0,1},2}),
  latcomb_consts({8,773.5,{-2,3,0,1},2}),
  latcomb_consts({9,773.5,{2,3,0,1},6}),
  latcomb_consts({10,-773.5,{2,3,0,-1},4}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,0,1},{1,1,0,1}}),std::vector<latcomb_consts>({
  latcomb_consts({1,-773.5,{2,3,0,-1},4}),
  latcomb_consts({2,773.5,{2,3,0,1},6}),
  latcomb_consts({3,773.5,{-2,3,0,1},2}),
  latcomb_consts({8,-773.5,{-2,3,0,1},2}),
  latcomb_consts({9,-773.5,{2,3,0,1},6}),
  latcomb_consts({10,773.5,{2,3,0,-1},4}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,0,-1},{1,1,0,1}}),std::vector<latcomb_consts>({
  latcomb_consts({1,-1547.0,{2,3,0,-1},4}),
  latcomb_consts({2,1547.0,{2,3,0,1},6}),
  latcomb_consts({3,1547.0,{-2,3,0,1},2}),
  latcomb_consts({8,-1547.0,{-2,3,0,1},2}),
  latcomb_consts({9,-1547.0,{2,3,0,1},6}),
  latcomb_consts({10,1547.0,{2,3,0,-1},4}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,0,-1},{1,1,0,1}}),std::vector<latcomb_consts>({
  latcomb_consts({1,773.5,{2,3,0,-1},4}),
  latcomb_consts({2,-773.5,{2,3,0,1},6}),
  latcomb_consts({3,-773.5,{-2,3,0,1},2}),
  latcomb_consts({8,773.5,{-2,3,0,1},2}),
  latcomb_consts({9,773.5,{2,3,0,1},6}),
  latcomb_consts({10,-773.5,{2,3,0,-1},4}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,2,1},{1,1,0,1}}),std::vector<latcomb_consts>({
  latcomb_consts({1,-773.5,{2,3,0,-1},4}),
  latcomb_consts({2,773.5,{2,3,0,1},6}),
  latcomb_consts({3,773.5,{-2,3,0,1},2}),
  latcomb_consts({8,-773.5,{-2,3,0,1},2}),
  latcomb_consts({9,-773.5,{2,3,0,1},6}),
  latcomb_consts({10,773.5,{2,3,0,-1},4}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,-2,-3},{1,1,0,1}}),std::vector<latcomb_consts>({
  latcomb_consts({1,-0.0,{2,3,0,-1},4}),
  latcomb_consts({2,0.0,{2,3,0,1},6}),
  latcomb_consts({3,0.0,{-2,3,0,1},2}),
  latcomb_consts({8,-0.0,{-2,3,0,1},2}),
  latcomb_consts({9,-0.0,{2,3,0,1},6}),
  latcomb_consts({10,0.0,{2,3,0,-1},4}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,2,3},{1,1,0,1}}),std::vector<latcomb_consts>({
  latcomb_consts({1,-0.0,{2,3,0,-1},4}),
  latcomb_consts({2,0.0,{2,3,0,1},6}),
  latcomb_consts({3,0.0,{-2,3,0,1},2}),
  latcomb_consts({8,-0.0,{-2,3,0,1},2}),
  latcomb_consts({9,-0.0,{2,3,0,1},6}),
  latcomb_consts({10,0.0,{2,3,0,-1},4}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,0,1},{1,1,1,0}}),std::vector<latcomb_consts>({
  latcomb_consts({1,773.5,{1,2,-1,0},2}),
  latcomb_consts({2,-773.5,{1,2,1,0},4}),
  latcomb_consts({3,-773.5,{-1,2,1,0},2}),
  latcomb_consts({4,773.5,{-1,2,1,0},2}),
  latcomb_consts({5,773.5,{1,2,1,0},4}),
  latcomb_consts({6,-773.5,{1,2,-1,0},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,-2,-1},{1,1,1,0}}),std::vector<latcomb_consts>({
  latcomb_consts({1,-0.0,{1,2,-1,0},2}),
  latcomb_consts({2,0.0,{1,2,1,0},4}),
  latcomb_consts({3,0.0,{-1,2,1,0},2}),
  latcomb_consts({4,-0.0,{-1,2,1,0},2}),
  latcomb_consts({5,-0.0,{1,2,1,0},4}),
  latcomb_consts({6,0.0,{1,2,-1,0},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,0,1},{1,1,1,0}}),std::vector<latcomb_consts>({
  latcomb_consts({1,-773.5,{1,2,-1,0},2}),
  latcomb_consts({2,773.5,{1,2,1,0},4}),
  latcomb_consts({3,773.5,{-1,2,1,0},2}),
  latcomb_consts({4,-773.5,{-1,2,1,0},2}),
  latcomb_consts({5,-773.5,{1,2,1,0},4}),
  latcomb_consts({6,773.5,{1,2,-1,0},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,0,-1},{1,1,1,0}}),std::vector<latcomb_consts>({
  latcomb_consts({1,-773.5,{1,2,-1,0},2}),
  latcomb_consts({2,773.5,{1,2,1,0},4}),
  latcomb_consts({3,773.5,{-1,2,1,0},2}),
  latcomb_consts({4,-773.5,{-1,2,1,0},2}),
  latcomb_consts({5,-773.5,{1,2,1,0},4}),
  latcomb_consts({6,773.5,{1,2,-1,0},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,0,-1},{1,1,1,0}}),std::vector<latcomb_consts>({
  latcomb_consts({1,773.5,{1,2,-1,0},2}),
  latcomb_consts({2,-773.5,{1,2,1,0},4}),
  latcomb_consts({3,-773.5,{-1,2,1,0},2}),
  latcomb_consts({4,773.5,{-1,2,1,0},2}),
  latcomb_consts({5,773.5,{1,2,1,0},4}),
  latcomb_consts({6,-773.5,{1,2,-1,0},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,2,1},{1,1,1,0}}),std::vector<latcomb_consts>({
  latcomb_consts({1,-0.0,{1,2,-1,0},2}),
  latcomb_consts({2,0.0,{1,2,1,0},4}),
  latcomb_consts({3,0.0,{-1,2,1,0},2}),
  latcomb_consts({4,-0.0,{-1,2,1,0},2}),
  latcomb_consts({5,-0.0,{1,2,1,0},4}),
  latcomb_consts({6,0.0,{1,2,-1,0},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,-1,-2,-3},{1,1,1,0}}),std::vector<latcomb_consts>({
  latcomb_consts({1,-0.0,{1,2,-1,0},2}),
  latcomb_consts({2,0.0,{1,2,1,0},4}),
  latcomb_consts({3,0.0,{-1,2,1,0},2}),
  latcomb_consts({4,-0.0,{-1,2,1,0},2}),
  latcomb_consts({5,-0.0,{1,2,1,0},4}),
  latcomb_consts({6,0.0,{1,2,-1,0},2}),
})));
LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.insert(std::make_pair(cell_valid({{0,1,2,3},{1,1,1,0}}),std::vector<latcomb_consts>({
  latcomb_consts({1,-0.0,{1,2,-1,0},2}),
  latcomb_consts({2,0.0,{1,2,1,0},4}),
  latcomb_consts({3,0.0,{-1,2,1,0},2}),
  latcomb_consts({4,-0.0,{-1,2,1,0},2}),
  latcomb_consts({5,-0.0,{1,2,1,0},4}),
  latcomb_consts({6,0.0,{1,2,-1,0},2}),
})));
}
