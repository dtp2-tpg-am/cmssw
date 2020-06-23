#ifndef L1Trigger_DTTriggerPhase2_MuonPathAnalyzerPerSL_h
#define L1Trigger_DTTriggerPhase2_MuonPathAnalyzerPerSL_h

#include "L1Trigger/DTTriggerPhase2/interface/MuonPathAnalyzer.h"

// ===============================================================================
// Previous definitions and declarations
// ===============================================================================
using namespace cmsdt;
// ===============================================================================
// Class declarations
// ===============================================================================

class MuonPathAnalyzerPerSL : public MuonPathAnalyzer {
public:
  // Constructors and destructor
  MuonPathAnalyzerPerSL(const edm::ParameterSet &pset, edm::ConsumesCollector &iC);
  virtual ~MuonPathAnalyzerPerSL();

  // Main methods
  void initialise(const edm::EventSetup &iEventSetup);
  void run(edm::Event &iEvent,
           const edm::EventSetup &iEventSetup,
           MuonPathPtrs &inMpath,
           std::vector<metaPrimitive> &metaPrimitives);
  void run(edm::Event &iEvent, const edm::EventSetup &iEventSetup, MuonPathPtrs &inMpath, MuonPathPtrs &outMPath){};

  void finish();

  // Other public methods
  void setBXTolerance(int t) { bxTolerance_ = t; };
  int bxTolerance(void) { return bxTolerance_; };

  void setChiSquareThreshold(float ch2Thr) { chiSquareThreshold_ = ch2Thr; };

  void setMinQuality(MP_QUALITY q) {
    if (minQuality_ >= LOWQGHOST)
      minQuality_ = q;
  };
  MP_QUALITY minQuality(void) { return minQuality_; };

  bool hasPosRF(int wh, int sec) { return wh > 0 || (wh == 0 && sec % 4 > 1); };

  // Public attributes
  DTGeometry const *dtGeo_;
  edm::ESGetToken<DTGeometry, MuonGeometryRecord> dtGeomH;

  //shift
  edm::FileInPath shift_filename_;
  std::map<int, float> shiftinfo_;

  int chosen_sl_;

private:
  // Private methods
  void analyze(MuonPathPtr &inMPath, std::vector<metaPrimitive> &metaPrimitives);

  void setCellLayout(const int layout[NUM_LAYERS]);
  void buildLateralities(void);
  bool isStraightPath(LATERAL_CASES sideComb[NUM_LAYERS]);

  void evaluatePathQuality(MuonPathPtr &mPath);
  void evaluateLateralQuality(int latIdx, MuonPathPtr &mPath, LATQ_TYPE *latQuality);
  void validate(LATERAL_CASES sideComb[3], int layerIndex[3], MuonPathPtr &mPath, PARTIAL_LATQ_TYPE *latq);

  int eqMainBXTerm(LATERAL_CASES sideComb[2], int layerIdx[2], MuonPathPtr &mPath);

  int eqMainTerm(LATERAL_CASES sideComb[2], int layerIdx[2], MuonPathPtr &mPath, int bxValue);

  void lateralCoeficients(LATERAL_CASES sideComb[2], int *coefs);
  bool sameBXValue(PARTIAL_LATQ_TYPE *latq);

  void calculatePathParameters(MuonPathPtr &mPath);
  void calcTanPhiXPosChamber(MuonPathPtr &mPath);
  void calcCellDriftAndXcoor(MuonPathPtr &mPath);
  void calcChiSquare(MuonPathPtr &mPath);

  void calcTanPhiXPosChamber3Hits(MuonPathPtr &mPath);
  void calcTanPhiXPosChamber4Hits(MuonPathPtr &mPath);

  int omittedHit(int idx);

  // Private attributes

  static const int LAYER_ARRANGEMENTS_[NUM_LAYERS][NUM_CELL_COMB];
  LATERAL_CASES lateralities_[NUM_LATERALITIES][NUM_LAYERS];
  LATQ_TYPE latQuality_[NUM_LATERALITIES];

  int totalNumValLateralities_;

  int bxTolerance_;
  MP_QUALITY minQuality_;
  float chiSquareThreshold_;
  bool debug_;
  double chi2Th_;
  double chi2corTh_;
  double tanPhiTh_;
  int cellLayout_[NUM_LAYERS];
  bool use_LSB_;
  double tanPsi_precision_;
  double x_precision_;
};

#endif
