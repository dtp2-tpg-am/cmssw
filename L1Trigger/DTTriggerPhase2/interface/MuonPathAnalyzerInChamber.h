#ifndef L1Trigger_DTTriggerPhase2_MuonPathAnalyzerInChamber_cc
#define L1Trigger_DTTriggerPhase2_MuonPathAnalyzerInChamber_cc

#include "L1Trigger/DTTriggerPhase2/interface/MuonPathAnalyzer.h"

// ===============================================================================
// Previous definitions and declarations
// ===============================================================================
namespace {
  constexpr int NLayers = 8;
  typedef std::array<LATERAL_CASES, NLayers> TLateralities;
}  // namespace
// ===============================================================================
// Class declarations
// ===============================================================================

class MuonPathAnalyzerInChamber : public MuonPathAnalyzer {
public:
  // Constructors and destructor
  MuonPathAnalyzerInChamber(const edm::ParameterSet &pset, edm::ConsumesCollector &iC);
  virtual ~MuonPathAnalyzerInChamber();

  // Main methods
  void initialise(const edm::EventSetup &iEventSetup);
  void run(edm::Event &iEvent,
           const edm::EventSetup &iEventSetup,
           MuonPathPtrs &inMpath,
           std::vector<metaPrimitive> &metaPrimitives) {}
  void run(edm::Event &iEvent, const edm::EventSetup &iEventSetup, MuonPathPtrs &inMpath, MuonPathPtrs &outMPath);

  void finish();

  // Other public methods
  void setBxTolerance(int t) { bxTolerance_ = t; };
  void setMinHits4Fit(int h) { minHits4Fit_ = h; };
  void setChiSquareThreshold(float ch2Thr) { chiSquareThreshold_ = ch2Thr; };
  void setMinimumQuality(MP_QUALITY q) {
    if (minQuality_ >= LOWQGHOST)
      minQuality_ = q;
  };

  int bxTolerance(void) { return bxTolerance_; };
  int minHits4Fit(void) { return minHits4Fit_; };
  MP_QUALITY minQuality(void) { return minQuality_; };

  bool hasPosRF(int wh, int sec) { return wh > 0 || (wh == 0 && sec % 4 > 1); };

  // Public attributes
  DTGeometry const *dtGeo_;
  edm::ESGetToken<DTGeometry, MuonGeometryRecord> dtGeomH;

  //shift
  std::map<int, float> shiftinfo_;

private:
  // Private methods
  void analyze(MuonPathPtr &inMPath, MuonPathPtrs &outMPaths);

  void setCellLayout(MuonPathPtr &mpath);
  void buildLateralities(MuonPathPtr &mpath);
  void setLateralitiesInMP(MuonPathPtr &mpath, TLateralities lat);
  void setWirePosAndTimeInMP(MuonPathPtr &mpath);
  void calculateFitParameters(MuonPathPtr &mpath, TLateralities lat, int present_layer[NLayers]);

  void evaluateQuality(MuonPathPtr &mPath);
  int totalNumValLateralities_;
  std::vector<TLateralities> lateralities_;
  std::vector<LATQ_TYPE> latQuality_;

  bool debug_;
  double chi2Th_;
  edm::FileInPath shift_filename_;
  int bxTolerance_;
  MP_QUALITY minQuality_;
  float chiSquareThreshold_;
  short minHits4Fit_;
  int cellLayout_[NLayers];
};

#endif
