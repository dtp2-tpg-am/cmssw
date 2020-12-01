#ifndef L1Trigger_DTTriggerPhase2_MuonPathAnalyticAnalyzer_h
#define L1Trigger_DTTriggerPhase2_MuonPathAnalyticAnalyzer_h

#include "L1Trigger/DTTriggerPhase2/interface/MuonPathAnalyzer.h"

// ===============================================================================
// Previous definitions and declarations
// ===============================================================================

struct cell_valid {
  int cell_horiz_layout[4]; 
  int valid[4];

  bool operator =(const cell_valid &o) const {
    std::cout << "Equal " << std::endl; 
    for (int i = 0; i<4; i++) {
      std::cout << cell_horiz_layout[i] << " " <<   o.cell_horiz_layout[i] << " " << valid[i] << " " << o.valid[i] << std::endl; 
      if (cell_horiz_layout[i] != o.cell_horiz_layout[i] || valid[i] != o.valid[i] ) return false; 
    }
    return true;
  }

  bool operator <(const cell_valid &o) const {
    /*for (int i = 0; i<4; i++) {
      std::cout << cell_horiz_layout[i] << " " <<   o.cell_horiz_layout[i] << " " << valid[i] << " " << o.valid[i] << std::endl; 
      if ( !( (cell_horiz_layout[i] < o.cell_horiz_layout[i]) || ( (cell_horiz_layout[i] == o.cell_horiz_layout[i]) &&  valid[i] < o.valid[i]) ) ) return false; 
    }
    std::cout << "Smaller " << std::endl; */
    return true;
  } 
};

struct latcomb_consts {
  int latcomb;
  double numer_const; 
  int numer_coeff[4]; 
  int denom; 
};

// ===============================================================================
// Class declarations
// ===============================================================================

class MuonPathAnalyticAnalyzer : public MuonPathAnalyzer {
public:
  // Constructors and destructor
  MuonPathAnalyticAnalyzer(const edm::ParameterSet &pset, edm::ConsumesCollector &iC);
  ~MuonPathAnalyticAnalyzer() override;

  // Main methods
  void initialise(const edm::EventSetup &iEventSetup) override;
  void run(edm::Event &iEvent,
           const edm::EventSetup &iEventSetup,
           MuonPathPtrs &inMpath,
           std::vector<cmsdt::metaPrimitive> &metaPrimitives) override;
  void run(edm::Event &iEvent,
           const edm::EventSetup &iEventSetup,
           MuonPathPtrs &inMpath,
           MuonPathPtrs &outMPath) override{};

  void finish() override;

  // Other public methods
  void setBXTolerance(int t) { bxTolerance_ = t; };
  int bxTolerance(void) { return bxTolerance_; };

  void setChiSquareThreshold(float ch2Thr) { chiSquareThreshold_ = ch2Thr; };

  void setMinQuality(cmsdt::MP_QUALITY q) {
    if (minQuality_ >= cmsdt::LOWQGHOST)
      minQuality_ = q;
  };
  cmsdt::MP_QUALITY minQuality(void) { return minQuality_; };

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
  void analyze(MuonPathPtr &inMPath, std::vector<cmsdt::metaPrimitive> &metaPrimitives);

  void setCellLayout(const int layout[cmsdt::NUM_LAYERS]);
  void buildLateralities(void);
  bool isStraightPath(cmsdt::LATERAL_CASES sideComb[cmsdt::NUM_LAYERS]);
  
  int eqMainTerm(cmsdt::LATERAL_CASES sideComb[2], int layerIdx[2], MuonPathPtr &mPath, int bxValue);
  void lateralCoeficients(cmsdt::LATERAL_CASES sideComb[2], int *coefs);
  void calculatePathParameters(MuonPathPtr &mPath);
  void calcTanPhiXPosChamber(MuonPathPtr &mPath);
  void calcCellDriftAndXcoor(MuonPathPtr &mPath);
  void calcChiSquare(MuonPathPtr &mPath);

  void calcTanPhiXPosChamber3Hits(MuonPathPtr &mPath);
  void calcTanPhiXPosChamber4Hits(MuonPathPtr &mPath);

  int omittedHit(int idx);
  
  // NEW ANALYZER FUNCTIONS
 
  void fillLAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER();
  
  void superlayer_datapath(MuonPathPtr &mPath);
  void segment_composer(int valids[cmsdt::NUM_LAYERS], int wires[cmsdt::NUM_LAYERS], int t0s[cmsdt::NUM_LAYERS],
    int cell_horiz_layout[cmsdt::NUM_LAYERS], MuonPathPtr &mPath);
  int validate_4hits_one_latcomb_reduced(int valids[cmsdt::NUM_LAYERS], int t0s[cmsdt::NUM_LAYERS],
    int cell_horiz_layout[cmsdt::NUM_LAYERS], latcomb_consts my_latcomb_consts);
  int tuck_t0(int t, int t0s[cmsdt::NUM_LAYERS], int valid[cmsdt::NUM_LAYERS]);
  int changeLatToEmu (int lat);

  // Private attributes

  static const int LAYER_ARRANGEMENTS_[cmsdt::NUM_LAYERS][cmsdt::NUM_CELL_COMB];
  cmsdt::LATERAL_CASES lateralities_[cmsdt::NUM_LATERALITIES][cmsdt::NUM_LAYERS];
  cmsdt::LATQ_TYPE latQuality_[cmsdt::NUM_LATERALITIES];

  int totalNumValLateralities_;

  int bxTolerance_;
  cmsdt::MP_QUALITY minQuality_;
  float chiSquareThreshold_;
  bool debug_;
  double chi2Th_;
  double chi2corTh_;
  double tanPhiTh_;
  int cellLayout_[cmsdt::NUM_LAYERS];
  bool use_LSB_;
  double tanPsi_precision_;
  double x_precision_;
  std::map<cell_valid, std::vector<latcomb_consts>> LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER;

};

#endif
