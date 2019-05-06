#ifndef Phase2L1Trigger_DTTrigger_MuonPathAnalyzerPerSL_cc
#define Phase2L1Trigger_DTTrigger_MuonPathAnalyzerPerSL_cc

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/DTSuperLayerId.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"

#include "L1Trigger/DTTriggerPhase2/interface/muonpath.h"
#include "L1Trigger/DTTriggerPhase2/interface/analtypedefs.h"
#include "L1Trigger/DTTriggerPhase2/interface/constants.h"

#include "L1Trigger/DTTriggerPhase2/interface/MuonPathAnalyzer.h"
#include "L1Trigger/DTTriggerPhase2/interface/InitialGrouping.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"

#include <iostream>
#include <fstream>
#include <math.h>

// ===============================================================================
// Previous definitions and declarations
// ===============================================================================
//const Int_t CELL_HORIZONTAL_LAYOUTS[8][4] = {
//  {0, -1, -2, -3}, {0, -1, -2, -1}, {0, -1, 0, -1}, {0, -1, 0, 1},
//  {0,  1,  0, -1}, {0,  1,  0,  1}, {0,  1, 2,  1}, {0,  1, 2, 3}
//};

//const int MAX_VERT_ARRANG=4;
// ===============================================================================
// Class declarations
// ===============================================================================

class MuonPathAnalyzerPerSL : public MuonPathAnalyzer {
 public:
   public:
  typedef struct {
    bool latQValid;
    int  bxValue;
  } PARTIAL_LATQ_TYPE;
  
  typedef struct {
    bool valid;
    int bxValue;
    int invalidateHitIdx;
    MP_QUALITY quality;
  } LATQ_TYPE;
  
  // Constructors and destructor
  MuonPathAnalyzerPerSL(const edm::ParameterSet& pset);
  virtual ~MuonPathAnalyzerPerSL();
  
  // Main methods
  void initialise(const edm::EventSetup& iEventSetup) override;
  void run(edm::Event& iEvent, const edm::EventSetup& iEventSetup, std::vector<MuonPath*> &inMpath, std::vector<metaPrimitive> &metaPrimitives) override;
  void run(edm::Event& iEvent, const edm::EventSetup& iEventSetup, std::vector<MuonPath*> &inMpath, std::vector<MuonPath*> &outMPath) override;
  void finish() override;
  
  // Other public methods
  void setBXTolerance(int t) { bxTolerance = t; }
  int  getBXTolerance(void)  { return bxTolerance; }
  void setChiSquareThreshold(float ch2Thr) { chiSquareThreshold = ch2Thr; }
  void setMinimumQuality(MP_QUALITY q) { if (minQuality >= LOWQGHOST) minQuality = q; }
  MP_QUALITY getMinimumQuality(void) { return minQuality; }
  
  // Public attributes
  //ttrig
  std::string ttrig_filename;
  std::map<int,float> ttriginfo;
  
  //z
  std::string z_filename;
  std::map<int,float> zinfo;
  
  //shift
  std::string shift_filename;
  std::map<int,float> shiftinfo;

  edm::ESHandle<DTGeometry> dtGeo;
  
 private:
  
  // Private methods
  void analyze(MuonPath *inMPath, std::vector<metaPrimitive> &metaPrimitives); 
  void analyze(MuonPath *inMPath, std::vector<MuonPath*> &outMPath); 
  
  void setCellLayout(const int layout[4]);
  void buildLateralities(void);
  bool isStraightPath(LATERAL_CASES sideComb[4]);
  
  /* This determines whether the values of the 4 primitives make up a 
     trajectory. The values have to be placed in the layer order: 
     0    -> innermost layer 
     1, 2 -> next layers
     3    -> outermost layes  */
  void evaluatePathQuality(MuonPath *mPath);
  void evaluateLateralQuality(int latIdx, MuonPath *mPath, LATQ_TYPE *latQuality);
  void validate(LATERAL_CASES sideComb[3], int layerIndex[3], MuonPath* mPath, PARTIAL_LATQ_TYPE *latq);
  
  int eqMainBXTerm(LATERAL_CASES sideComb[2], int layerIdx[2],
		   MuonPath* mPath);
  
  int eqMainTerm(LATERAL_CASES sideComb[2], int layerIdx[2], MuonPath* mPath,
		 int bxValue);
  
  void getLateralCoeficients(LATERAL_CASES sideComb[2], int *coefs);
  bool sameBXValue(PARTIAL_LATQ_TYPE *latq);
  bool hasPosRF(int wh,int sec);
  
  void calculatePathParameters(MuonPath *mPath);
  void calcTanPhiXPosChamber  (MuonPath *mPath);
  void calcCellDriftAndXcoor  (MuonPath *mPath);
  void calcChiSquare          (MuonPath *mPath);
  
  void calcTanPhiXPosChamber3Hits(MuonPath *mPath);
  void calcTanPhiXPosChamber4Hits(MuonPath *mPath);
  
  int getOmittedHit(int idx);
  
  
  // Private attributes
  Bool_t debug;
  
  static const int LAYER_ARRANGEMENTS[MAX_VERT_ARRANG][3]; 

  LATERAL_CASES lateralities[16][4];
  LATQ_TYPE latQuality[16];
  
  int cellLayout[4];
  int bxTolerance;
  MP_QUALITY minQuality;
  
  float chiSquareThreshold;
  int totalNumValLateralities;
  
  int chosen_sl;
  
};


#endif
