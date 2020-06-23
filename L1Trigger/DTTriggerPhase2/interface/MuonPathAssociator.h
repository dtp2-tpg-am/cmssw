#ifndef Phase2L1Trigger_DTTrigger_MuonPathAssociator_cc
#define Phase2L1Trigger_DTTrigger_MuonPathAssociator_cc

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/FrameworkfwdMostUsed.h"

#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/DTSuperLayerId.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "DataFormats/DTDigi/interface/DTDigiCollection.h"

#include "L1Trigger/DTTriggerPhase2/interface/MuonPath.h"
#include "L1Trigger/DTTriggerPhase2/interface/constants.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"

#include <iostream>
#include <fstream>

// ===============================================================================
// Previous definitions and declarations
// ===============================================================================

// ===============================================================================
// Class declarations
// ===============================================================================

class MuonPathAssociator {
public:
  // Constructors and destructor
  MuonPathAssociator(const edm::ParameterSet &pset, edm::ConsumesCollector &iC);
  ~MuonPathAssociator();

  // Main methods
  void initialise(const edm::EventSetup &iEventSetup);
  void run(edm::Event &iEvent,
           const edm::EventSetup &iEventSetup,
           edm::Handle<DTDigiCollection> digis,
           std::vector<metaPrimitive> &inMPaths,
           std::vector<metaPrimitive> &outMPaths);

  void finish();

  // Other public methods

  bool shareFit(metaPrimitive first, metaPrimitive second);
  bool isNotAPrimo(metaPrimitive first, metaPrimitive second);
  void removeSharingFits(std::vector<metaPrimitive> &chamberMPaths, std::vector<metaPrimitive> &allMPaths);
  void removeSharingHits(std::vector<metaPrimitive> &firstMPaths,
                         std::vector<metaPrimitive> &secondMPaths,
                         std::vector<metaPrimitive> &allMPaths);
  void printmPC(metaPrimitive mP);

  // Public attributes
  DTGeometry const *dtGeo_;
  edm::ESGetToken<DTGeometry, MuonGeometryRecord> dtGeomH_;

private:
  // Private methods
  void correlateMPaths(edm::Handle<DTDigiCollection> digis,
                       std::vector<metaPrimitive> &inMPaths,
                       std::vector<metaPrimitive> &outMPaths);


  bool hasPosRF(int wh, int sec) { return wh > 0 || (wh == 0 && sec % 4 > 1); }

  // Private attributes
  bool debug_;
  bool clean_chi2_correlation_;
  bool useBX_correlation_;
  bool allow_confirmation_;
  double dT0_correlate_TP_;
  double dBX_correlate_TP_;
  double dTanPsi_correlate_TP_;
  double minx_match_2digis_;
  double chi2corTh_;
  bool use_LSB_;
  double tanPsi_precision_;
  double x_precision_;

  //shift
  edm::FileInPath shift_filename_;
  std::map<int, float> shiftinfo_;
};

#endif
