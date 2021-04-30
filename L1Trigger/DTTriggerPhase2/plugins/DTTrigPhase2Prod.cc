#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ESProducts.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"

#include "L1Trigger/DTTriggerPhase2/interface/MuonPath.h"
#include "L1Trigger/DTTriggerPhase2/interface/constants.h"

#include "L1Trigger/DTTriggerPhase2/interface/MotherGrouping.h"
#include "L1Trigger/DTTriggerPhase2/interface/InitialGrouping.h"
#include "L1Trigger/DTTriggerPhase2/interface/HoughGrouping.h"
#include "L1Trigger/DTTriggerPhase2/interface/PseudoBayesGrouping.h"
#include "L1Trigger/DTTriggerPhase2/interface/MuonPathAnalyzer.h"
#include "L1Trigger/DTTriggerPhase2/interface/MuonPathAnalyzerPerSL.h"
#include "L1Trigger/DTTriggerPhase2/interface/MuonPathAnalyzerInChamber.h"
#include "L1Trigger/DTTriggerPhase2/interface/MuonPathAssociator.h"
#include "L1Trigger/DTTriggerPhase2/interface/MPFilter.h"
#include "L1Trigger/DTTriggerPhase2/interface/MPQualityEnhancerFilter.h"
#include "L1Trigger/DTTriggerPhase2/interface/MPRedundantFilter.h"
#include "L1Trigger/DTTriggerPhase2/interface/MPCleanHitsFilter.h"

#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/DTSuperLayerId.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "DataFormats/DTDigi/interface/DTDigiCollection.h"
//#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTExtPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTExtPhDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTPhDigi.h"

// DT trigger GeomUtils
#include "DQM/DTMonitorModule/interface/DTTrigGeomUtils.h"

//RPC TP
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "L1Trigger/DTTriggerPhase2/interface/RPCIntegrator.h"

#include <fstream>
#include <iostream>
#include <queue>
#include <cmath>

using namespace edm;
using namespace std;
using namespace cmsdt;

class DTTrigPhase2Prod : public edm::stream::EDProducer<> {
  typedef std::map<DTChamberId, DTDigiCollection, std::less<DTChamberId>> DTDigiMap;
  typedef DTDigiMap::iterator DTDigiMap_iterator;
  typedef DTDigiMap::const_iterator DTDigiMap_const_iterator;

public:
  //! Constructor
  DTTrigPhase2Prod(const edm::ParameterSet& pset);

  //! Destructor
  ~DTTrigPhase2Prod() override;

  //! Create Trigger Units before starting event processing
  void beginRun(edm::Run const& iRun, const edm::EventSetup& iEventSetup) override;

  //! Producer: process every event and generates trigger data
  void produce(edm::Event& iEvent, const edm::EventSetup& iEventSetup) override;

  //! endRun: finish things
  void endRun(edm::Run const& iRun, const edm::EventSetup& iEventSetup) override;

  // Methods
  int rango(const metaPrimitive& mp) const;
  bool outer(const metaPrimitive& mp) const;
  bool inner(const metaPrimitive& mp) const;
  void printmP(const std::string& ss, const metaPrimitive& mP) const;
  void printmPC(const std::string& ss, const metaPrimitive& mP) const;
  bool hasPosRF(int wh, int sec) const;

  // Getter-methods
  MP_QUALITY getMinimumQuality(void);

  // Setter-methods
  void setChiSquareThreshold(float ch2Thr);
  void setMinimumQuality(MP_QUALITY q);

  // data-members
  DTGeometry const* dtGeo_;
  edm::ESGetToken<DTGeometry, MuonGeometryRecord> dtGeomH;
  std::vector<std::pair<int, MuonPath>> primitives_;

private:
  // Trigger Configuration Manager CCB validity flag
  bool my_CCBValid_;

  // BX offset used to correct DTTPG output
  int my_BXoffset_;

  // Debug Flag
  bool debug_;
  bool dump_;
  double dT0_correlate_TP_;
  bool do_correlation_;
  int scenario_;
  bool df_extended_;

  // shift
  edm::FileInPath shift_filename_;
  std::map<int, float> shiftinfo_;

  // ParameterSet
  edm::EDGetTokenT<DTDigiCollection> dtDigisToken_;
  edm::EDGetTokenT<RPCRecHitCollection> rpcRecHitsLabel_;

  // Grouping attributes and methods
  int algo_;  // Grouping code
  std::unique_ptr<MotherGrouping> grouping_obj_;
  std::unique_ptr<MuonPathAnalyzer> mpathanalyzer_;
  std::unique_ptr<MPFilter> mpathqualityenhancer_;
  std::unique_ptr<MPFilter> mpathredundantfilter_;
  std::unique_ptr<MPFilter> mpathhitsfilter_;
  std::unique_ptr<MuonPathAssociator> mpathassociator_;

  // Buffering
  bool activateBuffer_;
  int superCellhalfspacewidth_;
  float superCelltimewidth_;
  std::vector<DTDigiCollection*> distribDigis(std::queue<std::pair<DTLayerId, DTDigi>>& inQ);
  void processDigi(std::queue<std::pair<DTLayerId, DTDigi>>& inQ,
                   std::vector<std::queue<std::pair<DTLayerId, DTDigi>>*>& vec);

  // RPC
  std::unique_ptr<RPCIntegrator> rpc_integrator_;
  bool useRPC_;

  void assignIndex(std::vector<metaPrimitive>& inMPaths);
  void assignIndexPerBX(std::vector<metaPrimitive>& inMPaths);
  int assignQualityOrder(const metaPrimitive& mP) const;

  const std::unordered_map<int, int> qmap_;
};

namespace {
  struct {
    bool operator()(std::pair<DTLayerId, DTDigi> a, std::pair<DTLayerId, DTDigi> b) const {
      return (a.second.time() < b.second.time());
    }
  } DigiTimeOrdering;
}  // namespace

DTTrigPhase2Prod::DTTrigPhase2Prod(const ParameterSet& pset)
    : qmap_({{9, 9}, {8, 8}, {7, 6}, {6, 7}, {5, 3}, {4, 5}, {3, 4}, {2, 2}, {1, 1}}) {
  produces<L1Phase2MuDTPhContainer<L1Phase2MuDTPhDigi>>();
  produces<L1Phase2MuDTPhContainer<L1Phase2MuDTExtPhDigi>>();

  debug_ = pset.getUntrackedParameter<bool>("debug");
  dump_ = pset.getUntrackedParameter<bool>("dump");

  do_correlation_ = pset.getParameter<bool>("do_correlation");
  scenario_ = pset.getParameter<int>("scenario");

  df_extended_ = pset.getParameter<bool>("df_extended");

  dtDigisToken_ = consumes<DTDigiCollection>(pset.getParameter<edm::InputTag>("digiTag"));

  rpcRecHitsLabel_ = consumes<RPCRecHitCollection>(pset.getParameter<edm::InputTag>("rpcRecHits"));
  useRPC_ = pset.getParameter<bool>("useRPC");

  // Choosing grouping scheme:
  algo_ = pset.getParameter<int>("algo");

  edm::ConsumesCollector consumesColl(consumesCollector());
  
  if (algo_ == PseudoBayes || algo_ == BayesFit) {
    grouping_obj_ =
        std::make_unique<PseudoBayesGrouping>(pset.getParameter<edm::ParameterSet>("PseudoBayesPattern"), consumesColl);
  } else if (algo_ == HoughTrans) {
    grouping_obj_ =
        std::make_unique<HoughGrouping>(pset.getParameter<edm::ParameterSet>("HoughGrouping"), consumesColl);
  } else {
    grouping_obj_ = std::make_unique<InitialGrouping>(pset, consumesColl);
  }

  if (algo_ == Standard) {
    if (debug_)
      LogDebug("DTTrigPhase2Prod") << "DTp2:constructor: JM analyzer";
    mpathanalyzer_ = std::make_unique<MuonPathAnalyzerPerSL>(pset, consumesColl);
  } else {
    if (debug_)
      LogDebug("DTTrigPhase2Prod") << "DTp2:constructor: Full chamber analyzer";
    mpathanalyzer_ = std::make_unique<MuonPathAnalyzerInChamber>(pset, consumesColl);
  }

  // Getting buffer option
  activateBuffer_ = pset.getParameter<bool>("activateBuffer");
  superCellhalfspacewidth_ = pset.getParameter<int>("superCellspacewidth") / 2;
  superCelltimewidth_ = pset.getParameter<double>("superCelltimewidth");

  mpathqualityenhancer_ = std::make_unique<MPQualityEnhancerFilter>(pset);
  mpathredundantfilter_ = std::make_unique<MPRedundantFilter>(pset);
  mpathhitsfilter_ = std::make_unique<MPCleanHitsFilter>(pset);
  mpathassociator_ = std::make_unique<MuonPathAssociator>(pset, consumesColl);
  rpc_integrator_ = std::make_unique<RPCIntegrator>(pset, consumesColl);

  dtGeomH = esConsumes<DTGeometry, MuonGeometryRecord, edm::Transition::BeginRun>();
}

DTTrigPhase2Prod::~DTTrigPhase2Prod() {
  if (debug_)
    LogDebug("DTTrigPhase2Prod") << "DTp2: calling destructor" << std::endl;
}

void DTTrigPhase2Prod::beginRun(edm::Run const& iRun, const edm::EventSetup& iEventSetup) {
  if (debug_)
    LogDebug("DTTrigPhase2Prod") << "beginRun " << iRun.id().run();
  if (debug_)
    LogDebug("DTTrigPhase2Prod") << "beginRun: getting DT geometry";

  grouping_obj_->initialise(iEventSetup);          // Grouping object initialisation
  mpathanalyzer_->initialise(iEventSetup);         // Analyzer object initialisation
  mpathqualityenhancer_->initialise(iEventSetup);  // Filter object initialisation
  mpathredundantfilter_->initialise(iEventSetup);  // Filter object initialisation
  mpathhitsfilter_->initialise(iEventSetup);
  mpathassociator_->initialise(iEventSetup);       // Associator object initialisation

  const MuonGeometryRecord& geom = iEventSetup.get<MuonGeometryRecord>();
  dtGeo_ = &geom.get(dtGeomH);
}

void DTTrigPhase2Prod::produce(Event& iEvent, const EventSetup& iEventSetup) {
  if (debug_)
    LogDebug("DTTrigPhase2Prod") << "produce";
  edm::Handle<DTDigiCollection> dtdigis;
  iEvent.getByToken(dtDigisToken_, dtdigis);

  if (debug_)
    LogDebug("DTTrigPhase2Prod") << "\t Getting the RPC RecHits" << std::endl;
  edm::Handle<RPCRecHitCollection> rpcRecHits;
  iEvent.getByToken(rpcRecHitsLabel_, rpcRecHits);

  ////////////////////////////////
  // GROUPING CODE:
  ////////////////////////////////
  DTDigiMap digiMap;
  DTDigiCollection::DigiRangeIterator detUnitIt;
  for (const auto& detUnitIt : *dtdigis) {
    const DTLayerId& layId = detUnitIt.first;
    const DTChamberId chambId = layId.superlayerId().chamberId();
    const DTDigiCollection::Range& range = detUnitIt.second;
    digiMap[chambId].put(range, layId);
  }

  // generate a list muon paths for each event!!!
  if (debug_ && activateBuffer_)
    LogDebug("DTTrigPhase2Prod") << "produce - Getting and grouping digis per chamber using a buffer and super cells.";
  else if (debug_)
    LogDebug("DTTrigPhase2Prod") << "produce - Getting and grouping digis per chamber.";

  MuonPathPtrs muonpaths;
  for (const auto& ich : dtGeo_->chambers()) {
    // The code inside this for loop would ideally later fit inside a trigger unit (in principle, a DT station) of the future Phase 2 DT Trigger.
    const DTChamber* chamb = ich;
    DTChamberId chid = chamb->id();
    DTDigiMap_iterator dmit = digiMap.find(chid);

    if (dmit == digiMap.end())
      continue;

    if (activateBuffer_) {  // Use buffering (per chamber) or not
      // Import digis from the station
      std::vector<std::pair<DTLayerId, DTDigi>> tmpvec;
      tmpvec.clear();

      for (const auto& dtLayerIdIt : (*dmit).second) {
        for (DTDigiCollection::const_iterator digiIt = (dtLayerIdIt.second).first;
             digiIt != (dtLayerIdIt.second).second;
             digiIt++) {
          tmpvec.emplace_back(dtLayerIdIt.first, *digiIt);
        }
      }

      // Check to enhance CPU time usage
      if (tmpvec.empty())
        continue;

      // Order digis depending on TDC time and insert them into a queue (FIFO buffer). TODO: adapt for MC simulations.
      std::sort(tmpvec.begin(), tmpvec.end(), DigiTimeOrdering);
      std::queue<std::pair<DTLayerId, DTDigi>> timequeue;

      for (const auto& elem : tmpvec)
        timequeue.emplace(std::move(elem));
      tmpvec.clear();

      // Distribute the digis from the queue into supercells
      std::vector<DTDigiCollection*> superCells;
      superCells = distribDigis(timequeue);

      // Process each supercell & collect the resulting muonpaths (as the muonpaths std::vector is only enlarged each time
      // the groupings access it, it's not needed to "collect" the final products).
      while (!superCells.empty()) {
	grouping_obj_->run(iEvent, iEventSetup, *(superCells.back()), muonpaths);
        superCells.pop_back();
      }
    } 
    else {
      if (algo_ == BayesFit) continue;
      grouping_obj_->run(iEvent, iEventSetup, (*dmit).second, muonpaths);
    }
  }
  digiMap.clear();


  // if (dump_) {

  //   for (unsigned int i = 0; i < muonpaths.size(); i++) {
  //     stringstream ss;

  //     ss   << iEvent.id().event() << "      mpath " << i << ": ";

  //     for (int lay = 0; lay < muonpaths.at(i)->nprimitives(); lay++){
  //       ss   << muonpaths.at(i)->primitive(lay)->channelId() << " ";
  //       cout << muonpaths.at(i)->primitive(lay)->channelId() << " ";
  //     }
  //     cout << "" << endl;
  //     for (int lay = 0; lay < muonpaths.at(i)->nprimitives(); lay++){
  //       ss   << muonpaths.at(i)->primitive(lay)->tdcTimeStamp() << " ";
  //       cout << muonpaths.at(i)->primitive(lay)->tdcTimeStamp() << " ";
  //     }
  //     cout << "" << endl;
  //     for (int lay = 0; lay < muonpaths.at(i)->nprimitives(); lay++){
  //       ss   << muonpaths.at(i)->primitive(lay)->laterality() << " ";
  //       cout << muonpaths.at(i)->primitive(lay)->laterality() << " ";
  //     }
  //     cout << "" << endl;
  //     for (int lay = 0; lay < muonpaths.at(i)->nprimitives(); lay++){
  //       ss   << muonpaths.at(i)->primitive(lay)->layerId() << " ";
  //       cout << muonpaths.at(i)->primitive(lay)->layerId() << " ";
  //     }
  //     cout << "" << endl;
  //     for (int lay = 0; lay < muonpaths.at(i)->nprimitives(); lay++){
  //       ss   << muonpaths.at(i)->primitive(lay)->cameraId() << " ";
  //       cout << muonpaths.at(i)->primitive(lay)->cameraId() << " ";
  //     }
  //     cout << "" << endl;
  //     LogInfo("DTTrigPhase2Prod") << ss.str();
  //   }
  // }
  
  // FILTER GROUPING
  MuonPathPtrs filteredmuonpaths;
  if (algo_ == BayesFit){
    // Read muonpaths from txt file
    ifstream inFile;
    // Nicely hard-coded
    inFile.open("/afs/cern.ch/user/n/ntrevisa/work/DT/CMSSW_11_2_0_pre2/src/list_muons_Q1.txt");
    TString x;

    std::string str; 
    // while (std::getline(file, str))
    // std::vector<std::string> words;
    while (std::getline(inFile, str)){
      // while (inFile >> x) {
      std::istringstream iss(str);
      std::vector<std::string> words(std::istream_iterator<std::string>{iss},
					std::istream_iterator<std::string>());

      // int aaa = words.size();
      // for (int word = 0; word < aaa; word++)
      // 	cout << words[word] << endl;
      // cout << "-----------------------" << endl;
    
      DTPrimitivePtrs prims;
      int n_prim_up = 0;
      int n_prim_down = 0;
      for (int n_prim = 0; n_prim < 8; ++n_prim){
	DTPrimitivePtr prim = std::make_shared<DTPrimitive>();
	int prim_channelId = stoi(words.at(9 + n_prim)); 
	int prim_TDCtimeStamp = stoi(words.at(18 + n_prim)); 
	int prim_laterality = stoi(words.at(27 + n_prim)); 
	int prim_wheel = stoi(words.at(3));
	int prim_station = stoi(words.at(7));
	int prim_sector = stoi(words.at(5)) + 1;
	int prim_superLayer;
	int prim_layer;
	if (n_prim < 4){
	  prim -> setSuperLayerId(1); 
	  prim -> setLayerId(n_prim);
	  if (prim_TDCtimeStamp != -1) ++n_prim_down;
	  prim_superLayer = 1;
	  prim_layer = n_prim;
	}
	else{
	  prim -> setSuperLayerId(3); 
	  prim -> setLayerId(n_prim - 4);
	  if (prim_TDCtimeStamp != -1) ++n_prim_up;
	  prim_superLayer = 3;
	  prim_layer = n_prim - 4;
	}
	DTLayerId layerId(prim_wheel, prim_station, prim_sector, prim_superLayer, prim_layer);
	// DTChamberId chambId = layerId.superlayerId().chamberId();
	// DTChamberId chambId(prim_wheel, prim_station, prim_sector);
	// prim -> setCameraId(chambId);
	// cout << "Sector: STD = " << prim_sector << ", Bayes = " << layerId.sector() << endl;
	prim -> setCameraId(layerId);
	prim -> setChannelId(prim_channelId);
	prim -> setTDCTimeStamp(prim_TDCtimeStamp);
	if (prim_laterality == -1)
	  prim -> setLaterality(NONE);
	else if (prim_laterality == 0)
	  prim -> setLaterality(LEFT);
	else if (prim_laterality == 1)
	  prim -> setLaterality(RIGHT);
	else return;
	prims.push_back(prim);
      }
      // cout << "emplacing back.." << endl;
      filteredmuonpaths.emplace_back(std::make_shared<MuonPath>(prims, n_prim_up, n_prim_down));
      // cout << "emplaced!" << endl;
    }
  }
  else if (algo_ == Standard) {
    mpathredundantfilter_->run(iEvent, iEventSetup, muonpaths, filteredmuonpaths);
  }
  else {
    mpathhitsfilter_->run(iEvent, iEventSetup, muonpaths, filteredmuonpaths);
  }

  // if (dump_) {
  //   for (unsigned int i = 0; i < filteredmuonpaths.size(); i++) {
  //     stringstream ss;

  //     ss   << iEvent.id().event() << " filt. mpath " << i << ": ";
  //     cout << "" << endl;
  //     cout << iEvent.id().event() << " filt. mpath " << i << ": " << endl;

  //     for (int lay = 0; lay < filteredmuonpaths.at(i)->nprimitives(); lay++){
  //       ss   << filteredmuonpaths.at(i)->primitive(lay)->channelId() << " ";
  //       cout << filteredmuonpaths.at(i)->primitive(lay)->channelId() << " ";
  //     }
  //     cout << "" << endl;
  //     for (int lay = 0; lay < filteredmuonpaths.at(i)->nprimitives(); lay++){
  //       ss   << filteredmuonpaths.at(i)->primitive(lay)->tdcTimeStamp() << " ";
  //       cout << filteredmuonpaths.at(i)->primitive(lay)->tdcTimeStamp() << " ";
  //     }
  //     cout << "" << endl;
  //     for (int lay = 0; lay < filteredmuonpaths.at(i)->nprimitives(); lay++){
  //       ss   << filteredmuonpaths.at(i)->primitive(lay)->laterality() << " ";
  //       cout << "    " << filteredmuonpaths.at(i)->primitive(lay)->laterality() << " ";
  //     }
  //     cout << "" << endl;
  //     for (int lay = 0; lay < filteredmuonpaths.at(i)->nprimitives(); lay++){
  //       ss   << filteredmuonpaths.at(i)->primitive(lay)->cameraId() << " ";
  //       cout << filteredmuonpaths.at(i)->primitive(lay)->cameraId() << " ";
  //     }
  //     cout << "" << endl;
  //     ss   << " | quality: " << filteredmuonpaths.at(i)->quality();
  //     cout << "quality: " << filteredmuonpaths.at(i)->quality() << endl;
  //     LogInfo("DTTrigPhase2Prod") << ss.str();
  //   }
  // }

  ///////////////////////////////////////////
  /// FITTING SECTION;
  ///////////////////////////////////////////
  if (debug_)
    LogDebug("DTTrigPhase2Prod") << "MUON PATHS found: " << muonpaths.size() << " (" << filteredmuonpaths.size()
                                 << ") in event " << iEvent.id().event();

  cout << "" << endl;
  cout << "Starting fitting step" << endl;
  cout << "" << endl;
  cout << "MUON PATHS found using algorithm " << algo_ << ": " << muonpaths.size() << endl;
  cout << "MUON PATHS found using algorithm " << algo_ << " after filtering: " << filteredmuonpaths.size()
       << " in event " << iEvent.id().event() << endl;
  cout << "" << endl;

  if (debug_)
    LogDebug("DTTrigPhase2Prod") << "filling NmetaPrimtives" << std::endl;
  std::vector<metaPrimitive> metaPrimitives;
  MuonPathPtrs outmpaths;
  if (algo_ == Standard) {
    if (debug_)
      LogDebug("DTTrigPhase2Prod") << "Fitting 1SL ";
    // cout << "Fitting 1SL " << endl;
    mpathanalyzer_->run(iEvent, iEventSetup, filteredmuonpaths, metaPrimitives);
  } else {
    // implementation for advanced (2SL) grouping, no filter required..
    if (debug_)
      LogDebug("DTTrigPhase2Prod") << "Fitting 2SL at once ";
    mpathanalyzer_->run(iEvent, iEventSetup, filteredmuonpaths, outmpaths);
  }

  // if (dump_) {
  //   for (unsigned int i = 0; i < outmpaths.size(); i++) {
  //     LogInfo("DTTrigPhase2Prod") << iEvent.id().event() << " mp " << i << ": " 
  // 				  << outmpaths.at(i)->bxTimeValue() << " "
  //                                 << outmpaths.at(i)->horizPos() << " " 
  // 				  << outmpaths.at(i)->tanPhi() << " "
  //                                 << outmpaths.at(i)->phi() << " " 
  // 				  << outmpaths.at(i)->phiB() << " "
  //                                 << outmpaths.at(i)->quality() << " " 
  // 				  << outmpaths.at(i)->chiSquare();

  //     cout << iEvent.id().event() << " mp " << i << ": " 
  // 	   << outmpaths.at(i)->bxTimeValue() << " "
  // 	   << outmpaths.at(i)->horizPos() << " " 
  // 	   << outmpaths.at(i)->tanPhi() << " "
  // 	   << outmpaths.at(i)->phi() << " " 
  // 	   << outmpaths.at(i)->phiB() << " "
  // 	   << outmpaths.at(i)->quality() << " " 
  // 	   << outmpaths.at(i)->chiSquare() << endl;
      
  //     cout << iEvent.id().event() << " filt. mpath "
  //     	   << i << " passed fitting step: " << endl;
  //     for (int lay = 0; lay < outmpaths.at(i)->nprimitives(); lay++){
  //       cout << outmpaths.at(i)->primitive(lay)->channelId() << " ";
  //     }
  //     	cout << "" << endl;
  //     for (int lay = 0; lay < outmpaths.at(i)->nprimitives(); lay++){
  //       cout << outmpaths.at(i)->primitive(lay)->tdcTimeStamp() << " ";
  //     }
  //     cout << "" << endl;
  //     for (int lay = 0; lay < outmpaths.at(i)->nprimitives(); lay++){
  //       cout << "    " << outmpaths.at(i)->primitive(lay)->laterality() << " ";
  //     }
  //     cout << "" << endl;
  //     cout << "quality: " << outmpaths.at(i)->quality() << endl;
    
  //     cout << iEvent.id().event() << " mp " << i << ": " 
  //     	   << outmpaths.at(i)->bxTimeValue() << " "
  //     	   << outmpaths.at(i)->horizPos() << " " 
  //     	   << outmpaths.at(i)->tanPhi() << " "
  //     	   << outmpaths.at(i)->phi() << " " 
  //     	   << outmpaths.at(i)->phiB() << " "
  //     	   << outmpaths.at(i)->quality() << " " 
  //     	   << outmpaths.at(i)->chiSquare() << endl;

  //     cout << "" << endl;

  //   }
    // for (unsigned int i = 0; i < metaPrimitives.size(); i++) {
    //   stringstream ss;
    //   ss << iEvent.id().event() << " mp " << i << ": ";
    //   printmP(ss.str(), metaPrimitives.at(i));
    // }
    // if (algo_ == Standard) {
    //   cout << "" << endl;
    //   cout << "Muon paths using algorithm " << algo_ << " after fitting: " << metaPrimitives.size() << " in event " << iEvent.id().event() << endl;
    // }
    // else{
    //   cout << "" << endl;
    //   cout << "Muon paths using algorithm " << algo_ << " after fitting: " << outmpaths.size() << " in event " << iEvent.id().event() << endl;
    // }

  // }

  muonpaths.clear();
  filteredmuonpaths.clear();

  /////////////////////////////////////
  //  FILTER SECTIONS:
  ////////////////////////////////////
  if (debug_)
    LogDebug("DTTrigPhase2Prod") << "declaring new vector for filtered" << std::endl;

  std::vector<metaPrimitive> filteredMetaPrimitives;
  if (algo_ == Standard)
    mpathqualityenhancer_->run(iEvent, iEventSetup, metaPrimitives, filteredMetaPrimitives);

  if (dump_) {
    for (unsigned int i = 0; i < filteredMetaPrimitives.size(); i++) {
      stringstream ss;
      ss   << iEvent.id().event() << " filtered mp " << i << ": ";
      cout << iEvent.id().event() << " filtered mp " << i << ": " << endl;
      printmP(ss.str(), filteredMetaPrimitives.at(i));
    }
  }

  metaPrimitives.clear();
  metaPrimitives.erase(metaPrimitives.begin(), metaPrimitives.end());

  if (debug_)
    LogDebug("DTTrigPhase2Prod") << "DTp2 in event:" << iEvent.id().event() << " we found "
                                 << filteredMetaPrimitives.size() << " filteredMetaPrimitives (superlayer)"
                                 << std::endl;
  if (debug_)
    LogDebug("DTTrigPhase2Prod") << "filteredMetaPrimitives: starting correlations" << std::endl;

  /////////////////////////////////////
  //// CORRELATION:
  /////////////////////////////////////
  std::vector<metaPrimitive> correlatedMetaPrimitives;
  if (algo_ == Standard)
    mpathassociator_->run(iEvent, iEventSetup, dtdigis, filteredMetaPrimitives, correlatedMetaPrimitives);
  else {
    for (const auto& muonpath : outmpaths) {
      correlatedMetaPrimitives.emplace_back(muonpath->rawId(),
                                            (double)muonpath->bxTimeValue(),
                                            muonpath->horizPos(),
                                            muonpath->tanPhi(),
                                            muonpath->phi(),
                                            muonpath->phiB(),
                                            muonpath->chiSquare(),
                                            (int)muonpath->quality(),
					    muonpath->primitive(0)->channelId(), 
					    muonpath->primitive(0)->tdcTimeStamp(), 
					    muonpath->primitive(0)->laterality(), 
					    muonpath->primitive(1)->channelId(), 
					    muonpath->primitive(1)->tdcTimeStamp(), 
					    muonpath->primitive(1)->laterality(), 
					    muonpath->primitive(2)->channelId(), 
					    muonpath->primitive(2)->tdcTimeStamp(), 
					    muonpath->primitive(2)->laterality(), 
					    muonpath->primitive(3)->channelId(), 
					    muonpath->primitive(3)->tdcTimeStamp(), 
					    muonpath->primitive(3)->laterality(), 
					    muonpath->primitive(4)->channelId(), 
					    muonpath->primitive(4)->tdcTimeStamp(), 
					    muonpath->primitive(4)->laterality(), 
					    muonpath->primitive(5)->channelId(), 
					    muonpath->primitive(5)->tdcTimeStamp(), 
					    muonpath->primitive(5)->laterality(), 
					    muonpath->primitive(6)->channelId(), 
					    muonpath->primitive(6)->tdcTimeStamp(), 
					    muonpath->primitive(6)->laterality(), 
					    muonpath->primitive(7)->channelId(), 
					    muonpath->primitive(7)->tdcTimeStamp(), 
					    muonpath->primitive(7)->laterality());
    }
  }
  filteredMetaPrimitives.clear();
  


  if (debug_)
    LogDebug("DTTrigPhase2Prod") << "DTp2 in event:" << iEvent.id().event() << " we found "
                                 << correlatedMetaPrimitives.size() << " correlatedMetPrimitives (chamber)";

  if (dump_) {
    LogInfo("DTTrigPhase2Prod") << "DTp2 in event:" << iEvent.id().event() << " we found "
                                << correlatedMetaPrimitives.size() << " correlatedMetPrimitives (chamber)";

    cout << "" << endl;
    cout << "Correlated metaPrimitive found using algorithm " << algo_ << " are " << correlatedMetaPrimitives.size() << " in event " << iEvent.id().event() << endl;
    cout << "" << endl;

    for (unsigned int i = 0; i < correlatedMetaPrimitives.size(); i++) {
      stringstream ss;
      ss   << iEvent.id().event() << " correlated mp " << i << ": ";
      cout << iEvent.id().event() << " correlated mp " << i << ": " << endl;
      printmPC(ss.str(), correlatedMetaPrimitives.at(i));
      
    }
  }

  double shift_back = 0;
  if (scenario_ == MC)  //scope for MC
    shift_back = 400;
  else if (scenario_ == DATA)  //scope for data
    shift_back = 0;
  else if (scenario_ == SLICE_TEST)  //scope for slice test
    shift_back = 0;

  // RPC integration
  if (useRPC_) {
    rpc_integrator_->initialise(iEventSetup, shift_back);
    rpc_integrator_->prepareMetaPrimitives(rpcRecHits);
    rpc_integrator_->matchWithDTAndUseRPCTime(correlatedMetaPrimitives);
    rpc_integrator_->makeRPCOnlySegments();
    rpc_integrator_->storeRPCSingleHits();
    rpc_integrator_->removeRPCHitsUsed();
  }

  /// STORING RESULTs
  vector<L1Phase2MuDTPhDigi> outP2Ph;
  vector<L1Phase2MuDTExtPhDigi> outExtP2Ph;
  
  // Assigning index value
  assignIndex(correlatedMetaPrimitives);
  for (const auto& metaPrimitiveIt : correlatedMetaPrimitives) {
    DTChamberId chId(metaPrimitiveIt.rawId);
    if (debug_)
      LogDebug("DTTrigPhase2Prod") << "looping in final vector: SuperLayerId" << chId << " x=" << metaPrimitiveIt.x
                                   << " quality=" << metaPrimitiveIt.quality
                                   << " BX=" << round(metaPrimitiveIt.t0 / 25.) << " index=" << metaPrimitiveIt.index;

    int sectorTP = chId.sector();
    //sectors 13 and 14 exist only for the outermost stations for sectors 4 and 10 respectively
    //due to the larger MB4 that are divided into two.
    if (sectorTP == 13)
      sectorTP = 4;
    if (sectorTP == 14)
      sectorTP = 10;
    sectorTP = sectorTP - 1;
    int sl = 0;
    // if (metaPrimitiveIt.quality < LOWLOWQ || metaPrimitiveIt.quality == CHIGHQ) {
    // If the metaPrimitive is not correlated, give it a reference SL
    if (metaPrimitiveIt.quality < H3PLUS3) {
      if (inner(metaPrimitiveIt))
        sl = 1;
      else
        sl = 3;
    }

    if (debug_)
      LogDebug("DTTrigPhase2Prod") << "pushing back phase-2 dataformat carlo-federica dataformat";

   
    if (df_extended_ == true){

      int pathWireId[8] = {metaPrimitiveIt.wi1,metaPrimitiveIt.wi2,metaPrimitiveIt.wi3,metaPrimitiveIt.wi4,
			   metaPrimitiveIt.wi5,metaPrimitiveIt.wi6,metaPrimitiveIt.wi7,metaPrimitiveIt.wi8};

      int pathTDC[8] = {metaPrimitiveIt.tdc1,metaPrimitiveIt.tdc2,metaPrimitiveIt.tdc3,metaPrimitiveIt.tdc4,
			metaPrimitiveIt.tdc5,metaPrimitiveIt.tdc6,metaPrimitiveIt.tdc7,metaPrimitiveIt.tdc8};
    
      int pathLat[8] = {metaPrimitiveIt.lat1,metaPrimitiveIt.lat2,metaPrimitiveIt.lat3,metaPrimitiveIt.lat4,
			metaPrimitiveIt.lat5,metaPrimitiveIt.lat6,metaPrimitiveIt.lat7,metaPrimitiveIt.lat8};
    
      outExtP2Ph.emplace_back(L1Phase2MuDTExtPhDigi(
					      (int)round(metaPrimitiveIt.t0 / 25.) - shift_back,  // ubx (m_bx) //bx en la orbita
					      chId.wheel(),    // uwh (m_wheel)     // FIXME: It is not clear who provides this?
					      sectorTP,        // usc (m_sector)    // FIXME: It is not clear who provides this?
					      chId.station(),  // ust (m_station)
					      sl,              // ust (m_station)
					      (int)round(metaPrimitiveIt.phi * 65536. / 0.8),    // uphi (_phiAngle)
					      (int)round(metaPrimitiveIt.phiB * 2048. / 1.4),    // uphib (m_phiBending)
					      metaPrimitiveIt.quality,                           // uqua (m_qualityCode)
					      metaPrimitiveIt.index,                             // uind (m_segmentIndex)
					      (int)round(metaPrimitiveIt.t0) - shift_back * 25,  // ut0 (m_t0Segment)
					      (int)round(metaPrimitiveIt.chi2 * 1000000),        // uchi2 (m_chi2Segment)
					      (int)round(metaPrimitiveIt.x * 1000),              // ux (m_xLocal)
					      (int)round(metaPrimitiveIt.tanPhi * 1000),         // utan (m_tanPsi)
					      (int)round(metaPrimitiveIt.phi * 65536. / 0.8),    // uphi (_phiAngleCMSSW)     <-- in theory. For the moment, same as phi
					      (int)round(metaPrimitiveIt.phiB * 2048. / 1.4),    // uphib (m_phiBendingCMSSW) <-- in theory. For the moment, same as phiB
					      metaPrimitiveIt.rpcFlag,                           // urpc (m_rpcFlag)
					      pathWireId,
					      pathTDC,
					      pathLat
					      ));
      // cout << "Local x in muon path = " << outP2Ph.xLocal() << endl; 
    }
    else{
    
      outP2Ph.emplace_back(L1Phase2MuDTPhDigi(
					   (int)round(metaPrimitiveIt.t0 / (float)LHC_CLK_FREQ) - shift_back,
					   chId.wheel(),                                                // uwh (m_wheel)
					   sectorTP,                                                    // usc (m_sector)
					   chId.station(),                                              // ust (m_station)
					   sl,                                                          // ust (m_station)
					   (int)round(metaPrimitiveIt.phi * PHIRES_CONV),               // uphi (_phiAngle)
					   (int)round(metaPrimitiveIt.phiB * PHIBRES_CONV),             // uphib (m_phiBending)
					   metaPrimitiveIt.quality,                                     // uqua (m_qualityCode)
					   metaPrimitiveIt.index,                                       // uind (m_segmentIndex)
					   (int)round(metaPrimitiveIt.t0) - shift_back * LHC_CLK_FREQ,  // ut0 (m_t0Segment)
					   (int)round(metaPrimitiveIt.chi2 * CHI2RES_CONV),             // uchi2 (m_chi2Segment)
					   metaPrimitiveIt.rpcFlag                                      // urpc (m_rpcFlag)
					   ));
    }
  }

  // Storing RPC hits that were not used elsewhere
  //  if (useRPC_) {
  //    for (auto rpc_dt_digi = rpc_integrator_->rpcRecHits_translated_.begin();
  //         rpc_dt_digi != rpc_integrator_->rpcRecHits_translated_.end();
  //         rpc_dt_digi++) {
  //      outP2Ph.push_back(*rpc_dt_digi);
  //    }
  //  }

  if (df_extended_) {
    auto resultP2Ph = std::make_unique<L1Phase2MuDTPhContainer<L1Phase2MuDTExtPhDigi>>();
    resultP2Ph->setContainer(outExtP2Ph);
    iEvent.put(std::move(resultP2Ph));
  }
  else {
    auto resultP2Ph = std::make_unique<L1Phase2MuDTPhContainer<L1Phase2MuDTPhDigi>>();
    resultP2Ph->setContainer(outP2Ph);
    iEvent.put(std::move(resultP2Ph));
  }
  outExtP2Ph.clear();
  outExtP2Ph.erase(outExtP2Ph.begin(), outExtP2Ph.end());
  outP2Ph.clear();
  outP2Ph.erase(outP2Ph.begin(), outP2Ph.end());

  // if (df_extended_ == true){
  //   auto resultP2Ph = std::make_unique<L1Phase2MuDTExtPhContainer>();
  //   resultP2Ph->setContainer(outP2Ph);
  //   iEvent.put(std::move(resultP2Ph));
  //   outP2Ph.clear();
  //   outP2Ph.erase(outP2Ph.begin(), outP2Ph.end());
  // }
  // else{
  
  // auto resultP2Ph = std::make_unique<L1Phase2MuDTExtPhContainer>();
  // auto resultP2Ph = std::make_unique<L1Phase2MuDTPhContainer>();
  // resultP2Ph->setContainer(outP2Ph);
  // iEvent.put(std::move(resultP2Ph));
  // outP2Ph.clear();
  // outP2Ph.erase(outP2Ph.begin(), outP2Ph.end());
  // }
}

void DTTrigPhase2Prod::endRun(edm::Run const& iRun, const edm::EventSetup& iEventSetup) {
  grouping_obj_->finish();
  mpathanalyzer_->finish();
  mpathqualityenhancer_->finish();
  mpathredundantfilter_->finish();
  mpathhitsfilter_->finish();
  mpathassociator_->finish();
  rpc_integrator_->finish();
};

bool DTTrigPhase2Prod::outer(const metaPrimitive& mp) const {
  int counter = (mp.wi5 != -1) + (mp.wi6 != -1) + (mp.wi7 != -1) + (mp.wi8 != -1);
  return (counter > 2);
}

bool DTTrigPhase2Prod::inner(const metaPrimitive& mp) const {
  int counter = (mp.wi1 != -1) + (mp.wi2 != -1) + (mp.wi3 != -1) + (mp.wi4 != -1);
  return (counter > 2);
}

bool DTTrigPhase2Prod::hasPosRF(int wh, int sec) const { return wh > 0 || (wh == 0 && sec % 4 > 1); }

void DTTrigPhase2Prod::printmP(const string& ss, const metaPrimitive& mP) const {
  DTSuperLayerId slId(mP.rawId);
  LogInfo("DTTrigPhase2Prod") << ss << (int)slId << "\t " << setw(2) << left << mP.wi1 << " " << setw(2) << left
                              << mP.wi2 << " " << setw(2) << left << mP.wi3 << " " << setw(2) << left << mP.wi4 << " "
                              << setw(5) << left << mP.tdc1 << " " << setw(5) << left << mP.tdc2 << " " << setw(5)
                              << left << mP.tdc3 << " " << setw(5) << left << mP.tdc4 << " " << setw(10) << right
                              << mP.x << " " << setw(9) << left << mP.tanPhi << " " << setw(5) << left << mP.t0 << " "
                              << setw(13) << left << mP.chi2 << " r:" << rango(mP);
}

void DTTrigPhase2Prod::printmPC(const string& ss, const metaPrimitive& mP) const {
  DTChamberId ChId(mP.rawId);
  LogInfo("DTTrigPhase2Prod") << ss << (int)ChId << "\t  " << setw(2) << left << mP.wi1 << " " << setw(2) << left
                              << mP.wi2 << " " << setw(2) << left << mP.wi3 << " " << setw(2) << left << mP.wi4 << " "
                              << setw(2) << left << mP.wi5 << " " << setw(2) << left << mP.wi6 << " " << setw(2) << left
                              << mP.wi7 << " " << setw(2) << left << mP.wi8 << " " << setw(5) << left << mP.tdc1 << " "
                              << setw(5) << left << mP.tdc2 << " " << setw(5) << left << mP.tdc3 << " " << setw(5)
                              << left << mP.tdc4 << " " << setw(5) << left << mP.tdc5 << " " << setw(5) << left
                              << mP.tdc6 << " " << setw(5) << left << mP.tdc7 << " " << setw(5) << left << mP.tdc8
                              << " " << setw(2) << left << mP.lat1 << " " << setw(2) << left << mP.lat2 << " "
                              << setw(2) << left << mP.lat3 << " " << setw(2) << left << mP.lat4 << " " << setw(2)
                              << left << mP.lat5 << " " << setw(2) << left << mP.lat6 << " " << setw(2) << left
                              << mP.lat7 << " " << setw(2) << left << mP.lat8 << " " << setw(10) << right << mP.x << " "
                              << setw(9) << left << mP.tanPhi << " " << setw(5) << left << mP.t0 << " " << setw(13)
                              << left << mP.chi2 << " r:" << rango(mP);


   // cout << ss << (int)ChId << "\t  " << setw(2) << left << mP.wi1 << " " << setw(2) << left
   // 	<< mP.wi2 << " " << setw(2) << left << mP.wi3 << " " << setw(2) << left << mP.wi4 << " "
   // 	<< setw(2) << left << mP.wi5 << " " << setw(2) << left << mP.wi6 << " " << setw(2) << left
   // 	<< mP.wi7 << " " << setw(2) << left << mP.wi8 << " " << setw(5) << left << mP.tdc1 << " "
   // 	<< setw(5) << left << mP.tdc2 << " " << setw(5) << left << mP.tdc3 << " " << setw(5)
   // 	<< left << mP.tdc4 << " " << setw(5) << left << mP.tdc5 << " " << setw(5) << left
   // 	<< mP.tdc6 << " " << setw(5) << left << mP.tdc7 << " " << setw(5) << left << mP.tdc8
   // 	<< " " << setw(2) << left << mP.lat1 << " " << setw(2) << left << mP.lat2 << " "
   // 	<< setw(2) << left << mP.lat3 << " " << setw(2) << left << mP.lat4 << " " << setw(2)
   // 	<< left << mP.lat5 << " " << setw(2) << left << mP.lat6 << " " << setw(2) << left
   // 	<< mP.lat7 << " " << setw(2) << left << mP.lat8 << " " << setw(10) << right << mP.x << " "
   // 	<< setw(9) << left << mP.tanPhi << " " << setw(5) << left << mP.t0 << " " << setw(13)
   // 	<< left << mP.chi2 << " r:" << rango(mP) << endl;
}

int DTTrigPhase2Prod::rango(const metaPrimitive& mp) const {
  if (mp.quality == 1 or mp.quality == 2)
    return 3;
  if (mp.quality == 3 or mp.quality == 4)
    return 4;
  return mp.quality;
}

void DTTrigPhase2Prod::assignIndex(std::vector<metaPrimitive>& inMPaths) {
  std::map<int, std::vector<metaPrimitive>> primsPerBX;
  for (const auto& metaPrimitive : inMPaths) {
    int BX = round(metaPrimitive.t0 / 25.);
    primsPerBX[BX].push_back(metaPrimitive);
  }
  inMPaths.clear();
  for (auto& prims : primsPerBX) {
    assignIndexPerBX(prims.second);
    for (const auto& primitive : prims.second)
      inMPaths.push_back(primitive);
  }
}

void DTTrigPhase2Prod::assignIndexPerBX(std::vector<metaPrimitive>& inMPaths) {
  // First we asociate a new index to the metaprimitive depending on quality or phiB;
  uint32_t rawId = -1;
  int numP = -1;
  for (auto& metaPrimitiveIt : inMPaths) {
    numP++;
    rawId = metaPrimitiveIt.rawId;
    int iOrder = assignQualityOrder(metaPrimitiveIt);
    int inf = 0;
    int numP2 = -1;
    for (auto& metaPrimitiveItN : inMPaths) {
      int nOrder = assignQualityOrder(metaPrimitiveItN);
      numP2++;
      if (rawId != metaPrimitiveItN.rawId)
        continue;
      if (numP2 == numP) {
        metaPrimitiveIt.index = inf;
        break;
      } else if (iOrder < nOrder) {
        inf++;
      } else if (iOrder > nOrder) {
        metaPrimitiveItN.index++;
      } else if (iOrder == nOrder) {
        if (std::abs(metaPrimitiveIt.phiB) >= std::abs(metaPrimitiveItN.phiB)) {
          inf++;
        } else if (std::abs(metaPrimitiveIt.phiB) < std::abs(metaPrimitiveItN.phiB)) {
          metaPrimitiveItN.index++;
        }
      }
    }  // ending second for
  }    // ending first for
}

int DTTrigPhase2Prod::assignQualityOrder(const metaPrimitive& mP) const {
  if (mP.quality > 9 || mP.quality < 1)
    return -1;

  return qmap_.find(mP.quality)->second;
}

std::vector<DTDigiCollection*> DTTrigPhase2Prod::distribDigis(std::queue<std::pair<DTLayerId, DTDigi>>& inQ) {
  std::vector<std::queue<std::pair<DTLayerId, DTDigi>>*> tmpVector;
  tmpVector.clear();
  std::vector<DTDigiCollection*> collVector;
  collVector.clear();
  while (!inQ.empty()) {
    processDigi(inQ, tmpVector);
  }
  for (auto& sQ : tmpVector) {
    DTDigiCollection tmpColl;
    while (!sQ->empty()) {
      tmpColl.insertDigi((sQ->front().first), (sQ->front().second));
      sQ->pop();
    }
    collVector.push_back(&tmpColl);
  }
  return collVector;
}

void DTTrigPhase2Prod::processDigi(std::queue<std::pair<DTLayerId, DTDigi>>& inQ,
                                   std::vector<std::queue<std::pair<DTLayerId, DTDigi>>*>& vec) {
  bool classified = false;
  if (!vec.empty()) {
    for (auto& sC : vec) {  // Conditions for entering a super cell.
      if ((sC->front().second.time() + superCelltimewidth_) > inQ.front().second.time()) {
        // Time requirement
        if (TMath::Abs(sC->front().second.wire() - inQ.front().second.wire()) <= superCellhalfspacewidth_) {
          // Spatial requirement
          sC->push(std::move(inQ.front()));
          classified = true;
        }
      }
    }
  }
  if (classified) {
    inQ.pop();
    return;
  }

  std::queue<std::pair<DTLayerId, DTDigi>> newQueue;

  std::pair<DTLayerId, DTDigi> tmpPair;
  tmpPair = std::move(inQ.front());
  newQueue.push(tmpPair);
  inQ.pop();

  vec.push_back(&newQueue);
}

DEFINE_FWK_MODULE(DTTrigPhase2Prod);
