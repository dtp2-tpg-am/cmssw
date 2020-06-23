#ifndef Phase2L1Trigger_DTTrigger_MPRedundantFilter_cc
#define Phase2L1Trigger_DTTrigger_MPRedundantFilter_cc

#include "L1Trigger/DTTriggerPhase2/interface/MPFilter.h"


#include <iostream>
#include <fstream>
#include <deque>

// ===============================================================================
// Previous definitions and declarations
// ===============================================================================

// ===============================================================================
// Class declarations
// ===============================================================================

class MPRedundantFilter : public MPFilter {
public:
  // Constructors and destructor
  MPRedundantFilter(const edm::ParameterSet& pset);
  virtual ~MPRedundantFilter();

  // Main methods
  void initialise(const edm::EventSetup& iEventSetup);
  void run(edm::Event& iEvent,
           const edm::EventSetup& iEventSetup,
           std::vector<metaPrimitive>& inMPath,
           std::vector<metaPrimitive>& outMPath){};
  void run(edm::Event& iEvent, const edm::EventSetup& iEventSetup, MuonPathPtrs& inMPath, MuonPathPtrs& outMPath);
  void finish() { buffer_.clear(); };

  // Other public methods

private:
  void filter(MuonPathPtr& mpath, MuonPathPtrs& outMPaths);
  bool isInBuffer(MuonPathPtr& mpath);

  // Private attributes
  bool debug_;
  unsigned int maxBufferSize_;
  std::deque<MuonPathPtr> buffer_;
};

#endif
