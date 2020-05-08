#include "L1Trigger/DTTriggerPhase2/interface/MotherGrouping.h"

using namespace edm;
using namespace std;

// ============================================================================
// Constructors and destructor
// ============================================================================
MotherGrouping::MotherGrouping(const ParameterSet& pset, edm::ConsumesCollector& iC) {
  // Obtention of parameters
  debug = pset.getUntrackedParameter<bool>("debug");
  if (debug)
    cout << "MotherGrouping: constructor" << endl;
}

MotherGrouping::~MotherGrouping() {
  if (debug)
    cout << "MotherGrouping: destructor" << endl;
}

// ============================================================================
// Main methods (initialise, run, finish)
// ============================================================================
void MotherGrouping::initialise(const edm::EventSetup& iEventSetup) {
  if (debug)
    cout << "MotherGrouping::initialiase" << endl;
}

void MotherGrouping::run(Event& iEvent,
                         const EventSetup& iEventSetup,
                         const DTDigiCollection& digis,
                         std::vector<MuonPath*>* mpaths) {
  if (debug)
    cout << "MotherGrouping: run" << endl;
}

void MotherGrouping::finish() {
  if (debug)
    cout << "MotherGrouping: finish" << endl;
};
