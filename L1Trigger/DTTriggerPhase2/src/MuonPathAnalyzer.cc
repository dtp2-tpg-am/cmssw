#include "L1Trigger/DTTriggerPhase2/interface/MuonPathAnalyzer.h"

using namespace edm;
using namespace std;

// ============================================================================
// Constructors and destructor
// ============================================================================
MuonPathAnalyzer::MuonPathAnalyzer(const ParameterSet& pset, edm::ConsumesCollector& iC) {
  // Obtention of parameters
  debug = pset.getUntrackedParameter<bool>("debug");
  if (debug)
    cout << "MuonPathAnalyzer: constructor" << endl;
}

MuonPathAnalyzer::~MuonPathAnalyzer() {
  if (debug)
    cout << "MuonPathAnalyzer: destructor" << endl;
}

// ============================================================================
// Main methods (initialise, run, finish)
// ============================================================================
void MuonPathAnalyzer::initialise(const edm::EventSetup& iEventSetup) {
  if (debug)
    cout << "MuonPathAnalyzer::initialiase" << endl;
}

/*
void MuonPathAnalyzer::run(edm::Event& iEvent, const edm::EventSetup& iEventSetup, std::vector<MuonPath*> *inMpath, std::vector<MuonPath*> &outMpath) {
  if (debug) cout <<"MuonPathAnalyzer: run" << endl;  

}
*/
void MuonPathAnalyzer::finish() {
  if (debug)
    cout << "MuonPathAnalyzer: finish" << endl;
};
