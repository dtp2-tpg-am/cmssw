//-------------------------------------------------
//
//   Class L1MuDTChambContainer
//
//   Description: trigger primtive data for the
//                muon barrel Phase2 trigger
//--------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTEtaContainer.h"

//-------------------------------
// Collaborating Class Headers --
//-------------------------------

//---------------
// C++ Headers --
//---------------

//-------------------
// Initializations --
//-------------------

//----------------
// Constructors --
//----------------
L1Phase2MuDTEtaContainer::L1Phase2MuDTEtaContainer() {}

//--------------
// Operations --
//--------------
void L1Phase2MuDTEtaContainer::setContainer(const Segment_Container& inputSegments) { m_segments = inputSegments; }

L1Phase2MuDTEtaContainer::Segment_Container const* L1Phase2MuDTEtaContainer::getContainer() const { return &m_segments; }
