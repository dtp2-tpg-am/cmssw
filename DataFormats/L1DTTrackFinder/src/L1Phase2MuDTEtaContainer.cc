#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTEtaContainer.h"

L1Phase2MuDTEtaContainer::L1Phase2MuDTEtaContainer() {}

void L1Phase2MuDTEtaContainer::setContainer(const Segment_Container& inputSegments) { m_segments = inputSegments; }

L1Phase2MuDTEtaContainer::Segment_Container const* L1Phase2MuDTEtaContainer::getContainer() const { return &m_segments; }
