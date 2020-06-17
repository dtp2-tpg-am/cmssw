#include "L1Trigger/DTTriggerPhase2/interface/InitialGrouping.h"

using namespace edm;
using namespace std;
using namespace cmsdt;
// ============================================================================
// Constructors and destructor
// ============================================================================
InitialGrouping::InitialGrouping(const ParameterSet &pset, edm::ConsumesCollector &iC)
    : MotherGrouping(pset, iC), currentBaseChannel(-1) {
  // Obtention of parameters
  debug = pset.getUntrackedParameter<bool>("debug");
  if (debug)
    cout << "InitialGrouping: constructor" << endl;

  chInDummy.push_back(DTPrimitivePtr(new DTPrimitive()));
  // Initialisation of channelIn array
  //  for (int lay = 0; lay < NUM_LAYERS; lay++) {
  //    for (int ch = 0; ch < NUM_CH_PER_LAYER; ch++) {
  //      channelIn[lay][ch] = chInDummy;
  //      channelIn[lay][ch].clear();
  //    }
  //  }
}

InitialGrouping::~InitialGrouping() {
  if (debug)
    cout << "InitialGrouping: destructor" << endl;
}

// ============================================================================
// Main methods (initialise, run, finish)
// ============================================================================
void InitialGrouping::initialise(const edm::EventSetup &iEventSetup) {
  if (debug)
    cout << "InitialGrouping::initialiase" << endl;
}

void InitialGrouping::run(Event &iEvent,
                          const EventSetup &iEventSetup,
                          const DTDigiCollection &digis,
                          MuonPathPtrs &mpaths) {
  if (debug)
    cout << "InitialGrouping: run" << endl;

  //   This function returns the analyzable mpath collection back to the the main function
  //   so it can be fitted. This is in fact doing the so-called grouping.

  for (int supLayer = 0; supLayer < NUM_SUPERLAYERS; supLayer++) {  // for each SL:
    if (debug)
      cout << "InitialGrouping::run Reading SL" << supLayer << endl;
    setInChannels(&digis, supLayer);

    for (int baseCh = 0; baseCh < TOTAL_BTI; baseCh++) {
      currentBaseChannel = baseCh;
      selectInChannels(currentBaseChannel);  //map a number of wires for a given base channel
      if (notEnoughDataInChannels())
        continue;

      if (debug)
        cout << "InitialGrouping::run --> now check pathId" << endl;
      for (int pathId = 0; pathId < 8; pathId++) {
        resetPrvTDCTStamp();
        if (debug)
          cout << "[InitialGrouping::run] mixChannels calling" << endl;
        mixChannels(supLayer, pathId, mpaths);
        if (debug)
          cout << "[InitialGrouping::run] mixChannels end" << endl;
      }
    }
  }
  if (debug)
    cout << "[InitialGrouping::run] end" << endl;
}

void InitialGrouping::finish() { return; };

// ============================================================================
// Other methods
// ============================================================================
void InitialGrouping::setInChannels(const DTDigiCollection *digis, int sl) {
  //   before setting channels we need to clear
  for (int lay = 0; lay < NUM_LAYERS; lay++) {
    for (int ch = 0; ch < NUM_CH_PER_LAYER; ch++) {
      channelIn[lay][ch].clear();
    }
  }

  // now fill with those primitives that makes sense:
  DTDigiCollection::DigiRangeIterator dtLayerId_It;
  for (dtLayerId_It = digis->begin(); dtLayerId_It != digis->end(); ++dtLayerId_It) {
    const DTLayerId dtLId = (*dtLayerId_It).first;
    if (dtLId.superlayer() != sl + 1)
      continue;  //skip digis not in SL...

    for (DTDigiCollection::const_iterator digiIt = ((*dtLayerId_It).second).first;
         digiIt != ((*dtLayerId_It).second).second;
         ++digiIt) {
      int layer = dtLId.layer() - 1;
      int wire = (*digiIt).wire() - 1;
      int digiTIME = (*digiIt).time();
      int digiTIMEPhase2 = digiTIME;

      auto dtpAux = DTPrimitivePtr(new DTPrimitive());
      dtpAux->setTDCTimeStamp(digiTIMEPhase2);
      dtpAux->setChannelId(wire);
      dtpAux->setLayerId(layer);    //  L=0, 1, 2, 3
      dtpAux->setSuperLayerId(sl);  // SL=0,1,2
      dtpAux->setCameraId(dtLId.rawId());
      channelIn[layer][wire].push_back(std::move(dtpAux));
    }
  }
}

void InitialGrouping::selectInChannels(int baseChannel) {
  // Channels are labeled following next schema:
  // Input Muxer Indexes
  // ---------------------------------
  // |   6   |   7   |   8   |   9   |
  // ---------------------------------
  // |   3   |   4   |   5   |
  // -------------------------
  // |   1   |   2   |
  // -----------------
  // |   0   |
  // ---------

  // ****** LAYER 0 ******
  muxInChannels[0] = channelIn[0][baseChannel];

  // ****** LAYER 1 ******
  muxInChannels[1] = channelIn[1][baseChannel];

  if (baseChannel + 1 < NUM_CH_PER_LAYER)
    muxInChannels[2] = channelIn[1][baseChannel + 1];
  else
    muxInChannels[2] = chInDummy;

  // ****** LAYER 2 ******
  if (baseChannel - 1 >= 0)
    muxInChannels[3] = channelIn[2][baseChannel - 1];
  else
    muxInChannels[3] = chInDummy;

  muxInChannels[4] = channelIn[2][baseChannel];

  if (baseChannel + 1 < NUM_CH_PER_LAYER)
    muxInChannels[5] = channelIn[2][baseChannel + 1];
  else
    muxInChannels[5] = chInDummy;

  // ****** LAYER 3 ******
  if (baseChannel - 1 >= 0)
    muxInChannels[6] = channelIn[3][baseChannel - 1];
  else
    muxInChannels[6] = chInDummy;

  muxInChannels[7] = channelIn[3][baseChannel];

  if (baseChannel + 1 < NUM_CH_PER_LAYER)
    muxInChannels[8] = channelIn[3][baseChannel + 1];
  else
    muxInChannels[8] = chInDummy;

  if (baseChannel + 2 < NUM_CH_PER_LAYER)
    muxInChannels[9] = channelIn[3][baseChannel + 2];
  else
    muxInChannels[9] = chInDummy;
}

bool InitialGrouping::notEnoughDataInChannels(void) {
  // Empty layer indicators
  bool lEmpty[4];

  lEmpty[0] = muxInChannels[0].empty();

  lEmpty[1] = muxInChannels[1].empty() && muxInChannels[2].empty();

  lEmpty[2] = muxInChannels[3].empty() && muxInChannels[4].empty() && muxInChannels[5].empty();

  lEmpty[3] =
      muxInChannels[6].empty() && muxInChannels[7].empty() && muxInChannels[8].empty() && muxInChannels[9].empty();

  // If there are at least two empty layers, you cannot link it to a possible trace
  if ((lEmpty[0] && lEmpty[1]) or (lEmpty[0] && lEmpty[2]) or (lEmpty[0] && lEmpty[3]) or (lEmpty[1] && lEmpty[2]) or
      (lEmpty[1] && lEmpty[3]) or (lEmpty[2] && lEmpty[3])) {
    return true;
  } else {
    return false;
  }
}

void InitialGrouping::resetPrvTDCTStamp(void) {
  for (int i = 0; i <= 3; i++)
    prevTDCTimeStamps[i] = -1;
}

bool InitialGrouping::isEqualComb2Previous(DTPrimitivePtrs dtPrims) {
  bool answer = true;

  for (int i = 0; i < (int)dtPrims.size(); i++) {
    if (dtPrims.at(i).get() == nullptr)
      continue;
    if (prevTDCTimeStamps[i] != dtPrims.at(i)->tdcTimeStamp()) {
      answer = false;
      for (int j = 0; j < (int)dtPrims.size(); j++) {
        prevTDCTimeStamps[j] = dtPrims.at(j)->tdcTimeStamp();
      }
      break;
    }
  }
  return answer;
}

void InitialGrouping::mixChannels(int supLayer, int pathId, MuonPathPtrs &outMuonPath) {
  if (debug)
    cout << "[InitialGrouping::mixChannel] begin" << endl;
  DTPrimitivePtrs data[4];

  int horizLayout[4];
  memcpy(horizLayout, CELL_HORIZONTAL_LAYOUTS[pathId], 4 * sizeof(int));

  int chIdxForPath[4];
  memcpy(chIdxForPath, CHANNELS_PATH_ARRANGEMENTS[pathId], 4 * sizeof(int));

  // Real amount of values extracted from each channel.
  int numPrimsPerLayer[4] = {0, 0, 0, 0};
  unsigned int canal;
  int channelEmptyCnt = 0;
  for (int layer = 0; layer <= 3; layer++) {
    canal = CHANNELS_PATH_ARRANGEMENTS[pathId][layer];
    if (muxInChannels[canal].empty())
      channelEmptyCnt++;
  }

  if (channelEmptyCnt >= 2)
    return;
  //

  // We extract the number of elements necesary from each channel as the combination requires
  for (int layer = 0; layer <= 3; layer++) {
    canal = CHANNELS_PATH_ARRANGEMENTS[pathId][layer];
    unsigned int maxPrimsToBeRetrieved = muxInChannels[canal].size();
    /*
    If the number of primitives is zero, in order to avoid that only one
    empty channel avoids mixing data from the other three, we, at least,
    consider one dummy element from this channel.
    In other cases, where two or more channels has zero elements, the final
    combination will be not analyzable (the condition for being analyzable is
    that it has at least three good TDC time values, not dummy), so it will
    be discarded and not sent to the analyzer.
  */
    if (maxPrimsToBeRetrieved == 0)
      maxPrimsToBeRetrieved = 1;

    for (unsigned int items = 0; items < maxPrimsToBeRetrieved; items++) {
      if (muxInChannels[canal].size() == 0)
        continue;

      auto dtpAux = DTPrimitivePtr(std::move(muxInChannels[canal].at(items)));
      /*
        I won't allow a whole loop cycle. When a DTPrimitive has an invalid
        time-stamp (TDC value = -1) it means that the buffer is empty or the
        buffer has reached the last element within the configurable time window.
        In this case the loop is broken, but only if there is, at least, one
        DTPrim (even invalid) on the outgoing array. This is mandatory to cope
        with the idea explained in the previous comment block
      */
      if (!dtpAux)
        break;
      if (dtpAux->tdcTimeStamp() < 0 && items > 0)
        break;

      // In this new schema, if the hit corresponds with the SL over which
      // you are doing the mixings, it is sent to the intermediate mixing
      // buffer. In the opposite case, a blank and invalid copy is sent to
      // allow them mixing to be complete, as it was done in the one SL case.

      // This is a kind of quick solution in which there will be no few cases
      // where you will have invalid mixings. Because of that, the verification
      // that is done later, where the segment is analysed to check whether it
      // can be analysed is essential.
      if (dtpAux->superLayerId() == supLayer)
        data[layer].push_back(std::move(dtpAux));  // values are 0, 1, 2
      else
        data[layer].push_back(DTPrimitivePtr(new DTPrimitive()));
      numPrimsPerLayer[layer]++;
    }
  }

  // Here we do the different combinations and send them to the output FIFO.
  int chIdx[4];
  for (chIdx[0] = 0; chIdx[0] < numPrimsPerLayer[0]; chIdx[0]++) {
    for (chIdx[1] = 0; chIdx[1] < numPrimsPerLayer[1]; chIdx[1]++) {
      for (chIdx[2] = 0; chIdx[2] < numPrimsPerLayer[2]; chIdx[2]++) {
        for (chIdx[3] = 0; chIdx[3] < numPrimsPerLayer[3]; chIdx[3]++) {
          // We build a copy of the object so that we can manipulate each one
          // in each thread of the process independently, allowing us also to
          // delete them whenever it is necessary, without relying upon a
          // unique reference all over the code.

          DTPrimitivePtrs ptrPrimitive;
          for (int i = 0; i <= 3; i++) {
            ptrPrimitive.push_back(std::move((data[i])[chIdx[i]]));
          }

          auto ptrMuonPath = MuonPathPtr(new MuonPath(ptrPrimitive));
          ptrMuonPath->setCellHorizontalLayout(horizLayout);

          /*
            This new version of this code is redundant with PathAnalyzer code,
            where every MuonPath not analyzable is discarded.
            I insert this discarding mechanism here, as well, to avoid inserting
            not-analyzable MuonPath into the candidate FIFO.
            Equivalent code must be removed in the future from PathAnalyzer, but
            it the mean time, at least during the testing state, I'll preserve
            both.
            Code in the PathAnalyzer should be doing nothing now.
          */
          if (ptrMuonPath->isAnalyzable()) {
            /*
            This is a very simple filter because, during the tests, it has been
            detected that many consecutive MuonPaths are duplicated mainly due
            to buffers empty (or dummy) that give a TDC time-stamp = -1
            With this filter, I'm removing those consecutive identical
            combinations.
            
            If duplicated combinations are not consecutive, they won't be
            detected here
          */
            ptrMuonPath->setBaseChannelId(currentBaseChannel);
            outMuonPath.push_back(std::move(ptrMuonPath));
            ptrPrimitive.clear();
          }
        }
      }
    }
  }
  for (int layer = 0; layer <= 3; layer++) {
    data[layer].clear();
  }
}
