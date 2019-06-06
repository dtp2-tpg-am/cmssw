// Modified version of the DTAB7RawToDigi unpacker by Oscar Gonzalez to improve
// readability and maintenance.
// Started on (2019_01_22)

//-------------------------------------------------
//
//   Class: DTAB7RawToDigi
//
//   L1 DT AB7 Raw-to-Digi
//
//
//
//   Author :
//   C. Heidemann - RWTH Aachen
//   J. Troconiz  - UAM
//   y el eterno retorno
//
//   Modification
//    M.C Fouz - CIEMAT. 14 June 2017
//    Interchanging Phi1 & Phi2 assigned wrongly
//    Forcing the chamber to be in Wheel-1
//        Since is the good type for the demonstrator
//
//--------------------------------------------------

#include "EventFilter/DTRawToDigi/plugins/OglezDTAB7RawToDigi.h"

// CMS Classes

#include "FWCore/Framework/interface/LuminosityBlock.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CondFormats/DTObjects/interface/DTReadOutMapping.h"
#include "CondFormats/DataRecord/interface/DTReadOutMappingRcd.h"

// DT-related classes

#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "EventFilter/DTRawToDigi/interface/DTROChainCoding.h"

//OLD #include "EventFilter/DTRawToDigi/plugins/BunchIDfromTDC.h"

#include <iostream>

//-----------------------------------------------------------------------
OglezDTAB7RawToDigi::OglezDTAB7RawToDigi(const edm::ParameterSet& pset)
// Constructor of the class
{
  // Collections to be produced for the EventRecord:
  produces<DTDigiCollection>();
  //  produces<L1MuDTChambPhContainer>();
  produces<L1Phase2MuDTPhContainer>();

  // Reading parameters

  DTAB7InputTag_ = pset.getParameter<edm::InputTag>("DTAB7_FED_Source");

  debug_ = pset.getUntrackedParameter<bool>("debug", false);

  doHexDumping_ = pset.getUntrackedParameter<bool>("doHexDumping", false);

  feds_ = pset.getUntrackedParameter<std::vector<int> >("feds", std::vector<int>());

//OLD  nfeds_ = feds_.size();

  rawToken_ = consumes<FEDRawDataCollection>(DTAB7InputTag_);
}

//-----------------------------------------------------------------------
OglezDTAB7RawToDigi::~OglezDTAB7RawToDigi(){}
// Destructor of the class.

//-----------------------------------------------------------------------
void OglezDTAB7RawToDigi::produce(edm::Event& e, const edm::EventSetup& c)
// Main routine to process in every event.
{
  digis_=new DTDigiCollection();
  primitives_.clear();

//OLD  DTDigiCollection digis;
//OLD  std::vector<L1MuDTChambPhDigi> primitives;

  if (!fillRawData(e, c)) return;

  auto AB7DTDigi_product = std::make_unique<DTDigiCollection>(*digis_);


  // Original structure of the trigger primitives
  //L1MuDTChambPhContainer primContainer;
  //primContainer.setContainer(primitives_);
  //auto AB7DTPrim_product = std::make_unique<L1MuDTChambPhContainer>(primContainer);

  // Phase2 structure of the trigger primitives (by Federica Primavera)
  L1Phase2MuDTPhContainer primContainer;
  primContainer.setContainer(primitives_);
  auto AB7DTPrim_product = std::make_unique<L1Phase2MuDTPhContainer>(primContainer);

//  e.put(std::move(AB7DTDigi_product), "DTAB7Digis");
//  e.put(std::move(AB7DTPrim_product), "DTAB7Primitives");
  e.put(std::move(AB7DTDigi_product));
  e.put(std::move(AB7DTPrim_product));

  delete digis_;
}

//-----------------------------------------------------------------------
bool OglezDTAB7RawToDigi::fillRawData(edm::Event& e,
                                 const edm::EventSetup& c
                                 //OLD DTDigiCollection& digis,
                                 //OLD std::vector<L1MuDTChambPhDigi>& primitives
                                 ) {

  edm::Handle<FEDRawDataCollection> data;
  e.getByToken(rawToken_, data);

  edm::ESHandle<DTReadOutMapping> mapping;
  c.get<DTReadOutMappingRcd>().get( mapping );

  //std::cout << " New Event.............................................................................." << std::endl;
//OLD  for (int w_i = 0; w_i < nfeds_; ++w_i) {
//OLD    process(feds_[w_i], data, mapping); //OLD digis, primitives);
//OLD  }
  for (auto &xval : feds_) {
    if (doHexDumping_)
      std::cout << "OGINFO-INFO: Dump ++++++ Information about FED: "<<xval<<" (in event: "
                <<e.run()<<"-"<<e.getLuminosityBlock().luminosityBlock()
                <<"-"<<e.eventAuxiliary().event()<<")"<<std::endl;

    process(xval, data, mapping);
  }

  return true;
}

//-----------------------------------------------------------------------
void OglezDTAB7RawToDigi::process(int DTAB7FED,
                                  edm::Handle<FEDRawDataCollection> data,
                                  edm::ESHandle<DTReadOutMapping> mapping
                                  //OLD DTDigiCollection& digis,
                                  //OLD std::vector<L1MuDTChambPhDigi>& primitives
                                  ) {

  // Reading the data for the requested FED:

  FEDRawData dturosdata = data->FEDData(DTAB7FED);
  if ( dturosdata.size() == 0 ) {
    if ( debug_ ) edm::LogWarning("oglez_dtab7_unpacker") << "No data for the requested FED "<<DTAB7FED <<std::endl;
    return;
  }

  fedLinePointer_=dturosdata.data();
  long dataWord=0;
  lineCounter_=0;

  hitOrder_.clear();

  newCRC_ = 0xFFFF;

  // Reading the headers:

  readLine(&dataWord);  // Reading the word
  calcCRC(dataWord);

  // Bits 60-63 should be 0x05
  if ( ((dataWord>>60)&0x0F) != 0x5) {
    if ( debug_ ) edm::LogWarning("oglez_dtab7_unpacker")
                    << "Not a DTAB7 header code: "
                    << std::hex << dataWord << std::dec;
    return;
  }

  // Bits 20-31 is the bunch crossing as given by the FED
  bxCounter_=((dataWord>>20)&0xFFF);

  // Bits 8-19 should be the FED number (source ID, it is called)
  if ( ((dataWord>>8)&0x0FFF) != DTAB7FED) {
    if ( debug_ ) edm::LogWarning("oglez_dtab7_unpacker")
                    << "Data does not correspond to the FED "
                    << DTAB7FED<<" "
                    << std::hex << dataWord << std::dec;
    return;
  }

  if (0 || doHexDumping_) std::cout << "OGINFO-INFO: Dump + Header/Block, BX: "<<bxCounter_<<std::endl;

  readLine(&dataWord);   // Second word of the header
  calcCRC(dataWord);

      int nslots = (dataWord>>52) & 0xF; // bits 52-55 Number of AMC/slots

//  std::cout<<"OG-INFO: Reading fed "<<DTAB7FED<<" with slots: "<<nslots<<std::endl;
//OLD      std::cout<<"+++++++++++++++++++++++++++++++++++++++++++++ Starting AMC "<<nslots<<std::endl;

  if (doHexDumping_) std::cout << "OGINFO-INFO: Dump + Header/muFOV word, nAMCs: "<<nslots<<" orbit number: "<<((dataWord>>4)&0xFFFFFFFF)<<std::endl;

  // Reading the information for the slots: one word per AMC
  std::map<int,int> slot_size;
  for (int j=0;j<nslots;++j) {
    readLine(&dataWord);   // AMC word for the slot
    calcCRC(dataWord);

    int slot=(dataWord>>16)&0xF; // Bits 16-19:
    if (slot<1 || slot>12) {
      if ( debug_ ) edm::LogWarning("oglez_dtab7_unpacker")
                      << "AMCnumber "
                      << slot << " out of range (1-12)";
      return;
    }
    slot_size[slot] = (dataWord>>32)&0xFFFFFF; // bits 32-55: n words for the AMC
    if (0 || doHexDumping_) std::cout << "OGINFO-INFO: Dump + AMC: "<<slot<<" size: "<<slot_size[slot]<<std::endl;
  }

  // Reading all the payloads for the AB7... each AMC
  for (int j=0;j<nslots;++j) {
    readLine(&dataWord);  // First word header AMC
    calcCRC(dataWord);

    int slot = (dataWord>>56)&0xF;      // bit 56-59: slot number
    if (slot<1 || slot>12) {
      if ( debug_ ) edm::LogWarning("oglez_dtab7_unpacker")
                      << "AMCnumber "
                      << slot << " for the information is out of range (1-12)";
      return;
    }

    // Collision information (from Carsten's code in the original(?))

    //int dataLenght = (dataWord & 0xFFFFF);         // positions 0 -> 19
    int bx = (dataWord >> 20 ) & 0xFFF;    // positions 20 -> 31
   int event      = (dataWord >> 32 ) & 0xFFFFFF; // positions 32 -> 55
//   int AMC_ID     = (dataWord >> 56 ) & 0xF;      // positions 56 -> 59
//   int control    = (dataWord >> 60 ) & 0xF;      // positions 59 -> 63
//   int wheel      = AB7Wheel;

   // CHECK!
   if (bx!=bxCounter_) {
     edm::LogWarning("oglez_dtab7_unpacker")
       << "AMC Slot has too many digis/TP: "<<slot_size[slot]<<std::endl;
   }

   if (0 || doHexDumping_) std::cout << "OGINFO-INFO: Dump + AMC Header 1, BX: "<<bx<<" event: "<<event<<std::endl;

    // Second word header AMC: nothing relevant
    readLine(&dataWord);
    calcCRC(dataWord);

    // Check on the number of slots... if large, we exclude the entry.
    if (slot_size[slot]>200) {
      if ( debug_ ) edm::LogWarning("oglez_dtab7_unpacker")
                      << "AMC Slot has too many digis/TP: "<<slot_size[slot]<<std::endl;
      if (slot_size[slot]>1000) return;  // Problematic events... code crashes due to crappy/crazy values?
    }

    // Reading the hits or trigger primitives (words!)
    for (int k=slot_size[slot]-1;k>2;--k) { // Just a counter of how many!
      readLine(&dataWord);
      calcCRC(dataWord);

      // Reading the AB7PayLoad:
      int type=(dataWord>>60)&0xF; // Bits 60-63 gives us the type of information:

      //OLD std::cout<<"OGDT-INFO: "<<std::hex<<type<<" "<<dataWord<<" "<<std::dec<<std::endl;
      if (doHexDumping_) std::cout<<"OGINFO-INFO: Dump + word info for digi/tp: "<<std::hex<<type<<std::dec<<std::endl;

      if (type==1) {  // Hit information: 2 hits of 30 bits!
        readAB7PayLoad_hitWord(dataWord,DTAB7FED,slot);
      }
      else if (((type>>2)&0x3)==2) {  // First trigger word... The other should come after wards!
        long firstWord=dataWord;
        readLine(&dataWord);  // Reading the second word
        calcCRC(dataWord);
        --k;

        // Checking it is fine...
        if (((dataWord>>62)&0x3)==3)  {  // It is ok!
          readAB7PayLoad_triggerPrimitive(firstWord,dataWord,DTAB7FED,slot);
        }
        else if ( debug_ ) edm::LogWarning("oglez_dtab7_unpacker")
                         << "Expected second trigger word that is not there "
                         <<DTAB7FED<<" "<<slot<<std::hex<<type<<" "<<((dataWord>>62)&0x3)<<" "<<dataWord<<std::dec;
      }
      else {
        // Problems with the value of the word
        if ( debug_ ) edm::LogWarning("oglez_dtab7_unpacker")
                         << "Error word for the ros information "
                         <<DTAB7FED<<" "<<slot<<" "<<std::hex<<type<<std::dec;
      }
    }

    // Trailer word of the AMC information
    readLine(&dataWord);
    calcCRC(dataWord);
  }

  // Trailer words for checks

  readLine(&dataWord);   // First trailer word (for the block) is not used, but it has some check information
  calcCRC(dataWord);

  if (doHexDumping_) std::cout << "OGINFO-INFO: Dump + Trailer/muFOV, BX: "<<(dataWord&0xFFF)<<std::endl;

  readLine(&dataWord);   // Second trailer word (final one for the FED)
  calcCRC(dataWord&0xFFFFFFFF0000FFFF);   // EXCLUYENDO LOS VALORES DE CRC

      if ( ((dataWord>>60)&0xF)!=0xA) {   // Bits 60-63 are 0xA (control)
        if ( debug_ )  edm::LogWarning("oglez_dtab7_unpacker")
                         << "Trailer word "
                         << std::hex << dataWord << std::dec
                         << " does not start with 0xA";
        return;
      }

      int evtLgth = ( dataWord >> 32 ) & 0xFFFFFF; // Bits 32-55 is the number of lines
      int CRC     = ( dataWord >> 16 ) & 0xFFFF;   // Bits 16-31 is the expected CRC

      if ( newCRC_ != CRC ) {
        if ( debug_ ) edm::LogWarning("oglez_dtab7_unpacker")
                        << "Calculated CRC "
                        << std::hex << newCRC_
                        << " differs from CRC in trailer "
                        << CRC << std::dec;
        std::cout<<"OG-ERROR: CRC does not match!!! "<<evtLgth<<" "<<lineCounter_<<" "<<std::hex<<CRC<<" "<<newCRC_<<" "<<std::dec<<std::endl;
        return;
      }

      if ( lineCounter_ != evtLgth ) {
        if ( debug_ ) edm::LogWarning("oglez_dtab7_unpacker")
                        << "Number of words read != event lenght "
                        << lineCounter_ << " " << evtLgth;
        std::cout<<"OG-ERROR: line counting does not match!!! "<<evtLgth<<" "<<lineCounter_<<" "<<std::hex<<CRC<<" "<<newCRC_<<" "<<std::dec<<std::endl;
        return;
      }

//      std::cout<<"OG-CHECK: "<<evtLgth<<" "<<lineCounter_<<" "
//               <<std::hex<<CRC<<" "<<newCRC_<<" "<<std::dec<<std::endl;

  // Everything was read here!
  return;
}

//-----------------------------------------------------------------------
void OglezDTAB7RawToDigi::readAB7PayLoad_hitWord (long dataWord,int fedno, int slot)
// Read the HIT information from the Payload of the AB7.
{
  int ioff=0;
  while (ioff<31) {  // Two hits in the word
    int hitinfo=(dataWord>>ioff);   // Information for the hit (bits 0-29)
    ioff+=30;

    // ORIG?    int tdcTime    = ( dataWord >> 32 ) & 0x3FFF; // positions  32 -> 45
    // ORIG?    int tdcChannel = ( dataWord >> 46 ) & 0x1F;   // positions  46 -> 50
//OLD    int tdcChannel = (hitinfo>>16)&0x1F;  // FIXME: ni puta idea de a qué corresponde esto en la estructura.
    // ORIG?    int tdcId      = ( dataWord >> 51 ) & 0x3;    // positions  51 -> 52
//OLD    int tdcId = (hitinfo>>21)&0x3;  // FIXME: ni puta idea de a qué corresponde esto en la estructura.
    // ORIG?    int link       = ( dataWord >> 53 ) & 0x7F;   // positions  53 -> 59
//OLD    int link = (hitinfo>>23)&0x7F;  // FIXME: ni puta idea de a qué corresponde esto en la estructura.
    // ORIG?    int kchannel   = ( dataWord >> 46 ) & 0x3FFF; // positions  46 -> 59
    // int kchannel = (hitinfo>>16)&0x3FFF;  // FIXME: ni puta idea de a qué corresponde esto en la estructura.
    // This is used to get the drift time for the digi...
    // Como me casca el codigo porque el kchannel parece ser enorme... intento
    // usar el CH_ID de la nueva version de Alvaro (ver abajo)

    // ORIG? No idea to what it corresponds???? int offset_t    =    ( dataWord >> 37 ) & 0x1FF;         // positions   5 -> 13
    //int offset_t = 0;  // FIXME...
    // I understand from Bilal's code that what it is called "offset" is in fact the BX in V3:
    int bx = (hitinfo>>5)&0xFFF; // positions   5 -> 16

    //OLD    int corrected_bx=bx-bxCounter_+421;
    //OLD if (corrected_bx>3563) corrected_bx-=3564;  // Is like this?

    if (doHexDumping_) std::cout<<"OGINFO-INFO: Dump + DIGI Info summary: BX: "<<bx<<std::endl;

    if (bx==0xFFF) {  // invalid (empty, not-needed because of odd-number) hits
      //if ( debug_ ) edm::LogWarning("oglez_dtab7_unpacker")
      //              << "Empty/invalid hit information because of an odd number of hits in the slot ";
      continue;
    }

    // ORIG? int tdc_hit_t   =    ( dataWord >> 32    & 0x1F );        // positions   0 ->  4
    // I understand from Bilal's code that this is the time within the BX, using 25ns/30 units in V3
    int tdc_hit_t = hitinfo&0x1F;   // Bits 0-4: TDC value

    // Identifying the channel

    int channelId = (hitinfo>>17)&0x1FF;   // Bits 17-25: Channel index
    int slId = (hitinfo>>26)&0x3;   // Bits 26-27: SL number
    int stationId = (hitinfo>>28)&0x3;   // Bits 28-29: Station number

    //std::cout<<"Information: "<<channelId<<" "<<stationId<<" "<<slId<<std::endl;

    // Getting the channel index.
//    int dduId = theDDU(fedno, slot, link);
//    int rosId = theROS(fedno, slot, link);
//    int robId = theROB(fedno, slot, link);
//    DTROChainCoding channelIndex(dduId, rosId, robId, tdcId, tdcChannel);
//    uint32_t chCode=channelIndex.getCode();

    // Using the one for the Phase-2 SX5 proptype
    uint32_t chCode = channelId;  //

    if (hitOrder_.find(chCode) == hitOrder_.end()) hitOrder_[chCode] = 0;
    else hitOrder_[chCode]++;

//Estupidez, ver abajo formula de relacion    int layerId[4] = {4,2,3,1};

    // Forcing the chamber to be in an specific wheel & sector
    // if needed change also in the if (selector2==0) part of code
    //OLD int wheelId=0, sectorId=4; //, stationId=1, sectorId=4, slId=1;  // MB1 negative
    int wheelId=2, sectorId=12; //, stationId=1, sectorId=4, slId=1;  // MB1 negative
    stationId=2;
    //OLD    slId=3;
    //OLD int wheelId=0, stationId=1, sectorId=10, slId=1; // MB1 positive

    // Getting the values from the current mapping
    int layerId=-999, wire=-999;
    OglezDTAB7RawToDigi::currentSX5ChannelMapping(channelId,&slId,&layerId,&wire);  // slId needs to be overwritten!

//          //if (slot == 11) slId = 4-slId;
//      //NI IDEA    if (slot == 9) slId = 4-slId;  // slot =9 ==> SL3 ==> SL Phi2
//
//          // if needed change also in the if (selector2==0) part of code
//          int firstConnector=1;
//          int cellOffset=4*(firstConnector-1);
//          // 0 if first cable on the first FE connector (starting from the left)
//          // 4 if first cable on secord connector
//          // ...
//          // 4x(n-1)
// ALL THE CRAP BEFORE IS USELESS... CONBSIDERING THE ID structure of the cell is completely "to be decided"
//    DTWireId detId = DTWireId(wheelId, stationId, sectorId, slId, layerId[channelId%4], (channelId/4)+1+cellOffset);

    // Clean way:
    DTWireId detId = DTWireId(wheelId, stationId, sectorId, slId, layerId, wire);
//Useless... stored as value!    int wire = detId.wire();

//    std::cout<<"        Estamos en el canal "<<channelId<<" que se corresponde con "<<slId<<" "<<layerId<<" "
//             <<wire<<" [wire="<<detId.wire()<<"]"<<std::endl;

    // Storing the digi:
    // This is the way Mara and Andrea saved the time... but it seems very strange.
    // DTDigi digi(wire,  offset_t*30+tdc_hit_t-1, hitOrder_[channelIndex.getCode()]);
    // No idea how to interpret the "time-based" constructor.

    // This is what I think is likely correct:
    //    DTDigi digi(wire, corrected_bx*30+tdc_hit_t-1, hitOrder_[channelIndex.getCode()]);

    // For tests... storing information (in 1/25 ns units... according to Cristina's equation)
    // Although I discovered that to use the legacy we should store "legacy TDC counts"
    // There is a *25 to convert from "BX units" to ns
    // There is a *(32/25) to convert from ns to TDCCounts  (and it needs to be positive!)
    // We need to subtract 1 in the tdc_hit, because it goes from 1-30, due to some
    //            convention (Alvaro indicated so)
    int tdccounts = int(32*(bx+(tdc_hit_t-1)/30.)+0.5)-32*bxCounter_;
    while (tdccounts<0) tdccounts+=114048;// 32*3564;

    DTDigi digi(wire,tdccounts, hitOrder_[chCode]);
    //OLD DTDigi digi(wire,int(32*(bx+tdc_hit_t/30.)+0.5)-32*bxCounter_, hitOrder_[chCode]);

    //std::cout<<"            OGINFO: Inserting DIGI! "<<detId.layerId()<<" "<<wire<<" "<<channelIndex.getCode()<<" "<<hitOrder_[channelIndex.getCode()]<<std::endl;
    //std::cout<<"     DIGI TIME: "<<bx<<" "<<(bx-bxCounter_)+tdc_hit_t/30.<<" "<<digi.time()<<std::endl;
    // To convert digi.time() en "BX" hay que dividir por 25, sumarle el BX_fed
    // y si lo que sale es mayor que 3564, restarle ese número.

    digis_->insertDigi(detId.layerId(),digi);
  }  // While for the two hits!
}

//-----------------------------------------------------------------------
void OglezDTAB7RawToDigi::readAB7PayLoad_triggerPrimitive (long firstWord,long secondWord,int fedno, int slot)
// Read the Trigger Primitive information from the Payload of the AB7.
{
  int station = ((firstWord>>60)&0x3);   // Bits 60-61 (first word) is the station (-1 (?))
  int superlayer = ((firstWord>>58)&0x3);   // Bits 58-59 (first word) is the superlayer (1-3) or a phi-primitive (0)

  int quality = ((firstWord>>35)&0x3F);   // Bits 35-40 (first word) is the quality of the TP

//OLD not sure what to do with it  int t_fine = ((firstWord>>41)&0x1F);   // Bits 41-45 (first word) is the time fine (according to Alvaro's notation)
  int time = ((firstWord>>41)&0x1FFFF) ; // Bits 41-57 (first word) is the time (in nanoseconds)
  int bx = time/25 - bxCounter_;  // Bunch crossing in LHC notation (as I understand the "0")
//NOT SURE WHAT TO DO WITH THIS  int t_fine = time-25*bx;

  //std::cout<<"HOLA-BX: "<<bx<<" (TP)"<<std::endl;

  int slope = ((firstWord>>16)&0x7FFF);   // Bits 16-30 (first word) is the slope (phi or theta depending on SL)
  int position = ((firstWord)&0xFFFF);   // Bits 0-15 (first word) is the position (phi or theta depending on SL)

  int chi2 = (secondWord&0x1F);             // Bits 0-4 (second word) is the chi2 of the fit
  int tpindex = ((secondWord>>5)&0x07);     // Bits 5-7 (second word) is the index of this trigger primitive in the event

  if (0 || doHexDumping_) std::cout<<"OGINFO-INFO: Dump + TP Info summary: BX: "<<bx<<std::endl; // TEST:" "<<32*bx/25<<std::endl;

//  std::cout<<"            PRUEBA-TP(1): "<<bx<<" "<<station<<" "<<superlayer<<" "<<quality<<" "<<slope<<" "<<position<<std::endl;
//  std::cout<<"            PRUEBA-TP(2): "<<chi2<<" "<<tpindex<<std::endl;

  // Getting the hits per layer:
  for (int i=0;i<4;++i) {  // Scanning layer number (inverted with respect to usual? Alvaro claims 4 is bottom...)
    int lat = ((secondWord>>(11-i))&1);  // Laterality (1 is right, 0 is left)
    int use = ((secondWord>>(15-i))&1);  // Was digi used in the primitive?
    int driftime = ((secondWord>>(28-4*i))&0xF);  // 1/32 of Drift time in ns (?), 15 means greater than 14, calculated as hit_time-primitive time
    int chan = ((secondWord>>(53-7*i))&0x3F);  // channel number inside its layer
//    std::cout<<"                 - Layer: "<<(i+1)<<" "<<lat<<" "<<use<<" "<<driftime<<" "<<chan<<std::endl;
  }

  // Some values seem to be "superlayer dependent"... now it is only phi.
  int phiAngle=position*250;    // Values as Camilo requested: position in um (original is in mm/4)

  int phiBending=slope;    // Acording to Jose Manuel (email 2019_04_03) this is 4096*tan(phi)
  // With Camilo we agree on passing this 4096*tan(phi) BUT we need to account for the sign.
  if ( ((slope>>14)&0x01)==0x01) { // Negative value!
    phiBending=(slope-32768);   // (1<<15)    // -1*(((~slope)&0x7FFF)+1);
  }

  //  std::cout<<"HOLA "<<std::hex<<slope<<" "<<phiBending<<" "<<std::dec<<" "<<phiBending<<" "<<(-1*(((~slope)&0x7FFF)+1))<<std::endl;

  int zCoordinate=position;   // I agree to store the original values in this fields.
  int zSlope=slope;

  // Adding the trigger primitive
  L1Phase2MuDTPhDigi trigprim(bx,   // ubx (m_bx)
			      2,   // uwh (m_wheel)     // FIXME: It is not clear who provides this?
			      12,   // usc (m_sector)    // FIXME: It is not clear who provides this?
			      2,
			      // CB line above is a little hack station,   // ust (m_station)
			      0, // the SL

			      phiAngle,   // uphi (_phiAngle)
			      phiBending,   // uphib (m_phiBending)
			      
			      quality,  // uqua (m_qualityCode)
			      tpindex,  // uind (m_segmentIndex)
			      time,  // ut0 (m_t0Segment)
			      chi2,  // uchi2 (m_chi2Segment)
			      -10);  // urpc (m_rpcFlag)
  primitives_.push_back(trigprim);
}

//-----------------------------------------------------------------------
int OglezDTAB7RawToDigi::theDDU(int crate, int slot, int link) {

     if (crate == 1368) {
       if (slot < 7) return 770;
       return 771;
     }

     if (crate == 1370) {
       if (slot > 6) return 773;
       return 774;
     }

     return 772;
}

//-----------------------------------------------------------------------
int OglezDTAB7RawToDigi::theROS(int crate, int slot, int link) {

  if (slot%6 == 5) return link+1;

  int ros = (link/24) + 3*(slot%6) - 2;
  return ros;
}

//-----------------------------------------------------------------------
int OglezDTAB7RawToDigi::theROB(int crate, int slot, int link) {

  if (slot%6 == 5) return 23;

  int rob = link%24;
  if (rob < 15) return rob;
  if (rob == 15) return 24;
  return rob-1;
}

//-----------------------------------------------------------------------
void OglezDTAB7RawToDigi::calcCRC(long word) {
// Routine to compute the CRC using a new word.

  int myCRC[16], D[64], C[16];

  for ( int i = 0; i < 64; ++i ) { D[i]    = (word >> i) & 0x1; }
  for ( int i = 0; i < 16; ++i ) { C[i]    = (newCRC_>>i)  & 0x1; }

  myCRC[0] = ( D[63] + D[62] + D[61] + D[60] + D[55] + D[54] +
               D[53] + D[52] + D[51] + D[50] + D[49] + D[48] +
               D[47] + D[46] + D[45] + D[43] + D[41] + D[40] +
               D[39] + D[38] + D[37] + D[36] + D[35] + D[34] +
               D[33] + D[32] + D[31] + D[30] + D[27] + D[26] +
               D[25] + D[24] + D[23] + D[22] + D[21] + D[20] +
               D[19] + D[18] + D[17] + D[16] + D[15] + D[13] +
               D[12] + D[11] + D[10] + D[9]  + D[8]  + D[7]  +
               D[6]  + D[5]  + D[4]  + D[3]  + D[2]  + D[1]  +
               D[0]  + C[0]  + C[1]  + C[2]  + C[3]  + C[4]  +
               C[5]  + C[6]  + C[7]  + C[12] + C[13] + C[14] +
               C[15] )%2;

  myCRC[1] = ( D[63] + D[62] + D[61] + D[56] + D[55] + D[54] +
               D[53] + D[52] + D[51] + D[50] + D[49] + D[48] +
               D[47] + D[46] + D[44] + D[42] + D[41] + D[40] +
               D[39] + D[38] + D[37] + D[36] + D[35] + D[34] +
               D[33] + D[32] + D[31] + D[28] + D[27] + D[26] +
               D[25] + D[24] + D[23] + D[22] + D[21] + D[20] +
               D[19] + D[18] + D[17] + D[16] + D[14] + D[13] +
               D[12] + D[11] + D[10] + D[9]  + D[8]  + D[7]  +
               D[6]  + D[5]  + D[4]  + D[3]  + D[2]  + D[1]  +
               C[0]  + C[1]  + C[2]  + C[3]  + C[4]  + C[5]  +
               C[6]  + C[7]  + C[8]  + C[13] + C[14] + C[15] )%2;

  myCRC[2] = ( D[61] + D[60] + D[57] + D[56] + D[46] + D[42] +
               D[31] + D[30] + D[29] + D[28] + D[16] + D[14] +
               D[1]  + D[0]  + C[8]  + C[9]  + C[12] + C[13] )%2;

  myCRC[3] = ( D[62] + D[61] + D[58] + D[57] + D[47] + D[43] +
               D[32] + D[31] + D[30] + D[29] + D[17] + D[15] +
               D[2]  + D[1]  + C[9]  + C[10] + C[13] + C[14] )%2;

  myCRC[4] = ( D[63] + D[62] + D[59] + D[58] + D[48] + D[44] +
               D[33] + D[32] + D[31] + D[30] + D[18] + D[16] +
               D[3]  + D[2]  + C[0]  + C[10] + C[11] + C[14] +
               C[15] )%2;

  myCRC[5] = ( D[63] + D[60] + D[59] + D[49] + D[45] + D[34] +
               D[33] + D[32] + D[31] + D[19] + D[17] + D[4]  +
               D[3]  + C[1]  + C[11] + C[12] + C[15] )%2;

  myCRC[6] = ( D[61] + D[60] + D[50] + D[46] + D[35] + D[34] +
               D[33] + D[32] + D[20] + D[18] + D[5]  + D[4]  +
               C[2]  + C[12] + C[13] )%2;

  myCRC[7] = ( D[62] + D[61] + D[51] + D[47] + D[36] + D[35] +
               D[34] + D[33] + D[21] + D[19] + D[6]  + D[5]  +
               C[3]  + C[13] + C[14] )%2;

  myCRC[8] = ( D[63] + D[62] + D[52] + D[48] + D[37] + D[36] +
               D[35] + D[34] + D[22] + D[20] + D[7]  + D[6]  +
               C[0]  + C[4]  + C[14] + C[15] )%2;

  myCRC[9] = ( D[63] + D[53] + D[49] + D[38] + D[37] + D[36] +
               D[35] + D[23] + D[21] + D[8]  + D[7]  + C[1]  +
               C[5]  + C[15] )%2;

  myCRC[10] = ( D[54] + D[50] + D[39] + D[38] + D[37] + D[36] +
                D[24] + D[22] + D[9]  + D[8]  + C[2]  + C[6] )%2;

  myCRC[11] = ( D[55] + D[51] + D[40] + D[39] + D[38] + D[37] +
                D[25] + D[23] + D[10] + D[9]  + C[3]  + C[7] )%2;

  myCRC[12] = ( D[56] + D[52] + D[41] + D[40] + D[39] + D[38] +
                D[26] + D[24] + D[11] + D[10] + C[4]  + C[8] )%2;

  myCRC[13] = ( D[57] + D[53] + D[42] + D[41] + D[40] + D[39] +
                D[27] + D[25] + D[12] + D[11] + C[5]  + C[9] )%2;

  myCRC[14] = ( D[58] + D[54] + D[43] + D[42] + D[41] + D[40] +
                D[28] + D[26] + D[13] + D[12] + C[6]  + C[10] )%2;

  myCRC[15] = ( D[63] + D[62] + D[61] + D[60] + D[59] + D[54] +
                D[53] + D[52] + D[51] + D[50] + D[49] + D[48] +
                D[47] + D[46] + D[45] + D[44] + D[42] + D[40] +
                D[39] + D[38] + D[37] + D[36] + D[35] + D[34] +
                D[33] + D[32] + D[31] + D[30] + D[29] + D[26] +
                D[25] + D[24] + D[23] + D[22] + D[21] + D[20] +
                D[19] + D[18] + D[17] + D[16] + D[15] + D[14] +
                D[12] + D[11] + D[10] + D[9]  + D[8]  + D[7]  +
                D[6]  + D[5]  + D[4]  + D[3]  + D[2]  + D[1]  +
                D[0]  + C[0]  + C[1]  + C[2]  + C[3]  + C[4]  +
                C[5]  + C[6]  + C[11] + C[12] + C[13] + C[14] +
                C[15] )%2;

  int tempC = 0x0;
  for ( int i = 0; i < 16 ; ++i) { tempC = tempC + ( myCRC[i] << i ); }

  newCRC_ = tempC;
  return;
}

//-----------------------------------------------------------------------
void OglezDTAB7RawToDigi::currentSX5ChannelMapping (int chanid, int *sl, int *layer, int *wire)
  // Static routine to process the current mapping between the channel ID and
  // the wire, Layer and SL of the digi for the Phase-2 prototype..
{
  // Not clear the logic)
  if (chanid>=64 && chanid<=79) {
    (*sl) = 1;
    (*wire) = 1 + ((chanid)-64)/4;  // wire 1 to 4
  }
  else if (chanid>=80 && chanid<=95) {
    (*sl) = 3;
    (*wire) = 17 + ((chanid)-80)/4;  // wire 17 to 20
  }
  else if (chanid>=112 && chanid<=127) {
    (*sl) = 1;
    (*wire) = 17 + ((chanid)-112)/4;  // wire 17 to 20
  }
  else {
    std::cout<<"OG-ERROR: Channel-ID does not seem to be correct in the currentSX5ChannelMapping: "<<chanid<<std::endl;
    (*sl) = -1;
    (*wire) = -999;
  }

  // Layer is simpler (If I got it right)
  int i = (chanid%4);
  (*layer) = 4 - 2*(i%2) - (i>=2);
}

//-----------------------------------------------------------------------
// #include "FWCore/Framework/interface/MakerMacros.h"
// DEFINE_FWK_MODULE(OglezDTAB7RawToDigi);
//=======================================================================
