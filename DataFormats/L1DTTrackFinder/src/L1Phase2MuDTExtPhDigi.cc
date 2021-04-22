//-------------------------------------------------
//
//   Class L1MuDTChambPhDigi
//
//   Description: trigger primtive data for the
//                muon barrel Phase2 trigger
//
//
//   Author List: Federica Primavera  Bologna INFN
//
//
//--------------------------------------------------

//-----------------------
// This Class's Header --
//-----------------------
#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTExtPhDigi.h"

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
L1Phase2MuDTExtPhDigi::L1Phase2MuDTExtPhDigi():

  // m_bx(-100),
  // m_wheel(0),
  // m_sector(0),
  // m_station(0),
  // m_superlayer(0),
  // m_phiAngle(0),
  // m_phiBending(0),
  // m_qualityCode(-1),
  // m_index(0),
  // m_t0(0),
  // m_chi2(0),
  // m_rpcFlag(-10){} 

  L1Phase2MuDTPhDigi(),

  m_xLocal(0),
  m_tanPsi(0) {
  
  for (int i=0; i<8; i++) {
    m_pathWireId[i] = -1;
    m_pathTDC[i] = -1;
    m_pathLat[i] = 2;
  }
}

L1Phase2MuDTExtPhDigi::L1Phase2MuDTExtPhDigi(int bx, int wh, int sc, int st, int sl, int phi, int phib, int qual, int idx, int t0, int chi2, int x, int tan, int rpc, int wireId[8], int tdc[8],int lat[8]): 

  // m_bx(bx),
  // m_wheel(wh),
  // m_sector(sc),
  // m_station(st),
  // m_superlayer(sl),
  // m_phiAngle(phi),
  // m_phiBending(phib),
  // m_qualityCode(qual),
  // m_index(idx),
  // m_t0(t0),
  // m_chi2(chi2),
  // m_rpcFlag(rpc), 

  L1Phase2MuDTPhDigi(bx, wh, sc, st, sl, phi, phib, qual, idx, t0, chi2, rpc),
  
  m_xLocal(x),
  m_tanPsi(tan) {
  
  for (int i=0; i<8; i++) {
    m_pathWireId[i] = wireId[i];
    m_pathTDC[i] = tdc[i];
    m_pathLat[i] = lat[i];
  }
  }

//--------------
// Operations --
//--------------
// int L1Phase2MuDTExtPhDigi::bxNum() const { return m_bx; }

// int L1Phase2MuDTExtPhDigi::whNum() const { return m_wheel; }

// int L1Phase2MuDTExtPhDigi::scNum() const { return m_sector; }

// int L1Phase2MuDTExtPhDigi::stNum() const { return m_station; }

// int L1Phase2MuDTExtPhDigi::slNum() const { return m_superlayer; }

// int L1Phase2MuDTExtPhDigi::phi() const { return m_phiAngle; }

// int L1Phase2MuDTExtPhDigi::phiBend() const { return m_phiBending; }

// int L1Phase2MuDTExtPhDigi::quality() const { return m_qualityCode; }

// int L1Phase2MuDTExtPhDigi::index() const { return m_index; }

// int L1Phase2MuDTExtPhDigi::t0() const { return m_t0; }

// int L1Phase2MuDTExtPhDigi::chi2() const { return m_chi2; }

int L1Phase2MuDTExtPhDigi::xLocal() const { return m_xLocal; }

int L1Phase2MuDTExtPhDigi::tanPsi() const { return m_tanPsi; } 

// int L1Phase2MuDTExtPhDigi::rpcFlag() const { return m_rpcFlag; }

int L1Phase2MuDTExtPhDigi::pathWireId(int i) const { return m_pathWireId[i]; }

int L1Phase2MuDTExtPhDigi::pathTDC(int i) const { return m_pathTDC[i]; }

int L1Phase2MuDTExtPhDigi::pathLat(int i) const { return m_pathLat[i]; }
