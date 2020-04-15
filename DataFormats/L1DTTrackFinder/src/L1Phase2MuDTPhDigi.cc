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
#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTPhDigi.h"

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
L1Phase2MuDTPhDigi::L1Phase2MuDTPhDigi()
    : m_bx(-100),
      m_wheel(0),
      m_sector(0),
      m_station(0),
      m_superlayer(0),
      m_phiAngle(0),
      m_phiBending(0),
      m_qualityCode(-1),
      m_index(0),
      m_t0(0),
      m_chi2(0),
      m_rpcFlag(-10) {
  for (int i=0; i<8; i++) {
    m_pathWireId[i] = -1;
    m_pathTDC[i] = -1;
    m_pathLat[i] = 2;
  }
}

L1Phase2MuDTPhDigi::L1Phase2MuDTPhDigi(int bx, int wh, int sc, int st, int sl, int phi, int phib, int qual, int idx, int t0, int chi2, int rpc, int wireId[8], int tdc[8],int lat[8])
    : m_bx(bx),
      m_wheel(wh),
      m_sector(sc),
      m_station(st),
      m_superlayer(sl),
      m_phiAngle(phi),
      m_phiBending(phib),
      m_qualityCode(qual),
      m_index(idx),
      m_t0(t0),
      m_chi2(chi2),
      m_rpcFlag(rpc) {

  for (int i=0; i<8; i++) {
    m_pathWireId[i] = wireId[i];
    m_pathTDC[i] = tdc[i];
    m_pathLat[i] = lat[i];
  }

}

//--------------
// Operations --
//--------------
int L1Phase2MuDTPhDigi::bxNum() const { return m_bx; }

int L1Phase2MuDTPhDigi::whNum() const { return m_wheel; }

int L1Phase2MuDTPhDigi::scNum() const { return m_sector; }

int L1Phase2MuDTPhDigi::stNum() const { return m_station; }

int L1Phase2MuDTPhDigi::slNum() const { return m_superlayer; }

int L1Phase2MuDTPhDigi::phi() const { return m_phiAngle; }

int L1Phase2MuDTPhDigi::phiBend() const { return m_phiBending; }

int L1Phase2MuDTPhDigi::quality() const { return m_qualityCode; }

int L1Phase2MuDTPhDigi::index() const { return m_index; }

int L1Phase2MuDTPhDigi::t0() const { return m_t0; }

int L1Phase2MuDTPhDigi::chi2() const { return m_chi2; }

int L1Phase2MuDTPhDigi::rpcFlag() const { return m_rpcFlag; }

int L1Phase2MuDTPhDigi::pathWireId(int i) const { return m_pathWireId[i]; }

int L1Phase2MuDTPhDigi::pathTDC(int i) const { return m_pathTDC[i]; }

int L1Phase2MuDTPhDigi::pathLat(int i) const { return m_pathLat[i]; }
