#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTEtaDigi.h"

L1Phase2MuDTEtaDigi::L1Phase2MuDTEtaDigi()
    : m_bx(-100),
      m_wheel(0),
      m_sector(0),
      m_station(0),
      m_superlayer(0),
      m_etaAngle(0),
      m_etaBending(0),
      m_qualityCode(-1),
      m_index(0),
      m_t0(0),
      m_chi2(0),
      m_rpcFlag(-10) {}

L1Phase2MuDTEtaDigi::L1Phase2MuDTEtaDigi(
    int bx, int wh, int sc, int st, int sl, int eta, int etab, int qual, int idx, int t0, int chi2, int rpc)
    : m_bx(bx),
      m_wheel(wh),
      m_sector(sc),
      m_station(st),
      m_superlayer(sl),
      m_etaAngle(eta),
      m_etaBending(etab),
      m_qualityCode(qual),
      m_index(idx),
      m_t0(t0),
      m_chi2(chi2),
      m_rpcFlag(rpc) {}

int L1Phase2MuDTEtaDigi::bxNum() const { return m_bx; }

int L1Phase2MuDTEtaDigi::whNum() const { return m_wheel; }

int L1Phase2MuDTEtaDigi::scNum() const { return m_sector; }

int L1Phase2MuDTEtaDigi::stNum() const { return m_station; }

int L1Phase2MuDTEtaDigi::slNum() const { return m_superlayer; }

int L1Phase2MuDTEtaDigi::eta() const { return m_etaAngle; }

int L1Phase2MuDTEtaDigi::etaBend() const { return m_etaBending; }

int L1Phase2MuDTEtaDigi::quality() const { return m_qualityCode; }

int L1Phase2MuDTEtaDigi::index() const { return m_index; }

int L1Phase2MuDTEtaDigi::t0() const { return m_t0; }

int L1Phase2MuDTEtaDigi::chi2() const { return m_chi2; }

int L1Phase2MuDTEtaDigi::rpcFlag() const { return m_rpcFlag; }
