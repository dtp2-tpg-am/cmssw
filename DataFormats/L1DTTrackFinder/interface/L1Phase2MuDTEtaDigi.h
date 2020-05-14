#ifndef L1Phase2MuDTEtaDigi_H
#define L1Phase2MuDTEtaDigi_H

class L1Phase2MuDTEtaDigi {
public:
  //  Constructors
  L1Phase2MuDTEtaDigi();

  L1Phase2MuDTEtaDigi(
      int bx, int wh, int sc, int st, int sl, int eta, int etab, int qual, int idx, int t0, int chi2, int rpc = -10);

  // Operations
  int bxNum() const;

  int whNum() const;
  int scNum() const;
  int stNum() const;
  int slNum() const;

  int eta() const;
  int etaBend() const;

  int quality() const;
  int index() const;

  int t0() const;
  int chi2() const;

  int rpcFlag() const;

private:
  int m_bx;
  int m_wheel;
  int m_sector;
  int m_station;
  int m_superlayer;

  int m_etaAngle;
  int m_etaBending;

  int m_qualityCode;
  int m_index;

  int m_t0;
  int m_chi2;

  int m_rpcFlag;
};

#endif
