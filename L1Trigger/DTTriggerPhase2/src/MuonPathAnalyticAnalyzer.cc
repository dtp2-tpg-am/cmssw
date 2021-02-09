#include "L1Trigger/DTTriggerPhase2/interface/MuonPathAnalyticAnalyzer.h"
#include <cmath>
#include <memory>

using namespace edm;
using namespace std;
using namespace cmsdt;
// ============================================================================
// Constructors and destructor
// ============================================================================
MuonPathAnalyticAnalyzer::MuonPathAnalyticAnalyzer(const ParameterSet &pset, edm::ConsumesCollector &iC)
    : MuonPathAnalyzer(pset, iC),
      bxTolerance_(30),
      minQuality_(LOWQGHOST),
      chiSquareThreshold_(50),
      debug_(pset.getUntrackedParameter<bool>("debug")),
      chi2Th_(pset.getUntrackedParameter<double>("chi2Th")),
      tanPhiTh_(pset.getUntrackedParameter<double>("tanPhiTh")),
      use_LSB_(pset.getUntrackedParameter<bool>("use_LSB")),
      tanPsi_precision_(pset.getUntrackedParameter<double>("tanPsi_precision")),
      x_precision_(pset.getUntrackedParameter<double>("x_precision")) {
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "MuonPathAnalyzer: constructor";

  // setChiSquareThreshold(chi2Th_ * 100.);
  setChiSquareThreshold(chi2Th_);
  fillLAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER();
  fillLAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL();

  //shift
  int rawId;
  shift_filename_ = pset.getParameter<edm::FileInPath>("shift_filename");
  std::ifstream ifin3(shift_filename_.fullPath());
  double shift;
  if (ifin3.fail()) {
    throw cms::Exception("Missing Input File")
        << "MuonPathAnalyticAnalyzer::MuonPathAnalyticAnalyzer() -  Cannot find " << shift_filename_.fullPath();
  }
  while (ifin3.good()) {
    ifin3 >> rawId >> shift;
    shiftinfo_[rawId] = shift;
  }

  chosen_sl_ = pset.getUntrackedParameter<int>("trigger_with_sl");

  if (chosen_sl_ != 1 && chosen_sl_ != 3 && chosen_sl_ != 4) {
    LogDebug("MuonPathAnalyticAnalyzer") << "chosen sl must be 1,3 or 4(both superlayers)";
    assert(chosen_sl_ != 1 && chosen_sl_ != 3 && chosen_sl_ != 4);  //4 means run using the two superlayers
  }

  dtGeomH = iC.esConsumes<DTGeometry, MuonGeometryRecord, edm::Transition::BeginRun>();
}

MuonPathAnalyticAnalyzer::~MuonPathAnalyticAnalyzer() {
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "MuonPathAnalyzer: destructor";
}

// ============================================================================
// Main methods (initialise, run, finish)
// ============================================================================
void MuonPathAnalyticAnalyzer::initialise(const edm::EventSetup &iEventSetup) {
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "MuonPathAnalyticAnalyzer::initialiase";

  const MuonGeometryRecord &geom = iEventSetup.get<MuonGeometryRecord>();
  dtGeo_ = &geom.get(dtGeomH);
}

void MuonPathAnalyticAnalyzer::run(edm::Event &iEvent,
                                const edm::EventSetup &iEventSetup,
                                MuonPathPtrs &muonpaths,
                                std::vector<metaPrimitive> &metaPrimitives) {
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "MuonPathAnalyticAnalyzer: run";

  // fit per SL (need to allow for multiple outputs for a single mpath)
  for (auto &muonpath : muonpaths) {
    analyze(muonpath, metaPrimitives);
  }
}

void MuonPathAnalyticAnalyzer::finish() {
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "MuonPathAnalyzer: finish";
};

constexpr int MuonPathAnalyticAnalyzer::LAYER_ARRANGEMENTS_[NUM_LAYERS][NUM_CELL_COMB] = {
    {0, 1, 2},
    {1, 2, 3},  // Consecutive groups
    {0, 1, 3},
    {0, 2, 3}  // Non-consecutive groups
};

//------------------------------------------------------------------
//--- Métodos privados
//------------------------------------------------------------------

void MuonPathAnalyticAnalyzer::analyze(MuonPathPtr &inMPath, std::vector<metaPrimitive> &metaPrimitives) {
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "DTp2:analyze \t\t\t\t starts";

  // LOCATE MPATH
  int selected_Id = 0;
  if (inMPath->primitive(0)->tdcTimeStamp() != -1)
    selected_Id = inMPath->primitive(0)->cameraId();
  else if (inMPath->primitive(1)->tdcTimeStamp() != -1)
    selected_Id = inMPath->primitive(1)->cameraId();
  else if (inMPath->primitive(2)->tdcTimeStamp() != -1)
    selected_Id = inMPath->primitive(2)->cameraId();
  else if (inMPath->primitive(3)->tdcTimeStamp() != -1)
    selected_Id = inMPath->primitive(3)->cameraId();

  DTLayerId thisLId(selected_Id);
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "Building up MuonPathSLId from rawId in the Primitive";
  DTSuperLayerId MuonPathSLId(thisLId.wheel(), thisLId.station(), thisLId.sector(), thisLId.superLayer());
  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer") << "The MuonPathSLId is" << MuonPathSLId;

  if (debug_)
    LogDebug("MuonPathAnalyticAnalyzer")
        << "DTp2:analyze \t\t\t\t In analyze function checking if inMPath->isAnalyzable() " << inMPath->isAnalyzable();


  if (chosen_sl_ < 4 && thisLId.superLayer() != chosen_sl_)
    return;  // avoid running when mpath not in chosen SL (for 1SL fitting)

  auto mPath = std::make_shared<MuonPath>(inMPath);
  mPath->setQuality(NOPATH);

  int wi[4], wires[4], t0s[4], valids[4];
  // bool is_four_hit = true;
  for (int j = 0; j < NUM_LAYERS; j++) {
    if (mPath->primitive(j)->isValidTime()){
      wi[j] = mPath->primitive(j)->channelId();
      wires[j] = mPath->primitive(j)->channelId();
      t0s[j] = mPath->primitive(j)->tdcTimeStamp();
      valids[j] = 1; 
    } else {
      wi[j] = -1;
      wires[j] = -1;
      t0s[j] = -1;
      valids[j] = 0;
      // is_four_hit = false;
    }
  }

  if (wi[0] < 0) wi[0] = wi[1];
  else if (wi[1] < 0) wi[1] = wi[0];
  else if (wi[2] < 0) wi[2] = wi[1] - 1;
  else if (wi[3] < 0) wi[3] = wi[2];

  int cell_horiz_layout[4];
  for (int lay = 0; lay < NUM_LAYERS; lay++) {
    cell_horiz_layout[lay] = (wi[lay] - wi[0]) * 2;
    if (lay % 2 != 0) cell_horiz_layout[lay]--;
  }

  // calculate the coarse offset position
  int tmp = 1;
  if (valids[1] == 0) tmp = 3;
  int coarse_pos = (wi[tmp] * 2 - cell_horiz_layout[tmp]) * 21 * std::pow(2, 4);

  //calculate the relative position of wires in mm wrt layer 0's cell wire
  int xwire_mm[4];
  for (int lay = 0; lay < NUM_LAYERS; lay++) {
    xwire_mm[lay] = 21 * cell_horiz_layout[lay];
  }

  // divide the timestamps in coarse + reduced part
  int valid_coarse_times[4], min_coarse_time = 999999, max_coarse_time = -999999;
  for (int lay = 0; lay < NUM_LAYERS; lay++) {
    if (valids[lay] == 1) {
      valid_coarse_times[lay] = (t0s[lay] >> (TDCTIME_REDUCED_SIZE - 1));
      if (valid_coarse_times[lay] < min_coarse_time) {
        min_coarse_time = valid_coarse_times[lay];
      }
      if (valid_coarse_times[lay] > max_coarse_time) {
        max_coarse_time = valid_coarse_times[lay];
      }
    } else {
      valid_coarse_times[lay] = -1;
    }
  }
  
  if (max_coarse_time - min_coarse_time >= 2) return; 
  int coarse_offset = max_coarse_time - 1;

  int reduced_times[4];
  for (int lay = 0; lay < NUM_LAYERS; lay++) {
    reduced_times[lay] = ( ( 1 - ((max_coarse_time & 1) ^ ((t0s[lay] >> (TDCTIME_REDUCED_SIZE - 1)) & 1))	) << (TDCTIME_REDUCED_SIZE - 1));
    reduced_times[lay] += (t0s[lay] & std::stoi(std::string(TDCTIME_REDUCED_SIZE - 1, '1'), nullptr, 2) );
  }
  std::vector <LATCOMB_CONSTANTS> latcomb_consts_arr;
  std::vector <LATCOMB_CONSTANTS_F> latcomb_consts_arr_real;
  // for (auto & elem : LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER) 
  for (size_t i = 0; i < LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.size(); i++) { 
    auto elem = LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER[i];
    if (elem.cell_valid_layout.valid[0] == valids[0] 
        && elem.cell_valid_layout.valid[1] == valids[1]
        && elem.cell_valid_layout.valid[2] == valids[2]
        && elem.cell_valid_layout.valid[3] == valids[3]
        && elem.cell_valid_layout.cell_horiz_layout[0] == cell_horiz_layout[0]
        && elem.cell_valid_layout.cell_horiz_layout[1] == cell_horiz_layout[1]
        && elem.cell_valid_layout.cell_horiz_layout[2] == cell_horiz_layout[2]
        && elem.cell_valid_layout.cell_horiz_layout[3] == cell_horiz_layout[3]) {
      for (auto & ind_latcomb_consts : elem.latcomb_constants)
        latcomb_consts_arr.push_back(ind_latcomb_consts);
      auto elem_real = LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL[i];
      for (auto & ind_latcomb_consts : elem_real.latcomb_constants)
        latcomb_consts_arr_real.push_back(ind_latcomb_consts);
    }
  }
  // for (auto & latcomb_consts : latcomb_consts_arr) {
  for (size_t i = 0; i < latcomb_consts_arr.size(); i++) {
    auto latcomb_consts = latcomb_consts_arr[i];
    auto latcomb_consts_real = latcomb_consts_arr_real[i];
    segment_fitter(MuonPathSLId, wires, t0s, valids, reduced_times, cell_horiz_layout, latcomb_consts,
      latcomb_consts_real, xwire_mm, coarse_pos, coarse_offset, metaPrimitives);
  }
}

int MuonPathAnalyticAnalyzer::compute_parameter(MAGNITUDE constants, int t0s[4], int DIV_SHR_BITS, int INCREASED_RES) {
  long int result = 0; 
  for (int lay = 0; lay < NUM_LAYERS; lay++){
    result += constants.coeff[lay] * t0s[lay];
  }
  result = ((result * int(std::pow(2, INCREASED_RES)) + constants.add) * constants.mult) >> DIV_SHR_BITS;
  
  return result;
}

float MuonPathAnalyticAnalyzer::compute_exact_parameter(MAGNITUDE_F constants, int t0s[4]) {
  float result = 0; 
  for (int lay = 0; lay < NUM_LAYERS; lay++){
    result += constants.coeff[lay] * t0s[lay];
  }
  if (constants.mult != 0) {
    result = (result  + constants.add) / constants.mult;
  } else result = 0;
  return result;
}

void MuonPathAnalyticAnalyzer::segment_fitter(DTSuperLayerId MuonPathSLId, int wires[4], int t0s[4], int valid[4],  int reduced_times[4],
    int cell_horiz_layout[4], LATCOMB_CONSTANTS latcomb_consts, LATCOMB_CONSTANTS_F latcomb_consts_exact, int xwire_mm[4], int coarse_pos, int coarse_offset,
    std::vector<cmsdt::metaPrimitive> &metaPrimitives) {
  auto latcomb = latcomb_consts.latcomb;
  auto constants = latcomb_consts.constants;
  auto exact_constants = latcomb_consts_exact.constants;
  bool is_four_hit = true;

  if (latcomb == 0) return;
  
  int lat_array[4];
  for (int lay = 0; lay < NUM_LAYERS; lay++) {
    if (((latcomb >> lay) & 1) != 0) {
      lat_array[lay] = 1;
    } else lat_array[lay] = -1;
  }
   
  int time = compute_parameter(constants.t0, reduced_times, DIV_SHR_BITS_T0, INCREASED_RES_T0);
  int pos = compute_parameter(constants.pos, reduced_times, DIV_SHR_BITS_POS, INCREASED_RES_POS);
  int slope = compute_parameter(constants.slope, reduced_times, DIV_SHR_BITS_SLOPE, INCREASED_RES_SLOPE);
  int slope_xhh = compute_parameter(constants.slope_xhh, reduced_times, DIV_SHR_BITS_SLOPE_XHH, INCREASED_RES_SLOPE_XHH);
  
  float exact_time = compute_exact_parameter(exact_constants.t0, reduced_times);
  float exact_pos = compute_exact_parameter(exact_constants.pos, reduced_times);
  // float exact_slope = compute_exact_parameter(exact_constants.slope, reduced_times);
  float exact_slope_xhh = compute_exact_parameter(exact_constants.slope_xhh, reduced_times);
  
  int bx_time = time + (coarse_offset << (TDCTIME_REDUCED_SIZE - 1));
  // float bx_time_perfect = exact_time + (coarse_offset << (TDCTIME_REDUCED_SIZE - 1));
  
  // cout << "Time: " << time << " " << exact_time << " ";
  // cout << "Pos: " << pos / (10 * pow(2, INCREASED_RES_POS)) << " " << exact_pos << " ";
  // cout << "Slope: " << slope / std::pow(2, INCREASED_RES_SLOPE) << " " << exact_slope << " ";
  // cout << "Slope_xhh: " << slope_xhh << " " << exact_slope_xhh << endl;

  pos += coarse_pos;
  exact_pos += coarse_pos/16.;  

  int chi2_mm2_p = 0;
  float chi2_mm2_perfect = 0;

  for (int lay = 0; lay < NUM_LAYERS; lay++) {
    int drift_time = reduced_times[lay] - time;
    if (valid[lay] == 1 && (drift_time < 0 || drift_time > MAXDRIFT)) return;

    int drift_dist = (( (drift_time * int(pow(2, 4)) + 5) * 445 ) >> 13);
    int xdist = xwire_mm[lay] * pow(2, 4) - (pos - coarse_pos) + lat_array[lay] * drift_dist;
    xdist -= (3 - 2 * (3 - lay)) * slope_xhh;
    
    float xdist_perfect = xwire_mm[lay] - (exact_pos - coarse_pos / 16.) + lat_array[lay] * (reduced_times[lay] - exact_time) * DRIFT_SPEED_F;
    xdist_perfect -= (3 - 2 * (3 - lay)) * exact_slope_xhh;
    
    int res = xdist;
    float res_perfect = xdist_perfect;
    if (valid[lay] == 0) {
      res = 0;
      res_perfect = 0.;
      is_four_hit = false;
    }
    chi2_mm2_p += res * res * 4;
    chi2_mm2_perfect += res_perfect * res_perfect; 
  }
  
  int quality = 4;
  if (!is_four_hit) quality = 2;
  
  // Obtain coordinate values in floating point
  double pos_f, slope_f, chi2_f;
  DTWireId wireId(MuonPathSLId, 2, 1);
  pos_f = double(pos) + int(10 * shiftinfo_[wireId.rawId()] * std::pow(2, INCREASED_RES_POS)); // position in mm * precision in JM RF
  pos_f /=  (10 * pow(2, INCREASED_RES_POS)); // position in cm in JM RF 
  slope_f = - (double(slope) / std::pow(2, INCREASED_RES_SLOPE));
  // chi2_f = double(chi2_mm2_p) / (16. * 64. * 100.);
  cout << chi2_mm2_perfect / 100. << endl;
  chi2_f = chi2_mm2_perfect / 100.;
  
  // Impose the thresholds
  if (std::abs(slope_f) > tanPhiTh_) return;
  if (chi2_f > (chiSquareThreshold_)) return;

  // Compute phi and phib
  // Implemented using cmssw geometry as of now, will implemented fw-like in the near future
  DTChamberId ChId(MuonPathSLId.wheel(), MuonPathSLId.station(), MuonPathSLId.sector());
  double z = 0;
  double z1 = Z_POS_SL;
  double z3 = -1. * z1;
  if (ChId.station() == 3 or ChId.station() == 4) {
	z1 = z1 + Z_SHIFT_MB4;
	z3 = z3 + Z_SHIFT_MB4;
  } else if (MuonPathSLId.superLayer() == 1)
	z = z1;
  else if (MuonPathSLId.superLayer() == 3)
	z = z3;
  GlobalPoint jm_x_cmssw_global = dtGeo_->chamber(ChId)->toGlobal(LocalPoint(pos_f, 0., z));
  int thisec = MuonPathSLId.sector();
  if (thisec == 13)
	thisec = 4;
  if (thisec == 14)
	thisec = 10;
  double phi = jm_x_cmssw_global.phi() - PHI_CONV * (thisec - 1);
  double psi = atan(slope_f);
  double phiB = hasPosRF(MuonPathSLId.wheel(), MuonPathSLId.sector()) ? psi - phi : -psi - phi;

  // get the lateralities (in reverse order) in order to fill the metaprimitive
  std::vector <int> lateralities = getLateralityCombination(latcomb);
  for (int lay = 0; lay < NUM_LAYERS; lay++) {
    if (valid[lay] == 0)
      lateralities[lay] = -1;
  }

  metaPrimitives.emplace_back(metaPrimitive({MuonPathSLId.rawId(),
                        double(bx_time),
                        pos_f,
                        slope_f,
                        phi,
                        phiB,
                        chi2_f,
                        quality,
                        wires[0],
                        t0s[0],
                        lateralities[0],
                        wires[1],
                        t0s[1],
                        lateralities[1],
                        wires[2],
                        t0s[2],
                        lateralities[2],
                        wires[3],
                        t0s[3],
                        lateralities[3],
                        -1,
                        -1,
                        -1,
                        -1,
                        -1,
                        -1,
                        -1,
                        -1,
                        -1,
                        -1,
                        -1,
                        -1,
                        -1}));
}


std::vector <int> MuonPathAnalyticAnalyzer::getLateralityCombination (int latcomb) {
  // returns the latcomb as a binary number represented in a vector of integers
  // careful, the output is in reverse order
  std::vector <int> binaryNum = {};
  while (latcomb > 1) {
    binaryNum.push_back(latcomb % 2);
    latcomb = latcomb / 2;
  }
  binaryNum.push_back(latcomb);
  while (binaryNum.size() < 4) binaryNum.push_back(0);
  return binaryNum;
}


void MuonPathAnalyticAnalyzer::fillLAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER(){
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,0,-1},{1,1,0,1}},
    {
      {1, { {-6170, {1,0,0,-1}, 56936}, {239, {0,1,0,-1}, 4380}, {37, {0,1,0,-1}, 3559}, {776, {2,3,0,-1}, 16384}, }},
      {2, { {-30885, {-1,3,0,-2}, 18979}, {-1583769, {1,0,0,-1}, 2920}, {-6133, {1,0,0,-1}, 2372}, {-771, {2,3,0,1}, 10923}, }},
      {3, { {-6170, {1,0,0,-1}, 56936}, {-1584008, {-1,1,0,0}, 8759}, {-6170, {-1,1,0,0}, 7117}, {-773, {-2,3,0,1}, 32768}, }},
      {8, { {-6170, {-1,0,0,1}, 56936}, {-1584008, {1,-1,0,0}, 8759}, {-6170, {1,-1,0,0}, 7117}, {775, {-2,3,0,1}, 32768}, }},
      {9, { {-30885, {1,-3,0,2}, 18979}, {-1583769, {-1,0,0,1}, 2920}, {-6133, {-1,0,0,1}, 2372}, {777, {2,3,0,1}, 10923}, }},
      {10, { {-6170, {-1,0,0,1}, 56936}, {239, {0,-1,0,1}, 4380}, {37, {0,-1,0,1}, 3559}, {-772, {2,3,0,-1}, 16384}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,0,1},{0,1,1,1}},
    {
      {2, { {-6170, {0,1,-1,0}, 56936}, {1584248, {0,0,1,-1}, 8759}, {6206, {0,0,1,-1}, 7117}, {1, {0,1,2,-1}, 32768}, }},
      {4, { {-6170, {0,-1,1,0}, 56936}, {3168495, {0,1,0,-1}, 4380}, {12413, {0,1,0,-1}, 3559}, {2, {0,1,2,1}, 16384}, }},
      {6, { {-6170, {0,2,-1,-1}, 56936}, {1584248, {0,-1,1,0}, 8759}, {6206, {0,-1,1,0}, 7117}, {1, {0,-1,2,1}, 32768}, }},
      {8, { {-6170, {0,-2,1,1}, 56936}, {1584248, {0,1,-1,0}, 8759}, {6206, {0,1,-1,0}, 7117}, {1, {0,-1,2,1}, 32768}, }},
      {10, { {-6170, {0,1,-1,0}, 56936}, {3168495, {0,-1,0,1}, 4380}, {12413, {0,-1,0,1}, 3559}, {2, {0,1,2,1}, 16384}, }},
      {12, { {-6170, {0,-1,1,0}, 56936}, {1584248, {0,0,-1,1}, 8759}, {6206, {0,0,-1,1}, 7117}, {1, {0,1,2,-1}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,-2,-3},{1,1,0,1}},
    {
      {1, { {-18546, {1,0,0,-1}, 56936}, {-3168017, {0,1,0,-1}, 4380}, {-12339, {0,1,0,-1}, 3559}, {2, {2,3,0,-1}, 16384}, }},
      {2, { {-55637, {-1,3,0,-2}, 18979}, {-4752025, {1,0,0,-1}, 2920}, {-18509, {1,0,0,-1}, 2372}, {3, {2,3,0,1}, 10923}, }},
      {3, { {-18546, {1,0,0,-1}, 56936}, {-1584008, {-1,1,0,0}, 8759}, {-6170, {-1,1,0,0}, 7117}, {1, {-2,3,0,1}, 32768}, }},
      {8, { {-18546, {-1,0,0,1}, 56936}, {-1584008, {1,-1,0,0}, 8759}, {-6170, {1,-1,0,0}, 7117}, {1, {-2,3,0,1}, 32768}, }},
      {9, { {-55637, {1,-3,0,2}, 18979}, {-4752025, {-1,0,0,1}, 2920}, {-18509, {-1,0,0,1}, 2372}, {3, {2,3,0,1}, 10923}, }},
      {10, { {-18546, {-1,0,0,1}, 56936}, {-3168017, {0,-1,0,1}, 4380}, {-12339, {0,-1,0,1}, 3559}, {2, {2,3,0,-1}, 16384}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,0,1},{0,1,1,1}},
    {
      {2, { {6206, {0,1,-1,0}, 56936}, {1584248, {0,0,1,-1}, 8759}, {6206, {0,0,1,-1}, 7117}, {775, {0,1,2,-1}, 32768}, }},
      {4, { {6206, {0,-1,1,0}, 56936}, {239, {0,1,0,-1}, 4380}, {37, {0,1,0,-1}, 3559}, {-772, {0,1,2,1}, 16384}, }},
      {6, { {18582, {0,2,-1,-1}, 56936}, {-1584008, {0,-1,1,0}, 8759}, {-6170, {0,-1,1,0}, 7117}, {-773, {0,-1,2,1}, 32768}, }},
      {8, { {18582, {0,-2,1,1}, 56936}, {-1584008, {0,1,-1,0}, 8759}, {-6170, {0,1,-1,0}, 7117}, {775, {0,-1,2,1}, 32768}, }},
      {10, { {6206, {0,1,-1,0}, 56936}, {239, {0,-1,0,1}, 4380}, {37, {0,-1,0,1}, 3559}, {776, {0,1,2,1}, 16384}, }},
      {12, { {6206, {0,-1,1,0}, 56936}, {1584248, {0,0,-1,1}, 8759}, {6206, {0,0,-1,1}, 7117}, {-773, {0,1,2,-1}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,2,1},{1,1,1,0}},
    {
      {1, { {18582, {1,1,-2,0}, 56936}, {1584248, {0,1,-1,0}, 8759}, {6206, {0,1,-1,0}, 7117}, {1, {1,2,-1,0}, 32768}, }},
      {2, { {18582, {0,1,-1,0}, 56936}, {3168495, {1,0,-1,0}, 4380}, {12413, {1,0,-1,0}, 3559}, {2, {1,2,1,0}, 16384}, }},
      {3, { {18582, {0,1,-1,0}, 56936}, {1584248, {-1,1,0,0}, 8759}, {6206, {-1,1,0,0}, 7117}, {1, {-1,2,1,0}, 32768}, }},
      {4, { {18582, {0,-1,1,0}, 56936}, {1584248, {1,-1,0,0}, 8759}, {6206, {1,-1,0,0}, 7117}, {1, {-1,2,1,0}, 32768}, }},
      {5, { {18582, {0,-1,1,0}, 56936}, {3168495, {-1,0,1,0}, 4380}, {12413, {-1,0,1,0}, 3559}, {2, {1,2,1,0}, 16384}, }},
      {6, { {18582, {-1,-1,2,0}, 56936}, {1584248, {0,-1,1,0}, 8759}, {6206, {0,-1,1,0}, 7117}, {1, {1,2,-1,0}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,0,-1},{1,0,1,1}},
    {
      {1, { {-6170, {1,0,0,-1}, 56936}, {-1584008, {0,0,1,-1}, 8759}, {-6170, {0,0,1,-1}, 7117}, {-773, {1,0,3,-2}, 32768}, }},
      {4, { {-6133, {-2,0,3,-1}, 18979}, {-1583769, {1,0,0,-1}, 2920}, {-6133, {1,0,0,-1}, 2372}, {777, {1,0,3,2}, 10923}, }},
      {5, { {-6170, {1,0,0,-1}, 56936}, {239, {-1,0,1,0}, 4380}, {37, {-1,0,1,0}, 3559}, {776, {-1,0,3,2}, 16384}, }},
      {8, { {-6170, {-1,0,0,1}, 56936}, {239, {1,0,-1,0}, 4380}, {37, {1,0,-1,0}, 3559}, {-772, {-1,0,3,2}, 16384}, }},
      {9, { {-6133, {2,0,-3,1}, 18979}, {-1583769, {-1,0,0,1}, 2920}, {-6133, {-1,0,0,1}, 2372}, {-771, {1,0,3,2}, 10923}, }},
      {12, { {-6170, {-1,0,0,1}, 56936}, {-1584008, {0,0,-1,1}, 8759}, {-6170, {0,0,-1,1}, 7117}, {775, {1,0,3,-2}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,-2,-1},{1,1,1,0}},
    {
      {1, { {-18546, {1,1,-2,0}, 56936}, {-1584008, {0,1,-1,0}, 8759}, {-6170, {0,1,-1,0}, 7117}, {1, {1,2,-1,0}, 32768}, }},
      {2, { {-18546, {0,1,-1,0}, 56936}, {-3168017, {1,0,-1,0}, 4380}, {-12339, {1,0,-1,0}, 3559}, {2, {1,2,1,0}, 16384}, }},
      {3, { {-18546, {0,1,-1,0}, 56936}, {-1584008, {-1,1,0,0}, 8759}, {-6170, {-1,1,0,0}, 7117}, {1, {-1,2,1,0}, 32768}, }},
      {4, { {-18546, {0,-1,1,0}, 56936}, {-1584008, {1,-1,0,0}, 8759}, {-6170, {1,-1,0,0}, 7117}, {1, {-1,2,1,0}, 32768}, }},
      {5, { {-18546, {0,-1,1,0}, 56936}, {-3168017, {-1,0,1,0}, 4380}, {-12339, {-1,0,1,0}, 3559}, {2, {1,2,1,0}, 16384}, }},
      {6, { {-18546, {-1,-1,2,0}, 56936}, {-1584008, {0,-1,1,0}, 8759}, {-6170, {0,-1,1,0}, 7117}, {1, {1,2,-1,0}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,-2,-3},{0,1,1,1}},
    {
      {2, { {-18546, {0,1,-1,0}, 56936}, {-1584008, {0,0,1,-1}, 8759}, {-6170, {0,0,1,-1}, 7117}, {1, {0,1,2,-1}, 32768}, }},
      {4, { {-18546, {0,-1,1,0}, 56936}, {-3168017, {0,1,0,-1}, 4380}, {-12339, {0,1,0,-1}, 3559}, {2, {0,1,2,1}, 16384}, }},
      {6, { {-18546, {0,2,-1,-1}, 56936}, {-1584008, {0,-1,1,0}, 8759}, {-6170, {0,-1,1,0}, 7117}, {1, {0,-1,2,1}, 32768}, }},
      {8, { {-18546, {0,-2,1,1}, 56936}, {-1584008, {0,1,-1,0}, 8759}, {-6170, {0,1,-1,0}, 7117}, {1, {0,-1,2,1}, 32768}, }},
      {10, { {-18546, {0,1,-1,0}, 56936}, {-3168017, {0,-1,0,1}, 4380}, {-12339, {0,-1,0,1}, 3559}, {2, {0,1,2,1}, 16384}, }},
      {12, { {-18546, {0,-1,1,0}, 56936}, {-1584008, {0,0,-1,1}, 8759}, {-6170, {0,0,-1,1}, 7117}, {1, {0,1,2,-1}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,-2,-1},{0,1,1,1}},
    {
      {2, { {-18546, {0,1,-1,0}, 56936}, {1584248, {0,0,1,-1}, 8759}, {6206, {0,0,1,-1}, 7117}, {775, {0,1,2,-1}, 32768}, }},
      {4, { {-18546, {0,-1,1,0}, 56936}, {239, {0,1,0,-1}, 4380}, {37, {0,1,0,-1}, 3559}, {-772, {0,1,2,1}, 16384}, }},
      {6, { {-6170, {0,2,-1,-1}, 56936}, {-1584008, {0,-1,1,0}, 8759}, {-6170, {0,-1,1,0}, 7117}, {-773, {0,-1,2,1}, 32768}, }},
      {8, { {-6170, {0,-2,1,1}, 56936}, {-1584008, {0,1,-1,0}, 8759}, {-6170, {0,1,-1,0}, 7117}, {775, {0,-1,2,1}, 32768}, }},
      {10, { {-18546, {0,1,-1,0}, 56936}, {239, {0,-1,0,1}, 4380}, {37, {0,-1,0,1}, 3559}, {776, {0,1,2,1}, 16384}, }},
      {12, { {-18546, {0,-1,1,0}, 56936}, {1584248, {0,0,-1,1}, 8759}, {6206, {0,0,-1,1}, 7117}, {-773, {0,1,2,-1}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,-2,-3},{1,1,1,0}},
    {
      {1, { {-18546, {1,1,-2,0}, 56936}, {-1584008, {0,1,-1,0}, 8759}, {-6170, {0,1,-1,0}, 7117}, {1, {1,2,-1,0}, 32768}, }},
      {2, { {-18546, {0,1,-1,0}, 56936}, {-3168017, {1,0,-1,0}, 4380}, {-12339, {1,0,-1,0}, 3559}, {2, {1,2,1,0}, 16384}, }},
      {3, { {-18546, {0,1,-1,0}, 56936}, {-1584008, {-1,1,0,0}, 8759}, {-6170, {-1,1,0,0}, 7117}, {1, {-1,2,1,0}, 32768}, }},
      {4, { {-18546, {0,-1,1,0}, 56936}, {-1584008, {1,-1,0,0}, 8759}, {-6170, {1,-1,0,0}, 7117}, {1, {-1,2,1,0}, 32768}, }},
      {5, { {-18546, {0,-1,1,0}, 56936}, {-3168017, {-1,0,1,0}, 4380}, {-12339, {-1,0,1,0}, 3559}, {2, {1,2,1,0}, 16384}, }},
      {6, { {-18546, {-1,-1,2,0}, 56936}, {-1584008, {0,-1,1,0}, 8759}, {-6170, {0,-1,1,0}, 7117}, {1, {1,2,-1,0}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,2,1},{0,1,1,1}},
    {
      {2, { {18582, {0,1,-1,0}, 56936}, {-1584008, {0,0,1,-1}, 8759}, {-6170, {0,0,1,-1}, 7117}, {-773, {0,1,2,-1}, 32768}, }},
      {4, { {18582, {0,-1,1,0}, 56936}, {239, {0,1,0,-1}, 4380}, {37, {0,1,0,-1}, 3559}, {776, {0,1,2,1}, 16384}, }},
      {6, { {6206, {0,2,-1,-1}, 56936}, {1584248, {0,-1,1,0}, 8759}, {6206, {0,-1,1,0}, 7117}, {775, {0,-1,2,1}, 32768}, }},
      {8, { {6206, {0,-2,1,1}, 56936}, {1584248, {0,1,-1,0}, 8759}, {6206, {0,1,-1,0}, 7117}, {-773, {0,-1,2,1}, 32768}, }},
      {10, { {18582, {0,1,-1,0}, 56936}, {239, {0,-1,0,1}, 4380}, {37, {0,-1,0,1}, 3559}, {-772, {0,1,2,1}, 16384}, }},
      {12, { {18582, {0,-1,1,0}, 56936}, {-1584008, {0,0,-1,1}, 8759}, {-6170, {0,0,-1,1}, 7117}, {775, {0,1,2,-1}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,-2,-1},{1,1,1,1}},
    {
      {4, { {-222510, {-6,-5,14,-3}, 4067}, {-6334836, {4,1,0,-5}, 626}, {-24494, {4,1,0,-5}, 508}, {-3087, {1,2,7,4}, 4681}, }},
      {6, { {-24715, {-1,1,1,-1}, 28468}, {-6335315, {3,-1,1,-3}, 876}, {-24568, {3,-1,1,-3}, 712}, {-772, {1,1,1,1}, 16384}, }},
      {7, { {-37018, {5,2,-1,-6}, 9489}, {-3168017, {-1,0,1,0}, 4380}, {-12339, {-1,0,1,0}, 3559}, {-2318, {-2,1,4,3}, 10923}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,2,3},{0,1,1,1}},
    {
      {2, { {18582, {0,1,-1,0}, 56936}, {1584248, {0,0,1,-1}, 8759}, {6206, {0,0,1,-1}, 7117}, {1, {0,1,2,-1}, 32768}, }},
      {4, { {18582, {0,-1,1,0}, 56936}, {3168495, {0,1,0,-1}, 4380}, {12413, {0,1,0,-1}, 3559}, {2, {0,1,2,1}, 16384}, }},
      {6, { {18582, {0,2,-1,-1}, 56936}, {1584248, {0,-1,1,0}, 8759}, {6206, {0,-1,1,0}, 7117}, {1, {0,-1,2,1}, 32768}, }},
      {8, { {18582, {0,-2,1,1}, 56936}, {1584248, {0,1,-1,0}, 8759}, {6206, {0,1,-1,0}, 7117}, {1, {0,-1,2,1}, 32768}, }},
      {10, { {18582, {0,1,-1,0}, 56936}, {3168495, {0,-1,0,1}, 4380}, {12413, {0,-1,0,1}, 3559}, {2, {0,1,2,1}, 16384}, }},
      {12, { {18582, {0,-1,1,0}, 56936}, {1584248, {0,0,-1,1}, 8759}, {6206, {0,0,-1,1}, 7117}, {1, {0,1,2,-1}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,0,-1},{1,1,1,1}},
    {
      {1, { {-37018, {6,1,-2,-5}, 9489}, {-3168017, {0,1,0,-1}, 4380}, {-12339, {0,1,0,-1}, 3559}, {-2318, {3,4,1,-2}, 10923}, }},
      {9, { {37, {1,-1,-1,1}, 28468}, {-6335315, {-3,1,-1,3}, 876}, {-24568, {-3,1,-1,3}, 712}, {-772, {1,1,1,1}, 16384}, }},
      {13, { {49762, {3,-14,5,6}, 4067}, {-6334836, {-5,0,1,4}, 626}, {-24494, {-5,0,1,4}, 508}, {-3087, {4,7,2,1}, 4681}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,-2,-1},{1,1,0,1}},
    {
      {1, { {-6170, {1,0,0,-1}, 56936}, {239, {0,1,0,-1}, 4380}, {37, {0,1,0,-1}, 3559}, {776, {2,3,0,-1}, 16384}, }},
      {2, { {-30885, {-1,3,0,-2}, 18979}, {-1583769, {1,0,0,-1}, 2920}, {-6133, {1,0,0,-1}, 2372}, {-771, {2,3,0,1}, 10923}, }},
      {3, { {-6170, {1,0,0,-1}, 56936}, {-1584008, {-1,1,0,0}, 8759}, {-6170, {-1,1,0,0}, 7117}, {-773, {-2,3,0,1}, 32768}, }},
      {8, { {-6170, {-1,0,0,1}, 56936}, {-1584008, {1,-1,0,0}, 8759}, {-6170, {1,-1,0,0}, 7117}, {775, {-2,3,0,1}, 32768}, }},
      {9, { {-30885, {1,-3,0,2}, 18979}, {-1583769, {-1,0,0,1}, 2920}, {-6133, {-1,0,0,1}, 2372}, {777, {2,3,0,1}, 10923}, }},
      {10, { {-6170, {-1,0,0,1}, 56936}, {239, {0,-1,0,1}, 4380}, {37, {0,-1,0,1}, 3559}, {-772, {2,3,0,-1}, 16384}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,2,3},{1,1,0,1}},
    {
      {1, { {18582, {1,0,0,-1}, 56936}, {3168495, {0,1,0,-1}, 4380}, {12413, {0,1,0,-1}, 3559}, {2, {2,3,0,-1}, 16384}, }},
      {2, { {55747, {-1,3,0,-2}, 18979}, {4752743, {1,0,0,-1}, 2920}, {18619, {1,0,0,-1}, 2372}, {3, {2,3,0,1}, 10923}, }},
      {3, { {18582, {1,0,0,-1}, 56936}, {1584248, {-1,1,0,0}, 8759}, {6206, {-1,1,0,0}, 7117}, {1, {-2,3,0,1}, 32768}, }},
      {8, { {18582, {-1,0,0,1}, 56936}, {1584248, {1,-1,0,0}, 8759}, {6206, {1,-1,0,0}, 7117}, {1, {-2,3,0,1}, 32768}, }},
      {9, { {55747, {1,-3,0,2}, 18979}, {4752743, {-1,0,0,1}, 2920}, {18619, {-1,0,0,1}, 2372}, {3, {2,3,0,1}, 10923}, }},
      {10, { {18582, {-1,0,0,1}, 56936}, {3168495, {0,-1,0,1}, 4380}, {12413, {0,-1,0,1}, 3559}, {2, {2,3,0,-1}, 16384}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,0,1},{1,1,1,0}},
    {
      {1, { {6206, {1,1,-2,0}, 56936}, {1584248, {0,1,-1,0}, 8759}, {6206, {0,1,-1,0}, 7117}, {775, {1,2,-1,0}, 32768}, }},
      {2, { {-6170, {0,1,-1,0}, 56936}, {239, {1,0,-1,0}, 4380}, {37, {1,0,-1,0}, 3559}, {-772, {1,2,1,0}, 16384}, }},
      {3, { {-6170, {0,1,-1,0}, 56936}, {-1584008, {-1,1,0,0}, 8759}, {-6170, {-1,1,0,0}, 7117}, {-773, {-1,2,1,0}, 32768}, }},
      {4, { {-6170, {0,-1,1,0}, 56936}, {-1584008, {1,-1,0,0}, 8759}, {-6170, {1,-1,0,0}, 7117}, {775, {-1,2,1,0}, 32768}, }},
      {5, { {-6170, {0,-1,1,0}, 56936}, {239, {-1,0,1,0}, 4380}, {37, {-1,0,1,0}, 3559}, {776, {1,2,1,0}, 16384}, }},
      {6, { {6206, {-1,-1,2,0}, 56936}, {1584248, {0,-1,1,0}, 8759}, {6206, {0,-1,1,0}, 7117}, {-773, {1,2,-1,0}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,0,-1},{0,1,1,1}},
    {
      {2, { {6206, {0,1,-1,0}, 56936}, {-1584008, {0,0,1,-1}, 8759}, {-6170, {0,0,1,-1}, 7117}, {1, {0,1,2,-1}, 32768}, }},
      {4, { {6206, {0,-1,1,0}, 56936}, {-3168017, {0,1,0,-1}, 4380}, {-12339, {0,1,0,-1}, 3559}, {2, {0,1,2,1}, 16384}, }},
      {6, { {6206, {0,2,-1,-1}, 56936}, {-1584008, {0,-1,1,0}, 8759}, {-6170, {0,-1,1,0}, 7117}, {1, {0,-1,2,1}, 32768}, }},
      {8, { {6206, {0,-2,1,1}, 56936}, {-1584008, {0,1,-1,0}, 8759}, {-6170, {0,1,-1,0}, 7117}, {1, {0,-1,2,1}, 32768}, }},
      {10, { {6206, {0,1,-1,0}, 56936}, {-3168017, {0,-1,0,1}, 4380}, {-12339, {0,-1,0,1}, 3559}, {2, {0,1,2,1}, 16384}, }},
      {12, { {6206, {0,-1,1,0}, 56936}, {-1584008, {0,0,-1,1}, 8759}, {-6170, {0,0,-1,1}, 7117}, {1, {0,1,2,-1}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,0,-1},{1,1,1,1}},
    {
      {2, { {-123502, {-3,14,-5,-6}, 4067}, {-6334836, {5,0,-1,-4}, 626}, {-24494, {5,0,-1,-4}, 508}, {-2314, {4,7,2,1}, 4681}, }},
      {10, { {-12339, {-1,1,-1,1}, 28468}, {479, {1,-1,-1,1}, 2190}, {74, {1,-1,-1,1}, 1779}, {-1543, {1,3,3,1}, 8192}, }},
      {3, { {-12339, {1,1,-1,-1}, 28468}, {-3168017, {-1,1,1,-1}, 4380}, {-12339, {-1,1,1,-1}, 3559}, {-1545, {-1,3,3,-1}, 16384}, }},
      {11, { {-49246, {6,5,-14,3}, 4067}, {-6334836, {-4,-1,0,5}, 626}, {-24494, {-4,-1,0,5}, 508}, {-2314, {1,2,7,4}, 4681}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,0,-1},{0,1,1,1}},
    {
      {2, { {-6170, {0,1,-1,0}, 56936}, {-1584008, {0,0,1,-1}, 8759}, {-6170, {0,0,1,-1}, 7117}, {-773, {0,1,2,-1}, 32768}, }},
      {4, { {-6170, {0,-1,1,0}, 56936}, {239, {0,1,0,-1}, 4380}, {37, {0,1,0,-1}, 3559}, {776, {0,1,2,1}, 16384}, }},
      {6, { {-18546, {0,2,-1,-1}, 56936}, {1584248, {0,-1,1,0}, 8759}, {6206, {0,-1,1,0}, 7117}, {775, {0,-1,2,1}, 32768}, }},
      {8, { {-18546, {0,-2,1,1}, 56936}, {1584248, {0,1,-1,0}, 8759}, {6206, {0,1,-1,0}, 7117}, {-773, {0,-1,2,1}, 32768}, }},
      {10, { {-6170, {0,1,-1,0}, 56936}, {239, {0,-1,0,1}, 4380}, {37, {0,-1,0,1}, 3559}, {-772, {0,1,2,1}, 16384}, }},
      {12, { {-6170, {0,-1,1,0}, 56936}, {-1584008, {0,0,-1,1}, 8759}, {-6170, {0,0,-1,1}, 7117}, {775, {0,1,2,-1}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,2,3},{1,1,1,1}},
    {
      {8, { {111495, {-5,-2,1,6}, 9489}, {3168495, {1,0,-1,0}, 4380}, {12413, {1,0,-1,0}, 3559}, {3, {-2,1,4,3}, 10923}, }},
      {12, { {37165, {-1,-1,1,1}, 28468}, {3168495, {1,-1,-1,1}, 4380}, {12413, {1,-1,-1,1}, 3559}, {2, {-1,3,3,-1}, 16384}, }},
      {14, { {111495, {-6,-1,2,5}, 9489}, {3168495, {0,-1,0,1}, 4380}, {12413, {0,-1,0,1}, 3559}, {3, {3,4,1,-2}, 10923}, }},
      {1, { {111495, {6,1,-2,-5}, 9489}, {3168495, {0,1,0,-1}, 4380}, {12413, {0,1,0,-1}, 3559}, {3, {3,4,1,-2}, 10923}, }},
      {3, { {37165, {1,1,-1,-1}, 28468}, {3168495, {-1,1,1,-1}, 4380}, {12413, {-1,1,1,-1}, 3559}, {2, {-1,3,3,-1}, 16384}, }},
      {7, { {111495, {5,2,-1,-6}, 9489}, {3168495, {-1,0,1,0}, 4380}, {12413, {-1,0,1,0}, 3559}, {3, {-2,1,4,3}, 10923}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,0,1},{1,0,1,1}},
    {
      {1, { {6206, {1,0,0,-1}, 56936}, {1584248, {0,0,1,-1}, 8759}, {6206, {0,0,1,-1}, 7117}, {775, {1,0,3,-2}, 32768}, }},
      {4, { {6243, {-2,0,3,-1}, 18979}, {1584487, {1,0,0,-1}, 2920}, {6243, {1,0,0,-1}, 2372}, {-771, {1,0,3,2}, 10923}, }},
      {5, { {6206, {1,0,0,-1}, 56936}, {239, {-1,0,1,0}, 4380}, {37, {-1,0,1,0}, 3559}, {-772, {-1,0,3,2}, 16384}, }},
      {8, { {6206, {-1,0,0,1}, 56936}, {239, {1,0,-1,0}, 4380}, {37, {1,0,-1,0}, 3559}, {776, {-1,0,3,2}, 16384}, }},
      {9, { {6243, {2,0,-3,1}, 18979}, {1584487, {-1,0,0,1}, 2920}, {6243, {-1,0,0,1}, 2372}, {777, {1,0,3,2}, 10923}, }},
      {12, { {6206, {-1,0,0,1}, 56936}, {1584248, {0,0,-1,1}, 8759}, {6206, {0,0,-1,1}, 7117}, {-773, {1,0,3,-2}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,2,1},{1,1,0,1}},
    {
      {1, { {6206, {1,0,0,-1}, 56936}, {239, {0,1,0,-1}, 4380}, {37, {0,1,0,-1}, 3559}, {-772, {2,3,0,-1}, 16384}, }},
      {2, { {30995, {-1,3,0,-2}, 18979}, {1584487, {1,0,0,-1}, 2920}, {6243, {1,0,0,-1}, 2372}, {777, {2,3,0,1}, 10923}, }},
      {3, { {6206, {1,0,0,-1}, 56936}, {1584248, {-1,1,0,0}, 8759}, {6206, {-1,1,0,0}, 7117}, {775, {-2,3,0,1}, 32768}, }},
      {8, { {6206, {-1,0,0,1}, 56936}, {1584248, {1,-1,0,0}, 8759}, {6206, {1,-1,0,0}, 7117}, {-773, {-2,3,0,1}, 32768}, }},
      {9, { {30995, {1,-3,0,2}, 18979}, {1584487, {-1,0,0,1}, 2920}, {6243, {-1,0,0,1}, 2372}, {-771, {2,3,0,1}, 10923}, }},
      {10, { {6206, {-1,0,0,1}, 56936}, {239, {0,-1,0,1}, 4380}, {37, {0,-1,0,1}, 3559}, {776, {2,3,0,-1}, 16384}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,0,-1},{1,1,1,0}},
    {
      {1, { {6206, {1,1,-2,0}, 56936}, {1584248, {0,1,-1,0}, 8759}, {6206, {0,1,-1,0}, 7117}, {775, {1,2,-1,0}, 32768}, }},
      {2, { {-6170, {0,1,-1,0}, 56936}, {239, {1,0,-1,0}, 4380}, {37, {1,0,-1,0}, 3559}, {-772, {1,2,1,0}, 16384}, }},
      {3, { {-6170, {0,1,-1,0}, 56936}, {-1584008, {-1,1,0,0}, 8759}, {-6170, {-1,1,0,0}, 7117}, {-773, {-1,2,1,0}, 32768}, }},
      {4, { {-6170, {0,-1,1,0}, 56936}, {-1584008, {1,-1,0,0}, 8759}, {-6170, {1,-1,0,0}, 7117}, {775, {-1,2,1,0}, 32768}, }},
      {5, { {-6170, {0,-1,1,0}, 56936}, {239, {-1,0,1,0}, 4380}, {37, {-1,0,1,0}, 3559}, {776, {1,2,1,0}, 16384}, }},
      {6, { {6206, {-1,-1,2,0}, 56936}, {1584248, {0,-1,1,0}, 8759}, {6206, {0,-1,1,0}, 7117}, {-773, {1,2,-1,0}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,2,1},{1,0,1,1}},
    {
      {1, { {6206, {1,0,0,-1}, 56936}, {-1584008, {0,0,1,-1}, 8759}, {-6170, {0,0,1,-1}, 7117}, {-1546, {1,0,3,-2}, 32768}, }},
      {4, { {43371, {-2,0,3,-1}, 18979}, {1584487, {1,0,0,-1}, 2920}, {6243, {1,0,0,-1}, 2372}, {1550, {1,0,3,2}, 10923}, }},
      {5, { {6206, {1,0,0,-1}, 56936}, {3168495, {-1,0,1,0}, 4380}, {12413, {-1,0,1,0}, 3559}, {1549, {-1,0,3,2}, 16384}, }},
      {8, { {6206, {-1,0,0,1}, 56936}, {3168495, {1,0,-1,0}, 4380}, {12413, {1,0,-1,0}, 3559}, {-1545, {-1,0,3,2}, 16384}, }},
      {9, { {43371, {2,0,-3,1}, 18979}, {1584487, {-1,0,0,1}, 2920}, {6243, {-1,0,0,1}, 2372}, {-1544, {1,0,3,2}, 10923}, }},
      {12, { {6206, {-1,0,0,1}, 56936}, {-1584008, {0,0,-1,1}, 8759}, {-6170, {0,0,-1,1}, 7117}, {1548, {1,0,3,-2}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,0,-1},{1,1,1,0}},
    {
      {1, { {-6170, {1,1,-2,0}, 56936}, {-1584008, {0,1,-1,0}, 8759}, {-6170, {0,1,-1,0}, 7117}, {-773, {1,2,-1,0}, 32768}, }},
      {2, { {6206, {0,1,-1,0}, 56936}, {239, {1,0,-1,0}, 4380}, {37, {1,0,-1,0}, 3559}, {776, {1,2,1,0}, 16384}, }},
      {3, { {6206, {0,1,-1,0}, 56936}, {1584248, {-1,1,0,0}, 8759}, {6206, {-1,1,0,0}, 7117}, {775, {-1,2,1,0}, 32768}, }},
      {4, { {6206, {0,-1,1,0}, 56936}, {1584248, {1,-1,0,0}, 8759}, {6206, {1,-1,0,0}, 7117}, {-773, {-1,2,1,0}, 32768}, }},
      {5, { {6206, {0,-1,1,0}, 56936}, {239, {-1,0,1,0}, 4380}, {37, {-1,0,1,0}, 3559}, {-772, {1,2,1,0}, 16384}, }},
      {6, { {-6170, {-1,-1,2,0}, 56936}, {-1584008, {0,-1,1,0}, 8759}, {-6170, {0,-1,1,0}, 7117}, {775, {1,2,-1,0}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,2,3},{1,1,1,0}},
    {
      {1, { {18582, {1,1,-2,0}, 56936}, {1584248, {0,1,-1,0}, 8759}, {6206, {0,1,-1,0}, 7117}, {1, {1,2,-1,0}, 32768}, }},
      {2, { {18582, {0,1,-1,0}, 56936}, {3168495, {1,0,-1,0}, 4380}, {12413, {1,0,-1,0}, 3559}, {2, {1,2,1,0}, 16384}, }},
      {3, { {18582, {0,1,-1,0}, 56936}, {1584248, {-1,1,0,0}, 8759}, {6206, {-1,1,0,0}, 7117}, {1, {-1,2,1,0}, 32768}, }},
      {4, { {18582, {0,-1,1,0}, 56936}, {1584248, {1,-1,0,0}, 8759}, {6206, {1,-1,0,0}, 7117}, {1, {-1,2,1,0}, 32768}, }},
      {5, { {18582, {0,-1,1,0}, 56936}, {3168495, {-1,0,1,0}, 4380}, {12413, {-1,0,1,0}, 3559}, {2, {1,2,1,0}, 16384}, }},
      {6, { {18582, {-1,-1,2,0}, 56936}, {1584248, {0,-1,1,0}, 8759}, {6206, {0,-1,1,0}, 7117}, {1, {1,2,-1,0}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,0,1},{1,1,1,0}},
    {
      {1, { {-6170, {1,1,-2,0}, 56936}, {-1584008, {0,1,-1,0}, 8759}, {-6170, {0,1,-1,0}, 7117}, {-773, {1,2,-1,0}, 32768}, }},
      {2, { {6206, {0,1,-1,0}, 56936}, {239, {1,0,-1,0}, 4380}, {37, {1,0,-1,0}, 3559}, {776, {1,2,1,0}, 16384}, }},
      {3, { {6206, {0,1,-1,0}, 56936}, {1584248, {-1,1,0,0}, 8759}, {6206, {-1,1,0,0}, 7117}, {775, {-1,2,1,0}, 32768}, }},
      {4, { {6206, {0,-1,1,0}, 56936}, {1584248, {1,-1,0,0}, 8759}, {6206, {1,-1,0,0}, 7117}, {-773, {-1,2,1,0}, 32768}, }},
      {5, { {6206, {0,-1,1,0}, 56936}, {239, {-1,0,1,0}, 4380}, {37, {-1,0,1,0}, 3559}, {-772, {1,2,1,0}, 16384}, }},
      {6, { {-6170, {-1,-1,2,0}, 56936}, {-1584008, {0,-1,1,0}, 8759}, {-6170, {0,-1,1,0}, 7117}, {775, {1,2,-1,0}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,0,1},{1,1,1,1}},
    {
      {2, { {-49246, {-3,14,-5,-6}, 4067}, {6338188, {5,0,-1,-4}, 626}, {25010, {5,0,-1,-4}, 508}, {-3087, {4,7,2,1}, 4681}, }},
      {6, { {37, {-1,1,1,-1}, 28468}, {6337709, {3,-1,1,-3}, 876}, {24936, {3,-1,1,-3}, 712}, {-772, {1,1,1,1}, 16384}, }},
      {14, { {37239, {-6,-1,2,5}, 9489}, {3168495, {0,-1,0,1}, 4380}, {12413, {0,-1,0,1}, 3559}, {-2318, {3,4,1,-2}, 10923}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,-2,-3},{1,0,1,1}},
    {
      {1, { {-18546, {1,0,0,-1}, 56936}, {-1584008, {0,0,1,-1}, 8759}, {-6170, {0,0,1,-1}, 7117}, {1, {1,0,3,-2}, 32768}, }},
      {4, { {-55637, {-2,0,3,-1}, 18979}, {-4752025, {1,0,0,-1}, 2920}, {-18509, {1,0,0,-1}, 2372}, {3, {1,0,3,2}, 10923}, }},
      {5, { {-18546, {1,0,0,-1}, 56936}, {-3168017, {-1,0,1,0}, 4380}, {-12339, {-1,0,1,0}, 3559}, {2, {-1,0,3,2}, 16384}, }},
      {8, { {-18546, {-1,0,0,1}, 56936}, {-3168017, {1,0,-1,0}, 4380}, {-12339, {1,0,-1,0}, 3559}, {2, {-1,0,3,2}, 16384}, }},
      {9, { {-55637, {2,0,-3,1}, 18979}, {-4752025, {-1,0,0,1}, 2920}, {-18509, {-1,0,0,1}, 2372}, {3, {1,0,3,2}, 10923}, }},
      {12, { {-18546, {-1,0,0,1}, 56936}, {-1584008, {0,0,-1,1}, 8759}, {-6170, {0,0,-1,1}, 7117}, {1, {1,0,3,-2}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,0,1},{1,1,1,1}},
    {
      {4, { {49762, {-6,-5,14,-3}, 4067}, {6338188, {4,1,0,-5}, 626}, {25010, {4,1,0,-5}, 508}, {-2314, {1,2,7,4}, 4681}, }},
      {12, { {12413, {-1,-1,1,1}, 28468}, {3168495, {1,-1,-1,1}, 4380}, {12413, {1,-1,-1,1}, 3559}, {-1545, {-1,3,3,-1}, 16384}, }},
      {5, { {12413, {1,-1,1,-1}, 28468}, {479, {-1,1,1,-1}, 2190}, {74, {-1,1,1,-1}, 1779}, {-1543, {1,3,3,1}, 8192}, }},
      {13, { {124018, {3,-14,5,6}, 4067}, {6338188, {-5,0,1,4}, 626}, {25010, {-5,0,1,4}, 508}, {-2314, {4,7,2,1}, 4681}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,0,1},{1,0,1,1}},
    {
      {1, { {6206, {1,0,0,-1}, 56936}, {1584248, {0,0,1,-1}, 8759}, {6206, {0,0,1,-1}, 7117}, {775, {1,0,3,-2}, 32768}, }},
      {4, { {6243, {-2,0,3,-1}, 18979}, {1584487, {1,0,0,-1}, 2920}, {6243, {1,0,0,-1}, 2372}, {-771, {1,0,3,2}, 10923}, }},
      {5, { {6206, {1,0,0,-1}, 56936}, {239, {-1,0,1,0}, 4380}, {37, {-1,0,1,0}, 3559}, {-772, {-1,0,3,2}, 16384}, }},
      {8, { {6206, {-1,0,0,1}, 56936}, {239, {1,0,-1,0}, 4380}, {37, {1,0,-1,0}, 3559}, {776, {-1,0,3,2}, 16384}, }},
      {9, { {6243, {2,0,-3,1}, 18979}, {1584487, {-1,0,0,1}, 2920}, {6243, {-1,0,0,1}, 2372}, {777, {1,0,3,2}, 10923}, }},
      {12, { {6206, {-1,0,0,1}, 56936}, {1584248, {0,0,-1,1}, 8759}, {6206, {0,0,-1,1}, 7117}, {-773, {1,0,3,-2}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,0,-1},{1,0,1,1}},
    {
      {1, { {-6170, {1,0,0,-1}, 56936}, {-1584008, {0,0,1,-1}, 8759}, {-6170, {0,0,1,-1}, 7117}, {-773, {1,0,3,-2}, 32768}, }},
      {4, { {-6133, {-2,0,3,-1}, 18979}, {-1583769, {1,0,0,-1}, 2920}, {-6133, {1,0,0,-1}, 2372}, {777, {1,0,3,2}, 10923}, }},
      {5, { {-6170, {1,0,0,-1}, 56936}, {239, {-1,0,1,0}, 4380}, {37, {-1,0,1,0}, 3559}, {776, {-1,0,3,2}, 16384}, }},
      {8, { {-6170, {-1,0,0,1}, 56936}, {239, {1,0,-1,0}, 4380}, {37, {1,0,-1,0}, 3559}, {-772, {-1,0,3,2}, 16384}, }},
      {9, { {-6133, {2,0,-3,1}, 18979}, {-1583769, {-1,0,0,1}, 2920}, {-6133, {-1,0,0,1}, 2372}, {-771, {1,0,3,2}, 10923}, }},
      {12, { {-6170, {-1,0,0,1}, 56936}, {-1584008, {0,0,-1,1}, 8759}, {-6170, {0,0,-1,1}, 7117}, {775, {1,0,3,-2}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,2,1},{1,1,1,1}},
    {
      {8, { {37239, {-5,-2,1,6}, 9489}, {3168495, {1,0,-1,0}, 4380}, {12413, {1,0,-1,0}, 3559}, {-2318, {-2,1,4,3}, 10923}, }},
      {9, { {24789, {1,-1,-1,1}, 28468}, {6337709, {-3,1,-1,3}, 876}, {24936, {-3,1,-1,3}, 712}, {-772, {1,1,1,1}, 16384}, }},
      {11, { {223026, {6,5,-14,3}, 4067}, {6338188, {-4,-1,0,5}, 626}, {25010, {-4,-1,0,5}, 508}, {-3087, {1,2,7,4}, 4681}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,-2,-1},{1,0,1,1}},
    {
      {1, { {-6170, {1,0,0,-1}, 56936}, {1584248, {0,0,1,-1}, 8759}, {6206, {0,0,1,-1}, 7117}, {1548, {1,0,3,-2}, 32768}, }},
      {4, { {-43261, {-2,0,3,-1}, 18979}, {-1583769, {1,0,0,-1}, 2920}, {-6133, {1,0,0,-1}, 2372}, {-1544, {1,0,3,2}, 10923}, }},
      {5, { {-6170, {1,0,0,-1}, 56936}, {-3168017, {-1,0,1,0}, 4380}, {-12339, {-1,0,1,0}, 3559}, {-1545, {-1,0,3,2}, 16384}, }},
      {8, { {-6170, {-1,0,0,1}, 56936}, {-3168017, {1,0,-1,0}, 4380}, {-12339, {1,0,-1,0}, 3559}, {1549, {-1,0,3,2}, 16384}, }},
      {9, { {-43261, {2,0,-3,1}, 18979}, {-1583769, {-1,0,0,1}, 2920}, {-6133, {-1,0,0,1}, 2372}, {1550, {1,0,3,2}, 10923}, }},
      {12, { {-6170, {-1,0,0,1}, 56936}, {1584248, {0,0,-1,1}, 8759}, {6206, {0,0,-1,1}, 7117}, {-1546, {1,0,3,-2}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,0,-1},{1,1,0,1}},
    {
      {1, { {-6170, {1,0,0,-1}, 56936}, {-3168017, {0,1,0,-1}, 4380}, {-12339, {0,1,0,-1}, 3559}, {-1545, {2,3,0,-1}, 16384}, }},
      {2, { {6243, {-1,3,0,-2}, 18979}, {-1583769, {1,0,0,-1}, 2920}, {-6133, {1,0,0,-1}, 2372}, {1550, {2,3,0,1}, 10923}, }},
      {3, { {-6170, {1,0,0,-1}, 56936}, {1584248, {-1,1,0,0}, 8759}, {6206, {-1,1,0,0}, 7117}, {1548, {-2,3,0,1}, 32768}, }},
      {8, { {-6170, {-1,0,0,1}, 56936}, {1584248, {1,-1,0,0}, 8759}, {6206, {1,-1,0,0}, 7117}, {-1546, {-2,3,0,1}, 32768}, }},
      {9, { {6243, {1,-3,0,2}, 18979}, {-1583769, {-1,0,0,1}, 2920}, {-6133, {-1,0,0,1}, 2372}, {-1544, {2,3,0,1}, 10923}, }},
      {10, { {-6170, {-1,0,0,1}, 56936}, {-3168017, {0,-1,0,1}, 4380}, {-12339, {0,-1,0,1}, 3559}, {1549, {2,3,0,-1}, 16384}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,-2,-3},{1,1,1,1}},
    {
      {8, { {-111274, {-5,-2,1,6}, 9489}, {-3168017, {1,0,-1,0}, 4380}, {-12339, {1,0,-1,0}, 3559}, {3, {-2,1,4,3}, 10923}, }},
      {12, { {-37091, {-1,-1,1,1}, 28468}, {-3168017, {1,-1,-1,1}, 4380}, {-12339, {1,-1,-1,1}, 3559}, {2, {-1,3,3,-1}, 16384}, }},
      {14, { {-111274, {-6,-1,2,5}, 9489}, {-3168017, {0,-1,0,1}, 4380}, {-12339, {0,-1,0,1}, 3559}, {3, {3,4,1,-2}, 10923}, }},
      {1, { {-111274, {6,1,-2,-5}, 9489}, {-3168017, {0,1,0,-1}, 4380}, {-12339, {0,1,0,-1}, 3559}, {3, {3,4,1,-2}, 10923}, }},
      {3, { {-37091, {1,1,-1,-1}, 28468}, {-3168017, {-1,1,1,-1}, 4380}, {-12339, {-1,1,1,-1}, 3559}, {2, {-1,3,3,-1}, 16384}, }},
      {7, { {-111274, {5,2,-1,-6}, 9489}, {-3168017, {-1,0,1,0}, 4380}, {-12339, {-1,0,1,0}, 3559}, {3, {-2,1,4,3}, 10923}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,0,1},{1,1,0,1}},
    {
      {1, { {6206, {1,0,0,-1}, 56936}, {239, {0,1,0,-1}, 4380}, {37, {0,1,0,-1}, 3559}, {-772, {2,3,0,-1}, 16384}, }},
      {2, { {30995, {-1,3,0,-2}, 18979}, {1584487, {1,0,0,-1}, 2920}, {6243, {1,0,0,-1}, 2372}, {777, {2,3,0,1}, 10923}, }},
      {3, { {6206, {1,0,0,-1}, 56936}, {1584248, {-1,1,0,0}, 8759}, {6206, {-1,1,0,0}, 7117}, {775, {-2,3,0,1}, 32768}, }},
      {8, { {6206, {-1,0,0,1}, 56936}, {1584248, {1,-1,0,0}, 8759}, {6206, {1,-1,0,0}, 7117}, {-773, {-2,3,0,1}, 32768}, }},
      {9, { {30995, {1,-3,0,2}, 18979}, {1584487, {-1,0,0,1}, 2920}, {6243, {-1,0,0,1}, 2372}, {-771, {2,3,0,1}, 10923}, }},
      {10, { {6206, {-1,0,0,1}, 56936}, {239, {0,-1,0,1}, 4380}, {37, {0,-1,0,1}, 3559}, {776, {2,3,0,-1}, 16384}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,1,2,3},{1,0,1,1}},
    {
      {1, { {18582, {1,0,0,-1}, 56936}, {1584248, {0,0,1,-1}, 8759}, {6206, {0,0,1,-1}, 7117}, {1, {1,0,3,-2}, 32768}, }},
      {4, { {55747, {-2,0,3,-1}, 18979}, {4752743, {1,0,0,-1}, 2920}, {18619, {1,0,0,-1}, 2372}, {3, {1,0,3,2}, 10923}, }},
      {5, { {18582, {1,0,0,-1}, 56936}, {3168495, {-1,0,1,0}, 4380}, {12413, {-1,0,1,0}, 3559}, {2, {-1,0,3,2}, 16384}, }},
      {8, { {18582, {-1,0,0,1}, 56936}, {3168495, {1,0,-1,0}, 4380}, {12413, {1,0,-1,0}, 3559}, {2, {-1,0,3,2}, 16384}, }},
      {9, { {55747, {2,0,-3,1}, 18979}, {4752743, {-1,0,0,1}, 2920}, {18619, {-1,0,0,1}, 2372}, {3, {1,0,3,2}, 10923}, }},
      {12, { {18582, {-1,0,0,1}, 56936}, {1584248, {0,0,-1,1}, 8759}, {6206, {0,0,-1,1}, 7117}, {1, {1,0,3,-2}, 32768}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER.push_back({
    {{0,-1,0,1},{1,1,0,1}},
    {
      {1, { {6206, {1,0,0,-1}, 56936}, {3168495, {0,1,0,-1}, 4380}, {12413, {0,1,0,-1}, 3559}, {1549, {2,3,0,-1}, 16384}, }},
      {2, { {-6133, {-1,3,0,-2}, 18979}, {1584487, {1,0,0,-1}, 2920}, {6243, {1,0,0,-1}, 2372}, {-1544, {2,3,0,1}, 10923}, }},
      {3, { {6206, {1,0,0,-1}, 56936}, {-1584008, {-1,1,0,0}, 8759}, {-6170, {-1,1,0,0}, 7117}, {-1546, {-2,3,0,1}, 32768}, }},
      {8, { {6206, {-1,0,0,1}, 56936}, {-1584008, {1,-1,0,0}, 8759}, {-6170, {1,-1,0,0}, 7117}, {1548, {-2,3,0,1}, 32768}, }},
      {9, { {-6133, {1,-3,0,2}, 18979}, {1584487, {-1,0,0,1}, 2920}, {6243, {-1,0,0,1}, 2372}, {1550, {2,3,0,1}, 10923}, }},
      {10, { {6206, {-1,0,0,1}, 56936}, {3168495, {0,-1,0,1}, 4380}, {12413, {0,-1,0,1}, 3559}, {-1545, {2,3,0,-1}, 16384}, }},
    }});
}

void MuonPathAnalyticAnalyzer::fillLAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL(){
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,0,-1},{1,1,0,1}},
    {
      {1, { {-386.75, {1,0,0,-1}, 36.8333333333}, {0.0, {0,1,0,-1}, 478.833333333}, {0.0, {0,1,0,-1}, 73.6666666667}, {773.5, {2,3,0,-1}, 4.0}, }},
      {2, { {-1933.75, {-1,3,0,-2}, 110.5}, {-386.75, {1,0,0,-1}, 718.25}, {-386.75, {1,0,0,-1}, 110.5}, {-773.5, {2,3,0,1}, 6.0}, }},
      {3, { {-386.75, {1,0,0,-1}, 36.8333333333}, {-386.75, {-1,1,0,0}, 239.416666667}, {-386.75, {-1,1,0,0}, 36.8333333333}, {-773.5, {-2,3,0,1}, 2.0}, }},
      {8, { {-386.75, {-1,0,0,1}, 36.8333333333}, {-386.75, {1,-1,0,0}, 239.416666667}, {-386.75, {1,-1,0,0}, 36.8333333333}, {773.5, {-2,3,0,1}, 2.0}, }},
      {9, { {-1933.75, {1,-3,0,2}, 110.5}, {-386.75, {-1,0,0,1}, 718.25}, {-386.75, {-1,0,0,1}, 110.5}, {773.5, {2,3,0,1}, 6.0}, }},
      {10, { {-386.75, {-1,0,0,1}, 36.8333333333}, {0.0, {0,-1,0,1}, 478.833333333}, {0.0, {0,-1,0,1}, 73.6666666667}, {-773.5, {2,3,0,-1}, 4.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,0,1},{0,1,1,1}},
    {
      {2, { {-386.75, {0,1,-1,0}, 36.8333333333}, {386.75, {0,0,1,-1}, 239.416666667}, {386.75, {0,0,1,-1}, 36.8333333333}, {0.0, {0,1,2,-1}, 2.0}, }},
      {4, { {-386.75, {0,-1,1,0}, 36.8333333333}, {773.5, {0,1,0,-1}, 478.833333333}, {773.5, {0,1,0,-1}, 73.6666666667}, {0.0, {0,1,2,1}, 4.0}, }},
      {6, { {-386.75, {0,2,-1,-1}, 36.8333333333}, {386.75, {0,-1,1,0}, 239.416666667}, {386.75, {0,-1,1,0}, 36.8333333333}, {0.0, {0,-1,2,1}, 2.0}, }},
      {8, { {-386.75, {0,-2,1,1}, 36.8333333333}, {386.75, {0,1,-1,0}, 239.416666667}, {386.75, {0,1,-1,0}, 36.8333333333}, {0.0, {0,-1,2,1}, 2.0}, }},
      {10, { {-386.75, {0,1,-1,0}, 36.8333333333}, {773.5, {0,-1,0,1}, 478.833333333}, {773.5, {0,-1,0,1}, 73.6666666667}, {0.0, {0,1,2,1}, 4.0}, }},
      {12, { {-386.75, {0,-1,1,0}, 36.8333333333}, {386.75, {0,0,-1,1}, 239.416666667}, {386.75, {0,0,-1,1}, 36.8333333333}, {0.0, {0,1,2,-1}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,-2,-3},{1,1,0,1}},
    {
      {1, { {-1160.25, {1,0,0,-1}, 36.8333333333}, {-773.5, {0,1,0,-1}, 478.833333333}, {-773.5, {0,1,0,-1}, 73.6666666667}, {0.0, {2,3,0,-1}, 4.0}, }},
      {2, { {-3480.75, {-1,3,0,-2}, 110.5}, {-1160.25, {1,0,0,-1}, 718.25}, {-1160.25, {1,0,0,-1}, 110.5}, {-0.0, {2,3,0,1}, 6.0}, }},
      {3, { {-1160.25, {1,0,0,-1}, 36.8333333333}, {-386.75, {-1,1,0,0}, 239.416666667}, {-386.75, {-1,1,0,0}, 36.8333333333}, {0.0, {-2,3,0,1}, 2.0}, }},
      {8, { {-1160.25, {-1,0,0,1}, 36.8333333333}, {-386.75, {1,-1,0,0}, 239.416666667}, {-386.75, {1,-1,0,0}, 36.8333333333}, {-0.0, {-2,3,0,1}, 2.0}, }},
      {9, { {-3480.75, {1,-3,0,2}, 110.5}, {-1160.25, {-1,0,0,1}, 718.25}, {-1160.25, {-1,0,0,1}, 110.5}, {0.0, {2,3,0,1}, 6.0}, }},
      {10, { {-1160.25, {-1,0,0,1}, 36.8333333333}, {-773.5, {0,-1,0,1}, 478.833333333}, {-773.5, {0,-1,0,1}, 73.6666666667}, {0.0, {2,3,0,-1}, 4.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,0,1},{0,1,1,1}},
    {
      {2, { {386.75, {0,1,-1,0}, 36.8333333333}, {386.75, {0,0,1,-1}, 239.416666667}, {386.75, {0,0,1,-1}, 36.8333333333}, {773.5, {0,1,2,-1}, 2.0}, }},
      {4, { {386.75, {0,-1,1,0}, 36.8333333333}, {0.0, {0,1,0,-1}, 478.833333333}, {0.0, {0,1,0,-1}, 73.6666666667}, {-773.5, {0,1,2,1}, 4.0}, }},
      {6, { {1160.25, {0,2,-1,-1}, 36.8333333333}, {-386.75, {0,-1,1,0}, 239.416666667}, {-386.75, {0,-1,1,0}, 36.8333333333}, {-773.5, {0,-1,2,1}, 2.0}, }},
      {8, { {1160.25, {0,-2,1,1}, 36.8333333333}, {-386.75, {0,1,-1,0}, 239.416666667}, {-386.75, {0,1,-1,0}, 36.8333333333}, {773.5, {0,-1,2,1}, 2.0}, }},
      {10, { {386.75, {0,1,-1,0}, 36.8333333333}, {0.0, {0,-1,0,1}, 478.833333333}, {0.0, {0,-1,0,1}, 73.6666666667}, {773.5, {0,1,2,1}, 4.0}, }},
      {12, { {386.75, {0,-1,1,0}, 36.8333333333}, {386.75, {0,0,-1,1}, 239.416666667}, {386.75, {0,0,-1,1}, 36.8333333333}, {-773.5, {0,1,2,-1}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,2,1},{1,1,1,0}},
    {
      {1, { {1160.25, {1,1,-2,0}, 36.8333333333}, {386.75, {0,1,-1,0}, 239.416666667}, {386.75, {0,1,-1,0}, 36.8333333333}, {0.0, {1,2,-1,0}, 2.0}, }},
      {2, { {1160.25, {0,1,-1,0}, 36.8333333333}, {773.5, {1,0,-1,0}, 478.833333333}, {773.5, {1,0,-1,0}, 73.6666666667}, {-0.0, {1,2,1,0}, 4.0}, }},
      {3, { {1160.25, {0,1,-1,0}, 36.8333333333}, {386.75, {-1,1,0,0}, 239.416666667}, {386.75, {-1,1,0,0}, 36.8333333333}, {-0.0, {-1,2,1,0}, 2.0}, }},
      {4, { {1160.25, {0,-1,1,0}, 36.8333333333}, {386.75, {1,-1,0,0}, 239.416666667}, {386.75, {1,-1,0,0}, 36.8333333333}, {0.0, {-1,2,1,0}, 2.0}, }},
      {5, { {1160.25, {0,-1,1,0}, 36.8333333333}, {773.5, {-1,0,1,0}, 478.833333333}, {773.5, {-1,0,1,0}, 73.6666666667}, {0.0, {1,2,1,0}, 4.0}, }},
      {6, { {1160.25, {-1,-1,2,0}, 36.8333333333}, {386.75, {0,-1,1,0}, 239.416666667}, {386.75, {0,-1,1,0}, 36.8333333333}, {-0.0, {1,2,-1,0}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,0,-1},{1,0,1,1}},
    {
      {1, { {-386.75, {1,0,0,-1}, 36.8333333333}, {-386.75, {0,0,1,-1}, 239.416666667}, {-386.75, {0,0,1,-1}, 36.8333333333}, {-773.5, {1,0,3,-2}, 2.0}, }},
      {4, { {-386.75, {-2,0,3,-1}, 110.5}, {-386.75, {1,0,0,-1}, 718.25}, {-386.75, {1,0,0,-1}, 110.5}, {773.5, {1,0,3,2}, 6.0}, }},
      {5, { {-386.75, {1,0,0,-1}, 36.8333333333}, {0.0, {-1,0,1,0}, 478.833333333}, {0.0, {-1,0,1,0}, 73.6666666667}, {773.5, {-1,0,3,2}, 4.0}, }},
      {8, { {-386.75, {-1,0,0,1}, 36.8333333333}, {0.0, {1,0,-1,0}, 478.833333333}, {0.0, {1,0,-1,0}, 73.6666666667}, {-773.5, {-1,0,3,2}, 4.0}, }},
      {9, { {-386.75, {2,0,-3,1}, 110.5}, {-386.75, {-1,0,0,1}, 718.25}, {-386.75, {-1,0,0,1}, 110.5}, {-773.5, {1,0,3,2}, 6.0}, }},
      {12, { {-386.75, {-1,0,0,1}, 36.8333333333}, {-386.75, {0,0,-1,1}, 239.416666667}, {-386.75, {0,0,-1,1}, 36.8333333333}, {773.5, {1,0,3,-2}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,-2,-1},{1,1,1,0}},
    {
      {1, { {-1160.25, {1,1,-2,0}, 36.8333333333}, {-386.75, {0,1,-1,0}, 239.416666667}, {-386.75, {0,1,-1,0}, 36.8333333333}, {-0.0, {1,2,-1,0}, 2.0}, }},
      {2, { {-1160.25, {0,1,-1,0}, 36.8333333333}, {-773.5, {1,0,-1,0}, 478.833333333}, {-773.5, {1,0,-1,0}, 73.6666666667}, {0.0, {1,2,1,0}, 4.0}, }},
      {3, { {-1160.25, {0,1,-1,0}, 36.8333333333}, {-386.75, {-1,1,0,0}, 239.416666667}, {-386.75, {-1,1,0,0}, 36.8333333333}, {0.0, {-1,2,1,0}, 2.0}, }},
      {4, { {-1160.25, {0,-1,1,0}, 36.8333333333}, {-386.75, {1,-1,0,0}, 239.416666667}, {-386.75, {1,-1,0,0}, 36.8333333333}, {-0.0, {-1,2,1,0}, 2.0}, }},
      {5, { {-1160.25, {0,-1,1,0}, 36.8333333333}, {-773.5, {-1,0,1,0}, 478.833333333}, {-773.5, {-1,0,1,0}, 73.6666666667}, {-0.0, {1,2,1,0}, 4.0}, }},
      {6, { {-1160.25, {-1,-1,2,0}, 36.8333333333}, {-386.75, {0,-1,1,0}, 239.416666667}, {-386.75, {0,-1,1,0}, 36.8333333333}, {0.0, {1,2,-1,0}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,-2,-3},{0,1,1,1}},
    {
      {2, { {-1160.25, {0,1,-1,0}, 36.8333333333}, {-386.75, {0,0,1,-1}, 239.416666667}, {-386.75, {0,0,1,-1}, 36.8333333333}, {-0.0, {0,1,2,-1}, 2.0}, }},
      {4, { {-1160.25, {0,-1,1,0}, 36.8333333333}, {-773.5, {0,1,0,-1}, 478.833333333}, {-773.5, {0,1,0,-1}, 73.6666666667}, {0.0, {0,1,2,1}, 4.0}, }},
      {6, { {-1160.25, {0,2,-1,-1}, 36.8333333333}, {-386.75, {0,-1,1,0}, 239.416666667}, {-386.75, {0,-1,1,0}, 36.8333333333}, {0.0, {0,-1,2,1}, 2.0}, }},
      {8, { {-1160.25, {0,-2,1,1}, 36.8333333333}, {-386.75, {0,1,-1,0}, 239.416666667}, {-386.75, {0,1,-1,0}, 36.8333333333}, {-0.0, {0,-1,2,1}, 2.0}, }},
      {10, { {-1160.25, {0,1,-1,0}, 36.8333333333}, {-773.5, {0,-1,0,1}, 478.833333333}, {-773.5, {0,-1,0,1}, 73.6666666667}, {-0.0, {0,1,2,1}, 4.0}, }},
      {12, { {-1160.25, {0,-1,1,0}, 36.8333333333}, {-386.75, {0,0,-1,1}, 239.416666667}, {-386.75, {0,0,-1,1}, 36.8333333333}, {0.0, {0,1,2,-1}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,-2,-1},{0,1,1,1}},
    {
      {2, { {-1160.25, {0,1,-1,0}, 36.8333333333}, {386.75, {0,0,1,-1}, 239.416666667}, {386.75, {0,0,1,-1}, 36.8333333333}, {773.5, {0,1,2,-1}, 2.0}, }},
      {4, { {-1160.25, {0,-1,1,0}, 36.8333333333}, {0.0, {0,1,0,-1}, 478.833333333}, {0.0, {0,1,0,-1}, 73.6666666667}, {-773.5, {0,1,2,1}, 4.0}, }},
      {6, { {-386.75, {0,2,-1,-1}, 36.8333333333}, {-386.75, {0,-1,1,0}, 239.416666667}, {-386.75, {0,-1,1,0}, 36.8333333333}, {-773.5, {0,-1,2,1}, 2.0}, }},
      {8, { {-386.75, {0,-2,1,1}, 36.8333333333}, {-386.75, {0,1,-1,0}, 239.416666667}, {-386.75, {0,1,-1,0}, 36.8333333333}, {773.5, {0,-1,2,1}, 2.0}, }},
      {10, { {-1160.25, {0,1,-1,0}, 36.8333333333}, {0.0, {0,-1,0,1}, 478.833333333}, {0.0, {0,-1,0,1}, 73.6666666667}, {773.5, {0,1,2,1}, 4.0}, }},
      {12, { {-1160.25, {0,-1,1,0}, 36.8333333333}, {386.75, {0,0,-1,1}, 239.416666667}, {386.75, {0,0,-1,1}, 36.8333333333}, {-773.5, {0,1,2,-1}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,-2,-3},{1,1,1,0}},
    {
      {1, { {-1160.25, {1,1,-2,0}, 36.8333333333}, {-386.75, {0,1,-1,0}, 239.416666667}, {-386.75, {0,1,-1,0}, 36.8333333333}, {-0.0, {1,2,-1,0}, 2.0}, }},
      {2, { {-1160.25, {0,1,-1,0}, 36.8333333333}, {-773.5, {1,0,-1,0}, 478.833333333}, {-773.5, {1,0,-1,0}, 73.6666666667}, {0.0, {1,2,1,0}, 4.0}, }},
      {3, { {-1160.25, {0,1,-1,0}, 36.8333333333}, {-386.75, {-1,1,0,0}, 239.416666667}, {-386.75, {-1,1,0,0}, 36.8333333333}, {0.0, {-1,2,1,0}, 2.0}, }},
      {4, { {-1160.25, {0,-1,1,0}, 36.8333333333}, {-386.75, {1,-1,0,0}, 239.416666667}, {-386.75, {1,-1,0,0}, 36.8333333333}, {-0.0, {-1,2,1,0}, 2.0}, }},
      {5, { {-1160.25, {0,-1,1,0}, 36.8333333333}, {-773.5, {-1,0,1,0}, 478.833333333}, {-773.5, {-1,0,1,0}, 73.6666666667}, {-0.0, {1,2,1,0}, 4.0}, }},
      {6, { {-1160.25, {-1,-1,2,0}, 36.8333333333}, {-386.75, {0,-1,1,0}, 239.416666667}, {-386.75, {0,-1,1,0}, 36.8333333333}, {0.0, {1,2,-1,0}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,2,1},{0,1,1,1}},
    {
      {2, { {1160.25, {0,1,-1,0}, 36.8333333333}, {-386.75, {0,0,1,-1}, 239.416666667}, {-386.75, {0,0,1,-1}, 36.8333333333}, {-773.5, {0,1,2,-1}, 2.0}, }},
      {4, { {1160.25, {0,-1,1,0}, 36.8333333333}, {0.0, {0,1,0,-1}, 478.833333333}, {0.0, {0,1,0,-1}, 73.6666666667}, {773.5, {0,1,2,1}, 4.0}, }},
      {6, { {386.75, {0,2,-1,-1}, 36.8333333333}, {386.75, {0,-1,1,0}, 239.416666667}, {386.75, {0,-1,1,0}, 36.8333333333}, {773.5, {0,-1,2,1}, 2.0}, }},
      {8, { {386.75, {0,-2,1,1}, 36.8333333333}, {386.75, {0,1,-1,0}, 239.416666667}, {386.75, {0,1,-1,0}, 36.8333333333}, {-773.5, {0,-1,2,1}, 2.0}, }},
      {10, { {1160.25, {0,1,-1,0}, 36.8333333333}, {0.0, {0,-1,0,1}, 478.833333333}, {0.0, {0,-1,0,1}, 73.6666666667}, {-773.5, {0,1,2,1}, 4.0}, }},
      {12, { {1160.25, {0,-1,1,0}, 36.8333333333}, {-386.75, {0,0,-1,1}, 239.416666667}, {-386.75, {0,0,-1,1}, 36.8333333333}, {773.5, {0,1,2,-1}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,-2,-1},{1,1,1,1}},
    {
      {4, { {-13923.0, {-6,-5,14,-3}, 515.666666667}, {-1547.0, {4,1,0,-5}, 3351.83333333}, {-1547.0, {4,1,0,-5}, 515.666666667}, {-3094.0, {1,2,7,4}, 14.0}, }},
      {6, { {-1547.0, {-1,1,1,-1}, 73.6666666667}, {-1547.0, {3,-1,1,-3}, 2394.16666667}, {-1547.0, {3,-1,1,-3}, 368.333333333}, {-773.5, {1,1,1,1}, 4.0}, }},
      {7, { {-2320.5, {5,2,-1,-6}, 221.0}, {-773.5, {-1,0,1,0}, 478.833333333}, {-773.5, {-1,0,1,0}, 73.6666666667}, {-2320.5, {-2,1,4,3}, 6.0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,2,3},{0,1,1,1}},
    {
      {2, { {1160.25, {0,1,-1,0}, 36.8333333333}, {386.75, {0,0,1,-1}, 239.416666667}, {386.75, {0,0,1,-1}, 36.8333333333}, {0.0, {0,1,2,-1}, 2.0}, }},
      {4, { {1160.25, {0,-1,1,0}, 36.8333333333}, {773.5, {0,1,0,-1}, 478.833333333}, {773.5, {0,1,0,-1}, 73.6666666667}, {-0.0, {0,1,2,1}, 4.0}, }},
      {6, { {1160.25, {0,2,-1,-1}, 36.8333333333}, {386.75, {0,-1,1,0}, 239.416666667}, {386.75, {0,-1,1,0}, 36.8333333333}, {-0.0, {0,-1,2,1}, 2.0}, }},
      {8, { {1160.25, {0,-2,1,1}, 36.8333333333}, {386.75, {0,1,-1,0}, 239.416666667}, {386.75, {0,1,-1,0}, 36.8333333333}, {0.0, {0,-1,2,1}, 2.0}, }},
      {10, { {1160.25, {0,1,-1,0}, 36.8333333333}, {773.5, {0,-1,0,1}, 478.833333333}, {773.5, {0,-1,0,1}, 73.6666666667}, {0.0, {0,1,2,1}, 4.0}, }},
      {12, { {1160.25, {0,-1,1,0}, 36.8333333333}, {386.75, {0,0,-1,1}, 239.416666667}, {386.75, {0,0,-1,1}, 36.8333333333}, {-0.0, {0,1,2,-1}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,0,-1},{1,1,1,1}},
    {
      {1, { {-2320.5, {6,1,-2,-5}, 221.0}, {-773.5, {0,1,0,-1}, 478.833333333}, {-773.5, {0,1,0,-1}, 73.6666666667}, {-2320.5, {3,4,1,-2}, 6.0}, }},
      {9, { {0.0, {1,-1,-1,1}, 73.6666666667}, {-1547.0, {-3,1,-1,3}, 2394.16666667}, {-1547.0, {-3,1,-1,3}, 368.333333333}, {-773.5, {1,1,1,1}, 4.0}, }},
      {13, { {3094.0, {3,-14,5,6}, 515.666666667}, {-1547.0, {-5,0,1,4}, 3351.83333333}, {-1547.0, {-5,0,1,4}, 515.666666667}, {-3094.0, {4,7,2,1}, 14.0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,-2,-1},{1,1,0,1}},
    {
      {1, { {-386.75, {1,0,0,-1}, 36.8333333333}, {0.0, {0,1,0,-1}, 478.833333333}, {0.0, {0,1,0,-1}, 73.6666666667}, {773.5, {2,3,0,-1}, 4.0}, }},
      {2, { {-1933.75, {-1,3,0,-2}, 110.5}, {-386.75, {1,0,0,-1}, 718.25}, {-386.75, {1,0,0,-1}, 110.5}, {-773.5, {2,3,0,1}, 6.0}, }},
      {3, { {-386.75, {1,0,0,-1}, 36.8333333333}, {-386.75, {-1,1,0,0}, 239.416666667}, {-386.75, {-1,1,0,0}, 36.8333333333}, {-773.5, {-2,3,0,1}, 2.0}, }},
      {8, { {-386.75, {-1,0,0,1}, 36.8333333333}, {-386.75, {1,-1,0,0}, 239.416666667}, {-386.75, {1,-1,0,0}, 36.8333333333}, {773.5, {-2,3,0,1}, 2.0}, }},
      {9, { {-1933.75, {1,-3,0,2}, 110.5}, {-386.75, {-1,0,0,1}, 718.25}, {-386.75, {-1,0,0,1}, 110.5}, {773.5, {2,3,0,1}, 6.0}, }},
      {10, { {-386.75, {-1,0,0,1}, 36.8333333333}, {0.0, {0,-1,0,1}, 478.833333333}, {0.0, {0,-1,0,1}, 73.6666666667}, {-773.5, {2,3,0,-1}, 4.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,2,3},{1,1,0,1}},
    {
      {1, { {1160.25, {1,0,0,-1}, 36.8333333333}, {773.5, {0,1,0,-1}, 478.833333333}, {773.5, {0,1,0,-1}, 73.6666666667}, {0.0, {2,3,0,-1}, 4.0}, }},
      {2, { {3480.75, {-1,3,0,-2}, 110.5}, {1160.25, {1,0,0,-1}, 718.25}, {1160.25, {1,0,0,-1}, 110.5}, {0.0, {2,3,0,1}, 6.0}, }},
      {3, { {1160.25, {1,0,0,-1}, 36.8333333333}, {386.75, {-1,1,0,0}, 239.416666667}, {386.75, {-1,1,0,0}, 36.8333333333}, {-0.0, {-2,3,0,1}, 2.0}, }},
      {8, { {1160.25, {-1,0,0,1}, 36.8333333333}, {386.75, {1,-1,0,0}, 239.416666667}, {386.75, {1,-1,0,0}, 36.8333333333}, {0.0, {-2,3,0,1}, 2.0}, }},
      {9, { {3480.75, {1,-3,0,2}, 110.5}, {1160.25, {-1,0,0,1}, 718.25}, {1160.25, {-1,0,0,1}, 110.5}, {-0.0, {2,3,0,1}, 6.0}, }},
      {10, { {1160.25, {-1,0,0,1}, 36.8333333333}, {773.5, {0,-1,0,1}, 478.833333333}, {773.5, {0,-1,0,1}, 73.6666666667}, {0.0, {2,3,0,-1}, 4.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,0,1},{1,1,1,0}},
    {
      {1, { {386.75, {1,1,-2,0}, 36.8333333333}, {386.75, {0,1,-1,0}, 239.416666667}, {386.75, {0,1,-1,0}, 36.8333333333}, {773.5, {1,2,-1,0}, 2.0}, }},
      {2, { {-386.75, {0,1,-1,0}, 36.8333333333}, {0.0, {1,0,-1,0}, 478.833333333}, {0.0, {1,0,-1,0}, 73.6666666667}, {-773.5, {1,2,1,0}, 4.0}, }},
      {3, { {-386.75, {0,1,-1,0}, 36.8333333333}, {-386.75, {-1,1,0,0}, 239.416666667}, {-386.75, {-1,1,0,0}, 36.8333333333}, {-773.5, {-1,2,1,0}, 2.0}, }},
      {4, { {-386.75, {0,-1,1,0}, 36.8333333333}, {-386.75, {1,-1,0,0}, 239.416666667}, {-386.75, {1,-1,0,0}, 36.8333333333}, {773.5, {-1,2,1,0}, 2.0}, }},
      {5, { {-386.75, {0,-1,1,0}, 36.8333333333}, {0.0, {-1,0,1,0}, 478.833333333}, {0.0, {-1,0,1,0}, 73.6666666667}, {773.5, {1,2,1,0}, 4.0}, }},
      {6, { {386.75, {-1,-1,2,0}, 36.8333333333}, {386.75, {0,-1,1,0}, 239.416666667}, {386.75, {0,-1,1,0}, 36.8333333333}, {-773.5, {1,2,-1,0}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,0,-1},{0,1,1,1}},
    {
      {2, { {386.75, {0,1,-1,0}, 36.8333333333}, {-386.75, {0,0,1,-1}, 239.416666667}, {-386.75, {0,0,1,-1}, 36.8333333333}, {0.0, {0,1,2,-1}, 2.0}, }},
      {4, { {386.75, {0,-1,1,0}, 36.8333333333}, {-773.5, {0,1,0,-1}, 478.833333333}, {-773.5, {0,1,0,-1}, 73.6666666667}, {0.0, {0,1,2,1}, 4.0}, }},
      {6, { {386.75, {0,2,-1,-1}, 36.8333333333}, {-386.75, {0,-1,1,0}, 239.416666667}, {-386.75, {0,-1,1,0}, 36.8333333333}, {0.0, {0,-1,2,1}, 2.0}, }},
      {8, { {386.75, {0,-2,1,1}, 36.8333333333}, {-386.75, {0,1,-1,0}, 239.416666667}, {-386.75, {0,1,-1,0}, 36.8333333333}, {0.0, {0,-1,2,1}, 2.0}, }},
      {10, { {386.75, {0,1,-1,0}, 36.8333333333}, {-773.5, {0,-1,0,1}, 478.833333333}, {-773.5, {0,-1,0,1}, 73.6666666667}, {0.0, {0,1,2,1}, 4.0}, }},
      {12, { {386.75, {0,-1,1,0}, 36.8333333333}, {-386.75, {0,0,-1,1}, 239.416666667}, {-386.75, {0,0,-1,1}, 36.8333333333}, {0.0, {0,1,2,-1}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,0,-1},{1,1,1,1}},
    {
      {2, { {-7735.0, {-3,14,-5,-6}, 515.666666667}, {-1547.0, {5,0,-1,-4}, 3351.83333333}, {-1547.0, {5,0,-1,-4}, 515.666666667}, {-2320.5, {4,7,2,1}, 14.0}, }},
      {10, { {-773.5, {-1,1,-1,1}, 73.6666666667}, {0.0, {1,-1,-1,1}, 957.666666667}, {0.0, {1,-1,-1,1}, 147.333333333}, {-1547.0, {1,3,3,1}, 8.0}, }},
      {3, { {-773.5, {1,1,-1,-1}, 73.6666666667}, {-773.5, {-1,1,1,-1}, 478.833333333}, {-773.5, {-1,1,1,-1}, 73.6666666667}, {-1547.0, {-1,3,3,-1}, 4.0}, }},
      {11, { {-3094.0, {6,5,-14,3}, 515.666666667}, {-1547.0, {-4,-1,0,5}, 3351.83333333}, {-1547.0, {-4,-1,0,5}, 515.666666667}, {-2320.5, {1,2,7,4}, 14.0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,0,-1},{0,1,1,1}},
    {
      {2, { {-386.75, {0,1,-1,0}, 36.8333333333}, {-386.75, {0,0,1,-1}, 239.416666667}, {-386.75, {0,0,1,-1}, 36.8333333333}, {-773.5, {0,1,2,-1}, 2.0}, }},
      {4, { {-386.75, {0,-1,1,0}, 36.8333333333}, {0.0, {0,1,0,-1}, 478.833333333}, {0.0, {0,1,0,-1}, 73.6666666667}, {773.5, {0,1,2,1}, 4.0}, }},
      {6, { {-1160.25, {0,2,-1,-1}, 36.8333333333}, {386.75, {0,-1,1,0}, 239.416666667}, {386.75, {0,-1,1,0}, 36.8333333333}, {773.5, {0,-1,2,1}, 2.0}, }},
      {8, { {-1160.25, {0,-2,1,1}, 36.8333333333}, {386.75, {0,1,-1,0}, 239.416666667}, {386.75, {0,1,-1,0}, 36.8333333333}, {-773.5, {0,-1,2,1}, 2.0}, }},
      {10, { {-386.75, {0,1,-1,0}, 36.8333333333}, {0.0, {0,-1,0,1}, 478.833333333}, {0.0, {0,-1,0,1}, 73.6666666667}, {-773.5, {0,1,2,1}, 4.0}, }},
      {12, { {-386.75, {0,-1,1,0}, 36.8333333333}, {-386.75, {0,0,-1,1}, 239.416666667}, {-386.75, {0,0,-1,1}, 36.8333333333}, {773.5, {0,1,2,-1}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,2,3},{1,1,1,1}},
    {
      {8, { {6961.5, {-5,-2,1,6}, 221.0}, {773.5, {1,0,-1,0}, 478.833333333}, {773.5, {1,0,-1,0}, 73.6666666667}, {0.0, {-2,1,4,3}, 6.0}, }},
      {12, { {2320.5, {-1,-1,1,1}, 73.6666666667}, {773.5, {1,-1,-1,1}, 478.833333333}, {773.5, {1,-1,-1,1}, 73.6666666667}, {-0.0, {-1,3,3,-1}, 4.0}, }},
      {14, { {6961.5, {-6,-1,2,5}, 221.0}, {773.5, {0,-1,0,1}, 478.833333333}, {773.5, {0,-1,0,1}, 73.6666666667}, {0.0, {3,4,1,-2}, 6.0}, }},
      {1, { {6961.5, {6,1,-2,-5}, 221.0}, {773.5, {0,1,0,-1}, 478.833333333}, {773.5, {0,1,0,-1}, 73.6666666667}, {-0.0, {3,4,1,-2}, 6.0}, }},
      {3, { {2320.5, {1,1,-1,-1}, 73.6666666667}, {773.5, {-1,1,1,-1}, 478.833333333}, {773.5, {-1,1,1,-1}, 73.6666666667}, {0.0, {-1,3,3,-1}, 4.0}, }},
      {7, { {6961.5, {5,2,-1,-6}, 221.0}, {773.5, {-1,0,1,0}, 478.833333333}, {773.5, {-1,0,1,0}, 73.6666666667}, {0.0, {-2,1,4,3}, 6.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,0,1},{1,0,1,1}},
    {
      {1, { {386.75, {1,0,0,-1}, 36.8333333333}, {386.75, {0,0,1,-1}, 239.416666667}, {386.75, {0,0,1,-1}, 36.8333333333}, {773.5, {1,0,3,-2}, 2.0}, }},
      {4, { {386.75, {-2,0,3,-1}, 110.5}, {386.75, {1,0,0,-1}, 718.25}, {386.75, {1,0,0,-1}, 110.5}, {-773.5, {1,0,3,2}, 6.0}, }},
      {5, { {386.75, {1,0,0,-1}, 36.8333333333}, {0.0, {-1,0,1,0}, 478.833333333}, {0.0, {-1,0,1,0}, 73.6666666667}, {-773.5, {-1,0,3,2}, 4.0}, }},
      {8, { {386.75, {-1,0,0,1}, 36.8333333333}, {0.0, {1,0,-1,0}, 478.833333333}, {0.0, {1,0,-1,0}, 73.6666666667}, {773.5, {-1,0,3,2}, 4.0}, }},
      {9, { {386.75, {2,0,-3,1}, 110.5}, {386.75, {-1,0,0,1}, 718.25}, {386.75, {-1,0,0,1}, 110.5}, {773.5, {1,0,3,2}, 6.0}, }},
      {12, { {386.75, {-1,0,0,1}, 36.8333333333}, {386.75, {0,0,-1,1}, 239.416666667}, {386.75, {0,0,-1,1}, 36.8333333333}, {-773.5, {1,0,3,-2}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,2,1},{1,1,0,1}},
    {
      {1, { {386.75, {1,0,0,-1}, 36.8333333333}, {0.0, {0,1,0,-1}, 478.833333333}, {0.0, {0,1,0,-1}, 73.6666666667}, {-773.5, {2,3,0,-1}, 4.0}, }},
      {2, { {1933.75, {-1,3,0,-2}, 110.5}, {386.75, {1,0,0,-1}, 718.25}, {386.75, {1,0,0,-1}, 110.5}, {773.5, {2,3,0,1}, 6.0}, }},
      {3, { {386.75, {1,0,0,-1}, 36.8333333333}, {386.75, {-1,1,0,0}, 239.416666667}, {386.75, {-1,1,0,0}, 36.8333333333}, {773.5, {-2,3,0,1}, 2.0}, }},
      {8, { {386.75, {-1,0,0,1}, 36.8333333333}, {386.75, {1,-1,0,0}, 239.416666667}, {386.75, {1,-1,0,0}, 36.8333333333}, {-773.5, {-2,3,0,1}, 2.0}, }},
      {9, { {1933.75, {1,-3,0,2}, 110.5}, {386.75, {-1,0,0,1}, 718.25}, {386.75, {-1,0,0,1}, 110.5}, {-773.5, {2,3,0,1}, 6.0}, }},
      {10, { {386.75, {-1,0,0,1}, 36.8333333333}, {0.0, {0,-1,0,1}, 478.833333333}, {0.0, {0,-1,0,1}, 73.6666666667}, {773.5, {2,3,0,-1}, 4.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,0,-1},{1,1,1,0}},
    {
      {1, { {386.75, {1,1,-2,0}, 36.8333333333}, {386.75, {0,1,-1,0}, 239.416666667}, {386.75, {0,1,-1,0}, 36.8333333333}, {773.5, {1,2,-1,0}, 2.0}, }},
      {2, { {-386.75, {0,1,-1,0}, 36.8333333333}, {0.0, {1,0,-1,0}, 478.833333333}, {0.0, {1,0,-1,0}, 73.6666666667}, {-773.5, {1,2,1,0}, 4.0}, }},
      {3, { {-386.75, {0,1,-1,0}, 36.8333333333}, {-386.75, {-1,1,0,0}, 239.416666667}, {-386.75, {-1,1,0,0}, 36.8333333333}, {-773.5, {-1,2,1,0}, 2.0}, }},
      {4, { {-386.75, {0,-1,1,0}, 36.8333333333}, {-386.75, {1,-1,0,0}, 239.416666667}, {-386.75, {1,-1,0,0}, 36.8333333333}, {773.5, {-1,2,1,0}, 2.0}, }},
      {5, { {-386.75, {0,-1,1,0}, 36.8333333333}, {0.0, {-1,0,1,0}, 478.833333333}, {0.0, {-1,0,1,0}, 73.6666666667}, {773.5, {1,2,1,0}, 4.0}, }},
      {6, { {386.75, {-1,-1,2,0}, 36.8333333333}, {386.75, {0,-1,1,0}, 239.416666667}, {386.75, {0,-1,1,0}, 36.8333333333}, {-773.5, {1,2,-1,0}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,2,1},{1,0,1,1}},
    {
      {1, { {386.75, {1,0,0,-1}, 36.8333333333}, {-386.75, {0,0,1,-1}, 239.416666667}, {-386.75, {0,0,1,-1}, 36.8333333333}, {-1547.0, {1,0,3,-2}, 2.0}, }},
      {4, { {2707.25, {-2,0,3,-1}, 110.5}, {386.75, {1,0,0,-1}, 718.25}, {386.75, {1,0,0,-1}, 110.5}, {1547.0, {1,0,3,2}, 6.0}, }},
      {5, { {386.75, {1,0,0,-1}, 36.8333333333}, {773.5, {-1,0,1,0}, 478.833333333}, {773.5, {-1,0,1,0}, 73.6666666667}, {1547.0, {-1,0,3,2}, 4.0}, }},
      {8, { {386.75, {-1,0,0,1}, 36.8333333333}, {773.5, {1,0,-1,0}, 478.833333333}, {773.5, {1,0,-1,0}, 73.6666666667}, {-1547.0, {-1,0,3,2}, 4.0}, }},
      {9, { {2707.25, {2,0,-3,1}, 110.5}, {386.75, {-1,0,0,1}, 718.25}, {386.75, {-1,0,0,1}, 110.5}, {-1547.0, {1,0,3,2}, 6.0}, }},
      {12, { {386.75, {-1,0,0,1}, 36.8333333333}, {-386.75, {0,0,-1,1}, 239.416666667}, {-386.75, {0,0,-1,1}, 36.8333333333}, {1547.0, {1,0,3,-2}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,0,-1},{1,1,1,0}},
    {
      {1, { {-386.75, {1,1,-2,0}, 36.8333333333}, {-386.75, {0,1,-1,0}, 239.416666667}, {-386.75, {0,1,-1,0}, 36.8333333333}, {-773.5, {1,2,-1,0}, 2.0}, }},
      {2, { {386.75, {0,1,-1,0}, 36.8333333333}, {0.0, {1,0,-1,0}, 478.833333333}, {0.0, {1,0,-1,0}, 73.6666666667}, {773.5, {1,2,1,0}, 4.0}, }},
      {3, { {386.75, {0,1,-1,0}, 36.8333333333}, {386.75, {-1,1,0,0}, 239.416666667}, {386.75, {-1,1,0,0}, 36.8333333333}, {773.5, {-1,2,1,0}, 2.0}, }},
      {4, { {386.75, {0,-1,1,0}, 36.8333333333}, {386.75, {1,-1,0,0}, 239.416666667}, {386.75, {1,-1,0,0}, 36.8333333333}, {-773.5, {-1,2,1,0}, 2.0}, }},
      {5, { {386.75, {0,-1,1,0}, 36.8333333333}, {0.0, {-1,0,1,0}, 478.833333333}, {0.0, {-1,0,1,0}, 73.6666666667}, {-773.5, {1,2,1,0}, 4.0}, }},
      {6, { {-386.75, {-1,-1,2,0}, 36.8333333333}, {-386.75, {0,-1,1,0}, 239.416666667}, {-386.75, {0,-1,1,0}, 36.8333333333}, {773.5, {1,2,-1,0}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,2,3},{1,1,1,0}},
    {
      {1, { {1160.25, {1,1,-2,0}, 36.8333333333}, {386.75, {0,1,-1,0}, 239.416666667}, {386.75, {0,1,-1,0}, 36.8333333333}, {0.0, {1,2,-1,0}, 2.0}, }},
      {2, { {1160.25, {0,1,-1,0}, 36.8333333333}, {773.5, {1,0,-1,0}, 478.833333333}, {773.5, {1,0,-1,0}, 73.6666666667}, {-0.0, {1,2,1,0}, 4.0}, }},
      {3, { {1160.25, {0,1,-1,0}, 36.8333333333}, {386.75, {-1,1,0,0}, 239.416666667}, {386.75, {-1,1,0,0}, 36.8333333333}, {-0.0, {-1,2,1,0}, 2.0}, }},
      {4, { {1160.25, {0,-1,1,0}, 36.8333333333}, {386.75, {1,-1,0,0}, 239.416666667}, {386.75, {1,-1,0,0}, 36.8333333333}, {0.0, {-1,2,1,0}, 2.0}, }},
      {5, { {1160.25, {0,-1,1,0}, 36.8333333333}, {773.5, {-1,0,1,0}, 478.833333333}, {773.5, {-1,0,1,0}, 73.6666666667}, {0.0, {1,2,1,0}, 4.0}, }},
      {6, { {1160.25, {-1,-1,2,0}, 36.8333333333}, {386.75, {0,-1,1,0}, 239.416666667}, {386.75, {0,-1,1,0}, 36.8333333333}, {-0.0, {1,2,-1,0}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,0,1},{1,1,1,0}},
    {
      {1, { {-386.75, {1,1,-2,0}, 36.8333333333}, {-386.75, {0,1,-1,0}, 239.416666667}, {-386.75, {0,1,-1,0}, 36.8333333333}, {-773.5, {1,2,-1,0}, 2.0}, }},
      {2, { {386.75, {0,1,-1,0}, 36.8333333333}, {0.0, {1,0,-1,0}, 478.833333333}, {0.0, {1,0,-1,0}, 73.6666666667}, {773.5, {1,2,1,0}, 4.0}, }},
      {3, { {386.75, {0,1,-1,0}, 36.8333333333}, {386.75, {-1,1,0,0}, 239.416666667}, {386.75, {-1,1,0,0}, 36.8333333333}, {773.5, {-1,2,1,0}, 2.0}, }},
      {4, { {386.75, {0,-1,1,0}, 36.8333333333}, {386.75, {1,-1,0,0}, 239.416666667}, {386.75, {1,-1,0,0}, 36.8333333333}, {-773.5, {-1,2,1,0}, 2.0}, }},
      {5, { {386.75, {0,-1,1,0}, 36.8333333333}, {0.0, {-1,0,1,0}, 478.833333333}, {0.0, {-1,0,1,0}, 73.6666666667}, {-773.5, {1,2,1,0}, 4.0}, }},
      {6, { {-386.75, {-1,-1,2,0}, 36.8333333333}, {-386.75, {0,-1,1,0}, 239.416666667}, {-386.75, {0,-1,1,0}, 36.8333333333}, {773.5, {1,2,-1,0}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,0,1},{1,1,1,1}},
    {
      {2, { {-3094.0, {-3,14,-5,-6}, 515.666666667}, {1547.0, {5,0,-1,-4}, 3351.83333333}, {1547.0, {5,0,-1,-4}, 515.666666667}, {-3094.0, {4,7,2,1}, 14.0}, }},
      {6, { {0.0, {-1,1,1,-1}, 73.6666666667}, {1547.0, {3,-1,1,-3}, 2394.16666667}, {1547.0, {3,-1,1,-3}, 368.333333333}, {-773.5, {1,1,1,1}, 4.0}, }},
      {14, { {2320.5, {-6,-1,2,5}, 221.0}, {773.5, {0,-1,0,1}, 478.833333333}, {773.5, {0,-1,0,1}, 73.6666666667}, {-2320.5, {3,4,1,-2}, 6.0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,-2,-3},{1,0,1,1}},
    {
      {1, { {-1160.25, {1,0,0,-1}, 36.8333333333}, {-386.75, {0,0,1,-1}, 239.416666667}, {-386.75, {0,0,1,-1}, 36.8333333333}, {-0.0, {1,0,3,-2}, 2.0}, }},
      {4, { {-3480.75, {-2,0,3,-1}, 110.5}, {-1160.25, {1,0,0,-1}, 718.25}, {-1160.25, {1,0,0,-1}, 110.5}, {-0.0, {1,0,3,2}, 6.0}, }},
      {5, { {-1160.25, {1,0,0,-1}, 36.8333333333}, {-773.5, {-1,0,1,0}, 478.833333333}, {-773.5, {-1,0,1,0}, 73.6666666667}, {0.0, {-1,0,3,2}, 4.0}, }},
      {8, { {-1160.25, {-1,0,0,1}, 36.8333333333}, {-773.5, {1,0,-1,0}, 478.833333333}, {-773.5, {1,0,-1,0}, 73.6666666667}, {-0.0, {-1,0,3,2}, 4.0}, }},
      {9, { {-3480.75, {2,0,-3,1}, 110.5}, {-1160.25, {-1,0,0,1}, 718.25}, {-1160.25, {-1,0,0,1}, 110.5}, {0.0, {1,0,3,2}, 6.0}, }},
      {12, { {-1160.25, {-1,0,0,1}, 36.8333333333}, {-386.75, {0,0,-1,1}, 239.416666667}, {-386.75, {0,0,-1,1}, 36.8333333333}, {0.0, {1,0,3,-2}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,0,1},{1,1,1,1}},
    {
      {4, { {3094.0, {-6,-5,14,-3}, 515.666666667}, {1547.0, {4,1,0,-5}, 3351.83333333}, {1547.0, {4,1,0,-5}, 515.666666667}, {-2320.5, {1,2,7,4}, 14.0}, }},
      {12, { {773.5, {-1,-1,1,1}, 73.6666666667}, {773.5, {1,-1,-1,1}, 478.833333333}, {773.5, {1,-1,-1,1}, 73.6666666667}, {-1547.0, {-1,3,3,-1}, 4.0}, }},
      {5, { {773.5, {1,-1,1,-1}, 73.6666666667}, {0.0, {-1,1,1,-1}, 957.666666667}, {0.0, {-1,1,1,-1}, 147.333333333}, {-1547.0, {1,3,3,1}, 8.0}, }},
      {13, { {7735.0, {3,-14,5,6}, 515.666666667}, {1547.0, {-5,0,1,4}, 3351.83333333}, {1547.0, {-5,0,1,4}, 515.666666667}, {-2320.5, {4,7,2,1}, 14.0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,0,1},{1,0,1,1}},
    {
      {1, { {386.75, {1,0,0,-1}, 36.8333333333}, {386.75, {0,0,1,-1}, 239.416666667}, {386.75, {0,0,1,-1}, 36.8333333333}, {773.5, {1,0,3,-2}, 2.0}, }},
      {4, { {386.75, {-2,0,3,-1}, 110.5}, {386.75, {1,0,0,-1}, 718.25}, {386.75, {1,0,0,-1}, 110.5}, {-773.5, {1,0,3,2}, 6.0}, }},
      {5, { {386.75, {1,0,0,-1}, 36.8333333333}, {0.0, {-1,0,1,0}, 478.833333333}, {0.0, {-1,0,1,0}, 73.6666666667}, {-773.5, {-1,0,3,2}, 4.0}, }},
      {8, { {386.75, {-1,0,0,1}, 36.8333333333}, {0.0, {1,0,-1,0}, 478.833333333}, {0.0, {1,0,-1,0}, 73.6666666667}, {773.5, {-1,0,3,2}, 4.0}, }},
      {9, { {386.75, {2,0,-3,1}, 110.5}, {386.75, {-1,0,0,1}, 718.25}, {386.75, {-1,0,0,1}, 110.5}, {773.5, {1,0,3,2}, 6.0}, }},
      {12, { {386.75, {-1,0,0,1}, 36.8333333333}, {386.75, {0,0,-1,1}, 239.416666667}, {386.75, {0,0,-1,1}, 36.8333333333}, {-773.5, {1,0,3,-2}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,0,-1},{1,0,1,1}},
    {
      {1, { {-386.75, {1,0,0,-1}, 36.8333333333}, {-386.75, {0,0,1,-1}, 239.416666667}, {-386.75, {0,0,1,-1}, 36.8333333333}, {-773.5, {1,0,3,-2}, 2.0}, }},
      {4, { {-386.75, {-2,0,3,-1}, 110.5}, {-386.75, {1,0,0,-1}, 718.25}, {-386.75, {1,0,0,-1}, 110.5}, {773.5, {1,0,3,2}, 6.0}, }},
      {5, { {-386.75, {1,0,0,-1}, 36.8333333333}, {0.0, {-1,0,1,0}, 478.833333333}, {0.0, {-1,0,1,0}, 73.6666666667}, {773.5, {-1,0,3,2}, 4.0}, }},
      {8, { {-386.75, {-1,0,0,1}, 36.8333333333}, {0.0, {1,0,-1,0}, 478.833333333}, {0.0, {1,0,-1,0}, 73.6666666667}, {-773.5, {-1,0,3,2}, 4.0}, }},
      {9, { {-386.75, {2,0,-3,1}, 110.5}, {-386.75, {-1,0,0,1}, 718.25}, {-386.75, {-1,0,0,1}, 110.5}, {-773.5, {1,0,3,2}, 6.0}, }},
      {12, { {-386.75, {-1,0,0,1}, 36.8333333333}, {-386.75, {0,0,-1,1}, 239.416666667}, {-386.75, {0,0,-1,1}, 36.8333333333}, {773.5, {1,0,3,-2}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,2,1},{1,1,1,1}},
    {
      {8, { {2320.5, {-5,-2,1,6}, 221.0}, {773.5, {1,0,-1,0}, 478.833333333}, {773.5, {1,0,-1,0}, 73.6666666667}, {-2320.5, {-2,1,4,3}, 6.0}, }},
      {9, { {1547.0, {1,-1,-1,1}, 73.6666666667}, {1547.0, {-3,1,-1,3}, 2394.16666667}, {1547.0, {-3,1,-1,3}, 368.333333333}, {-773.5, {1,1,1,1}, 4.0}, }},
      {11, { {13923.0, {6,5,-14,3}, 515.666666667}, {1547.0, {-4,-1,0,5}, 3351.83333333}, {1547.0, {-4,-1,0,5}, 515.666666667}, {-3094.0, {1,2,7,4}, 14.0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
      {0, { {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, {0, {0,0,0,0}, 0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,-2,-1},{1,0,1,1}},
    {
      {1, { {-386.75, {1,0,0,-1}, 36.8333333333}, {386.75, {0,0,1,-1}, 239.416666667}, {386.75, {0,0,1,-1}, 36.8333333333}, {1547.0, {1,0,3,-2}, 2.0}, }},
      {4, { {-2707.25, {-2,0,3,-1}, 110.5}, {-386.75, {1,0,0,-1}, 718.25}, {-386.75, {1,0,0,-1}, 110.5}, {-1547.0, {1,0,3,2}, 6.0}, }},
      {5, { {-386.75, {1,0,0,-1}, 36.8333333333}, {-773.5, {-1,0,1,0}, 478.833333333}, {-773.5, {-1,0,1,0}, 73.6666666667}, {-1547.0, {-1,0,3,2}, 4.0}, }},
      {8, { {-386.75, {-1,0,0,1}, 36.8333333333}, {-773.5, {1,0,-1,0}, 478.833333333}, {-773.5, {1,0,-1,0}, 73.6666666667}, {1547.0, {-1,0,3,2}, 4.0}, }},
      {9, { {-2707.25, {2,0,-3,1}, 110.5}, {-386.75, {-1,0,0,1}, 718.25}, {-386.75, {-1,0,0,1}, 110.5}, {1547.0, {1,0,3,2}, 6.0}, }},
      {12, { {-386.75, {-1,0,0,1}, 36.8333333333}, {386.75, {0,0,-1,1}, 239.416666667}, {386.75, {0,0,-1,1}, 36.8333333333}, {-1547.0, {1,0,3,-2}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,0,-1},{1,1,0,1}},
    {
      {1, { {-386.75, {1,0,0,-1}, 36.8333333333}, {-773.5, {0,1,0,-1}, 478.833333333}, {-773.5, {0,1,0,-1}, 73.6666666667}, {-1547.0, {2,3,0,-1}, 4.0}, }},
      {2, { {386.75, {-1,3,0,-2}, 110.5}, {-386.75, {1,0,0,-1}, 718.25}, {-386.75, {1,0,0,-1}, 110.5}, {1547.0, {2,3,0,1}, 6.0}, }},
      {3, { {-386.75, {1,0,0,-1}, 36.8333333333}, {386.75, {-1,1,0,0}, 239.416666667}, {386.75, {-1,1,0,0}, 36.8333333333}, {1547.0, {-2,3,0,1}, 2.0}, }},
      {8, { {-386.75, {-1,0,0,1}, 36.8333333333}, {386.75, {1,-1,0,0}, 239.416666667}, {386.75, {1,-1,0,0}, 36.8333333333}, {-1547.0, {-2,3,0,1}, 2.0}, }},
      {9, { {386.75, {1,-3,0,2}, 110.5}, {-386.75, {-1,0,0,1}, 718.25}, {-386.75, {-1,0,0,1}, 110.5}, {-1547.0, {2,3,0,1}, 6.0}, }},
      {10, { {-386.75, {-1,0,0,1}, 36.8333333333}, {-773.5, {0,-1,0,1}, 478.833333333}, {-773.5, {0,-1,0,1}, 73.6666666667}, {1547.0, {2,3,0,-1}, 4.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,-2,-3},{1,1,1,1}},
    {
      {8, { {-6961.5, {-5,-2,1,6}, 221.0}, {-773.5, {1,0,-1,0}, 478.833333333}, {-773.5, {1,0,-1,0}, 73.6666666667}, {0.0, {-2,1,4,3}, 6.0}, }},
      {12, { {-2320.5, {-1,-1,1,1}, 73.6666666667}, {-773.5, {1,-1,-1,1}, 478.833333333}, {-773.5, {1,-1,-1,1}, 73.6666666667}, {0.0, {-1,3,3,-1}, 4.0}, }},
      {14, { {-6961.5, {-6,-1,2,5}, 221.0}, {-773.5, {0,-1,0,1}, 478.833333333}, {-773.5, {0,-1,0,1}, 73.6666666667}, {-0.0, {3,4,1,-2}, 6.0}, }},
      {1, { {-6961.5, {6,1,-2,-5}, 221.0}, {-773.5, {0,1,0,-1}, 478.833333333}, {-773.5, {0,1,0,-1}, 73.6666666667}, {0.0, {3,4,1,-2}, 6.0}, }},
      {3, { {-2320.5, {1,1,-1,-1}, 73.6666666667}, {-773.5, {-1,1,1,-1}, 478.833333333}, {-773.5, {-1,1,1,-1}, 73.6666666667}, {-0.0, {-1,3,3,-1}, 4.0}, }},
      {7, { {-6961.5, {5,2,-1,-6}, 221.0}, {-773.5, {-1,0,1,0}, 478.833333333}, {-773.5, {-1,0,1,0}, 73.6666666667}, {0.0, {-2,1,4,3}, 6.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,0,1},{1,1,0,1}},
    {
      {1, { {386.75, {1,0,0,-1}, 36.8333333333}, {0.0, {0,1,0,-1}, 478.833333333}, {0.0, {0,1,0,-1}, 73.6666666667}, {-773.5, {2,3,0,-1}, 4.0}, }},
      {2, { {1933.75, {-1,3,0,-2}, 110.5}, {386.75, {1,0,0,-1}, 718.25}, {386.75, {1,0,0,-1}, 110.5}, {773.5, {2,3,0,1}, 6.0}, }},
      {3, { {386.75, {1,0,0,-1}, 36.8333333333}, {386.75, {-1,1,0,0}, 239.416666667}, {386.75, {-1,1,0,0}, 36.8333333333}, {773.5, {-2,3,0,1}, 2.0}, }},
      {8, { {386.75, {-1,0,0,1}, 36.8333333333}, {386.75, {1,-1,0,0}, 239.416666667}, {386.75, {1,-1,0,0}, 36.8333333333}, {-773.5, {-2,3,0,1}, 2.0}, }},
      {9, { {1933.75, {1,-3,0,2}, 110.5}, {386.75, {-1,0,0,1}, 718.25}, {386.75, {-1,0,0,1}, 110.5}, {-773.5, {2,3,0,1}, 6.0}, }},
      {10, { {386.75, {-1,0,0,1}, 36.8333333333}, {0.0, {0,-1,0,1}, 478.833333333}, {0.0, {0,-1,0,1}, 73.6666666667}, {773.5, {2,3,0,-1}, 4.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,1,2,3},{1,0,1,1}},
    {
      {1, { {1160.25, {1,0,0,-1}, 36.8333333333}, {386.75, {0,0,1,-1}, 239.416666667}, {386.75, {0,0,1,-1}, 36.8333333333}, {0.0, {1,0,3,-2}, 2.0}, }},
      {4, { {3480.75, {-2,0,3,-1}, 110.5}, {1160.25, {1,0,0,-1}, 718.25}, {1160.25, {1,0,0,-1}, 110.5}, {0.0, {1,0,3,2}, 6.0}, }},
      {5, { {1160.25, {1,0,0,-1}, 36.8333333333}, {773.5, {-1,0,1,0}, 478.833333333}, {773.5, {-1,0,1,0}, 73.6666666667}, {-0.0, {-1,0,3,2}, 4.0}, }},
      {8, { {1160.25, {-1,0,0,1}, 36.8333333333}, {773.5, {1,0,-1,0}, 478.833333333}, {773.5, {1,0,-1,0}, 73.6666666667}, {0.0, {-1,0,3,2}, 4.0}, }},
      {9, { {3480.75, {2,0,-3,1}, 110.5}, {1160.25, {-1,0,0,1}, 718.25}, {1160.25, {-1,0,0,1}, 110.5}, {-0.0, {1,0,3,2}, 6.0}, }},
      {12, { {1160.25, {-1,0,0,1}, 36.8333333333}, {386.75, {0,0,-1,1}, 239.416666667}, {386.75, {0,0,-1,1}, 36.8333333333}, {-0.0, {1,0,3,-2}, 2.0}, }},
    }});
  LAYOUT_VALID_TO_LATCOMB_CONSTS_ENCODER_REAL.push_back({
    {{0,-1,0,1},{1,1,0,1}},
    {
      {1, { {386.75, {1,0,0,-1}, 36.8333333333}, {773.5, {0,1,0,-1}, 478.833333333}, {773.5, {0,1,0,-1}, 73.6666666667}, {1547.0, {2,3,0,-1}, 4.0}, }},
      {2, { {-386.75, {-1,3,0,-2}, 110.5}, {386.75, {1,0,0,-1}, 718.25}, {386.75, {1,0,0,-1}, 110.5}, {-1547.0, {2,3,0,1}, 6.0}, }},
      {3, { {386.75, {1,0,0,-1}, 36.8333333333}, {-386.75, {-1,1,0,0}, 239.416666667}, {-386.75, {-1,1,0,0}, 36.8333333333}, {-1547.0, {-2,3,0,1}, 2.0}, }},
      {8, { {386.75, {-1,0,0,1}, 36.8333333333}, {-386.75, {1,-1,0,0}, 239.416666667}, {-386.75, {1,-1,0,0}, 36.8333333333}, {1547.0, {-2,3,0,1}, 2.0}, }},
      {9, { {-386.75, {1,-3,0,2}, 110.5}, {386.75, {-1,0,0,1}, 718.25}, {386.75, {-1,0,0,1}, 110.5}, {1547.0, {2,3,0,1}, 6.0}, }},
      {10, { {386.75, {-1,0,0,1}, 36.8333333333}, {773.5, {0,-1,0,1}, 478.833333333}, {773.5, {0,-1,0,1}, 73.6666666667}, {-1547.0, {2,3,0,-1}, 4.0}, }},
    }});
}
