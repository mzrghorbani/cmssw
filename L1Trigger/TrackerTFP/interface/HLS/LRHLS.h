/**
 * Author Maziar Ghorbani - Brunel University London
 * Based on implementation of Linear-Regression algorithm by Dr Thomas Schuh
 *
 * LR-HLS class declaration for
 * Processing stub parameters
 * calculating Linear-Regression fit
 * Find largest residuals to the fit
 * Remove the stubs with largest residuals
 */

#ifndef L1Trigger_TrackerTFP_LRHLS_H_
#define L1Trigger_TrackerTFP_LRHLS_H_

#ifdef CMSSW_GIT_HASH
#include "L1Trigger/TrackerTFP/interface/HLS/LRHLS_types.h"
#else

#include "LRHLS_types.h"

#endif

#ifdef CMSSW_GIT_HASH
namespace trackerTFP {

namespace trackerHLS {
#endif

class LRHLS {
public:
  LRHLS() {}
  ~LRHLS() {}

  LRStub *fit(LRStub stubs[lrhlsNumStubs], LRTrack &track);
  void initFit();
  uint1_t exit_t();
  void countStubsPS(const LRStub *stubs);
  void countStubsAll(const LRStub *stubs);
  void countLayers(const LRStub *stubs);
  uint1_t calcHelix(const LRStub *stubs);
  uint1_t findLargestResid(const LRStub *stubs);
  void killLargestResid(LRStub *stubs);

public:
  ap_fixed<32, 4> qOverPt_;
  ap_fixed<32, 9> phiT_;
  ap_fixed<26, 11> cot_;
  ap_fixed<26, 16> zT_;
  ap_uint<4> largestResid_;
  ap_uint<4> allStubs_;
  ap_uint<4> psStubs_;
  ap_uint<1> valid_;
  ap_uint<3> layerPopulation_[lrhlsNumLayers];
};

#ifdef CMSSW_GIT_HASH
}

}
#endif

#endif