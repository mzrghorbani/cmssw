/**
 * Author Maziar Ghorbani - Brunel University London
 * Based on implementation of Linear-Regression algorithm by Dr Thomas Schuh
 *
 * All-Types and
 * Class 'info' and 'CheckCalc' are based on Dr Ian Tomalin implementation of Kalman filter
 */

#ifndef L1Trigger_TrackerTFP_LRHLS_TYPES_H_
#define L1Trigger_TrackerTFP_LRHLS_TYPES_H_

#ifdef CMSSW_GIT_HASH
#include <iostream>
#include <ap_fixed.h>
#include <hls_stream.h>
#include "LRHLS_utility.h"
#else
#include <ap_fixed.h>
#include <hls_stream.h>
#include "LRHLS_utility.h"
#endif

#ifdef CMSSW_GIT_HASH
namespace trackerTFP {

namespace trackerHLS {
#endif

enum {
  B9 = 9, B1 = 1
};

// Parameters for LR-HLS algorithm
#define lrhlsNumStubs 12
#define lrhlsNumLayers 7

typedef ap_uint<B9> uint9_t;
typedef ap_uint<B1> uint1_t;

// Physical parameters and reciprocals
static const unsigned int lrhlsMinLayers = 4;
static const unsigned int lrhlsMinLayersPS = 2;
static const double lrhlsChosenRofZ = 50.0;
static const double lrhlsChosenRofPhi = 67.24;
static const double lrhlsResidPhi = 1. / 0.001;
static const double lrhlsResidZPS = 1. / 0.07;
static const double lrhlsResidZ2S = 1. / 2.5;

// Struct to hold stub parameters
struct LRStub {
  ap_fixed<18, 7> r = 0;
  ap_fixed<18, 3> phi = 0;
  ap_fixed<18, 7> z = 0;
  ap_uint<3> layer = 0;
  ap_uint<1> barrel = 0;
  ap_uint<1> psModule = 0;
  ap_uint<1> valid = 0;
};

typedef hls::stream<LRStub> strm_t;

// Struct to hold track parameters
struct LRTrack {
  ap_fixed<18, 1> qOverPt = 0;
  ap_fixed<18, 3> phiT = 0;
  ap_fixed<18, 1> cot = 0;
  ap_fixed<18, 3> zT = 0;
};

#ifdef CMSSW_GIT_HASH
}

}
#endif

#endif
