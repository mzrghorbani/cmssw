/**
 * Author Maziar Ghorbani - Brunel University London
 * Based on implementation of Linear-Regression algorithm by Dr Thomas Schuh
 *
 * Top-level function declaration for LR-HLS class
 */

#ifndef L1Trigger_TrackerTFP_LRHLS_TOP_H_
#define L1Trigger_TrackerTFP_LRHLS_TOP_H_

#ifdef CMSSW_GIT_HASH
#include "L1Trigger/TrackerTFP/interface/HLS/LRHLS_utility.h"
#include "L1Trigger/TrackerTFP/interface/HLS/LRHLS_types.h"
#include "L1Trigger/TrackerTFP/interface/HLS/LRHLS.h"
#else
#include "LRHLS_utility.h"
#include "LRHLS_types.h"
#include "LRHLS.h"
#endif

#ifdef CMSSW_GIT_HASH
namespace trackerTFP {

namespace trackerHLS {
#endif

// Top-level function declaration for LR-HLS class
LRTrack LRHLS_top(strm_t& strm_in, strm_t& strm_out, const uint9_t& strm_len);

#ifdef CMSSW_GIT_HASH
}

}
#endif

#endif