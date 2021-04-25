/**
 * Author Maziar Ghorbani - Brunel University London
 * Based on implementation of Linear-Regression algorithm by Dr Thomas Schuh
 *
 * Top-level function definition for LR-HLS class
 */

#ifdef CMSSW_GIT_HASH
#include "L1Trigger/TrackerTFP/interface/HLS/LRHLS_top.h"
#else

#include "LRHLS_top.h"

#endif

#ifdef CMSSW_GIT_HASH
namespace trackerTFP {

namespace trackerHLS {
#endif

// Function to read data from hls::stream and store them in an array
void dummy_process_fe(strm_t& strm_in, LRStub *stubs, const uint9_t& strm_len) {
  uint9_t i;
  for (i = 0; i < strm_len; i++)
    strm_in >> stubs[i];
}

// Function to write processed data form array to hls::stream
void dummy_process_be(LRStub *stubs, strm_t& strm_out, const uint9_t& strm_len) {
  uint9_t i;
  for (i = 0; i < strm_len; i++)
    strm_out << stubs[i];
}

// Top-level function and interfaces necessary for LR-HLS
LRTrack LRHLS_top(strm_t& strm_in, strm_t& strm_out, const uint9_t& strm_len) {

  LRTrack track;

  // The size of array is fixed to constant integral value
  LRStub stubs_in[lrhlsNumStubs];

  // Pointer to first element in output array
  LRStub *stubs_out = &stubs_in[0];

  // Reading data from input hls::stream
  dummy_process_fe(strm_in, stubs_in, strm_len);

  // LR-HLS class object
  LRHLS lrhls;

  // LR-HLS member function generates the LR fit and returns stubs and track parameters
  stubs_out = lrhls.fit(stubs_in, track);

  // Writing data to output hls::stream
  dummy_process_be(stubs_out, strm_out, strm_len);

#ifdef PRINT_SUMMARY
  double qOverPt_top_ld = track.qOverPt;
  double phiT_top_ld = track.phiT;
  double cot_top_ld = track.cot;
  double zT_top_ld = track.zT;
  CHECK_AP::checkCalc("qOverPt_top", track.qOverPt, qOverPt_top_ld, 0.05);
  CHECK_AP::checkCalc("phiT_top", track.phiT, phiT_top_ld, 0.05);
  CHECK_AP::checkCalc("cot_top", track.cot, cot_top_ld, 0.05);
  CHECK_AP::checkCalc("zT_top", track.zT, zT_top_ld, 0.05);
#endif

  return track;
}

#ifdef CMSSW_GIT_HASH
}

}
#endif