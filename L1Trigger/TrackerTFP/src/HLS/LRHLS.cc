/**
 * Author Maziar Ghorbani - Brunel University London
 * Class LR-HLS implementation
 * Based on implementation of Linear-Regression algorithm by Dr Thomas Schuh
 *
 * Class LR-HLS definitions
 */

#ifdef CMSSW_GIT_HASH
#include "L1Trigger/TrackerTFP/interface/HLS/LRHLS.h"
#else

#include "LRHLS.h"

#endif

#ifdef CMSSW_GIT_HASH
namespace trackerTFP {

namespace trackerHLS {
#endif

// Main loop to process stub data and calculate LR-HLS fit parameters
LRStub *LRHLS::fit(LRStub *stubs, LRTrack &track) {
#ifdef PRINT_SUMMARY
  double lrhlsChosenRofZ_ld = 50.0;
  double lrhlsChosenRofPhi_ld = 67.24;
  double lrhlsResidPhi_ld = 1. / 0.001;
  double lrhlsResidZPS_ld = 1. / 0.07;
  double lrhlsResidZ2S_ld = 1. / 2.5;
  CHECK_AP::checkCalc("lrhlsChosenRofZ", ap_ufixed<16, 7>(lrhlsChosenRofZ), lrhlsChosenRofZ_ld, 0.05);
  CHECK_AP::checkCalc("lrhlsChosenRofPhi", ap_ufixed<19, 7>(lrhlsChosenRofPhi), lrhlsChosenRofPhi_ld, 0.0001);
  CHECK_AP::checkCalc("lrhlsResidPhi", ap_ufixed<16, 10>(lrhlsResidPhi), lrhlsResidPhi_ld, 0.0001);
  CHECK_AP::checkCalc("lrhlsResidZPS", ap_ufixed<16, 5>(lrhlsResidZPS), lrhlsResidZPS_ld, 0.0001);
  CHECK_AP::checkCalc("lrhlsResidZ2S", ap_ufixed<19, 1>(lrhlsResidZ2S), lrhlsResidZ2S_ld, 0.0001);
#endif

  for (uint9_t k = 0; k < lrhlsNumStubs; k++) {
    initFit();
    countStubsAll(stubs);
    countStubsPS(stubs);
    countLayers(stubs);
    valid_ = calcHelix(stubs);
    valid_ = findLargestResid(stubs);
    if (exit_t() == 1)
      break;
    killLargestResid(stubs);
  }

  track.qOverPt = qOverPt_;
  track.phiT = phiT_;
  track.cot = cot_;
  track.zT = zT_;

  return stubs;
}

// Clearing all class parameters
void LRHLS::initFit() {
  allStubs_ = 0;
  psStubs_ = 0;
  largestResid_ = 0;
  for (uint9_t i = 0; i < lrhlsNumLayers; i++)
    layerPopulation_[i] = 0;
  valid_ = 1;
}

// Finding the total number-of-stubs in PS regions
void LRHLS::countStubsPS(const LRStub *stubs) {
  for (uint9_t i = 0; i < lrhlsNumStubs; i++) {
    LRStub stub = stubs[i];
    if (stub.valid == 1)
      if (stub.psModule == 1)
        psStubs_++;
  }
}

// Finding the total number-of-stubs in all regions
void LRHLS::countStubsAll(const LRStub *stubs) {
  for (uint9_t i = 0; i < lrhlsNumStubs; i++) {
    LRStub stub = stubs[i];
    if (stub.valid == 1)
      allStubs_++;
  }
}

// Finding layer-populations for all stubs
void LRHLS::countLayers(const LRStub *stubs) {
  for (uint9_t i = 0; i < lrhlsNumStubs; i++) {
    LRStub stub = stubs[i];
    if (stub.valid == 1)
      layerPopulation_[stub.layer]++;
  }
}

// Exit condition for terminating the main-loop
uint1_t LRHLS::exit_t() {
  uint1_t single = 1;
  for (uint9_t i = 0; i < lrhlsNumLayers; i++)
    if (layerPopulation_[i] > 1)
      single = 0;
  if (valid_ == 0)
    return 1;
  if (single == 1)
    if (allStubs_ == ap_uint<3>(lrhlsMinLayers) || psStubs_ == ap_uint<2>(lrhlsMinLayersPS))
      return 1;
  return 0;
}

// Calculating track-parameters generated from a fit for PS and all stubs
uint1_t LRHLS::calcHelix(const LRStub *stubs) {
  ap_uint<4> rphi_n = 0;
  ap_uint<4> rz_n = 0;
  ap_fixed<36, 10> rphi_xy = 0;
  ap_fixed<32, 10> rphi_x = 0;
  ap_fixed<36, 14> rphi_y = 0;
  ap_fixed<32, 14> rphi_xx = 0;
  ap_fixed<36, 14> rz_xy = 0;
  ap_fixed<32, 10> rz_x = 0;
  ap_fixed<36, 10> rz_y = 0;
  ap_fixed<32, 14> rz_xx = 0;
#ifdef PRINT_SUMMARY
  unsigned int rphi_n_ld = 0;
  unsigned int rz_n_ld = 0;
  double rphi_xy_ld = 0;
  double rphi_x_ld = 0;
  double rphi_y_ld = 0;
  double rphi_xx_ld = 0;
  double rz_xy_ld = 0;
  double rz_x_ld = 0;
  double rz_y_ld = 0;
  double rz_xx_ld = 0;
#endif

  for (uint9_t i = 0; i < lrhlsNumLayers; i++) {
    if (layerPopulation_[i] != 0) {
      ap_fixed<18, 8> RPhi_layerMinPos = 127;
      ap_fixed<26, 8> Phi_layerMinPos = 127;
      ap_fixed<18, 8> RZ_layerMinPos = 127;
      ap_fixed<18, 8> Z_layerMinPos = 127;
      ap_fixed<18, 8> RPhi_layerMaxPos = -127;
      ap_fixed<26, 8> Phi_layerMaxPos = -127;
      ap_fixed<18, 8> RZ_layerMaxPos = -127;
      ap_fixed<18, 8> Z_layerMaxPos = -127;
#ifdef PRINT_SUMMARY
      double RPhi_layerMinPos_ld = 127;
      double Phi_layerMinPos_ld = 127;
      double RZ_layerMinPos_ld = 127;
      double Z_layerMinPos_ld = 127;
      double RPhi_layerMaxPos_ld = -127;
      double Phi_layerMaxPos_ld = -127;
      double RZ_layerMaxPos_ld = -127;
      double Z_layerMaxPos_ld = -127;
      CHECK_AP::checkCalc("RPhi_layerMinPos", RPhi_layerMinPos, RPhi_layerMinPos_ld, 0.05);
      CHECK_AP::checkCalc("Phi_layerMinPos", Phi_layerMinPos, Phi_layerMinPos_ld, 0.05);
      CHECK_AP::checkCalc("RZ_layerMinPos", RZ_layerMinPos, RZ_layerMinPos_ld, 0.05);
      CHECK_AP::checkCalc("Z_layerMinPos", Z_layerMinPos, Z_layerMinPos_ld, 0.05);
      CHECK_AP::checkCalc("RPhi_layerMaxPos", RPhi_layerMaxPos, RPhi_layerMaxPos_ld, 0.05);
      CHECK_AP::checkCalc("Phi_layerMaxPos", Phi_layerMaxPos, Phi_layerMaxPos_ld, 0.05);
      CHECK_AP::checkCalc("RZ_layerMaxPos", RZ_layerMaxPos, RZ_layerMaxPos_ld, 0.05);
      CHECK_AP::checkCalc("Z_layerMaxPos", Z_layerMaxPos, Z_layerMaxPos_ld, 0.05);
#endif
      uint1_t ps = 0;
      for (uint9_t j = 0; j < lrhlsNumStubs; j++) {
        LRStub stub = stubs[j];
        if (stub.valid == 1) {
          if (stub.layer == i) {
            ap_fixed<22, 10> RZ_pos = (stub.r + ap_ufixed<19, 7>(lrhlsChosenRofPhi)) - ap_ufixed<16, 7>(lrhlsChosenRofZ);
#ifdef PRINT_SUMMARY
            double RZ_pos_ld = double((stub.r + ap_ufixed<19, 7>(lrhlsChosenRofPhi)) - ap_ufixed<16, 7>(lrhlsChosenRofZ));
            CHECK_AP::checkCalc("RZ_pos", RZ_pos, RZ_pos_ld, 0.05);
#endif
            if (stub.psModule == 1) {
              ps = 1;
              if (RPhi_layerMinPos > stub.r)
                RPhi_layerMinPos = stub.r;
              if (Phi_layerMinPos > stub.phi)
                Phi_layerMinPos = stub.phi;
              if (RZ_layerMinPos > RZ_pos)
                RZ_layerMinPos = RZ_pos;
              if (Z_layerMinPos > stub.z)
                Z_layerMinPos = stub.z;
              if (RPhi_layerMaxPos < stub.r)
                RPhi_layerMaxPos = stub.r;
              if (Phi_layerMaxPos < stub.phi)
                Phi_layerMaxPos = stub.phi;
              if (RZ_layerMaxPos < RZ_pos)
                RZ_layerMaxPos = RZ_pos;
              if (Z_layerMaxPos < stub.z)
                Z_layerMaxPos = stub.z;
            } else {
              if (RPhi_layerMinPos > stub.r)
                RPhi_layerMinPos = stub.r;
              if (Phi_layerMinPos > stub.phi)
                Phi_layerMinPos = stub.phi;
              if (RPhi_layerMaxPos < stub.r)
                RPhi_layerMaxPos = stub.r;
              if (Phi_layerMaxPos < stub.phi)
                Phi_layerMaxPos = stub.phi;
            }
          }
        }
      }
      ap_fixed<18, 8> RPhi_pos = (RPhi_layerMinPos + RPhi_layerMaxPos);
      ap_fixed<26, 8> Phi_pos = (Phi_layerMinPos + Phi_layerMaxPos);
      ap_fixed<18, 8> RZ_pos = (RZ_layerMinPos + RZ_layerMaxPos);
      ap_fixed<18, 8> Z_pos = (Z_layerMinPos + Z_layerMaxPos);
#ifdef PRINT_SUMMARY
      double RPhi_pos_ld = double(RPhi_layerMinPos + RPhi_layerMaxPos);
      double Phi_pos_ld = double(Phi_layerMinPos + Phi_layerMaxPos);
      double RZ_pos_ld = double(RZ_layerMinPos + RZ_layerMaxPos);
      double Z_pos_ld = double(Z_layerMinPos + Z_layerMaxPos);
      CHECK_AP::checkCalc("RPhi_pos", RPhi_pos, RPhi_pos_ld, 0.05);
      CHECK_AP::checkCalc("Phi_pos", Phi_pos, Phi_pos_ld, 0.05);
      CHECK_AP::checkCalc("RZ_pos", RZ_pos, RZ_pos_ld, 0.05);
      CHECK_AP::checkCalc("Z_pos", Z_pos, Z_pos_ld, 0.05);
#endif

      ap_fixed<18, 8> RPhi_layerPos = ap_fixed<19, 9>(RPhi_pos) >> 1;
      ap_fixed<32, 8> Phi_layerPos = ap_fixed<32, 9>(Phi_pos) >> 1;
      ap_fixed<18, 8> RZ_layerPos = ap_fixed<19, 9>(RZ_pos) >> 1;
      ap_fixed<18, 8> Z_layerPos = ap_fixed<19, 9>(Z_pos) >> 1;
#ifdef PRINT_SUMMARY
      double RPhi_layerPos_ld = double(RPhi_layerMinPos + RPhi_layerMaxPos) / 2.;
      double Phi_layerPos_ld = double(Phi_layerMinPos + Phi_layerMaxPos) / 2.;
      double RZ_layerPos_ld = double(RZ_layerMinPos + RZ_layerMaxPos) / 2.;
      double Z_layerPos_ld = double(Z_layerMinPos + Z_layerMaxPos) / 2.;
      CHECK_AP::checkCalc("RPhi_layerPos", RPhi_layerPos, RPhi_layerPos_ld, 0.05);
      CHECK_AP::checkCalc("Phi_layerPos", Phi_layerPos, Phi_layerPos_ld, 0.05);
      CHECK_AP::checkCalc("RZ_layerPos", RZ_layerPos, RZ_layerPos_ld, 0.05);
      CHECK_AP::checkCalc("Z_layerPos", Z_layerPos, Z_layerPos_ld, 0.05);
#endif
      rphi_n++;
      rphi_xy += (RPhi_layerPos * Phi_layerPos);
      rphi_x += RPhi_layerPos;
      rphi_y += Phi_layerPos;
      rphi_xx += (RPhi_layerPos * RPhi_layerPos);
      if (ps == 1) {
        rz_n++;
        rz_xy += (RZ_layerPos * Z_layerPos);
        rz_x += RZ_layerPos;
        rz_y += Z_layerPos;
        rz_xx += (RZ_layerPos * RZ_layerPos);
      }
#ifdef PRINT_SUMMARY
      rphi_n_ld++;
      rphi_xy_ld += double(RPhi_layerPos * Phi_layerPos);
      rphi_x_ld += double(RPhi_layerPos);
      rphi_y_ld += double(Phi_layerPos);
      rphi_xx_ld += double(RPhi_layerPos * RPhi_layerPos);
      if (ps == 1) {
        rz_n_ld++;
        rz_xy_ld += double(RZ_layerPos * Z_layerPos);
        rz_x_ld += double(RZ_layerPos);
        rz_y_ld += double(Z_layerPos);
        rz_xx_ld += double(RZ_layerPos * RZ_layerPos);
      }
      CHECK_AP::checkCalc("rphi_xy", rphi_xy, rphi_xy_ld, 0.05);
      CHECK_AP::checkCalc("rphi_x", rphi_x, rphi_x_ld, 0.05);
      CHECK_AP::checkCalc("rphi_y", rphi_y, rphi_y_ld, 0.05);
      CHECK_AP::checkCalc("rphi_xx", rphi_xx, rphi_xx_ld, 0.05);
      CHECK_AP::checkCalc("rz_xy", rphi_xy, rphi_xy_ld, 0.05);
      CHECK_AP::checkCalc("rz_x", rphi_x, rphi_x_ld, 0.05);
      CHECK_AP::checkCalc("rz_y", rphi_y, rphi_y_ld, 0.05);
      CHECK_AP::checkCalc("rz_xx", rphi_xx, rphi_xx_ld, 0.05);
#endif
    }
  }
  ap_fixed<38, 17> rphi_den = (rphi_n * rphi_xx) - (rphi_x * rphi_x);
  if(rphi_den == 0)
    return 0;
  qOverPt_ = ((rphi_n * rphi_xy) - (rphi_x * rphi_y)) / rphi_den;
  phiT_ = ((rphi_y * rphi_xx) - (rphi_x * rphi_xy)) / rphi_den;
  ap_fixed<38, 17> rz_den = (rz_n * rz_xx) - (rz_x * rz_x);
  if(rz_den == 0)
    return 0;
  cot_ = ((rz_n * rz_xy) - (rz_x * rz_y)) / rz_den;
  zT_ = ((rz_y * rz_xx) - (rz_x * rz_xy)) / rz_den;
#ifdef PRINT_SUMMARY
  double rphi_den_ld = (rphi_n * rphi_xx) - (rphi_x * rphi_x);
  double qOverPt_ld = ((rphi_n * rphi_xy) - (rphi_x * rphi_y)) / rphi_den;
  double phiT_ld = ((rphi_y * rphi_xx) - (rphi_x * rphi_xy)) / rphi_den;
  double rz_den_ld = (rz_n * rz_xx) - (rz_x * rz_x);
  double cot_ld = ((rz_n * rz_xy) - (rz_x * rz_y)) / rz_den;
  double zT_ld = ((rz_y * rz_xx) - (rz_x * rz_xy)) / rz_den;
  CHECK_AP::checkCalc("rphi_den", rphi_den, rphi_den_ld, 0.05);
  CHECK_AP::checkCalc("qOverPt", qOverPt_, qOverPt_ld, 0.05);
  CHECK_AP::checkCalc("phiT", phiT_, phiT_ld, 0.05);
  CHECK_AP::checkCalc("rz_den", rz_den, rz_den_ld, 0.05);
  CHECK_AP::checkCalc("cot", cot_, cot_ld, 0.05);
  CHECK_AP::checkCalc("zT", zT_, zT_ld, 0.05);
#endif
  return 1;
}

// Finding the residuals of all stubs to the calculated LR-HLS fit
uint1_t LRHLS::findLargestResid(const LRStub *stubs) {
  ap_fixed<26, 21> phi_largestResid = -1;
  ap_fixed<26, 21> z_largestResid = -1;
  largestResid_ = 0;
  for (uint9_t i = 0; i < lrhlsNumStubs; i++) {
    LRStub stub = stubs[i];
    if (stub.valid == 1) {
      ap_fixed<36, 9> phi_resid = stub.phi - (phiT_ - qOverPt_ * stub.r);
#ifdef PRINT_SUMMARY
      double phi_resid_ld = stub.phi - (phiT_ - qOverPt_ * stub.r);
      CHECK_AP::checkCalc("phi_resid", phi_resid, phi_resid_ld, 0.05);
#endif
      if (phi_resid < 0)
        phi_resid = -phi_resid;

      ap_fixed<26, 14> z_resid = stub.z - (zT_ + cot_ * (stub.r + ap_ufixed<19, 7>(lrhlsChosenRofPhi)) - ap_ufixed<16, 7>(lrhlsChosenRofZ));
#ifdef PRINT_SUMMARY
      double z_resid_ld = stub.z - (zT_ + cot_ * (stub.r + ap_ufixed<19, 7>(lrhlsChosenRofPhi)) - ap_ufixed<16, 7>(lrhlsChosenRofZ));
      CHECK_AP::checkCalc("z_resid", z_resid, z_resid_ld, 0.05);
#endif
      if (z_resid < 0)
        z_resid = -z_resid;
      ap_fixed<32, 18> phi_resid_mult = (phi_resid * ap_ufixed<16, 10>(lrhlsResidPhi));
#ifdef PRINT_SUMMARY
      double phi_resid_mult_ld = (phi_resid * ap_ufixed<16, 10>(lrhlsResidPhi));
      CHECK_AP::checkCalc("phi_resid_mult", phi_resid_mult, phi_resid_mult_ld, 0.05);
#endif
      if(cot_ == 0)
        return 0;
      ap_fixed<26, 11> abs_cot = cot_;
      if (abs_cot < 0)
        abs_cot = -abs_cot;
      ap_fixed<27, 16> z_resid_div;
      if (stub.barrel == 0) {
        z_resid_div = ap_fixed<27, 16>(z_resid) / abs_cot;
#ifdef PRINT_SUMMARY
        double z_resid_div_ld = z_resid / abs_cot;
        CHECK_AP::checkCalc("z_resid_div", z_resid_div, z_resid_div_ld, 0.05);
#endif
      }
      ap_fixed<36, 20> z_resid_mult;
#ifdef PRINT_SUMMARY
      double z_resid_mult_ld;
#endif
      if (stub.psModule == 1) {
        z_resid_mult = (z_resid_div * ap_ufixed<16, 5>(lrhlsResidZPS));
#ifdef PRINT_SUMMARY
        z_resid_mult_ld = (z_resid_div * ap_ufixed<16, 5>(lrhlsResidZPS));
        CHECK_AP::checkCalc("z_resid_mult", z_resid_mult, z_resid_mult_ld, 0.05);
#endif
      } else {
        z_resid_mult = (z_resid_div * ap_ufixed<19, 1>(lrhlsResidZ2S));
#ifdef PRINT_SUMMARY
        z_resid_mult_ld = (z_resid_div * ap_ufixed<19, 1>(lrhlsResidZ2S));
        CHECK_AP::checkCalc("z_resid_mult", z_resid_mult, z_resid_mult_ld, 0.05);
#endif
      }
      ap_fixed<36, 20> phiz_resid = phi_resid_mult + z_resid_mult;
      ap_fixed<36, 20> phiz_largestResid = phi_largestResid + z_largestResid;
      if (phiz_resid > phiz_largestResid) {
#ifdef PRINT_SUMMARY
        double phiz_resid_ld = phi_resid_mult + z_resid_mult;
        double phiz_largest_ld = phi_largestResid + z_largestResid;
        CHECK_AP::checkCalc("phiz_resid", phiz_resid, phiz_resid_ld, 0.05);
        CHECK_AP::checkCalc("phiz_largest", phiz_largestResid, phiz_largest_ld, 0.05);
#endif
        phi_largestResid = phi_resid_mult;
        z_largestResid = z_resid_mult;
        largestResid_ = i;
      }
    }
  }
  return 1;
}

// Killing the stub with the largest-residual to the calculated fit
void LRHLS::killLargestResid(LRStub *stubs) {
  stubs[largestResid_].valid = 0;
}

#ifdef CMSSW_GIT_HASH
}

}
#endif




