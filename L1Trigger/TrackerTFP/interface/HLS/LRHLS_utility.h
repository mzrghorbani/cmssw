/**
 * Author Maziar Ghorbani - Brunel University London
 * Based on implementation of Linear-Regression algorithm by Dr Thomas Schuh
 *
 * Utility functions for analyzing fixed-point variables
 * Class 'info' and 'CheckCalc' are based on Dr Ian Tomalin implementation of Kalman filter
 */

#ifndef L1Trigger_TrackerTFP_LRHLS_UTILITY_H_
#define L1Trigger_TrackerTFP_LRHLS_UTILITY_H_

#ifdef CMSSW_GIT_HASH
#include <iostream>
#else
#include <iostream>
#endif

#ifndef __SYNTHESIS__
#define PRINT
#endif

#ifdef PRINT
#define PRINT_SUMMARY
#endif

#ifdef PRINT_SUMMARY
#include <iostream>
#include <string>
#include <typeinfo>
#include <map>
#include <cmath>
#include <utility>
#endif

#ifdef CMSSW_GIT_HASH
namespace trackerTFP {

namespace trackerHLS {
#endif

#ifdef PRINT_SUMMARY
// Utility to compare real and fixed-point variables
namespace CHECK_AP {

struct INFO {
  INFO() {}
  ~INFO() {}
  INFO(std::string className, int intBitsCfg, int intBitsSeenHigh, int intBitsSeenLow) :
      className_(std::move(className)), intBitsCfg_(intBitsCfg), intBitsSeenHigh_(intBitsSeenHigh),
      intBitsSeenLow_(intBitsSeenLow) {}
  std::string className_;
  int intBitsCfg_;
  int intBitsSeenHigh_;
  int intBitsSeenLow_;
};

// Function checkCalc prints Loss of Precision and word-length errors
template<class T>
bool checkCalc(const std::string &varName, T res_fix, double res_float, double reltol = 0.1, double tol = 0.0) {
  bool OK = true;
  static unsigned int nErrors = 0;
  std::map<std::string, CHECK_AP::INFO> apCheckMap_;
  double res_float_abs = fabs(res_float);
  int intBitsSeen = (res_float_abs > 0) ? std::ceil(log(res_float_abs) / log(2.)) : -99;
  std::string cNameUgly = typeid(T).name();
  std::string cName = "unknown";
  int intBitsCfg = 0;
  if (cNameUgly.find("ap_int") != std::string::npos) {
    cName = "ap_int";
    intBitsCfg = res_fix.width - 1;
  } else if (cNameUgly.find("ap_uint") != std::string::npos) {
    cName = "ap_uint";
    intBitsCfg = res_fix.width;
  } else if (cNameUgly.find("ap_fixed") != std::string::npos) {
    cName = "ap_fixed";
    intBitsCfg = res_fix.iwidth - 1;
  } else if (cNameUgly.find("ap_ufixed") != std::string::npos) {
    cName = "ap_ufixed";
    intBitsCfg = res_fix.iwidth;
  }
  if (apCheckMap_.find(varName) == apCheckMap_.end()) {
    apCheckMap_[varName] = CHECK_AP::INFO(cName, intBitsCfg, intBitsSeen, intBitsSeen);
  } else {
    int intHighOld = apCheckMap_[varName].intBitsSeenHigh_;
    int intLowOld = apCheckMap_[varName].intBitsSeenLow_;
    if (intHighOld < intBitsSeen) {
      apCheckMap_[varName] = CHECK_AP::INFO(cName, intBitsCfg, intBitsSeen, intLowOld);
    } else if (intLowOld > intBitsSeen) {
      apCheckMap_[varName] = CHECK_AP::INFO(cName, intBitsCfg, intHighOld, intBitsSeen);
    }
  }
  if (res_float < 0) {
    if (cName == "ap_uint" || cName == "ap_ufixed") {
      OK = false;
      if (nErrors < 999)
        std::cout << "checkCalc SIGN ERROR (" << cName << "): " << varName << " " << res_float << std::endl;
    }
  }
  if (res_float != 0.) {
    double err = fabs(double(res_fix) - res_float);
    double relerr = err / fabs(res_float);
    if (relerr > reltol && err > tol) {
      OK = false;
      if (nErrors < 999) {
        if (intBitsSeen > intBitsCfg) {
          std::cout << "checkCalc TOO FEW INTEGER BITS (" << cName << "): " <<
                    varName << " " << intBitsSeen << " > " << intBitsCfg << std::endl;
        } else {
          std::cout << "checkCalc TOO FEW FRACTIONAL BITS (" << cName << "): " <<
                    varName << " relerr=" << relerr << " err=" << err << " fix=" << res_fix << " flt="
                    << res_float << std::endl;
        }
      }
    }
  }
  if (not OK)
    nErrors++;
  return OK;
}

}
#endif

#ifdef CMSSW_GIT_HASH
}

}
#endif

#endif