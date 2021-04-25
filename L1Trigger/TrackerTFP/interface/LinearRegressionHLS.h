#ifndef L1Trigger_TrackerTFP_LinearRegressionHLS_h
#define L1Trigger_TrackerTFP_LinearRegressionHLS_h

#include "L1Trigger/TrackerDTC/interface/Setup.h"
#include "L1Trigger/TrackerTFP/interface/DataFormats.h"
#include "L1Trigger/TrackerTFP/interface/HLS/LRHLS_types.h"
#include "L1Trigger/TrackerTFP/interface/HLS/LRHLS_top.h"

namespace trackerTFP {

// Class to find in a region finer rough candidates in r-phi
class LinearRegressionHLS {
public:
  LinearRegressionHLS(const edm::ParameterSet& iConfig, const trackerDTC::Setup* setup, const DataFormats* dataFormats, int region);
  ~LinearRegressionHLS(){}

  // read in and organize input product
  void consume(const TTDTC::Streams& streams);
  // fill output products
  void produce(TTDTC::Streams& stubs, TTDTC::Streams& tracks);

private:
  //
  const trackerDTC::Setup* setup_;
  //
  const DataFormats* dataFormats_;
  //
  int region_;
  //
  std::vector<std::vector<StubMHT*>> input_;
  //
  std::vector<StubMHT> stubsMHT_;
  //
  std::vector<StubLRHLS> stubsLRHLS_;
  //
  bool valid_;
  //
  double phiT_;
  //
  double qOverPt_;
  //
  double zT_;
  //
  double cot_;
  //
};

}

#endif