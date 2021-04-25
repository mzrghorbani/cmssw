#include "L1Trigger/TrackerTFP/interface/LinearRegressionHLS.h"

using namespace std;
using namespace edm;
using namespace trackerDTC;

#include <numeric>
#include <vector>
#include <iterator>
#include <algorithm>

namespace trackerTFP {

LinearRegressionHLS::LinearRegressionHLS(const ParameterSet& iConfig, const Setup* setup, const DataFormats* dataFormats, int region) :
    setup_(setup),
    dataFormats_(dataFormats),
    region_(region),
    input_(dataFormats_->numChannel(Process::mht)) {}

void LinearRegressionHLS::consume(const TTDTC::Streams& streams) {
  auto valid = [](int& sum, const TTDTC::Frame& frame){ return sum += (frame.first.isNonnull() ? 1 : 0); };
  int nStubsMHT(0);
  for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
    const TTDTC::Stream& stream = streams[region_ * dataFormats_->numChannel(Process::mht) + channel];
    nStubsMHT += accumulate(stream.begin(), stream.end(), 0, valid);
  }
  stubsMHT_.reserve(nStubsMHT);
  //stubsLRHLS_.reserve(nStubsMHT);
  for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
    const TTDTC::Stream& stream = streams[region_ * dataFormats_->numChannel(Process::mht) + channel];
    vector<StubMHT*>& input = input_[channel];
    input.reserve(stream.size());
    for (const TTDTC::Frame& frame : stream) {
      StubMHT* stub = nullptr;
      if (frame.first.isNonnull()) {
        stubsMHT_.emplace_back(frame, dataFormats_);
        stub = & stubsMHT_.back();
      }
      input.push_back(stub);
    }
  }
}

void LinearRegressionHLS::produce(TTDTC::Streams& stubs, TTDTC::Streams& tracks) {
  for (int channel = 0; channel < dataFormats_->numChannel(Process::mht); channel++) {
    const vector<StubMHT*>& input = input_[channel];
    TTDTC::Stream& streamStubs = stubs[region_ * dataFormats_->numChannel(Process::mht) + channel];
    TTDTC::Stream& streamTracks = tracks[region_ * dataFormats_->numChannel(Process::mht) + channel];
    streamStubs.reserve(input.size());
    streamTracks.reserve(input.size());
    for (auto it = input.begin(); it != input.end();) {
      if (!*it) {
        streamStubs.emplace_back(TTDTC::Frame());
        streamTracks.emplace_back(TTDTC::Frame());
        it++;
        continue;
      }
      const auto start = it;
      const int id = (*it)->trackId();
      auto different = [id](StubMHT* stub){ return !stub || id != stub->trackId(); };
      it = find_if(it, input.end(), different);
      vector<StubMHT*> stubsMHT;
      vector<StubLRHLS> stubsLRHLS;
      StubLRHLS* track = nullptr;
      stubsMHT.reserve(distance(start, it));
      copy(start, it, back_inserter(stubsMHT));

      // hls::stream declarations for LR-HLS input/output data
      trackerHLS::strm_t strm_in;
      trackerHLS::strm_t strm_out;

      // Reading data from input stream and writing to LR-HLS hls::stream
      for(StubMHT* stubMHT : stubsMHT) {
        trackerHLS::LRStub stubLRHLS;
        stubLRHLS.r = stubMHT->r();
        stubLRHLS.phi = stubMHT->phi();
        stubLRHLS.z = stubMHT->z();
        stubLRHLS.barrel = stubMHT->barrel();
        stubLRHLS.psModule = stubMHT->psModule();
        stubLRHLS.layer = stubMHT->layer();
        stubLRHLS.valid = 1;
        strm_in.write(stubLRHLS);
#ifdef PRINT_SUMMARY
        trackerHLS::CHECK_AP::checkCalc("stub.r", stubLRHLS.r, stubMHT->r(), 0.05);
        trackerHLS::CHECK_AP::checkCalc("stub.phi", stubLRHLS.phi, stubMHT->phi(), 0.05);
        trackerHLS::CHECK_AP::checkCalc("stub.z", stubLRHLS.z, stubMHT->z(), 0.05);
#endif
      }

      // Storing hls::stream length for variable-sized storage
      trackerHLS::uint9_t strm_len = strm_in.size();

      // Calling LR-HLS Top-Level function
      trackerHLS::LRTrack trackLR = trackerHLS::LRHLS_top(strm_in, strm_out, strm_len);

      // Reading processed data from LR-HLS and writing it into LinearRegressionHLS module
      for(trackerHLS::uint9_t i = 0; i < strm_len; i++) {
        trackerHLS::LRStub stub = strm_out.read();
        double phiT = trackLR.phiT;
        double qOverPt = trackLR.qOverPt;
        double zT = trackLR.zT;
        double cot = trackLR.cot;
        StubLRHLS stubLRHLS(*stubsMHT.at(i), phiT, qOverPt, zT, cot);
        stubsLRHLS.push_back(stubLRHLS);
        track = &stubLRHLS;
      }

      // Storing data into output stream
      for (StubLRHLS stubLRHLS : stubsLRHLS)
        streamStubs.emplace_back(stubLRHLS ? stubLRHLS.frame() : TTDTC::Frame());
      streamTracks.insert(streamTracks.end(), stubsLRHLS.size(), track ? track->frame() : TTDTC::Frame());
    }
  }
}

} // namespace trackerTFP