/*
Part of the Fluid Corpus Manipulation Project (http://www.flucoma.org/)
Copyright 2017-2019 University of Huddersfield.
Licensed under the BSD-3 License.
See license.md file in the project root for full license information.
This project has received funding from the European Research Council (ERC)
under the European Union’s Horizon 2020 research and innovation programme
(grant agreement No 725899).
*/
#pragma once

#include "algorithms/GraphPlay.hpp"
#include "algorithms/public/MelBands.hpp"
#include "clients/common/BufferedProcess.hpp"
#include "clients/common/FluidBaseClient.hpp"
#include "clients/nrt/NRTClient.hpp"
#include "clients/common/ParameterConstraints.hpp"
#include "clients/common/ParameterSet.hpp"
#include "clients/common/ParameterTrackChanges.hpp"
#include "clients/common/ParameterTypes.hpp"
#include "clients/common/BufferAdaptor.hpp"
#include <clients/common/Result.hpp>

namespace fluid {
namespace client {
namespace graphplay {
  enum GraphPlayParamIndex {
    kSourceBuf,
    kNumBands,
    kThreshold,
    kMinDur,
    kMinDist,
    kForget,
    kStart,
    kOutputBuffer,
    kFFT,
    kMaxFFTSize
  };

  constexpr auto GraphPlayParams = defineParameters(
                  InputBufferParam("source", "Source Buffer"),
                  LongParam("numBands", "Number of Mel bands", 64),
                  FloatParam("threshold", "Threshold", 0.3, Min(0), Max(1.0)),
                  LongParam("minDur", "Min duration (frames)", 10, Min(1)),
                  LongParam("minDist", "Min distance (frames)", 10, Min(1)),
                  LongParam("forget", "Forget time (frames)", 1, Min(1)),
                  FloatParam("start", "Start point", 0, Min(0), Max(1)),
                  BufferParam("outputBuffer","Actual start/end points"),
                  FFTParam<kMaxFFTSize>("fftSettings", "FFT Settings",
                                             2048, 512, -1),
                  LongParam<Fixed<true>>("maxFFTSize", "Maxiumm FFT Size",
                                              16384, Min(4), PowerOfTwo{})
  );

  constexpr auto STFTParams = defineParameters(
      FFTParam<kMaxFFTSize>("fftSettings", "FFT Settings", 1024, -1, -1),
      LongParam<Fixed<true>>("maxFFTSize", "Maxiumm FFT Size", 16384, Min(4),
                             PowerOfTwo{})
  );


class GraphPlayClient : public FluidBaseClient, public AudioOut, ModelObject {

public:
  using ParamDescType = decltype(GraphPlayParams);

  using ParamSetViewType = ParameterSetView<ParamDescType>;
  std::reference_wrapper<ParamSetViewType> mParams;

  using STFTParamDescType = decltype(STFTParams);
  using STFTParamSetType = ParameterSet<STFTParamDescType>;
  STFTParamSetType mSTFTParams{STFTParams};


  void setParams(ParamSetViewType& p) { mParams = p; }

  template <size_t N>
  auto& get() const
  {
    return mParams.get().template get<N>();
  }

  static constexpr auto& getParameterDescriptors()
  {
    return GraphPlayParams;
  }


  GraphPlayClient(ParamSetViewType &p)
      : mParams{p}, mSTFTProcessor{get<kMaxFFTSize>(), 0, 1} {
    audioChannelsIn(0);
    audioChannelsOut(1);
  }

  index latency() { return get<kFFT>().winSize(); }

  void reset() { mSTFTProcessor.reset();}

  MessageResult<void> analyze(){
    using namespace algorithm;
    algorithm::GraphPlay newAlgoritm;
    auto source = BufferAdaptor::ReadAccess(get<kSourceBuf>().get());
    double sampleRate = source.sampleRate();
    auto fftParams = get<kFFT>();
    if(!source.exists())
      return {Result::Status::kError, "Source Buffer Supplied But Invalid"};
    index srcFrames = source.numFrames();
    if (srcFrames <= 0)
      return {Result::Status::kError, "Empty source buffer"};
    RealVector srcTmp{source.samps(0, srcFrames, 0)};

    RealVector outputData(1);

    newAlgoritm.init(srcTmp, sampleRate,
                get<kFFT>().winSize(),
                get<kFFT>().fftSize(),
                get<kFFT>().hopSize(),
                get<kNumBands>(),
                7,
                get<kThreshold>(),
                outputData
    );
    mNewAlgorithm = newAlgoritm;
    mNewAlgorithmReady = true;
    return OK();
  }


  template <typename T>
  void process(std::vector<HostVector<T>> &,
               std::vector<HostVector<T>> &output, FluidContext &c) {
    assert(audioChannelsOut() && "No control channels");
    assert(output.size() >= asUnsigned(audioChannelsOut()) &&
           "Too few output channels");
    if (mNewAlgorithmReady) {
      std::swap(mAlgorithm, mNewAlgorithm);
      mNewAlgorithmReady = false;
      mSTFTParams.template get<0>() = FFTParams(
        mAlgorithm.mWindowSize, mAlgorithm.mHopSize, mAlgorithm.mFFTSize);
    }
    RealVector outputData(2);
    auto outBuf = BufferAdaptor::Access(get<kOutputBuffer>().get());
    bool validOutput = (outBuf.exists() && outBuf.numFrames() == 2);
    mSTFTProcessor.processOutput(
          mSTFTParams, output, c,
          [&](ComplexMatrixView out) {
            if(mAlgorithm.initialized()){
              mAlgorithm.processFrame(out.row(0), get<kStart>(),
              get<kThreshold>(), get<kMinDur>(), get<kMinDist>(),
              get<kForget>(), outputData);
              if(validOutput) outBuf.samps(0) = outputData;
            }
          });
    }

    static auto getMessageDescriptors()
    {
      return defineMessages(
        makeMessage("analyze", &GraphPlayClient::analyze)
      );
    }

private:
  STFTBufferedProcess<STFTParamSetType, 0, true> mSTFTProcessor;
  algorithm::GraphPlay mAlgorithm;
  algorithm::GraphPlay mNewAlgorithm;
  bool mNewAlgorithmReady{false};


};
}
using RTGraphPlayClient = ClientWrapper<graphplay::GraphPlayClient>;
} // namespace client
} // namespace fluid
