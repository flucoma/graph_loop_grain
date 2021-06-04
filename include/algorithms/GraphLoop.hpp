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

#include "algorithms/util/RTPGHI.hpp"
#include "algorithms/util/PeakDetection.hpp"
#include "algorithms/public/DataSetIdSequence.hpp"
#include "algorithms/util/DistanceFuncs.hpp"
#include "algorithms/public/STFT.hpp"
#include "algorithms/public/KDTree.hpp"
#include "algorithms/public/MelBands.hpp"
#include "algorithms/util/AlgorithmUtils.hpp"
#include "algorithms/util/FluidEigenMappings.hpp"
#include "algorithms/GraphPlayUtils.hpp"
#include "data/TensorTypes.hpp"
#include "data/FluidDataSet.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <fstream>

namespace fluid {
namespace algorithm {

class GraphLoop {

public:
  using  MatrixXd = Eigen::MatrixXd;
  using DataSet = FluidDataSet<std::string, double, 1>;

  void init(RealVectorView audio, index sampleRate,
            index windowSize, index fftSize, index hopSize, index numBands,
            index distance, double threshold, bool quantize, RealVectorView output) {
    using namespace Eigen;
    using namespace _impl;
    using namespace std;

    mWindowSize = windowSize;
    mFFTSize = fftSize;
    mHopSize = hopSize;
    mFrameSize = (mFFTSize / 2) + 1;
    mThreshold = threshold;
    mBeat = 0;


    STFT stft = STFT(mWindowSize, mFFTSize, mHopSize);
    mLength = std::floor((audio.size() + mHopSize) / mHopSize);
    mSpectrogram = ComplexMatrix(mLength, mFrameSize);
    stft.process(audio, mSpectrogram);
    mMagnitude = RealMatrix(mLength, mFrameSize);
    stft.magnitude(mSpectrogram, mMagnitude);
    mDM = mUtils.computeDM(mMagnitude, numBands, sampleRate, windowSize,
                           fftSize, distance);
    mDM.diagonal().setZero();

    MatrixXd sim = 1 - mDM.array();
    ArrayXd beatSpectrum = ArrayXd::Zero(mLength);
    for(index i = 0; i < mLength; i++){
      beatSpectrum(i) = sim.diagonal(i).sum() / (mLength - i);
    }
    PeakDetection pd;
    auto bsPeaks = pd.process(beatSpectrum.segment(1,lrint(beatSpectrum.size()/2)), 3, 0, false, true);
    mBeat = bsPeaks[0].first;
    if(bsPeaks.size() > 1 && bsPeaks[1].first < mBeat)mBeat = bsPeaks[1].first;
    if(bsPeaks.size() > 2 && bsPeaks[2].first < mBeat)mBeat = bsPeaks[2].first;
    mFilter.init(5);
    ArrayXd odf = mDM.diagonal(1).array();
    for(index i = 0; i < odf.size(); i++){
      odf(i) = odf(i) - mFilter.processSample(odf(i));
    }
    auto onsets = mPD.process(odf, 0, 0.1, false, false);
    mOnsets = Eigen::VectorXi::Zero(mLength);
    for(index i = 0; i < onsets.size(); i++){
      mOnsets(onsets[i].first) = 1;
    }
    mLoop = RealVector{0, static_cast<double>(mLength)};
    fit(threshold, quantize);
    output(0)  = mLoop(0);
    output(1)  = mLoop(1);
    output(2)  = mBeat;
    output(3)  = mNumLinks;
    mInitialized = true;
  }

  void fit(double threshold, bool quantize){
    index stride = quantize?mBeat:1;
    algorithm::DataSetIdSequence seq("", 0, 0);
    mDataSet = DataSet(2);
    for(index i = 0; i <mLength; i++){
      for(index j = i + stride; j < mDM.rows(); j+=stride){
        if(mDM(i,j) < threshold || (quantize && mOnsets(i) > 0)){
          RealVector tmp{
            static_cast<double>(i),
            static_cast<double>(j)
          };
          mDataSet.add(seq.next(),tmp);
        }
      }
    }
    mTree = KDTree(mDataSet);
    mNumLinks = mTree.size();
  }

  void findLoop(){
    RealVector tmpPoint(2);
    auto query = RealVector { static_cast<double>(mStartFrame), static_cast<double>(mEndFrame)};
    auto nearest = mTree.kNearest(query, 1);
    auto nearestIds = nearest.getIds();
    if(nearestIds.size() > 0){
      mDataSet.get(nearestIds(0), mLoop);
    }

  }

  void processFrame(ComplexVectorView out, double start, double end, RealVectorView output) {
    using namespace Eigen;
    using namespace _impl;
    index startFrame = lrint(start * mSpectrogram.rows());
    index endFrame = lrint(end * mSpectrogram.rows());
    if(startFrame != mStartFrame || endFrame != mEndFrame){
      mStartFrame = startFrame;
      mEndFrame = endFrame;
      findLoop();
    }
    out = mSpectrogram.row(mPos);
    mPos = (mPos + 1) % mSpectrogram.rows();
    if(mPos >= mLoop(1))mPos = mLoop(0);
    output(0)  = mLoop(0);
    output(1)  = mLoop(1);
    output(2)  = mBeat;
    output(3)  = mNumLinks;
  }

  bool initialized(){
    return mInitialized;
  }

  index mWindowSize;
  index mHopSize;
  index mFFTSize;

private:
  index mFrameSize;
  GraphPlayUtils mUtils;
  RealVector mLoop;
  ComplexMatrix mSpectrogram;
  RealMatrix mMagnitude;
  RealMatrix mMelSpectrogram;
  Eigen::VectorXi mOnsets;
  KDTree mTree;
  MatrixXd mDM;
  DataSet mDataSet;
  bool mInitialized{false};
  int mPos{0};
  index mLength;
  index mBeat;
  index mStartFrame;
  index mEndFrame;
  double mThreshold;
  index mNumLinks;
  MedianFilter mFilter;
  PeakDetection mPD;
};
} // namespace algorithm
} // namespace fluid
