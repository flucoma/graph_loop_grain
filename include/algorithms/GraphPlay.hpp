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

#include "algorithms/public/STFT.hpp"
#include "algorithms/util/AlgorithmUtils.hpp"
#include "algorithms/util/FluidEigenMappings.hpp"
#include "algorithms/GraphPlayUtils.hpp"
#include "data/TensorTypes.hpp"
#include "data/FluidDataSet.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <random>

namespace fluid {
namespace algorithm {

class GraphPlay {

public:
  using  MatrixXd = Eigen::MatrixXd;
  using  VectorXd = Eigen::VectorXd;
  using  DataSet = FluidDataSet<std::string, double, 1>;

  void init(RealVectorView audio, index sampleRate,
            index windowSize, index fftSize, index hopSize, index numBands,
            index distance, double threshold, RealVectorView output) {
    using namespace Eigen;
    using namespace _impl;
    using namespace std;
    mWindowSize = windowSize;
    mFFTSize = fftSize;
    mHopSize = hopSize;
    mFrameSize = (mFFTSize / 2) + 1;
    mThreshold = threshold;
    STFT stft = STFT(mWindowSize, mFFTSize, mHopSize);
    mLength = std::floor((audio.size() + mHopSize) / mHopSize);
    mSpectrogram = ComplexMatrix(mLength, mFrameSize);
    stft.process(audio, mSpectrogram);
    RealMatrix magnitude(mLength, mFrameSize);
    stft.magnitude(mSpectrogram, magnitude);
    mDM = mUtils.computeDM(magnitude, numBands, sampleRate, windowSize, fftSize, distance);
    mDM.diagonal().setZero();
    mForbidden = ArrayXXd::Ones(mDM.rows(), mDM.cols());
    mVisited = ArrayXXd::Zero(mDM.rows(), mDM.cols());
    ArrayXd odf = mDM.diagonal(1).array();
    mRP = (mDM.array() < threshold).cast<double>();
    mRP = mRP.array() * mForbidden.array();
    mInitialized = true;
  }


  void processFrame(ComplexVectorView out, double start, double threshold,
    index minLength, index minDist, index forget, RealVectorView output) {
    using namespace Eigen;
    using namespace _impl;
    using namespace std;
    if(mThreshold != threshold){
      mRP =  (mDM.array() < threshold).cast<double>();
      mRP = mRP.array() * mForbidden.array();
      mThreshold = threshold;
    }
    mVisited = (mVisited.array() - 1).cwiseMax(0);
    index startFrame = lrint(start * (mSpectrogram.rows() - 1));
    if(startFrame != mStartFrame ){
      mStartFrame = startFrame;
      mPos = mStartFrame;
      index nNeighbors = mRP.col(mPos).sum();
      while(nNeighbors == 0){
        nNeighbors = mRP.col(mPos).sum();
        mPos = (mPos + 1) % mSpectrogram.rows();
      }
      mCount = 0;
    }
    else if (mCount < minLength){
      mPos = (mPos + 1) % mSpectrogram.rows();
      mCount++;
    }
    else{
        index nNeighbors = mRP.col(mPos).sum();
        index nForbidden = mForbidden.col(mPos).sum();
        index prevPos = mPos;
        if(nNeighbors > 0){
          std::vector<index> candidates(nNeighbors);
          index nCandidates = 0;
          for(index i = 0; i < mRP.rows(); i++)
            if(abs(i - mPos) > minDist && mRP(mPos, i) > 0 && mVisited(mPos, i) <= 0){
              candidates[nCandidates++] = i;
            }
            index next = mUtils.randInt(nCandidates);
            mPos = candidates[next];
            mCount = 0;
        }
        if (mPos == prevPos){
          mPos = (mPos + 1) % mSpectrogram.rows();
        }
        mVisited(prevPos, mPos) = forget;
    }
    out = mSpectrogram.row(mPos);
    output(0)  = mPos;
  }

  bool initialized(){
    return mInitialized;
  }

  index num{0};

  index mWindowSize;
  index mHopSize;
  index mFFTSize;

private:
  GraphPlayUtils mUtils;
  index mFrameSize;
  ComplexMatrix mSpectrogram;
  MatrixXd mDM;
  MatrixXd mRP;
  MatrixXd mForbidden;
  MatrixXd mVisited;
  VectorXd mDeg;
  bool mInitialized{false};
  int mPos{0};
  index mLength;
  index mStartFrame{-1};
  index mEndFrame;
  double mThreshold;
  index mCount{0};
};
} // namespace algorithm
} // namespace fluid
