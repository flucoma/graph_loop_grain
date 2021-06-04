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

#include "algorithms/GraphPlayUtils.hpp"
#include "algorithms/util/RTPGHI.hpp"
#include "algorithms/public/STFT.hpp"
#include "algorithms/util/AlgorithmUtils.hpp"
#include "algorithms/util/FluidEigenMappings.hpp"
#include "data/FluidDataSet.hpp"
#include "data/TensorTypes.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <fstream>
#include <random>
#include <vector>

namespace fluid {
namespace algorithm {

class GraphGrain {

public:
  using MatrixXd = Eigen::MatrixXd;
  using VectorXd = Eigen::VectorXd;
  using DataSet = FluidDataSet<std::string, double, 1>;

  void init(RealVectorView audio, index sampleRate, index windowSize,
            index fftSize, index hopSize, index numBands, index distance,
            double threshold, index nClusters, RealVectorView output) {
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
    mMagnitude = RealMatrix(mLength, mFrameSize);
    stft.magnitude(mSpectrogram, mMagnitude);
    mDM = mUtils.computeDM(mMagnitude, numBands, sampleRate, windowSize,
                           fftSize, distance);
    mDM.diagonal().setZero();
    mForbidden = ArrayXXd::Ones(mDM.rows(), mDM.cols());
    mVisited = ArrayXXd::Zero(mDM.rows(), mDM.cols());
    ArrayXd odf = mDM.diagonal(1).array();
    mClusters = FluidTensor<index, 1>(mLength);
    mUtils.onsetDetection(odf, mForbidden);
    if (nClusters != 1) {
      mClusters = mUtils.spectralClustering(mDM, nClusters);
      for (index i = 0; i < mLength; i++) {
        for (index j = 0; j < mLength; j++) {
          if (mClusters(i) != mClusters(j)) {
            mForbidden(i, j) = 0;
            mForbidden(j, i) = 0;
          }
        }
      }
    }
    mRP = (mDM.array() < threshold).cast<double>();
    mRP = mRP.array() * mForbidden.array();
    mRTPGHI.init(fftSize);
    mPrevMag = RealVector(mFrameSize);
    mPrevMag = mMagnitude.row(0);
    mInitialized = true;
  }

  index selectProb() {
    index nNeighbors = mRP.col(mPos).sum();
    if (nNeighbors == 0)
      return nextInCluster(mPos);
    std::vector<index> candidates(nNeighbors);
    std::vector<double> acumProbs(nNeighbors);
    index nCandidates = 0;
    double prob = 0;
    for (index i = 0; i < mRP.rows(); i++) {
      if (mRP(mPos, i) > 0 && mVisited(mPos, i) <= 0) {
        prob = prob + (1 - mDM(mPos, i));
        acumProbs[nCandidates] = prob;
        candidates[nCandidates++] = i;
      }
    }
    double rnd = prob * mUtils.rand();
    index selected = 0;
    while (acumProbs[selected] < rnd)
      selected++;
    return candidates[selected];
  }

  index nextInCluster(index current) {
    for (index i = current + 1; i < mLength + current; i++) {
      index pos = i % mLength;
      if (mClusters(pos) == mClusters(current))
        return pos;
    }
  }

  index selectRand(double randomness) {
    index nNeighbors = mRP.col(mPos).sum();
    if (nNeighbors == 0) {
      mVisited.row(mPos).setZero();
      return nextInCluster(mPos);
    }
    std::vector<index> candidates(nNeighbors);
    index nCandidates = 0;
    for (index i = 0; i < mRP.rows(); i++)
      if (i != mPos && mRP(mPos, i) > 0 && mVisited(mPos, i) <= 0) {
        candidates[nCandidates++] = i;
      }
    std::sort(candidates.begin(), candidates.end(),
              [this](index a, index b) { return mDM(mPos, a) < mDM(mPos, b); });
    if (randomness == 0)
      return candidates[0];
    index k = lrint(randomness * nCandidates);
    index next = mUtils.randInt(k);
    return candidates[next];
  }

  index selectNearest() {
    index nNeighbors = mRP.col(mPos).sum();
    if (nNeighbors == 0) {
      return nextInCluster(mPos);
    }
    double minDist = infinity;
    index selected = 0;
    for (index i = 0; i < mRP.rows(); i++) {
      if ((i != mPos) && (mRP(mPos, i) > 0) && (mDM(mPos, i) < minDist) &&
          (mVisited(mPos, i) == 0)) {
        selected = i;
        minDist = mDM(mPos, i);
      }
    }
    if (selected == 0)
      return nextInCluster(mPos);
    return selected;
  }

  void processFrame(ComplexVectorView out, double start, double threshold,
                    index forget, double rand, index phaseGen,
                    RealVectorView output) {
    using namespace Eigen;
    using namespace _impl;

    if (mThreshold != threshold) {
      mRP = (mDM.array() < threshold).cast<double>();
      mRP = mRP.array() * mForbidden.array();
      mThreshold = threshold;
    }

    index next = (mPos + 1) % mSpectrogram.rows();
    mVisited = (mVisited.array() - 1).cwiseMax(0);
    index startFrame = lrint(start * (mSpectrogram.rows() - 1));
    if (startFrame != mStartFrame) {
      mStartFrame = startFrame;
      mPos = mStartFrame;
      index nNeighbors = mRP.col(mPos).sum();
      while (nNeighbors == 0) {
        nNeighbors = mRP.col(mPos).sum();
        mPos = (mPos + 1) % mSpectrogram.rows();
      }
      mCount = 0;
    } else {
      index prevPos = mPos;
      mPos = selectRand(rand);
      mVisited(prevPos, mPos) = forget;
    }
    RealVectorView frame = mMagnitude.row(mPos);
    if (phaseGen > 0)
      mRTPGHI.processFrame(frame, out, mWindowSize, mFFTSize, mHopSize, 1e-5);
    else
      out = mSpectrogram.row(mPos);
    output(0) = mPos;
    output(1) = mClusters(mPos);
  }

  bool initialized() { return mInitialized; }

  index mWindowSize;
  index mHopSize;
  index mFFTSize;

private:
  GraphPlayUtils mUtils;
  ComplexMatrix mSpectrogram;
  RealMatrix mMagnitude;
  RealVector mPrevMag;
  index mFrameSize;
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
  RTPGHI mRTPGHI;
  FluidTensor<index, 1> mClusters;
  double mPrevGain{0};
};
} // namespace algorithm
} // namespace fluid
