#pragma once

#include "algorithms/util/PeakDetection.hpp"
#include "algorithms/public/DataSetIdSequence.hpp"
#include "algorithms/util/DistanceFuncs.hpp"
#include "algorithms/public/KDTree.hpp"
#include "algorithms/public/KMeans.hpp"
#include "algorithms/public/MelBands.hpp"
#include "algorithms/util/AlgorithmUtils.hpp"
#include "algorithms/util/FluidEigenMappings.hpp"
#include "algorithms/util/MedianFilter.hpp"
#include "algorithms/util/SpectralEmbedding.hpp"
#include "data/TensorTypes.hpp"
#include "data/FluidDataSet.hpp"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <random>

namespace fluid {
namespace algorithm {

class GraphPlayUtils {

public:
  using  MatrixXd = Eigen::MatrixXd;
  using  VectorXd = Eigen::VectorXd;
  using DataSet = FluidDataSet<std::string, double, 1>;

  GraphPlayUtils(){
    using namespace std;
    random_device rd;
    mGen = mt19937(rd());
    mDis = uniform_real_distribution<> (0.0, 1.0);
    mFilter.init(5);
  }

  index randInt(index N){
    return static_cast<index>(mDis(mGen) * N);
  }

  double rand(){
    return mDis(mGen);
  }

  Eigen::ArrayXXd computeDM(RealMatrixView mag, index numBands,
    double sampleRate, index windowSize, index fftSize, index dist){
    using namespace Eigen;
    using namespace _impl;
    MelBands melBands = MelBands(numBands, fftSize);
    melBands.init(20, 5000, numBands, mag.cols(), sampleRate, windowSize);
    RealMatrix melSpec = RealMatrix(mag.rows(), numBands);
    for(index i = 0; i < mag.rows(); i++){
      melBands.processFrame(mag.row(i), melSpec.row(i), true, false, false);
    }
    MatrixXd tmp = asEigen<Matrix>(melSpec);
    return DistanceMatrix(tmp, dist);
  }

  void onsetDetection(Eigen::Ref<Eigen::ArrayXd> odf,
                      Eigen::Ref<Eigen::ArrayXXd> transitions, index offset = 2){
    for(index i = 0; i < odf.size(); i++){
      odf(i) = odf(i) - mFilter.processSample(odf(i));
    }
    auto onsets = mPD.process(odf, 0, 0.1, false, false);
    for(index i = 0; i < onsets.size(); i++){
      index pos = onsets[i].first;
      index start = std::max(index(0), pos - offset);
      index end = std::min(transitions.rows() - 1, pos + offset + 1);
      transitions.block(start, 0, end - start, transitions.cols()).setZero();
      transitions.block(0, start, transitions.rows(), end - start).setZero();
    }
  }

  FluidTensor<index, 1> kmeans(RealMatrixView data, index nClusters){
      algorithm::DataSetIdSequence seq("", 0, 0);
      index minClusterSize = 10;
      DataSet tmpDS = DataSet(data.cols());
      for(index i = 0; i < data.rows(); i++){
        tmpDS.add(seq.next(), data.row(i));
      }
      mKMeans.clear();
      mKMeans.train(tmpDS, nClusters, 100);
      FluidTensor<index, 1> clusters(data.rows());
      mKMeans.getAssignments(clusters);
      RealMatrix means(mKMeans.size(), mKMeans.dims());
      mKMeans.getMeans(means);
      for(index i = 0; i < nClusters; i++){
        index cSize = mKMeans.getClusterSize(i);
        //std::cout<<cSize<<std::endl;
        if(cSize > 0 && cSize < minClusterSize){
          RealMatrix distances(1, mKMeans.size());
          RealMatrix mean(1, mKMeans.dims());
          mean.row(0) =  means.col(i);
          mKMeans.getDistances(mean, distances);
          index closest = std::min_element(distances.begin(), distances.end()) - distances.begin();
          for(index j = 0; j < clusters.size(); j++)
            if(clusters(j) == i) clusters(j) = closest;
        }
      }
      //std::ofstream ofs ("clusters.mat", std::ofstream::out);ofs << clusters;ofs.close();
      return clusters;
  }

  void writeMatrix(Eigen::Ref<Eigen::MatrixXd> mat, std::string name){
    //std::ofstream ofs (name+".mat", std::ofstream::out);ofs << mat;ofs.close();
  }

  FluidTensor<index, 1>  spectralClustering(Eigen::Ref<Eigen::ArrayXXd> dm, index numClusters = 0){
    using namespace Eigen;
    index maxClusters = numClusters > 0? numClusters : std::min(index(50), dm.rows());
    index nPoints = dm.rows();
    MatrixXd tmpRp = (dm.array() < 0.25).cast<double>();
    MatrixXd weightedGraph = (1 - dm) * tmpRp.array();
    SpectralEmbedding spectralEmbedding;
    spectralEmbedding.train(weightedGraph.sparseView(), maxClusters);
    if(numClusters == 0 ){
      VectorXd eigenValues  =  spectralEmbedding.eigenValues();
      ArrayXd diff = (
        eigenValues.segment(1, eigenValues.size() - 2) -
        eigenValues.segment(0, eigenValues.size() - 2));
        VectorXd::Index maxIndex;
        double maxVal = diff.maxCoeff(&maxIndex);
        numClusters = std::min(2*(maxIndex + 1), maxClusters);
    }
    MatrixXd eigenVectors = spectralEmbedding.eigenVectors().block(0,0, nPoints, numClusters);
    eigenVectors.rowwise().normalize();
    return kmeans( _impl::asFluid(eigenVectors), numClusters);
  }


private:
  MedianFilter mFilter;
  PeakDetection mPD;
  KMeans mKMeans;
  std::mt19937 mGen;
  std::uniform_real_distribution<> mDis;
};
} // namespace algorithm
} // namespace fluid
