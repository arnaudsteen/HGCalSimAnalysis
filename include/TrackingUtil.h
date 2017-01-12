#ifndef HGCALTRACKINGUTIL_HH
#define HGCALTRACKINGUTIL_HH

#include <iostream>
#include <vector>
#include <cmath>
//#include <Eigen/Dense>

#include "CLHEP/Vector/ThreeVector.h"

#include "HGCalCluster.h"
#include "HGCalTrack.h"
#include "HGCalHit.h"

#include "Distance.h"

template <typename T>
class LeastSquare
{
 public : 
  LeastSquare(){;}
  ~LeastSquare(){;}
    
  void run(std::vector<T> &t, std::vector<float> &par, std::vector<float> &parerror)
  { 
    par.clear();
    parerror.clear();
    float xsum = 0.0;
    float ysum = 0.0;
    float zsum = 0.0;
    float zzsum = 0.0;
    float xzsum = 0.0;
    float yzsum = 0.0;
    for( typename std::vector<T>::iterator it=t.begin(); it!=t.end(); ++it ){  
      xsum = xsum + (*it).x();
      ysum = ysum + (*it).y();
      zsum = zsum + (*it).z();
      xzsum = xzsum + (*it).x()*(*it).z();
      yzsum = yzsum + (*it).y()*(*it).z();
      zzsum = zzsum + (*it).z()*(*it).z();
    }

    par.push_back( (zzsum*xsum-xzsum*zsum)/(t.size()*zzsum-zsum*zsum) );
    par.push_back( (xzsum*t.size()-xsum*zsum)/(t.size()*zzsum-zsum*zsum) );
    par.push_back( (zzsum*ysum-yzsum*zsum)/(t.size()*zzsum-zsum*zsum) );
    par.push_back( (yzsum*t.size()-ysum*zsum)/(t.size()*zzsum-zsum*zsum) );
  
    parerror.push_back( std::sqrt( zzsum/(t.size()*zzsum-zsum*zsum) ) );
    parerror.push_back( std::sqrt( t.size()/(t.size()*zzsum-zsum*zsum) ) );
    parerror.push_back( std::sqrt( zzsum/(t.size()*zzsum-zsum*zsum) ) );
    parerror.push_back( std::sqrt( t.size()/(t.size()*zzsum-zsum*zsum) ) );
  }

  float chi2(std::vector<T> &t, std::vector<float> &par)
  {
    float _chi2=0.0;
    for( typename std::vector<T>::iterator it=t.begin(); it!=t.end(); ++it ){
      CLHEP::Hep3Vector p_track(  par[0]+par[1]*(*it).z(), par[2]+par[3]*(*it).z(), (*it).z() );
      _chi2 += (p_track-(*it).position()).mag2()*12;
    }
    return _chi2/t.size();
  }
};

template <typename T>
class WeightedLeastSquare
{
 public : 
  WeightedLeastSquare(){;}
  ~WeightedLeastSquare(){;}
    
  void run(std::vector<T> &t, std::vector<float> &par, std::vector<float> &parerror)
  { 
    par.clear();
    parerror.clear();
    float xsum = 0.0;
    float ysum = 0.0;
    float zsum = 0.0;
    float zzsum = 0.0;
    float xzsum = 0.0;
    float yzsum = 0.0;
    float esum = 0.0;
    for( typename std::vector<T>::iterator it=t.begin(); it!=t.end(); ++it ){
      xsum = xsum + (*it).x()*(*it).energy();
      ysum = ysum + (*it).y()*(*it).energy();
      zsum = zsum + (*it).z()*(*it).energy();
      xzsum = xzsum + (*it).x()*(*it).z()*(*it).energy();
      yzsum = yzsum + (*it).y()*(*it).z()*(*it).energy();
      zzsum = zzsum + (*it).z()*(*it).z()*(*it).energy();
      esum = esum + (*it).energy();
    }

    par.push_back( (zzsum*xsum-xzsum*zsum)/(esum*zzsum-zsum*zsum) );
    par.push_back( (xzsum*esum-xsum*zsum)/(esum*zzsum-zsum*zsum) );
    par.push_back( (zzsum*ysum-yzsum*zsum)/(esum*zzsum-zsum*zsum) );
    par.push_back( (yzsum*esum-ysum*zsum)/(esum*zzsum-zsum*zsum) );
  
    parerror.push_back( std::sqrt( zzsum/(esum*zzsum-zsum*zsum) ) );
    parerror.push_back( std::sqrt( esum/(esum*zzsum-zsum*zsum) ) );
    parerror.push_back( std::sqrt( zzsum/(esum*zzsum-zsum*zsum) ) );
    parerror.push_back( std::sqrt( esum/(esum*zzsum-zsum*zsum) ) );
  }

  float chi2(std::vector<T> &t, std::vector<float> &par)
  {
    float _chi2=0.0;
    for( typename std::vector<T>::iterator it=t.begin(); it!=t.end(); ++it ){
      CLHEP::Hep3Vector p_track(  par[0]+par[1]*(*it).z(), par[2]+par[3]*(*it).z(), (*it).z() );
      _chi2 += (p_track-(*it).position()).mag2()*12;
    }
    return _chi2;
  }
};

template <typename T>
class TrackCleaner
{
 public:
  TrackCleaner(){;}
  ~TrackCleaner(){;}
  void clean( std::vector<T> &col, std::vector<T> &cleancol, HGCalTrack &track, double maxDistance )
  {
    Distance<T,HGCalTrack> dist;
    for( typename std::vector<T>::iterator it=col.begin(); it!=col.end(); ++it )
      if( dist.distance( (*it),track )<maxDistance )
	cleancol.push_back(*it);
  }
};

//template <typename T>
//class PrincipalComponentAnalysis
//{
// public : 
//  PrincipalComponentAnalysis(){;}
//  ~PrincipalComponentAnalysis(){;}
//      
//  Eigen::MatrixXd covarianceMatrix() const {return _covarianceMatrix;};
//  Eigen::MatrixXd eigenVectors() const {return _eigenVectors;}
//  std::vector<double> &eigenValues() const {return _eigenValues; }
// private :
//  Eigen::MatrixXd _covarianceMatrix;
//  Eigen::MatrixXd _eigenVectors;
//  std::vector<double> _eigenValues;
// public:
//  void runPCA(T t)
//  {
//    covarianceMatrix = Eigen::MatrixXd(t.size(),3);
//    float x,y,z; x=y=z=0.0;
//    for( auto p : t ){
//      x+=p.x();
//      y+=p.y();
//      z+=p.z();
//    }
//    math::XYZPoint mean=math::XYZPoint(x,y,z)/t.size();
//    int count=0;
//    for( auto p : t ){
//      covarianceMatrix(count,0) = p.x()-mean.x();
//      covarianceMatrix(count,1) = p.y()-mean.y();
//      covarianceMatrix(count,2) = p.z()-mean.z();
//      count++;
//    }    
//    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver( _covarianceMatrix.transpose()*_covarianceMatrix );
//    _eigenVectors=eigensolver.eigenvectors();
//    for(unsigned int i=0; i<eigensolver.eigenvalues().size(); i++){
//      _eigenValues.push_back(eigensolver.eigenvalues()[i]);
//    }
//  }
//};

#endif
