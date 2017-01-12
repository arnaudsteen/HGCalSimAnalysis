#ifndef HGCALCALOTRACK_HH
#define HGCALCALOTRACK_HH 1

#include "HGCalHit.h"
#include "CLHEP/Vector/ThreeVector.h"

class HGCalTrack 
{
 public:
  HGCalTrack( ){ isNull_=true;}
  HGCalTrack( double _chi2, double _ndof, const CLHEP::Hep3Vector &_vertex,
	      const CLHEP::Hep3Vector &_momentum );
  //std::vector<HGCalHit> &_hits);
      
  bool isNull() const { return isNull_; }
  double chi2() const { return chi2_; }
  double ndof() const { return ndof_; }
  double normalisedChi2() const { return chi2_/ndof_; }

  /* HGCalHit hit(int i) const { return hits_.at(i); } */
  /* std::vector<HGCalHit> &hits() { return hits_; } */

  CLHEP::Hep3Vector vertex() const { return vertex_; }
  CLHEP::Hep3Vector momentum() const { return momentum_; }
  //CovarianceMatrix covarianceMatrix() const { return cov_; }; 

  CLHEP::Hep3Vector expectedTrackProjection(float z) const { return CLHEP::Hep3Vector( vertex_.x()+momentum_.x()*z,vertex_.y()+momentum_.y()*z,z); }

 private:
  bool isNull_;
  double chi2_;
  double ndof_; 
      
  CLHEP::Hep3Vector vertex_;
  CLHEP::Hep3Vector momentum_; 
  //std::vector<HGCalHit> hits_;
};

#endif
