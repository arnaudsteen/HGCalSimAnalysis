#include "HGCalTrack.h"
#include <iostream>

HGCalTrack::HGCalTrack(double _chi2, double _ndof, const CLHEP::Hep3Vector &_vertex,
		       const CLHEP::Hep3Vector &_momentum)
		       //		       std::vector<HGCalHit> &_hits) 
{
  if( _chi2!=_chi2 )
    return;

  isNull_ = false;
  chi2_ = _chi2;
  ndof_ = _ndof; 
    
  vertex_ = _vertex;
  momentum_ = _momentum; 

  //  hits_ = _hits;
}
