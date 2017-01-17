#ifndef DISTANCE_HH
#define DISTANCE_HH

#include <iostream>
#include <cmath>

#include "HGCalCluster.h"
#include "HGCalHit.h"
#include "HGCalTrack.h"

#include "CLHEP/Vector/ThreeVector.h"

template <typename T, typename S>
  class Distance
{
 public : 
  Distance(){;}
  ~Distance(){;}
    
  float distance(T t, S s){ return (t.position()-s.position()).mag(); }
};

template <typename T>
class Distance<T, CLHEP::Hep3Vector >
{
 public : 
  Distance(){;}
  ~Distance(){;}
    
  float distance(T t, CLHEP::Hep3Vector p){ return (t.position()-p).mag(); }
};

template <typename S> 
class Distance<S, HGCalTrack>
{ 
 public :  
  Distance(){;} 
  ~Distance(){;} 
      
  float distance(S s, HGCalTrack t)
  {
    /*
      point H(x,y,z) (S s)
      track T, orientation vector u(tx,ty,tz)
      d(H,T) = || vec(BH) * u || / || u || where B is a point from the track
    */

    CLHEP::Hep3Vector H=s.position();
    
    //B : a track point 
    CLHEP::Hep3Vector B=t.vertex();
    CLHEP::Hep3Vector u=t.momentum();
    CLHEP::Hep3Vector v=B-H;
    return (v.cross(u)).mag()/u.mag();
  }
  
  float distanceInLayer(S s,HGCalTrack t)
  {
    CLHEP::Hep3Vector impact=t.expectedTrackProjection(s.z());
    return (impact-s.position()).mag();
  }
}; 

class DistanceBetweenTrackAndPoint
{ 
 public :  
  DistanceBetweenTrackAndPoint(){;} 
  ~DistanceBetweenTrackAndPoint(){;} 
      
  float distance(HGCalTrack track, CLHEP::Hep3Vector xyz)
  {
    CLHEP::Hep3Vector B=track.vertex();

    CLHEP::Hep3Vector u=track.momentum();
    CLHEP::Hep3Vector v=B-xyz;
	
    return (v.cross(u)).mag()/u.mag();
  }

  float distanceInLayer(HGCalTrack track, CLHEP::Hep3Vector xyz)
  {
    CLHEP::Hep3Vector impact=track.expectedTrackProjection(xyz.z());
    return (impact-xyz).mag();
  }
}; 
#endif
