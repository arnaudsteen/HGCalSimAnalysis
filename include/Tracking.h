#ifndef TRACKING_HH 
#define TRACKING_HH 1

#include "HGCalTrack.h" 
#include "HGCalCluster.h"
#include "CLHEP/Vector/ThreeVector.h"

struct TrackingParameters{
  bool doTrackCleaning;
  int maxTransverse;
  int minTouchedLayers;
  float minEnergy;
  float maxEnergy;
TrackingParameters() :   doTrackCleaning(false),
    maxTransverse(13),
    minTouchedLayers(4),
    minEnergy(25.0e-6),
    maxEnergy(10000)
  {;}
}; 

class Tracking
{
 public:
  Tracking(){;}
  ~Tracking(){;}
  void run( HGCalTrack &track, std::vector<HGCalCluster> &clusters );  
  inline void setParameters( TrackingParameters par ){settings=par;}
  
 private:
  TrackingParameters settings;
};

#endif
