#include "Tracking.h"
#include "TrackingUtil.h"

#include <iostream>
#include <set>
#include <cmath>

void Tracking::run( HGCalTrack &track, std::vector<HGCalCluster> &clusters )
{
  std::vector<HGCalCluster> tmp;
  std::vector<HGCalHit> hits;
  std::set<int> touchedLayers;
  for( auto cluster : clusters ){
    if( cluster.energy()<settings.minEnergy || cluster.energy()>settings.maxEnergy ) 
      continue;
    tmp.push_back(cluster);
    touchedLayers.insert( cluster.layer() );
    for( auto hit : cluster.hits() )
      hits.push_back(hit);
  }
  if( (int)touchedLayers.size()>settings.minTouchedLayers ){
    WeightedLeastSquare<HGCalCluster> wls;
    std::vector<float> trackPar;
    std::vector<float> trackParError;
    wls.run( tmp, trackPar, trackParError );
    float chi2 = wls.chi2( tmp, trackPar );
    int ndof = tmp.size();
    CLHEP::Hep3Vector momentum = CLHEP::Hep3Vector(-1., 0., trackPar[1]).cross( CLHEP::Hep3Vector(0., -1., trackPar[3]) );
    CLHEP::Hep3Vector vertex = CLHEP::Hep3Vector( trackPar[0], trackPar[2], 0.0 );
    track=HGCalTrack( chi2, ndof, vertex, momentum);//, hits);
    if( settings.doTrackCleaning ){
      std::vector<HGCalCluster> cleancol;
      TrackCleaner<HGCalCluster> cleaner;
      cleaner.clean( tmp, cleancol, track, settings.maxTransverse );
      touchedLayers.clear();
      hits.clear();
      for( auto cluster : cleancol ){
      	touchedLayers.insert( cluster.layer() );
      	for( auto hit : cluster.hits() )
      	  hits.push_back(hit);
      }
      if( (int)touchedLayers.size()>settings.minTouchedLayers ){
	wls.run( cleancol, trackPar, trackParError);
	chi2 = wls.chi2( cleancol, trackPar);
	ndof = cleancol.size();
	momentum = CLHEP::Hep3Vector(-1., 0., trackPar[1]).cross( CLHEP::Hep3Vector(0., -1., trackPar[3]) );
	vertex = CLHEP::Hep3Vector( trackPar[0], trackPar[2], 0.0 );
	track = HGCalTrack( chi2, ndof, vertex, momentum);//, hits);
      }
    }
  }
}
