#ifndef CLUSTERING_HH 
#define CLUSTERING_HH 1

#include "HGCalCluster.h" 
#include "HGCalHit.h"
#include "HexGeom.h"
#include "CLHEP/Vector/ThreeVector.h"
#include <set>
#include <string>

struct ClusteringParameters{
  int maxTransverse;
  std::string weightOption;
  std::vector<float> logParams;
  std::vector<float> layerZPositions;
ClusteringParameters() : maxTransverse(1) ,
    weightOption("linear")
  {
    static const float par[]={5,1};
    logParams=std::vector<float>( par,par+sizeof(par)/sizeof(float));
    static const float pos[]={0.0, 5.35, 10.52, 14.44, 18.52, 19.67, 23.78, 25.92};
    layerZPositions=std::vector<float>( pos,pos+sizeof(pos)/sizeof(float));
  }
}; 

class Clustering
{
 public:
  Clustering(HexGeom g);
  ~Clustering(){;}
  void run( std::vector<HGCalCluster> &clusters, std::vector<HGCalHit> &hits );  
  inline void setParameters( ClusteringParameters par ){settings=par;}
  
 private:
  void buildCluster( std::vector<HGCalHit> &hits,
		     std::vector<HGCalHit> &temp,
		     std::vector<HGCalHit> &clusterHits);
  CLHEP::Hep3Vector clusterPosition(HGCalCluster cluster);
  
  ClusteringParameters settings;
  HexGeom geom;
};

#endif
