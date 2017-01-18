#include "Clustering.h"
#include <cmath>
#include <algorithm>

Clustering::Clustering(HexGeom g)
{
  geom=g;
}

void Clustering::run( std::vector<HGCalCluster> &clusters, std::vector<HGCalHit> &hits )
{
  clusters.clear();
  std::vector<HGCalHit> temp;
  int clId=0;
  for( std::vector<HGCalHit>::iterator it=hits.begin(); it!=hits.end(); ++it ){
    if( std::find(temp.begin(),temp.end(),(*it))!=temp.end() ) continue;
    temp.push_back( (*it) );
    std::vector<HGCalHit> clusterHits;
    clusterHits.push_back( (*it) );
    buildCluster( hits,temp,clusterHits );
    HGCalCluster cluster( 100*(*it).id().layer() + clId );
    cluster.setLayer( (*it).id().layer() );
    float energy=0.;
    for( std::vector<HGCalHit>::iterator jt=clusterHits.begin(); jt!=clusterHits.end(); ++jt){
      energy+=(*jt).energy();
      cluster.addHit( (*jt) );
    }
    cluster.setEnergy(energy);
    cluster.setPosition( clusterPosition(cluster) );
    clusters.push_back(cluster);
  }
}

void Clustering::buildCluster(  std::vector<HGCalHit> &hits,
				std::vector<HGCalHit> &temp,
				std::vector<HGCalHit> &clusterHits
				)
{
  HGCalHit hit=clusterHits.back();
  for( std::vector<HGCalHit>::const_iterator it=hits.begin(); it!=hits.end(); ++it){
    if( std::find(temp.begin(), temp.end(), (*it) )!=temp.end() ) continue;
    int deltaU=hit.id().iu()-(*it).id().iu();
    int deltaV=hit.id().iv()-(*it).id().iv();
    if( abs(deltaU+deltaV)>settings.maxTransverse ||
	abs(deltaU)>settings.maxTransverse ||
	abs(deltaV)>settings.maxTransverse ) continue;
    temp.push_back( (*it) );
    clusterHits.push_back( (*it) );
    buildCluster(hits, temp, clusterHits);
  }
}

CLHEP::Hep3Vector Clustering::clusterPosition(HGCalCluster cluster)
{
  double sumweight=0.;
  double x=0.;
  double y=0.;
  CLHEP::Hep3Vector xyz;
  for( std::vector<HGCalHit>::const_iterator it=cluster.hits().begin(); it!=cluster.hits().end(); ++it){
    float weight;
    if( settings.weightOption==std::string("logarithmic") )
      weight=settings.logParams[0]+settings.logParams[1]*std::log( (*it).energy()/cluster.energy() );
    else 
      weight=(*it).energy();
    if( weight<0 ) continue;
    xyz=(*it).position();
    x+=xyz.x()*weight;
    y+=xyz.y()*weight;
    sumweight+=weight;
  }
  x/=sumweight;
  y/=sumweight;
  return CLHEP::Hep3Vector( x,y,settings.layerZPositions.at( cluster.layer()-1 ) );
}
