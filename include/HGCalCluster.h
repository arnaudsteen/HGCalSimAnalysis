#ifndef HGCALCLUSTER_HH
#define HGCALCLUSTER_HH 1

#include "HGCalHit.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "vector"

class HGCalCluster
{
 public:
  HGCalCluster(int id);
  HGCalCluster(int id, int layer, float energy );
	
  inline int   layer() const { return _layer; }
  inline float energy() const {	return _energy; }
  inline int id() const {return _id;}
  
  inline std::vector<HGCalHit> const &hits() { return _hits; }
  
  inline void setLayer(int val){ _layer=val; }
  inline void setEnergy(float val){ _energy=val; }
  inline void setPosition( CLHEP::Hep3Vector vec ){ _position=vec; }

  inline CLHEP::Hep3Vector position() const { return _position; }
  inline double x() const { return _position.x(); }
  inline double y() const { return _position.y(); }
  inline double z() const { return _position.z(); }
  
  void addHit( HGCalHit hit );
  
 protected :
  int _id;
  int _layer;
  float _energy;
  std::vector<HGCalHit> _hits;
  CLHEP::Hep3Vector _position;
};

bool operator==(HGCalCluster const& c0,HGCalCluster const& c1);

#endif
