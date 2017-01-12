#ifndef HGCALHIT_HH
#define HGCALHIT_HH 1

#include "CLHEP/Vector/ThreeVector.h"
#include "HGCalDetId.h"

#include <iostream>

class HGCalHit
{
 public:
  HGCalHit(){;}
  HGCalHit( unsigned int cellId, float energy, HexGeom g );
  ~HGCalHit(){;}

  CLHEP::Hep3Vector position() const;
  double x() const;
  double y() const;
  double z() const;
  
  inline HGCalDetId id() const {return _id;}
  inline float energy() const {return _energy;}

 private:
  HGCalDetId _id;
  unsigned int _cellId;
  float _energy;
  HexGeom geom;
};

bool operator==(const HGCalHit& h0, const HGCalHit& h1);

#endif
