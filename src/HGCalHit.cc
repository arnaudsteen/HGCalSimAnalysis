#include "HGCalHit.h"

HGCalHit::HGCalHit( unsigned int cellId, float energy, HexGeom g )
{
  _cellId=cellId;
  _energy=energy;
  geom=g;
  _id = HGCalDetId( _cellId,geom ); 
}
  
CLHEP::Hep3Vector HGCalHit::position() const
{
  return CLHEP::Hep3Vector( geom.position2D(_id.cellnum()).first,
			    geom.position2D(_id.cellnum()).second,
			    geom.layerPosition(_id.layer()) );
}

double HGCalHit::x() const
{
  return geom.position2D(_id.cellnum()).first;
}
double HGCalHit::y() const
{
  return geom.position2D(_id.cellnum()).second;
}
double HGCalHit::z() const
{
  return geom.layerPosition(_id.layer());
}

bool operator==(const HGCalHit& h0, const HGCalHit& h1)
{
  return h0.id()==h1.id();
}
