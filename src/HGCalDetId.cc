#include "HGCalDetId.h"

HGCalDetId::HGCalDetId( unsigned int cellId,HexGeom geom ) : _cellId( cellId )
{
  _iu = geom.iu( cellnum() );
  _iv = geom.iv( cellnum() );
}


bool operator==(HGCalDetId const& id0, HGCalDetId const& id1)
{
  return id0.rawId()==id1.rawId();
}
