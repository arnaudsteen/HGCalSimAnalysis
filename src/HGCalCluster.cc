#include "HGCalCluster.h"

HGCalCluster::HGCalCluster(int id) : _id(id)
{
  _layer=0;
  _energy=0.;
  _hits.clear();
}

HGCalCluster::HGCalCluster(int id, int layer, float energy ) : _id(id),
							       _layer(layer),
							       _energy(energy)
{
  _hits.clear();
}

void HGCalCluster::addHit( HGCalHit hit )
{
  _hits.push_back(hit);
}

bool operator==(HGCalCluster const& c0,HGCalCluster const& c1)
{
  return c0.id()==c1.id();
}
