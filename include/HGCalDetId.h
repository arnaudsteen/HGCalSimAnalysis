#ifndef HGCALDETID_HH
#define HGCALDETID_HH 1

#include "HexGeom.h"

class HGCalDetId
{
 public:
  HGCalDetId( unsigned int cellId,HexGeom geom );
  HGCalDetId(){;}
  ~HGCalDetId(){;}

  static const int kHGCalCellHOffset      = 0;
  static const int kHGCalCellHMask        = 0xFF;
  static const int kHGCalCellTypHOffset   = 8;
  static const int kHGCalCellTypHMask     = 0x1;
  static const int kHGCalWaferHOffset     = 9;
  static const int kHGCalWaferHMask       = 0x3FF;
  static const int kHGCalLayerHOffset     = 19;
  static const int kHGCalLayerHMask       = 0x7F;
  static const int kHGCalZsideHOffset     = 26;
  static const int kHGCalZsideHMask       = 0x1;
  static const int kHGCalSubdetHOffset    = 27;
  static const int kHGCalSubdetHMask      = 0x7;

  inline unsigned int rawId() const { return _cellId; }
  inline int layer() const { return (_cellId>>kHGCalLayerHOffset)&kHGCalLayerHMask; }
  inline int cellnum() const { return (_cellId>>kHGCalCellHOffset)&kHGCalCellHMask; }
  inline int celltype() const{ return (_cellId>>kHGCalCellTypHOffset)&kHGCalCellTypHMask; }
  inline int iu() const {return _iu;}
  inline int iv() const {return _iv;}
 private:
  unsigned int _cellId;
  int _iu;
  int _iv;
};

bool operator==(HGCalDetId const& id0, HGCalDetId const& id1);

#endif
