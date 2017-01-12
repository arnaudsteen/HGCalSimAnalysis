#ifndef HEXGEOM_HH
#define HEXGEOM_HH 1

#include <iostream>
#include <vector>

class HexGeom {

public :
  HexGeom(){;}
  HexGeom(bool fine);
  virtual ~HexGeom() {}

  inline void setZPositions(std::vector<double> pos){positionZ=pos;}

  inline double layerPosition(int i) const { return (i>0&&i<positionZ.size()) ? positionZ.at(i) : -1000.;}
  std::pair<double,double> position2D(int cell) const;
  int iu(int cell) const { return cell >= 0 && cell < (int)(iuiv.size()) ? iuiv[cell].first : 10; }
  int iv(int cell) const { return cell >= 0 && cell < (int)(iuiv.size()) ? iuiv[cell].second : 10; }
  
private :
  std::vector<std::pair<double,double> > xypos;
  std::vector<std::pair<int,int> > iuiv;
  std::vector<double> positionZ;

};

#endif
