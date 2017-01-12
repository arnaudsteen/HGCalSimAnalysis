#include "HexGeom.h"
#include <cmath>

HexGeom::HexGeom(bool fine)
{
  const int nC(15), nF(20);
  int nCoarse(11), nyCoarse(-42), nFine(15), nyFine(-56);
  int cellCoarse[nC] = {2,5,8,11,12,11,12,11,12,11,12,11,8,5,2};
  int cellFine[nF] = {3,6,9,12,15,16,15,16,15,16,15,16,15,16,15,14,11,8,5,2};
  double wafer(123.7);
  int    rows = (fine) ? nF : nC;
  double cell = (fine) ? wafer/nFine : wafer/nCoarse;
  double dx   = 0.5*cell;
  double dy   = 0.5*dx*tan(30.0*M_PI/180.0);
  int    ny   = (fine) ? nyFine : nyCoarse;
  for (int ir = 0; ir < rows; ++ir) {
    int    column = (fine) ? cellFine[ir] : cellCoarse[ir];
    int    nx     = 1 - column;
    double ypos   = dy*ny;
    for (int ic = 0; ic<column; ++ic) {
      double xpos = dx*nx;
      nx += 2;
      xypos.push_back(std::pair<double,double>(xpos/10.0,ypos/10.0));
    }
    ny += 6;
  }
  //iu iv not yet implemented for fine geom
  if( !fine ){
    int minIU[nC] = {3,1,-1,-3,-4,-4,-5,-5,-6,-6,-7,-7,-6,-5,-4};
    for(int i=0; i<nC; i++){
      int maxIU=minIU[i]+cellCoarse[i];
      for(int j=minIU[i]; j<maxIU; j++)
	iuiv.push_back( std::pair<int,int>(j,-nC/2+i) );
    }
  }

  // std::cout << "Initialize HexGeom for " << xypos.size() << " cells"
  // 	    << std::endl;
  double pos[]={0.0, 5.35, 10.52, 14.44, 18.52, 19.67, 23.78, 25.92};
  positionZ=std::vector<double>( pos,pos+sizeof(pos)/sizeof(float));
}

std::pair<double,double> HexGeom::position2D(const int cell) const{
  std::pair<double,double> xy;
  if (cell >= 0 && cell < (int)(xypos.size())) {
    xy = xypos[cell];
  } else {
    xy = std::pair<double,double>(0,0);
  }
  return xy;
}
