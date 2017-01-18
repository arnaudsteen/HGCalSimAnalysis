#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>
#include <map>
#include <cmath>
#include "HGCalHit.h"
#include "HGCalDetId.h"
#include "HGCalCluster.h"
#include "Clustering.h"
#include "Tracking.h"
#include "Distance.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"

int main(int argc, char **argv)
{
  std::ostringstream os;
  os.str("");
  os << argv[1];
  TFile file(os.str().c_str(),"READ");
  if( file.IsOpen() )
    file.Print();
  else
    std::cout << "can not open " << os.str() << std::endl;
  TDirectory * dir = (TDirectory*)file.Get("HGCalTBAnalyzer");
  TTree *tree = (TTree*)dir->Get("HGCTB");
  if (!tree){
    std::cout << " -- Error, tree cannot be opened. Exiting..." << std::endl;
    return 0;
  }
  tree->Print();

  double xBeam, yBeam, zBeam, pBeam;
  std::vector<unsigned int> *simHitCellIdE=0;
  std::vector<double> *simHitCellEnE=0;
  tree->SetBranchAddress("xBeam",&xBeam);
  tree->SetBranchAddress("yBeam",&yBeam);
  tree->SetBranchAddress("zBeam",&zBeam);
  tree->SetBranchAddress("pBeam",&pBeam);
  tree->SetBranchAddress("simHitCellIdE",&simHitCellIdE);
  tree->SetBranchAddress("simHitCellEnE",&simHitCellEnE);

  std::map< int,std::vector<HGCalHit> > hitMap;
  HexGeom geom(false);
  const unsigned nEvts = tree->GetEntries();

  ClusteringParameters m_ClusteringParameters;
  m_ClusteringParameters.maxTransverse=1;
  m_ClusteringParameters.weightOption=std::string("logarithmic");
  float par[]={5,1};
  m_ClusteringParameters.logParams=std::vector<float>( par,par+sizeof(par)/sizeof(float));
  float pos[]={0.0, 5.35, 10.52, 14.44, 18.52, 19.67, 23.78, 25.92};
  m_ClusteringParameters.layerZPositions=std::vector<float>( pos,pos+sizeof(pos)/sizeof(float));
  Clustering clustering_algo( geom );
  clustering_algo.setParameters( m_ClusteringParameters );

  float minEnergy=25.e-6;//around 0.5 MIP
  float hgcalFullCellSide=0.6493;//cm
  int maxTransverseProfile=250;
  int nlayers=8;
  int maxEvents=2000;
  
  TrackingParameters m_TrackingParameters;
  m_TrackingParameters.doTrackCleaning=false;
  m_TrackingParameters.minTouchedLayers=4;
  m_TrackingParameters.minEnergy=minEnergy;
  m_TrackingParameters.maxEnergy=100000.0;
  Tracking tracking_algo;
  tracking_algo.setParameters( m_TrackingParameters );
  
  TFile outFile("Output.root","RECREATE");
  TDirectory * outdir = new TDirectory("hgcaltbshower","");
  TTree *outtree = new TTree("tree","HGCAL TB variables tree");
  int _evtID=0;
  float _theta;
  float _phi;
  int _nhit;
  float _x0;
  float _y0;
  float _ax;
  float _ay;
  float _dx;
  float _dy;
  std::vector<double> _transverseprofile;
  std::vector<double> _energylayer;
  outtree->Branch( "evtID",&_evtID ); 
  outtree->Branch( "theta",&_theta );
  outtree->Branch( "phi",&_phi );
  outtree->Branch( "nhit",&_nhit );
  outtree->Branch( "vx0",&_x0 );
  outtree->Branch( "vy0",&_y0 );
  outtree->Branch( "ax",&_ax );
  outtree->Branch( "ay",&_ay );
  outtree->Branch( "dx",&_dx );
  outtree->Branch( "dy",&_dy );
  outtree->Branch( "transverseprofile","std::vector<double>",&_transverseprofile );
  outtree->Branch( "energylayer","std::vector<double>",&_energylayer );

  for( unsigned ievt(0); ievt<nEvts; ++ievt ){
    if( ievt>maxEvents )
      break;
    _evtID++;
    _transverseprofile.clear();
    for( int ir=0; ir<maxTransverseProfile; ir++ )
      _transverseprofile.push_back(0);
    _energylayer.clear();
    for( int i=0; i<nlayers; i++ )
      _energylayer.push_back(0);
    
    hitMap.clear();
    _nhit=0;
    tree->GetEntry(ievt);
    if( simHitCellIdE->size()!=simHitCellEnE->size() ){
      std::cout << "problem with vector size : \t"
		<< "simHitCellIdE->size() = " << simHitCellIdE->size() << "\t"
		<< "simHitCellEnE->size() = " << simHitCellEnE->size() << std::endl;
      return 0;
    }
    for( unsigned int i=0; i<simHitCellIdE->size(); i++ ){
      if( simHitCellEnE->at(i)<minEnergy )
	continue;
      _nhit++;
      HGCalHit hit( simHitCellIdE->at(i),
		    simHitCellEnE->at(i),
		    geom );
      if( hitMap.find( hit.id().layer() )!=hitMap.end() )
	hitMap[ hit.id().layer() ].push_back(hit);
      else{
	std::vector<HGCalHit> vec;
	vec.push_back( hit );
	hitMap.insert( std::pair< int,std::vector<HGCalHit> >(hit.id().layer(),vec) );
      }
    }
    std::vector<HGCalCluster> clusters;
    for( std::map<int,std::vector<HGCalHit> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it ){
      std::vector<HGCalCluster> tmp;
      clustering_algo.run( tmp, it->second);
      clusters.insert(clusters.end(), tmp.begin(), tmp.end());
    }
    HGCalTrack track;
    tracking_algo.run(track, clusters);
    if( track.isNull()==false ){
      _theta = track.momentum().theta()*180/M_PI;
      _phi = track.momentum().phi()*180/M_PI;
      _x0 = track.vertex().x();
      _y0 = track.vertex().y();
      _ax = track.momentum().x();
      _ay = track.momentum().y();
      _dx = xBeam-track.vertex().x();
      _dy = yBeam-track.vertex().y();
      Distance<HGCalHit,HGCalTrack> dist;
      for( std::map<int,std::vector<HGCalHit> >::iterator it=hitMap.begin(); it!=hitMap.end(); ++it ){
	for( std::vector<HGCalHit>::iterator jt=it->second.begin(); jt!=it->second.end(); ++jt ){
	  HGCalHit hit=(*jt);
	  _energylayer[ hit.id().layer()-1 ] += hit.energy();
	  int ring = (int)( 10*dist.distance(hit,track)/hgcalFullCellSide );
	  if( ring<maxTransverseProfile )
	    _transverseprofile[ring]+=hit.energy();
	}
      }
    }
    outFile.cd();
    outdir->cd();
    outtree->Fill();
    std::cout << "Process " << ievt+1 << " events" << std::endl;
  }
  outdir->Write();
  outFile.Write();
  return 0;
}
