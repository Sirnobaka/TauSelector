
#include <memory>
#include <string>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
//
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
//
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"


#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include <cstdio>
#include <iostream>
#include <math.h>

//
// class decleration
//

class MySimpleParticle {
  public:
    MySimpleParticle() { _px=_py=_pz=_e=0; _pdgid=0; }
    MySimpleParticle(double x, double y, double z, double e, int id) { _px=x; _py=y; _pz=z; _e=e; _pdgid=id; }
  
    double px()    { return _px;    }
    double py()    { return _py;    }
    double pz()    { return _pz;    }
    double e ()    { return _e;     }
    double pt()    { return sqrt(px()*px()+py()*py()); }
    int    pdgid() { return _pdgid; }
  
    void setPx(double x ) { _px = x;     }
    void setPy(double y ) { _py = y;     }
    void setPz(double z ) { _pz = z;     }
    void setE (double e ) { _e  = e;     }
    void setPdgid(int id) { _pdgid = id; }

    double getAnglePhi();
    double getAngleEta();
 
  private:
    double _px,_py,_pz,_e;
    int    _pdgid;
};

inline double MySimpleParticle::getAnglePhi() {
  // conventions as in ANGFI(X,Y) of tauola.f and PHOAN1 of photos.f
  // but now 'phi' in name define that it is rotation in px py

  double buf = 0.;
  
  if(fabs(py())<fabs(px())) {
    buf = atan( fabs(py()/px()) );
    if(px()<0.) buf = M_PI-buf;
  }
  else buf = acos( px()/sqrt(px()*px()+py()*py()) );
   
  if(py()<0.)     buf = -buf;
 
  return buf;
}

inline double MySimpleParticle::getAngleEta()
 {
   
  double buf = 0.;
/*
  if(fabs(px())<fabs(pz())) {
    buf = atan( fabs(px()/pz()) );
    if(pz()<0.) buf = M_PI-buf;
  }
  else buf = acos( pz()/sqrt(pz()*pz()+px()*px()) );
*/

  buf = 0.5*log( (sqrt(px()*px()+py()*py()+pz()*pz()) + pz())/(sqrt(px()*px()+py()*py()+pz()*pz()) - pz()) );

  return buf;
}
