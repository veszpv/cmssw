#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"

#include <climits>
#include <iostream>

class Phase1GeometryTopologyAnalyzer : public edm::one::EDAnalyzer<> {
public:
  explicit Phase1GeometryTopologyAnalyzer( const edm::ParameterSet& ) {};
  ~Phase1GeometryTopologyAnalyzer() {};
  
  void beginJob() override {}
  void analyze(edm::Event const& iEvent, edm::EventSetup const&) override;
  void endJob() override {}
  std::string print(const TrackerTopology* tTopo, DetId id) const;
  std::string print(const TrackerTopology* tTopo, const TrackerGeometry* tkGeom, DetId id) const;
  std::string print(const GlobalVector& gvx, const GlobalVector& gvy, const GlobalVector& gvz, const GlobalPoint& gp) const;
};

std::string Phase1GeometryTopologyAnalyzer::print(const TrackerTopology* tTopo, DetId id) const {
  return tTopo->print(id);
  /*
  uint32_t subdet=id.subdetId();
  std::stringstream strstr;

  if ( subdet == PixelSubdetector::PixelBarrel ) {
    unsigned int theLayer  = tTopo->pxbLayer(id);
    unsigned int theLadder = tTopo->pxbLadder(id);
    unsigned int theModule = tTopo->pxbModule(id);
    strstr << "PixelBarrel" 
           << " Layer " << theLayer
           << " Ladder " << theLadder
           << " Module " << theModule ;
    strstr << " (" << id.rawId() << ")";
    return strstr.str();
  }

  if ( subdet == PixelSubdetector::PixelEndcap ) {
    unsigned int theSide   = tTopo->pxfSide(id);
    unsigned int theDisk   = tTopo->pxfDisk(id);
    unsigned int theBlade  = tTopo->pxfBlade(id);
    unsigned int thePanel  = tTopo->pxfPanel(id);
    unsigned int theModule = tTopo->pxfModule(id);
    std::string side  = (tTopo->pxfSide(id) == 1 ) ? "-" : "+";
    strstr << "PixelEndcap" 
           << " Side   " << theSide << side
           << " Disk   " << theDisk
           << " Blade  " << theBlade
           << " Panel  " << thePanel
           << " Module " << theModule ;
    strstr << " (" << id.rawId() << ")";
    return strstr.str();
  }
  
  return strstr.str();*/
}


std::string Phase1GeometryTopologyAnalyzer::print(const TrackerTopology* tTopo, const TrackerGeometry* tkGeom, DetId id) const {

  uint32_t subdet=id.subdetId();
  if ( subdet != PixelSubdetector::PixelBarrel && subdet != PixelSubdetector::PixelEndcap ) return "";

  assert((const void*) tkGeom->idToDetUnit(id)==(const void*) tkGeom->idToDet(id));

  std::stringstream strstr;

  const Surface& surface = tkGeom->idToDet(id)->surface();
  // Center of surface
  LocalPoint lp(0., 0., 0.);
  GlobalPoint gp = surface.toGlobal(lp);

  // Rotation of surface -- not really working in FreeCAD
  GlobalVector rot(surface.rotation().zy()-surface.rotation().yz(),
		   surface.rotation().zx()-surface.rotation().xz(),
		   surface.rotation().yx()-surface.rotation().xy());
  float angle= (acos((surface.rotation().xx()+surface.rotation().yy()+surface.rotation().zz()-1)/2)/atan(1))*45;
  
  // Printout for CMSSW checks
  if ( subdet == PixelSubdetector::PixelBarrel ) {
    // width is short edge along local x, length is long edge along local y
    strstr << "Dimensions "<< "(" << ", "<< surface.bounds().width() << surface.bounds().length() << ", "<< surface.bounds().thickness() << "); ";
    strstr << "Center "<< "(" << gp.x() << ", "<< gp.y() << ", "<< gp.z() << "); ";
    LocalVector lvx(1., 0., 0.);
    GlobalVector gvx = surface.toGlobal(lvx);
    float R_x_gvx=gvx.x()*gp.x()+gvx.y()*gp.y()+gvx.z()*gp.z();
    strstr << "R="<< sqrt(gp.x()*gp.x()+gp.y()*gp.y()) << "; phi=" << (atan2(gp.y(), gp.x())/atan(1))*45 << "; R_x_gvx =" << ((R_x_gvx < 1e-6) ? 0 : R_x_gvx);
  } 
  else if ( subdet == PixelSubdetector::PixelEndcap ) {
    // width is short edge along local x, length is long edge along local y
    strstr << "Dimensions "<< "(" << ", "<< surface.bounds().width() << surface.bounds().length() << ", "<< surface.bounds().thickness() << "); ";
    strstr << "Center "<< "(" << gp.x() << ", "<< gp.y() << ", "<< gp.z() << "); ";
    LocalVector lvx(1., 0., 0.);
    GlobalVector gvx = surface.toGlobal(lvx);
    float R_x_gvx=gvx.x()*gp.x()+gvx.y()*gp.y()+gvx.z()*gp.z();
    strstr << "R="<< sqrt(gp.x()*gp.x()+gp.y()*gp.y()) << "; phi=" << (atan2(gp.y(), gp.x())/atan(1))*45 << "; R_x_gvx =" << ((R_x_gvx < 1e-6) ? 0 : R_x_gvx);
  }
  
  // Printout for FreeCAD import
//   if ( subdet == PixelSubdetector::PixelEndcap && tTopo->pxfDisk(id)==1 && tTopo->pxfSide(id)==2 && 
//        ((tTopo->pxfBlade(id)>=1&&tTopo->pxfBlade(id)<=5) || (tTopo->pxfBlade(id)>=17&&tTopo->pxfBlade(id)<=22) || (tTopo->pxfBlade(id)>=23&&tTopo->pxfBlade(id)<=30) || (tTopo->pxfBlade(id)>=48&&tTopo->pxfBlade(id)<=56))
//      ) {
  if ( subdet == PixelSubdetector::PixelEndcap && tTopo->pxfDisk(id)==1 && tTopo->pxfSide(id)==2 && 
       ((tTopo->pxfBlade(id)>=1&&tTopo->pxfBlade(id)<=5) || (tTopo->pxfBlade(id)>=17&&tTopo->pxfBlade(id)<=22) || (tTopo->pxfBlade(id)>=23&&tTopo->pxfBlade(id)<=30) || (tTopo->pxfBlade(id)>=48&&tTopo->pxfBlade(id)<=56))
       ) {
//  if ( subdet == PixelSubdetector::PixelBarrel && tTopo->pxbLayer(id) == 2 && tTopo->pxbLadder(id)>=8 && tTopo->pxbLadder(id)<=21) {
//  if ( subdet == PixelSubdetector::PixelBarrel && tTopo->pxbLayer(id) == 4) {

    if (0) { // Building module from box using affin matrix for positioning
      float length = surface.bounds().width();
      float width =  surface.bounds().length();
      float height =  surface.bounds().thickness();
      
      // Rotated principle axes
      LocalVector lvx(1., 0., 0.);
      GlobalVector gvx = surface.toGlobal(lvx);
      LocalVector lvy(0., 1., 0.);
      GlobalVector gvy = surface.toGlobal(lvy);
      LocalVector lvz(0., 0., 1.);
      GlobalVector gvz = surface.toGlobal(lvz);
      LocalPoint lp0(length/-2., width/-2., height/-2.);
      GlobalPoint gp0 = surface.toGlobal(lp0);
      LocalPoint lp1(length/-2., width/-2., height/-1.);
      GlobalPoint gp1 = surface.toGlobal(lp1);

      printf("App.ActiveDocument.addObject(\"Part::Box\",\"Box%d\")\n", id.rawId());
      printf("FreeCAD.getDocument(\"Unnamed\").getObject(\"Box%d\").Length = %f\n", id.rawId(), length*10.);
      printf("FreeCAD.getDocument(\"Unnamed\").getObject(\"Box%d\").Width = %f\n", id.rawId(),  width*10.);
      printf("FreeCAD.getDocument(\"Unnamed\").getObject(\"Box%d\").Height = %f\n", id.rawId(), height*10.);
      //     printf("FreeCAD.getDocument(\"Unnamed\").getObject(\"Box%d\").Placement = FreeCAD.Placement(FreeCAD.Vector(%f,%f,%f), FreeCAD.Rotation(FreeCAD.Vector(%f,%f,%f), %f),FreeCAD.Vector(%f,%f,%f))\n", 
      // 	   id.rawId(), gp.x()-length/2., gp.y()-width/2., gp.z()-height/2., rot.x(), rot.y(), rot.z(), angle, length/2., width/2., height/2.);
      printf("FreeCAD.getDocument(\"Unnamed\").getObject(\"Box%d\").Placement = %s\n", id.rawId(), print(gvx, gvy, gvz, gp0).c_str());
      printf("Gui.getDocument(\"Unnamed\").getObject(\"Box%d\").ShapeColor=(0.0,1.0,0.0)\n", id.rawId());

      printf("App.ActiveDocument.addObject(\"Part::Box\",\"Box%d_dir\")\n", id.rawId());
      printf("FreeCAD.getDocument(\"Unnamed\").getObject(\"Box%d_dir\").Length = %f\n", id.rawId(), height*50.);
      printf("FreeCAD.getDocument(\"Unnamed\").getObject(\"Box%d_dir\").Width = %f\n", id.rawId(),  height*50.);
      printf("FreeCAD.getDocument(\"Unnamed\").getObject(\"Box%d_dir\").Height = %f\n", id.rawId(), height*100.);
      printf("FreeCAD.getDocument(\"Unnamed\").getObject(\"Box%d_dir\").Placement = %s\n", id.rawId(), print(gvx, gvy, gvz, gp1).c_str());
      printf("Gui.getDocument(\"Unnamed\").getObject(\"Box%d_dir\").ShapeColor=(1.0,0.0,0.0)\n", id.rawId());
    }

    if (0) { // Building module from meshes using 4 points (in local x-y-z, x points right, y points up): bottom-lower-right, bottom-upper-right, bottom-upper-left, top-upper-left
      LocalPoint lp1(surface.bounds().width()/2.0, surface.bounds().length()/-2.0, surface.bounds().thickness()/-2.0); // width is short edge along local x, length is long edge along local y
      LocalPoint lp2(surface.bounds().width()/2.0, surface.bounds().length()/2.0, surface.bounds().thickness()/-2.0);
      LocalPoint lp3(surface.bounds().width()/-2.0, surface.bounds().length()/2.0, surface.bounds().thickness()/-2.0);
      LocalPoint lp4(surface.bounds().width()/-2.0, surface.bounds().length()/2.0, surface.bounds().thickness()/2.0);
      GlobalPoint gp1 = surface.toGlobal(lp1);
      GlobalPoint gp2 = surface.toGlobal(lp2);
      GlobalPoint gp3 = surface.toGlobal(lp3);
      GlobalPoint gp4 = surface.toGlobal(lp4);
      printf("module%d=module(\"module%d\", %f,%f,%f, %f,%f,%f, %f,%f,%f, %f,%f,%f)\n",  id.rawId(),  id.rawId(), gp1.x(),gp1.y(),gp1.z(), gp2.x(),gp2.y(),gp2.z(), gp3.x(),gp3.y(),gp3.z(), gp4.x(),gp4.y(),gp4.z());
    }

  }

  return strstr.str();
}

std::string Phase1GeometryTopologyAnalyzer::print(const GlobalVector& gvx, const GlobalVector& gvy, const GlobalVector& gvz, const GlobalPoint& gp) const {
    std::stringstream strstr;
    strstr << "FreeCAD.Matrix(";
//     strstr << gvx.x()<<", "<<gvx.y()<<", "<<gvx.z()<<", "<<gp.x()<<",";
//     strstr << gvy.x()<<", "<<gvy.y()<<", "<<gvy.z()<<", "<<gp.y()<<",";
//     strstr << gvz.x()<<", "<<gvz.y()<<", "<<gvz.z()<<", "<<gp.z()<<",";
    strstr << gvx.x()<<", "<<gvy.x()<<", "<<gvz.x()<<", "<<gp.x()*10.<<",";
    strstr << gvx.y()<<", "<<gvy.y()<<", "<<gvz.y()<<", "<<gp.y()*10.<<",";
    strstr << gvx.z()<<", "<<gvy.z()<<", "<<gvz.z()<<", "<<gp.z()*10.<<",";
    strstr << "0.,0.,0.,1.)";
    return strstr.str();
}

void Phase1GeometryTopologyAnalyzer::analyze( const edm::Event &iEvent, const edm::EventSetup& iSetup) {

  edm::ESHandle<TrackerTopology> tTopoHandle;
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
  const TrackerTopology* const tTopo = tTopoHandle.product();

  typedef std::vector<DetId>                 DetIdContainer;

  edm::ESHandle<TrackerGeometry> tracker;
  iSetup.get<TrackerDigiGeometryRecord>().get(tracker);    
  const TrackerGeometry *tkGeom = &(*tracker);

  DetIdContainer allIds=tkGeom->detIds();

  for( DetIdContainer::const_iterator id = allIds.begin(), detUnitIdEnd = allIds.end(); id != detUnitIdEnd; ++id ) {
    if ( id->det()==DetId::Tracker ) {
      unsigned int subdet=id->subdetId();
      
      std::string DetIdPrint;
      std::string GeomPrint;

      if (subdet == PixelSubdetector::PixelBarrel) {
	DetIdPrint = print(tTopo, *id);
	GeomPrint = print(tTopo,tkGeom, *id);
	std::cout << DetIdPrint << " " << GeomPrint << std::endl;
      }
      else if (subdet == PixelSubdetector::PixelEndcap) {
	DetIdPrint = print(tTopo, *id);
	GeomPrint = print(tTopo,tkGeom, *id);
	std::cout << DetIdPrint << " " << GeomPrint << std::endl;
      }
    }
  }

}

DEFINE_FWK_MODULE(Phase1GeometryTopologyAnalyzer);

