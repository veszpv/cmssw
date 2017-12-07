///////////////////////////////////////////////////////////////////////////////
// File: DDTelescopePlanesAlgo.cc
// Description:  
///////////////////////////////////////////////////////////////////////////////


#include <cmath>
#include <algorithm>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DetectorDescription/Core/interface/DDCurrentNamespace.h"
#include "DetectorDescription/Core/interface/DDSplit.h"
#include "Geometry/PixelTelescope/plugins/DDTelescopePlanesAlgo.h"
#include "DetectorDescription/Core/interface/DDRotationMatrix.h"
#include "DetectorDescription/Core/interface/DDTransform.h"
#include "CLHEP/Units/GlobalPhysicalConstants.h"
#include "CLHEP/Units/GlobalSystemOfUnits.h"


DDTelescopePlanesAlgo::DDTelescopePlanesAlgo() {
  LogDebug("TrackerGeom") << "DDTelescopePlanesAlgo info: Creating an instance";
}


DDTelescopePlanesAlgo::~DDTelescopePlanesAlgo() {}


void DDTelescopePlanesAlgo::initialize(const DDNumericArguments & nArgs,
				  const DDVectorArguments & vArgs,
				  const DDMapArguments & ,
				  const DDStringArguments & sArgs,
				  const DDStringVectorArguments & ) {

  n             = int(nArgs["N"]);
  tiltAngle     = nArgs["tiltAngle"];
  skewAngle     = nArgs["skewAngle"];
  deltaZ        = nArgs["deltaZ"];
  
  LogDebug("TrackerGeom") << "DDTelescopePlanesAlgo debug: Parameters for position"
			  << "ing:: n " << n << " Phase 1 modules with deltaZ "
			  << deltaZ << ", tiltAngle "
			  << tiltAngle/CLHEP::deg << ", skew angle " 
			  << skewAngle/CLHEP::deg;

  idNameSpace = DDCurrentNamespace::ns();
  childName   = sArgs["ChildName"];

  DDName parentName = parent().name();
  LogDebug("TrackerGeom") << "DDTelescopePlanesAlgo debug: Parent " << parentName
			  << "\tChild " << childName << " NameSpace "
			  << idNameSpace;
}


void DDTelescopePlanesAlgo::execute(DDCompactView& cpv) {

  DDRotation prepaRot, tiltRot, skewRot, globalRot; // Identity
  DDRotationMatrix prepaMatrix, tiltMatrix, skewMatrix, globalRotMatrix; // Identity matrix
  std::string rotstr = "RTelescopePlanesAlgo";

  // prepaMatrix calculus
  std::string prepaRotstr = rotstr + "Prepa";
  prepaRot = DDRotation(DDName(prepaRotstr, idNameSpace));
  if (!prepaRot) {
    LogDebug("TrackerGeom") << "DDTelescopePlanesAlgo test: Creating a new rotation: " << prepaRotstr
			    << "\t180., 0., "
			    << "90., 90., "
			    << "90., 0.";
    prepaRot = DDrot(DDName(prepaRotstr, idNameSpace), 
		     180.*CLHEP::deg, 0., 90.*CLHEP::deg, 90.*CLHEP::deg, 90.*CLHEP::deg, 0.);
  }
  prepaMatrix = *prepaRot.matrix();

  // tiltMatrix calculus
  std::string tiltRotstr = rotstr + "Tilt" + std::to_string(tiltAngle/CLHEP::deg);
  tiltRot = DDRotation(DDName(tiltRotstr, idNameSpace));
  if (!tiltRot) {
    LogDebug("TrackerGeom") << "DDTelescopePlanesAlgo test: Creating a new rotation: " << tiltRotstr
			    << "\t90., 0., "
			    << 90. - tiltAngle/CLHEP::deg << ", 90., "
			    << tiltAngle/CLHEP::deg << ", 270.";
    tiltRot = DDrot(DDName(tiltRotstr, idNameSpace), 
		    90.*CLHEP::deg, 0., 90.*CLHEP::deg - tiltAngle, 90.*CLHEP::deg, tiltAngle, 270.*CLHEP::deg);
  }
  tiltMatrix = *tiltRot.matrix();
  tiltMatrix *= prepaMatrix;

  // skewMatrix calculus
  std::string skewRotstr = rotstr + "Skew" + std::to_string(skewAngle/CLHEP::deg);
  skewRot = DDRotation(DDName(skewRotstr, idNameSpace));
  if (!skewRot) {
    LogDebug("TrackerGeom") << "DDTelescopePlanesAlgo test: Creating a new rotation: " << skewRotstr
			    << "\t" << 90. + skewAngle/CLHEP::deg << ", 0., "
			    << "90., 90., "
			    << skewAngle/CLHEP::deg << ", 0.";
    skewRot = DDrot(DDName(skewRotstr, idNameSpace), 
		    90.*CLHEP::deg + skewAngle, 0., 90.*CLHEP::deg, 90.*CLHEP::deg, skewAngle, 0.);
  }
  skewMatrix = *skewRot.matrix();
  skewMatrix *= tiltMatrix;

  // globalRot def
  std::string globalRotstr = rotstr + "Global";
  globalRot = DDRotation(DDName(globalRotstr, idNameSpace));
  if (!globalRot) {
    LogDebug("TrackerGeom") << "DDTelescopePlanesAlgo test: Creating a new "
			    << "rotation: " << globalRotstr;
    globalRotMatrix = skewMatrix;
    globalRot = DDrot(DDName(globalRotstr, idNameSpace), new DDRotationMatrix(globalRotMatrix));
  }

  // Loops for all n modules
  DDName mother = parent().name();
  DDName child(DDSplit(childName).first, DDSplit(childName).second);
  int    copy   = 1;

  for (int i=0; i<n; i++) {
  
    // translation def
    double xpos = 0.;
    double ypos = 0.;
    double zpos = (-1.5 + i) * deltaZ;
    DDTranslation tran(xpos, ypos, zpos);
  
    // Positions child with respect to parent
    cpv.position(child, mother, copy, tran, globalRot);
    LogDebug("TrackerGeom") << "DDTelescopePlanesAlgo test " << child << " number "
			    << copy << " positioned in " << mother << " at "
			    << tran  << " with " << globalRot;

    copy += 1;
  }
}
