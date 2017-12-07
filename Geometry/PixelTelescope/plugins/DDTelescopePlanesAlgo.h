#ifndef DD_TelescopePlanesAlgo_h
#define DD_TelescopePlanesAlgo_h

/*
  
*/

#include <map>
#include <string>
#include <vector>
#include "DetectorDescription/Core/interface/DDTypes.h"
#include "DetectorDescription/Core/interface/DDAlgorithm.h"

class DDTelescopePlanesAlgo : public DDAlgorithm {
 
public:
  // Constructor and Destructor
  DDTelescopePlanesAlgo(); 
  ~DDTelescopePlanesAlgo() override;
  
  void initialize(const DDNumericArguments & nArgs,
		  const DDVectorArguments & vArgs,
		  const DDMapArguments & mArgs,
		  const DDStringArguments & sArgs,
		  const DDStringVectorArguments & vsArgs) override;

  void execute(DDCompactView& cpv) override;

private:

  int           n;              //Number of copies
  double        tiltAngle;      //
  double        skewAngle;      //
  double        deltaZ;         //

  std::string   idNameSpace;    //Namespace of this and ALL sub-parts
  std::string   childName;      //Child name
};

#endif
