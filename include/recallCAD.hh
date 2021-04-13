// Class for cell elements
// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch

#ifndef recallCAD_h
#define recallCAD_h 1

#include <iostream>
#include <string>
#include <fstream>
#include <utility>

#include "ConfigFile.hh"
#include "MyMaterials.hh"
#include "LedFiberTiming.hh"
// #include "DetectorParameterisation.hh"

#include "G4Material.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4UniformMagField.hh"



// class for Module Cells
class recallCAD
{
public:

  //! ctor
  recallCAD ();

  //! dtor
  ~recallCAD ();

  //! construct method

  G4LogicalVolume* returnMeshedLV(std::string OBJfile, G4Material* Material, const char*);

} ;

#endif /*CADrecall_h*/
