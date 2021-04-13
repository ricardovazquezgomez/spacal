// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch

#include "recallCAD.hh"

#include "CreateTree.hh"

#include <algorithm>
#include <string>
#include <sstream>

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4ExplicitEuler.hh"
#include "G4ChordFinder.hh"
#include "G4EqMagElectricField.hh"
#include "G4PropagatorInField.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"
#include "G4PVParameterised.hh"
#include <G4TransportationManager.hh>
#include <G4MagneticField.hh>
#include <G4UniformMagField.hh>
#include <G4FieldManager.hh>
#include <G4Cons.hh>
#include "CADMesh.hh"

using namespace CLHEP;

recallCAD::recallCAD()
{
  // constructor
  // doesn't actually do anything, all will be set by the methods
}

recallCAD::~recallCAD(){}


G4LogicalVolume* recallCAD::returnMeshedLV(std::string OBJfile, G4Material* Material, const char*)
{
    auto vol_mesh = CADMesh::TessellatedMesh::FromOBJ(OBJfile);
    G4LogicalVolume* LV = new G4LogicalVolume (vol_mesh->GetSolid(), Material, "test_to_be_changed");
    return LV;
}