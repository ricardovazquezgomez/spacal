// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch


#include "Absorber.hh"
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

using namespace CLHEP;

// absorber_optical_surface = NULL;
//
// absorber_optical_surface = MyMaterials::ABS_SURF(absorber_reflectivity,absorber_specularLobe, absorber_specularSpike, absorber_backScatter);
// if(absorber_reflectivity != -1)
// {
//   G4cout << "Absorber_optical_surface Reflectivity set to " << absorber_reflectivity << G4endl ;
// }
//
// SeparationMaterial = NULL;
// if      ( separation_material == 1 ) SeparationMaterial = MyMaterials::Brass () ;
// else if ( separation_material == 2 ) SeparationMaterial = MyMaterials::Tungsten () ;
// else if ( separation_material == 3 ) SeparationMaterial = MyMaterials::Lead () ;
// else if ( separation_material == 4 ) SeparationMaterial = MyMaterials::Iron () ;
// else if ( separation_material == 5 ) SeparationMaterial = MyMaterials::Aluminium () ;
// else if ( separation_material == 6 ) SeparationMaterial = MyMaterials::CopperTungstenAlloy(W_fraction) ;
// else if ( separation_material == 7 ) SeparationMaterial = MyMaterials::PureTungsten1() ;
// else if ( separation_material == 8 ) SeparationMaterial = MyMaterials::PureTungsten2() ;
// else if ( separation_material == 9 ) SeparationMaterial = MyMaterials::AirKiller() ;
// else if ( separation_material == 10 ) SeparationMaterial = MyMaterials::Air() ;
// else
// {
//   G4cerr << "<Absorber>: Invalid separation material specifier " << separation_material << G4endl ;
//   exit (-1) ;
// }

Absorber::Absorber()
{
  // constructor
  // doesn't actually do anything, all will be set by the methods, called in DetectorConstruction

}

Absorber::~Absorber(){}


void Absorber::AddCell(Cell cell)
{
  cells.push_back(cell);
}

Cell Absorber::GetCell(int id)
{
  return cells[id];
}

G4int Absorber::GetNumberOfCells()
{
  return cells.size();
}


void Absorber::SetMaterialIdentifier(G4int mat)
{
  MaterialIdentifier = mat;
}

void Absorber::SetMaterial(G4int mat,G4double W_fraction)
{
  SetMaterialIdentifier(mat);
  AbMaterial = NULL ;
  if      ( mat == 1 ) AbMaterial = MyMaterials::Brass () ;
  else if ( mat == 2 ) AbMaterial = MyMaterials::Tungsten () ;
  else if ( mat == 3 ) AbMaterial = MyMaterials::Lead () ;
  else if ( mat == 4 ) AbMaterial = MyMaterials::Iron () ;
  else if ( mat == 5 ) AbMaterial = MyMaterials::Aluminium () ;
  else if ( mat == 6 ) AbMaterial = MyMaterials::CopperTungstenAlloy(W_fraction) ;
  else if ( mat == 7 ) AbMaterial = MyMaterials::PureTungsten1() ;
  else if ( mat == 8 ) AbMaterial = MyMaterials::PureTungsten2() ;
  else if ( mat == 9 ) AbMaterial = MyMaterials::AirKiller() ;
  else if ( mat == 10 ) AbMaterial = MyMaterials::Air() ;
  else if ( mat == 11 ) AbMaterial = MyMaterials::LAPPD_Average() ;
  else if ( mat == 12 ) AbMaterial = MyMaterials::StainlessSteel() ;
  else if ( mat == 13 ) AbMaterial = MyMaterials::GarthTypographicAlloy() ;
  else
  {
    G4cerr << "<Absorber>: Invalid absorber material specifier " << mat << G4endl ;
    exit (-1) ;
  }
  // G4cout << "Absorber " << name << " material: "<< AbMaterial << G4endl ;
}

void Absorber::SetSeparationMaterial(G4int mat,G4double W_fraction)
{
  SeparationMaterial = NULL;
  if(SingleSeparationLayer)
  {
    if      ( mat == 1 )  SeparationMaterial = MyMaterials::Brass () ;
    else if ( mat == 2 )  SeparationMaterial = MyMaterials::Tungsten () ;
    else if ( mat == 3 )  SeparationMaterial = MyMaterials::Lead () ;
    else if ( mat == 4 )  SeparationMaterial = MyMaterials::Iron () ;
    else if ( mat == 5 )  SeparationMaterial = MyMaterials::Aluminium () ;
    else if ( mat == 6 )  SeparationMaterial = MyMaterials::CopperTungstenAlloy(W_fraction) ;
    else if ( mat == 7 )  SeparationMaterial = MyMaterials::PureTungsten1() ;
    else if ( mat == 8 )  SeparationMaterial = MyMaterials::PureTungsten2() ;
    else if ( mat == 9 )  SeparationMaterial = MyMaterials::AirKiller() ;
    else if ( mat == 10 ) SeparationMaterial = MyMaterials::Air() ;
    else if ( mat == 11 ) SeparationMaterial = MyMaterials::LAPPD_Average() ;
    else if ( mat == 12 ) SeparationMaterial = MyMaterials::StainlessSteel() ;
    else if ( mat == 13 ) SeparationMaterial = MyMaterials::GarthTypographicAlloy() ;
    else
    {
      G4cerr << "<Separation>: Invalid separation material specifier " << mat << G4endl ;
      exit (-1) ;
    }
  }
  else // in multi separation layer mode, set a container volume of air 
  {
    SeparationMaterial = MyMaterials::Air() ;
  }
  
  G4cout << "Separation " << name << " material: "<< SeparationMaterial << G4endl ;
}

void Absorber::SetLAPPDMaterials(std::vector<G4int> mat,G4double W_fraction)
{
  LAPPD_layers_materials.clear();
  for (auto & elements : mat) {
    SeparationMaterial = NULL;
    if      ( elements == 1 )  LAPPD_layers_materials.push_back( MyMaterials::Air () );
    else if ( elements == 2 )  LAPPD_layers_materials.push_back( MyMaterials::LAPPD_Window() );
    else if ( elements == 3 )  LAPPD_layers_materials.push_back( MyMaterials::Vacuum () );
    else if ( elements == 4 )  LAPPD_layers_materials.push_back( MyMaterials::LAPPD_MCP () );
    else if ( elements == 5 )  LAPPD_layers_materials.push_back( MyMaterials::LAPPD_PCB () );
    else if ( elements == 6 )  LAPPD_layers_materials.push_back( MyMaterials::Aluminium () );
    else if ( elements == 7 )  LAPPD_layers_materials.push_back( MyMaterials::StainlessSteel () );
    else if ( elements == 8 )  LAPPD_layers_materials.push_back( MyMaterials::ESR_Vikuiti () );
   else
    {
      G4cerr << "<LAPPD>: Invalid LAPPD material specifier " << elements << G4endl ;
      exit (-1) ;
    }
  }
}

void Absorber::SetOpticalSurface(double reflectivity,double specularLobe,double specularSpike,double backScatter,double sigma_alpha_v)
{
  absorber_optical_surface = MyMaterials::ABS_SURF(reflectivity,specularLobe,specularSpike,backScatter);
  sigma_alpha = sigma_alpha_v;
  if(reflectivity != -1)
  {
    G4cout << "Absorber_optical_surface Reflectivity set to " << reflectivity << G4endl ;
  }
}
