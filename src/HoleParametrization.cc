#include "HoleParametrization.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HoleParametrization::HoleParametrization(
                                            std::vector<G4double> hole_x,
                                            std::vector<G4double> hole_y,
                                            std::vector<G4double> hole_z,
                                            G4ThreeVector hole_dimensions,
                                            std::vector<G4double> crystal_x,
                                            std::vector<G4double> crystal_y,
                                            std::vector<G4double> crystal_z,
                                            G4ThreeVector crystal_dimensions,
                                             int moduleType,
                                             int material
                                              )
 : G4VPVParameterisation()
{
  // fnModules = nModules;
  fHole_x = hole_x;
  fHole_y = hole_y;
  fHole_z = hole_z;
  fHoleDimensions = hole_dimensions;
  fCrystal_x = crystal_x;
  fCrystal_y = crystal_y;
  fCrystal_z = crystal_z;
  fCrystalDimensions = crystal_dimensions;
  fType = moduleType;
  fMaterial = material;
  // fSeparation_z = module_separation_z;
  // fSections = module_sections;


  // call method to save modules in struct
  SaveElements();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

HoleParametrization::~HoleParametrization()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HoleParametrization::ComputeTransformation
(const G4int copyNo,
 G4VPhysicalVolume* physVol) const
{
  // // Note: copyNo will start with zero!
  // G4double Zposition = fStartZ + copyNo * fSpacing;
  G4double x = fHole_x[copyNo];
  G4double y = fHole_y[copyNo];
  G4double z = fHole_z[copyNo];
  G4ThreeVector origin(x,y,z);
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void HoleParametrization::ComputeDimensions
(G4Box &module,
 const G4int copyNo,
 const G4VPhysicalVolume* physVol) const
{
  module.SetXHalfLength (fHoleDimensions.getX()/2.0);
  module.SetYHalfLength (fHoleDimensions.getY()/2.0);
  module.SetZHalfLength (fHoleDimensions.getZ()/2.0);
  // // Note: copyNo will start with zero!
}


void HoleParametrization::SaveElements()
{
  for(int i = 0; i < fHole_x.size(); i++)
  {
    if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveStructure)
    {

      // save the hole TTree in output file
      CreateTree::Instance()->holeID           = i;
      CreateTree::Instance()->holeMaterial     = 0;
      CreateTree::Instance()->holeType         = fType;
      CreateTree::Instance()->holeX            = fHole_x[i];
      CreateTree::Instance()->holeY            = fHole_y[i];
      CreateTree::Instance()->holeZ            = fHole_z[i];
      CreateTree::Instance()->holeDX           = fHoleDimensions.getX();
      CreateTree::Instance()->holeDY           = fHoleDimensions.getY();
      CreateTree::Instance()->holeDZ           = fHoleDimensions.getZ();
      CreateTree::Instance()->holeSeparationZ  = 0;
      CreateTree::Instance()->holeSections     = 1;
      CreateTree::Instance()->holes->Fill();

      // save the crystal in the fibre TTree in output file
      CreateTree::Instance()->fibreID          = i;
      CreateTree::Instance()->fibreMaterial    = (Int_t) fMaterial;
      CreateTree::Instance()->fibreType        = fType;
      CreateTree::Instance()->fibreX           = fCrystal_x[i];
      CreateTree::Instance()->fibreY           = fCrystal_y[i];
      CreateTree::Instance()->fibreZ           = fCrystal_z[i];
      CreateTree::Instance()->fibreDX          = fCrystalDimensions.getX();
      CreateTree::Instance()->fibreDY          = fCrystalDimensions.getY();
      CreateTree::Instance()->fibreDZ          = fCrystalDimensions.getZ();
      CreateTree::Instance()->fibreSeparationZ = 0;
      CreateTree::Instance()->fibreSections    = 1;
      CreateTree::Instance()->fibres->Fill();
      
    }
  }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
