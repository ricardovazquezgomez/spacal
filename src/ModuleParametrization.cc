#include "ModuleParametrization.hh"

#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4SystemOfUnits.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModuleParametrization::ModuleParametrization(G4int    nModules,
                                             std::vector<G4double> x,
                                             std::vector<G4double> y,
                                             std::vector<G4double> z,
                                             G4ThreeVector dimensions,
                                             int moduleType,
                                             double module_separation_z,
                                             int module_sections  )
 : G4VPVParameterisation()
{
  fnModules = nModules;
  fx = x;
  fy = y;
  fz = z;
  fDimensions = dimensions;
  fType = moduleType;
  fSeparation_z = module_separation_z;
  fSections = module_sections;


  // call method to save modules in struct
  SaveModules();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ModuleParametrization::~ModuleParametrization()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModuleParametrization::ComputeTransformation
(const G4int copyNo,
 G4VPhysicalVolume* physVol) const
{
  // // Note: copyNo will start with zero!
  // G4double Zposition = fStartZ + copyNo * fSpacing;
  G4double x = fx[copyNo];
  G4double y = fy[copyNo];
  G4double z = fz[copyNo];
  G4ThreeVector origin(x,y,z);
  physVol->SetTranslation(origin);
  physVol->SetRotation(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModuleParametrization::ComputeDimensions
(G4Box &module,
 const G4int copyNo,
 const G4VPhysicalVolume* physVol) const
{
  module.SetXHalfLength (fDimensions.getX()/2.0);
  module.SetYHalfLength (fDimensions.getY()/2.0);
  module.SetZHalfLength (fDimensions.getZ()/2.0);
  // // Note: copyNo will start with zero!
  // G4double rmax = fRmaxFirst + copyNo * fRmaxIncr;
  // trackerChamber.SetInnerRadius(0);
  // trackerChamber.SetOuterRadius(rmax);
  // trackerChamber.SetZHalfLength(fHalfWidth);
  // trackerChamber.SetStartPhiAngle(0.*deg);
  // trackerChamber.SetDeltaPhiAngle(360.*deg);
}


void ModuleParametrization::SaveModules()
{
  for(int i = 0; i < fx.size(); i++)
  {
    if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveStructure)
    {
      CreateTree::Instance()->moduleID                = i;
      CreateTree::Instance()->moduleMaterial          = 0;
      CreateTree::Instance()->moduleType              = fType;
      CreateTree::Instance()->moduleX                 = fx[i];
      CreateTree::Instance()->moduleY                 = fy[i];
      CreateTree::Instance()->moduleZ                 = fz[i];
      CreateTree::Instance()->moduleDX                = fDimensions.getX();
      CreateTree::Instance()->moduleDY                = fDimensions.getY();
      CreateTree::Instance()->moduleDZ                = fDimensions.getZ();
      CreateTree::Instance()->moduleSeparationZ       = fSeparation_z;
      CreateTree::Instance()->moduleSections          = fSections;
      CreateTree::Instance()->modules->Fill();
    }
  }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
