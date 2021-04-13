// Det const
// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch

//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include <iostream>
#include <string>
#include <fstream>
#include <utility>

#include "ConfigFile.hh"
#include "MyMaterials.hh"
#include "LedFiberTiming.hh"
#include "Parameters.hh"
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



#include "Absorber.hh"
#include "Spacal.hh"


class DetectorConstruction : public G4VUserDetectorConstruction
{
public:

  //! ctor
  DetectorConstruction  () ;
  DetectorConstruction  (const string& configFileName) ;

  //! dtor
  ~DetectorConstruction () ;

  //! construct method
  G4VPhysicalVolume* Construct () ;

  // build a module
  // void BuildModule(G4LogicalVolume * moduleLV) ;

  // module positioning with replicas
  // void PositionModulesWithReplica(G4LogicalVolume * moduleLV,G4LogicalVolume * calorimeterLV);


  // ! other methods
  G4double GetModule_z () const { return 0.0 ; } ; // FIXME hardocoded now because of multi modules. What was its purpose?

  void initializeMaterials () ;
  void ConstructField () ;
  // spacal 2018 test beam readout
  // void ConstructLightGuides(G4LogicalVolume *Interface_LV, G4PVPlacement *Interface_PV, std::string moduleName,G4double lguide_edge,G4RotationMatrix *rot);
  // // void JustAirGap(G4LogicalVolume *Interface_LV, G4PVPlacement *Interface_PV, std::string moduleName,G4double lguide_edge,G4RotationMatrix *rot);
  // void ConstructPMTs(G4LogicalVolume *Readout_LV, std::string moduleName,G4double lg_pitch,G4RotationMatrix *rot);
  //
  // void ConstructOneLargePMT(G4LogicalVolume *Readout_LV, std::string moduleName,G4double PMT_radius,G4RotationMatrix *rot);
  // // spacal 2019 test beam readout
  // // void ConstructLightGuides2019(G4LogicalVolume *Interface_LV, G4PVPlacement *Interface_PV, std::string moduleName,G4RotationMatrix *rot);
  // // void ConstructPMTSforTB2019(G4LogicalVolume *Readout_LV, std::string moduleName,G4RotationMatrix *rot);
  //
  // void ConstructClearFiber(G4LogicalVolume *Interface_LV, G4PVPlacement *Interface_PV, G4PVPlacement *Gap_Interface_PV, G4PVPlacement *Gap_Interface_Readout_PV, std::string moduleName);
  // // include clear fiber from kuraray, 10 cm long clear_PSM, 1mm diameter


  Fiber* GetFiber() { return &fib ; } ;


private:
  G4bool    checkOverlaps ;

  bool saveAll            ;
  bool saveTree           ;
  bool saveShower         ;
  bool saveStructure      ;
  bool savePrimaries      ;
  bool savePhotons        ;
  bool savePhotonAbsPoint ;
  bool saveLAPPD          ;

  G4double  expHall_x ;
  G4double  expHall_y ;
  G4double  expHall_z ;

  G4int    world_material ;    // world material
  G4double W_fraction ;      // fraction of Tungsten in the alloy

  G4int surface_lg;
  G4int glue_interface;
  G4int cone_material;
  bool esr_on_cones;

  int pipe_modules_nx;
  int pipe_modules_ny;


  G4double scaleFactor;



  G4double fibre_cladRIndex;
  G4double fibre_radius ;
  G4double fibre_length ;
  G4double fibre_absLength ;   // absorption length in the fiber

  G4int crystalsVisibility;
  G4int worldVisibility;
  G4int caloVisibility;
  G4int moduleVisibility;
  G4int interfaceVisibility;
  G4int holeVisibility;
  G4int gapsVisibility;
  G4int readoutVisibility;
  G4int absorberVisibility;
  G4int esrVisibility;
  G4int lgVisibility;


  G4int wireFrame;
  G4double crystal_lateral_depolishing;
  G4double crystal_exit_depolishing;
  G4int cell_separation_type;
  G4double cell_separator_position;
  G4double separation_thickness;
  G4int separation_material;
  G4double esrTransmittance;

  std::vector<G4double> attLengths;

  G4UniformMagField * B_field ;
  G4bool   B_field_IsInitialized ;
  G4double B_field_intensity ;     // magnetic field, in units of Tesla

  Fiber fib ;

  //Materials
  G4Material* WoMaterial ;
  G4Material* AbMaterial ;
  G4Material* AbMaterial2;
  G4Material* CoMaterial ;
  G4Material* ClMaterial ;
  G4Material* ClSSMaterial;
  G4Material* ClSSSMaterial;
  G4Material* PLEXMaterial;
  G4Material* PVCMaterial;
  G4Material* PMTMaterial;
  G4Material* WiresMaterial;
  G4Material* PlaneMaterial;

  G4Material* AirMaterial;
  G4Material* AirKillerMaterial;
  G4Material* InterfaceContainerMaterial;

  G4Material* LAPPD_average;

  G4Material* GlueMaterial ;
  G4Material* InterfaceMaterial ;

  G4Material* GapAbsToInterfaceMaterial;
  G4Material* GapInterfaceToReadoutMaterial;

  G4Material* Cl3Material ;
  G4Material* Cl4SSMaterial;
  G4Material* Cl43SSMaterial;

  G4Material* GaMaterial ;
  G4Material* DeMaterial ;

  G4Material* SeparationMaterial;

  G4MaterialPropertiesTable* absorber_optical_surface;
  G4double absorber_reflectivity;
  G4double absorber_specularLobe;
  G4double absorber_specularSpike;
  G4double absorber_backScatter ;

  Absorber absorber;

  G4int modules_nx;
  G4int modules_ny;
  std::vector<G4int>      ForbiddenModules;

  G4double LGrelevantDimension;

  G4double module_pos_x;
  G4double module_pos_y;
  G4double module_pos_z;

  G4double module_size_x;
  G4double module_size_y;
  G4double module_size_z;

  G4double calorimeter_size_x;
  G4double calorimeter_size_y;
  G4double calorimeter_size_z;

  G4int calorimeter_position;
  G4double calorimeter_pos_x;
  G4double calorimeter_pos_y;
  G4double calorimeter_pos_z;

  G4double gapSize ;
  G4double InterfaceSizeX;
  G4double InterfaceSizeY;
  G4double InterfaceSizeZ;
  G4double ReadoutSizeX;
  G4double ReadoutSizeY;
  G4double ReadoutSizeZ;

  G4String   AbsName;
  G4double   AbsSizeX;
  G4double   AbsSizeY;
  G4double   AbsSizeZ;
  G4double   AbsPositionX;
  G4double   AbsPositionY;
  G4double   AbsPositionZ;
  G4int      AbsMaterial;

  G4int gap_abs_interface_material;
  G4int gap_interface_readout_material;

  G4int readoutType;


  G4Colour  white   ;
  G4Colour  grey    ;
  G4Colour  black   ;
  G4Colour  red     ;
  G4Colour  green   ;
  G4Colour  palegreen   ;
  G4Colour  blue    ;
  G4Colour  cyan    ;
  G4Colour  air     ;
  G4Colour  magenta ;
  G4Colour  yellow  ;
  G4Colour  brass   ;
  G4Colour  brown   ;
  G4Colour  orange  ;
  G4Colour  paleorange  ;

  std::vector<G4Colour> color_vector;
  bool containerVisibility;

  std::vector<Spacal*> modules_list;
  G4double   module_block_size_x;
  G4double   module_block_size_y;
  G4double   module_block_size_z;

} ;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*DetectorConstruction_h*/
