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

#ifndef Spacal_h
#define Spacal_h 1

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


#include "Absorber.hh"
#include "Parameters.hh"

// class of Spacal module
class Spacal
{
  public:

    //ctor(s)
    Spacal();
    Spacal(Data_t input_parameters);

    //dtor
    ~Spacal();

    void ComputeVolumes();
    // absorber
    void ComputeAbsorber();

    // cells
    void ComputeCells();

    // interface
    void ComputeInterface();

    // readout
    void ComputeReadout();

    void ComputeModule();

    void BuildModule(G4LogicalVolume * moduleLV);
    // void BuildInterfaces(G4LogicalVolume * moduleLV);
    // void BuildReadouts(G4LogicalVolume * moduleLV);

    

    G4LogicalVolume* MakeHole(Cell cell );
    void InitializeMaterials();
    // void SaveCrystalMaterials();

    void ConstructLightGuides(G4LogicalVolume *Interface_LV, G4PVPlacement *Interface_PV, std::string moduleName,G4double lguide_edge,G4RotationMatrix *rot, G4double shiftCAD);
    void ConstructPMTs(G4LogicalVolume *Readout_LV, std::string moduleName,G4double lg_pitch,G4RotationMatrix *rot);
    void ConstructOneLargePMT(G4LogicalVolume *Readout_LV, std::string moduleName,G4double PMT_radius,G4RotationMatrix *rot);
    void ConstructClearFiber(G4LogicalVolume *Interface_LV, G4PVPlacement *Interface_PV,  G4PVPlacement *Gap_Interface_PV,  G4PVPlacement *Gap_Interface_Readout_PV, std::string moduleName);

    G4double GetModulePositionX(){return module_pos_x;};
    G4double GetModulePositionY(){return module_pos_y;};
    G4double GetModulePositionZ(){return module_pos_z;};
    G4double GetModuleSizeX()    {return module_size_x;};
    G4double GetModuleSizeY()    {return module_size_y;};
    G4double GetModuleSizeZ()    {return module_size_z;};

    G4double GetXmapMin(){return xMapMin;};
    G4double GetXmapMax(){return xMapMax;};
    G4double GetYmapMin(){return yMapMin;};
    G4double GetYmapMax(){return yMapMax;};

    G4LogicalVolume *GetModuleLV(){return moduleLV;};

    G4double GetModuleZshift(){return moduleZshift;};

    int GetECALposition(){return ecal_position;};
    double GetSeparationPositionZ(){return separation_z_position;};
    int  GetNofSections(){return sections;};
    bool GetModuleVisibility(){return module_visibility;};
    bool GetModuleWireFrame(){return module_wireframe;};
    std::string GetContainerVolumeFileName(){return volume_file;};
    std::vector<Parameters::Point> GetCoordinates(){return coordinates;};

    void ScanMap(TH2D *mapEcal);



  private:

    std::stringstream sname;
    std::string moduleName;

    //colors
    G4Colour  white   ;
    G4Colour  grey    ;
    G4Colour  black   ;
    G4Colour  red     ;
    G4Colour  green   ;
    G4Colour  blue    ;
    G4Colour  cyan    ;
    G4Colour  air     ;
    G4Colour  magenta ;
    G4Colour  yellow  ;
    G4Colour  brass   ;
    G4Colour  brown   ;
    G4Colour  orange  ;

    // data from config file
    Data_t parameters;

    // specific data saved for this module
    Parameters::Module_t module;

    // Absorber
    Absorber absorber;

    // dimensions
    G4double module_pos_x;
    G4double module_pos_y;
    G4double module_pos_z;
    G4double module_size_x;
    G4double module_size_y;
    G4double module_size_z;
    G4double AbsSizeX;
    G4double AbsSizeY;
    G4double AbsSizeZ;
    G4double InterfaceSizeX;
    G4double InterfaceSizeY;
    G4double InterfaceSizeZ;
    G4double ReadoutSizeX;
    G4double ReadoutSizeY;
    G4double ReadoutSizeZ;
    G4double LGrelevantDimension;

    G4double xMapMin;
    G4double xMapMax;
    G4double yMapMin;
    G4double yMapMax;

    // position type
    int ecal_position;

    // relative z position of cell sepation
    double separation_z_position;

    // number of longitudinal sections
    int sections;

    // container volume
    std::string volume_file;

    // module volume visibility
    bool module_visibility;
    bool module_wireframe;

    std::vector<Parameters::Point> coordinates;

    // solid and logic volumes
    G4VSolid * moduleS;
    G4LogicalVolume * moduleLV;

    //Materials
    G4Material* AirMaterial;
    G4Material* AirKillerMaterial;
    G4Material* WoMaterial;
    G4Material* InterfaceContainerMaterial;
    G4Material* InterfaceMaterial;
    G4Material* GapAbsToInterfaceMaterial;
    G4Material* GapInterfaceToReadoutMaterial;
    G4Material* PLEXMaterial;
    G4Material* PVCMaterial;


    // mandatory volumes
    // prepare pointers
    G4VSolid        *Interface_Positive_Z_S;
    G4LogicalVolume *Interface_Positive_Z_LV;
    G4PVPlacement   *Interface_Positive_Z_PV;

    G4VSolid        *Gap_Abs_Interface_Positive_Z_S;
    G4LogicalVolume *Gap_Abs_Interface_Positive_Z_LV;
    G4PVPlacement   *Gap_Abs_Interface_Positive_Z_PV;

    G4VSolid        *Gap_Interface_Readout_Positive_Z_S;
    G4LogicalVolume *Gap_Interface_Readout_Positive_Z_LV;
    G4PVPlacement   *Gap_Interface_Readout_Positive_Z_PV;

    G4VSolid        *Interface_Negative_Z_S;
    G4LogicalVolume *Interface_Negative_Z_LV;
    G4PVPlacement   *Interface_Negative_Z_PV;

    G4VSolid        *Gap_Abs_Interface_Negative_Z_S;
    G4LogicalVolume *Gap_Abs_Interface_Negative_Z_LV;
    G4PVPlacement   *Gap_Abs_Interface_Negative_Z_PV;

    G4VSolid        *Gap_Interface_Readout_Negative_Z_S;
    G4LogicalVolume *Gap_Interface_Readout_Negative_Z_LV;
    G4PVPlacement   *Gap_Interface_Readout_Negative_Z_PV;

    G4VSolid        *Readout_Positive_Z_S;
    G4LogicalVolume *Readout_Positive_Z_LV;
    G4PVPlacement   *Readout_Positive_Z_PV;

    G4VSolid        *Readout_Negative_Z_S;
    G4LogicalVolume *Readout_Negative_Z_LV;
    G4PVPlacement   *Readout_Negative_Z_PV;

    G4VSolid        *Absorber_S;
    G4LogicalVolume *Absorber_LV;
    G4PVPlacement   *Absorber_PV;

    G4VSolid        **abs_section_S;
    G4LogicalVolume **abs_section_LV;
    G4PVPlacement   **abs_section_PV;

    G4VSolid        *esr_S;
    G4LogicalVolume *esr_LV;
    G4PVPlacement   *esr_PV;

    std::vector<G4VSolid*>        LAPPD_S;
    std::vector<G4LogicalVolume*> LAPPD_LV;
    std::vector<G4PVPlacement*>   LAPPD_PV;
    std::vector<G4Material*>      LAPPD_materials_vec;

    G4double moduleZshift;
    G4double abs_interface_extraGap;



};




#endif /*Spacal_h*/
