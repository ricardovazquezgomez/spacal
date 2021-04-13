// Parameters
// Marco Pizzichemi 20.04.2020 marco.pizzichemi@cern.ch

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

#ifndef Parameters_h
#define Parameters_h 1

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
#include "G4ParticleDefinition.hh"


#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
// #include "Absorber.hh"

#include "Data.hh"

// class for reading multiple input parameters
class Parameters
{


public:

  // ctor
  Parameters();

  //! dtor
  ~Parameters ();

  int calorimeter_method;

  

  // a struct to hold relevant data on the different module types generated in the simulation
  struct Module_t
  {
    std::string fileName; // config file name for this module type
    std::vector<int> crystalMaterialList;
    int ecal_type;
  };

  struct Point
  {
    G4double x, y , z;
    Point(G4double a, G4double b, int c) { this->x = a; this->y = b; this->z = c;}
  };

  std::vector<Point> coordinates;
  G4double module_block_size_x;
  G4double module_block_size_y;
  G4double xMapMin;
  G4double xMapMax;
  G4double yMapMin;
  G4double yMapMax;
  std::vector<int> ecal_type_number;
  std::vector<G4double> ecal_type_xMapMin;
  std::vector<G4double> ecal_type_xMapMax;
  std::vector<G4double> ecal_type_yMapMin;
  std::vector<G4double> ecal_type_yMapMax;


  // structures for main config file and for module config files
  Data_t main_config;
  std::vector<Data_t> module_config;

 
  
  // vector of Module_t
  std::vector<Module_t> module;

  // read main config file
  void ReadMainConfig(std::string fileName);

  // read a module config file
  void ReadModuleConfig(std::string fileName,int ecal_type,std::string container_volume);

  // routine that actually reads config files
  struct Data_t FillStructure(std::string fileName);

  // routine to extract positional info from ROOT file
  void ScanMap(TH2D *histo2d);
  TH2D *mapEcal;

  // print a config structure data
  void PrintConfig(Data_t data);

  // primary particle definition
  std::string particleDefinition;

  // global material list
  std::map<int, G4Material*> CrystalMaterialMap;
  std::vector<int>           CrystalMaterialList;
  void               AddCrystalMaterial(int num, G4Material *aMaterial);
  G4Material*        GetCrystalMaterial(int num);
  std::vector<int>   GetCrystalMaterialList(){return CrystalMaterialList;};
  void               AddTimeHisto(TH1F* aHisto){tHisto.push_back(aHisto);};
  void               AddEnergyHisto(TH1F* aHisto){enHisto.push_back(aHisto);};
  void               AddAbsHisto(TH1F* aHisto){absHisto.push_back(aHisto);};
  std::vector<TH1F*> enHisto;
  std::vector<TH1F*> tHisto;
  std::vector<TH1F*> absHisto;

  // these members are static, which means that
  // no matter how many times this class will be called
  // they will be set only once. Combined with the setting
  // of fInstance in the constructor, and the return if
  // fInstance is already defined, this ensures that parameters
  // will be created only once, and makes it accessible from everywhere
  // in the program, provided that this header is included. Sounds
  // like a dirty trick to me, but it works so let's use it
  static Parameters* Instance () { return fInstance ; } ;
  static Parameters* fInstance ;

  void WriteParameters(TFile *outfile);
  void doWriteParameters(TFile *outfile, Data_t parameters);
  void CreateParamTree();

  TTree* paramTTree;
  // struct to save params in output file 
  // just a service struct, this struct will not be saved 
  // like it is in the ttree, but its values will 
  Data_t paramStruct;
  
};




#endif /*Parameters_h*/
