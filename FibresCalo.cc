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
// $Id: exampleN06.cc,v 1.18 2010-10-23 19:33:55 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
// Description: Test of Continuous Process G4Cerenkov
//              and RestDiscrete Process G4Scintillation
//              -- Generation Cerenkov Photons --
//              -- Generation Scintillation Photons --
//              -- Transport of optical Photons --
// Version:     5.0
// Created:     1996-04-30
// Author:      Juliet Armstrong
// mail:        gum@triumf.ca
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <iostream>
#include <string>
#include <vector>
#include <ctime>

#include "TString.h"
#include "TRandom3.h"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4PhysListFactory.hh"
#include "G4EmUserPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4VModularPhysicsList.hh"
#include "G4GeometryTolerance.hh"
#include "QGSP_BERT.hh"
#include "G4GeometryManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"
#include "EventAction.hh"
#include "SteppingVerbose.hh"
#include "CreateTree.hh"
#include "StackingAction.hh"
#include "Parameters.hh"

// #ifdef G4VIS_USE
#include "G4VisExecutive.hh"
// #endif

// #ifdef G4UI_USE
#include "G4UIExecutive.hh"
// #endif


# include <unistd.h>
# include <sys/time.h>
# include <sys/resource.h>



using namespace CLHEP;



long int CreateSeed();



int main(int argc,char** argv)
{
  //gInterpreter -> GenerateDictionary("vector<float>","vector");
  //gInterpreter -> GenerateDictionary("map<int,vector<float> >","map,vector");


  if (argc < 2 || argc > 6)
  {
    cout << "Syntax for exec with gps, seed, and flux: FibresCalo <configuration file> <output file without extension> <gps file> <random seed> <flux file>" << endl;
    cout << "Syntax for exec with gps and seed: FibresCalo <configuration file> <output file without extension> <gps file> <random seed>" << endl;
    cout << "Syntax for exec with gps: FibresCalo <configuration file> <output file without extension> <gps file> " << endl;
    cout << "Syntax for exec: FibresCalo <configuration file> <output file>" << endl;
    cout << "Syntax for viz: FibresCalo <configuration file>" << endl;
    return 0;
  }

  string file;
  string filename;
  TFile* outfile = NULL;
  // create a fake outfile
  outfile = new TFile(".temp.root","RECREATE");
  bool fakeOutput = true;
  bool overriding_flux_file = false;
  // float x,y,z;

  if(argc == 4 || argc == 5 || argc == 6)
  {
    cout << "Starting special mode..." << endl;
    file = argv[2];
    filename = file + ".root";
    G4cout << "Writing data to file '" << filename << "' ..." << G4endl;

    if(argc == 5 || argc == 6)
    {
      cout << "Random seed fix by user at " << argv[4] << endl;
    }

    if(argc == 6)
    {
      overriding_flux_file = true;
      cout << "Flux source file fixed by user to " << argv[5] << endl;
    }

    fakeOutput = false;
    outfile = new TFile((TString)filename,"RECREATE");
    outfile -> cd();
  }

  if(argc == 3)
  {
    cout << "Starting exec mode..." << endl;
    file = argv[2];
    filename = file + ".root";
    G4cout << "Writing data to file '" << filename << "' ..." << G4endl;

    fakeOutput = false;
    outfile = new TFile((TString)filename,"RECREATE");
    outfile -> cd();
  }



  if (argc == 2)
  {
    cout<<"Starting viz mode..."<<endl;
  }
  
  // G4double worldLength = std::max({expHall_x, expHall_y, expHall_z});
  
  
  // G4double worldLength = 1000 * m;
  // G4GeometryManager :: GetInstance ()
  //                      -> SetWorldMaximumExtent ( worldLength );
  // G4cout << " Computed tolerance = "
  //        << G4GeometryTolerance :: GetInstance ()
  //                                  -> GetSurfaceTolerance ()/ nm
  //        << " nm " << G4endl ;
  
  
  // std::cout << "G4GeometryTolerance::GetInstance()->GetSurfaceTolerance() = " <<  G4GeometryTolerance::GetInstance()->GetSurfaceTolerance() << std::endl;
  // std::cout << "G4GeometryTolerance::GetInstance()->GetAngularTolerance() = " <<  G4GeometryTolerance::GetInstance()->GetAngularTolerance() << std::endl;
  // std::cout << "G4GeometryTolerance::GetInstance()->GetRadialTolerance()  = " <<  G4GeometryTolerance::GetInstance()->GetRadialTolerance() << std::endl;

  // std::cout << G4GeometryTolerance::GetInstance()->GetSurfaceTolerance() << std::endl;
  // std::cout << G4GeometryTolerance::GetInstance()->GetAngularTolerance() << std::endl;
  // std::cout << G4GeometryTolerance::GetInstance()->GetRadialTolerance() << std::endl;

  // G4GeometryTolerance::GetInstance()->SetSurfaceTolerance(1e-5) ;

  
  
  cout<<"=====>   C O N F I G U R A T I O N   <====\n"<<endl;

  // GET INFO ON PC AND FILES
  char hostname[1024];
  hostname[1023] = '\0';
  gethostname(hostname, 1023);
  std::string HostNameString(hostname);
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::string PWDstring(cwd);
  //PUT THE ENTIRE CONFIG FILE IN A STRINGSTREAM, TO SAVE IT LATER IN THE OUTPUT FILE
  std::vector<std::string> sConfigFiles;
  std::string ConfigFileName ;
  ConfigFileName = argv[1];
  std::stringstream streamConfigFile;
  std::string strConfig;
  std::ifstream inConfig;
  inConfig.open (ConfigFileName.c_str(), std::ifstream::in);
  while (std::getline(inConfig, strConfig))
  {
    streamConfigFile << strConfig << std::endl;
  }
  sConfigFiles.push_back(streamConfigFile.str());
  inConfig.close();

  G4cout << "Main configuration file: '" << argv[1] << "'" << G4endl;
  ConfigFile config(argv[1]);



  // new config reading
  Parameters* parameters = new Parameters();
  // read main config
  Parameters::Instance()->ReadMainConfig(std::string(argv[1]));

  // override the flux input file, if passed by command line arg 
  // and force use_gps to false
  if(overriding_flux_file)
  {
    Parameters::Instance()->main_config.input_filename = argv[5];
    Parameters::Instance()->main_config.use_gps = false;
  }

  // now check if main config specifies other config files, and in that case read them
  // and prepare to save them in output file
  if(Parameters::Instance()->main_config.module_config_list.size() > 0)
  {
    Parameters::Instance()->calorimeter_method = 2;
    //
    for(int iC = 0 ; iC < Parameters::Instance()->main_config.module_config_list.size() ; iC++ )
    {
      G4cout << "Module configuration file  '" << Parameters::Instance()->main_config.module_config_list[iC] << "'" << G4endl;
      Parameters::Instance()->ReadModuleConfig(Parameters::Instance()->main_config.module_config_list[iC],
                                               Parameters::Instance()->main_config.ecal_type[iC],
                                               Parameters::Instance()->main_config.container_volume[iC]);
      // prepare string for saving module config files into output file
      std::string ModuleConfigFileName ;
      ModuleConfigFileName = Parameters::Instance()->main_config.module_config_list[iC];
      std::stringstream streamModuleConfigFile;
      std::string strModuleConfig;
      std::ifstream inModuleConfig;
      inModuleConfig.open (ModuleConfigFileName.c_str(), std::ifstream::in);
      while (std::getline(inModuleConfig, strModuleConfig))
      {
        streamModuleConfigFile << strModuleConfig << std::endl;
      }
      sConfigFiles.push_back(streamModuleConfigFile.str());
      inModuleConfig.close();
    }
    // read ecal map
    TFile *_file0 = TFile::Open(Parameters::Instance()->main_config.ecal_map_file.c_str());
    TH2D *mapEcal = (TH2D*) _file0->Get(Parameters::Instance()->main_config.ecal_map_histo.c_str());
    G4cout << Parameters::Instance()->main_config.ecal_map_histo.c_str() << " " << mapEcal << G4endl;
    Parameters::Instance()->ScanMap(mapEcal);

    Parameters::Instance()->mapEcal = mapEcal;
    // set the directory back to output file, or things are going to get crazy
    // thanks to one of the nicest features of ROOT: always create objects in pwd
    outfile->cd();
  }
  else
  {
    Parameters::Instance()->calorimeter_method = 1;
  }

  // Seed the random number generator manually
  //
  G4long myseed = Parameters::Instance()->main_config.seed;
  if(argc == 5)
  {
    myseed = atoi(argv[4]);
  }

  if( myseed == -1 )
  {
    G4cout << "Creating random seed..." << G4endl;
    myseed = CreateSeed();
  }
  std::stringstream streamSeed;
  streamSeed << myseed;

  G4cout << "Random seed : " << myseed << G4endl;
  CLHEP::HepRandom::setTheSeed(myseed);

  std::vector<float> attLengths;
  // config.readIntoVect(attLengths,"attLengths");

  // RUN OPTIONS
  // Parameters::Instance()->simulation_flags.primaries               = config.read<int>("primaries",1);
  // Parameters::Instance()->simulation_flags.saveAll                 = config.read<bool>("saveAll",1);
  // Parameters::Instance()->simulation_flags.saveTree                = config.read<bool>("saveTree",0);
  // Parameters::Instance()->simulation_flags.saveShower              = config.read<bool>("saveShower",0);
  // Parameters::Instance()->simulation_flags.saveStructure           = config.read<bool>("saveStructure",1);
  // Parameters::Instance()->simulation_flags.savePrimaries           = config.read<bool>("savePrimaries",0);
  // Parameters::Instance()->simulation_flags.savePhotons             = config.read<bool>("savePhotons",0);
  // Parameters::Instance()->simulation_flags.savePhotonAbsPoint      = config.read<bool>("savePhotonAbsPoint",0);
  // Parameters::Instance()->simulation_flags.saveSimple              = config.read<bool>("saveSimple",0);
  // Parameters::Instance()->simulation_flags.saveSummary             = config.read<bool>("saveSummary",0);
  // Parameters::Instance()->simulation_flags.savePhotonGen           = config.read<bool>("savePhotonGen",0);
  // Parameters::Instance()->simulation_flags.saveEnergyPerModule     = config.read<bool>("saveEnergyPerModule",0);
  // Parameters::Instance()->simulation_flags.switchOnScintillation   = config.read<int> ("switchOnScintillation");
  // Parameters::Instance()->simulation_flags.switchOnCerenkov        = config.read<int> ("switchOnCerenkov");
  // Parameters::Instance()->simulation_flags.propagateScintillation  = config.read<int> ("propagateScintillation");
  // Parameters::Instance()->simulation_flags.propagateCerenkov       = config.read<int> ("propagateCerenkov");

  // read specific runs flags
  // Parameters::Instance()->simulationType          = config.read<int>("simulationType",-1);
  // simulation types:
  // not specified        = -1
  // only energy depo     = 0
  // full ray tracing     = 1
  // hybrid               = 2
  // optical calibration  = 3
  // if user do not specify in config file, the simulation type is "free"
  // i.e. the save flags above are left untouched.
  // otherwise, those flags will be overwritten according to the simulation type
  if(Parameters::Instance()->main_config.simulationType == -1)
  {
    cout<<"=====>   S I M U L A T I O N      T Y P E   <===="<<endl;
    cout<<"=====>           FREE SIMULATION            <===="<<endl;
    cout<<"=====>  ----------------------------------  <===="<<endl;
    // do nothing
  }
  else
  {
    if(Parameters::Instance()->main_config.simulationType == 0)
    {
      cout<<"=====>   S I M U L A T I O N      T Y P E   <===="<<endl;
      cout<<"=====>        ONLY ENERGY DEPOSITION        <===="<<endl;
      cout<<"=====>  ----------------------------------  <===="<<endl;
      Parameters::Instance()->main_config.saveAll                = 0;
      Parameters::Instance()->main_config.saveTree               = 0;
      Parameters::Instance()->main_config.saveShower             = 1;
      Parameters::Instance()->main_config.saveStructure          = 1;
      Parameters::Instance()->main_config.savePrimaries          = 1;
      Parameters::Instance()->main_config.savePhotons            = 0;
      Parameters::Instance()->main_config.savePhotonAbsPoint     = 0;
      Parameters::Instance()->main_config.saveSimple             = 0;
      Parameters::Instance()->main_config.saveSummary            = 0;
      Parameters::Instance()->main_config.savePhotonGen          = 0;
      Parameters::Instance()->main_config.switchOnScintillation  = 0;
      Parameters::Instance()->main_config.switchOnCerenkov       = 0;
      Parameters::Instance()->main_config.propagateScintillation = 0;
      Parameters::Instance()->main_config.propagateCerenkov      = 0;
      Parameters::Instance()->main_config.primaries              = 1;
    }
    if(Parameters::Instance()->main_config.simulationType == 1)
    {
      cout<<"=====>   S I M U L A T I O N      T Y P E   <===="<<endl;
      cout<<"=====>           FULL RAY TRACING           <===="<<endl;
      cout<<"=====>  ----------------------------------  <===="<<endl;
      Parameters::Instance()->main_config.saveAll                = 0;
      Parameters::Instance()->main_config.saveTree               = 0;
      Parameters::Instance()->main_config.saveShower             = 1;
      Parameters::Instance()->main_config.saveStructure          = 1;
      Parameters::Instance()->main_config.savePrimaries          = 1;
      Parameters::Instance()->main_config.savePhotons            = 1;
      Parameters::Instance()->main_config.savePhotonAbsPoint     = 0;
      Parameters::Instance()->main_config.saveSimple             = 0;
      Parameters::Instance()->main_config.saveSummary            = 0;
      Parameters::Instance()->main_config.savePhotonGen          = 0;  // in real full ray tracing this is not useful, but it would be if one wants to compare also to the "photon gen" approach  (which anyway works bad, as it does not take the cerenkov into account correcty)
      Parameters::Instance()->main_config.switchOnScintillation  = 1;
      Parameters::Instance()->main_config.switchOnCerenkov       = 1;
      Parameters::Instance()->main_config.propagateScintillation = 1;
      Parameters::Instance()->main_config.propagateCerenkov      = 1;
      Parameters::Instance()->main_config.primaries              = 1;
    }
    if(Parameters::Instance()->main_config.simulationType == 2)
    {
      cout<<"=====>   S I M U L A T I O N      T Y P E   <===="<<endl;
      cout<<"=====>          HYBRID SIMULATION           <===="<<endl;
      cout<<"=====>  ----------------------------------  <===="<<endl;
      Parameters::Instance()->main_config.saveAll                = 0;
      Parameters::Instance()->main_config.saveTree               = 0;
      Parameters::Instance()->main_config.saveShower             = 1;
      Parameters::Instance()->main_config.saveStructure          = 1;
      Parameters::Instance()->main_config.savePrimaries          = 1;
      Parameters::Instance()->main_config.savePhotons            = 1;
      Parameters::Instance()->main_config.savePhotonAbsPoint     = 0;
      Parameters::Instance()->main_config.saveSimple             = 0;
      Parameters::Instance()->main_config.saveSummary            = 0;
      Parameters::Instance()->main_config.savePhotonGen          = 0;
      Parameters::Instance()->main_config.switchOnScintillation  = 0;
      Parameters::Instance()->main_config.switchOnCerenkov       = 1;
      Parameters::Instance()->main_config.propagateScintillation = 0;
      Parameters::Instance()->main_config.propagateCerenkov      = 1;
      Parameters::Instance()->main_config.primaries              = 1;
    }
    if(Parameters::Instance()->main_config.simulationType == 3)
    {
      cout<<"=====>   S I M U L A T I O N      T Y P E   <===="<<endl;
      cout<<"=====>         OPTICAL CALIBRATION          <===="<<endl;
      cout<<"=====>  ----------------------------------  <===="<<endl;
      Parameters::Instance()->main_config.saveAll                = 0;
      Parameters::Instance()->main_config.saveTree               = 0;
      Parameters::Instance()->main_config.saveShower             = 0;
      Parameters::Instance()->main_config.saveStructure          = 1;
      Parameters::Instance()->main_config.savePrimaries          = 0;
      Parameters::Instance()->main_config.savePhotons            = 1;
      Parameters::Instance()->main_config.savePhotonAbsPoint     = 0;
      Parameters::Instance()->main_config.saveSimple             = 0;
      Parameters::Instance()->main_config.saveSummary            = 0;
      Parameters::Instance()->main_config.savePhotonGen          = 0;
      Parameters::Instance()->main_config.switchOnScintillation  = 1;
      Parameters::Instance()->main_config.switchOnCerenkov       = 1;
      Parameters::Instance()->main_config.propagateScintillation = 1;
      Parameters::Instance()->main_config.propagateCerenkov      = 1;
      // different set for primaries. by default they should be 100k. if the user has put 1 or any value lower than 100k, we assume it's a mistake and correct it.
      // otherwise, we use the user choice
      if(Parameters::Instance()->main_config.primaries != 1)
      {
        cout << "WARNING: you chose an Optical Calibration run, but you set the primaries to " << Parameters::Instance()->main_config.primaries << std::endl;
        cout << "The number of primaries in Optical Calibration mode had to be set to 1, because the actual number will be controlled by GPS file" << std::endl;
        Parameters::Instance()->main_config.primaries = 1;
      }
      else
      {
        // leave the choice of the user untouched
      }

      // if(Parameters::Instance()->main_config.optical_calibration_module_type < 0)
      // {
      //   cout << "ERROR: Invalid choice of module type! This is set to be an optical calibration run, but the optical_calibration_module_type " << std::endl;
      //   cout << "       has not be set in the config file! Aborting... " << std::endl;
      //   return -1;
      // }
      // else // set in CreateTree
      // {
      //   CreateTree::Instance()->optType.push_back(Parameters::Instance()->main_config.optical_calibration_module_type);
      // }
    }
  }

  std::stringstream streamPrimaries;
  streamPrimaries << Parameters::Instance()->main_config.primaries;

  // int opticalCalibrationRun = 0;
  // if(Parameters::Instance()->main_config.simulationType == 3)
  // {
    // opticalCalibrationRun = 1;
  // }

  CreateTree* mytree = new CreateTree ("tree") ;

  // if(Parameters::Instance()->main_config.simulationType == 3)
  // {
  //   CreateTree::Instance()->SaveOpticalRunVector = true;
  // }
  // std::string logging_volume = config.read<std::string>("logging_volume","PMT");
  // G4cout << "Logging volume contains : " << logging_volume << G4endl;

  // std::string pre_volume = config.read<std::string>("pre_volume","");
  // if(pre_volume == "")
  // {
  //   G4cout << "Not checking for pre volume " << G4endl;
  // }
  // else
  // {
  //   G4cout << "Pre volume contains : " << pre_volume << G4endl;
  // }


  // // YANXI MOD
  // int modules_nx      = config.read<int>("modules_nx",1);
  // int modules_ny      = config.read<int>("modules_ny",1);
  // int pipe_modules_nx = config.read<int>("pipe_modules_nx",0);
  // int pipe_modules_ny = config.read<int>("pipe_modules_ny",0);
  // CreateTree::Instance()->module_array_x = modules_nx;
  // CreateTree::Instance()->module_array_y = modules_ny;
  // CreateTree::Instance()->pipe_modules_nx = pipe_modules_nx;
  // CreateTree::Instance()->pipe_modules_ny = pipe_modules_ny;
  //
  // CreateTree::Instance()->module_total = modules_nx*modules_ny;
  // G4cout<<"Number of modules "<< modules_nx<<"x"<<modules_ny<<G4endl;
  // // these numbers define the 4 or 1 regions of the calorimeter. We need to save this info
  // // as usual, we use CreateTree as storing place
  // if((pipe_modules_nx == 0) && (pipe_modules_ny == 0))
  // {
  //   CreateTree::Instance()->pipe = 0;
  // }
  // else
  // {
  //   CreateTree::Instance()->pipe = 1;
  // }
  // CreateTree::Instance()->depositedEnergyByModuleFront = new float[CreateTree::Instance()->module_total];
  // CreateTree::Instance()->depositedEnergyByModuleRear = new float[CreateTree::Instance()->module_total];
  // CreateTree::Instance()->depositionTimingByModuleFront = new float[CreateTree::Instance()->module_total];
  // CreateTree::Instance()->depositionTimingByModuleRear = new float[CreateTree::Instance()->module_total];
  // if(Parameters::Instance()->main_config.saveAll || Parameters::Instance()->main_config.saveEnergyPerModule)
  // {
  //   CreateTree::Instance()->CreateTreeEnergyDepositByModule();
  // }

  // END OF YANXI MOD

  // User Verbose output class
  //
  G4VSteppingVerbose* verbosity = new SteppingVerbose;
  G4VSteppingVerbose::SetInstance(verbosity);

  // Run manager
  //
  G4RunManager* runManager = new G4RunManager();
  runManager->SetVerboseLevel(0);
  runManager->SetPrintProgress(-1);


  //Physics list defined using PhysListFactory
  //
  std::string physName("");

  G4PhysListFactory factory;
  const std::vector<G4String>& names = factory.AvailablePhysLists();
  for(unsigned n = 0; n != names.size(); n++)
    G4cout << "PhysicsList: " << names[n] << G4endl;

  if( physName == "")
  {
    char* path = getenv("PHYSLIST");
    if( path ) physName = G4String(path);
  }

  if ( physName == "" || factory.IsReferencePhysList(physName))
  {
    physName = "FTFP_BERT";
    // physName = "FTFP_BERT_EMZ"; // high accuracy but for low energy
    // physName = "FTFP_BERT_EMV"; // less precise, but supposed to be faster
                                // it might be the one used by CMS
  }



  std::cout << "Using physics list: " << physName << std::endl;


  // UserInitialization classes - mandatory
  //

  G4cout << ">>> Define physics list::begin <<<" << G4endl;
  double defaultCutValue = Parameters::Instance()->main_config.defaultCut;
  G4double defaultCut = defaultCutValue * mm;
  G4cout << "> Default production cut set to " << defaultCut << " mm" << G4endl;
  G4VModularPhysicsList* physics = factory.GetReferencePhysList(physName);
  physics->SetVerboseLevel(0);
  physics->SetDefaultCutValue(defaultCut);
  physics->RegisterPhysics(new G4EmUserPhysics(Parameters::Instance()->main_config.switchOnScintillation,Parameters::Instance()->main_config.switchOnCerenkov));

  runManager-> SetUserInitialization(physics);
  G4cout << ">>> Define physics list::end <<<" << G4endl;

  G4cout << ">>> Define DetectorConstruction::begin <<<" << G4endl;
  DetectorConstruction* detector = new DetectorConstruction(argv[1]);
  runManager-> SetUserInitialization(detector);
  G4cout << ">>> Define DetectorConstruction::end <<<" << G4endl;

  G4cout << ">>> Define PrimaryGeneratorAction::begin <<<" << G4endl;

  // G4ThreeVector posCentre(
  //     0. * m,
  //     0. * m,
  //     -1. * (detector->GetModule_z () / m) / 2. * m
  //   ) ;


  PrimaryGeneratorAction* gen_action = new PrimaryGeneratorAction(config);
  // G4bool use_gps;
  // config.readInto(use_gps,"useGPS",false);
  // G4bool use_flux_fromTree= !use_gps;
  runManager->SetUserAction(gen_action);
  G4cout << ">>> Define PrimaryGeneratorAction::end <<<" << G4endl;

  // UserAction classes
  //

  G4cout << ">>> Define RunAction::begin <<<" << G4endl;
  G4UserRunAction* run_action = new RunAction;
  runManager->SetUserAction(run_action);
  G4cout << ">>> Define RunAction::end <<<" << G4endl;

  G4cout << ">>> Define EventAction::begin <<<" << G4endl;
  G4UserEventAction* event_action = new EventAction();
  runManager->SetUserAction(event_action);
  G4cout << ">>> Define EventAction::end <<<" << G4endl;

  G4cout << ">>> Define TrackingAction::begin <<<" << G4endl;
  TrackingAction* tracking_action = new TrackingAction;
  runManager->SetUserAction(tracking_action);
  G4cout << ">>> Define TrackingAction::end <<<" << G4endl;

  //
  G4cout << ">>> Define StackingAction::begin <<<" << G4endl;
  G4UserStackingAction* stacking_action = new StackingAction();
  runManager->SetUserAction(stacking_action);
  G4cout << ">>> Define StackingAction::end <<<" << G4endl;


  G4cout << ">>> Define SteppingAction::begin <<<" << G4endl;

  SteppingAction* stepping_action = new SteppingAction(gen_action);
  runManager->SetUserAction(stepping_action);
  G4cout << ">>> Define SteppingAction::end <<<" << G4endl;

  string gps_instructions_file = "" ;

  if (argc == 2)   // Define UI session for interactive mode
  {
    // Initialize G4 kernel
    //
    runManager -> Initialize();
    G4UIExecutive* ui = 0;
    ui = new G4UIExecutive(argc, argv);
    // Initialize visualization
    //
    G4VisManager* visManager = new G4VisExecutive;
    // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
    // G4VisManager* visManager = new G4VisExecutive("Quiet");
    visManager->Initialize();
    // Get the pointer to the User Interface manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
  }
  else if(argc == 3)
  {
    runManager -> Initialize();
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    gps_instructions_file = Parameters::Instance()->main_config.gps_instructions_file;
    // config.readInto (gps_instructions_file, "gps_instructions_file") ;
    UImanager -> ApplyCommand("/control/execute " + gps_instructions_file);
  }
  else if(argc == 4 || argc == 5 || argc == 6)
  {
    runManager -> Initialize();
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    gps_instructions_file = argv[3];
    G4cout << gps_instructions_file << G4endl;
    // config.readInto (gps_instructions_file, "gps_instructions_file") ;
    UImanager -> ApplyCommand("/control/execute " + gps_instructions_file);
  }
  else
  {
    G4cout << "ERROR: invalid number of command line arguments" << G4endl;
    return 0;
  }

  // delete runManager;
  // delete verbosity;

  // save gps file into string
  //PUT THE ENTIRE CONFIG FILE IN A STRINGSTREAM, TO SAVE IT LATER IN THE OUTPUT FILE
  std::string GpsFileName ;
  GpsFileName = gps_instructions_file;
  std::stringstream streamGpsFile;
  std::string strGps;
  std::ifstream inGps;
  inGps.open (GpsFileName.c_str(), std::ifstream::in);
  while (std::getline(inGps, strGps))
  {
    streamGpsFile << strGps << std::endl;
  }
  inGps.close();

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not
  // be deleted in the main() program !

  delete runManager;
  delete verbosity;

  if(fakeOutput)
  {
    if( remove( ".temp.root" ) != 0 )
      G4cout << "Error deleting file .temp.root" << G4endl;
    else
      G4cout << ".temp.root file successfully deleted" << G4endl;

  }


  if(argc == 3 || argc == 4 || argc == 5 || argc == 6)
  {


    G4cout << "Writing tree to file " << filename << " ..." << G4endl;



    mytree->Write (outfile) ;
    // std::vector<int> modules_n;
    // std::vector<int> pipe_modules_n;
    // modules_n.push_back(modules_nx);
    // modules_n.push_back(modules_ny);
    // pipe_modules_n.push_back(pipe_modules_nx);
    // pipe_modules_n.push_back(pipe_modules_ny);
    // gDirectory->WriteObject(&modules_n, "modules_n");
    // gDirectory->WriteObject(&pipe_modules_n, "pipe_modules_n");
    // // int module_nx = config.read<float>("module_nx");
    // // int module_ny = config.read<float>("module_ny");
    // float absorber_size_x = config.read<float>("absorber_size_x");
    // float absorber_size_y = config.read<float>("absorber_size_y");
    // float absorber_size_z = config.read<float>("absorber_size_z");
    // std::vector<int> module_n;
    // std::vector<float> absorber_size;
    // absorber_size.push_back(absorber_size_x);
    // absorber_size.push_back(absorber_size_y);
    // absorber_size.push_back(absorber_size_z);
    // gDirectory->WriteObject(&absorber_size, "absorber_size");
    // std::vector<float> separator_z_position;
    // separator_z_position.push_back(CreateTree::Instance()->zSeparationAbsolutePosition);
    // gDirectory->WriteObject(&separator_z_position, "separator_z_position");

    // save type of module positioning

    // std::vector<int> calorimeter_method;
    // calorimeter_method.push_back(Parameters::Instance()->calorimeter_method);
    // gDirectory->WriteObject(&calorimeter_method, "calorimeter_method");
    
    Parameters::Instance()->WriteParameters(outfile);

    //save info on sim configuration
    TDirectory *configDir = outfile->mkdir("Configuration");
    configDir->cd();
    TNamed SeedNameD("Seed",streamSeed.str().c_str());
    TNamed HostNameD("Hostname",HostNameString.c_str());
    TNamed PWDNameD("PWD",PWDstring.c_str());
    for(unsigned int i = 0; i < sConfigFiles.size(); i++)
    {
      std::stringstream sConf;
      sConf << "ConfigFile";
      if(i != 0) sConf << i; // for back compatibility
      TNamed ConfigNameD(sConf.str().c_str(),sConfigFiles[i].c_str());
      ConfigNameD.Write();
    }
    TNamed GpsNameD("GpsFile",streamGpsFile.str().c_str());
    TNamed PrimariesNameD("Primaries",streamPrimaries.str().c_str());
    SeedNameD.Write();
    HostNameD.Write();
    PWDNameD.Write();
    GpsNameD.Write();
    PrimariesNameD.Write();
    configDir->cd("..");
    outfile->Close () ;
  }

  return 0;
}



long int CreateSeed()
{
  TRandom3 rangen;

  long int sec = time(0);
  G4cout << "Time : " << sec << G4endl;

  sec += getpid();
  G4cout << "PID  : " << getpid() << G4endl;

  FILE* fp = fopen ("/proc/uptime", "r");
  int upsecs = 0;
  if( fp != NULL )
  {
    char buf[BUFSIZ];
    char *b = fgets(buf,BUFSIZ,fp);
    if( b == buf )
    {
      /* The following sscanf must use the C locale.  */
      setlocale(LC_NUMERIC, "C");
      setlocale(LC_NUMERIC, "");
    }
    fclose(fp);
  }
  G4cout << "Upsecs: " << upsecs << G4endl;
  sec += upsecs;

  G4cout << "Seed for srand: " << sec << G4endl;
  srand(sec);
  rangen.SetSeed(rand());
  long int seed = round(1000000*rangen.Uniform());
  return seed;
}
