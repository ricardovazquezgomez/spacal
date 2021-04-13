// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch
//

#ifndef SteppingAction_H
#define SteppingAction_H 1

#include "G4UserSteppingAction.hh"
#include "G4SteppingManager.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "EventAction.hh"
#include "G4ProcessManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4UnitsTable.hh"
#include "ConfigFile.hh"

#include "CreateTree.hh"
#include "DetectorConstruction.hh"
#include "TrackInformation.hh"
#include "MyMaterials.hh"
#include "LedFiberTiming.hh"

#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TRandom3.h"
// #include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
// #include "G4VPVParameterisation.hh"



class SteppingAction : public G4UserSteppingAction
{
public:

  SteppingAction  (const string& configFileName) ;

  SteppingAction(PrimaryGeneratorAction* gen_action);
  ~SteppingAction();
  float        propagate(float zstart, float zend, float refractive=1.8);
  virtual void UserSteppingAction(const G4Step*);


private:
  // DetectorConstruction* fDetectorConstruction;
  PrimaryGeneratorAction *fPrimaryGeneratorAction;

  // G4int propagateScintillation;
  // G4int propagateCerenkov;
  // G4double absorber_x ;     //size of rectangle containing fibres inside module
  // G4double absorber_y ;
  // G4double module_xy;
  // G4double module_yx;

  std::string logging_volume;
  std::string pre_volume;

  std::size_t modBegin;
  std::string modString ;

  bool saveAll            ;
  bool saveTree           ;
  bool saveShower         ;
  bool saveStructure      ;
  bool savePrimaries      ;
  bool savePhotons        ;
  bool savePhotonAbsPoint ;
  bool saveSummary ;
  bool saveLAPPD          ;
  // bool opticalCalibrationRun;
};

#endif
