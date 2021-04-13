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
// $Id: SteppingVerbose.cc,v 1.4 2006-06-29 17:54:37 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingVerbose.hh"

#include "G4SteppingManager.hh"
#include "G4UnitsTable.hh"

#include "G4OpBoundaryProcess.hh"


std::string OPprocessName[] =  {  "Undefined",
                                  "Transmission",
                                  "FresnelRefraction",
                                  "FresnelReflection",
                                  "TotalInternalReflection",
                                  "LambertianReflection",
                                  "LobeReflection",
                                  "SpikeReflection",
                                  "BackScattering",
                                  "Absorption",
                                  "Detection",
                                  "NotAtBoundary",
                                  "SameMaterial",
                                  "StepTooSmall",
                                  "NoRINDEX",
                                  "PolishedLumirrorAirReflection",
                                  "PolishedLumirrorGlueReflection",
                                  "PolishedAirReflection",
                                  "PolishedTeflonAirReflection",
                                  "PolishedTiOAirReflection",
                                  "PolishedTyvekAirReflection",
                                  "PolishedVM2000AirReflection",
                                  "PolishedVM2000GlueReflection",
                                  "EtchedLumirrorAirReflection",
                                  "EtchedLumirrorGlueReflection",
                                  "EtchedAirReflection",
                                  "EtchedTeflonAirReflection",
                                  "EtchedTiOAirReflection",
                                  "EtchedTyvekAirReflection",
                                  "EtchedVM2000AirReflection",
                                  "EtchedVM2000GlueReflection",
                                  "GroundLumirrorAirReflection",
                                  "GroundLumirrorGlueReflection",
                                  "GroundAirReflection",
                                  "GroundTeflonAirReflection",
                                  "GroundTiOAirReflection",
                                  "GroundTyvekAirReflection",
                                  "GroundVM2000AirReflection",
                                  "GroundVM2000GlueReflection",
                                  "Dichroic"
                                  };

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingVerbose::SteppingVerbose()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingVerbose::~SteppingVerbose()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingVerbose::StepInfo()
{
  CopyState();

  G4int prec = G4cout.precision(3);

  if( verboseLevel >= 1 ){
    if( verboseLevel >= 4 ) VerboseTrack();
    if( verboseLevel >= 3 ){
      G4cout << G4endl;
      G4cout << std::setw(15) << "#Step#"
	           << std::setw(15) << "X"
	           << std::setw(15) << "Y"
	           << std::setw(15) << "Z"
	           << std::setw(15) << "KineE"
	           << std::setw(15) << "dEStep"
	           << std::setw(15) << "StepLeng"
	           << std::setw(15) << "TrakLeng"
             << std::setw(15) << "BoundSt"
	           << std::setw(25) << "Volume"
	           << std::setw(15) << "Process"   << G4endl;
    }

    G4cout << std::setw(15) << fTrack->GetCurrentStepNumber()
           << std::setw(15) << G4BestUnit(fTrack->GetPosition().x(),"Length")
           << std::setw(15) << G4BestUnit(fTrack->GetPosition().y(),"Length")
           << std::setw(15) << G4BestUnit(fTrack->GetPosition().z(),"Length")
           << std::setw(15) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
           << std::setw(15) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
           << std::setw(15) << G4BestUnit(fStep->GetStepLength(),"Length")
           << std::setw(15) << G4BestUnit(fTrack->GetTrackLength(),"Length");
    // boundary process
    G4OpBoundaryProcessStatus boundaryStatus=Undefined;
    static G4ThreadLocal G4OpBoundaryProcess* boundary=NULL;

    //find the boundary process only once
    if(!boundary)
    {
      G4ProcessManager* pm = fTrack->GetDefinition()->GetProcessManager();
      G4int nprocesses = pm->GetProcessListLength();
      G4ProcessVector* pv = pm->GetProcessList();
      G4int i;
      for( i=0;i<nprocesses;i++)
      {
        if((*pv)[i]->GetProcessName()=="OpBoundary")
        {
          boundary = (G4OpBoundaryProcess*)(*pv)[i];
          break;
        }
      }
    }
    if(boundary)
    {
      boundaryStatus = boundary->GetStatus();
      G4cout << std::setw(15) << OPprocessName[boundaryStatus];
    }
    else
    {
      G4cout << std::setw(15) << "Not Boundary";
    }

    G4StepPoint * thePrePoint  = fStep->GetPreStepPoint () ;
    G4StepPoint * thePostPoint = fStep->GetPostStepPoint () ;
    // const G4ThreeVector & thePrePosition  = thePrePoint->GetPosition () ;
    // G4VPhysicalVolume * thePrePV  = thePrePoint->GetPhysicalVolume () ;
    // G4VPhysicalVolume * thePostPV = thePostPoint->GetPhysicalVolume () ;
    G4String thePrePVName  = "" ; if ( thePrePoint->GetPhysicalVolume () )  thePrePVName  = thePrePoint->GetPhysicalVolume ()  -> GetName () ;
    G4String thePostPVName = "" ; if ( thePostPoint->GetPhysicalVolume () ) thePostPVName = thePostPoint->GetPhysicalVolume () -> GetName () ;
    G4cout << std::setw(40) << thePrePVName << std::setw(40) << thePostPVName ;


    // if( fStepStatus != fWorldBoundary){
    if( fTrack->GetNextVolume() != 0 ) {
      G4cout << std::setw(25) << fTrack->GetVolume()->GetName();
    } else {
      G4cout << std::setw(15) << "OutOfWorld";
    }

    if(fStep->GetPostStepPoint()->GetProcessDefinedStep() != 0){
      G4cout << "  "
             << std::setw(15)
	     << fStep->GetPostStepPoint()->GetProcessDefinedStep()
	                                 ->GetProcessName();
    } else {
      G4cout << "   UserLimit";
    }

    G4cout << G4endl;

    if( verboseLevel == 2 ){
      G4int tN2ndariesTot = fN2ndariesAtRestDoIt +
	                    fN2ndariesAlongStepDoIt +
	                    fN2ndariesPostStepDoIt;
      if(tN2ndariesTot>0){
	G4cout << "    :----- List of 2ndaries - "
	       << "#SpawnInStep=" << std::setw(3) << tN2ndariesTot
	       << "(Rest="  << std::setw(2) << fN2ndariesAtRestDoIt
	       << ",Along=" << std::setw(2) << fN2ndariesAlongStepDoIt
	       << ",Post="  << std::setw(2) << fN2ndariesPostStepDoIt
	       << "), "
	       << "#SpawnTotal=" << std::setw(3) << (*fSecondary).size()
	       << " ---------------"
	       << G4endl;

	for(size_t lp1=(*fSecondary).size()-tN2ndariesTot;
                        lp1<(*fSecondary).size(); lp1++){
	  G4cout << "    : "
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().x(),"Length")
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().y(),"Length")
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetPosition().z(),"Length")
		 << std::setw(6)
		 << G4BestUnit((*fSecondary)[lp1]->GetKineticEnergy(),"Energy")
		 << std::setw(10)
		 << (*fSecondary)[lp1]->GetDefinition()->GetParticleName();
	  G4cout << G4endl;
	}

	G4cout << "    :-----------------------------"
	       << "----------------------------------"
	       << "-- EndOf2ndaries Info ---------------"
	       << G4endl;
      }
    }

  }
  G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingVerbose::TrackingStarted()
{

  CopyState();
G4int prec = G4cout.precision(3);
  if( verboseLevel > 0 ){
    // G4cout << "-----------------------> " << fTrack->GetCurrentStepNumber() << G4endl;

    // boundary process
    G4OpBoundaryProcessStatus boundaryStatus=Undefined;
    static G4ThreadLocal G4OpBoundaryProcess* boundary=NULL;


    //find the boundary process only once
    if(!boundary)
    {

      G4ProcessManager* pm = fTrack->GetDefinition()->GetProcessManager();
      G4int nprocesses = pm->GetProcessListLength();
      G4ProcessVector* pv = pm->GetProcessList();
      G4int i;
      for( i=0;i<nprocesses;i++)
      {
        if((*pv)[i]->GetProcessName()=="OpBoundary")
        {
          boundary = (G4OpBoundaryProcess*)(*pv)[i];
          break;
        }
      }
    }


    G4cout << std::setw(15) << "Step#"
           << std::setw(18) << "X"
	         << std::setw(18) << "Y"
	         << std::setw(19) << "Z"
	         << std::setw(21) << "KineE"
	         << std::setw(18) << "dEStep"
	         << std::setw(18) << "StepLeng"
	         << std::setw(18) << "TrakLeng";
    G4cout << std::setw(15) << "BoundSt";
	  G4cout << std::setw(25) << "Volume"
	         << std::setw(15) << "Process"    << G4endl;

   G4cout << std::setw(15) << fTrack->GetCurrentStepNumber()
          << std::setw(15) << G4BestUnit(fTrack->GetPosition().x(),"Length")
          << std::setw(15) << G4BestUnit(fTrack->GetPosition().y(),"Length")
          << std::setw(15) << G4BestUnit(fTrack->GetPosition().z(),"Length")
          << std::setw(15) << G4BestUnit(fTrack->GetKineticEnergy(),"Energy")
          << std::setw(15) << G4BestUnit(fStep->GetTotalEnergyDeposit(),"Energy")
          << std::setw(15) << G4BestUnit(fStep->GetStepLength(),"Length")
          << std::setw(15) << G4BestUnit(fTrack->GetTrackLength(),"Length");

    if(boundary)
    {
      boundaryStatus = boundary->GetStatus();
      G4cout << std::setw(15) << OPprocessName[boundaryStatus];
    }
    else
    {
      G4cout << std::setw(15) << "Not Boundary";
    }

    if(fTrack->GetNextVolume()){
      G4cout << std::setw(25) << fTrack->GetVolume()->GetName();
    } else {
      G4cout << "OutOfWorld";
    }
    G4cout << "    initStep" << G4endl;
    // G4cout << "-----------------------> " << fTrack->GetCurrentStepNumber() << G4endl;


  }
  G4cout.precision(prec);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
