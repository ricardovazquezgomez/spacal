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
// $Id: RunAction.cc,v 1.10 2006-06-29 17:54:31 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Make this appear first!

#include "RunAction.hh"
#include "Parameters.hh"

#include "G4Timer.hh"
#include "G4Run.hh"
#include "CreateTree.hh"
#include "TRandom3.h"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
  timer = new G4Timer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete timer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  if(Parameters::Instance()->main_config.verbosity > 0) G4cout << "### Run :: " << aRun->GetRunID() << " started ..." << G4endl;
  // set event number in photons tree
  CreateTree::Instance ()->Run = aRun->GetRunID();
  CreateTree::Instance ()->primaryFirstEntranceFound = false;

  timer->Start();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{

  // CreateTree::Instance()->events_per_run = aRun->GetNumberOfEvent();
  if(Parameters::Instance()->main_config.simulationType == 3)
  {
    CreateTree::Instance()->optN.push_back(aRun->GetNumberOfEvent());
    // CreateTree::Instance()->optY.push_back(y);
    // CreateTree::Instance()->optZ.push_back(z);
    // CreateTree::Instance()->optE.push_back(InitEnergy/CLHEP::eV);
  }
  // CreateTree::Instance()->optCaliInfo->Fill();
  //
  timer->Stop();
  if(Parameters::Instance()->main_config.verbosity > 0) G4cout << "number of event = " << aRun->GetNumberOfEvent()
                                                               << " " << *timer << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
