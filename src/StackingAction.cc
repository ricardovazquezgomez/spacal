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

#include "StackingAction.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
// #include "StackingActionMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
// #include "G4String.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction()
 : G4UserStackingAction(),
   fStage(0)
{
  // fMessenger = new StackingActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  // start with all new tracks in waiting
  G4ClassificationOfNewTrack classification = fWaiting;
  // get particle definition and step number
  G4ParticleDefinition * particleType = aTrack->GetDefinition();
  // G4int nStep   = aTrack -> GetCurrentStepNumber();
  if(aTrack->GetParentID()!=0) // non primary, because in optical photons as primary need to survive (for optical calibrations!)
  {
    if(particleType == G4OpticalPhoton::OpticalPhotonDefinition())
    {
      G4String processName = aTrack->GetCreatorProcess()->GetProcessName();
      if(processName == "Scintillation")
      {
        if(!Parameters::Instance()->main_config.propagateScintillation)
        {
          return fKill;
        }
      }
      if(processName == "Cerenkov")
      {
        if(!Parameters::Instance()->main_config.propagateCerenkov)
        {
          return fKill;
        }
      }
    }
  }
  
  return fUrgent;
  
  
  // // start with all new tracks in waiting
  // G4ClassificationOfNewTrack classification = fWaiting;
  // // get particle definition and step number
  // G4ParticleDefinition * particleType = aTrack->GetDefinition();
  // // G4int nStep   = aTrack -> GetCurrentStepNumber();

  // // decide what stack, based on stage
  // // sim divided in 2+ stages
  // // 1. propagate only the prymary, puh the rest to fWaiting
  // // 2. propagate all the fWaiting, but not the optical photons
  // // 3. propagate the optical photons, only if user wants to. And save the gen position, if user wants to
  // switch(fStage)
  // {
  // case 0: // Stage 0 : Primary only
  //   // G4cout << "Case 0, aTrack->GetParentID() = " <<  aTrack->GetParentID() << G4endl;
  //   if(aTrack->GetParentID()==0)
  //   {
  //     classification = fUrgent;
  //     // special case:
  //     // if primary is actually an optical photon
  //     if(particleType == G4OpticalPhoton::OpticalPhotonDefinition() )
  //     {
  //       G4String processName = "Manual";
  //       // if this is:
  //       // 2. User wants to save optical photons at generation
  //       // 3. OptPhoton was created in a crystal (to avoid logging air cherenkov etc)
  //       if ( Parameters::Instance()->main_config.savePhotonGen  )
  //       {
  //         CreateTree::Instance()->gen_x            = aTrack->GetPosition().x();
  //         CreateTree::Instance()->gen_y            = aTrack->GetPosition().y();
  //         CreateTree::Instance()->gen_z            = aTrack->GetPosition().z();
  //         CreateTree::Instance()->gen_ux           = aTrack->GetMomentumDirection().x();
  //         CreateTree::Instance()->gen_uy           = aTrack->GetMomentumDirection().y();
  //         CreateTree::Instance()->gen_uz           = aTrack->GetMomentumDirection().z();
  //         CreateTree::Instance()->gen_energy       = aTrack->GetDynamicParticle()->GetTotalEnergy()/CLHEP::eV;
  //         CreateTree::Instance()->gen_globalTime   = aTrack->GetGlobalTime()/CLHEP::ns;
  //         CreateTree::Instance()->gen_localTime    = aTrack->GetLocalTime()/CLHEP::ns;
  //         CreateTree::Instance()->gen_process      = processName;
  //         CreateTree::Instance()->opticalPhotonGen->Fill();
  //       }
  //       // break; // so remain fWaiting
  //     }
  //   }
  //   break;

  // case 1: // Stage 1 : secondaries not opticalphoton
  //   // G4cout << "Case 1, aTrack->GetParentID() = " <<  aTrack->GetParentID() << G4endl;
  //   if(particleType == G4OpticalPhoton::OpticalPhotonDefinition() )
  //   {
  //     break; // so remain fWaiting
  //   }
  //   // otherwise...
  //   classification = fUrgent;
  //   break;

  // default: // Stage 2 : propagate all the rest (i.e. optical photons)
  //          // or kill if requested by user
  //          // or kill and save if requested
  //   // if not Optical Photon, something was wrong
  //   // G4cout << "default, aTrack->GetParentID() = " <<  aTrack->GetParentID() << G4endl;
  //   if(particleType != G4OpticalPhoton::OpticalPhotonDefinition() )
  //   {
  //     G4cout << "ERROR! A non-optical photon eneded up in stage > 1 in StackingAction??? Please kill sim and check... " << G4endl;// something was wrong!
  //     break;
  //   }
  //   else // it's optical photon
  //   {
  //     // G4cout << "in opt phot " << nStep << G4endl;
  //     // set all to fUrgent, but then check
  //     classification = fUrgent;
  //     // get creator process name
  //     G4String processName = aTrack->GetCreatorProcess()->GetProcessName();
  //     // if this is:
  //     // 2. User wants to save optical photons at generation
  //     // 3. OptPhoton was created in a crystal (to avoid logging air cherenkov etc)
  //     if ( (Parameters::Instance()->main_config.savePhotonGen) && (aTrack->GetVolume()->GetLogicalVolume()->GetName().contains("crystal") ) )
  //     {
  //       CreateTree::Instance()->gen_x            = aTrack->GetPosition().x();
  //       CreateTree::Instance()->gen_y            = aTrack->GetPosition().y();
  //       CreateTree::Instance()->gen_z            = aTrack->GetPosition().z();
  //       CreateTree::Instance()->gen_ux           = aTrack->GetMomentumDirection().x();
  //       CreateTree::Instance()->gen_uy           = aTrack->GetMomentumDirection().y();
  //       CreateTree::Instance()->gen_uz           = aTrack->GetMomentumDirection().z();
  //       CreateTree::Instance()->gen_energy       = aTrack->GetDynamicParticle()->GetTotalEnergy()/CLHEP::eV;
  //       CreateTree::Instance()->gen_globalTime   = aTrack->GetGlobalTime()/CLHEP::ns;
  //       CreateTree::Instance()->gen_localTime    = aTrack->GetLocalTime()/CLHEP::ns;
  //       CreateTree::Instance()->gen_process      = processName;
  //       CreateTree::Instance()->opticalPhotonGen->Fill();
  //     }
  //     // then stop propagating the optical photons if propagate flags are off
  //     if(processName == "Scintillation")
  //     {
  //       if(!Parameters::Instance()->main_config.propagateScintillation)
  //       {
  //         classification = fKill;
  //       }
  //       break;
  //     }
  //     if(processName == "Cerenkov")
  //     {
  //       if(!Parameters::Instance()->main_config.propagateCerenkov)
  //       {
  //         classification = fKill;
  //       }
  //       break;
  //     }
  //   } // end of optical photon, stage > 1
  //   // classification = fUrgent;

  // } // end of switch

  // return classification;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// G4bool StackingAction::InsideRoI(const G4Track * aTrack,G4double ang)
// {
//   if(!fMuonHits)
//   { fMuonHits = (RE05MuonHitsCollection*)GetCollection("muonCollection"); }
//   if(!fMuonHits)
//   { G4cerr << "muonCollection NOT FOUND" << G4endl;
//     return true; }
//
//   G4int nhits = fMuonHits->entries();
//
//   const G4ThreeVector trPos = aTrack->GetPosition();
//   for(G4int i=0;i<nhits;i++)
//   {
//     G4ThreeVector muHitPos = (*fMuonHits)[i]->GetPos();
//     G4double angl = muHitPos.angle(trPos);
//     if(angl<ang) { return true; }
//   }
//
//   return false;
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// G4VHitsCollection* StackingAction::GetCollection(G4String colName)
// {
//   G4SDManager* SDMan = G4SDManager::GetSDMpointer();
//   G4RunManager* runMan = G4RunManager::GetRunManager();
//   int colID = SDMan->GetCollectionID(colName);
//   if(colID>=0)
//   {
//     const G4Event* currentEvent = runMan->GetCurrentEvent();
//     G4HCofThisEvent* HCE = currentEvent->GetHCofThisEvent();
//     return HCE->GetHC(colID);
//   }
//   return 0;
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StackingAction::NewStage()
{
  // do nothing 
  fStage++;
  stackManager->ReClassify();
//   fStage++;
//   G4int nhits;
//   if(fStage==1)
//   {
// //   // Stage 0->1 : Inform the user
//     // G4cout << "Finished propagating Primary particle(s). Moving on." << G4endl;
//     stackManager->ReClassify();
//     return;
//   }

//   else if(fStage==2)
//   {
// //   // Stage 1->2 : check the isolation of muon tracks
//     // G4cout << "Finished propagating all particle that are not Optical Photons. Moving to Optical Photons." << G4endl;
//     stackManager->ReClassify();
//     return;
//   }

//   else
//   {
//   // Other fStage change : just re-classify
//     stackManager->ReClassify();
//   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StackingAction::PrepareNewEvent()
{
  fStage = 0;
  // fTrkHits = 0;
  // fMuonHits = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
