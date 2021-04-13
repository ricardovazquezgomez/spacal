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
// $Id: PrimaryGeneratorAction.hh,v 1.6 2006-06-29 17:54:04 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"
#include "TH1F.h"
#include "ConfigFile.hh"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "G4RunManager.hh"
#include "Parameters.hh"

class G4ParticleGun;
class G4Event;
class G4GeneralParticleSource;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction(ConfigFile config);
  virtual ~PrimaryGeneratorAction();

public:
  virtual void GeneratePrimaries(G4Event*);
  // virtual void MakeEnergyHisto_GaGG_Ce_Mg();
  // virtual void MakeTimeHisto_GaGG_Ce_Mg();
  // virtual void MakeEnergyHisto_YAG_Ce();
  // virtual void MakeTimeHisto_YAG_Ce();
  // virtual void MakeEnergyHisto_Ideal(float monoEnergy);
  // virtual void MakeTimeHisto_Ideal();
  // virtual void MakeAbsHisto_GaGG_Ce_Mg();
  // virtual void MakeAbsHisto_GaGG_Ce_Mg_old();
  // virtual void MakeAbsHisto_YAG_Ce();
  TH1F*      MakeAbsHisto(int num,G4Material *aMaterial);
  TH1F*        MakeEnergyHisto(int num,G4Material *aMaterial);
  TH1F*        MakeTimeHisto  (int num,G4Material *aMaterial);
  G4ParticleDefinition* GetParticleDefinition(){return particleDefinition;};
  G4double              sample_time(G4double tau1, G4double tau2);
  inline G4double       single_exp(G4double t, G4double tau2)
  {
    return std::exp(-1.0*t/tau2)/tau2;
  }
  inline G4double bi_exp(G4double t, G4double tau1, G4double tau2)
  {
    return std::exp(-1.0*t/tau2)*(1-std::exp(-1.0*t/tau1))/tau2/tau2*(tau1+tau2);
  }



private:
  // G4ParticleGun* fParticleGun;
  G4GeneralParticleSource* gun=NULL;
  G4ParticleGun * fParticleGun=NULL;
  G4ParticleDefinition* particleDefinition;
  G4ParticleTable* particleTable ;
  std::vector<TH1F*> enHisto;
  std::vector<TH1F*> tHisto;
  std::vector<TH1F*> aHisto;
  TH1F* indeal_enHisto;
  TH1F* indeal_tHisto;
  TGraph *abs_GaGG_Ce_Mg;
  TGraph *abs_GaGG_Ce_Mg_old;
  TGraph *abs_YAG_Ce;
  int primaries;
  TH1F* eHistoGen ;
  TH1F* tHistoGen ;

  int opticalMaterial;
  // float monoEnergy;
  int simulationType;
  //-----------
  G4bool use_gps;
  G4bool use_flux_fromTree;
  G4int SkipEvents;
  TFile *finput;
  // TTree *input_tree;
  TTree* source_tree;
  G4double Zshift; //shifted from LHCb coordinate to module cooridate

  G4double prod_x;
  G4double prod_y;
  G4double prod_z;
  G4double entry_x;
  G4double entry_y;
  G4double entry_z;
  G4double px;
  G4double py;
  G4double pz;
  G4double eTot;
  G4double eKinetic;
  G4int pdgID;
  G4int evtID;
  G4double timing;

  G4int partTot;
  G4double energyTot;

  G4int firstEvtID;

  //-----------
  // GAGG_Ce_Mg





};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*PrimaryGeneratorAction_h*/
