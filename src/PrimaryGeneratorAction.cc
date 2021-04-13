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
// $Id: PrimaryGeneratorAction.cc,v 1.6 2006-06-29 17:54:27 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorAction.hh"

//#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4GeneralParticleSource.hh"
#include "G4SystemOfUnits.hh"
#include "TH1F.h"
#include "TGraph.h"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "TRandom3.h"
#include "CreateTree.hh"
#include "G4Material.hh"
#include "G4MaterialPropertiesTable.hh"

#include "G4MaterialPropertyVector.hh"
#include "Parameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(ConfigFile config)
{
  //----------
  SkipEvents = 0;
  finput = NULL;
  source_tree = NULL;
  // config.readInto(use_gps,"useGPS",false);
  use_flux_fromTree = !Parameters::Instance()->main_config.use_gps;
  //----------
  // read specific runs flags
  // gRandom->SetSeed(myseed);
  // primaries = config.read("primaries",1);
  // opticalMaterial = config.read("opticalMaterial",-1);
  // simulationType      = config.read<int>("simulationType",-1); // default to -1
  // simulation types:
  // not specified        = -1
  // only energy depo     = 0
  // full ray tracing     = 1
  // hybrid               = 2
  // optical calibration  = 3
  // if user do not specify in config file, the simulation type is "free"
  // i.e. the save flags above are left untouched.
  // otherwise, those flags will be overwritten according to the simulation type
  // here, only the primaries and opticalMaterial flags need to be written again
  if(Parameters::Instance()->main_config.simulationType == -1)
  {
    // do nothing
  }
  else
  {
    if(Parameters::Instance()->main_config.simulationType == 0)
    {
      if(Parameters::Instance()->main_config.use_gps) Parameters::Instance()->main_config.primaries              = 1;
      // while if fParticleGun it's irrelevant
      Parameters::Instance()->main_config.opticalMaterial        = -1; // irrelevant anyway in this modality
    }
    if(Parameters::Instance()->main_config.simulationType == 1)
    {
      if(Parameters::Instance()->main_config.use_gps) Parameters::Instance()->main_config.primaries              = 1;
    }
    if(Parameters::Instance()->main_config.simulationType == 2)
    {
      if(Parameters::Instance()->main_config.use_gps) Parameters::Instance()->main_config.primaries              = 1;
      Parameters::Instance()->main_config.opticalMaterial        = -1; // irrelevant anyway in this modality
    }
    if(Parameters::Instance()->main_config.simulationType == 3)
    {
      if(use_flux_fromTree)
      {
        std::cout << "WARNING: if simulationType == 3 (optical calibration) you CANNOT use use particle flux from Tree. Setting useGPS = true"<< std::endl;
        Parameters::Instance()->main_config.use_gps = true;
        use_flux_fromTree = false;
      }
      // in opticalCalibrationRun (i.e. simulationType == 3) the user sets in gps the number of events per run N (with /run/beamOn N). Each event is set to 1 photon
      if(Parameters::Instance()->main_config.primaries != 1)
      {
        // cout << "WARNING: you chose an Optical Calibration run, but you set the primaries to " << Parameters::Instance()->main_config.primaries << std::endl;
        // cout << "The number of primaries in Optical Calibration mode had to be set to 1, because the actual number will be controlled by GPS file" << std::endl;
        Parameters::Instance()->main_config.primaries = 1;
      }
      else
      {
        // leave the choice of the user untouched
      }
      Parameters::Instance()->main_config.opticalMaterial        = 0; // instantaneous monocromatic
    }
  }
  //----------

  if(use_flux_fromTree)
  {

    SkipEvents  = Parameters::Instance()->main_config.skipEvents;
    TString input_filename("");
    TString input_treename("DecayTree");
    // config.readInto (input_treename, "input_tree",TString("DecayTree")) ;
    // config.readInto (input_filename, "input_file",TString("LHCbGaussSimulation_Option/input_fakeGamma_fixedE50GeV.root")) ;
    input_filename = Parameters::Instance()->main_config.input_filename;
    input_treename = Parameters::Instance()->main_config.input_treename;
    G4cout<<"Will use LHCb particle flux to generate primary particles, from file "<<input_filename.Data()<<G4endl;
    finput = new TFile(input_filename.Data());
    source_tree = (TTree *) finput->Get(input_treename.Data());
    // Zshift      = config.read<double>("deltaZ",-12500.); //shifted from LHCb coordinate to module coordinate
    Zshift = Parameters::Instance()->main_config.Zshift;
    source_tree->SetBranchAddress("prod_vertex_x", &prod_x);   //in mm
    source_tree->SetBranchAddress("prod_vertex_y", &prod_y);
    source_tree->SetBranchAddress("prod_vertex_z", &prod_z);
    source_tree->SetBranchAddress("entry_x", &entry_x);
    source_tree->SetBranchAddress("entry_y", &entry_y);
    source_tree->SetBranchAddress("entry_z", &entry_z);
    source_tree->SetBranchAddress("px", &px);       //in GeV
    source_tree->SetBranchAddress("py", &py);       //in GeV
    source_tree->SetBranchAddress("pz", &pz);       //in GeV
    source_tree->SetBranchAddress("eTot", &eTot);   //in GeV
    source_tree->SetBranchAddress("eKinetic", &eKinetic); // in GeV
    source_tree->SetBranchAddress("pdgID", &pdgID);
    source_tree->SetBranchAddress("evtIndex", &evtID);
    source_tree->SetBranchAddress("timing", &timing);
    source_tree->GetEntry(0);
    firstEvtID = max(evtID,0);
    fParticleGun = new G4ParticleGun();
  }
  else
  {
    gun = new G4GeneralParticleSource();
  }

  particleTable = G4ParticleTable::GetParticleTable();
  particleDefinition = NULL;




  //----------
  // G4GeneralParticleSource* gps = new G4GeneralParticleSource();
  // gun = gps;
  // find material optical distribution
  // particleDefinition = gun->GetParticleDefinition();
  // G4cout << particleDefinition->GetParticleName() << G4endl;







  // write the primaries number
  CreateTree::Instance()->primaries_per_event = Parameters::Instance()->main_config.primaries;
  // int events_per_run = gun->GetNumberOfParticles();
  // CreateTree::Instance()->events_per_run = events_per_run;
  // CreateTree::Instance()->abs_GaGG_Ce_Mg = NULL;
  // CreateTree::Instance()->abs_GaGG_Ce_Mg_old = NULL;
  // CreateTree::Instance()->abs_YAG_Ce = NULL;
  // // make anyway the two abs histo for GaGG_Ce_Mg and YAG_Ce
  // MakeAbsHisto_GaGG_Ce_Mg();
  // MakeAbsHisto_GaGG_Ce_Mg_old();
  // MakeAbsHisto_YAG_Ce();
  // CreateTree::Instance()->abs_GaGG_Ce_Mg = abs_GaGG_Ce_Mg;
  // CreateTree::Instance()->abs_GaGG_Ce_Mg_old = abs_GaGG_Ce_Mg_old;
  // CreateTree::Instance()->abs_YAG_Ce = abs_YAG_Ce;


  std::vector<int> CrystalMaterialList = Parameters::Instance()->GetCrystalMaterialList();
  std::cout << "Found " << CrystalMaterialList.size() << " crystal materials "<< std::endl;
  for(unsigned int i = 0 ; i < CrystalMaterialList.size(); i++)
  {
    G4Material *aMaterial = Parameters::Instance()->GetCrystalMaterial(CrystalMaterialList[i]);
    std::cout << "Crystal Material name " << aMaterial->GetName() << std::endl;
    TH1F* eHistoTemp = MakeEnergyHisto(CrystalMaterialList[i],aMaterial);
    TH1F* tHistoTemp = MakeTimeHisto(CrystalMaterialList[i],aMaterial);
    TH1F* aHistoTemp = MakeAbsHisto(CrystalMaterialList[i],aMaterial);
    enHisto.push_back(eHistoTemp);
    tHisto.push_back(tHistoTemp);
    aHisto.push_back(aHistoTemp);

    Parameters::Instance()->AddTimeHisto(tHistoTemp);
    Parameters::Instance()->AddEnergyHisto(eHistoTemp);
    Parameters::Instance()->AddAbsHisto(aHistoTemp);

  }

}

TH1F* PrimaryGeneratorAction::MakeAbsHisto(int num,G4Material *aMaterial)
{
  G4MaterialPropertiesTable* aMaterialPropertiesTable =
                               aMaterial->GetMaterialPropertiesTable();
  G4MaterialPropertyVector* absSpectrumVector = aMaterialPropertiesTable->GetProperty("ABSLENGTH");
  G4double minEnergy = absSpectrumVector->GetMinLowEdgeEnergy()/eV;
  G4double maxEnergy = absSpectrumVector->GetMaxEnergy()/eV;
  int bins = 1000;
  std::stringstream ssname;
  ssname << "absHisto_" << num;
  TH1F *absHistogram = new TH1F(ssname.str().c_str(),ssname.str().c_str(),bins,minEnergy,maxEnergy);
  // TH1F *absHistogram = new TH1F(ssname.str().c_str(),ssname.str().c_str(),bins,minEnergy,maxEnergy);
  for(int i = 1;i< bins;i++)
  {
    G4double binCenter = (absHistogram->GetBinCenter(i))*eV; // mandatory to use G4double in function Value, with correct units
    G4double value = absSpectrumVector->Value(binCenter);
    // std::cout << binCenter << " " << value << std::endl;
    absHistogram->SetBinContent(i,value);
  }
  ssname.str("");
  ssname << "Abs length of " << aMaterial->GetName();
  absHistogram->SetTitle(ssname.str().c_str());
  absHistogram->GetXaxis()->SetTitle("eV");
  absHistogram->GetYaxis()->SetTitle("mm");
  return absHistogram;

}

TH1F* PrimaryGeneratorAction::MakeEnergyHisto(int num,G4Material *aMaterial)
{
  G4MaterialPropertiesTable* aMaterialPropertiesTable =
                               aMaterial->GetMaterialPropertiesTable();
  G4MaterialPropertyVector* emissionSpectrumVector = aMaterialPropertiesTable->GetProperty("FASTCOMPONENT");
  G4double minEnergy = emissionSpectrumVector->GetMinLowEdgeEnergy()/eV;
  G4double maxEnergy = emissionSpectrumVector->GetMaxEnergy()/eV;
  // sample a 1000 bins
  int bins = 1000;
  std::stringstream ssname;
  ssname << "enHisto_" << num;
  TH1F *energyHistogram = new TH1F(ssname.str().c_str(),ssname.str().c_str(),bins,minEnergy,maxEnergy);
  for(int i = 1;i< bins;i++)
  {
    G4double binCenter = (energyHistogram->GetBinCenter(i))*eV; // mandatory to use G4double in function Value, with correct units
    G4double value = emissionSpectrumVector->Value(binCenter);
    // std::cout << binCenter << " " << value << std::endl;
    energyHistogram->SetBinContent(i,value);
  }
  ssname.str("");
  ssname << "Emission spectrum of " << aMaterial->GetName();
  energyHistogram->SetTitle(ssname.str().c_str());
  energyHistogram->GetXaxis()->SetTitle("eV");
  energyHistogram->GetYaxis()->SetTitle("Counts");
  return energyHistogram;

}

G4double PrimaryGeneratorAction::sample_time(G4double tau1, G4double tau2)
{
// tau1: rise time and tau2: decay time

        while(1) {
          // two random numbers
          G4double ran1 = G4UniformRand();
          G4double ran2 = G4UniformRand();
          //
          // exponential distribution as envelope function: very efficient
          //
          G4double d = (tau1+tau2)/tau2;
          // make sure the envelope function is
          // always larger than the bi-exponential
          G4double t = -1.0*tau2*std::log(1-ran1);
          G4double gg = d*single_exp(t,tau2);
          if (ran2 <= bi_exp(t,tau1,tau2)/gg) return t;
        }
        return -1.0;
}

TH1F* PrimaryGeneratorAction::MakeTimeHisto(int num,G4Material *aMaterial)
{
  G4MaterialPropertiesTable* aMaterialPropertiesTable =
                               aMaterial->GetMaterialPropertiesTable();
  G4MaterialPropertyVector* Fast_Intensity = aMaterialPropertiesTable->GetProperty("FASTCOMPONENT");
  G4MaterialPropertyVector* Slow_Intensity = aMaterialPropertiesTable->GetProperty("SLOWCOMPONENT");
  G4int nscnt = 1;
  if (Fast_Intensity && Slow_Intensity)
  {
    nscnt = 2;
  }

  // fast scintillation
  G4double FastScintillationDecayTime = 0.*ns;
  G4double FastScintillationRiseTime = 0.*ns;

  // slow scintillation
  G4double SlowScintillationDecayTime = 0.*ns;
  G4double SlowScintillationRiseTime = 0.*ns;

  // assign
  if(Fast_Intensity) // fast
  {
    FastScintillationRiseTime  = aMaterialPropertiesTable->GetConstProperty("FASTSCINTILLATIONRISETIME");;
    FastScintillationDecayTime = aMaterialPropertiesTable->GetConstProperty("FASTTIMECONSTANT");;
  }
  if(Slow_Intensity) // slow
  {
    SlowScintillationRiseTime  = aMaterialPropertiesTable->GetConstProperty("SLOWSCINTILLATIONRISETIME");;
    SlowScintillationDecayTime = aMaterialPropertiesTable->GetConstProperty("SLOWTIMECONSTANT");;
  }

  G4double yieldRatio = 1.0;
  if(nscnt == 2)
  {
    yieldRatio = aMaterialPropertiesTable->GetConstProperty("YIELDRATIO");
  }

  // take max decay time
  G4double maxDecayTime = 0.*ns;
  if(Fast_Intensity && !Slow_Intensity) // only fast
  {
    maxDecayTime = FastScintillationDecayTime;
  }
  else
  {
    if(!Fast_Intensity && Slow_Intensity) // only slow
    {
      maxDecayTime = SlowScintillationDecayTime;
    }
    else // both
    {
      maxDecayTime = std::max(SlowScintillationDecayTime,FastScintillationDecayTime); // should always be the slow one, but you never know...
    }
  }


  std::stringstream ssname;
  ssname << "tHisto_" << num;
  // N points

  // t histo from 0 to 10 slow decay times...
  // bins of 1 ps, t max is in ns, so tmax *1000 bins
  int bins = (int) 10000*(maxDecayTime/ns);
  std::cout << SlowScintillationDecayTime << " " << FastScintillationDecayTime << " " << maxDecayTime << " " << SlowScintillationRiseTime << " "<< FastScintillationRiseTime << " " << bins << std::endl;
  TH1F *timeHistogram = new TH1F(ssname.str().c_str(),ssname.str().c_str(),bins,0,maxDecayTime*10.0);
  ssname.str("");
  ssname << "Time distribution for scintillation, " << aMaterial->GetName();
  timeHistogram->SetTitle(ssname.str().c_str());
  timeHistogram->GetXaxis()->SetTitle("ns");
  timeHistogram->GetYaxis()->SetTitle("Counts");

  for (int scnt = 1; scnt <= nscnt; scnt++)
  {
    G4double riseTime = 0 *ns;
    G4double decayTime = 0*ns;
    // for each scintillation
    if (scnt == 1) { // first component
      if (nscnt == 1) { // and there is only one component
        if(Fast_Intensity && !Slow_Intensity) // it's fast
        {
          riseTime  = FastScintillationRiseTime;
          decayTime = FastScintillationDecayTime;
        }
        else
        {
          riseTime  = SlowScintillationRiseTime;
          decayTime = SlowScintillationDecayTime;
        }
        for(int i = 1 ; i < bins; i++)
        {
          // fill histo in position t given by geant4 strategy
          timeHistogram->SetBinContent(i,bi_exp(timeHistogram->GetBinCenter(i),riseTime,decayTime));
        }
      }
      else // first component, but there are 2 components
      {
        // take first component as fast, fill with weight
        riseTime  = FastScintillationRiseTime;
        decayTime = FastScintillationDecayTime;
        for(int i = 1 ; i < bins; i++)
        {
          // fill histo in position t given by geant4 strategy
          timeHistogram->SetBinContent(i,yieldRatio*bi_exp(timeHistogram->GetBinCenter(i),riseTime,decayTime));
        }
      }
    }
    else // second component , take as slow one
    {
      // take second component as slow, fill with weight (1-yieldRatio)
      riseTime  = SlowScintillationRiseTime;
      decayTime = SlowScintillationDecayTime;
      for(int i = 1 ; i < bins; i++)
      {
        // fill histo in position t given by geant4 strategy
        timeHistogram->SetBinContent(i, timeHistogram->GetBinContent(i) + (1.0-yieldRatio)*bi_exp(timeHistogram->GetBinCenter(i),riseTime,decayTime));
      }
    }
  }

  return timeHistogram;

}



// void PrimaryGeneratorAction::MakeAbsHisto_GaGG_Ce_Mg()
// {
//
//   const G4int nEntries_ABS = 140;
//   G4double PhotonEnergy_ABS[nEntries_ABS] =
//     { 1.23746, 1.2465, 1.25567, 1.26498, 1.27442, 1.28401, 1.29375, 1.30363, 1.31366, 1.32385, 1.3342,
//     1.34471, 1.35539, 1.36624, 1.37727, 1.38847, 1.39986, 1.41144, 1.42321, 1.43518, 1.44735, 1.45973,
//     1.47232, 1.48513, 1.49817, 1.51144, 1.52494, 1.53869, 1.55269, 1.56695, 1.58147, 1.59626, 1.61133,
//     1.62669, 1.64234, 1.6583, 1.67457, 1.69116, 1.70809, 1.72536, 1.74298, 1.76096, 1.77932, 1.79807,
//     1.81721, 1.83677, 1.84173, 1.85172, 1.86182, 1.87203, 1.88235, 1.89279, 1.90335, 1.91402, 1.92481,
//     1.93573, 1.94677, 1.95793, 1.96923, 1.98066, 1.99222, 2.00391, 2.01575, 2.02772, 2.03984, 2.0521,
//     2.06452, 2.07708, 2.0898, 2.10267, 2.1157, 2.1289, 2.14226, 2.15579, 2.16949, 2.18337, 2.19742,
//     2.21166, 2.2407, 2.25551, 2.27051, 2.28571, 2.30112, 2.31674, 2.33257, 2.34862, 2.3649, 2.3814,
//     2.39813, 2.41509, 2.4323, 2.44976, 2.46747, 2.48544, 2.50367, 2.52217, 2.54094, 2.56, 2.57935,
//     2.59898, 2.61893, 2.63918, 2.65974, 2.68063, 2.70185, 2.7234, 2.74531, 2.76757, 2.79019, 2.81319,
//     2.83657, 2.86034, 2.88451, 2.90909, 2.95954, 2.98542, 3.03858, 3.06587, 3.12195, 3.24051, 3.27157,
//     3.30323, 3.3355, 3.36842, 3.40199, 3.43624, 3.47119, 3.50685, 3.54325, 3.58042, 3.61837, 3.65714,
//     3.69675, 3.73723, 3.8209, 3.86415, 3.95367, 4.096, 5.62461};
//   G4double Absorption[nEntries_ABS] =
//     {955.954000, 668.943000, 614.032000, 467.353000, 559.162000, 579.555000, 650.626000, 535.162000, 524.755000,
//     540.255000, 600.531000, 519.526000, 526.569000, 516.447000, 502.306000, 518.467000, 522.762000, 503.951000,
//     537.711000, 520.839000, 553.001000, 488.909000, 521.548000, 471.766000, 450.077000, 404.852000, 441.221000,
//     440.960000, 411.248000, 380.195000, 399.852000, 380.692000, 372.110000, 362.841000, 370.628000, 355.130000,
//     371.095000, 389.190000, 390.557000, 376.662000, 394.540000, 419.485000, 423.900000, 424.207000, 459.364000,
//     432.117000, 450.723000, 444.802000, 451.629000, 415.643000, 458.704000, 473.048000, 388.097000, 445.714000,
//     356.658000, 362.018000, 323.960000, 249.410000, 222.770000, 173.466000, 126.428000, 94.944500, 64.388500,
//     45.198300, 30.697900, 20.705400, 14.095600, 9.563980, 6.555540, 4.590660, 3.251770, 2.439640, 1.899230,
//     1.561530, 1.644290, 1.454600, 1.344050, 1.392430, 1.474630, 1.365870, 1.418190, 1.449410, 1.429420, 1.588940,
//     1.456210, 1.313570, 1.626080, 1.611680, 1.429990, 1.472460, 1.501370, 1.565330, 1.556520, 1.413260, 1.651750,
//     1.532130, 1.908330, 2.906480, 4.388340, 6.664950, 9.720680, 12.470900, 13.226400, 11.510300, 8.747250, 6.425540,
//     4.638230, 3.356940, 2.389330, 1.730910, 1.621300, 1.349290, 1.371340, 1.454750, 1.353480, 1.092090, 1.486040,
//     1.375890, 1.153950, 1.623740, 1.443300, 1.605000, 1.705960, 1.510680, 1.525980, 1.586010, 1.616290, 1.489230,
//     1.483450, 1.411050, 1.685210, 1.409220, 1.232600, 1.707510, 1.382160, 1.322320, 1.441060, 1.478580, 0 };
//   //intrinsic absorption spectrum
//    // const G4int nEntries_ABS = 122;
//    // G4double PhotonEnergy_ABS[nEntries_ABS] =
//  //     {1.55, 1.55975, 1.56962, 1.57962, 1.58974, 1.6, 1.61039, 1.62092,
//  // 1.63158, 1.64238, 1.65333, 1.66443, 1.67568, 1.68707, 1.69863, 1.71034,
//  // 1.72222, 1.73427, 1.74648, 1.75887, 1.77143, 1.78417, 1.7971, 1.81022,
//  // 1.82353, 1.83704, 1.85075, 1.86466, 1.87879, 1.89313, 1.90769, 1.92248,
//  // 1.9375, 1.95276, 1.96825, 1.984, 2, 2.01626, 2.03279, 2.04959, 2.06667,
//  // 2.08403, 2.10169, 2.11966, 2.13793, 2.15652, 2.17544, 2.19469, 2.21429,
//  // 2.23423, 2.25455, 2.27523, 2.2963, 2.31776, 2.33962, 2.3619, 2.38462,
//  // 2.40777, 2.43137, 2.45545, 2.48, 2.50505, 2.53061, 2.5567, 2.58333,
//  // 2.61053, 2.6383, 2.66667, 2.69565, 2.72527, 2.75556, 2.78652, 2.81818,
//  // 2.85057, 2.88372, 2.91765, 2.95238, 2.98795, 3.02439, 3.06173, 3.1,
//  // 3.13924, 3.17949, 3.22078, 3.26316, 3.30667, 3.35135, 3.39726, 3.44444,
//  // 3.49296, 3.54286, 3.60465, 3.65782, 3.72372, 3.78049, 3.83901, 3.89937,
//  // 3.96166, 4.02597, 4.09241, 4.16107, 4.23208, 4.30556, 4.38163, 4.46043,
//  // 4.55882, 4.64419, 4.75096, 4.84375, 4.94024, 5.04065, 5.14523, 5.25424,
//  // 5.36797, 5.48673, 5.61086, 5.74074, 5.87678, 6.01942, 6.16915, 6.32653,
//  // 6.49215};
//  //
//  //
//  //
//  //   G4double Absorption[nEntries_ABS] =
//  //     {1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000.,
//  // 1000., 950., 900., 850., 800., 750., 700., 650., 600., 600., 550.,
//  // 525., 500., 495.113, 493.925, 490.185, 485.183, 480.656, 475.745, 474.191,
//  // 470.37, 465.699, 460.696, 455.925, 450.141, 445.712, 440.798, 435.2,
//  // 430.659, 425.727, 420.328, 415.261, 410.212, 405.257, 402.131, 401.244,
//  // 400.65, 395.193, 390.745, 385.259, 379.265, 373.136, 374.806, 326.031,
//  // 311.609, 281.734, 217.074, 156.221, 101.25, 59.9635, 32.9009, 17.6943,
//  // 9.62468, 5.33133, 3.06042, 1.86737, 0.931539, 1.17313, 1.33421, 0.5, 0.5,
//  // 1.0692, 0.5, 0.5, 0.5, 0.5, 0.994237, 1.11628, 0.5, 1.38984, 2.55123,
//  // 5.06262, 10.5534, 20.0187, 24.5102, 17.6148, 10.163, 5.66707, 3.08815,
//  // 1.63268, 0.5, 1.0365, 0.5, 0.5, 1.06797, 0.5, 0.5, 0.5, 0.926825,
//  // 0.5, 0.5, 0.926322, 0.925969, 0.5, 0.5, 1.19071, 0.5, 0.5, 0.5, 0.5,
//  // 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.25215, 1.07403,
//  // 1.16123, 1.22498};
//
//  abs_GaGG_Ce_Mg = new TGraph(nEntries_ABS,&PhotonEnergy_ABS[0],&Absorption[0]);
//  abs_GaGG_Ce_Mg->SetName("GaGG_Ce_Mg abs length");
//  abs_GaGG_Ce_Mg->SetTitle("GaGG_Ce_Mg abs length");
//  abs_GaGG_Ce_Mg->GetXaxis()->SetTitle("Energy [eV]");
//  abs_GaGG_Ce_Mg->GetYaxis()->SetTitle("Length [mm]");
//
// }

// void PrimaryGeneratorAction::MakeAbsHisto_GaGG_Ce_Mg_old()
// {
//   //intrinsic absorption spectrum
//    const G4int nEntries_ABS = 122;
//    G4double PhotonEnergy_ABS[nEntries_ABS] =
//      {1.55, 1.55975, 1.56962, 1.57962, 1.58974, 1.6, 1.61039, 1.62092,
//  1.63158, 1.64238, 1.65333, 1.66443, 1.67568, 1.68707, 1.69863, 1.71034,
//  1.72222, 1.73427, 1.74648, 1.75887, 1.77143, 1.78417, 1.7971, 1.81022,
//  1.82353, 1.83704, 1.85075, 1.86466, 1.87879, 1.89313, 1.90769, 1.92248,
//  1.9375, 1.95276, 1.96825, 1.984, 2, 2.01626, 2.03279, 2.04959, 2.06667,
//  2.08403, 2.10169, 2.11966, 2.13793, 2.15652, 2.17544, 2.19469, 2.21429,
//  2.23423, 2.25455, 2.27523, 2.2963, 2.31776, 2.33962, 2.3619, 2.38462,
//  2.40777, 2.43137, 2.45545, 2.48, 2.50505, 2.53061, 2.5567, 2.58333,
//  2.61053, 2.6383, 2.66667, 2.69565, 2.72527, 2.75556, 2.78652, 2.81818,
//  2.85057, 2.88372, 2.91765, 2.95238, 2.98795, 3.02439, 3.06173, 3.1,
//  3.13924, 3.17949, 3.22078, 3.26316, 3.30667, 3.35135, 3.39726, 3.44444,
//  3.49296, 3.54286, 3.60465, 3.65782, 3.72372, 3.78049, 3.83901, 3.89937,
//  3.96166, 4.02597, 4.09241, 4.16107, 4.23208, 4.30556, 4.38163, 4.46043,
//  4.55882, 4.64419, 4.75096, 4.84375, 4.94024, 5.04065, 5.14523, 5.25424,
//  5.36797, 5.48673, 5.61086, 5.74074, 5.87678, 6.01942, 6.16915, 6.32653,
//  6.49215};
//
//
//
//    G4double Absorption[nEntries_ABS] =
//      {1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000.,
//  1000., 950., 900., 850., 800., 750., 700., 650., 600., 600., 550.,
//  525., 500., 495.113, 493.925, 490.185, 485.183, 480.656, 475.745, 474.191,
//  470.37, 465.699, 460.696, 455.925, 450.141, 445.712, 440.798, 435.2,
//  430.659, 425.727, 420.328, 415.261, 410.212, 405.257, 402.131, 401.244,
//  400.65, 395.193, 390.745, 385.259, 379.265, 373.136, 374.806, 326.031,
//  311.609, 281.734, 217.074, 156.221, 101.25, 59.9635, 32.9009, 17.6943,
//  9.62468, 5.33133, 3.06042, 1.86737, 0.931539, 1.17313, 1.33421, 0.5, 0.5,
//  1.0692, 0.5, 0.5, 0.5, 0.5, 0.994237, 1.11628, 0.5, 1.38984, 2.55123,
//  5.06262, 10.5534, 20.0187, 24.5102, 17.6148, 10.163, 5.66707, 3.08815,
//  1.63268, 0.5, 1.0365, 0.5, 0.5, 1.06797, 0.5, 0.5, 0.5, 0.926825,
//  0.5, 0.5, 0.926322, 0.925969, 0.5, 0.5, 1.19071, 0.5, 0.5, 0.5, 0.5,
//  0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.25215, 1.07403,
//  1.16123, 1.22498};
//
//  abs_GaGG_Ce_Mg_old = new TGraph(nEntries_ABS,&PhotonEnergy_ABS[0],&Absorption[0]);
//  abs_GaGG_Ce_Mg_old->SetName("GaGG_Ce_Mg_old abs length");
//  abs_GaGG_Ce_Mg_old->SetTitle("GaGG_Ce_Mg_old abs length");
//  abs_GaGG_Ce_Mg_old->GetXaxis()->SetTitle("Energy [eV]");
//  abs_GaGG_Ce_Mg_old->GetYaxis()->SetTitle("Length [mm]");
//
// }


// void PrimaryGeneratorAction::MakeAbsHisto_YAG_Ce()
// {
//   const G4int nEntries_ABS = 84;
//   G4double PhotonEnergy_ABS[nEntries_ABS] =
//     {4.42857, 4.35088, 4.27586, 4.20339, 4.13333, 4.06557, 4, 3.93651,
//   3.875,
//     3.81538, 3.75758, 3.70149, 3.64706, 3.5942, 3.54286, 3.49296, 3.44444,
//     3.39726, 3.35135, 3.30667, 3.26316, 3.22078, 3.17949, 3.13924, 3.1,
//   3.06173,
//     3.02439, 2.98795, 2.95238, 2.91765, 2.88372, 2.85057, 2.81818, 2.78652,
//     2.75556, 2.72527, 2.69565, 2.66667, 2.6383, 2.61053, 2.58333, 2.5567,
//     2.53061, 2.50505, 2.48, 2.45545, 2.43137, 2.40777, 2.38462, 2.3619,
//     2.33962, 2.31776, 2.2963, 2.27523, 2.25455, 2.23423, 2.21429, 2.19469,
//     2.17544, 2.15652, 2.13793, 2.11966, 2.10169, 2.08403, 2.06667, 2.04959,
//     2.03279, 2.01626, 2, 1.984, 1.96825, 1.95276, 1.9375, 1.92248, 1.90769,
//     1.89313, 1.87879, 1.86466, 1.85075, 1.83704, 1.82353, 1.81022, 1.7971,
//   1.78417};
//
//   G4double Absorption[nEntries_ABS] =
//     {12.7389, 7.96854, 9.79472, 11.3642, 13.9544, 14.078, 9.48293, 4.60862,
//     0.45205, 0.14906, 0.09215, 0.05898, 0.09349, 0.09448, 4.00614, 14.055,
//     59.4712, 177.765, 322.779, 407.883, 423.384, 250.626, 107.986, 40.2593,
//     16.1867, 7.08721, 3.3827, 1.85802, 0.42111, 0.348, 0.36044, 0.37582,
//     0.39199, 0.33268, 0.27998, 0.33348, 0.21815, 0.28273, 0.37124, 0.4064,
//     1.09597, 1.45447, 2.03141, 3.62354, 6.92578, 13.9067, 29.0131, 60.2374,
//     129.476, 226.564, 553.51, 516.632, 600., 700., 800., 900.,
//     1000., 1000., 1000., 1000.2, 1000., 1000., 1000., 1000.,
//     1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000.,
//     1000., 1000., 1000., 1000., 1000., 1000., 1000., 1000.,
//     1000., 1000., 1000., 1000.};
//
//  abs_YAG_Ce = new TGraph(nEntries_ABS,&PhotonEnergy_ABS[0],&Absorption[0]);
//  abs_YAG_Ce->SetName("YAG_Ce_Mg abs length");
//  abs_YAG_Ce->SetTitle("YAG_Ce_Mg abs length");
//  abs_YAG_Ce->GetXaxis()->SetTitle("Energy [eV]");
//  abs_YAG_Ce->GetYaxis()->SetTitle("Length [mm]");
// }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  if(gun) delete gun;
  if(fParticleGun) delete fParticleGun;
}



void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  this->partTot = 0;
  this->energyTot = 0.;
  G4double EMenergyTot = 0.;

  if(fParticleGun)
  {
    G4cout<<"WARNING: can crash if particle vertex not setup properly due to mis-alignment between LHCb and proto-type coordinate system!!"<<G4endl;
    for(int ievt =0;ievt<source_tree->GetEntries(); ievt++)
    {
      source_tree->GetEntry(ievt);
      if(this->evtID-this->firstEvtID < anEvent->GetEventID()) continue;
      if(this->evtID-this->firstEvtID > anEvent->GetEventID()) break;

      particleDefinition = particleTable->FindParticle(this->pdgID);
      // Parameters::Instance()->particleDefinition = particleDefinition->GetParticleName();

      if(!particleDefinition)
      {
        G4cout<<"Particle with ID "<<pdgID<<" not found. Skip!"<<G4endl;
        continue;
      }

      if(particleDefinition->GetParticleName() == "opticalphoton" )
      {
        if(Parameters::Instance()->main_config.opticalMaterial == -1)
        {
          G4cout << "ERROR in src/PrimaryGeneratorAction.cc: Optical simulation chosen, but no optical material specified! "<< G4endl;
          G4cout << "Use the opticalMaterial key in .cfg file"<< G4endl;
          exit(-1);
        }
      }
      G4double pp = sqrt(this->px*this->px + this->py*this->py + this->pz*this->pz);
      if(pp<0.000001 || eKinetic<0.000001) continue; //unexpected zero energy

      fParticleGun->SetParticleDefinition(particleDefinition);
      fParticleGun->SetParticleEnergy(this->eKinetic*CLHEP::GeV/CLHEP::MeV);
      fParticleGun->SetParticleMomentumDirection(G4ThreeVector(px/pp, py/pp, pz/pp));
      fParticleGun->SetParticlePosition(G4ThreeVector(this->entry_x/CLHEP::mm, this->entry_y/CLHEP::mm, this->entry_z/CLHEP::mm + Zshift));
      fParticleGun->SetParticleTime(this->timing);

      // G4cout << "id = " << this->evtID << G4endl;
      fParticleGun->GeneratePrimaryVertex(anEvent);
      // info
      // std::cout << "a vertex" << std::endl;
      // std::cout << particleDefinition->GetParticleName() << std::endl;
      // std::cout << this->eKinetic*CLHEP::GeV/CLHEP::MeV << std::endl;
      // std::cout << this->entry_x/CLHEP::mm << " " << this->entry_y/CLHEP::mm << " " << this->entry_z/CLHEP::mm + Zshift << std::endl;
      // std::cout << px/pp << " " << py/pp << " " << pz/pp << std::endl;

      this->partTot +=1;
      this->energyTot +=  this->eKinetic;
      if(pdgID==22 || abs(pdgID)==11) EMenergyTot += this->eKinetic;
    }

    G4cout<<"Total primary particles to simulate, "<<this->partTot<<", with total energy "<<this->energyTot<<" GeV, and EM energy "<<EMenergyTot<<" GeV."<<G4endl;
    if(this->partTot==0)
    {
      G4cout<<G4endl<<"#######: No more events in the loop"<<G4endl;
      //G4RunManager::GetRunManager()->SetVerboseLevel(1); //not working
      //anEvent->SetEventID( 100000-1);
      //G4RunManager::GetRunManager()->TerminateEventLoop(); //not working
      G4RunManager::GetRunManager()->RunTermination();  //kind of work, continues to run all empty events
      //G4RunManager::GetRunManager()->AbortRun(); //not working
    }
    // the following line to hack the SteppingAction.cc
    particleDefinition = particleTable->FindParticle(13);
  }
  else if(gun)
  {
    // std::cout << "New Primary" << std::endl;
    particleDefinition = gun->GetParticleDefinition();
    // G4cout << particleDefinition->GetParticleName() << G4endl;

    if(Parameters::Instance()->main_config.simulationType == 3)
    {
      // if optical calibration run
      // first check if the user actually used the proper gps file...
      if(particleDefinition->GetParticleName() == "opticalphoton" )
      {

      }
      else
      {
        G4cout << "ERROR!! You chose simulationType = 3 , i.e. Optical Calibration, but you provided a GPS file that specifies primary particle different from opticalphoton!!!" << G4endl;
        G4cout << "Please use a proper GPS file... Aborting." << G4endl;
        exit(-1);
      }
    }


    // check if it's a simulation with optical primaries, and in this case
    // if an optical material has been defined
    if(particleDefinition->GetParticleName() == "opticalphoton" )
    {
      if(Parameters::Instance()->main_config.opticalMaterial > 0)
      {
        std::stringstream ssname;
        ssname << "enHisto_" << Parameters::Instance()->main_config.opticalMaterial;
        bool foundMaterial = false;
        for(int iMat = 0; iMat < enHisto.size();iMat++)
        {
          if( ssname.str() == enHisto[iMat]->GetName() )
          {
            foundMaterial = true;
            tHistoGen  = tHisto[iMat];
            eHistoGen = enHisto[iMat];
          }
        }
        if(foundMaterial == false)
        {
          G4cout << "ERROR! Wrong optical material specified! User specified opticalMaterial = " <<  Parameters::Instance()->main_config.opticalMaterial << ", but this material does not exist in the simulation!"<< " Aborting..." << std::endl;
          exit(-1);
        }
      }
    }

    // check if it's a simulation with optical primaries, and in this case
    // if an optical material has been defined
    if(particleDefinition->GetParticleName() == "opticalphoton" )
    {
      if(Parameters::Instance()->main_config.opticalMaterial == -1)
      {
        G4cout << "ERROR in src/PrimaryGeneratorAction.cc: Optical simulation chosen, but no optical material specified! "<< G4endl;
        G4cout << "Use the opticalMaterial key in .cfg file"<< G4endl;
        exit(-1);
      }
    }
    for(int p = 0 ; p < Parameters::Instance()->main_config.primaries ; p++)
    {
      if(particleDefinition->GetParticleName() == "opticalphoton" )
      {

        // time and energy
        // time is defined by scintillation time profile if user requires "real" scintillation center, otherwise 0 for insta source
        // energy is defined by scintillation emission profile, or set by macro (so no need to over write)
        if(Parameters::Instance()->main_config.opticalMaterial != 0)
        {
          gun->SetParticleTime((tHistoGen->GetRandom())*ns);
          gun->GetCurrentSource()->GetEneDist()->SetMonoEnergy((eHistoGen->GetRandom())*eV);
        }
        else // for opticalMaterial == 0, instantaneous monocromatic source, t=0 and en decided by gps
        {
          gun->SetParticleTime(0*ns);
          // gun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(monoEnergy*eV);
        }
      }
      gun->GeneratePrimaryVertex(anEvent);
    }
  }



}
