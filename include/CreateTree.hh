// Class for creating output ttrees
// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch

#ifndef CreateTree_H
#define CreateTree_H 1

#include <iostream>
#include <vector>
#include <map>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TGraph.h"
#include "G4Material.hh"
#include "Parameters.hh"


class CreateTree
{
private:

  TTree*  ftree ;
  TString fname ;

public:

  CreateTree (TString name);

  ~CreateTree () ;

  TTree*             GetTree  () const { return ftree ; } ;
  TString            GetName  () const { return fname ; } ;
  void               AddEnergyDeposit (int index, float deposit, std::map<int,float>& depositAtt) ;
  // void               AddEnergyDepositABS (int Aindex, float depositA, std::map<int,float>) ;
 // void               AddEnergyDeposit_1st_Section (int indexion, float deposit) ;
 // void               AddEnergyDeposit_2nd_Section (int index2, float deposit) ;
  void               AddScintillationPhoton (int index) ;
  void               AddCerenkovPhoton (int index) ;
  int                Fill     () ;
  bool               Write    (TFile *) ;
  bool               WriteStructure  (TFile *) ;
  void               Clear    () ;
  static CreateTree* Instance () { return fInstance ; } ;
  static CreateTree* fInstance ;
  // void               AddCrystalMaterial(int num, G4Material *aMaterial);
  // G4Material*        GetCrystalMaterial(int num);
  // std::vector<int>   GetCrystalMaterialList(){return CrystalMaterialList;};
  void               CreateTreeEnergyDepositByModule() ;
  // void               AddTimeHisto(TH1F* aHisto){tHisto.push_back(aHisto);};
  // void               AddEnergyHisto(TH1F* aHisto){enHisto.push_back(aHisto);};
  // void               AddAbsHisto(TH1F* aHisto){absHisto.push_back(aHisto);};

  // std::vector<TH1F*> enHisto;
  // std::vector<TH1F*> tHisto;
  // std::vector<TH1F*> absHisto;

  int   Run;
  int   Event ;

  std::vector<float> * energyPerCrystal ;

  std::vector<double> * inputMomentum ; // Px Py Pz E
  std::vector<double> * inputInitialPosition ; // x, y, z





  float depositedEnergyTotal ;
  float depositedEnergy_1st_Section ;
  float depositedEnergy_2nd_Section ;
  float depositedEnergyFibres ;
  float depositedEnergyFibres_1st_Section ;
  float depositedEnergyFibres_2nd_Section ;
  float depositedEnergyFibresCross;
  float depositedEnergyFibresCenter;
  float depositedEnergyFibresCorners;
  float depositedEnergy_2nd_Sect_FibresCross;
  float depositedEnergy_2nd_Sect_FibresCenter;
  float depositedEnergy_2nd_Sect_FibresCorners;
  float depositedEnergyAbsorber ;
  float depositedEnergyCell1 ;
  float depositedEnergyCell2 ;
  float depositedEnergyCell3 ;
  float depositedEnergyCell4 ;
  float depositedEnergyCell5 ;
  float depositedEnergyCell6 ;
  float depositedEnergyCell7 ;
  float depositedEnergyCell8 ;
  float depositedEnergyCell9 ;
  float depositedEnergyCell10 ;
  float depositedEnergyCell11 ;
  float depositedEnergyCell12 ;
  float depositedEnergyCell13 ;
  float depositedEnergyCell14 ;
  float depositedEnergyCell15 ;
  float depositedEnergyCell16 ;
  float depositedEnergyCell17 ;
  float depositedEnergyCell18 ;
  float depositedEnergyAbsorber_1st_Section ;
  float depositedEnergyAbsorber_2nd_Section ;
  float depositedEnergyFibres_post ;
  float depositedEnergyAbsorber_post ;
  float depositedEnergyWorld ;
  std::vector<float> * depositedEnergyFibresAtt ;

  int module_array_x;
  int module_array_y;
  int module_total;
  float * depositedEnergyByModuleFront;
  float * depositedEnergyByModuleRear;
  float * depositionTimingByModuleFront;
  float * depositionTimingByModuleRear;

  float totalTrackLengthFibres ;
  float totalTrackLengthOverThFibres ;

  int tot_phot_cer;
  int tot_phot_cer_post;
  int tot_gap_phot_cer;
  int tot_det_phot_cer;
  std::vector<float> * tot_gap_photFast_cer;
  std::vector<float> * tot_det_photFast_cer;

  // energy deposited in each fibre of a tower
  std::vector<float> * depositedEnergies ;

  std::vector<std::vector<float> > * depositedEnergiesAtt ;
  // index of the fibre where the deposit happens
  std::vector<int> * depositFibres ;

  // scintillation photons produced in each fibre of a tower
  std::vector<int> * scintillationPhotons ;
  // index of the fibre where the deposit happens
  std::vector<int> * scintillationFibres ;

  // cerenkov photons produced in each fibre of a tower
  std::vector<int> * cerenkovPhotons ;
  // index of the fibre where the deposit happens
  std::vector<int> * cerenkovFibres ;

  std::map<int, G4Material*> CrystalMaterialMap;
  std::vector<int> CrystalMaterialList;

  float Radial_stepLength;
  float Longitudinal_stepLength;
  float Radial_ion_energy_absorber[5000];
  float Longitudinal_ion_energy_absorber[5000];

  float PrimaryParticleX[1000];
  float PrimaryParticleY[1000];
  float PrimaryParticleZ[1000];
  float PrimaryParticleE[1000];

  // to be filled at the beginning of the event generation only

  TTree       * EnergyPerModuleFront;
  TTree       * EnergyPerModuleRear;

  TTree       * fibres ;
  Int_t         fibreID;
  Int_t         fibreMaterial;
  Int_t         fibreType;
  Float_t       fibreX;
  Float_t       fibreY;
  Float_t       fibreZ;
  Float_t       fibreDX;
  Float_t       fibreDY;
  Float_t       fibreDZ;
  Float_t       fibreSeparationZ;
  Int_t         fibreSections;

  TTree       * holes ;
  Int_t         holeID;
  Int_t         holeMaterial;
  Int_t         holeType;
  Float_t       holeX;
  Float_t       holeY;
  Float_t       holeZ;
  Float_t       holeDX;
  Float_t       holeDY;
  Float_t       holeDZ;
  Float_t       holeSeparationZ;
  Int_t         holeSections;

  TTree       * cells ;
  Int_t         cellID;
  Int_t         cellMaterial;
  Int_t         cellType;
  Float_t       cellX;
  Float_t       cellY;
  Float_t       cellZ;
  Float_t       cellDX;
  Float_t       cellDY;
  Float_t       cellDZ;
  Float_t       cellSeparationZ;
  Int_t         cellSections;

  TTree       * absorbers ;
  Int_t         absorberID;
  Int_t         absorberMaterial;
  Int_t         absorberType;
  Float_t       absorberX;
  Float_t       absorberY;
  Float_t       absorberZ;
  Float_t       absorberDX;
  Float_t       absorberDY;
  Float_t       absorberDZ;
  Float_t       absorberSeparationZ;
  Int_t         absorberSections;

  TTree       * modules ;
  Int_t         moduleID;
  Int_t         moduleType;
  Int_t         moduleMaterial;
  Float_t       moduleX;
  Float_t       moduleY;
  Float_t       moduleZ;
  Float_t       moduleDX;
  Float_t       moduleDY;
  Float_t       moduleDZ;
  Float_t       moduleSeparationZ;
  Int_t         moduleSections;

  TTree       * calorimeter ;
  Int_t         calorimeterID;
  Int_t         calorimeterType;
  Int_t         calorimeterMaterial;
  Float_t       calorimeterX;
  Float_t       calorimeterY;
  Float_t       calorimeterZ;
  Float_t       calorimeterDX;
  Float_t       calorimeterDY;
  Float_t       calorimeterDZ;
  Float_t       calorimeterSeparationZ;
  Int_t         calorimeterSections;

  std::vector<int>   optN;
  std::vector<float> optX;
  std::vector<float> optY;
  std::vector<float> optZ;
  std::vector<float> optE;
  // std::vector<int>   optType;

  // float run_optX;
  // float run_optY;
  // float run_optZ;
  // float run_optE;

  bool SaveOpticalRunVector;

  TTree           *shower ;
  TTree           *LAPPD ;

  Int_t           showerIsInCrystal;
  // Int_t           crystalID;
  // Int_t           shower_cellID;
  // Int_t           shower_moduleID;
  Int_t           shower_moduleType;
  // Float_t         shower_module_x;
  // Float_t         shower_module_y;
  // Float_t         shower_module_z;
  Int_t           pdgID;
  Int_t           trackID;
  Int_t           primaryID;
  Int_t           primaryPDGID;
  Float_t         primaryEnergy;
  Float_t         showerX    ;
  Float_t         showerY    ;
  Float_t         showerZ    ;
  Float_t         showerT    ;
  Float_t         showerPx   ;
  Float_t         showerPy   ;
  Float_t         showerPz   ;
  Float_t         showerTotalEnDep ;
  Float_t         showerIonizingEnDep;
  Float_t         showerNonIonizingEnDep;
  std::string     showerProcessName;
  std::string     showerMaterialName;
  std::string     showerMaterialNamePre;
  std::string     showerMaterialNamePost;
  Int_t           showerMaterialNumber;
  Float_t         showerLocalX    ;
  Float_t         showerLocalY    ;
  Float_t         showerLocalZ    ;


  TTree       * primaries ;


  Float_t         primaryPositionAtVertexX    ;
  Float_t         primaryPositionAtVertexY    ;
  Float_t         primaryPositionAtVertexZ    ;
  Float_t         primaryMomentumAtVertexX    ;
  Float_t         primaryMomentumAtVertexY    ;
  Float_t         primaryMomentumAtVertexZ    ;
  Float_t         primaryEnergyAtVertex       ;

  Float_t         primaryPositionOnAbsorberX    ;
  Float_t         primaryPositionOnAbsorberY    ;
  Float_t         primaryPositionOnAbsorberZ    ;
  Float_t         primaryMomentumOnAbsorberX    ;
  Float_t         primaryMomentumOnAbsorberY    ;
  Float_t         primaryMomentumOnAbsorberZ    ;
  Float_t         primaryPositionOnAbsorberT   ;
  Float_t         primaryEnergyOnAbsorber       ;

  TTree       *simple_time;
  std::vector<double> listOfTimestamps_front_back_0;
  std::vector<double> listOfTimestamps_front_back_1;
  double timestamp_front_back_0;
  double timestamp_front_back_1;

  // TTree       * optCaliInfo ;

  TTree * summary;
  Int_t primary_type   ;
  Float_t total_energy_deposited;

  bool primaryFirstEntranceFound;


  TTree * attenuationLengths ;

  TTree * photons;
  // int   eventNumber      ;
  int   front_back      ;
  int   pmt_number      ;
  int   module_number       ;
  double pde;

  // int   phPerPMT[2][9];
  int primaries_per_event;
  int events_per_run;
  float vertX         ;
  float vertY         ;
  float vertZ         ;
  float vertMomentumX  ;
  float vertMomentumY  ;
  float vertMomentumZ  ;
  float PositionX     ;
  float PositionY     ;
  float PositionZ     ;
  float PreMomentumX  ;
  float PreMomentumY  ;
  float PreMomentumZ  ;
  float PostMomentumX ;
  float PostMomentumY ;
  float PostMomentumZ ;
  float globalTime    ;
  float localTime    ;
  float PhotonEnergy  ;
  std::string processName;
  int photonDetected;
  Int_t           photon_cellID;
  Int_t           photon_moduleID;
  Int_t           photon_moduleType;
  Float_t         photon_module_x;
  Float_t         photon_module_y;
  Float_t         photon_module_z;

  // Float_t         photons_primaryPositionOnAbsorberX    ;
  // Float_t         photons_primaryPositionOnAbsorberY    ;
  // Float_t         photons_primaryPositionOnAbsorberZ    ;
  // Float_t         photons_primaryMomentumOnAbsorberX    ;
  // Float_t         photons_primaryMomentumOnAbsorberY    ;
  // Float_t         photons_primaryMomentumOnAbsorberZ    ;
  // Float_t         photons_primaryEnergyOnAbsorber       ;

  float zSeparationAbsolutePosition;
  TTree* opticalPhotonGen ;
  float  gen_x            ;
  float  gen_y            ;
  float  gen_z            ;
  float  gen_ux           ;
  float  gen_uy           ;
  float  gen_uz           ;
  float  gen_energy       ;
  float  gen_globalTime   ;
  float  gen_localTime    ;
  std::string gen_process ;
  // Float_t         gen_primaryPositionOnAbsorberX    ;
  // Float_t         gen_primaryPositionOnAbsorberY    ;
  // Float_t         gen_primaryPositionOnAbsorberZ    ;
  // Float_t         gen_primaryMomentumOnAbsorberX    ;
  // Float_t         gen_primaryMomentumOnAbsorberY    ;
  // Float_t         gen_primaryMomentumOnAbsorberZ    ;
  // Float_t         gen_primaryEnergyOnAbsorber       ;

  TTree * photonsAbsPoint;
  // int abs_eventNumber;
  float abs_x         ;
  float abs_y          ;
  float abs_z           ;
  float abs_globalTime   ;
  float abs_PhotonEnergy  ;

  std::vector<float> * attLengths;

  // TH1F *enHisto;
  // TH1F  *tHisto;
  TGraph* abs_GaGG_Ce_Mg;
  TGraph* abs_GaGG_Ce_Mg_old;
  TGraph* abs_YAG_Ce;

  int pipe;
  int pipe_modules_nx;
  int pipe_modules_ny;
  float absorberLength;
};

#endif
