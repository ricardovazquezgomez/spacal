// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch

// compile with
// g++ -o ../build/readHybrid readHybrid.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore

#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3I.h"
#include "TString.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TTree.h"
#include "TF2.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TTreeFormula.h"
#include "TMath.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TError.h"
#include "TEllipse.h"
#include "TFormula.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TMultiGraph.h"
#include "TCutG.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TVector.h"
#include "TRandom3.h"
// #include "TRatioPlot.h"
#include "TLegendEntry.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <getopt.h>
#include <algorithm>    // std::sort

#include <sys/stat.h>
#include <dirent.h>
#include "regions.h"

void usage();

struct detector_t
{
  std::vector<Float_t> listOfTimestamps;
  std::vector<Float_t> listOft_dep     ;
  std::vector<Float_t> listOft_gen     ;
  std::vector<Float_t> listOft_prop    ;
  std::vector<Float_t> listOft_dep_gen ;
  std::vector<Float_t> listOft_dep_prop;
  std::vector<Float_t> listOft_gen_prop;
  std::vector<Float_t> u_listOfTimestamps;
  std::vector<Float_t> u_listOft_dep     ;
  std::vector<Float_t> u_listOft_gen     ;
  std::vector<Float_t> u_listOft_prop    ;
  std::vector<Float_t> u_listOft_dep_gen ;
  std::vector<Float_t> u_listOft_dep_prop;
  std::vector<Float_t> u_listOft_gen_prop;
  Float_t timestamp;
  Float_t t_dep     ;
  Float_t t_gen     ;
  Float_t t_prop    ;
  Float_t t_dep_gen ;
  Float_t t_dep_prop;
  Float_t t_gen_prop;
  Float_t u_timestamp;
  Float_t u_t_dep     ;
  Float_t u_t_gen     ;
  Float_t u_t_prop    ;
  Float_t u_t_dep_gen ;
  Float_t u_t_dep_prop;
  Float_t u_t_gen_prop;
  long int counts;
};

struct module_t
{
  int id;
  detector_t** detector;
};

struct event_t
{
  module_t* module;

  bool primarySet;
  float primaryPositionOnAbsorberX    ;
  float primaryPositionOnAbsorberY    ;
  float primaryPositionOnAbsorberZ    ;
  float primaryMomentumOnAbsorberX    ;
  float primaryMomentumOnAbsorberY    ;
  float primaryMomentumOnAbsorberZ    ;
  float primaryEnergyOnAbsorber       ;

};

int main(int argc, char** argv)
{

  // check args
  if(argc < 2)
  {
    std::cout << argv[0] << std::endl;
    usage();
    return 1;
  }

  // save command line string
  std::stringstream streamCommandLine;
  for (int i = 0; i < argc ; i++)
  {
    streamCommandLine << argv[i] << " ";
  }
  std::string CommandLineString(streamCommandLine.str());

  // save info on pc
  char hostname[1024];
  hostname[1023] = '\0';
  gethostname(hostname, 1023);
  std::string HostNameString(hostname);
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::string PWDstring(cwd);

  // input and output file names
  std::string inputFileName = "";
  std::string outputFileName = "output.root";
  UInt_t seed = 0;
  float qe = 0.07;
  int Nevents = -1;
  int Ndet = -1;
  int Nmod = -1;
  float sigmaSPTR = 0.087; // in ns
  int numb_of_phot_for_time_average = 5;
  int verbose = 0;

  // parse command line arguments
  static struct option longOptions[] =
  {
    { "input", required_argument, 0, 0 },
    { "output", required_argument, 0, 0 },
    { "seed", required_argument, 0, 0 },
    { "qe", required_argument, 0, 0 },
    { "events", required_argument, 0, 0 },
    // { "detectors", required_argument, 0, 0 },
    { "photons", required_argument, 0, 0 },
    { "sptr", required_argument, 0, 0 },
    { "verbose", required_argument, 0, 0 },
    // { "modules", required_argument, 0, 0 },
    { NULL, 0, 0, 0 }
  };
  while(1) {
    int optionIndex = 0;
    int c = getopt_long(argc, argv, "i:o:s:", longOptions, &optionIndex);
    if (c == -1) {
      break;
    }
    else if (c == 'i'){
      inputFileName = (char *)optarg;
    }
    else if (c == 'o'){
      outputFileName = (char *)optarg;
    }
    else if (c == 's'){
      seed = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 0){
      inputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 1){
      outputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 2){
      seed = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 3){
      qe = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 4){
      Nevents = atoi((char *)optarg);
    }
    // else if (c == 0 && optionIndex == 5){
      // Ndet = atoi((char *)optarg);
    // }
    else if (c == 0 && optionIndex == 5){
      numb_of_phot_for_time_average = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 6){
      sigmaSPTR = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 7){
      verbose = atoi((char *)optarg);
    }
    // else if (c == 0 && optionIndex == 9){
      // Nmod = atoi((char *)optarg);
    // }
    else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
      usage();
      return 1;
    }
  }

  // check if input file name is given
  if(inputFileName == "")
  {
    std::cout << std::endl;
    std::cout << "ERROR! You need to an input file name!" << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }

  // random seed
  TRandom3 *rand = new TRandom3(seed);
  UInt_t randomSeed = rand->GetSeed();

  // feedback
  std::stringstream feedbackString;
  feedbackString << std::endl;
  feedbackString << CommandLineString << std::endl;
  feedbackString << std::endl;
  feedbackString << "|***************************************************************|" << std::endl;
  feedbackString << "| HYBRID MONTE CARLO FOR SPACAL-RD                              |" << std::endl;
  feedbackString << "| READING OUTPUT FILE                                           |" << std::endl;
  feedbackString << "|***************************************************************|" << std::endl;
  feedbackString << "" << std::endl;
  feedbackString << "COMMAND LINE PARAMETERS:" << std::endl;
  feedbackString << "Input file                  = " << inputFileName <<std::endl;
  feedbackString << "Output file                 = " << outputFileName <<std::endl;
  feedbackString << "QE                          = " << qe << std::endl;
  feedbackString << "Number of events            = " << Nevents << "\t (-1 means not given)" <<std::endl;
  feedbackString << "Number of detectors         = " << Ndet    << "\t (-1 means not given)" <<std::endl;
  feedbackString << "N of photons for time avg   = " << numb_of_phot_for_time_average << std::endl;
  feedbackString << "SPTR sigma [ns]             = " << sigmaSPTR << std::endl;
  feedbackString << "Verbose level               = " << verbose << std::endl;

  feedbackString << std::endl;


  //----------------------------------//
  // INPUT                            //
  //----------------------------------//
  // open input file
  TFile *fileIn = TFile::Open(inputFileName.c_str());
  // get TTrees
  // hybrid
  // declare variables
  TTree* hybrid = (TTree*)  fileIn->Get("hybrid");

  Int_t       hb_event;
  Float_t     hb_x;
  Float_t     hb_y;
  Float_t     hb_z;
  Int_t       hb_detector;
  Int_t       hb_module;
  Int_t       hb_moduleType;
  Float_t     hb_timestamp;
  Float_t     hb_t_deposition;
  Float_t     hb_t_generation;
  Float_t     hb_t_propagation;
  Float_t     hb_photon_energy;
  std::string *hb_processName = 0;
  Int_t       hb_processNumber;
  float      hb_primaryPositionOnAbsorberX    ;
  float      hb_primaryPositionOnAbsorberY    ;
  float      hb_primaryPositionOnAbsorberZ    ;
  float      hb_primaryMomentumOnAbsorberX    ;
  float      hb_primaryMomentumOnAbsorberY    ;
  float      hb_primaryMomentumOnAbsorberZ    ;
  float      hb_primaryEnergyOnAbsorber       ;
  TBranch *b_hb_event;
  TBranch *b_hb_x;
  TBranch *b_hb_y;
  TBranch *b_hb_z;
  TBranch *b_hb_detector;
  TBranch *b_hb_module;
  TBranch *b_hb_moduleType;
  TBranch *b_hb_timestamp;
  TBranch *b_hb_t_deposition;
  TBranch *b_hb_t_generation;
  TBranch *b_hb_t_propagation;
  TBranch *b_hb_photon_energy;
  TBranch *b_hb_processNumber;
  TBranch *b_hb_primaryPositionOnAbsorberX    ;
  TBranch *b_hb_primaryPositionOnAbsorberY    ;
  TBranch *b_hb_primaryPositionOnAbsorberZ    ;
  TBranch *b_hb_primaryMomentumOnAbsorberX    ;
  TBranch *b_hb_primaryMomentumOnAbsorberY    ;
  TBranch *b_hb_primaryMomentumOnAbsorberZ    ;
  TBranch *b_hb_primaryEnergyOnAbsorber       ;
  // set adreesses
  hybrid->SetBranchAddress("event"                       , &hb_event                     ,&b_hb_event);
  hybrid->SetBranchAddress("x"                           , &hb_x                         ,&b_hb_x);
  hybrid->SetBranchAddress("y"                           , &hb_y                         ,&b_hb_y);
  hybrid->SetBranchAddress("z"                           , &hb_z                         ,&b_hb_z);
  hybrid->SetBranchAddress("detector"                    , &hb_detector                  ,&b_hb_detector);
  hybrid->SetBranchAddress("module"                      , &hb_module                    ,&b_hb_module);
  hybrid->SetBranchAddress("moduleType"                  , &hb_moduleType                ,&b_hb_moduleType);
  hybrid->SetBranchAddress("timestamp"                   , &hb_timestamp                 ,&b_hb_timestamp);
  hybrid->SetBranchAddress("t_deposition"                , &hb_t_deposition              ,&b_hb_t_deposition);
  hybrid->SetBranchAddress("t_generation"                , &hb_t_generation              ,&b_hb_t_generation);
  hybrid->SetBranchAddress("t_propagation"               , &hb_t_propagation             ,&b_hb_t_propagation);
  hybrid->SetBranchAddress("PhotonEnergy"                , &hb_photon_energy             ,&b_hb_photon_energy);
  hybrid->SetBranchAddress("processName"                 , &hb_processName);
  hybrid->SetBranchAddress("processNumber"               , &hb_processNumber             ,&b_hb_processNumber);
  hybrid->SetBranchAddress("primary_PositionOnAbsorberX" , &hb_primaryPositionOnAbsorberX,&b_hb_primaryPositionOnAbsorberX );
  hybrid->SetBranchAddress("primary_PositionOnAbsorberY" , &hb_primaryPositionOnAbsorberY,&b_hb_primaryPositionOnAbsorberY );
  hybrid->SetBranchAddress("primary_PositionOnAbsorberZ" , &hb_primaryPositionOnAbsorberZ,&b_hb_primaryPositionOnAbsorberZ );
  hybrid->SetBranchAddress("primary_MomentumOnAbsorberX" , &hb_primaryMomentumOnAbsorberX,&b_hb_primaryMomentumOnAbsorberX );
  hybrid->SetBranchAddress("primary_MomentumOnAbsorberY" , &hb_primaryMomentumOnAbsorberY,&b_hb_primaryMomentumOnAbsorberY );
  hybrid->SetBranchAddress("primary_MomentumOnAbsorberZ" , &hb_primaryMomentumOnAbsorberZ,&b_hb_primaryMomentumOnAbsorberZ );
  hybrid->SetBranchAddress("primary_EnergyOnAbsorber"    , &hb_primaryEnergyOnAbsorber   ,&b_hb_primaryEnergyOnAbsorber    );

  // shower
  TTree* shower = (TTree*) fileIn->Get("shower");
  int    showerRun;
  int    showerEvent           ;
  int    showerIsInCrystal     ;
  int    showerCrystalID       ;
  float showerX               ;
  float showerY               ;
  float showerZ               ;
  float showerT               ;
  float showerTotalEnDep      ;
  float showerIonizingEnDep   ;
  float showerNonIonizingEnDep;
  std::string      *showerProcessName = 0;
  std::string      *showerMaterialName = 0;
  int    showerMaterialNumber;
  float shower_primaryPositionOnAbsorberX    ;
  float shower_primaryPositionOnAbsorberY    ;
  float shower_primaryPositionOnAbsorberZ    ;
  float shower_primaryMomentumOnAbsorberX    ;
  float shower_primaryMomentumOnAbsorberY    ;
  float shower_primaryMomentumOnAbsorberZ    ;
  float shower_primaryEnergyOnAbsorber       ;
  TBranch *b_showerRun;
  TBranch *b_showerEvent;
  TBranch *b_showerIsInCrystal;
  TBranch *b_showerCrystalID;
  TBranch *b_showerX    ;
  TBranch *b_showerY    ;
  TBranch *b_showerZ    ;
  TBranch *b_showerT    ;
  TBranch *b_showerTotalEnDep ;
  TBranch *b_showerIonizingEnDep;
  TBranch *b_showerNonIonizingEnDep;
  TBranch *b_showerMaterialNumber;
  TBranch *b_shower_primaryPositionOnAbsorberX    ;
  TBranch *b_shower_primaryPositionOnAbsorberY    ;
  TBranch *b_shower_primaryPositionOnAbsorberZ    ;
  TBranch *b_shower_primaryMomentumOnAbsorberX    ;
  TBranch *b_shower_primaryMomentumOnAbsorberY    ;
  TBranch *b_shower_primaryMomentumOnAbsorberZ    ;
  TBranch *b_shower_primaryEnergyOnAbsorber       ;
  shower->SetBranchAddress("run"                ,&showerRun             ,&b_showerRun             );
  shower->SetBranchAddress("event"              ,&showerEvent           ,&b_showerEvent           );
  shower->SetBranchAddress("isInCrystal"        ,&showerIsInCrystal     ,&b_showerIsInCrystal     );
  shower->SetBranchAddress("crystalID"          ,&showerCrystalID       ,&b_showerCrystalID       );
  shower->SetBranchAddress("x"                  ,&showerX               ,&b_showerX               );
  shower->SetBranchAddress("y"                  ,&showerY               ,&b_showerY               );
  shower->SetBranchAddress("z"                  ,&showerZ               ,&b_showerZ               );
  shower->SetBranchAddress("t"                  ,&showerT               ,&b_showerT               );
  shower->SetBranchAddress("totalEnDep"         ,&showerTotalEnDep      ,&b_showerTotalEnDep      );
  shower->SetBranchAddress("ionizingEnDep"      ,&showerIonizingEnDep   ,&b_showerIonizingEnDep   );
  shower->SetBranchAddress("nonIonizingEnDep"   ,&showerNonIonizingEnDep,&b_showerNonIonizingEnDep);
  shower->SetBranchAddress("processName"        ,&showerProcessName);
  shower->SetBranchAddress("materialName"       ,&showerMaterialName);
  shower->SetBranchAddress("materialNumber"     ,&showerMaterialNumber  ,&b_showerMaterialNumber             );
  shower->SetBranchAddress("primary_PositionOnAbsorberX" ,&shower_primaryPositionOnAbsorberX    ,&b_shower_primaryPositionOnAbsorberX);
  shower->SetBranchAddress("primary_PositionOnAbsorberY" ,&shower_primaryPositionOnAbsorberY    ,&b_shower_primaryPositionOnAbsorberY);
  shower->SetBranchAddress("primary_PositionOnAbsorberZ" ,&shower_primaryPositionOnAbsorberZ    ,&b_shower_primaryPositionOnAbsorberZ);
  shower->SetBranchAddress("primary_MomentumOnAbsorberX" ,&shower_primaryMomentumOnAbsorberX    ,&b_shower_primaryMomentumOnAbsorberX);
  shower->SetBranchAddress("primary_MomentumOnAbsorberY" ,&shower_primaryMomentumOnAbsorberY    ,&b_shower_primaryMomentumOnAbsorberY);
  shower->SetBranchAddress("primary_MomentumOnAbsorberZ" ,&shower_primaryMomentumOnAbsorberZ    ,&b_shower_primaryMomentumOnAbsorberZ);
  shower->SetBranchAddress("primary_EnergyOnAbsorber"    ,&shower_primaryEnergyOnAbsorber       ,&b_shower_primaryEnergyOnAbsorber   );
  //----------------------------------//





  //----------------------------------//
  // COUNT EVENTS                     //
  //----------------------------------//
  // quick count of events...
  long int Nhybrid = hybrid->GetEntries();
  long int counter = 0;
  std::cout << "Hybrid entries = " << Nhybrid << std::endl;




  if(Nevents == -1)
  {
    std::cout << "Number of events not provided, need to check in the input file ..." << std::endl;
    std::cout << "(if you want to provide it, run with --events)" << std::endl;
    std::cout << "Counting events ..." << std::endl;

    std::vector<int> eventsNumbers;
    // std::vector<int> detectorNumbers;
    for(int i = 0 ; i < Nhybrid; i++)
    {
      hybrid->GetEntry(i);
      bool eventIsThere= false;
      for(int iEv = 0; iEv < eventsNumbers.size(); iEv++)
      {
        if(eventsNumbers[iEv] == hb_event)
        {
          eventIsThere = true;
        }
      }
      if(!eventIsThere)
      {
        eventsNumbers.push_back(hb_event);
      }

      bool detectorIsThere= false;
      // for(int iEv = 0; iEv < detectorNumbers.size(); iEv++)
      // {
      //   if(detectorNumbers[iEv] == hb_detector)
      //   {
      //     detectorIsThere = true;
      //   }
      // }
      // if(!detectorIsThere)
      // {
      //   detectorNumbers.push_back(hb_detector);
      // }


      // counter++;
      // if((Nevents/10) != 0)
      // {
      //   if((counter % (Nevents/10))==0)
      //   {
      //     if(Nevents != 0)
      //     {
      //       std::cout << ((100*counter)/Nevents) << "% done... " << std::endl;
      //     }
      //   }
      // }

    }
    std::cout << std::endl;

    Nevents = eventsNumbers.size();
    std::cout << "Found " << Nevents <<  " events." << std::endl;
    // Ndet = detectorNumbers.size();
    // std::cout << "Found " << Ndet <<  " detectors." << std::endl;
  }



  // find modules and detectors per module
  // get modules ttree, count modules and types
  TTree *modules = (TTree*) fileIn->Get("modules");
  std::vector<region_t> modules_regions = CalculateRegions(modules);
  std::vector<int> moduleTypes;
  // find how many module types are involved
  for(int iMod = 0 ; iMod < modules_regions.size(); iMod++)
  {
    bool moduleTypeIsAlreadyThere = false;
    for(unsigned int im = 0; im < moduleTypes.size(); im++)
    {
      if(modules_regions[iMod].type == moduleTypes[im])
      {
        moduleTypeIsAlreadyThere = true;
      }
    }
    if(!moduleTypeIsAlreadyThere)
    {
      moduleTypes.push_back(modules_regions[iMod].type);
    }
  }


  int Nmodules  = modules_regions.size();
  int NmodTypes = moduleTypes.size();

  int* nModulePerType = new int[NmodTypes];

  int *nFrontDetPerType = new int[NmodTypes];
  int *nBackDetPerType  = new int[NmodTypes];
  for(int i = 0; i < NmodTypes; i++)
  {
    std::stringstream snameFront;
    std::stringstream snameBack;
    snameFront  << "cell_map_front_" << moduleTypes[i];
    snameBack   << "cell_map_back_"  << moduleTypes[i];
    nFrontDetPerType[i] = ((TH2I*) fileIn->Get(snameFront.str().c_str()))->GetEntries();
    nBackDetPerType [i] = ((TH2I*) fileIn->Get(snameBack .str().c_str()))->GetEntries();
  }

  feedbackString << std::endl;
  feedbackString << "EVENTS AND DETECTORS:" << std::endl;
  feedbackString << "Events                      = " << Nevents   << std::endl;
  feedbackString << "Modules                     = " << Nmodules  << std::endl;
  feedbackString << "Module types                = " << NmodTypes << std::endl;
  for(unsigned int im = 0; im < moduleTypes.size(); im++)
  {
    feedbackString << "Module Type n " << moduleTypes[im] << " has " << nFrontDetPerType[im] << " modules front" << std::endl;
    feedbackString << "Module Type n " << moduleTypes[im] << " has " << nBackDetPerType[im] << " modules back" << std::endl;
  }

  // feedbackString << "Detectors                   = " << Ndet <<std::endl;
  feedbackString << std::endl;
  std::cout << feedbackString.str();

  // array of events
  event_t *event = new event_t[Nevents];

  for(int iEv = 0; iEv < Nevents; iEv++)
  {
    // for each event
    event[iEv].primarySet = false;
    // array of modules
    event[iEv].module = new module_t[Nmodules];
    for(int iMod = 0 ; iMod < Nmodules; iMod++)
    {
      int typeID = -1;
      for(unsigned int im = 0; im < moduleTypes.size(); im++)
      {
        if(modules_regions[iMod].type == moduleTypes[im])
        {
          typeID = im;
        }
      }
      // prepare det arrays
      event[iEv].module[iMod].detector = new detector_t*[2];

      event[iEv].module[iMod].detector[0] = new detector_t[nFrontDetPerType[typeID]];
      event[iEv].module[iMod].detector[1] = new detector_t[nBackDetPerType[typeID]];
      for(int iDet = 0; iDet < nFrontDetPerType[typeID]; iDet++)
      {
        event[iEv].module[iMod].detector[0][iDet].timestamp = 0.;
        event[iEv].module[iMod].detector[0][iDet].counts = 0;
      }
      for(int iDet = 0; iDet < nBackDetPerType[typeID]; iDet++)
      {
        event[iEv].module[iMod].detector[1][iDet].timestamp = 0.;
        event[iEv].module[iMod].detector[1][iDet].counts = 0;
      }
    }
  }

  // output file needs to be opened before creating output ttree...
  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  //-------------------//
  // Output TTree      //
  //-------------------//
  ULong64_t *counts;
  counts = new ULong64_t[Ndet];
  Float_t *timestamp;
  Float_t *t_dep;
  Float_t *t_gen;
  Float_t *t_prop;
  Float_t *t_dep_gen;
  Float_t *t_dep_prop;
  Float_t *t_gen_prop;
  Float_t *u_timestamp;
  Float_t *u_t_dep;
  Float_t *u_t_gen;
  Float_t *u_t_prop;
  Float_t *u_t_dep_gen;
  Float_t *u_t_dep_prop;
  Float_t *u_t_gen_prop;

  timestamp  = new Float_t[Ndet];
  t_dep      = new Float_t[Ndet];
  t_gen      = new Float_t[Ndet];
  t_prop     = new Float_t[Ndet];
  t_dep_gen  = new Float_t[Ndet];
  t_dep_prop = new Float_t[Ndet];
  t_gen_prop = new Float_t[Ndet];
  u_timestamp  = new Float_t[Ndet];
  u_t_dep      = new Float_t[Ndet];
  u_t_gen      = new Float_t[Ndet];
  u_t_prop     = new Float_t[Ndet];
  u_t_dep_gen  = new Float_t[Ndet];
  u_t_dep_prop = new Float_t[Ndet];
  u_t_gen_prop = new Float_t[Ndet];

  TTree* outTree = new TTree("tree","tree");
  //branches of the channels data
  std::stringstream snames,stypes;
  for (int i = 0 ; i < Ndet ; i++)
  {
    //empty the stringstreams
    snames.str("");
    stypes.str("");
    counts[i] = 0;
    snames << "ch" << i;
    stypes << "ch" << i << "/l";
    outTree->Branch(snames.str().c_str(),&counts[i],stypes.str().c_str());
    snames.str("");
    stypes.str("");
  }
  for (int i = 0 ; i < Ndet ; i++)
  {
    // smeared
    //empty the stringstreams
    snames.str("");
    stypes.str("");
    timestamp[i] = 0;
    snames << "t" << i;
    stypes << "t" << i << "/F";
    outTree->Branch(snames.str().c_str(),&timestamp[i],stypes.str().c_str());
    snames.str("");
    stypes.str("");

    //empty the stringstreams
    snames.str("");
    stypes.str("");
    t_dep[i] = 0;
    snames << "t_dep" << i;
    stypes << "t_dep" << i << "/F";
    outTree->Branch(snames.str().c_str(),&t_dep[i],stypes.str().c_str());
    snames.str("");
    stypes.str("");

    //empty the stringstreams
    snames.str("");
    stypes.str("");
    t_gen[i] = 0;
    snames << "t_gen" << i;
    stypes << "t_gen" << i << "/F";
    outTree->Branch(snames.str().c_str(),&t_gen[i],stypes.str().c_str());
    snames.str("");
    stypes.str("");

    //empty the stringstreams
    snames.str("");
    stypes.str("");
    t_prop[i] = 0;
    snames << "t_prop" << i;
    stypes << "t_prop" << i << "/F";
    outTree->Branch(snames.str().c_str(),&t_prop[i],stypes.str().c_str());
    snames.str("");
    stypes.str("");

    //empty the stringstreams
    snames.str("");
    stypes.str("");
    t_dep_gen[i] = 0;
    snames << "t_dep_gen" << i;
    stypes << "t_dep_gen" << i << "/F";
    outTree->Branch(snames.str().c_str(),&t_dep_gen[i],stypes.str().c_str());
    snames.str("");
    stypes.str("");

    //empty the stringstreams
    snames.str("");
    stypes.str("");
    t_dep_prop[i] = 0;
    snames << "t_dep_prop" << i;
    stypes << "t_dep_prop" << i << "/F";
    outTree->Branch(snames.str().c_str(),&t_dep_prop[i],stypes.str().c_str());
    snames.str("");
    stypes.str("");

    //empty the stringstreams
    snames.str("");
    stypes.str("");
    t_gen_prop[i] = 0;
    snames << "t_gen_prop" << i;
    stypes << "t_gen_prop" << i << "/F";
    outTree->Branch(snames.str().c_str(),&t_gen_prop[i],stypes.str().c_str());
    snames.str("");
    stypes.str("");

    // smeared
    //empty the stringstreams
    snames.str("");
    stypes.str("");
    u_timestamp[i] = 0;
    snames << "u_t" << i;
    stypes << "u_t" << i << "/F";
    outTree->Branch(snames.str().c_str(),&u_timestamp[i],stypes.str().c_str());
    snames.str("");
    stypes.str("");

    //empty the stringstreams
    snames.str("");
    stypes.str("");
    u_t_dep[i] = 0;
    snames << "u_t_dep" << i;
    stypes << "u_t_dep" << i << "/F";
    outTree->Branch(snames.str().c_str(),&u_t_dep[i],stypes.str().c_str());
    snames.str("");
    stypes.str("");

    //empty the stringstreams
    snames.str("");
    stypes.str("");
    u_t_gen[i] = 0;
    snames << "u_t_gen" << i;
    stypes << "u_t_gen" << i << "/F";
    outTree->Branch(snames.str().c_str(),&u_t_gen[i],stypes.str().c_str());
    snames.str("");
    stypes.str("");

    //empty the stringstreams
    snames.str("");
    stypes.str("");
    u_t_prop[i] = 0;
    snames << "u_t_prop" << i;
    stypes << "u_t_prop" << i << "/F";
    outTree->Branch(snames.str().c_str(),&u_t_prop[i],stypes.str().c_str());
    snames.str("");
    stypes.str("");

    //empty the stringstreams
    snames.str("");
    stypes.str("");
    u_t_dep_gen[i] = 0;
    snames << "u_t_dep_gen" << i;
    stypes << "u_t_dep_gen" << i << "/F";
    outTree->Branch(snames.str().c_str(),&u_t_dep_gen[i],stypes.str().c_str());
    snames.str("");
    stypes.str("");

    //empty the stringstreams
    snames.str("");
    stypes.str("");
    u_t_dep_prop[i] = 0;
    snames << "u_t_dep_prop" << i;
    stypes << "u_t_dep_prop" << i << "/F";
    outTree->Branch(snames.str().c_str(),&u_t_dep_prop[i],stypes.str().c_str());
    snames.str("");
    stypes.str("");

    //empty the stringstreams
    snames.str("");
    stypes.str("");
    u_t_gen_prop[i] = 0;
    snames << "u_t_gen_prop" << i;
    stypes << "u_t_gen_prop" << i << "/F";
    outTree->Branch(snames.str().c_str(),&u_t_gen_prop[i],stypes.str().c_str());
    snames.str("");
    stypes.str("");



  }
  snames.str("");
  stypes.str("");
  // save primary info
  //----------------------------------//
  float primaryPositionOnAbsorberX    ;
  float primaryPositionOnAbsorberY    ;
  float primaryPositionOnAbsorberZ    ;
  float primaryMomentumOnAbsorberX    ;
  float primaryMomentumOnAbsorberY    ;
  float primaryMomentumOnAbsorberZ    ;
  float primaryEnergyOnAbsorber       ;
  outTree->Branch("primary_PositionOnAbsorberX", &primaryPositionOnAbsorberX,"primary_PositionOnAbsorberX/F");
  outTree->Branch("primary_PositionOnAbsorberY", &primaryPositionOnAbsorberY,"primary_PositionOnAbsorberY/F");
  outTree->Branch("primary_PositionOnAbsorberZ", &primaryPositionOnAbsorberZ,"primary_PositionOnAbsorberZ/F");
  outTree->Branch("primary_MomentumOnAbsorberX", &primaryMomentumOnAbsorberX,"primary_MomentumOnAbsorberX/F");
  outTree->Branch("primary_MomentumOnAbsorberY", &primaryMomentumOnAbsorberY,"primary_MomentumOnAbsorberY/F");
  outTree->Branch("primary_MomentumOnAbsorberZ", &primaryMomentumOnAbsorberZ,"primary_MomentumOnAbsorberZ/F");
  outTree->Branch("primary_EnergyOnAbsorber"   , &primaryEnergyOnAbsorber   ,"primary_EnergyOnAbsorber/F");

  //----------------------------------//
  // LOOP ON HYBRID                   //
  //----------------------------------//
  counter = 0;
  std::cout << "Loop on events..." << std::endl;
  for(int i = 0 ; i < Nhybrid; i++)
  {
    hybrid->GetEntry(i);

    // write position in event

    if(event[hb_event].primarySet == false)
    {
      event[hb_event].primaryPositionOnAbsorberX = hb_primaryPositionOnAbsorberX   ;
      event[hb_event].primaryPositionOnAbsorberY = hb_primaryPositionOnAbsorberY   ;
      event[hb_event].primaryPositionOnAbsorberZ = hb_primaryPositionOnAbsorberZ   ;
      event[hb_event].primaryMomentumOnAbsorberX = hb_primaryMomentumOnAbsorberX   ;
      event[hb_event].primaryMomentumOnAbsorberY = hb_primaryMomentumOnAbsorberY   ;
      event[hb_event].primaryMomentumOnAbsorberZ = hb_primaryMomentumOnAbsorberZ   ;
      event[hb_event].primaryEnergyOnAbsorber    = hb_primaryEnergyOnAbsorber      ;
      event[hb_event].primarySet = true;
    }
    // add to the event and the detector written in ttree, but passing by QE and smearing by sptr
    double dice = rand->Uniform(1.0);
    // std::cout << i << " " << hb_event << " " << hb_detector << std::endl;
    if(dice < qe)
    {
      // find front_back (front = 0 back = 1)
      int fb = -1;
      if(hb_z > 0) fb = 0;
      if(hb_z < 0) fb = 1;
      event[hb_event].module[hb_module].detector[fb][hb_detector].counts++;
      // smear, unless sigmaSPTR = 0

      // full
      Float_t smearedTimeStamp;
      if(sigmaSPTR == 0)
      {
        smearedTimeStamp = hb_timestamp;
      }
      else
      {
        smearedTimeStamp = (Float_t) rand->Gaus(hb_timestamp,sigmaSPTR);
      }
      event[hb_event].module[hb_module].detector[fb][hb_detector].listOfTimestamps.push_back(smearedTimeStamp);

      // dep
      if(sigmaSPTR == 0)
      {
        smearedTimeStamp = hb_t_deposition;
      }
      else
      {
        smearedTimeStamp = (Float_t) rand->Gaus(hb_t_deposition,sigmaSPTR);
      }
      event[hb_event].module[hb_module].detector[fb][hb_detector].listOft_dep.push_back(smearedTimeStamp);

      // gen
      if(sigmaSPTR == 0)
      {
        smearedTimeStamp = hb_t_generation;
      }
      else
      {
        smearedTimeStamp = (Float_t) rand->Gaus(hb_t_generation,sigmaSPTR);
      }
      event[hb_event].module[hb_module].detector[fb][hb_detector].listOft_gen.push_back(smearedTimeStamp);

      // prop
      if(sigmaSPTR == 0)
      {
        smearedTimeStamp = hb_t_propagation;
      }
      else
      {
        smearedTimeStamp = (Float_t) rand->Gaus(hb_t_propagation,sigmaSPTR);
      }
      event[hb_event].module[hb_module].detector[fb][hb_detector].listOft_prop.push_back(smearedTimeStamp);

      // dep + gen
      if(sigmaSPTR == 0)
      {
        smearedTimeStamp = hb_t_deposition+hb_t_generation;
      }
      else
      {
        smearedTimeStamp = (Float_t) rand->Gaus(hb_t_deposition+hb_t_generation,sigmaSPTR);
      }
      event[hb_event].module[hb_module].detector[fb][hb_detector].listOft_dep_gen.push_back(smearedTimeStamp);

      // dep + prop
      if(sigmaSPTR == 0)
      {
        smearedTimeStamp = hb_t_deposition+hb_t_propagation;
      }
      else
      {
        smearedTimeStamp = (Float_t) rand->Gaus(hb_t_deposition+hb_t_propagation,sigmaSPTR);
      }
      event[hb_event].module[hb_module].detector[fb][hb_detector].listOft_dep_prop.push_back(smearedTimeStamp);

      // gen + prop
      if(sigmaSPTR == 0)
      {
        smearedTimeStamp = hb_t_generation+hb_t_propagation;
      }
      else
      {
        smearedTimeStamp = (Float_t) rand->Gaus(hb_t_generation+hb_t_propagation,sigmaSPTR);
      }
      event[hb_event].module[hb_module].detector[fb][hb_detector].listOft_gen_prop.push_back(smearedTimeStamp);

      // no smear
      smearedTimeStamp = hb_timestamp;
      event[hb_event].module[hb_module].detector[fb][hb_detector].u_listOfTimestamps.push_back(smearedTimeStamp);
      // dep
      smearedTimeStamp = (Float_t) hb_t_deposition;
      event[hb_event].module[hb_module].detector[fb][hb_detector].u_listOft_dep.push_back(smearedTimeStamp);
      // gen
      smearedTimeStamp = hb_t_generation;
      event[hb_event].module[hb_module].detector[fb][hb_detector].u_listOft_gen.push_back(smearedTimeStamp);
      // prop
      smearedTimeStamp = hb_t_propagation;
      event[hb_event].module[hb_module].detector[fb][hb_detector].u_listOft_prop.push_back(smearedTimeStamp);
      // dep + gen
      smearedTimeStamp = hb_t_deposition+hb_t_generation;
      event[hb_event].module[hb_module].detector[fb][hb_detector].u_listOft_dep_gen.push_back(smearedTimeStamp);
      // dep + prop
      smearedTimeStamp = hb_t_deposition+hb_t_propagation;
      event[hb_event].module[hb_module].detector[fb][hb_detector].u_listOft_dep_prop.push_back(smearedTimeStamp);
      // gen + prop
      smearedTimeStamp = hb_t_generation+hb_t_propagation;
      event[hb_event].module[hb_module].detector[fb][hb_detector].u_listOft_gen_prop.push_back(smearedTimeStamp);

    }


    // counter++;
    // if((Nevents/10) != 0)
    // {
    //   if((counter % (Nevents/10))==0)
    //   {
    //     if(Nevents != 0)
    //     {
    //       std::cout << ((100*counter)/Nevents) << "% done... " << std::endl;
    //     }
    //   }
    // }

  }
  std::cout <<  std::endl;

  std::cout << "Calculating timestamps..." << std::endl;
  counter = 0;
  // run on the events, calc average timestamps
  for(int iEv = 0; iEv < Nevents; iEv++)
  {
    if(verbose == 2) std::cout << "Event " << iEv << std::endl;
    // set primaries
    primaryPositionOnAbsorberX = event[iEv].primaryPositionOnAbsorberX;
    primaryPositionOnAbsorberY = event[iEv].primaryPositionOnAbsorberY;
    primaryPositionOnAbsorberZ = event[iEv].primaryPositionOnAbsorberZ;
    primaryMomentumOnAbsorberX = event[iEv].primaryMomentumOnAbsorberX;
    primaryMomentumOnAbsorberY = event[iEv].primaryMomentumOnAbsorberY;
    primaryMomentumOnAbsorberZ = event[iEv].primaryMomentumOnAbsorberZ;
    primaryEnergyOnAbsorber    = event[iEv].primaryEnergyOnAbsorber   ;
    // if(verbose == 2) std::cout << "------------------> "<< iEv << " ";
    // run on each detector

    for(int iMod = 0 ; iMod < Nmodules; iMod++)
    {
      int typeID = -1;
      for(unsigned int im = 0; im < moduleTypes.size(); im++)
      {
        if(modules_regions[iMod].type == moduleTypes[im])
        {
          typeID = im;
        }
      }
      for(int fb = 0 ; fb < 2; fb++)
      {
        for(int iDet = 0; iDet < nFrontDetPerType[typeID]; iDet++) // should be the same number of det front and back, come on..
        {


        }
      }



    }
    for(int iDet = 0; iDet < Ndet; iDet++)
    {
      // std::cout << iDet << " ";
      // store charge
      counts[iDet] = event[iEv].detector[iDet].counts;
      if(verbose == 1) std::cout << iDet << std::endl;
      //--------------------------------//
      // SMEARED
      //--------------------------------//
      // if timestamps are there..
      // event[hb_event].module[hb_module].detector[fb][hb_detector]

      if(event[iEv].detector[iDet].listOfTimestamps.size())
      {
        // ...calculate the detector timestamp
        // sort timestamps
        if(verbose == 2) std::cout << "timestamp" << std::endl;
        std::sort(event[iEv].detector[iDet].listOfTimestamps.begin(),event[iEv].detector[iDet].listOfTimestamps.end());
        // take the average of the first N
        // check if there are enough photons--
        int effectiveN = numb_of_phot_for_time_average;

        if(numb_of_phot_for_time_average > event[iEv].detector[iDet].listOfTimestamps.size())
        {
          effectiveN = event[iEv].detector[iDet].listOfTimestamps.size();
        }
        // std::cout << effectiveN << " ";
        // calc timestamp
        Float_t final_timestamp = 0.0;
        for(int j = 0 ; j < effectiveN; j++)
        {
          final_timestamp +=  (Float_t) ((event[iEv].detector[iDet].listOfTimestamps[j] / effectiveN)); // avg
          // std::cout << event[iEv].detector[iDet].listOfTimestamps[j] << " " ;
        }
        // std::cout << final_timestamp << std::endl;
        event[iEv].detector[iDet].timestamp = final_timestamp;
      }
      else
      {
        event[iEv].detector[iDet].timestamp = -1000;  // FIXME - for now put unreasonably in the past, so it will be discarded. find a more elegant and maintaniable way!
      }
      // and store it
      timestamp[iDet] = event[iEv].detector[iDet].timestamp;


      // if timestamps are there..
      if(event[iEv].detector[iDet].listOft_dep.size())
      {
        // ...calculate the detector timestamp
        // sort timestamps
        if(verbose == 2) std::cout << "t_dep" << std::endl;
        std::sort(event[iEv].detector[iDet].listOft_dep.begin(),event[iEv].detector[iDet].listOft_dep.end());
        // take the average of the first N
        // check if there are enough photons--
        int effectiveN = numb_of_phot_for_time_average;
        if(numb_of_phot_for_time_average > event[iEv].detector[iDet].listOft_dep.size())
        {
          effectiveN = event[iEv].detector[iDet].listOft_dep.size();
        }
        // calc timestamp
        Float_t final_timestamp = 0.0;
        for(int j = 0 ; j < effectiveN; j++)
        {
          final_timestamp +=  (Float_t) ((event[iEv].detector[iDet].listOft_dep[j] / effectiveN)); // avg
        }
        event[iEv].detector[iDet].t_dep = final_timestamp;
      }
      else
      {
        event[iEv].detector[iDet].t_dep = -1000;
      }
      // and store it
      t_dep[iDet] = event[iEv].detector[iDet].t_dep;
      // end ----------


      // if timestamps are there..
      if(event[iEv].detector[iDet].listOft_gen.size())
      {
        // ...calculate the detector timestamp
        // sort timestamps
        if(verbose == 2) std::cout << "t_gen" << std::endl;
        std::sort(event[iEv].detector[iDet].listOft_gen.begin(),event[iEv].detector[iDet].listOft_gen.end());
        // take the average of the first N
        // check if there are enough photons--
        int effectiveN = numb_of_phot_for_time_average;
        if(numb_of_phot_for_time_average > event[iEv].detector[iDet].listOft_gen.size())
        {
          effectiveN = event[iEv].detector[iDet].listOft_gen.size();
        }
        // calc timestamp
        Float_t final_timestamp = 0.0;
        for(int j = 0 ; j < effectiveN; j++)
        {
          final_timestamp +=  (Float_t) ((event[iEv].detector[iDet].listOft_gen[j] / effectiveN)); // avg
        }
        event[iEv].detector[iDet].t_gen = final_timestamp;
      }
      else
      {
        event[iEv].detector[iDet].t_gen = -1000;
      }
      // and store it
      t_gen[iDet] = event[iEv].detector[iDet].t_gen;
      // end ----------


      // if timestamps are there..
      if(event[iEv].detector[iDet].listOft_prop.size())
      {
        // ...calculate the detector timestamp
        // sort timestamps
        if(verbose == 2) std::cout << "t_prop" << std::endl;

        std::sort(event[iEv].detector[iDet].listOft_prop.begin(),event[iEv].detector[iDet].listOft_prop.end());
        // take the average of the first N
        // check if there are enough photons--
        int effectiveN = numb_of_phot_for_time_average;
        if(numb_of_phot_for_time_average > event[iEv].detector[iDet].listOft_prop.size())
        {
          effectiveN = event[iEv].detector[iDet].listOft_prop.size();
        }
        // calc timestamp
        Float_t final_timestamp = 0.0;
        for(int j = 0 ; j < effectiveN; j++)
        {
          final_timestamp +=  (Float_t) ((event[iEv].detector[iDet].listOft_prop[j] / effectiveN)); // avg
        }
        event[iEv].detector[iDet].t_prop = final_timestamp;
      }
      else
      {
        event[iEv].detector[iDet].t_prop = -1000;
      }
      // and store it
      t_prop[iDet] = event[iEv].detector[iDet].t_prop;
      // end ----------


      // if timestamps are there..
      if(event[iEv].detector[iDet].listOft_dep_gen.size())
      {
        // ...calculate the detector timestamp
        // sort timestamps
        if(verbose == 2) std::cout << "t_dep_gen" << std::endl;

        std::sort(event[iEv].detector[iDet].listOft_dep_gen.begin(),event[iEv].detector[iDet].listOft_dep_gen.end());
        // take the average of the first N
        // check if there are enough photons--
        int effectiveN = numb_of_phot_for_time_average;
        if(numb_of_phot_for_time_average > event[iEv].detector[iDet].listOft_dep_gen.size())
        {
          effectiveN = event[iEv].detector[iDet].listOft_dep_gen.size();
        }
        // calc timestamp
        Float_t final_timestamp = 0.0;
        for(int j = 0 ; j < effectiveN; j++)
        {
          final_timestamp +=  (Float_t) ((event[iEv].detector[iDet].listOft_dep_gen[j] / effectiveN)); // avg
        }
        event[iEv].detector[iDet].t_dep_gen = final_timestamp;
      }
      else
      {
        event[iEv].detector[iDet].t_dep_gen = -1000;
      }
      // and store it
      t_dep_gen[iDet] = event[iEv].detector[iDet].t_dep_gen;
      // end ----------


      // if timestamps are there..
      if(event[iEv].detector[iDet].listOft_dep_prop.size())
      {
        // ...calculate the detector timestamp
        // sort timestamps
        if(verbose == 2) std::cout << "t_dep_prop" << std::endl;
        std::sort(event[iEv].detector[iDet].listOft_dep_prop.begin(),event[iEv].detector[iDet].listOft_dep_prop.end());
        // take the average of the first N
        // check if there are enough photons--
        int effectiveN = numb_of_phot_for_time_average;
        if(numb_of_phot_for_time_average > event[iEv].detector[iDet].listOft_dep_prop.size())
        {
          effectiveN = event[iEv].detector[iDet].listOft_dep_prop.size();
        }
        // calc timestamp
        Float_t final_timestamp = 0.0;
        for(int j = 0 ; j < effectiveN; j++)
        {
          final_timestamp +=  (Float_t) ((event[iEv].detector[iDet].listOft_dep_prop[j] / effectiveN)); // avg
        }
        event[iEv].detector[iDet].t_dep_prop = final_timestamp;
      }
      else
      {
        event[iEv].detector[iDet].t_dep_prop = -1000;
      }
      // and store it
      t_dep_prop[iDet] = event[iEv].detector[iDet].t_dep_prop;
      // end ----------

      // if timestamps are there..
      if(event[iEv].detector[iDet].listOft_gen_prop.size())
      {
        // ...calculate the detector timestamp
        // sort timestamps
        if(verbose == 2) std::cout << "t_gen_prop" << std::endl;

        std::sort(event[iEv].detector[iDet].listOft_gen_prop.begin(),event[iEv].detector[iDet].listOft_gen_prop.end());
        // take the average of the first N
        // check if there are enough photons--
        int effectiveN = numb_of_phot_for_time_average;
        if(numb_of_phot_for_time_average > event[iEv].detector[iDet].listOft_gen_prop.size())
        {
          effectiveN = event[iEv].detector[iDet].listOft_gen_prop.size();
        }
        // calc timestamp
        Float_t final_timestamp = 0.0;
        for(int j = 0 ; j < effectiveN; j++)
        {
          final_timestamp +=  (Float_t) ((event[iEv].detector[iDet].listOft_gen_prop[j] / effectiveN)); // avg
        }
        event[iEv].detector[iDet].t_gen_prop = final_timestamp;
      }
      else
      {
        event[iEv].detector[iDet].t_gen_prop = -1000;
      }
      // and store it
      t_gen_prop[iDet] = event[iEv].detector[iDet].t_gen_prop;
      // end ----------

      //------------------------------//

      //--------------------------------//
      // UNSMEARED
      //--------------------------------//
      // if timestamps are there..
      if(event[iEv].detector[iDet].u_listOfTimestamps.size())
      {
        // ...calculate the detector timestamp
        // sort timestamps
        if(verbose == 2) std::cout << "u_timestamp" << std::endl;

        std::sort(event[iEv].detector[iDet].u_listOfTimestamps.begin(),event[iEv].detector[iDet].u_listOfTimestamps.end());
        // take the average of the first N
        // check if there are enough photons--
        int effectiveN = numb_of_phot_for_time_average;

        if(numb_of_phot_for_time_average > event[iEv].detector[iDet].u_listOfTimestamps.size())
        {
          effectiveN = event[iEv].detector[iDet].u_listOfTimestamps.size();
        }
        // std::cout << effectiveN << " ";
        // calc timestamp
        Float_t final_timestamp = 0.0;
        for(int j = 0 ; j < effectiveN; j++)
        {
          final_timestamp +=  (Float_t) ((event[iEv].detector[iDet].u_listOfTimestamps[j] / effectiveN)); // avg
          // std::cout << event[iEv].detector[iDet].listOfTimestamps[j] << " " ;
        }
        // std::cout << final_timestamp << std::endl;
        event[iEv].detector[iDet].u_timestamp = final_timestamp;
      }
      else
      {
        event[iEv].detector[iDet].u_timestamp = -1000;  // FIXME - for now put unreasonably in the past, so it will be discarded. find a more elegant and maintaniable way!
      }
      // and store it
      u_timestamp[iDet] = event[iEv].detector[iDet].u_timestamp;


      // if timestamps are there..
      if(event[iEv].detector[iDet].u_listOft_dep.size())
      {
        // ...calculate the detector timestamp
        // sort timestamps
        if(verbose == 2) std::cout << "u_t_dep" << std::endl;

        std::sort(event[iEv].detector[iDet].u_listOft_dep.begin(),event[iEv].detector[iDet].u_listOft_dep.end());
        // take the average of the first N
        // check if there are enough photons--
        int effectiveN = numb_of_phot_for_time_average;
        if(numb_of_phot_for_time_average > event[iEv].detector[iDet].u_listOft_dep.size())
        {
          effectiveN = event[iEv].detector[iDet].u_listOft_dep.size();
        }
        // calc timestamp
        Float_t final_timestamp = 0.0;
        for(int j = 0 ; j < effectiveN; j++)
        {
          final_timestamp +=  (Float_t) ((event[iEv].detector[iDet].u_listOft_dep[j] / effectiveN)); // avg
        }
        event[iEv].detector[iDet].u_t_dep = final_timestamp;
      }
      else
      {
        event[iEv].detector[iDet].u_t_dep = -1000;
      }
      // and store it
      u_t_dep[iDet] = event[iEv].detector[iDet].u_t_dep;
      // end ----------


      // if timestamps are there..
      if(event[iEv].detector[iDet].u_listOft_gen.size())
      {
        // ...calculate the detector timestamp
        // sort timestamps
        if(verbose == 2) std::cout << "u_t_gen" << std::endl;

        std::sort(event[iEv].detector[iDet].u_listOft_gen.begin(),event[iEv].detector[iDet].u_listOft_gen.end());
        // take the average of the first N
        // check if there are enough photons--
        int effectiveN = numb_of_phot_for_time_average;
        if(numb_of_phot_for_time_average > event[iEv].detector[iDet].u_listOft_gen.size())
        {
          effectiveN = event[iEv].detector[iDet].u_listOft_gen.size();
        }
        // calc timestamp
        Float_t final_timestamp = 0.0;
        for(int j = 0 ; j < effectiveN; j++)
        {
          final_timestamp +=  (Float_t) ((event[iEv].detector[iDet].u_listOft_gen[j] / effectiveN)); // avg
        }
        event[iEv].detector[iDet].u_t_gen = final_timestamp;
      }
      else
      {
        event[iEv].detector[iDet].u_t_gen = -1000;
      }
      // and store it
      u_t_gen[iDet] = event[iEv].detector[iDet].u_t_gen;
      // end ----------


      // if timestamps are there..
      if(event[iEv].detector[iDet].u_listOft_prop.size())
      {
        // ...calculate the detector timestamp
        // sort timestamps
        if(verbose == 2) std::cout << "u_t_prop" << std::endl;

        std::sort(event[iEv].detector[iDet].u_listOft_prop.begin(),event[iEv].detector[iDet].u_listOft_prop.end());
        // take the average of the first N
        // check if there are enough photons--
        int effectiveN = numb_of_phot_for_time_average;
        if(numb_of_phot_for_time_average > event[iEv].detector[iDet].u_listOft_prop.size())
        {
          effectiveN = event[iEv].detector[iDet].u_listOft_prop.size();
        }
        // calc timestamp
        Float_t final_timestamp = 0.0;
        for(int j = 0 ; j < effectiveN; j++)
        {
          final_timestamp +=  (Float_t) ((event[iEv].detector[iDet].u_listOft_prop[j] / effectiveN)); // avg
        }
        event[iEv].detector[iDet].u_t_prop = final_timestamp;
      }
      else
      {
        event[iEv].detector[iDet].u_t_prop = -1000;
      }
      // and store it
      u_t_prop[iDet] = event[iEv].detector[iDet].u_t_prop;
      // end ----------


      // if timestamps are there..
      if(event[iEv].detector[iDet].u_listOft_dep_gen.size())
      {
        // ...calculate the detector timestamp
        // sort timestamps
        if(verbose == 2) std::cout << "u_t_dep_gen" << std::endl;

        std::sort(event[iEv].detector[iDet].u_listOft_dep_gen.begin(),event[iEv].detector[iDet].u_listOft_dep_gen.end());
        // take the average of the first N
        // check if there are enough photons--
        int effectiveN = numb_of_phot_for_time_average;
        if(numb_of_phot_for_time_average > event[iEv].detector[iDet].u_listOft_dep_gen.size())
        {
          effectiveN = event[iEv].detector[iDet].u_listOft_dep_gen.size();
        }
        // calc timestamp
        Float_t final_timestamp = 0.0;
        for(int j = 0 ; j < effectiveN; j++)
        {
          final_timestamp +=  (Float_t) ((event[iEv].detector[iDet].u_listOft_dep_gen[j] / effectiveN)); // avg
        }
        event[iEv].detector[iDet].u_t_dep_gen = final_timestamp;
      }
      else
      {
        event[iEv].detector[iDet].u_t_dep_gen = -1000;
      }
      // and store it
      u_t_dep_gen[iDet] = event[iEv].detector[iDet].u_t_dep_gen;
      // end ----------


      // if timestamps are there..
      if(event[iEv].detector[iDet].u_listOft_dep_prop.size())
      {
        // ...calculate the detector timestamp
        // sort timestamps
        if(verbose == 2) std::cout << "u_t_dep_prop" << std::endl;

        std::sort(event[iEv].detector[iDet].u_listOft_dep_prop.begin(),event[iEv].detector[iDet].u_listOft_dep_prop.end());
        // take the average of the first N
        // check if there are enough photons--
        int effectiveN = numb_of_phot_for_time_average;
        if(numb_of_phot_for_time_average > event[iEv].detector[iDet].u_listOft_dep_prop.size())
        {
          effectiveN = event[iEv].detector[iDet].u_listOft_dep_prop.size();
        }
        // calc timestamp
        Float_t final_timestamp = 0.0;
        for(int j = 0 ; j < effectiveN; j++)
        {
          final_timestamp +=  (Float_t) ((event[iEv].detector[iDet].u_listOft_dep_prop[j] / effectiveN)); // avg
        }
        event[iEv].detector[iDet].u_t_dep_prop = final_timestamp;
      }
      else
      {
        event[iEv].detector[iDet].u_t_dep_prop = -1000;
      }
      // and store it
      u_t_dep_prop[iDet] = event[iEv].detector[iDet].u_t_dep_prop;
      // end ----------

      // if timestamps are there..
      if(event[iEv].detector[iDet].u_listOft_gen_prop.size())
      {
        // ...calculate the detector timestamp
        // sort timestamps
        if(verbose == 2) std::cout << "u_t_gen_prop" << std::endl;

        std::sort(event[iEv].detector[iDet].u_listOft_gen_prop.begin(),event[iEv].detector[iDet].u_listOft_gen_prop.end());
        // take the average of the first N
        // check if there are enough photons--
        int effectiveN = numb_of_phot_for_time_average;
        if(numb_of_phot_for_time_average > event[iEv].detector[iDet].u_listOft_gen_prop.size())
        {
          effectiveN = event[iEv].detector[iDet].u_listOft_gen_prop.size();
        }
        // calc timestamp
        Float_t final_timestamp = 0.0;
        for(int j = 0 ; j < effectiveN; j++)
        {
          final_timestamp +=  (Float_t) ((event[iEv].detector[iDet].u_listOft_gen_prop[j] / effectiveN)); // avg
        }
        event[iEv].detector[iDet].u_t_gen_prop = final_timestamp;
      }
      else
      {
        event[iEv].detector[iDet].u_t_gen_prop = -1000;
      }
      // and store it
      u_t_gen_prop[iDet] = event[iEv].detector[iDet].u_t_gen_prop;
      // end ----------

      //------------------------------//



    }
    if(verbose == 1) std::cout << "fill " << iEv << std::endl;
    // save detectors data in output ttree
    outTree->Fill();

    if(verbose == 1) std::cout << "fill done " << iEv << std::endl;


    // counter++;
    // if((Nevents/10) != 0)
    // {
    //   if((counter % (Nevents/10))==0)
    //   {
    //     if(Nevents != 0)
    //     {
    //       std::cout << ((100*counter)/Nevents) << "% done... " << std::endl;
    //     }
    //   }
    // }

  }
  std::cout <<  std::endl;

  if(verbose == 1) std::cout << "done loop " << std::endl;


  std::cout << std::endl;
  std::cout << "Saving data to file " << outputFileName << std::endl;
  outputFile->cd();
  outTree->Write();
  if(verbose == 1) std::cout << "tree wrote " << std::endl;
  TDirectory *configDir = outputFile->mkdir("Configuration");
  configDir->cd();
  TNamed FeedbackNameD("Parameters",feedbackString.str().c_str());
  TNamed CommandLineNameD("CommandLine",CommandLineString.c_str());
  TNamed HostNameNameD("Hostname",HostNameString.c_str());
  TNamed PWDNameD("PWD",PWDstring.c_str());
  FeedbackNameD.Write();
  if(verbose == 1) std::cout << "feedback wrote " << std::endl;
  CommandLineNameD.Write();
  if(verbose == 1) std::cout << "cmd wrote " << std::endl;
  HostNameNameD.Write();
  if(verbose == 1) std::cout << "host wrote " << std::endl;
  PWDNameD.Write();
  if(verbose == 1) std::cout << "pwd wrote " << std::endl;
  outputFile->Close();

  std::cout << "Done." << std::endl;
  return 0;
}



void usage()
{
  std::cout << "\t\t" << "[-i | --input]          <input file>   " << std::endl
            << "\t\t" << "[-o | --output]         <output file name>    " << std::endl
            << "\t\t" << "[-s | --seed]           <seed for TRandom3>                  - default = 0 , i.e. seed automatically computed via a TUUID object" << std::endl
            << "\t\t" << "[--qe]                  <quantum efficiency>                 - default = 0.07 "  << std::endl
            << "\t\t" << "[ --photons]            <average time on first N photons>    - default = 5]"     << std::endl
            << "\t\t" << "[ --sptr]               <sigma sptr of the detector>         - default = 0.087]" << std::endl
            << "\t\t" << "[ --events]             <primary events in input file>       - default = -1 (if you leave the default, the program will calculate it, but it's slower)]"     << std::endl
            // << "\t\t" << "[ --modules]            <number of modules>                  - default = -1 (if you leave the default, the program will calculate it, but it's slower)]" << std::endl
            // << "\t\t" << "[ --detectors]          <number of detectors per module>     - default = -1 (if you leave the default, the program will calculate it, but it's slower)]" << std::endl
            << "\t\t" << "[ --verbose]            <verbose level>                      - default = 0]" << std::endl
            << "\t\t" << std::endl;
}
