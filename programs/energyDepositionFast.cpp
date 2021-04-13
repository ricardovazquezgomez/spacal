// Marco Pizzichemi 29.04.2020 marco.pizzichemi@cern.ch

// compile with
// g++ -o ../build/energyDepositionFast energyDepositionFast.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore

// calculates energy deposited in the crystals everywhere in the simulation
// the program retrieves the simulation geometry from the input file, runs on the entries of the shower ttree
// and sums the energies deposited the crystals
// for help on command line parameters, run the program without arguments

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
#include "TGraphDelaunay.h"
#include "TVector.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <getopt.h>
#include <limits>
#include <dirent.h>
#include <numeric>

#include "TEnv.h"

#include "../parametrization/regions.h"


// struct cell_t
// {
//   bool save;
//   TH1F *h_EnDepInCrystal   ;
//   // TH1F *h_EnDepNotInCrystal;
//   // TH1F *h_EnDepTotal       ;
// };

struct module_t
{
  bool save;
  int type;
  int sections;
  int *nCellsPerSections;
  TH1F *h_EnDepInCrystal   ;
  TH1F *h_EnDepNotInCrystal;
  TH1F *h_EnDepTotal       ;
  bool **saveCell;
  TH1F ***h_Cell_EnDepInCrystal   ;
  // cell_t **cell;
};

void usage()
{
  std::cout << "\t\t" << "[-f | --folder]  <input folder>   " << std::endl
            << "\t\t" << "[-i | --input]   <input file prefix>   " << std::endl
            << "\t\t" << "[-o | --output]  <output file name>   " << std::endl
            << "\t\t" << "[--type]         <css list of module types - default = 0 >  " << std::endl
            << "\t\t" << "[--verbose]      <verbosity level - default = 0>  " << std::endl
            << "\t\t" << "[--fraction]     <csv list of sampling fractions per module type - default = 1>  " << std::endl
            << "\t\t" << "[--sum-all]      <sum all energy deposits, not just in crystal. False if not given.>  " << std::endl
            << "\t\t" << std::endl;
}



int main(int argc, char** argv)
{

  //------------------------------------------------//
  // CHECK CMD LINE INPUT                           //
  //------------------------------------------------//
  if(argc < 2)
  {
    std::cout << argv[0] << std::endl;
    usage();
    return 1;
  }
  // end of CHECK CMD LINE INPUT
  //------------------------------------------------//

  gEnv->GetValue("TFile.Recover", 0);



  //------------------------------------------------//
  // PARSE INPUT                                    //
  //------------------------------------------------//
  // prepare variable, giving defaults
  std::string inputFilePrefix = "";
  std::string inputFolderName = "./";
  std::string outputFileName = "";
  // bool singleEvent = false;
  int verbose = 0;
  // int moduleTypeNumber = 0;
  // bool singleModuleType = false;
  
  bool multiType = false;
  bool multiFrac = false;
  std::vector<int> typeNumber;
  std::string typeString = "";
  std::vector<float> fracValue;
  std::string fracString = "";

  bool sumAll = false;

  static struct option longOptions[] =
  {
			{ "input", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { "folder", required_argument, 0, 0 },
      { "verbose", required_argument, 0, 0 },
      { "type", required_argument, 0, 0 },
      { "fraction", required_argument, 0, 0 },
      { "sum-all", no_argument, 0, 0 },
			{ NULL, 0, 0, 0 }
	};
  // read input
  while(1) {
		int optionIndex = 0;
		int c = getopt_long(argc, argv, "i:o:f:", longOptions, &optionIndex);
		if (c == -1) {
			break;
		}
		if (c == 'i'){
			inputFilePrefix= (char *)optarg;
    }
		else if (c == 'o'){
      outputFileName = (char *)optarg;
    }
    else if (c == 'f'){
      inputFolderName = (char *)optarg;
    }
		else if (c == 0 && optionIndex == 0){
      inputFilePrefix = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 1){
      outputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 2){
      inputFolderName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 3){
      verbose = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 4){
      typeString = (char *)optarg;
      multiType = true;
    }
    else if (c == 0 && optionIndex == 5){
      fracString = (char *)optarg;
      multiFrac = true;
    }
    else if (c == 0 && optionIndex == 6){
      sumAll = true;
    }
		else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
	}
  // check mandatory input
  if(inputFilePrefix == "")
  {
    std::cout << std::endl;
    std::cout << "ERROR! You need to provide an input file name!" << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }
  if(outputFileName == "")
  {
    std::cout << std::endl;
    std::cout << "ERROR! You need to provide an output file name!" << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }
  // frac and type
  if(multiFrac)
  {
    if(!multiType)
    {
      std::cout << std::endl;
      std::cout << "ERROR! If you speficify fraction(s) you need to specify an equal number of type(s)" << std::endl;
      std::cout << "See program usage below..." << std::endl;
      std::cout << std::endl;
      std::cout << argv[0];
      usage();
     return 1;
    }
  }
  if(!multiFrac)
  {
    if(multiType)
    {
      std::cout << std::endl;
      std::cout << "ERROR! If you speficify type(s) you need to specify an equal number of fraction(s)" << std::endl;
      std::cout << "See program usage below..." << std::endl;
      std::cout << std::endl;
      std::cout << argv[0];
      usage();
      return 1;
    }
  }

  // parse frac and types 
  if(multiType)
  {
    std::stringstream type_stream(typeString); //create string stream from the string
    while(type_stream.good()) {
      std::string substr;
      getline(type_stream, substr, ','); //get first string delimited by comma
      int t_number;
      std::istringstream ( substr ) >> t_number;
      typeNumber.push_back(t_number);
    }
  }
  
  if(multiFrac)
  {
    std::stringstream frac_stream(fracString); //create string stream from the string
    while(frac_stream.good()) {
      std::string substr;
      getline(frac_stream, substr, ','); //get first string delimited by comma
      float t_number;
      std::istringstream ( substr ) >> t_number;
      fracValue.push_back(t_number);
    }
  }

  if(typeNumber.size() != fracValue.size())
  {
    std::cout << std::endl;
      std::cout << "ERROR! If you speficify type(s) you need to specify an equal number of fraction(s)" << std::endl;
      std::cout << "See program usage below..." << std::endl;
      std::cout << std::endl;
      std::cout << argv[0];
      usage();
      return 1;
  }

  // put in a map, from crystal type to fraction
  std::map<int, float> fractionMap;
  bool usingFractionMap = true;
  
  if((typeNumber.size() == 0) && (fracValue.size() == 0))
  {
    usingFractionMap = false;
  }
  else 
  {
    for(unsigned int i = 0 ; i < typeNumber.size(); i++)
    {
      fractionMap.insert(std::make_pair(typeNumber[i],fracValue[i]));
    }
  }

  for(unsigned int i = 0 ; i < typeNumber.size(); i++)
  {
    std::cout << " type " << typeNumber[i] << std::endl;
  }

  for(unsigned int i = 0 ; i < fracValue.size(); i++)
  {
    std::cout << " fraction " << fracValue[i] << std::endl;
  }
  
  


  //----------------------------------//
  // INPUT FILES                   //
  //----------------------------------//
  std::vector<std::string> v;
  read_directory(inputFolderName, v);
  // extract files with correct prefix
  std::vector<std::string> listInputFiles;
  for(unsigned int i = 0 ; i < v.size() ; i++)
  {
    if(!v[i].compare(0,inputFilePrefix.size(),inputFilePrefix))
    {
      listInputFiles.push_back(inputFolderName + "/" + v[i]);
    }
  }
  // check if it's empty
  if(listInputFiles.size() == 0)
  {
    std::cout << std::endl;
    std::cout << "ERROR! Some input files do not exists! Aborting." << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }
  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");

  // calc regions
  std::stringstream feedbackString;
  feedbackString << "|***************************************************************|" << std::endl;
  feedbackString << "| READING SPACAL STRUCTURE                                      |" << std::endl;
  feedbackString << "|***************************************************************|" << std::endl;
  feedbackString << "" << std::endl;
  TFile *_file0 = TFile::Open(listInputFiles[0].c_str());
  // TFile *inputFile = TFile::Open(inputFileName.c_str());
  TTree *modules = (TTree *)_file0->Get("modules");
  TTree *absorbers = (TTree *)_file0->Get("absorbers");
  TTree *cells = (TTree *)_file0->Get("cells");
  TTree *holes = (TTree *)_file0->Get("holes");
  TTree *fibres = (TTree *)_file0->Get("fibres");
  // TTree *shower = (TTree *)_file0->Get("shower");
  _file0->cd("Configuration");
  TNamed SeedNameD("Seed", ((TNamed *)gDirectory->Get("Seed"))->GetTitle());
  TNamed HostNameD("Hostname", ((TNamed *)gDirectory->Get("Hostname"))->GetTitle());
  TNamed PWDNameD("PWD", ((TNamed *)gDirectory->Get("PWD"))->GetTitle());
  TNamed ConfigNameD("ConfigFile", ((TNamed *)gDirectory->Get("ConfigFile"))->GetTitle());
  TNamed GpsNameD("GpsFile", ((TNamed *)gDirectory->Get("GpsFile"))->GetTitle());
  TNamed PrimariesNameD("Primaries", ((TNamed *)gDirectory->Get("Primaries"))->GetTitle());
  _file0->cd();

  // produce the module ID TH2I map
  int det_counter = 0;
  TH2I *module_ID_map = ComputeElementMap(modules, "module_ID_map", -1, -1, det_counter);
  // and the detector id map
  // TH2I *detector_ID_map_front = ComputeDetectorMap(modules,cells,0);
  // TH2I *detector_ID_map_back  = ComputeDetectorMap(modules,cells,1);

  std::vector<region_t> modules_regions = CalculateRegions(modules);
  std::vector<int> moduleTypes;
  std::vector<int> moduleSections;
  std::vector<float> moduleSeparationZ;
  std::vector<float> modulePositionZ;
  // find how many module types are involved
  for (int iMod = 0; iMod < modules_regions.size(); iMod++)
  {
    bool moduleTypeIsAlreadyThere = false;
    for (unsigned int im = 0; im < moduleTypes.size(); im++)
    {
      if (modules_regions[iMod].type == moduleTypes[im])
      {
        moduleTypeIsAlreadyThere = true;
      }
    }
    if (!moduleTypeIsAlreadyThere)
    {
      moduleTypes.push_back(modules_regions[iMod].type);
      // also get sections
      moduleSections.push_back(modules_regions[iMod].sections);
      moduleSeparationZ.push_back(modules_regions[iMod].separation_z);
      modulePositionZ.push_back(modules_regions[iMod].position_z);
    }
  }

  // lots of maps
  std::map<int, int> moduleNumberMap;         // map from module type to number of modules
  std::map<int, int> moduleNumberOfBeforeMap; // map from module type to number of other modules
  std::map<int, int> moduleSectionsMap;       // map from module type to number of segments
  std::map<int, float> moduleMinMap;          // map from module type to min z
  std::map<int, float> moduleMaxMap;          // map from module type to max z
  std::map<int, float> moduleSeparationMap;   // map from module type to z of separation

  for (int i = 0; i < moduleTypes.size(); i++)
  {
    int modCounter = 0;
    for (int iMod = 0; iMod < modules_regions.size(); iMod++)
    {
      if (modules_regions[iMod].type == moduleTypes[i])
      {
        modCounter++;
      }
    }
    moduleNumberMap.insert(std::make_pair(moduleTypes[i], modCounter));
  }

  // sum all modules
  // std::map<int,int>::iterator it_moduleNumberMap = moduleNumberMap.begin();

  int totModules = 0;
  for (std::map<int, int>::iterator it = moduleNumberMap.begin(); it != moduleNumberMap.end(); ++it)
  {
    moduleNumberOfBeforeMap.insert(std::make_pair(it->first, totModules));
    totModules += it->second;
  }
  feedbackString << "Total number of modules                 = " << totModules << std::endl;

  int NofModuleTypes = moduleTypes.size();
  // int NofSpacalSegments = moduleTypes.size();
  feedbackString << "Number of module types                 = " << NofModuleTypes << std::endl;

  std::map<int, int> mapOfCellMaps; // map from showerModuleType to i-index of cell_map_per_module (then j is 0 or 1)
  // produce 2 cell maps for each module type
  TH2I ***cell_map_per_module;
  cell_map_per_module = new TH2I **[NofModuleTypes];
  for (int i = 0; i < NofModuleTypes; i++)
  {
    cell_map_per_module[i] = new TH2I *[2];
    int det_counter = 0;
    for (int j = 0; j < 2; j++)
    {
      std::stringstream cname;
      cname << "cell_map_";
      if (j == 0) // front
      {
        cname << "front_";
      }
      else // back
      {
        cname << "back_";
      }
      cname << moduleTypes[i];
      // if sections == 1, then keep always front_back
      int fb = j;
      if (moduleSections[i] == 1)
      {
        fb = -1;
      }
      cell_map_per_module[i][j] = ComputeElementMap(cells, cname.str(), fb, moduleTypes[i], det_counter);
      mapOfCellMaps.insert(std::make_pair(moduleTypes[i], i));
    }
  }

  // vector of found
  std::vector<std::vector<bool>> foundLimit;

  std::vector<region_t> cells_regions = CalculateRegions(cells);

  for (int i = 0; i < moduleTypes.size(); i++)
  {
    // feedbackString << "Type " << moduleTypes[i] << " has " << moduleSections[i] << " sections " << std::endl;
    // feedbackString << "Type " << moduleTypes[i] << " has " << moduleNumberMap[moduleTypes[i]] << " modules " << std::endl;
    moduleSectionsMap.insert(std::make_pair(moduleTypes[i], moduleSections[i]));
    moduleSeparationMap.insert(std::make_pair(moduleTypes[i], moduleSeparationZ[i] + modulePositionZ[i]));
    std::vector<bool> foundLimitForModule;
    foundLimitForModule.push_back(false);
    foundLimitForModule.push_back(false);
    foundLimit.push_back(foundLimitForModule);
  }

  for (int i = 0; i < moduleTypes.size(); i++)
  {
    for (int iCell = 0; iCell < cells_regions.size(); iCell++)
    {
      if (foundLimit[i][0] && foundLimit[i][1])
      {
        break;
      }
      // look for a cell of this type
      if (cells_regions[iCell].type == moduleTypes[i])
      {
        if (moduleSectionsMap.at(moduleTypes[i]) == 1)
        {
          // only 1 sections, just take max and min
          foundLimit[i][1] = true;
          moduleMaxMap.insert(std::make_pair(moduleTypes[i], cells_regions[iCell].max[2] + modulePositionZ[i]));
          foundLimit[i][0] = true;
          moduleMinMap.insert(std::make_pair(moduleTypes[i], cells_regions[iCell].min[2] + modulePositionZ[i]));
        }
        else
        {
          // find mean
          float meanCellZ = (cells_regions[iCell].min[2] + cells_regions[iCell].max[2]) / 2.0;
          if (meanCellZ > 0)
          {
            if (foundLimit[i][1] == false)
            {
              foundLimit[i][1] = true;
              moduleMaxMap.insert(std::make_pair(moduleTypes[i], cells_regions[iCell].max[2] + modulePositionZ[i]));
            }
          }
          else
          {
            if (foundLimit[i][0] == false)
            {
              foundLimit[i][0] = true;
              moduleMinMap.insert(std::make_pair(moduleTypes[i], cells_regions[iCell].min[2] + modulePositionZ[i]));
            }
          }
        }
      }
    }
  }

  for (int i = 0; i < moduleTypes.size(); i++)
  {
    feedbackString << "Type " << moduleTypes[i] << " has " << moduleSectionsMap.at(moduleTypes[i]) << " sections " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " has " << moduleNumberMap.at(moduleTypes[i]) << " modules " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " has " << moduleNumberOfBeforeMap.at(moduleTypes[i]) << " modules before" << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " min " << moduleMinMap.at(moduleTypes[i]) << " mm " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " max " << moduleMaxMap.at(moduleTypes[i]) << " mm " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " separation at " << moduleSeparationMap.at(moduleTypes[i]) << " mm " << std::endl;
    if(usingFractionMap)
    {
      feedbackString << "Type " << moduleTypes[i] << " has sampling fracton " << fractionMap.at(moduleTypes[i])  << std::endl;
    }
  }
  
  if(usingFractionMap) feedbackString << std::endl << "==> Using different sampling fractions" << std::endl;
  if(sumAll)           feedbackString << std::endl << "==> Summing all energy depositions, not just in crystals" << std::endl;

  feedbackString << std::endl;
  std::cout << feedbackString.str() << std::endl;
  

  outputFile->cd();
  




  //----------------------------------//
  // LOOP ON INPUT                    //
  //----------------------------------//
  // TH1F *hTotal = new TH1F("hTotal","hTotal",100000,0,100000);
  // TH1F *hFront = new TH1F("hFront","hFront",100000,0,100000);
  // TH1F *hBack = new TH1F("hBack","hBack",100000,0,100000);

  std::vector<float> v_total;
  std::vector<float> v_front;
  std::vector<float> v_back;


  TFile *inputFile = ((TFile *)0);
  TTree *shower = ((TTree*) 0);
  // now run on all input files

  float totalEnergyPerEvent = 0;
  float totalEnergyPerEventFront = 0;
  float totalEnergyPerEventBack = 0;
  for(unsigned int i = 0 ; i < listInputFiles.size(); i++)
  {
    if (inputFile) delete inputFile; // just a precaution
    // if(i == 0)
    // {
    //   // keep the first file...
    //   inputFile = firstFile;
    // }


    inputFile = new TFile(listInputFiles[i].c_str());

    if (!inputFile || inputFile->IsZombie() || inputFile->TestBit(TFile::kRecovered))
    {
      std::cout << "Skipping file " << listInputFiles[i] << std::endl;
      continue;
    }
    std::cout << "Opening file " << listInputFiles[i] << std::endl;
    inputFile->cd();
    shower = (TTree*) inputFile->Get("shower");
    if (!shower) continue; // requested ROOT file does not exist or is unreadable



    // read input shower
    //----------------------------------//
    // INPUT SHOWER                     //
    //----------------------------------//
    int              showerRun              ;
    int              showerEvent            ;
    int              showerIsInCrystal      ;
    int              showerModuleType       ;
    int              shower_pdgID           ;
    int              shower_primaryID       ;
    int              shower_primaryPDGID    ;
    float            shower_primaryEnergy   ;
    float            showerX                ;
    float            showerY                ;
    float            showerZ                ;
    float            showerT                ;
    float            showerTotalEnDep       ;
    float            showerIonizingEnDep    ;
    float            showerNonIonizingEnDep ;
    std::string      *showerProcessName = 0 ;
    std::string      *showerMaterialName = 0;
    int              showerMaterialNumber   ;
    float            shower_local_X         ;
    float            shower_local_Y         ;
    float            shower_local_Z         ;
    float            shower_primaryPositionOnAbsorberX    ;
    float            shower_primaryPositionOnAbsorberY    ;
    float            shower_primaryPositionOnAbsorberZ    ;
    float            shower_primaryMomentumOnAbsorberX    ;
    float            shower_primaryMomentumOnAbsorberY    ;
    float            shower_primaryMomentumOnAbsorberZ    ;
    float            shower_primaryEnergyOnAbsorber       ;
    TBranch *b_showerRun;
    TBranch *b_showerEvent;
    TBranch *b_showerIsInCrystal;
    TBranch *b_showerModuleType       ;
    TBranch *b_shower_pdgID           ;
    TBranch *b_shower_primaryID       ;
    TBranch *b_shower_primaryPDGID    ;
    TBranch *b_shower_primaryEnergy   ;
    TBranch *b_showerX    ;
    TBranch *b_showerY    ;
    TBranch *b_showerZ    ;
    TBranch *b_showerT    ;
    TBranch *b_showerTotalEnDep ;
    TBranch *b_showerIonizingEnDep;
    TBranch *b_showerNonIonizingEnDep;
    TBranch *b_showerMaterialNumber;
    TBranch *b_shower_local_X    ;
    TBranch *b_shower_local_Y    ;
    TBranch *b_shower_local_Z    ;
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
    shower->SetBranchAddress("moduleType"         ,&showerModuleType      ,&b_showerModuleType      );
    shower->SetBranchAddress("pdgID"              ,&shower_pdgID          ,&b_shower_pdgID          );
    shower->SetBranchAddress("primaryID"          ,&shower_primaryID      ,&b_shower_primaryID      );
    shower->SetBranchAddress("primaryPDGID"       ,&shower_primaryPDGID   ,&b_shower_primaryPDGID   );
    shower->SetBranchAddress("primaryEnergy"      ,&shower_primaryEnergy  ,&b_shower_primaryEnergy  );
    shower->SetBranchAddress("x"                  ,&showerX               ,&b_showerX               );
    shower->SetBranchAddress("y"                  ,&showerY               ,&b_showerY               );
    shower->SetBranchAddress("z"                  ,&showerZ               ,&b_showerZ               );
    shower->SetBranchAddress("t"                  ,&showerT               ,&b_showerT               );
    shower->SetBranchAddress("totalEnDep"         ,&showerTotalEnDep      ,&b_showerTotalEnDep      );
    shower->SetBranchAddress("ionizingEnDep"      ,&showerIonizingEnDep   ,&b_showerIonizingEnDep   );
    shower->SetBranchAddress("nonIonizingEnDep"   ,&showerNonIonizingEnDep,&b_showerNonIonizingEnDep);
    shower->SetBranchAddress("processName"        ,&showerProcessName);
    shower->SetBranchAddress("materialName"       ,&showerMaterialName);
    shower->SetBranchAddress("materialNumber"     ,&showerMaterialNumber  ,&b_showerMaterialNumber  );
    shower->SetBranchAddress("localX"             ,&shower_local_X        ,&b_shower_local_X        );
    shower->SetBranchAddress("localY"             ,&shower_local_Y        ,&b_shower_local_Y        );
    shower->SetBranchAddress("localZ"             ,&shower_local_Z        ,&b_shower_local_Z        );
    shower->SetBranchAddress("primary_PositionOnAbsorberX" ,&shower_primaryPositionOnAbsorberX    ,&b_shower_primaryPositionOnAbsorberX);
    shower->SetBranchAddress("primary_PositionOnAbsorberY" ,&shower_primaryPositionOnAbsorberY    ,&b_shower_primaryPositionOnAbsorberY);
    shower->SetBranchAddress("primary_PositionOnAbsorberZ" ,&shower_primaryPositionOnAbsorberZ    ,&b_shower_primaryPositionOnAbsorberZ);
    shower->SetBranchAddress("primary_MomentumOnAbsorberX" ,&shower_primaryMomentumOnAbsorberX    ,&b_shower_primaryMomentumOnAbsorberX);
    shower->SetBranchAddress("primary_MomentumOnAbsorberY" ,&shower_primaryMomentumOnAbsorberY    ,&b_shower_primaryMomentumOnAbsorberY);
    shower->SetBranchAddress("primary_MomentumOnAbsorberZ" ,&shower_primaryMomentumOnAbsorberZ    ,&b_shower_primaryMomentumOnAbsorberZ);
    shower->SetBranchAddress("primary_EnergyOnAbsorber"    ,&shower_primaryEnergyOnAbsorber       ,&b_shower_primaryEnergyOnAbsorber   );
    // //------------------------------------------------//


    int eventID = -1;
    int feedbackDivision;
    long int nEntriesShower = shower->GetEntries();
    std::cout << "Energy deposition events = " << nEntriesShower <<  std::endl;
    std::cout << "Running on energy depositions..." << std::endl;
    long int counter = 0;
    int lostScintillation = 0;
    int lostShower = 0;
    bool sumEvent = false;
    feedbackDivision = (int) (nEntriesShower/20);
    if(feedbackDivision == 0) feedbackDivision = 1;

    
    if(verbose > 0)
    {
      // header
      std::cout << "event\t"
              << "front\t"
              << "back\t"
              << "f/(f+b)\t" 
              << std::endl;
    }
    for(int iEntry = 0 ; iEntry < nEntriesShower ; iEntry++)
    {
      // // progress feedback
      // if (  iEntry % feedbackDivision == 0 )
      // {
      //   std::cout << (int) ((double) iEntry/nEntriesShower * 100) << "% done..." << std::endl;
      // }
      shower->GetEvent(iEntry);

      // if(singleModuleType)
      // {
      //   if (showerModuleType != moduleTypeNumber)
      //   {
      //     continue;
      //   }
      // }
      // summation of energy needs to be performed only until a new event arrives in the list
      // an event corresponds to 1 primary shot. the info is stored in event
      // so we start be setting a event indicator, eventID, to -1. At the first loop step,
      // the first value of variable event will be 0, so what comes in the next if statement
      // (event saving and clean up) will not be executed, and instead the event accumulation starts
      // and the end of event accumulation the eventID becomes = event, then a new step of the cycle
      // starts. Now eventID is = 0, so "event saving and clean up" will not be performed until
      // the value of event becomes != 0. When that happens, the variables were the enegies were
      // summed are used to filled the relevant histograms, than they are cleared and a new accumulation
      // starts. Of course this means that after the end of the loop, the last event still needs to be saved
      //
      //
      //------------------------------//
      // event saving and clean up    //
      //------------------------------//
      // std::cout << eventID << " " << showerEvent;
      if((eventID != showerEvent) && (eventID != -1)) // new event, dump the previous
      {
        // fill
        // hTotal->Fill(totalEnergyPerEvent);
        // hFront->Fill(totalEnergyPerEventFront);
        // hBack->Fill(totalEnergyPerEventBack);
        v_total.push_back(totalEnergyPerEvent);
        v_front.push_back(totalEnergyPerEventFront);
        v_back .push_back(totalEnergyPerEventBack);
        
        if(verbose > 0)
        { 
          // command line output 
          std::cout << eventID << "\t" 
                    << totalEnergyPerEventFront << "\t" 
                    << totalEnergyPerEventBack << "\t";
          if( ( ( (float) totalEnergyPerEventFront) + ( (float) totalEnergyPerEventBack) ) != 0 )
          {
            std::cout << ((float) totalEnergyPerEventFront) / ( ( (float) totalEnergyPerEventFront) + ( (float) totalEnergyPerEventBack) );
          } 
          else 
          {
            std::cout << "0";
          }
          std::cout << std::endl;
        }        


        // cleanup
        totalEnergyPerEvent = 0;
        totalEnergyPerEventFront = 0;
        totalEnergyPerEventBack = 0;


      }

      // check if events has to be summed or not
      
      sumEvent = false;
      if(sumAll)
      {
        sumEvent = true;
      }
      else 
      {
        if(showerIsInCrystal)
        {
          sumEvent = true;
        }
      }

      if (sumEvent) 
      {
        // necessary distinction for depositions in absorber. 
        // we do not save the module type in these case (basically to save CPU during simulation runs)
        // so we need to ensure we don't check for module type
        // and set artificially the separation to 0, i.e. the front back distribution
        // when you use --sum-all will have no meaning, unless by coincidence the 
        // module was in fact separated in z = 0 
        float moduleSeparation = 0.;

        if(!sumAll)
        {
          //fetch separation point in z
          std::map<int, float> ::const_iterator iterModSep = moduleSeparationMap.find(showerModuleType);
          if(iterModSep != moduleSeparationMap.end())
          {
            // iter is item pair in the map. The value will be accessible as `iter->second`.
            moduleSeparation = (iterModSep->second);
          }
          else
          {
            if(verbose > 0)
            {
              std::cout << "WARNING: Module separation not found for showerModuleType " << showerModuleType << std::endl;
              std::cout << "         This event will be ignored! Check the input files provided." << std::endl;
              std::cout << std::endl;
            }
            //lostScintillation++;
            continue;
          }
        }

        if(usingFractionMap)
        {
          // get sampling fraction
          float samplingFraction;
          std::map<int, float> ::const_iterator iterModFrac = fractionMap.find(showerModuleType);
          if(iterModFrac != fractionMap.end())
          {
            // iter is item pair in the map. The value will be accessible as `iter->second`.
            samplingFraction = (iterModFrac->second);
          }
          else
          {
            if(verbose > 0)
            {
              std::cout << "WARNING: Sampling fraction not found for showerModuleType " << showerModuleType << std::endl;
              std::cout << "         This event will be ignored! Check the input files provided." << std::endl;
              std::cout << std::endl;
            }
            //lostScintillation++;
            continue;
          }
          showerTotalEnDep = showerTotalEnDep / samplingFraction;
        }

        // correct by sampling, if specified
        // accumulate total energy
        totalEnergyPerEvent += showerTotalEnDep;
        // std::cout << showerZ << " " << moduleSeparation << " " << showerTotalEnDep ;
        if(showerZ < moduleSeparation) // then it's front
        {
          totalEnergyPerEventFront += showerTotalEnDep;
          // std::cout << " 0" << std::endl;
        }
        else
        {
          totalEnergyPerEventBack += showerTotalEnDep;
          // std::cout << " 1" << std::endl;
        }
      }
      eventID = showerEvent;
    }
    // fill last event!
    // fill
    // hTotal->Fill(totalEnergyPerEvent);
    // hFront->Fill(totalEnergyPerEventFront);
    // hBack->Fill(totalEnergyPerEventBack);
    v_total.push_back(totalEnergyPerEvent);
    v_front.push_back(totalEnergyPerEventFront);
    v_back .push_back(totalEnergyPerEventBack);
    
    if(verbose > 0)
    {
      // command line output 
      std::cout << eventID << "\t" 
              << totalEnergyPerEventFront << "\t" 
              << totalEnergyPerEventBack << "\t";
      if( ( ( (float) totalEnergyPerEventFront) + ( (float) totalEnergyPerEventBack) ) != 0 )
      {
        std::cout << ((float) totalEnergyPerEventFront) / ( ( (float) totalEnergyPerEventFront) + ( (float) totalEnergyPerEventBack) );
      } 
      else 
      {
        std::cout << "0";
      }
      std::cout << std::endl;
    }
    // cleanup
    totalEnergyPerEvent = 0;
    totalEnergyPerEventFront = 0;
    totalEnergyPerEventBack = 0;
  }

  
  // fill histograms
  // first calc mean and stdevs
  float avg_total = 0. ;
  float avg_front = 0. ;
  float avg_back  = 0. ;
  float std_total = 0. ;
  float std_front = 0. ;
  float std_back  = 0. ;
  unsigned int nt = v_total.size(); 
  unsigned int nf = v_front.size(); 
  unsigned int nb = v_back.size();
  if ( nt != 0) {
     avg_total = accumulate( v_total.begin(), v_total.end(), 0.0) / nt; 
     double sq_sum = std::inner_product(v_total.begin(), v_total.end(), v_total.begin(), 0.0);
     std_total = std::sqrt(sq_sum / v_total.size() - avg_total * avg_total);
  }
  if ( nf != 0) {
     avg_front = accumulate( v_front.begin(), v_front.end(), 0.0) / nf; 
     double sq_sum = std::inner_product(v_front.begin(), v_front.end(), v_front.begin(), 0.0);
     std_front = std::sqrt(sq_sum / v_front.size() - avg_front * avg_front);
  }
  if ( nb != 0) {
     avg_back  = accumulate( v_back.begin(), v_back.end(), 0.0) / nb; 
     double sq_sum = std::inner_product(v_back.begin(), v_back.end(), v_back.begin(), 0.0);
     std_back = std::sqrt(sq_sum / v_back.size() - avg_back * avg_back);
  }

  // keep assuming the histogram to go until 100 GeV, but calculate appropriate bin size 
  // by assuming 100 bins in the +/- 5 stdevs range
  float max = 100000.;
  int bins = 100;
  float range_total = (avg_total + 4.0*std_total) - (avg_total - 4.0*std_total);
  float range_front = (avg_front + 4.0*std_front) - (avg_front - 4.0*std_front);
  float range_back  = (avg_back + 4.0*std_back) - (avg_back - 4.0*std_back);
  float bin_size_total = range_total / bins;
  float bin_size_front = range_front / bins;
  float bin_size_back  = range_back  / bins;
  int bins_total = (int) std::round(max / bin_size_total); 
  int bins_front = (int) std::round(max / bin_size_front); 
  int bins_back  = (int) std::round(max / bin_size_back ); 

  TH1F *hTotal = new TH1F("hTotal","hTotal",bins_total,0,max);
  TH1F *hFront = new TH1F("hFront","hFront",bins_front,0,max);
  TH1F *hBack  = new TH1F("hBack","hBack"  ,bins_back ,0,max);
  
  for(unsigned int i = 0; i < v_total.size(); i++)
  {
    hTotal->Fill(v_total[i]);
  }
  for(unsigned int i = 0; i < v_front.size(); i++)
  {
    hFront->Fill(v_front[i]);
  }
  for(unsigned int i = 0; i < v_back.size(); i++)
  {
    hBack->Fill(v_back[i]);
  }

  std::cout << "Writing output to file " << outputFileName << std::endl;
  outputFile->cd();
  hTotal->Write();
  hFront->Write();
  hBack->Write();
  module_ID_map->Write();

  // detector_ID_map_front->Write();
  // detector_ID_map_back->Write();
  for(int i = 0; i < NofModuleTypes; i++)
  {
    for(int j = 0; j < 2 ; j++)
    {
      cell_map_per_module[i][j]->Write();
    }
  }

  std::cout << "Done." << std::endl;
  std::cout << "Energy resolution:" << std::endl;
  std::cout << hTotal->GetRMS() / hTotal->GetMean() << " " << ( TMath::Sqrt(TMath::Power( (hTotal)->GetStdDevError() / (hTotal)->GetStdDev() ,2) + TMath::Power( (hTotal)->GetMeanError() / (hTotal)->GetMean(),2)  ) * (hTotal)->GetStdDev()/(hTotal)->GetMean()) << std::endl;
  outputFile->Close();
  return 0;
}