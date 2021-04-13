// Marco Pizzichemi 29.04.2020 marco.pizzichemi@cern.ch

// compile with
// g++ -o ../build/showerAnalysis showerAnalysis.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore

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
#include "ROOT/RDataFrame.hxx"


#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <getopt.h>
#include <limits>
#include <dirent.h>
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
            << "\t\t" << "[--type]         <module type - default = 0 >  " << std::endl
            << "\t\t" << "[--max]          <max distance from shower line to test [mm] - default = 30. >  " << std::endl
            << "\t\t" << "[--divs]         <divisions of max distance - default = 10 >  " << std::endl
            << "\t\t" << "[--stop]         <stop after N event, default = INFINITY >  " << std::endl
            << "\t\t" << "[--verbose]      <verbosity level - default = 0>  " << std::endl
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
  // ROOT::EnableImplicitMT(); // Enable ROOT's implicit multi-threading

  

  //------------------------------------------------//
  // PARSE INPUT                                    //
  //------------------------------------------------//
  // prepare variable, giving defaults
  std::string inputFilePrefix = "";
  std::string inputFolderName = "./";
  std::string outputFileName = "";
  // bool singleEvent = false;
  int verbose = 0;
  int moduleTypeNumber = 0;
  bool singleModuleType = false;
  float maxRadius_mm = 30.;
  int radius_divs = 10;
  int maxEvent = std::numeric_limits<int>::max();
  float th_moliere_mm = 10.0;
  static struct option longOptions[] =
  {
			{ "input", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { "folder", required_argument, 0, 0 },
      { "verbose", required_argument, 0, 0 },
      { "type", required_argument, 0, 0 },
      { "max", required_argument, 0, 0 },
      { "divs", required_argument, 0, 0 },
      { "stop", required_argument, 0, 0 },
      { "moliere", required_argument, 0, 0 },
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
      singleModuleType = true;
      moduleTypeNumber = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 5){
      maxRadius_mm = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 6){
      radius_divs = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 7){
      maxEvent = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 8){
      th_moliere_mm = atof((char *)optarg);
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
  }

  feedbackString << std::endl;
  std::cout << feedbackString.str() << std::endl;
  

  outputFile->cd();
  




  //----------------------------------//
  // LOOP ON INPUT                    //
  //----------------------------------//
  TH1F *hTotal = new TH1F("hTotal","hTotal",100000,0,100000);
  TH1F *hFront = new TH1F("hFront","hFront",100000,0,100000);
  TH1F *hBack = new TH1F("hBack","hBack",100000,0,100000);
  TH1F *hMoliere = new TH1F("hMoliere","Radius of cylinder containing 90 percent of energy",500,0,50);
  TH1F *hMoliere2 = new TH1F("hMoliere2","Fraction of Energy contained in 1 theoretical Moliere Radius",100,0,1);

  TFile *inputFile = ((TFile *)0);
  TTree *shower = ((TTree*) 0);
  // now run on all input files

  float totalEnergyPerEvent = 0;
  float totalEnergyPerEventFront = 0;
  float totalEnergyPerEventBack = 0;

  // 
  float d ;
  float t ; 
  float Sx;
  float Sy;
  float Sz;
  float distance ;

  float moduleSeparation;
  
  std::vector<float> gx,gy,gMol;
  std::vector<float> AVGgx,AVGgy,AVGgMol;
  std::vector<float> energy_depo_within;
  int tot_events = 0;
  // std::vector<float> rm_list;
 
  // create the radius steps 
  // and the vector of energy accumulators
  for(int i = 0; i < radius_divs; i++)
  {
    gy.push_back((maxRadius_mm/radius_divs)*(i+1));
    gx.push_back(0.);
    gMol.push_back(0.);
    AVGgy.push_back(0.);
    AVGgx.push_back(0.);
    AVGgMol.push_back(0.);
    energy_depo_within.push_back(0.);
  }

  // for(int i = 0; i < radius_divs; i++)
  // {
  //   std::cout << gy[i] << " " << gx[i] << " " << energy_depo_within[i] << std::endl;
  // }

  std::vector<TGraph*> graphs;
  std::vector<TGraph*> graphs_moliere;
  std::vector<TGraph*> graphs_mm;

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

  
  // ROOT::RDataFrame shower("shower",inputFileName);
  // auto isInCrystal = [](Float_t &PositionZ) { return PositionZ > 0; };

  // primaries, hardcoded values, only need to run once...
  // ROOT::RDataFrame primaries("primaries",listInputFiles[i].c_str());
  // auto v_meanXprimaries = primaries.Mean("PositionOnAbsorberX");
  // auto v_meanYprimaries = primaries.Mean("PositionOnAbsorberY");
  // auto v_meanZprimaries = primaries.Mean("PositionOnAbsorberZ");
  // auto v_meanMXprimaries = primaries.Mean("MomentumOnAbsorberX");
  // auto v_meanMYprimaries = primaries.Mean("MomentumOnAbsorberY");
  // auto v_meanMZprimaries = primaries.Mean("MomentumOnAbsorberZ");
   
  // float Ox = *v_meanXprimaries ;
  // float Oy = *v_meanYprimaries ;
  // float Oz = *v_meanZprimaries ;
  // float Vx = *v_meanMXprimaries;
  // float Vy = *v_meanMYprimaries;
  // float Vz = *v_meanMZprimaries;



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

      if(showerEvent == maxEvent )
        break;

      if(singleModuleType)
      {
        if (showerModuleType != moduleTypeNumber)
        {
          continue;
        }
      }


      
      
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
      if((eventID != showerEvent) && (eventID != -1)) // new event, dump the previous, and store the primaries info
      {
        // fill
        hTotal->Fill(totalEnergyPerEvent);
        hFront->Fill(totalEnergyPerEventFront);
        hBack->Fill(totalEnergyPerEventBack);
        
        if(verbose > 0)
        { 
          // command line output 
          std::cout << eventID << "\t" 
                    << totalEnergyPerEventFront << "\t" 
                    << totalEnergyPerEventBack << "\t"
                    << ((float) totalEnergyPerEventFront) / ( ( (float) totalEnergyPerEventFront) + ( (float) totalEnergyPerEventBack) )
                    << std::endl;
        }        

        // calculate moliere quantities 
        // and update averages
        tot_events++;
        for(int iRad = 0; iRad < gy.size(); iRad++)
        {
          
          gx[iRad] = ((float) energy_depo_within[iRad])/((float) totalEnergyPerEvent);
          gMol[iRad] = gy[iRad] / th_moliere_mm;
          AVGgy[iRad] += gy[iRad];
          AVGgx[iRad] += gx[iRad];
          AVGgMol[iRad] += gMol[iRad];
          // std::cout << gy[iRad] << " " << energy_depo_within[iRad] << " " <<  gx[iRad] << " ";
        }
        // std::cout << std::endl;
        // make the TGraph 
        TGraph* moliere = new TGraph(gx.size(),&gx[0],&gy[0]);
        std::stringstream sname;
        sname << "graph" << eventID;
        moliere->SetName(sname.str().c_str());
        moliere->SetTitle(sname.str().c_str());
        moliere->GetXaxis()->SetTitle("Fraction of energy contained in cylinder");
        moliere->GetYaxis()->SetTitle("Cylinder radius [mm]");
        TCanvas *cTemp = new TCanvas("cTemp","cTemp",1200,800);
        cTemp->cd();
        moliere->Draw("AL*");
        std::cout << eventID << " " << moliere->Eval(0.9) ;
        hMoliere->Fill(moliere->Eval(0.9));
        // rm_list.push_back(moliere->Eval(0.9));
        graphs.push_back(moliere);
        sname.str("");
        delete cTemp;

        
        TGraph* moliere2 = new TGraph(gMol.size(),&gMol[0],&gx[0]);
        sname << "graph_moliere" << eventID;
        moliere2->SetName(sname.str().c_str());
        moliere2->SetTitle(sname.str().c_str());
        moliere2->GetYaxis()->SetTitle("Fraction of energy contained in cylinder");
        moliere2->GetXaxis()->SetTitle("Units of Theoretical Moliere Radius");
        cTemp = new TCanvas("cTemp","cTemp",1200,800);
        cTemp->cd();
        moliere2->Draw("AL*");
        std::cout << " " << moliere2->Eval(1) << std::endl;
        hMoliere2->Fill(moliere2->Eval(1));
        // rm_list.push_back(moliere->Eval(0.9));
        graphs_moliere.push_back(moliere2);
        sname.str("");
        delete cTemp;

        TGraph* moliere_mm = new TGraph(gMol.size(),&gMol[0],&gx[0]);
        sname << "graph_moliere_mm" << eventID;
        moliere_mm->SetName(sname.str().c_str());
        moliere_mm->SetTitle(sname.str().c_str());
        moliere_mm->GetYaxis()->SetTitle("Fraction of energy contained in cylinder");
        moliere_mm->GetXaxis()->SetTitle("Cylinder Radius [mm]");
        cTemp = new TCanvas("cTemp","cTemp",1200,800);
        cTemp->cd();
        moliere_mm->Draw("AL*");
        // std::cout << " " << moliere_mm->Eval(1) << std::endl;
        // moliere_mm->Fill(moliere_mm->Eval(1));
        // rm_list.push_back(moliere->Eval(0.9));
        graphs_mm.push_back(moliere_mm);
        sname.str("");
        delete cTemp;





        // cleanup
        totalEnergyPerEvent = 0;
        totalEnergyPerEventFront = 0;
        totalEnergyPerEventBack = 0;
        for(int iRad = 0; iRad < gy.size(); iRad++)
        {
          energy_depo_within[iRad] = 0;
        }
      }

      // if (showerIsInCrystal && showerEvent == 0)
      if (showerIsInCrystal ) 
      {
        
        //fetch separation point in z
        moduleSeparation;
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

        // Ox = shower_primaryMomentumOnAbsorberX;
        // Oy = shower_primaryPositionOnAbsorberY;
        // Oz = shower_primaryPositionOnAbsorberZ;
        // Vx = shower_primaryMomentumOnAbsorberX;
        // Vy = shower_primaryMomentumOnAbsorberY;
        // Vz = shower_primaryMomentumOnAbsorberZ;

        // calculate distance to primary line 
        // find plane perpendicular to line that passes by point (x,y,z)
        // general equation of plane
        // ax + by + cz + d = 0 (1)
        // so generic set of planes perpendicular to vector (Vx,Vy,Vz) is just 
        // Vx*x + Vy*y + Vz*z + d = 0 (2)
        // then impose passage for point P and find d
        // Vx*Px + Vy*Py + Vz*Pz + d = 0 (3)
        d = -(shower_primaryMomentumOnAbsorberX*showerX + shower_primaryMomentumOnAbsorberY*showerY + shower_primaryMomentumOnAbsorberZ*showerZ);
        // find intersection of primary line from O in V direction 
        // with this plane. parametric line eq are 
        // x = Ox + Vx*t
        // y = Oy + Vy*t
        // z = Oz + Vz*t
        // substitute in eq (2)
        // Vx*(Ox + Vx*t) + Vy*(Oy + Vy*t) + Vz*(Oz + Vz*t) + d = 0 (4)
        // we get 
        // t = -(d + Vx*Ox +  Vy*Oy + Vz*Oz)/(Vx + Vy + Vz) (5)
        t = -(d + shower_primaryMomentumOnAbsorberX*shower_primaryPositionOnAbsorberX +  shower_primaryMomentumOnAbsorberY*shower_primaryPositionOnAbsorberY + shower_primaryMomentumOnAbsorberZ*shower_primaryPositionOnAbsorberZ)/(shower_primaryMomentumOnAbsorberX + shower_primaryMomentumOnAbsorberY + shower_primaryMomentumOnAbsorberZ);
        // and back in the line parametric eq to find Sx,Sy,Sz
        Sx = shower_primaryPositionOnAbsorberX + shower_primaryMomentumOnAbsorberX*t;
        Sy = shower_primaryPositionOnAbsorberY + shower_primaryMomentumOnAbsorberY*t;
        Sz = shower_primaryPositionOnAbsorberZ + shower_primaryMomentumOnAbsorberZ*t;
        // finally, calculate the distance between P and S
        distance = fabs(sqrt(pow((Sx-showerX),2) + pow((Sy-showerY),2) + pow((Sz-showerZ),2)));

        // std::cout << d << " "
        //           << t << " "
        //           // << Ox << " "
        //           // << Oy << " "
        //           // << Oz << " "
        //           // << Vx << " "
        //           // << Vy << " "
        //           // << Vz << " "
        //           << Sx << " "
        //           << Sy << " "
        //           << Sz << " "
        //           << shower_primaryPositionOnAbsorberX << " "
        //           << shower_primaryPositionOnAbsorberY << " "
        //           << shower_primaryPositionOnAbsorberZ << " "
        //           << shower_primaryMomentumOnAbsorberX << " "
        //           << shower_primaryMomentumOnAbsorberY << " "
        //           << shower_primaryMomentumOnAbsorberZ << " "
        //           << distance << std::endl;

        for(int iRad = 0; iRad < gy.size(); iRad++)
        {
          if(distance < gy[iRad])
          {
            energy_depo_within[iRad] += showerTotalEnDep;
          }
        }
      }
      eventID = showerEvent;
    }
    // fill last event!
    // fill
    hTotal->Fill(totalEnergyPerEvent);
    hFront->Fill(totalEnergyPerEventFront);
    hBack->Fill(totalEnergyPerEventBack);

    // calculate moliere quantities for last event
    // only if event limit was not set
    if(maxEvent == std::numeric_limits<int>::max())
    {
      tot_events++;
      for(int iRad = 0; iRad < gy.size(); iRad++)
      {
        gx[iRad] = ((float) energy_depo_within[iRad])/((float) totalEnergyPerEvent);
        gMol[iRad] = gy[iRad] / th_moliere_mm;
        AVGgy[iRad] += gy[iRad];
        AVGgx[iRad] += gx[iRad];
        AVGgMol[iRad] += gMol[iRad];
        // std::cout << gy[iRad] << " " << energy_depo_within[iRad] << " " <<  gx[iRad] << " ";
      }
      // std::cout << std::endl;
      // make the TGraph 
      TGraph* moliere = new TGraph(gx.size(),&gx[0],&gy[0]);
      std::stringstream sname;
      sname << "graph" << eventID;
      moliere->SetName(sname.str().c_str());
      moliere->SetTitle(sname.str().c_str());
      TCanvas *cTemp = new TCanvas("cTemp","cTemp",1200,800);
      cTemp->cd();
      moliere->Draw("AL*");
      std::cout << showerEvent << " " << moliere->Eval(0.9) << std::endl;
      hMoliere->Fill(moliere->Eval(0.9));
      graphs.push_back(moliere);
      sname.str("");
      delete cTemp;

      TGraph* moliere2 = new TGraph(gMol.size(),&gMol[0],&gx[0]);
      sname << "graph_moliere" << eventID;
      moliere2->SetName(sname.str().c_str());
      moliere2->SetTitle(sname.str().c_str());
      moliere2->GetYaxis()->SetTitle("Fraction of energy contained in cylinder");
      moliere2->GetXaxis()->SetTitle("Units of Theoretical Moliere Radius");
      cTemp = new TCanvas("cTemp","cTemp",1200,800);
      cTemp->cd();
      moliere2->Draw("AL*");
      std::cout << " " << moliere2->Eval(1) << std::endl;
      hMoliere2->Fill(moliere2->Eval(1));
      // rm_list.push_back(moliere->Eval(0.9));
      graphs_moliere.push_back(moliere2);
      sname.str("");
      delete cTemp;

      TGraph* moliere_mm = new TGraph(gMol.size(),&gMol[0],&gx[0]);
      sname << "graph_moliere_mm" << eventID;
      moliere_mm->SetName(sname.str().c_str());
      moliere_mm->SetTitle(sname.str().c_str());
      moliere_mm->GetYaxis()->SetTitle("Fraction of energy contained in cylinder");
      moliere_mm->GetXaxis()->SetTitle("Cylinder Radius [mm]");
      cTemp = new TCanvas("cTemp","cTemp",1200,800);
      cTemp->cd();
      moliere_mm->Draw("AL*");
      // std::cout << " " << moliere_mm->Eval(1) << std::endl;
      // moliere_mm->Fill(moliere_mm->Eval(1));
      // rm_list.push_back(moliere->Eval(0.9));
      graphs_mm.push_back(moliere_mm);
      sname.str("");
      delete cTemp;

    } 
    
    
    if(verbose > 0)
    {
      // command line output 
      std::cout << eventID << "\t" 
              << totalEnergyPerEventFront << "\t" 
              << totalEnergyPerEventBack << "\t"
              << ((float) totalEnergyPerEventFront) / ( ( (float) totalEnergyPerEventFront) + ( (float) totalEnergyPerEventBack) )
              << std::endl;
    }
    // cleanup
    totalEnergyPerEvent = 0;
    totalEnergyPerEventFront = 0;
    totalEnergyPerEventBack = 0;
    for(int iRad = 0; iRad < gy.size(); iRad++)
    {
      energy_depo_within[iRad] = 0;
    }
  }

  
  // do the averages and the graphs 
  for(int iRad = 0; iRad < AVGgx.size(); iRad++)
  {
    AVGgx[iRad] = AVGgx[iRad] / tot_events;
  }
  TGraph* AVGmoliere = new TGraph(AVGgx.size(),&AVGgx[0],&gy[0]);
  AVGmoliere->SetName("average_graph");
  AVGmoliere->SetTitle("average_graph");
  AVGmoliere->GetXaxis()->SetTitle("Average Fraction of energy contained in cylinder");
  AVGmoliere->GetYaxis()->SetTitle("Cylinder radius [mm]");
  
      
  TGraph* AVGmoliere2 = new TGraph(gMol.size(),&gMol[0],&AVGgx[0]);
  AVGmoliere2->SetName("average_graph_moliere");
  AVGmoliere2->SetTitle("average_graph_moliere");
  AVGmoliere2->GetYaxis()->SetTitle("Average Fraction of energy contained in cylinder");
  AVGmoliere2->GetXaxis()->SetTitle("Units of Theoretical Moliere Radius");
  
  TGraph* AVGmoliere_mm = new TGraph(gy.size(),&gy[0],&AVGgx[0]);
  AVGmoliere_mm->SetName("average_graph_moliere_per_mm");
  AVGmoliere_mm->SetTitle("average_graph_moliere_per_mm");
  AVGmoliere_mm->GetYaxis()->SetTitle("Average Fraction of energy contained in cylinder");
  AVGmoliere_mm->GetXaxis()->SetTitle("Cylinder Radius [mm]");
  

  std::cout << "Writing output to file " << outputFileName << std::endl;
  outputFile->cd();
  hTotal->Write();
  hFront->Write();
  hBack->Write();
  hMoliere->Write();
  hMoliere2->Write();
  module_ID_map->Write();
  AVGmoliere->Write();
  AVGmoliere2->Write();
  AVGmoliere_mm->Write();
  for(int i = 0; i < graphs.size(); i++)
  {
    graphs[i]->Write();
  }
  for(int i = 0; i < graphs_moliere.size(); i++)
  {
    graphs_moliere[i]->Write();
  }
  for(int i = 0; i < graphs_mm.size(); i++)
  {
    graphs_mm[i]->Write();
  }

  // detector_ID_map_front->Write();
  // detector_ID_map_back->Write();
  for(int i = 0; i < NofModuleTypes; i++)
  {
    for(int j = 0; j < 2 ; j++)
    {
      cell_map_per_module[i][j]->Write();
    }
  }

  // BROKEN FOR NOW
  // for(int iMod = 0 ; iMod < modules_regions.size(); iMod++)
  // {
  //   std::cout << "iMod " << iMod << std::endl;
  //   std::cout << "iMod.save" << module[iMod].save << std::endl;
  //   if(module[iMod].save)
  //   {
  //     module[iMod].h_EnDepInCrystal   ->Write();
  //     std::cout << "1" << std::endl;
  //     module[iMod].h_EnDepNotInCrystal->Write();
  //     std::cout << "2" << std::endl;
  //     module[iMod].h_EnDepTotal       ->Write();
  //     std::cout << "3" << std::endl;
  //     for(int iSec = 0; iSec < module[iMod].sections; iSec++)
  //     {
  //       std::cout << "iSec " << iSec << std::endl;
  //       for(int iCel = 0; iCel <  module[iMod].nCellsPerSections[iSec]; iCel++)
  //       {
  //         std::cout << "iCel " << iCel << std::endl;
  //         if(module[iMod].saveCell[iSec][iCel])
  //         {
  //           module[iMod].h_Cell_EnDepInCrystal[iSec][iCel]->Write();
  //           std::cout << "4" << std::endl;
  //         }
  //       }
  //     }
  //   }
  //
  // }
  //
  //
  // // manual delete
  // for(int iMod = 0 ; iMod < modules_regions.size(); iMod++)
  // {
  //
  //
  //     for(int iSec = 0; iSec < module[iMod].sections; iSec++)
  //     {
  //       for(int iCel = 0; iCel <  module[iMod].nCellsPerSections[iSec]; iCel++)
  //       {
  //         if(module[iMod].saveCell[iSec][iCel])
  //         {
  //           delete module[iMod].h_Cell_EnDepInCrystal[iSec][iCel];
  //         }
  //       }
  //     }
  //
  // }
  // BROKEN FOR NOW
  std::cout << "Done." << std::endl;
  // gROOT->ProcessLine(".qqqqqqqqqqqqqqqqqqqq");
  outputFile->Close();
  return 0;
}

  //
  //
  //
  // //------------------------------------------------//
  // // CALCULATE EN RES.                              //
  // //------------------------------------------------//
  // // prepare vars to sum energies
  // // one variable for each abs, cell and fiber , so dinamic arrays
  // float *absorberEnergyInCrystal;
  // float *absorberEnergyNotInCrystal;
  // float *cellEnergyInCrystal;
  // float *cellEnergyNotInCristal;
  // float *fibreEnergy;
  // // prepare debug variables
  // float *debug_absorberEnergyInCrystal;
  // float *debug_absorberEnergyNotInCrystal;
  // float *debug_cellEnergyInCrystal;
  // float *debug_cellEnergyNotInCristal;
  // float *debug_fibreEnergy;
  // // prepare histograms to plot energies
  // // one histo for each abs, cell and fiber, std::vectors of histos
  // std::stringstream sname;
  // std::stringstream stitle;
  // std::vector<TH1F*> abs_histoIn;
  // std::vector<TH1F*> abs_histoTot;
  // std::vector<TH1F*> cell_histoIn;
  // std::vector<TH1F*> cell_histoTot;
  // std::vector<TH1F*> fiber_histoIn;
  //
  //
  // std::vector<TH1F*> histo_pair_cell_Tot;
  // std::vector<TH1F*> histo_pair_cell_In;
  // std::vector<TH1F*> histo_ratio_pair_cell_Tot;
  // std::vector<TH1F*> histo_ratio_pair_cell_In;
  // float *pair_cellEnergyInCrystal;
  // float *pair_cellEnergyNotInCrystal;
  // float *ratio_pair_cellEnergyInCrystal;
  // float *ratio_pair_cellEnergyNotInCrystal;
  //
  // //hardcoded
  //
  // float energyInFrontInCry = 0.;
  // float energyInBackInCry  = 0.;
  // float energyInFrontNotInCry = 0.;
  // float energyInBackNotInCry  = 0.;
  // TH1F *h_energyInFrontInCry = new TH1F("EnergyInFrontInCry","EnergyInFrontInCry",bins,min,max);
  // h_energyInFrontInCry->GetXaxis()->SetTitle("MeV");
  // TH1F *h_energyInBackInCry  = new TH1F("EnergyInBackInCry" ,"EnergyInBackInCry" ,bins,min,max);
  // h_energyInBackInCry->GetXaxis()->SetTitle("MeV");
  // TH1F *h_energyInFrontInAbs = new TH1F("EnergyInFrontInAbs","EnergyInFrontInAbs",bins,min,max);
  // h_energyInFrontInAbs->GetXaxis()->SetTitle("MeV");
  // TH1F *h_energyInBackInAbs  = new TH1F("EnergyInBackInAbs" ,"EnergyInBackInAbs" ,bins,min,max);
  // h_energyInBackInAbs->GetXaxis()->SetTitle("MeV");
  // pair_cellEnergyInCrystal = new float[cells_regions.size()/2];
  // pair_cellEnergyNotInCrystal = new float[cells_regions.size()/2];
  // ratio_pair_cellEnergyInCrystal = new float[cells_regions.size()/2];
  // ratio_pair_cellEnergyNotInCrystal = new float[cells_regions.size()/2];
  // for(int i = 0 ; i < cells_regions.size()/2; i++)
  // {
  //   pair_cellEnergyInCrystal[i]          = 0.0;
  //   pair_cellEnergyNotInCrystal[i]       = 0.0;
  //   ratio_pair_cellEnergyInCrystal[i]    = 0.0;
  //   ratio_pair_cellEnergyNotInCrystal[i] = 0.0;
  // }
  // for(int i = 0 ; i < cells_regions.size()/2; i++)
  // {
  //   sname << "EnInCrystals_Cells_" << i << "_" << i+cells_regions.size()/2;
  //   stitle << "Energy deposited in crystals of Cells " << i << " and " << i+cells_regions.size()/2;
  //   TH1F *t_histoIn = new TH1F(sname.str().c_str(),stitle.str().c_str(),bins,min,max);
  //   t_histoIn->GetXaxis()->SetTitle("MeV");
  //   histo_pair_cell_In.push_back(t_histoIn);
  //   sname.str("");
  //   stitle.str("");
  // }
  // for(int i = 0 ; i < cells_regions.size()/2; i++)
  // {
  //   sname << "TotEn_Cells_" << i << "_" << i+cells_regions.size()/2;
  //   stitle << "Total energy deposited in Cells " << i << " and " << i+cells_regions.size()/2 << " (crystals + absorber + air)";
  //   TH1F *t_histoIn = new TH1F(sname.str().c_str(),stitle.str().c_str(),bins,min,max);
  //   t_histoIn->GetXaxis()->SetTitle("MeV");
  //   histo_pair_cell_Tot.push_back(t_histoIn);
  //   sname.str("");
  //   stitle.str("");
  // }
  // for(int i = 0 ; i < cells_regions.size()/2; i++)
  // {
  //   sname << "RatioBackFront_inCrystals_Cells" << i << "_" << i+cells_regions.size()/2;;
  //   stitle << "Ratio Back/Front, energy deposited in crystals of Cells " << i << " and " << i+cells_regions.size()/2;
  //   TH1F *t_histoIn = new TH1F(sname.str().c_str(),stitle.str().c_str(),bins,0,100);
  //   t_histoIn->GetXaxis()->SetTitle("Ratio back/front");
  //   histo_ratio_pair_cell_In.push_back(t_histoIn);
  //   sname.str("");
  //   stitle.str("");
  // }
  // for(int i = 0 ; i < cells_regions.size()/2; i++)
  // {
  //   sname << "RatioBackFront_Tot_Cells" << i << "_" << i+cells_regions.size()/2;;
  //   stitle << "Ratio Back/Front, total energy deposited in Cells " << i << " and " << i+cells_regions.size()/2 << " (crystals + absorber + air)";
  //   TH1F *t_histoIn = new TH1F(sname.str().c_str(),stitle.str().c_str(),bins,0,100);
  //   t_histoIn->GetXaxis()->SetTitle("Ratio back/front");
  //   histo_ratio_pair_cell_Tot.push_back(t_histoIn);
  //   sname.str("");
  //   stitle.str("");
  // }
  //
  //
  // if(sum_abs)
  // {
  //   absorberEnergyInCrystal = new float[absobers_regions.size()];
  //   absorberEnergyNotInCrystal = new float[absobers_regions.size()];
  //   if(debug) // some debug
  //   {
  //     debug_absorberEnergyInCrystal = new float[absobers_regions.size()];
  //     debug_absorberEnergyNotInCrystal = new float[absobers_regions.size()];
  //   }
  //   for(int i = 0 ; i < absobers_regions.size(); i++)
  //   {
  //     absorberEnergyInCrystal[i] = 0;
  //     absorberEnergyNotInCrystal[i] = 0;
  //     if(debug) // some debug
  //     {
  //       debug_absorberEnergyInCrystal[i] = 0;
  //       debug_absorberEnergyNotInCrystal[i] = 0;
  //     }
  //   }
  //   for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //   {
  //     sname << "EnInCrystals_Absorber" << iAbs;
  //     stitle << "Energy deposited in crystals of Absorber " << iAbs;
  //     TH1F *t_histoIn = new TH1F(sname.str().c_str(),stitle.str().c_str(),bins,min,max);
  //     t_histoIn->GetXaxis()->SetTitle("MeV");
  //     abs_histoIn.push_back(t_histoIn);
  //     sname.str("");
  //     stitle.str("");
  //   }
  //   for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //   {
  //     sname << "EnTot_Absorber" << iAbs;
  //     stitle << "Total Energy in Absorber " << iAbs;
  //     TH1F *t_histoIn = new TH1F(sname.str().c_str(),stitle.str().c_str(),bins,min,max);
  //     t_histoIn->GetXaxis()->SetTitle("MeV");
  //     abs_histoTot.push_back(t_histoIn);
  //     sname.str("");
  //     stitle.str("");
  //   }
  // }
  // if(sum_cell)
  // {
  //   cellEnergyInCrystal     = new float[cells_regions.size()];
  //   cellEnergyNotInCristal     = new float[cells_regions.size()];
  //   if(debug) // some debug
  //   {
  //     debug_cellEnergyInCrystal     = new float[cells_regions.size()];
  //     debug_cellEnergyNotInCristal     = new float[cells_regions.size()];
  //   }
  //   for(int i = 0 ; i < cells_regions.size(); i++)
  //   {
  //     cellEnergyInCrystal[i] = 0;
  //     cellEnergyNotInCristal[i] = 0;
  //     if(debug) // some debug
  //     {
  //       debug_cellEnergyInCrystal[i] = 0;
  //       debug_cellEnergyNotInCristal[i] = 0;
  //     }
  //   }
  //   for(int i = 0 ; i < cells_regions.size(); i++)
  //   {
  //     sname << "EnInCrystals_Cell" << i;
  //     stitle << "Energy deposited in crystals of Cell " << i;
  //     TH1F *t_histoIn = new TH1F(sname.str().c_str(),stitle.str().c_str(),bins,min,max);
  //     t_histoIn->GetXaxis()->SetTitle("MeV");
  //     cell_histoIn.push_back(t_histoIn);
  //     sname.str("");
  //     stitle.str("");
  //   }
  //   for(int i = 0 ; i < cells_regions.size(); i++)
  //   {
  //     sname << "EnTot_Cell" << i;
  //     stitle << "Total energy deposited in Cell " << i << " (crystals + absorber + air)";
  //     TH1F *t_histoTot = new TH1F(sname.str().c_str(),stitle.str().c_str(),bins,min,max);
  //     t_histoTot->GetXaxis()->SetTitle("MeV");
  //     cell_histoTot.push_back(t_histoTot);
  //     sname.str("");
  //     stitle.str("");
  //   }
  // }
  // if(sum_cry)
  // {
  //   fibreEnergy    = new float[fibres_regions.size()];
  //   if(debug)  // some debug
  //   {
  //     debug_fibreEnergy    = new float[fibres_regions.size()];
  //   }
  //   for(int i = 0 ; i < fibres_regions.size(); i++)
  //   {
  //     fibreEnergy[i] = 0;
  //     if(debug)  // some debug
  //     {
  //       debug_fibreEnergy = 0;
  //     }
  //   }
  //   for(int i = 0 ; i < fibres_regions.size(); i++)
  //   {
  //     sname << "EnIn_Crystal" << i;
  //     stitle << "Total energy deposited in Crystal " << i ;
  //     TH1F *t_histoTot = new TH1F(sname.str().c_str(),stitle.str().c_str(),bins,min,max);
  //     t_histoTot->GetXaxis()->SetTitle("MeV");
  //     fiber_histoIn.push_back(t_histoTot);
  //     sname.str("");
  //     stitle.str("");
  //   }
  // }
  //
  // // loop on events
  // long int nEvents = shower->GetEntries();
  // int eventID = -1;
  // long long int counter = 0;
  // float dumb_sum = 0.;
  // for(int i = 0 ; i < nEvents ; i++) // loop on ttree entries
  // {
  //   bool inAbs = false;
  //   bool inCry = false;
  //   bool inDumb = false;
  //   // get the entry
  //   shower->GetEvent(i);
  //
  //   // now summation of energy needs to be performed only until a new event arrives in the list
  //   // an event corresponds to 1 primary shot. the info is stored in event
  //   // so we start be setting a event indicator, eventID, to -1. At the first loop step,
  //   // the first value of variable event will be 0, so what comes in the next if statement
  //   // (event saving and clean up) will not be executed, and instead the event accumulation starts
  //   // and the end of event accumulation the eventID becomes = event, then a new step of the cycle
  //   // starts. Now eventID is = 0, so "event saving and clean up" will not be performed until
  //   // the value of event becomes != 0. When that happens, the variables were the enegies were
  //   // summed are used to filled the relevant histograms, than they are cleared and a new accumulation
  //   // starts. Of course this means that after the end of the loop, the last event still needs to be saved
  //
  //
  //   //------------------------------//
  //   // event saving and clean up    //
  //   //------------------------------//
  //   if((eventID != event) && (eventID != -1)) // new event, dump the previous
  //   {
  //     //fill histograms
  //     if(sum_abs)
  //     {
  //       for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //       {
  //         abs_histoIn[iAbs]->Fill(absorberEnergyInCrystal[iAbs]);
  //         abs_histoTot[iAbs]->Fill(absorberEnergyInCrystal[iAbs] + absorberEnergyNotInCrystal[iAbs]);
  //       }
  //     }
  //     if(sum_cell)
  //     {
  //       h_energyInFrontInCry->Fill(energyInFrontInCry);
  //       h_energyInFrontInAbs->Fill(energyInFrontInCry + energyInFrontNotInCry);
  //       h_energyInBackInCry->Fill(energyInBackInCry);
  //       h_energyInBackInAbs->Fill(energyInBackInCry + energyInBackNotInCry);
  //       for(int iCell = 0 ; iCell < cells_regions.size(); iCell++)
  //       {
  //         cell_histoIn[iCell]->Fill(cellEnergyInCrystal[iCell]);
  //         cell_histoTot[iCell]->Fill(cellEnergyInCrystal[iCell] + cellEnergyNotInCristal[iCell]);
  //       }
  //       for(int iCell = 0 ; iCell < cells_regions.size()/2; iCell++)
  //       {
  //         histo_pair_cell_In[iCell]->Fill(pair_cellEnergyInCrystal[iCell]);
  //         histo_pair_cell_Tot[iCell]->Fill(pair_cellEnergyInCrystal[iCell] + pair_cellEnergyNotInCrystal[iCell]);
  //         histo_ratio_pair_cell_In[iCell]->Fill(cellEnergyInCrystal[iCell+9]/cellEnergyInCrystal[iCell]);
  //         histo_ratio_pair_cell_Tot[iCell]->Fill( (cellEnergyInCrystal[iCell+9] + cellEnergyNotInCristal[iCell+9])/(cellEnergyInCrystal[iCell] + cellEnergyNotInCristal[iCell]) );
  //       }
  //
  //     }
  //     if(sum_cry)
  //     {
  //       for(int iCry = 0 ; iCry < fibres_regions.size(); iCry++)
  //       {
  //         fiber_histoIn[iCry]->Fill(fibreEnergy[iCry]);
  //       }
  //     }
  //
  //     if(debug)// DEBUG output sums
  //     {
  //       std::cout << std::endl;
  //       std::cout << "EVENT " << eventID << std::endl;
  //       if(sum_abs)
  //       {
  //         for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //         {
  //           std::cout << "Absorber " << iAbs << std::endl;
  //           std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
  //                     << absorberEnergyInCrystal[iAbs] << std::endl;
  //           std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
  //                     << absorberEnergyNotInCrystal[iAbs] << std::endl;
  //         }
  //       }
  //       if(sum_cell)
  //       {
  //
  //         for(int iCell = 0 ; iCell < cells_regions.size(); iCell++)
  //         {
  //           std::cout << "Cell " << iCell << std::endl;
  //           std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
  //                     << cellEnergyInCrystal[iCell] << std::endl;
  //           std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
  //                     << cellEnergyNotInCristal[iCell] << std::endl;
  //
  //         }
  //
  //       }
  //       if(sum_cry)
  //       {
  //         float sum_in_cry = 0;
  //         for(int iCry = 0 ; iCry < fibres_regions.size(); iCry++)
  //         {
  //           // if(fibreEnergy[iCry] != 0)
  //           // {
  //           //   std::cout << fibreEnergy[iCry] << std::endl;
  //           // }
  //           sum_in_cry += fibreEnergy[iCry];
  //         }
  //         std::cout << "SUM IN CRYSTALS = " << sum_in_cry << std::endl;
  //       }
  //       std::cout << "dumb sum = " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << dumb_sum << std::endl;
  //       std::cout << std::endl;
  //     } // end of debug
  //
  //
  //
  //     //clean variables
  //     dumb_sum = 0;
  //     energyInFrontInCry = 0.0;
  //     energyInFrontNotInCry = 0.0;
  //     energyInBackInCry = 0.0;
  //     energyInBackNotInCry = 0.0;
  //
  //     if(sum_abs)
  //     {
  //       for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //       {
  //         absorberEnergyInCrystal[iAbs] = 0;
  //         absorberEnergyNotInCrystal[iAbs] = 0;
  //       }
  //     }
  //     if(sum_cell)
  //     {
  //       for(int iCell = 0 ; iCell < cells_regions.size(); iCell++)
  //       {
  //         cellEnergyInCrystal[iCell] = 0;
  //         cellEnergyNotInCristal[iCell] = 0;
  //       }
  //       for(int iCell = 0 ; iCell < cells_regions.size()/2; iCell++)
  //       {
  //         pair_cellEnergyInCrystal[iCell] = 0;
  //         pair_cellEnergyNotInCrystal[iCell] = 0;
  //       }
  //     }
  //     if(sum_cry)
  //     {
  //       for(int iCry = 0 ; iCry < fibres_regions.size(); iCry++)
  //       {
  //         fibreEnergy[iCry] = 0;
  //       }
  //     }
  //   }
  //   // end of event saving and clean up //
  //   //------------------------------//
  //
  //
  //   //------------------------------//
  //   // event accumulation           //
  //   //------------------------------//
  //
  //   if(crystalID != -1)
  //   {
  //     inDumb = true;
  //     dumb_sum += showerIonizingEnDep;
  //   }
  //
  //   if(sum_abs)
  //   {
  //     //check if event inside absorber
  //     for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //     {
  //       if( (showerX >= absobers_regions[iAbs].min[0]) &&
  //           (showerX <  absobers_regions[iAbs].max[0]) &&
  //           (showerY >= absobers_regions[iAbs].min[1]) &&
  //           (showerY <  absobers_regions[iAbs].max[1]) &&
  //           (showerZ >= absobers_regions[iAbs].min[2]) &&
  //           (showerZ <  absobers_regions[iAbs].max[2])   ) // inside the abs region
  //       {
  //         if(showerIsInCrystal) /// if in crystal
  //         {
  //           // for(int iCell = 0 ; iCell < cells_regions.size(); iCell++) // find cell
  //           // {
  //           //   if( (showerX >= cells_regions[iCell].min[0]) &&
  //           //       (showerX <  cells_regions[iCell].max[0]) &&
  //           //       (showerY >= cells_regions[iCell].min[1]) &&
  //           //       (showerY <  cells_regions[iCell].max[1]) &&
  //           //       (showerZ >= cells_regions[iCell].min[2]) &&
  //           //       (showerZ <  cells_regions[iCell].max[2])   ) // inside this cell region
  //           //   {
  //           //
  //           //     // if(showerIsInCrystal)
  //           //     // {
  //           //     //   cellEnergyInCrystal[iCell] += showerIonizingEnDep;
  //           //     // }
  //           //     // else
  //           //     // {
  //           //     //   cellEnergyNotInCristal[iCell] += showerIonizingEnDep;
  //           //     // }
  //           //   }
  //           // }
  //
  //
  //
  //           absorberEnergyInCrystal[iAbs] += showerIonizingEnDep;
  //           inAbs = true;
  //         }
  //         else
  //         {
  //           absorberEnergyNotInCrystal[iAbs] += showerIonizingEnDep;
  //         }
  //       }
  //     }
  //   }
  //   if(sum_cell)
  //   {
  //     //check if event inside cell
  //     for(int iCell = 0 ; iCell < cells_regions.size(); iCell++)
  //     {
  //       if( (showerX >= cells_regions[iCell].min[0]) &&
  //           (showerX <  cells_regions[iCell].max[0]) &&
  //           (showerY >= cells_regions[iCell].min[1]) &&
  //           (showerY <  cells_regions[iCell].max[1]) &&
  //           (showerZ >= cells_regions[iCell].min[2]) &&
  //           (showerZ <  cells_regions[iCell].max[2])   ) // inside this cell region
  //       {
  //
  //         if(showerIsInCrystal)
  //         {
  //           cellEnergyInCrystal[iCell] += showerIonizingEnDep;
  //           if(iCell > 8)
  //           {
  //             pair_cellEnergyInCrystal[iCell - 9] += showerIonizingEnDep;
  //             energyInBackInCry  += showerIonizingEnDep;
  //           }
  //           else
  //           {
  //             pair_cellEnergyInCrystal[iCell] += showerIonizingEnDep;
  //             energyInFrontInCry  += showerIonizingEnDep;
  //           }
  //         }
  //         else
  //         {
  //           cellEnergyNotInCristal[iCell] += showerIonizingEnDep;
  //           if(iCell > 8)
  //           {
  //             pair_cellEnergyNotInCrystal[iCell - 9] += showerIonizingEnDep;
  //             energyInBackNotInCry  += showerIonizingEnDep;
  //           }
  //           else
  //           {
  //             pair_cellEnergyNotInCrystal[iCell] += showerIonizingEnDep;
  //             energyInFrontNotInCry  += showerIonizingEnDep;
  //           }
  //         }
  //       }
  //     }
  //   }
  //   if(sum_cry)
  //   {
  //     //check if event inside fiber
  //     for(int iCry = 0 ; iCry < fibres_regions.size(); iCry++)
  //     {
  //       if( (showerX >= fibres_regions[iCry].min[0]) &&
  //           (showerX <= fibres_regions[iCry].max[0]) &&
  //           (showerY >= fibres_regions[iCry].min[1]) &&
  //           (showerY <= fibres_regions[iCry].max[1]) &&
  //           (showerZ >= fibres_regions[iCry].min[2]) &&
  //           (showerZ <= fibres_regions[iCry].max[2])   ) // inside this fiber region
  //       {
  //         // no need to check if it is in a crystal...
  //         inCry = true;
  //         fibreEnergy[iCry] += showerIonizingEnDep;
  //       }
  //
  //     }
  //   }
  //   // end of event accumulation
  //   //------------------------------//
  //
  //   //DEBUG
  //   if(debug)  // some debug
  //   {
  //     if(inAbs == true)
  //     {
  //       if(inCry == false)
  //       {
  //         std::cout << "found " << event << " " << i << std::endl;
  //         std::cout << showerIsInCrystal << " "
  //                   << std::setprecision(std::numeric_limits<double>::digits10 + 1) << showerX << " "
  //                   << std::setprecision(std::numeric_limits<double>::digits10 + 1) << showerY << " "
  //                   << std::setprecision(std::numeric_limits<double>::digits10 + 1) << showerZ << " "
  //                   << std::endl;
  //       }
  //     }
  //
  //     if(inAbs == false)
  //     {
  //       if(inCry == true)
  //       {
  //         std::cout << "weird " << event << " " << i << std::endl;
  //         std::cout << showerIsInCrystal << " "
  //                   << std::setprecision(std::numeric_limits<double>::digits10 + 1) << showerX << " "
  //                   << std::setprecision(std::numeric_limits<double>::digits10 + 1) << showerY << " "
  //                   << std::setprecision(std::numeric_limits<double>::digits10 + 1) << showerZ << " "
  //                   << std::endl;
  //       }
  //     }
  //   }
  //
  //
  //
  //   // update event id
  //   eventID = event;
  //   // give feedback to user
  //   counter++;
  //   int perc = ((100*counter)/nEvents);
  //   if(!debug)
  //   {
  //     if( (perc % 10) == 0 )
  //     {
  //       std::cout << "\r";
  //       std::cout << perc << "% done... ";
  //     }
  //   }
  // } // end of loop on ttree entries
  //
  // //fill last event
  // if(sum_abs)
  // {
  //   for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //   {
  //     abs_histoIn[iAbs]->Fill(absorberEnergyInCrystal[iAbs]);
  //     abs_histoTot[iAbs]->Fill(absorberEnergyInCrystal[iAbs] + absorberEnergyNotInCrystal[iAbs]);
  //   }
  // }
  // if(sum_cell)
  // {
  //   for(int iCell = 0 ; iCell < cells_regions.size(); iCell++)
  //   {
  //     cell_histoIn[iCell]->Fill(cellEnergyInCrystal[iCell]);
  //     cell_histoTot[iCell]->Fill(cellEnergyInCrystal[iCell] + cellEnergyNotInCristal[iCell]);
  //   }
  // }
  // if(sum_cry)
  // {
  //   for(int iCry = 0 ; iCry < fibres_regions.size(); iCry++)
  //   {
  //     fiber_histoIn[iCry]->Fill(fibreEnergy[iCry]);
  //   }
  // }
  //
  // if(debug)// DEBUG output sums
  // {
  //   std::cout << std::endl;
  //   std::cout << "EVENT " << eventID << std::endl;
  //   if(sum_abs)
  //   {
  //     for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //     {
  //       std::cout << "Absorber " << iAbs << std::endl;
  //       std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << absorberEnergyInCrystal[iAbs] << std::endl;
  //       std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << absorberEnergyNotInCrystal[iAbs] << std::endl;
  //     }
  //   }
  //   if(sum_cell)
  //   {
  //     for(int iCell = 0 ; iCell < cells_regions.size(); iCell++)
  //     {
  //       std::cout << "Cell " << iCell << std::endl;
  //       std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << cellEnergyInCrystal[iCell] << std::endl;
  //       std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << cellEnergyNotInCristal[iCell] << std::endl;
  //     }
  //   }
  //   if(sum_cry)
  //   {
  //     float sum_in_cry = 0;
  //     for(int iCry = 0 ; iCry < fibres_regions.size(); iCry++)
  //     {
  //       // if(fibreEnergy[iCry] != 0)
  //       // {
  //       //   std::cout << fibreEnergy[iCry] << std::endl;
  //       // }
  //       sum_in_cry += fibreEnergy[iCry];
  //     }
  //     std::cout << "SUM IN CRYSTALS = " << sum_in_cry << std::endl;
  //   }
  //
  //   std::cout << "dumb sum = " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) <<  dumb_sum << std::endl;
  //   std::cout << std::endl;
  // } // end of debug
  // // end of CALCULATE EN RES
  // //------------------------------------------------//
  //
  //
  // //------------------------------------------------//
  // // SAVE TO OUTPUT                                 //
  // //------------------------------------------------//
  // std::string outputRootFileName = outputFileName + ".root";
  // std::string outputTextFileName = outputFileName + ".txt";
  // // save to text file
  // std::ofstream outputTextFile;
  // std::cout <<std::endl;
  // std::cout << "Saving results in "<< outputTextFileName.c_str() << std::endl ;
  // outputTextFile.open (outputTextFileName.c_str(),std::ofstream::out);
  // outputTextFile << "Abs Number\t FIT Fraction Tot En[MeV]\t FIT Energy Res. In Cry (sigma/E) \tFraction Tot En[MeV]\t Energy Res. In Cry (sigma/E)" << std::endl;
  // if(sum_abs)
  // {
  //   // fit total and crystal abs spectra
  //   for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //   {
  //     TF1 *gaussIn = new TF1("gaussIn","gaus",min,max);
  //     abs_histoIn[iAbs]->Fit(gaussIn,"Q");
  //
  //     TF1 *gaussTot = new TF1("gaussTot","gaus",min,max);
  //     abs_histoTot[iAbs]->Fit(gaussTot,"Q");
  //
  //     outputTextFile << iAbs << "\t"
  //                    << gaussTot->GetParameter(1) / max << "\t"
  //                    << gaussIn->GetParameter(2) / gaussIn->GetParameter(1) << "\t"
  //                    << abs_histoTot[iAbs]->GetMean() / max << "\t"
  //                    << abs_histoIn[iAbs]->GetRMS() / abs_histoIn[iAbs]->GetMean() << "\t"
  //                    << std::endl;
  //   }
  // }
  // outputTextFile.close();
  // std::cout <<std::endl;
  // // save to ROOT file
  // std::cout << "Saving histograms to "<< outputRootFileName.c_str() << " ..." ;
  // TFile *outputFile = new TFile(outputRootFileName.c_str(),"RECREATE");
  // outputFile->cd();
  // if(sum_abs)
  // {
  //   for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //   {
  //     abs_histoIn[iAbs]->Write();
  //     abs_histoTot[iAbs]->Write();
  //   }
  // }
  // if(sum_cell)
  // {
  //   h_energyInFrontInCry->Write();
  //   h_energyInFrontInAbs->Write();
  //   h_energyInBackInCry->Write();
  //   h_energyInBackInAbs->Write();
  //   for(int iCell = 0 ; iCell < cells_regions.size(); iCell++)
  //   {
  //     cell_histoIn[iCell]->Write();
  //     cell_histoTot[iCell]->Write();
  //   }
  //   for(int iCell = 0 ; iCell < cells_regions.size()/2; iCell++)
  //   {
  //     histo_pair_cell_In[iCell]->Write();
  //     histo_pair_cell_Tot[iCell]->Write();
  //     histo_ratio_pair_cell_In[iCell]->Write();
  //     histo_ratio_pair_cell_Tot[iCell]->Write();
  //   }
  //
  // }
  // if(sum_cry)
  // {
  //   for(int iCry = 0 ; iCry < fibres_regions.size(); iCry++)
  //   {
  //     fiber_histoIn[iCry]->Write();
  //   }
  // }
  // outputFile->Close();
  // // end of SAVE TO OUTPUT
  // //------------------------------------------------//
  //
  //
  // std::cout <<" done!" << std::endl;
