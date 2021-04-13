// Marco Pizzichemi 19.10.2020 marco.pizzichemi@cern.ch

// compile with
// g++ -o ../build/showerContainment showerContainment.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore

// calculates 

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
  int type;
  float minX;
  float maxX;
  float minY;
  float maxY;
  float minZ;
  float maxZ;
  float total_in_crystal;
  float total_in_module;
  
  float fraction;
  float fraction_in_front;
  float fraction_in_back;

  std::vector<float> v_fraction;
  std::vector<float> v_fraction_in_front;
  std::vector<float> v_fraction_in_back;

  std::vector<float> v_sampling_fraction;
  std::vector<float> v_total_in_module;
  std::vector<float> v_total_in_front;
  std::vector<float> v_total_in_back;
  std::vector<float> v_total_in_crystal;
  TH1F *hTotal   ;
  TH1F *hFront   ;
  TH1F *hBack   ;
  TH1F *hFraction;
  TH1F *hFractionFront;
  TH1F *hFractionBack;

  // TH1F *h_EnDepTotal       ;
};

void usage()
{
  std::cout << "\t\t" << "[-f | --folder]  <input folder>   " << std::endl
            << "\t\t" << "[-i | --input]   <input file prefix>   " << std::endl
            << "\t\t" << "[-o | --output]  <output file name>   " << std::endl
            << "\t\t" << "[-e | --energy]  <primary energy in MeV - default = 1000>   " << std::endl
            << "\t\t" << "[--type]         <module type - default = 0 >  " << std::endl            
            << "\t\t" << "[--stop]         <stop after N event, default = INFINITY >  " << std::endl
            << "\t\t" << "[--verbose]      <verbosity level - default = 0>  " << std::endl
            << "\t\t" << "[--limitx]      <user defined limit (abs) for energy summing - default = ignored>  " << std::endl
            << "\t\t" << "[--limity]      <user defined limit (abs) for energy summing - default = ignored>  " << std::endl
            << "\t\t" << "[--text]        <output is only text (no ROOT file)>  " << std::endl
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
  int maxEvent = std::numeric_limits<int>::max();
  float energy = 1000.;
  float limitx = std::numeric_limits<float>::max();
  float limity = std::numeric_limits<float>::max();
  bool useLimits = false;
  bool onlyText = false;
  static struct option longOptions[] =
  {
			{ "input", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { "folder", required_argument, 0, 0 },
      { "verbose", required_argument, 0, 0 },
      { "type", required_argument, 0, 0 },
      { "stop", required_argument, 0, 0 },
      { "energy", required_argument, 0, 0 },
      { "limitx", required_argument, 0, 0 },
      { "limity", required_argument, 0, 0 },
      { "text", no_argument, 0, 0 },
			{ NULL, 0, 0, 0 }
	};
  // read input
  while(1) {
		int optionIndex = 0;
		int c = getopt_long(argc, argv, "i:o:f:e:", longOptions, &optionIndex);
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
    else if (c == 'e'){
      energy = atof((char *)optarg);
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
      maxEvent = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 6){
      energy = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 7){
      useLimits = true;
      limitx = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 8){
      useLimits = true;
      limity = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 9){
      onlyText = true;
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
  if((outputFileName == "") && (onlyText == false) )
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
  TFile *outputFile;
  if(!onlyText) outputFile = new TFile(outputFileName.c_str(),"RECREATE");

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
  std::vector<float> modulePositionX;
  std::vector<float> modulePositionY;
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
      modulePositionX.push_back(modules_regions[iMod].position_x);
      modulePositionY.push_back(modules_regions[iMod].position_y);
      modulePositionZ.push_back(modules_regions[iMod].position_z);
    }
  }

  // lots of maps
  std::map<int, int> moduleNumberMap;         // map from module type to number of modules
  std::map<int, int> moduleNumberOfBeforeMap; // map from module type to number of other modules
  std::map<int, int> moduleSectionsMap;       // map from module type to number of segments
  std::map<int, float> moduleXMinMap;          // map from module type to min X
  std::map<int, float> moduleXMaxMap;          // map from module type to max X
  std::map<int, float> moduleYMinMap;          // map from module type to min Y
  std::map<int, float> moduleYMaxMap;          // map from module type to max Y
  std::map<int, float> moduleZMinMap;          // map from module type to min z
  std::map<int, float> moduleZMaxMap;          // map from module type to max z
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
        moduleXMaxMap.insert(std::make_pair(moduleTypes[i], cells_regions[iCell].max[0] + modulePositionX[i]));
        moduleYMaxMap.insert(std::make_pair(moduleTypes[i], cells_regions[iCell].max[1] + modulePositionY[i]));
        moduleXMinMap.insert(std::make_pair(moduleTypes[i], cells_regions[iCell].min[0] + modulePositionX[i]));
        moduleYMinMap.insert(std::make_pair(moduleTypes[i], cells_regions[iCell].min[1] + modulePositionY[i]));
        if (moduleSectionsMap.at(moduleTypes[i]) == 1)
        {
          // only 1 sections, just take max and min
          foundLimit[i][1] = true;
          moduleZMaxMap.insert(std::make_pair(moduleTypes[i], cells_regions[iCell].max[2] + modulePositionZ[i]));
          foundLimit[i][0] = true; 
          moduleZMinMap.insert(std::make_pair(moduleTypes[i], cells_regions[iCell].min[2] + modulePositionZ[i]));
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
              moduleZMaxMap.insert(std::make_pair(moduleTypes[i], cells_regions[iCell].max[2] + modulePositionZ[i]));
            }
          }
          else
          {
            if (foundLimit[i][0] == false)
            {
              foundLimit[i][0] = true;
              moduleZMinMap.insert(std::make_pair(moduleTypes[i], cells_regions[iCell].min[2] + modulePositionZ[i]));
            }
          }
        }
      }
    }
  }
  
  std::vector<float> all_x;
  std::vector<float> all_y;
  std::vector<float> all_z;

  for (int i = 0; i < moduleTypes.size(); i++)
  {
    feedbackString << "Type " << moduleTypes[i] << " has " << moduleSectionsMap.at(moduleTypes[i]) << " sections " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " has " << moduleNumberMap.at(moduleTypes[i]) << " modules " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " has " << moduleNumberOfBeforeMap.at(moduleTypes[i]) << " modules before" << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " X min " << moduleXMinMap.at(moduleTypes[i]) << " mm " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " X max " << moduleXMaxMap.at(moduleTypes[i]) << " mm " << std::endl;
    all_x.push_back(moduleXMinMap.at(moduleTypes[i]));
    all_x.push_back(moduleXMaxMap.at(moduleTypes[i]));
    feedbackString << "Type " << moduleTypes[i] << " Y min " << moduleYMinMap.at(moduleTypes[i]) << " mm " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " Y max " << moduleYMaxMap.at(moduleTypes[i]) << " mm " << std::endl;
    all_y.push_back(moduleYMinMap.at(moduleTypes[i]));
    all_y.push_back(moduleYMaxMap.at(moduleTypes[i]));
    feedbackString << "Type " << moduleTypes[i] << " Z min " << moduleZMinMap.at(moduleTypes[i]) << " mm " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " Z max " << moduleZMaxMap.at(moduleTypes[i]) << " mm " << std::endl;
    all_z.push_back(moduleZMinMap.at(moduleTypes[i]));
    all_z.push_back(moduleZMaxMap.at(moduleTypes[i]));
    feedbackString << "Type " << moduleTypes[i] << " separation at " << moduleSeparationMap.at(moduleTypes[i]) << " mm " << std::endl;
  }


  
  
  
  // // finally calc the whole module min and max x,y,z
  // C++17 way...
  // const auto [minX, maxX] = std::minmax_element(begin(all_x), end(all_x));
  // const auto [minY, maxY] = std::minmax_element(begin(all_y), end(all_y));
  // const auto [minZ, maxZ] = std::minmax_element(begin(all_z), end(all_z));
  // std::cout << "minX = " << *minX << ", maxX = " << *maxX << '\n';
  // std::cout << "minY = " << *minY << ", maxY = " << *maxY << '\n';
  // std::cout << "minZ = " << *minZ << ", maxZ = " << *maxZ << '\n';
  // more standard way...
  float minX = *min_element(all_x.begin(), all_x.end());
  float maxX = *max_element(all_x.begin(), all_x.end());
  float minY = *min_element(all_y.begin(), all_y.end());
  float maxY = *max_element(all_y.begin(), all_y.end());
  float minZ = *min_element(all_z.begin(), all_z.end());
  float maxZ = *max_element(all_z.begin(), all_z.end());
  feedbackString << "Module minX = " << minX << ", maxX = " << maxX << '\n';
  feedbackString << "Module minY = " << minY << ", maxY = " << maxY << '\n';
  feedbackString << "Module minZ = " << minZ << ", maxZ = " << maxZ << '\n';

  
  // prepare arrays etc for each module type (regardless of some module being actually the same)
  std::map<int, module_t> mapOfResultsArray;
  std::map<int, float> mapOfSamplingFraction;
  std::map<int, float> mapOfTotalEnergyInside;

  for (int i = 0; i < moduleTypes.size(); i++)
  {
    module_t tempModule;
    tempModule.type = moduleTypes[i];
    tempModule.minX = moduleXMinMap.at(moduleTypes[i]);
    tempModule.maxX = moduleXMaxMap.at(moduleTypes[i]);
    tempModule.minY = moduleYMinMap.at(moduleTypes[i]);
    tempModule.maxY = moduleYMaxMap.at(moduleTypes[i]);
    tempModule.minZ = moduleZMinMap.at(moduleTypes[i]);
    tempModule.maxZ = moduleZMaxMap.at(moduleTypes[i]);
    tempModule.total_in_crystal = 0.;
    tempModule.total_in_module = 0.;
    tempModule.fraction = 0.;
    tempModule.fraction_in_front = 0.;
    tempModule.fraction_in_back = 0.;
    mapOfResultsArray.insert(std::make_pair(moduleTypes[i], tempModule));
    mapOfSamplingFraction.insert(std::make_pair(moduleTypes[i],0.));
    mapOfSamplingFraction.insert(std::make_pair(moduleTypes[i],0.));
  }
  for (int i = 0; i < moduleTypes.size(); i++)
  {
    feedbackString << moduleTypes[i] << " "
                   << (mapOfResultsArray.at(moduleTypes[i])).type << " "
                   << (mapOfResultsArray.at(moduleTypes[i])).minX << " "
                   << (mapOfResultsArray.at(moduleTypes[i])).maxX << " "
                   << (mapOfResultsArray.at(moduleTypes[i])).minY << " "
                   << (mapOfResultsArray.at(moduleTypes[i])).maxY << " "
                   << (mapOfResultsArray.at(moduleTypes[i])).minZ << " "
                   << (mapOfResultsArray.at(moduleTypes[i])).maxZ << " "
                   << std::endl;
  }



  feedbackString << std::endl;
  std::cout << feedbackString.str() << std::endl;
  

  // outputFile->cd();
  




  //----------------------------------//
  // LOOP ON INPUT                    //
  //----------------------------------//
  TFile *inputFile = ((TFile *)0);
  TTree *shower = ((TTree*) 0);
  // now run on all input files
  float fractionInsideModule = 0;
  float totalEnergyInCrystals = 0;
  float moduleSeparation;
  
  std::vector<float> v_total_inside;
  std::vector<float> v_total_in_crystals;
  std::vector<float> v_fraction;

  for(unsigned int i = 0 ; i < listInputFiles.size(); i++)
  {
    if (inputFile) delete inputFile; // just a precaution

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
    feedbackDivision = (int) (nEntriesShower/20);
    if(feedbackDivision == 0) feedbackDivision = 1;
    // float en_at_surf = 0;
    
    for(int iEntry = 0 ; iEntry < nEntriesShower ; iEntry++)
    {
      
      shower->GetEvent(iEntry);

      // std::cout << iEntry  << " " 
      //           << eventID << " " 
      //           << showerEvent << std::endl;

      if(showerEvent == maxEvent )
        break;

      if(singleModuleType)
      {
        if (showerModuleType != moduleTypeNumber)
        {
          continue;
        }
      }
      
      // skip event if outside user limits
      if(useLimits)
      {
        if((fabs(showerX) > limitx) || (fabs(showerY) > limity))
        {
          continue;
        }
      }

      // if(iEntry == 0)
      // {
      //   en_at_surf = shower_primaryEnergyOnAbsorber;
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
      if((eventID != showerEvent) && (eventID != -1)) // new event, dump the previous, and store the primaries info
      {
        // fill
        // std::cout << "------------------------------------------ " ;
        if(fractionInsideModule > 1) std::cout <<  "WARNING - " << fractionInsideModule << " " <<eventID << std::endl;
        v_fraction.push_back(fractionInsideModule);
        v_total_in_crystals.push_back(totalEnergyInCrystals);

        // single modules 
        for (int i = 0; i < moduleTypes.size(); i++)
        {
          (mapOfResultsArray.at(moduleTypes[i])).v_total_in_module .push_back((mapOfResultsArray.at(moduleTypes[i])).total_in_module );
          (mapOfResultsArray.at(moduleTypes[i])).v_total_in_crystal.push_back((mapOfResultsArray.at(moduleTypes[i])).total_in_crystal);
          (mapOfResultsArray.at(moduleTypes[i])).v_fraction.push_back( (mapOfResultsArray.at(moduleTypes[i])).fraction);
          (mapOfResultsArray.at(moduleTypes[i])).v_fraction_in_front .push_back((mapOfResultsArray.at(moduleTypes[i])).fraction_in_front );
          (mapOfResultsArray.at(moduleTypes[i])).v_fraction_in_back  .push_back((mapOfResultsArray.at(moduleTypes[i])).fraction_in_back  );
        }
        
        

        // cleanup
        fractionInsideModule = 0;
        totalEnergyInCrystals = 0;
        
        for (int i = 0; i < moduleTypes.size(); i++)
        {
          (mapOfResultsArray.at(moduleTypes[i])).total_in_module = 0;
          (mapOfResultsArray.at(moduleTypes[i])).total_in_crystal = 0;
          (mapOfResultsArray.at(moduleTypes[i])).fraction = 0;
          (mapOfResultsArray.at(moduleTypes[i])).fraction_in_front = 0;
          (mapOfResultsArray.at(moduleTypes[i])).fraction_in_back = 0;
        }

      }
      // end of event saving 

      // event accumulating
      if( (showerX >= minX) && (showerX <= maxX) &&
          (showerY >= minY) && (showerY <= maxY) &&
          (showerZ >= minZ) && (showerZ <= maxZ) 
        )
      {
        
        

        float moduleSeparation = 0.;
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
        
        fractionInsideModule += showerTotalEnDep/energy;
        // totalInsideModule += showerTotalEnDep;
        (mapOfResultsArray.at(showerModuleType)).fraction += showerTotalEnDep/energy;
        if(showerZ <= moduleSeparation)
        {
          (mapOfResultsArray.at(showerModuleType)).fraction_in_front += showerTotalEnDep/energy;
        }
        else 
        {
          (mapOfResultsArray.at(showerModuleType)).fraction_in_back += showerTotalEnDep/energy;
        }
        

        // and now fill proper submodule...
        // avoid moduleType == 0, should be only on borders if multi modules are defined...
        if(moduleTypes.size() > 1)
        {
          if(!showerModuleType) continue;
        }

        (mapOfResultsArray.at(showerModuleType)).total_in_module+=showerTotalEnDep;
        
        if(showerIsInCrystal)
        {
          totalEnergyInCrystals+=showerTotalEnDep;
          (mapOfResultsArray.at(showerModuleType)).total_in_crystal+=showerTotalEnDep;
        }
      }

      
      
      eventID = showerEvent;
    }
    
    // fill last event!
    // v_total_inside.push_back(totalEnergyInsideModule);
    // v_energy_at_surface.push_back(en_at_surf);
    // std::cout << fractionInsideModule << std::endl;
    v_fraction.push_back(fractionInsideModule);
    v_total_in_crystals.push_back(totalEnergyInCrystals);
    for (int i = 0; i < moduleTypes.size(); i++)
    {
      (mapOfResultsArray.at(moduleTypes[i])).v_total_in_module .push_back((mapOfResultsArray.at(moduleTypes[i])).total_in_module );
      (mapOfResultsArray.at(moduleTypes[i])).v_total_in_crystal.push_back((mapOfResultsArray.at(moduleTypes[i])).total_in_crystal);
      (mapOfResultsArray.at(moduleTypes[i])).v_fraction.push_back((mapOfResultsArray.at(moduleTypes[i])).fraction);
      (mapOfResultsArray.at(moduleTypes[i])).v_fraction_in_front .push_back((mapOfResultsArray.at(moduleTypes[i])).fraction_in_front );
      (mapOfResultsArray.at(moduleTypes[i])).v_fraction_in_back  .push_back((mapOfResultsArray.at(moduleTypes[i])).fraction_in_back  );
    }

    // clean last event
    fractionInsideModule = 0;
    totalEnergyInCrystals = 0;
    for (int i = 0; i < moduleTypes.size(); i++)
    {
      (mapOfResultsArray.at(moduleTypes[i])).total_in_module = 0;
      (mapOfResultsArray.at(moduleTypes[i])).total_in_crystal = 0;
      (mapOfResultsArray.at(moduleTypes[i])).fraction = 0;
      (mapOfResultsArray.at(moduleTypes[i])).fraction_in_front = 0;
      (mapOfResultsArray.at(moduleTypes[i])).fraction_in_back = 0;
    }
    
  }
  
  
  // fractions in cells 
// and for modules...
  for (int i = 0; i < moduleTypes.size(); i++)
  {
    // from vector to value
    float avg_fraction = 0. ;
    float std_fraction  = 0. ;
    float sstd_fraction  = 0. ;
    unsigned int n = (mapOfResultsArray.at(moduleTypes[i])).v_fraction.size();
    if ( n != 0) {
       avg_fraction = accumulate( (mapOfResultsArray.at(moduleTypes[i])).v_fraction.begin(), (mapOfResultsArray.at(moduleTypes[i])).v_fraction.end(), 0.0) / n; 
       double sq_sum = std::inner_product((mapOfResultsArray.at(moduleTypes[i])).v_fraction.begin(), (mapOfResultsArray.at(moduleTypes[i])).v_fraction.end(), (mapOfResultsArray.at(moduleTypes[i])).v_fraction.begin(), 0.0);
       std_fraction = std::sqrt(sq_sum / (mapOfResultsArray.at(moduleTypes[i])).v_fraction.size() - avg_fraction * avg_fraction);
       sstd_fraction = std_fraction / sqrt(n-1);
    }

    float avg_fraction_front = 0. ;
    float std_fraction_front  = 0. ;
    float sstd_fraction_front  = 0. ;
    unsigned int n_front = (mapOfResultsArray.at(moduleTypes[i])).v_fraction_in_front.size();
    if ( n_front != 0) {
       avg_fraction_front = accumulate( (mapOfResultsArray.at(moduleTypes[i])).v_fraction_in_front.begin(), (mapOfResultsArray.at(moduleTypes[i])).v_fraction_in_front.end(), 0.0) / n_front; 
       double sq_sum = std::inner_product((mapOfResultsArray.at(moduleTypes[i])).v_fraction_in_front.begin(), (mapOfResultsArray.at(moduleTypes[i])).v_fraction_in_front.end(), (mapOfResultsArray.at(moduleTypes[i])).v_fraction_in_front.begin(), 0.0);
       std_fraction_front = std::sqrt(sq_sum / (mapOfResultsArray.at(moduleTypes[i])).v_fraction_in_front.size() - avg_fraction_front * avg_fraction_front);
       sstd_fraction_front = std_fraction_front / sqrt(n_front-1);
    }
    
    float avg_fraction_back = 0. ;
    float std_fraction_back  = 0. ;
    float sstd_fraction_back  = 0. ;
    unsigned int n_back = (mapOfResultsArray.at(moduleTypes[i])).v_fraction_in_back.size();
    if ( n_back != 0) {
       avg_fraction_back = accumulate( (mapOfResultsArray.at(moduleTypes[i])).v_fraction_in_back.begin(), (mapOfResultsArray.at(moduleTypes[i])).v_fraction_in_back.end(), 0.0) / n_back; 
       double sq_sum = std::inner_product((mapOfResultsArray.at(moduleTypes[i])).v_fraction_in_back.begin(), (mapOfResultsArray.at(moduleTypes[i])).v_fraction_in_back.end(), (mapOfResultsArray.at(moduleTypes[i])).v_fraction_in_back.begin(), 0.0);
       std_fraction_back = std::sqrt(sq_sum / (mapOfResultsArray.at(moduleTypes[i])).v_fraction_in_back.size() - avg_fraction_back * avg_fraction_back);
       sstd_fraction_back = std_fraction_back / sqrt(n_back-1);
    }



    std::cout << (mapOfResultsArray.at(moduleTypes[i])).type << "\t"
              << avg_fraction << "\t"
              << sstd_fraction << "\t"
              <<  avg_fraction_back << "\t"
              << sstd_fraction_back << "\t"
              <<  avg_fraction_front << "\t"
              << sstd_fraction_front << "\t"
              

              << std::endl;
  }

  // from vector to value
  float avg_fraction = 0. ;
  float std_fraction  = 0. ;
  float sstd_fraction  = 0. ;
  unsigned int n = v_fraction.size();
  if ( n != 0) {
     avg_fraction = accumulate( v_fraction.begin(), v_fraction.end(), 0.0) / n; 
     double sq_sum = std::inner_product(v_fraction.begin(), v_fraction.end(), v_fraction.begin(), 0.0);
     std_fraction = std::sqrt(sq_sum / v_fraction.size() - avg_fraction * avg_fraction);
     sstd_fraction = std_fraction / sqrt(n-1);
  }


  // from vector to value
  float avg_total = 0. ;
  float std_total  = 0. ;
  float sstd_total  = 0. ;
  unsigned int nt = v_total_in_crystals.size();
  if ( nt != 0) {
     avg_total = accumulate( v_total_in_crystals.begin(), v_total_in_crystals.end(), 0.0) / nt; 
     double sq_sum = std::inner_product(v_total_in_crystals.begin(), v_total_in_crystals.end(), v_total_in_crystals.begin(), 0.0);
     std_total = std::sqrt(sq_sum / v_total_in_crystals.size() - avg_total * avg_total);
     sstd_total = std_total / sqrt(nt-1);
  }
  
  // keep assuming the histogram to go until 100 GeV, but calculate appropriate bin size 
  // by assuming 100 bins in the +/- 5 stdevs range
  float max = 100000.;
  int bins = 100;
  float range_total = (avg_total + 4.0*std_total) - (avg_total - 4.0*std_total);
  float range_fraction = (avg_fraction + 4.0*std_fraction) - (avg_fraction - 4.0*std_fraction);
  float bin_size_total = range_total / bins;
  float bin_size_fraction = range_fraction / bins;
  int bins_total = (int) std::round(max / bin_size_total); 
  int bins_fraction = (int) std::round(1. / bin_size_fraction); 

  // safety 
  if(bins_total < 1)
  {
    bins_total = 100;
  }
  if(bins_fraction < 1)
  {
    bins_fraction = 100;
  }


  TH1F *hTotal = new TH1F("hTotal","hTotal",bins_total,0,max);
  TH1F *hFraction = new TH1F("hFraction","hFraction",bins_fraction,0,1);
  
  for(unsigned int i = 0; i < v_total_in_crystals.size(); i++)
  {
    hTotal->Fill(v_total_in_crystals[i]);
  }
  for(unsigned int i = 0; i < v_fraction.size(); i++)
  {
    hFraction->Fill(v_fraction[i]);
  }
  
  
  
  float en_res = hTotal->GetRMS() / hTotal->GetMean();
  float en_res_err = ( TMath::Sqrt(TMath::Power( (hTotal)->GetStdDevError() / (hTotal)->GetStdDev() ,2) + TMath::Power( (hTotal)->GetMeanError() / (hTotal)->GetMean(),2)  ) * (hTotal)->GetStdDev()/(hTotal)->GetMean());
  
  if(!onlyText)
  {
    std::cout << "Writing output to file " << outputFileName << std::endl;
    outputFile->cd();
    hTotal->Write();
    hFraction->Write();
    std::cout << "Done." << std::endl;
    outputFile->Close();
  }

  

  std::string outTxtName = listInputFiles[0].substr(0, listInputFiles[0].size()-5) + ".txt";
  std::ofstream ofs;
  ofs.open (outTxtName, std::ofstream::out);

  for (int i = 0; i < moduleTypes.size(); i++)
  {
    // calc all the sampling fractions per module
    for(int iVal = 0; iVal < (mapOfResultsArray.at(moduleTypes[i])).v_total_in_module.size(); iVal++)
    {
      (mapOfResultsArray.at(moduleTypes[i])).v_sampling_fraction.push_back((mapOfResultsArray.at(moduleTypes[i])).v_total_in_crystal[iVal]/(mapOfResultsArray.at(moduleTypes[i])).v_total_in_module[iVal]);
    }
    float avg_mod = 0. ;
    float std_mod  = 0. ;
    float sstd_mod  = 0. ;
    unsigned int n = (mapOfResultsArray.at(moduleTypes[i])).v_sampling_fraction.size();
    if ( n != 0) {
      avg_mod = accumulate( (mapOfResultsArray.at(moduleTypes[i])).v_sampling_fraction.begin(), (mapOfResultsArray.at(moduleTypes[i])).v_sampling_fraction.end(), 0.0) / n; 
      double sq_sum = std::inner_product( (mapOfResultsArray.at(moduleTypes[i])).v_sampling_fraction.begin(), (mapOfResultsArray.at(moduleTypes[i])).v_sampling_fraction.end(), (mapOfResultsArray.at(moduleTypes[i])).v_sampling_fraction.begin(), 0.0);
      std_mod = std::sqrt(sq_sum / (mapOfResultsArray.at(moduleTypes[i])).v_sampling_fraction.size() - avg_mod * avg_mod);
      sstd_mod = std_mod / sqrt(n-1);
    }
    std::cout << "Sampling fraction module type "
              << moduleTypes[i] << " [avg sstd] = "
              << avg_mod << " " 
              << sstd_mod << std::endl;
    ofs       << "Sampling fraction module type "
              << moduleTypes[i] << " [avg sstd] = "
              << avg_mod << " " 
              << sstd_mod << std::endl;
  }

  
  ofs.close();

  

  std::cout << avg_fraction << " " << sstd_fraction << std::endl;

  return 0;
}
