// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch

// compile with
// g++ -o ../build/propagateHybrid propagateHybrid.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore


// program to perform the hybrid propagation of optical photons


#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3I.h"
#include "TH1I.h"
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


#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <getopt.h>
#include <set>

#include <sys/stat.h>
#include <dirent.h>

#include "regions.h"
#include "Data.cpp"


// structs
struct calibration_point_t
{
  // float x;
  // float y;
  // float z;
  // int   s;
  // float energy;
  // int primaries;
  float eff;
  TH1I *delay;
};

struct calibration_t
{
  int type;
  std::string fileName;
  // int sections;
  // float zSeparation;
  // float region_min_z;
  // float region_maz_z;

  calibration_point_t *****calibration_point;
  float xmin  ;
  float xstep ;
  float ymin  ;
  float ystep ;
  float zmin  ;
  float zstep ;
  float emin  ;
  float estep ;
  int   nSides;
  int xn;
  int yn;
  int zn;
  int en;

  std::vector <float> x;
  std::vector <float> y;
  std::vector <float> z;
  std::vector <float> energy;

  bool footer;
  float moduleCenter;
  float showerShift;

  relevant_data_t parameters;
};


// methods
void usage();
void importCalibrationFile(std::string calibrationFileName,calibration_t &calibration);
void trim( std::string& st );


// find index of closest element in sorted vector
// taken from
// https://stackoverflow.com/questions/698520/search-for-nearest-value-in-an-array-of-doubles-in-c/701141#701141
template <typename BidirectionalIterator, typename T>
BidirectionalIterator getClosest(BidirectionalIterator first,
                                 BidirectionalIterator last,
                                 const T & value)
{
    BidirectionalIterator before = std::lower_bound(first, last, value);

    if (before == first) return first;
    if (before == last)  return --last; // iterator must be bidirectional

    BidirectionalIterator after = before;
    --before;

    return (*after - value) < (value - *before) ? after : before;
}

template <typename BidirectionalIterator, typename T>
std::size_t getClosestIndex(BidirectionalIterator first,
                            BidirectionalIterator last,
                            const T & value)
{
    return std::distance(first, getClosest(first, last, value));
}

int getPositionOfLevel(std::vector<float> myarray, float level)
{
    unsigned int array_length = myarray.size();
    return getClosestIndex(myarray.begin(), myarray.begin() + array_length, level);
}

// main
int main(int argc, char** argv)
{

  if(argc < 2)
  {
    std::cout << argv[0] << std::endl;
    usage();
    return 1;
  }

  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  std::string PWDstring(cwd);

  // save command line string
  std::stringstream streamCommandLine;
  for (int i = 0; i < argc ; i++)
  {
    streamCommandLine << argv[i] << " ";
  }
  std::string CommandLineString(streamCommandLine.str());


  std::string inputFileName = "";

  std::string outputFileName = "output.root";
  std::string externalFileName = "";
  float ScintillationYield = 0.0 ; // GaGG_Ce_Mg
  float ResolutionScale = 1.0;   // for now, but needs to be changed!

  UInt_t photonSeed = 0;

  int strategy = 0;
  bool useExternal = false;
  std::vector<int> typeNumber;
  std::string typeString = "";
  std::vector<std::string> calibrationFileName;
  std::string calibrationFileNameString = "";
  int verbose = 0;
  // bool repoFolderGiven = false;
  // std::string repoFolder = "";

  static struct option longOptions[] =
  {
    { "input", required_argument, 0, 0 },
    { "calibration", required_argument, 0, 0 },
    { "output", required_argument, 0, 0 },
    { "photonSeed", required_argument, 0, 0 },
    { "external", required_argument, 0, 0 },
    { "type", required_argument, 0, 0 },
    { "verbose", required_argument, 0, 0 },
    // { "repo", required_argument, 0, 0 },
    { NULL, 0, 0, 0 }
  };

  while(1) {
    int optionIndex = 0;
    int c = getopt_long(argc, argv, "i:c:o:p:", longOptions, &optionIndex);
    if (c == -1) {
      break;
    }
    // if (c == 'f'){
    //   inputFolderName = (char *)optarg;
    // }
    else if (c == 'i'){
      inputFileName = (char *)optarg;
    }
    else if (c == 'c'){
      calibrationFileNameString = (char *)optarg;
    }
    else if (c == 'o'){
      outputFileName = (char *)optarg;
    }
    else if (c == 'p'){
      photonSeed = atoi((char *)optarg);
    }
    // else if (c == 0 && optionIndex == 0){
    //   inputFolderName = (char *)optarg;
    // }
    else if (c == 0 && optionIndex == 0){
      inputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 1){
      calibrationFileNameString = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 2){
      outputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 3){
      photonSeed = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 4){
      useExternal = true;
      externalFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 5){
      typeString = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 6){
      verbose = atoi((char *)optarg);
    }
    // else if (c == 0 && optionIndex == 7){
    //   repoFolder = (char *)optarg;
    //   repoFolderGiven = true;
    // }
    // else if (c == 0 && optionIndex == 4){
    //   pdfBins = atoi((char *)optarg);
    // }
    // else if (c == 0 && optionIndex == 5){
    //   pdfEnd = atof((char *)optarg);
    // }
    else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
      usage();
      return 1;
    }
  }

  if(inputFileName == "")
  {
    std::cout << std::endl;
    std::cout << "ERROR! You need to provide an input file name!" << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }
  if(calibrationFileNameString == "")
  {
    std::cout << std::endl;
    std::cout << "ERROR! You need to provide at least one optical calibration file name!" << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }



  std::stringstream feedbackString;
  feedbackString << "|***************************************************************|" << std::endl;
  feedbackString << "| HYBRID MONTE CARLO FOR SPACAL-RD                              |" << std::endl;
  feedbackString << "| PROPAGATION OF SCINTILLATION AND CERENKOV                     |" << std::endl;
  feedbackString << "|***************************************************************|" << std::endl;
  feedbackString << "" << std::endl;
  feedbackString << "COMMAND LINE PARAMETERS:" << std::endl;
  feedbackString << "Input file                      = " << inputFileName <<std::endl;
  feedbackString << "Optical calibration file        = " << calibrationFileNameString <<std::endl;
  feedbackString << "Type                            = " << typeString <<std::endl;
  feedbackString << "Output file                     = " << outputFileName <<std::endl;
  feedbackString << "Random seed                     = " << photonSeed << " \t (if = 0, seed automatically computed via a TUUID object)" <<std::endl;
  if(useExternal)
  {
    feedbackString << "Using external material from    = " << externalFileName << std::endl;
  }
  // if(repoFolderGiven)
  // {
  //   feedbackString << "Scanning calibration library (if necessary) = " << repoFolder << std::endl;
  // }
  feedbackString << "" << std::endl;
  // trandom
  TRandom3 *rand = new TRandom3(photonSeed);
  UInt_t randomSeed = rand->GetSeed();


  // calc regions
  feedbackString << "|***************************************************************|" << std::endl;
  feedbackString << "| READING SPACAL STRUCTURE                                      |" << std::endl;
  feedbackString << "|***************************************************************|" << std::endl;
  feedbackString << "" << std::endl;
  TFile *inputFile = TFile::Open(inputFileName.c_str());
  TTree *modules = (TTree*)   inputFile->Get("modules");
  TTree *absorbers = (TTree*) inputFile->Get("absorbers");
  TTree *cells = (TTree*)     inputFile->Get("cells");
  TTree *holes = (TTree*)     inputFile->Get("holes");
  TTree *fibres = (TTree*)    inputFile->Get("fibres");
  TTree* shower = (TTree*)    inputFile->Get("shower");

  std::map<int, int> relevantParametersMap; // map from type to index of relevant parameters vector
  // load also the parameters ttree, if it's there
  TTree* parametersTree =NULL;
  parametersTree = (TTree*)    inputFile->Get("parameters");
  std::vector<relevant_data_t> parameters;
  bool parametersFound = false;
  // and read...
  if(parametersTree)
  {
    parametersFound = true;
    feedbackString << ">>>>> parameters found " << std::endl;
    Data_t temp_parameters;
    // read parameters and fill struct     
    readParametersFillStruc(parametersTree,temp_parameters);
    // extract relevant parameters
    int nOfParameters = parametersTree->GetEntries();
    for(int i = 0 ; i < nOfParameters ; i++)
    {
      relevant_data_t relevant;
      parametersTree->GetEntry(i);
      fillRelevant(relevant,temp_parameters);
      parameters.push_back(relevant);
    }
  }

  feedbackString << "Number of parameters found in energy deposition file = " << parameters.size() << std::endl;
  // build a quick map out of it
  for(int i = 0 ; i < parameters.size(); i++)
  {
    relevantParametersMap.insert(std::make_pair(parameters[i].ecal_position,i));
  }




  inputFile->cd("Configuration");
  TNamed SeedNameD     ("Seed"      ,((TNamed*) gDirectory->Get("Seed"))->GetTitle());
  TNamed HostNameD     ("Hostname"  ,((TNamed*) gDirectory->Get("Hostname"))->GetTitle());
  TNamed PWDNameD      ("PWD"       ,((TNamed*) gDirectory->Get("PWD"))->GetTitle());
  TNamed ConfigNameD   ("ConfigFile",((TNamed*) gDirectory->Get("ConfigFile"))->GetTitle());
  TNamed GpsNameD      ("GpsFile"   ,((TNamed*) gDirectory->Get("GpsFile"))->GetTitle());
  TNamed PrimariesNameD("Primaries" ,((TNamed*) gDirectory->Get("Primaries"))->GetTitle());
  inputFile->cd();

  // produce the module ID TH2I map
  int det_counter = 0;
  TH2I *module_ID_map = ComputeElementMap(modules,"module_ID_map",-1,-1,det_counter);
  // and the detector id map
  // TH2I *detector_ID_map_front = ComputeDetectorMap(modules,cells,0);
  // TH2I *detector_ID_map_back  = ComputeDetectorMap(modules,cells,1);



  std::vector<region_t> modules_regions = CalculateRegions(modules);
  std::vector<int> moduleTypes;
  std::vector<int> moduleSections;
  std::vector<float> moduleSeparationZ;
  std::vector<float> modulePositionZ;

  

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
      // also get sections
      moduleSections.push_back(modules_regions[iMod].sections);
      moduleSeparationZ.push_back(modules_regions[iMod].separation_z);
      modulePositionZ.push_back(modules_regions[iMod].position_z);
    }
  }

  std::vector<region_t> absorbers_regions = CalculateRegions(absorbers);
  std::vector<int>   absorbersTypes;
  std::vector<float> absorbersPositionZ;

  for(int iAbs = 0 ; iAbs < absorbers_regions.size(); iAbs++)
  {
    bool absTypeIsAlreadyThere = false;
    for(unsigned int im = 0; im < absorbersTypes.size(); im++)
    {
      if(absorbers_regions[iAbs].type == absorbersTypes[im])
      {
        absTypeIsAlreadyThere = true;
      }
    }
    if(!absTypeIsAlreadyThere)
    {
      absorbersTypes.push_back(absorbers_regions[iAbs].type);
      absorbersPositionZ.push_back(absorbers_regions[iAbs].position_z);
    }
  }


  // lots of maps
  std::map<int, int>   moduleNumberMap; // map from module type to number of modules
  std::map<int, int>   moduleNumberOfBeforeMap; // map from module type to number of other modules
  std::map<int, int>   moduleSectionsMap; // map from module type to number of segments
  std::map<int, float> moduleMinMap; // map from module type to min z
  std::map<int, float> moduleMaxMap; // map from module type to max z
  std::map<int, float> moduleSeparationMap; // map from module type to z of separatio
  std::map<int, float> modulePositionMap; // map from module type to z of separation
  
  std::map<int, int>   absNumberMap; // map from module type to number of modules
  std::map<int, float> absPositionMap; // map from type to z of abs in that module

  for(int i = 0 ; i < moduleTypes.size(); i++)
  {
    int modCounter = 0;
    for(int iMod = 0 ; iMod < modules_regions.size(); iMod++)
    {
      if(modules_regions[iMod].type == moduleTypes[i])
      {
        modCounter++;
      }
    }
    moduleNumberMap.insert(std::make_pair(moduleTypes[i],modCounter));
  }

  for(int i = 0 ; i < absorbersTypes.size(); i++)
  {
    int absCounter = 0;
    for(int iAbs = 0 ; iAbs < absorbers_regions.size(); iAbs++)
    {
      if(absorbers_regions[iAbs].type == absorbersTypes[i])
      {
        absCounter++;
      }
    }
    absNumberMap.insert(std::make_pair(absorbersTypes[i],absCounter));
    
  }
  
  for(int i = 0 ; i < absorbersTypes.size(); i++)
  {
    absPositionMap.insert(std::make_pair(absorbersTypes[i],absorbersPositionZ[i]));
  }

  // sum all modules
  // std::map<int,int>::iterator it_moduleNumberMap = moduleNumberMap.begin();

  int totModules = 0;
  for (std::map<int,int>::iterator it=moduleNumberMap.begin(); it!=moduleNumberMap.end(); ++it)
  {
    moduleNumberOfBeforeMap.insert(std::make_pair(it->first,totModules));
    totModules += it->second;
  }
  feedbackString << "Total number of modules                 = "<< totModules << std::endl;




  int NofModuleTypes    = moduleTypes.size();
  // int NofSpacalSegments = moduleTypes.size();
  feedbackString << "Number of module types                 = "<< NofModuleTypes << std::endl;

  std::map<int, int>   mapOfCellMaps; // map from showerModuleType to i-index of cell_map_per_module (then j is 0 or 1)
  // produce 2 cell maps for each module type
  TH2I *** cell_map_per_module;
  cell_map_per_module = new TH2I**[NofModuleTypes];
  for(int i = 0; i < NofModuleTypes; i++)
  {
    cell_map_per_module[i] = new TH2I*[2];
    int det_counter = 0;
    for(int j = 0; j < 2 ; j++)
    {
      std::stringstream cname;
      cname << "cell_map_";
      if(j == 0) // front
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
      if(moduleSections[i] == 1)
      {
        fb = -1;
      }
      cell_map_per_module[i][j] = ComputeElementMap(cells,cname.str(),fb,moduleTypes[i],det_counter);
      mapOfCellMaps.insert(std::make_pair(moduleTypes[i],i));
    }
  }

  // vector of found
  std::vector< std::vector<bool> > foundLimit ;

  std::vector<region_t> cells_regions = CalculateRegions(cells);


  for(int i = 0 ; i < moduleTypes.size(); i++)
  {
    // feedbackString << "Type " << moduleTypes[i] << " has " << moduleSections[i] << " sections " << std::endl;
    // feedbackString << "Type " << moduleTypes[i] << " has " << moduleNumberMap[moduleTypes[i]] << " modules " << std::endl;
    moduleSectionsMap.insert(std::make_pair(moduleTypes[i],moduleSections[i]));
    moduleSeparationMap.insert(std::make_pair(moduleTypes[i],moduleSeparationZ[i]+modulePositionZ[i]));
    modulePositionMap.insert(std::make_pair(moduleTypes[i],modulePositionZ[i]));
    // moduleAbsPositionMap.insert(std::make_pair(moduleTypes[i],modulePositionZ[i]));
    std::vector<bool> foundLimitForModule;
    foundLimitForModule.push_back(false);
    foundLimitForModule.push_back(false);
    foundLimit.push_back(foundLimitForModule);
  }

  for(int i = 0 ; i < moduleTypes.size(); i++)
  {
    for(int iCell = 0 ; iCell < cells_regions.size(); iCell++)
    {
      if(foundLimit[i][0] && foundLimit[i][1])
      {
        break;
      }
      // look for a cell of this type
      if(cells_regions[iCell].type == moduleTypes[i])
      {
        if(moduleSectionsMap.at(moduleTypes[i]) == 1)
        {
          // only 1 sections, just take max and min
          foundLimit[i][1] = true;
          moduleMaxMap.insert(std::make_pair(moduleTypes[i],cells_regions[iCell].max[2]+modulePositionZ[i]));
          foundLimit[i][0] = true;
          moduleMinMap.insert(std::make_pair(moduleTypes[i],cells_regions[iCell].min[2]+modulePositionZ[i]));
        }
        else
        {
          // find mean
          float meanCellZ = (cells_regions[iCell].min[2] + cells_regions[iCell].max[2])/2.0;
          if(meanCellZ > 0)
          {
            if(foundLimit[i][1] == false)
            {
              foundLimit[i][1] = true;
              moduleMaxMap.insert(std::make_pair(moduleTypes[i],cells_regions[iCell].max[2]+modulePositionZ[i]));
            }
          }
          else
          {
            if(foundLimit[i][0] == false)
            {
              foundLimit[i][0] = true;
              moduleMinMap.insert(std::make_pair(moduleTypes[i],cells_regions[iCell].min[2]+modulePositionZ[i]));
            }
          }
        }

      }
    }
  }

  for(int i = 0 ; i < moduleTypes.size(); i++)
  {
    feedbackString << "Type " << moduleTypes[i] << " has " << moduleSectionsMap.at(moduleTypes[i]) << " sections " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " has " << moduleNumberMap.at(moduleTypes[i]) << " modules " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " has " << moduleNumberOfBeforeMap.at(moduleTypes[i]) << " modules before" << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " min " << moduleMinMap.at(moduleTypes[i]) << " mm " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " max " << moduleMaxMap.at(moduleTypes[i]) << " mm " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " separation at " << moduleSeparationMap.at(moduleTypes[i]) << " mm " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " module position at " << modulePositionMap.at(moduleTypes[i]) << " mm " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " absorber position at " << absPositionMap.at(moduleTypes[i]) << " mm " << std::endl;
  }


  feedbackString << std::endl;
  std::cout << feedbackString.str() << std::endl;




  //----------------------------------//
  // READING CALIBRATIONS             //
  //----------------------------------//
  // tokenize and read
  // calibration string
  std::stringstream cali_stream(calibrationFileNameString); //create string stream from the string
  while(cali_stream.good()) {
      std::string substr;
      getline(cali_stream, substr, ','); //get first string delimited by comma
      trim(substr);
      calibrationFileName.push_back(substr);
  }

  // types string
  if(typeString != "")
  {
    std::stringstream type_stream(typeString); //create string stream from the string
    while(type_stream.good()) {
        std::string substr;
        getline(type_stream, substr, ','); //get first string delimited by comma
        trim(substr);
        int t_number;
        std::istringstream ( substr ) >> t_number;
        typeNumber.push_back(t_number);
    }
  }



  // if calibrationFileName is only 1, and typeNumber is empty, set by default the only typeNumber to 0
  // it means only one module type has been simulated, and if not type is specified then the user didn't
  // use the CAD module and the simulation would give 0 to module type
  if(calibrationFileName.size() == 1)
  {
    if(typeNumber.size() == 0)
    {
      typeNumber.push_back(0);
    }
    if(typeNumber.size() > 1)
    {
      std::cout << "ERROR: You need to provide an equal number of arguments for --calibration and --type! See usage..." << std::endl ;
      usage();
      return -1;
    }
  }
  else
  {
    if(calibrationFileName.size() != typeNumber.size())
    {
      std::cout << "ERROR: You need to provide an equal number of arguments for --calibration and --type! See usage..." << std::endl ;
      usage();
      return -1;
    }
  }

  for(int i = 0; i<calibrationFileName.size(); i++) {    //print all splitted strings
      std::cout << calibrationFileName.at(i) << std::endl;
  }
  for(int i = 0; i<typeNumber.size(); i++) {    //print all splitted strings
      std::cout << typeNumber.at(i) << std::endl;
  }

  // int module_types_final[2] = {5,6};
  // // hardcoded for now
  // std::vector<std::string> calibrationFileNameTemp;
  // calibrationFileNameTemp.push_back("calibrationYAG.data");
  // calibrationFileNameTemp.push_back("calibrationGAGG.data");

  // map from module type to calibrations std::vector index
  std::map<int, int> calibrationMap;
  // fill calibration structures
  std::vector<calibration_t> calibrations;

  // need to find a way not to repeat import and vector filling 
  // of the same calibration file. strategy 
  // 1. import calibration files only once 
  // 2. write the correspondance map 
  //
  // 1. Run on calibration file names, remove duplicates 
  // std::vector<std::string>
  std::set<std::string> s;
  unsigned sizeCali = calibrationFileName.size();
  for( unsigned i = 0; i < sizeCali; ++i ) s.insert( calibrationFileName[i] );
  std::vector<std::string> uniqueCali;
  uniqueCali.assign(s.begin(), s.end());
  // std::cout << "------------------------" << std::endl;
  // for(int i = 0 ; i < calibrationFileName.size(); i++)
  // {
  //   std::cout << calibrationFileName[i] << std::endl;
  // }
  // std::cout << "------------------------" << std::endl;
  // std::cout << "------------------------" << std::endl;
  // for(int i = 0 ; i < uniqueCali.size(); i++)
  // {
  //   std::cout << uniqueCali[i] << std::endl;
  // }
  // std::cout << "------------------------" << std::endl;

  // now import only these files 

  for(int i = 0 ; i < uniqueCali.size(); i++)
  {
    calibration_t calibration_temp;
    // calibration_temp.type = typeNumber[i];
    calibration_temp.type = i; // temp
    calibration_temp.footer = false; // temp
    importCalibrationFile(uniqueCali[i],calibration_temp);
    

    // calibrationMap.insert(std::make_pair(typeNumber[i],calibrations.size()));
    calibrations.push_back(calibration_temp);
  }

  // check that calibrations and typeNumber have the same size 
  if(calibrations.size() != NofModuleTypes)
  {
    std::cout << ">>>>>>>>>>>>>>>>>>>> WARNING! Incompatible number of input module types and calibrations found!" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>          1) Part of the energy deposition data could be lost in the hybrid files" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>          2) Also, this might cause a crash, unless you specified a proper --type when calling this program" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>          Energy deposition files has " << NofModuleTypes << " types" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>>>          Unique calibrations are     " << calibrations.size() << std::endl;
    std::cout << std::endl;
  }

  // now run on all calibration file names 
  for(int i = 0 ; i < calibrationFileName.size(); i++)
  {
    // for each, run on calibrations imported 
    for(int j = 0 ; j < calibrations.size(); j++)
    {
      // compare filenames 
      // std::cout << calibrationFileName[i] << " " << calibrations[j].fileName << std::endl;
      if(calibrationFileName[i].compare(calibrations[j].fileName) == 0)
      {
        calibrations[j].type = typeNumber[i];
        // get and calc shift 
        float caliFile_moduleCenter = calibrations[j].moduleCenter;
        // out module center 
        float showerShift = 0.;
        if(calibrations[j].footer) // calculate shift only if the calibration file has this info
        {
          float outFile_moduleCenter = modulePositionMap.at(calibrations[j].type) + absPositionMap.at(calibrations[j].type);
          showerShift = caliFile_moduleCenter - outFile_moduleCenter;
          calibrations[j].showerShift = showerShift;
          std::cout << "Type " << calibrations[j].type << " "
                    << "Cali file center " << caliFile_moduleCenter << " mm, " 
                    << "Out file center " << outFile_moduleCenter  << " mm, "
                    << "Shower shift " << showerShift  << " mm "
                    << std::endl;
        }
        else 
        {
          std::cout << "Old type calibration files, no (possible) shifts between calibration and module will be considered! " << std::endl;
        }
        calibrationMap.insert(std::make_pair(typeNumber[i],j));
        break;
      }
    }
  }

  // return 0;


  //----------------------------------//
  // COMPARE CALIBRATION TO EXPECTED  //
  //----------------------------------//

  // we have 
  // calibrationMap        -> map from type number to calibration struct (which includes the relevant struct)
  // relevantParametersMap -> map from type number to parameters map, i.e. the relevant struct 
  // we can loop on them and compare 
  // remember that when there is only 1 config file, the type number is set to 0 by the geant simulation 
  // (because only main config file is there, check method void Parameters::ReadMainConfig(std::string fileName) 
  // in Parameters.cc).
  // so it would be 0 in optical cali relevant struct, but in propagate the user specifies a type to assign 
  // to an optical calibration, and that's in the map first entry. remember that 
  // if calibrationFileName is only 1, and typeNumber is empty, set by default the only type number to 0
  std::cout << std::endl;
  // attempt comparison only if en depo file has parameters
  if(parametersFound)
  {
    // and do only if types > 1
    // if(typeNumber.size() > 1) // is this necessary?
  
    std::cout << "/////////////////////////////////////////////////////////////////" << std::endl;
    std::cout << "//COMPARISON OF MODULE AND CALIBRATION STRUCTURES              //" << std::endl;
    std::cout << "/////////////////////////////////////////////////////////////////" << std::endl;
    std::cout << std::endl;
    for(int i = 0; i<typeNumber.size(); i++) // loop on types 
    {
      std::cout << "/////////////////////////////////////////////////////////////////" << std::endl;
      std::cout << "// Type " <<  typeNumber[i] << " compatibility test" << std::endl;
      // find index of calibration vector 
      std::map<int, int> ::const_iterator iterCal = calibrationMap.find(typeNumber[i]);
      int calIndex ;
      if(iterCal != calibrationMap.end())
      {
        // iterCal is item pair in the map. The value will be accessible as `iter->second`.
        calIndex = iterCal->second;
      }
      else
      {
        std::cout << "// WARNING: No calibration file found for type " << typeNumber[i] << std::endl;
        std::cout << "// propagateHybrid will continue, we will assume you chose the correct calibration file." << std::endl;
        continue;
      }

      // go on only if chosen calibration has footer 
      if(calibrations[calIndex].footer)
      {
        // find index of relevant vector
        std::map<int, int> ::const_iterator iterRel = relevantParametersMap.find(typeNumber[i]);
        int relIndex ;
        if(iterRel != relevantParametersMap.end())
        {
          // iter is item pair in the map. The value will be accessible as `iter->second`.
          relIndex = iterRel->second;
        }
        else
        {
          std::cout << "// WARNING: No parameters in energy deposition file found for type " << typeNumber[i] << std::endl;
          std::cout << "// propagateHybrid will continue, we will assume you chose the correct calibration file." << std::endl;
          continue;
        }

        // and finally compare...
        int equals = compareRelevant(parameters[relIndex],calibrations[calIndex].parameters);
        if(equals == 0)
        {
          std::cout << "// Type " <<  typeNumber[i] << " test PASSED!" << std::endl;
        }
        else 
        {
          std::cout << ">>>>>>>>>>>>>>>>>>>> ERROR! Incompatible calibration(s) found!" << std::endl;
          std::cout << ">>>>>>>>>>>>>>>>>>>> Aborting..." << std::endl;
          // force root to quit..
          gROOT->ProcessLine(".qqqqqqqqqqqqqqqqqqqq");
          return -1;
        }
      }
      else 
      {
        std::cout << "// Test not possible, chosen calibration file has no footer info (old calibration version)!" << std::endl;
        std::cout << "// propagateHybrid will continue, we will assume you chose the correct calibration file." << std::endl;
      }
    }
    std::cout << "/////////////////////////////////////////////////////////////////" << std::endl;
    std::cout << std::endl;
  }
  else 
  {
    std::cout << "/////////////////////////////////////////////////////////////////" << std::endl;
    std::cout << "// Compatibility test not possible, energy deposition file has no parameters TTree (old out version)" << std::endl;
    std::cout << "// propagateHybrid will continue, we will assume you chose the correct calibration file." << std::endl;
    std::cout << "/////////////////////////////////////////////////////////////////" << std::endl;
  }
  std::cout << std::endl;
  
  


  //----------------------------------//
  // OUTPUT TTREE                     //
  //----------------------------------//
  // variables
  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");


  TTree* outTreeHybrid = new TTree ("hybrid", "hybrid") ;
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
  std::string hb_processName;
  Int_t       hb_processNumber;
  float      hb_primaryPositionOnAbsorberX    ;
  float      hb_primaryPositionOnAbsorberY    ;
  float      hb_primaryPositionOnAbsorberZ    ;
  float      hb_primaryMomentumOnAbsorberX    ;
  float      hb_primaryMomentumOnAbsorberY    ;
  float      hb_primaryMomentumOnAbsorberZ    ;
  float      hb_primaryEnergyOnAbsorber       ;

  outTreeHybrid->Branch("event"    ,   &hb_event        , "event/I");
  outTreeHybrid->Branch("x",   &hb_x    , "x/F");
  outTreeHybrid->Branch("y",   &hb_y    , "y/F");
  outTreeHybrid->Branch("z",   &hb_z    , "z/F");
  outTreeHybrid->Branch("detector" ,   &hb_detector     , "detector/I");
  outTreeHybrid->Branch("module" ,   &hb_module     , "module/I");
  outTreeHybrid->Branch("moduleType" ,   &hb_moduleType     , "moduleType/I");
  outTreeHybrid->Branch("timestamp",   &hb_timestamp    , "timestamp/F");
  outTreeHybrid->Branch("t_deposition",     &hb_t_deposition     , "t_deposition/F");
  outTreeHybrid->Branch("t_generation",     &hb_t_generation     , "t_generation/F");
  outTreeHybrid->Branch("t_propagation",    &hb_t_propagation    , "t_propagation/F");
  outTreeHybrid->Branch("PhotonEnergy",&hb_photon_energy, "PhotonEnergy/F");
  outTreeHybrid->Branch("processName" ,&hb_processName);
  outTreeHybrid->Branch("processNumber" ,   &hb_processNumber        , "processNumber/I");
  outTreeHybrid->Branch("primary_PositionOnAbsorberX" ,&hb_primaryPositionOnAbsorberX    ,"primaryPositionOnAbsorberX/F");
  outTreeHybrid->Branch("primary_PositionOnAbsorberY" ,&hb_primaryPositionOnAbsorberY    ,"primaryPositionOnAbsorberY/F");
  outTreeHybrid->Branch("primary_PositionOnAbsorberZ" ,&hb_primaryPositionOnAbsorberZ    ,"primaryPositionOnAbsorberZ/F");
  outTreeHybrid->Branch("primary_MomentumOnAbsorberX" ,&hb_primaryMomentumOnAbsorberX    ,"primaryMomentumOnAbsorberX/F");
  outTreeHybrid->Branch("primary_MomentumOnAbsorberY" ,&hb_primaryMomentumOnAbsorberY    ,"primaryMomentumOnAbsorberY/F");
  outTreeHybrid->Branch("primary_MomentumOnAbsorberZ" ,&hb_primaryMomentumOnAbsorberZ    ,"primaryMomentumOnAbsorberZ/F");
  outTreeHybrid->Branch("primary_EnergyOnAbsorber" ,&hb_primaryEnergyOnAbsorber       ,"primaryEnergyOnAbsorber/F");

  //----------------------------------//
  // OUTPUT SHOWER                    //
  //----------------------------------//
  TTree*           out_shower = new TTree ("shower" , "shower") ;
  Int_t            out_showerRun;
  Int_t            out_showerEvent;
  Int_t            out_showerIsInCrystal;
  Int_t            out_shower_moduleType;
  Int_t            out_pdgID;
  Int_t            out_primaryID;
  Int_t            out_primaryPDGID;
  Float_t          out_primaryEnergy;
  Float_t          out_showerX    ;
  Float_t          out_showerY    ;
  Float_t          out_showerZ    ;
  Float_t          out_showerT    ;
  Float_t          out_showerTotalEnDep ;
  Float_t          out_showerIonizingEnDep;
  Float_t          out_showerNonIonizingEnDep;
  std::string      out_showerProcessName;
  std::string      out_showerMaterialName;
  Int_t            out_showerMaterialNumber;
  Float_t          out_showerLocalX    ;
  Float_t          out_showerLocalY    ;
  Float_t          out_showerLocalZ    ;
  Float_t          out_shower_primaryPositionOnAbsorberX    ;
  Float_t          out_shower_primaryPositionOnAbsorberY    ;
  Float_t          out_shower_primaryPositionOnAbsorberZ    ;
  Float_t          out_shower_primaryMomentumOnAbsorberX    ;
  Float_t          out_shower_primaryMomentumOnAbsorberY    ;
  Float_t          out_shower_primaryMomentumOnAbsorberZ    ;
  Float_t          out_shower_primaryEnergyOnAbsorber       ;
  out_shower->Branch("run"                        , &out_showerRun                         , "run/I");
  out_shower->Branch("event"                      , &out_showerEvent                       , "event/I");
  out_shower->Branch("isInCrystal"                , &out_showerIsInCrystal                 , "isInCrystal/I");
  out_shower->Branch("moduleType"                 , &out_shower_moduleType                 , "moduleType/I");
  out_shower->Branch("pdgID"                      , &out_pdgID                             , "pdgID/I");
  out_shower->Branch("primaryID"                  , &out_primaryID                         , "primaryID/I");
  out_shower->Branch("primaryPDGID"               , &out_primaryPDGID                      , "primaryPDGID/I");
  out_shower->Branch("primaryEnergy"              , &out_primaryEnergy                     , "primaryEnergy/F");
  out_shower->Branch("x"                          , &out_showerX                           , "x/F");
  out_shower->Branch("y"                          , &out_showerY                           , "y/F");
  out_shower->Branch("z"                          , &out_showerZ                           , "z/F");
  out_shower->Branch("t"                          , &out_showerT                           , "t/F");
  out_shower->Branch("totalEnDep"                 , &out_showerTotalEnDep                  , "totalEnDep/F");
  out_shower->Branch("ionizingEnDep"              , &out_showerIonizingEnDep               , "ionizingEnDep/F");
  out_shower->Branch("nonIonizingEnDep"           , &out_showerNonIonizingEnDep            , "nonIonizingEnDep/F");
  out_shower->Branch("processName"                , &out_showerProcessName);
  out_shower->Branch("materialName"               , &out_showerMaterialName);
  out_shower->Branch("materialNumber"             , &out_showerMaterialNumber              , "materialNumber/I");
  out_shower->Branch("localX"                     , &out_showerLocalX                      , "localX/F");
  out_shower->Branch("localY"                     , &out_showerLocalY                      , "localY/F");
  out_shower->Branch("localZ"                     , &out_showerLocalZ                      , "localZ/F");
  out_shower->Branch("primary_PositionOnAbsorberX", &out_shower_primaryPositionOnAbsorberX , "primary_PositionOnAbsorberX/F");
  out_shower->Branch("primary_PositionOnAbsorberY", &out_shower_primaryPositionOnAbsorberY , "primary_PositionOnAbsorberY/F");
  out_shower->Branch("primary_PositionOnAbsorberZ", &out_shower_primaryPositionOnAbsorberZ , "primary_PositionOnAbsorberZ/F");
  out_shower->Branch("primary_MomentumOnAbsorberX", &out_shower_primaryMomentumOnAbsorberX , "primary_MomentumOnAbsorberX/F");
  out_shower->Branch("primary_MomentumOnAbsorberY", &out_shower_primaryMomentumOnAbsorberY , "primary_MomentumOnAbsorberY/F");
  out_shower->Branch("primary_MomentumOnAbsorberZ", &out_shower_primaryMomentumOnAbsorberZ , "primary_MomentumOnAbsorberZ/F");
  out_shower->Branch("primary_EnergyOnAbsorber"   , &out_shower_primaryEnergyOnAbsorber    , "primary_EnergyOnAbsorber/F");

  //----------------------------------//
  // OUTPUT STRUCTURES                //
  //----------------------------------//

  TTree* out_fibres = new TTree ("fibres", "fibres") ;
  TTree* out_cells = new TTree ("cells", "cells") ;
  TTree* out_holes = new TTree ("holes", "holes") ;
  TTree* out_modules = new TTree ("modules", "modules") ;
  TTree* out_absorbers = new TTree ("absorbers", "absorbers") ;
  CopyStructure(fibres,out_fibres);
  CopyStructure(cells,out_cells);
  CopyStructure(holes,out_holes);
  CopyStructure(absorbers,out_absorbers);
  CopyStructure(modules,out_modules);

  long long int counter = 0;

  //----------------------------------//
  // TIME AND ENERGY HISTOS           //
  //----------------------------------//

  // two cases :
  // 1. use material characteristics from out file
  // 2. import material from external file
  // it should all be down to pointing to the correct file...
  TFile* externalFile = NULL;
  int externalMaterialNumber = -1;
  if(useExternal)
  {
    externalFile = TFile::Open(externalFileName.c_str());
    externalFile->cd();
  }
  else
  {
    inputFile->cd();
  }

  std::map<int, TH1F*> eHistoMap;
  std::map<int, TH1F*> tHistoMap;
  // get relevant energy and time histos
  TList *listKeys = gDirectory->GetListOfKeys();
  int nKeys = listKeys->GetEntries();
  std::vector<std::string> keysName;
  // fill a vector with the leaves names
  std::string e_prefix("enHisto_");
  std::string t_prefix("tHisto_");
  for(int i = 0 ; i < nKeys ; i++){
    keysName.push_back(listKeys->At(i)->GetName());
  }
  // get e histo and t histo
  for(unsigned int i = 0 ; i < keysName.size() ; i++)
  {
    if (!keysName[i].compare(0, e_prefix.size(), e_prefix))
    {
      // en histo
      std::size_t found      = keysName[i].find_last_of(e_prefix);
      std::string mat_string = keysName[i].substr(found  +1);
      int crystalNumber      = atoi( mat_string.c_str() );
      if(useExternal)
      {
        // in this special case, only 1 material should be in the external file (if not, the user is crazy)
        // so we save the number in order to be able to overwrite the choice later
        externalMaterialNumber = crystalNumber;
      }
      TH1F* energyHistogram  = (TH1F*) gDirectory->Get(keysName[i].c_str());
      eHistoMap.insert(std::make_pair(crystalNumber,energyHistogram));
    }
    if (!keysName[i].compare(0, t_prefix.size(), t_prefix))
    {
      // en histo
      std::size_t found      = keysName[i].find_last_of(t_prefix);
      std::string mat_string = keysName[i].substr(found  +1);
      int crystalNumber      = atoi( mat_string.c_str() );
      TH1F* timeHistogram  = (TH1F*) gDirectory->Get(keysName[i].c_str());
      tHistoMap.insert(std::make_pair(crystalNumber,timeHistogram));
    }
  }

  // get crystal numbers, light yields and res scales
  std::vector<int>   *pCrystalNumber;
  std::vector<float> *pLightYield;
  std::vector<float> *pResolutionScale;
  gDirectory->GetObject("crystalMaterialList"   ,pCrystalNumber);
  gDirectory->GetObject("crystalLightYield"     ,pLightYield);
  gDirectory->GetObject("crystalResolutionScale",pResolutionScale);

  std::vector<int>   CrystalNumberList   = pCrystalNumber[0];
  std::vector<float> LightYieldList      = pLightYield[0];
  std::vector<float> ResolutionScaleList = pResolutionScale[0];

  // and put in maps..
  std::map<int, float> LightYieldMap;
  std::map<int, float> ResolutionScaleMap;
  for(int i = 0; i < CrystalNumberList.size();i++)
  {
    LightYieldMap.insert(std::make_pair(CrystalNumberList[i],LightYieldList[i]));
    ResolutionScaleMap.insert(std::make_pair(CrystalNumberList[i],ResolutionScaleList[i]));
  }

  for(int i = 0; i < CrystalNumberList.size();i++)
  {
    std::cout << CrystalNumberList[i] << " "
              << LightYieldList[i] << " "
              << ResolutionScaleList[i] << " "
              << std::endl;
  }

  // in any case, go back to input cd
  inputFile->cd();

  bool foundInAll = false;
  bool foundLY = false;
  bool foundRes = false;
  bool foundEHisto = false;
  bool foundTHisto = false;



  //----------------------------------//
  // INPUT SHOWER                     //
  //----------------------------------//
  int              showerRun              ;
  int              showerEvent            ;
  int              showerIsInCrystal      ;
  // int              showerCrystalID        ;
  // int              showerCellID           ;
  // int              showerModuleID         ;
  int              showerModuleType       ;
  // float            shower_module_x        ;
  // float            shower_module_y        ;
  // float            shower_module_z        ;
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
  // TBranch *b_showerCrystalID;
  // TBranch *b_showerCellID           ;
  // TBranch *b_showerModuleID           ;
  TBranch *b_showerModuleType       ;
  // TBranch *b_shower_module_x        ;
  // TBranch *b_shower_module_y        ;
  // TBranch *b_shower_module_z        ;
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
  // shower->SetBranchAddress("crystalID"          ,&showerCrystalID       ,&b_showerCrystalID       );
  // shower->SetBranchAddress("cellID"             ,&showerCellID          ,&b_showerCellID          );
  // shower->SetBranchAddress("moduleID"           ,&showerModuleID        ,&b_showerModuleID        );
  shower->SetBranchAddress("moduleType"         ,&showerModuleType      ,&b_showerModuleType      );
  // shower->SetBranchAddress("module_x"           ,&shower_module_x       ,&b_shower_module_x       );
  // shower->SetBranchAddress("module_y"           ,&shower_module_y       ,&b_shower_module_y       );
  // shower->SetBranchAddress("module_z"           ,&shower_module_z       ,&b_shower_module_z       );
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

  int feedbackDivision;
  long int nEntriesShower = shower->GetEntries();
  long long int lostScintillation = 0;
  long long int lostShower = 0;
  // long int nEntriesOptPhotGen = opticalPhotonGen->GetEntries();
  // std::cout << nEntriesOptPhotGen <<  std::endl;
  std::cout << "Energy deposition events = " << nEntriesShower <<  std::endl;
  std::cout << "Running on energy depositions..." << std::endl;
  long long int totalNumberOfOpticalPhotons = 0;
  long long int totalNumberOfOpticalPhotonsOut = 0;

  counter = 0;
  feedbackDivision = (int) (nEntriesShower/20);
  if(feedbackDivision == 0) feedbackDivision = 1;
  for(int iEntry = 0 ; iEntry < nEntriesShower ; iEntry++)
  // for(int iEntry = 5974 ; iEntry < 5975 ; iEntry++)
  {
    if (  iEntry % feedbackDivision == 0 ) { std::cout << (int) ((double) iEntry/nEntriesShower * 100) << "% done..." << std::endl;}
    shower->GetEvent(iEntry);

    // copy the shower data to output shower

    out_showerRun                         = showerRun                        ;
    out_showerEvent                       = showerEvent                      ;
    out_showerIsInCrystal                 = showerIsInCrystal                ;
    out_shower_moduleType                 = showerModuleType                 ;
    out_pdgID                             = shower_pdgID                     ;
    out_primaryID                         = shower_primaryID                 ;
    out_primaryPDGID                      = shower_primaryPDGID              ;
    out_primaryEnergy                     = shower_primaryEnergy             ;
    out_showerX                           = showerX                          ;
    out_showerY                           = showerY                          ;
    out_showerZ                           = showerZ                          ;
    out_showerT                           = showerT                          ;
    out_showerTotalEnDep                  = showerTotalEnDep                 ;
    out_showerIonizingEnDep               = showerIonizingEnDep              ;
    out_showerNonIonizingEnDep            = showerNonIonizingEnDep           ;
    out_showerProcessName                 = showerProcessName->c_str()       ;
    out_showerMaterialName                = showerMaterialName->c_str()      ;
    out_showerMaterialNumber              = showerMaterialNumber             ;
    out_showerLocalX                      = shower_local_X                   ;
    out_showerLocalY                      = shower_local_Y                   ;
    out_showerLocalZ                      = shower_local_Z                   ;
    out_shower_primaryPositionOnAbsorberX = shower_primaryPositionOnAbsorberX;
    out_shower_primaryPositionOnAbsorberY = shower_primaryPositionOnAbsorberY;
    out_shower_primaryPositionOnAbsorberZ = shower_primaryPositionOnAbsorberZ;
    out_shower_primaryMomentumOnAbsorberX = shower_primaryMomentumOnAbsorberX;
    out_shower_primaryMomentumOnAbsorberY = shower_primaryMomentumOnAbsorberY;
    out_shower_primaryMomentumOnAbsorberZ = shower_primaryMomentumOnAbsorberZ;
    out_shower_primaryEnergyOnAbsorber    = shower_primaryEnergyOnAbsorber   ;
    out_shower->Fill();

    // just overwrite the showerMaterialNumber in the special case of external material
    if(useExternal)
    {
      showerMaterialNumber = externalMaterialNumber;
    }




    // now find closest point on grid or interpolate between points ...
    // brutally closest x y z e
    // this means of course the same approach




    // for now xyz ok like this, but in reality at xy will change for diff crystal, if only optical scan on axis ignore, otherwise need to recenter...
    if(showerIsInCrystal)
    {


      // fetch index of calibrations vector
      // int calIndex = calibrationMap.find(showerModuleType);
      std::map<int, int> ::const_iterator iter = calibrationMap.find(showerModuleType);
      int calIndex ;
      if(iter != calibrationMap.end())
      {
        // iter is item pair in the map. The value will be accessible as `iter->second`.
        calIndex = iter->second;
      }
      else
      {
        if(verbose > 0)
        {
          std::cout << "WARNING: Calibration not found for shower event located in moduleType  = " << showerModuleType << std::endl;
          std::cout << "                                                        at position    = " << showerX  << " " << showerY << " " << showerZ << std::endl;
          std::cout << "         This energy deposition event will be ignored! Check the input calibration files provided." << std::endl;
          std::cout << std::endl;
        }
        lostShower++;
        continue;
      }


      int xIndex = getPositionOfLevel(calibrations[calIndex].x,shower_local_X);
      int yIndex = getPositionOfLevel(calibrations[calIndex].y,shower_local_Y);
      // shower z needs to be shifted by the (possible) difference in module centers, in world z coordinates
      int zIndex = getPositionOfLevel(calibrations[calIndex].z,showerZ+calibrations[calIndex].showerShift);
      // int zIndex = getPositionOfLevel(calibrations[calIndex].z,showerZ);
      // std::cout << showerZ << " " << zIndex << std::endl;
      // take  and x and y by local x and y, so it's always referred to the axis of the crystal and you can use just one crystal in optical calibration
      // int xIndex = (int) roundf((shower_local_X - calibrations[calIndex].xmin)/calibrations[calIndex].xstep);
      // int yIndex = (int) roundf((shower_local_Y - calibrations[calIndex].ymin)/calibrations[calIndex].ystep);
      // take z from global z
      // int zIndex = (int) roundf((showerZ - calibrations[calIndex].zmin)/calibrations[calIndex].zstep); // this is of course ok, but remember that the optical calibration needs to be consistent in z
      // i.e. the optical calibration needs to be performed with the same calorimeter_position key value, and you need to properly prepare the
      // optical calibration run, setting the zmin zmax values accordingly

      // HACK to avoid borders when x y z are exactly on boundaries
      // if(xIndex < 0) xIndex = 0;
      // if(yIndex < 0) yIndex = 0;
      // if(zIndex < 0) zIndex = 0;
      // if(xIndex >= calibrations[calIndex].xn) xIndex = calibrations[calIndex].xn -1;
      // if(yIndex >= calibrations[calIndex].yn) yIndex = calibrations[calIndex].yn -1;
      // if(zIndex >= calibrations[calIndex].zn) zIndex = calibrations[calIndex].zn -1;


      // int cZ_Index = getPositionOfLevel(calibrations[calIndex].z,showerZ);
      //
      // if(cZ_Index != zIndex)
      // {
      //   std::cout << showerZ << " " << zIndex << " " << cZ_Index << std::endl;
      // }



      // then need to decide how many photons to produce
      // look into g4 source code (see end of code *-*-*)

      // calc photons produced
      float MeanNumberOfPhotons;
      std::map<int, float> ::const_iterator iterLYmap = LightYieldMap.find(showerMaterialNumber);
      if(iterLYmap != LightYieldMap.end())
      {
        // iter is item pair in the map. The value will be accessible as `iter->second`.
        MeanNumberOfPhotons = (iterLYmap->second)*showerTotalEnDep;
      }
      else
      {
        if(verbose > 0)
        {
          std::cout << "WARNING: Light yield map not found for showerMaterialNumber " << showerMaterialNumber << std::endl;
          std::cout << "         This event will be ignored! Check the input files provided." << std::endl;
          std::cout << std::endl;
        }
        lostShower++;
        continue;
      }
       // = LightYieldMap.at(showerMaterialNumber)*showerTotalEnDep;
      int NumPhotons = 0;
      if (MeanNumberOfPhotons > 10.)
      {

        float sigma;
        std::map<int, float> ::const_iterator iterRSmap = ResolutionScaleMap.find(showerMaterialNumber);
        if(iterRSmap != ResolutionScaleMap.end())
        {
          // iter is item pair in the map. The value will be accessible as `iter->second`.
          sigma = (iterRSmap->second) * std::sqrt(MeanNumberOfPhotons);
        }
        else
        {
          if(verbose > 0)
          {
            std::cout << "WARNING: ResolutionScale map not found for showerMaterialNumber " << showerMaterialNumber << std::endl;
            std::cout << "         This event will be ignored! Check the input files provided." << std::endl;
            std::cout << std::endl;
          }
          lostShower++;
          continue;
        }
        NumPhotons = (int) (rand->Gaus(MeanNumberOfPhotons,sigma) + 0.5);
      }
      else
      {
        NumPhotons = (int) (rand->Poisson(MeanNumberOfPhotons));
      }
      totalNumberOfOpticalPhotons+=NumPhotons;

      // fetch histograms here
      TH1F *enHisto;
      std::map<int, TH1F*> ::const_iterator iterEnHistomap = eHistoMap.find(showerMaterialNumber);
      if(iterEnHistomap != eHistoMap.end())
      {
        // iter is item pair in the map. The value will be accessible as `iter->second`.
        enHisto = (iterEnHistomap->second);
      }
      else
      {
        if(verbose > 0)
        {
          std::cout << "WARNING: Energy Histogram map not found event for showerMaterialNumber " << showerMaterialNumber << std::endl;
          std::cout << "         This event will be ignored! Check the input files provided." << std::endl;
          std::cout << std::endl;
        }
        lostShower++;
        continue;
      }


      TH1F *tHisto;
      std::map<int, TH1F*> ::const_iterator iterThistomap = tHistoMap.find(showerMaterialNumber);
      if(iterThistomap != tHistoMap.end())
      {
        // iter is item pair in the map. The value will be accessible as `iter->second`.
        tHisto = (iterThistomap->second);
      }
      else
      {
        if(verbose > 0)
        {
          std::cout << "WARNING: Time Histogram map not found event for showerMaterialNumber " << showerMaterialNumber << std::endl;
          std::cout << "         This event will be ignored! Check the input files provided." << std::endl;
          std::cout << std::endl;
        }
        lostShower++;
        continue;
      }



      // TH1F *enHisto = eHistoMap.at(showerMaterialNumber);
      // TH1F *tHisto  = tHistoMap.at(showerMaterialNumber);

      if(verbose > 0)
      {
        std::cout << "---------------------------------------------------------" << std::endl;
        std::cout << "|                  SHOWER ENTRY N "              << iEntry << std::endl;
        std::cout << "---------------------------------------------------------" << std::endl;
        std::cout << "showerIsInCrystal " << showerIsInCrystal << std::endl;
        std::cout << "showerModuleType " << showerModuleType << std::endl;
        std::cout << "calIndex         " << calIndex << std::endl;
        std::cout << "calibrations[calIndex].type     " << calibrations[calIndex].type << std::endl;
        std::cout << "calibrations[calIndex].fileName " << calibrations[calIndex].fileName << std::endl;
        std::cout << "shower_local_X xIndex           " << shower_local_X << " " << xIndex << std::endl;
        std::cout << "shower_local_Y yIndex           " << shower_local_Y << " " << yIndex << std::endl;
        std::cout << "showerZ        zIndex           " << showerZ        << " " << zIndex << std::endl;
        std::cout << "MeanNumberOfPhotons " << MeanNumberOfPhotons << std::endl;
        std::cout << "NumPhotons " << NumPhotons << std::endl;
        std::cout << "enHisto       " << enHisto->GetName() << std::endl;
        std::cout << "tHisto        " << tHisto->GetName() << std::endl;
      }

      for(int iPhot = 0; iPhot < NumPhotons; iPhot++)// run on all produced
      {
        // assign a wavelenght
        // use same as primary gen
        float photon_energy = enHisto->GetRandom(); // in eV

        // find closest calibration point index
        int eIndex = getPositionOfLevel(calibrations[calIndex].energy,photon_energy);
        // int eIndex = (int) roundf((photon_energy - calibrations[calIndex].emin)/calibrations[calIndex].estep);
        // HACK to avoid borders when e is exactly on boundaries or far beyond
        // if(eIndex < 0) eIndex = 0;
        // if(eIndex >= calibrations[calIndex].en) eIndex = calibrations[calIndex].en -1;

        // roll the dice for detection
        double roll_the_dice =  rand->Uniform(1.0);
        // assign to a readout accordingly, or discard
        int sIndex = -1;
        float totProb = 0;
        bool foundIndex = false;

        if(verbose > 1)
        {
          std::cout << "-----------------------------" << std::endl;
          std::cout << "|   Photon " << iPhot << std::endl;
          std::cout << "-----------------------------" << std::endl;

          std::cout << "photon_energy " << photon_energy << std::endl;
          std::cout << "calibrations[calIndex].emin  " << calibrations[calIndex].emin  << std::endl;
          std::cout << "calibrations[calIndex].estep " << calibrations[calIndex].estep << std::endl;
          std::cout << "calibrations[calIndex].en    " << calibrations[calIndex].en << std::endl;
          std::cout << "eIndex        " << eIndex << std::endl;
          std::cout << "roll_the_dice        " << roll_the_dice << std::endl;
        }

        // look for side
        for(int iSide = 0 ; iSide < calibrations[calIndex].nSides; iSide++)
        {
          if(foundIndex)
          {
            break;
          }

          totProb += calibrations[calIndex].calibration_point[xIndex][yIndex][zIndex][eIndex][iSide].eff;


          if(roll_the_dice < totProb)
          {
            sIndex = iSide;
            foundIndex = true;
          }
          if(verbose > 1)
          {
            std::cout << "totProb iSide " << iSide << " = "<< totProb << std::endl;
            std::cout << "foundIndex " << foundIndex << std::endl;
          }
        }
        if(verbose > 1)
        {
          std::cout << "sIndex " << sIndex << std::endl;
        }


        // record or discard
        if(sIndex == -1)
        {
          // do nothing, discard
        }
        else
        {
          // find cell
          // int cell_number = showerCellID;
          // int cell_number = 0;

          // if(cell_number == -1 )
          // {
          //   std::cout << "WARNING: cell not found for Scintillation event located at position " << showerX  << " " << showerY << " " << showerZ << std::endl;
          //   std::cout << std::endl;
          //   lostScintillation++;
          // }



          // // int module_index = module_index_x + module_index_y*ModulesN[1];
          // // int cell_global_number = cell_number + module_index*cells_regions.size();
          // int module_number_before;
          // std::map<int, int> ::const_iterator iterModNumbBefore = moduleNumberOfBeforeMap.find(showerModuleType);
          // if(iterModNumbBefore != moduleNumberOfBeforeMap.end())
          // {
          //   // iter is item pair in the map. The value will be accessible as `iter->second`.
          //   module_number_before = (iterModNumbBefore->second);
          // }
          // else
          // {
          //   if(verbose > 0)
          //   {
          //     std::cout << "WARNING: Module number before not found event for showerModuleType " << showerModuleType << std::endl;
          //     std::cout << "         This event will be ignored! Check the input files provided." << std::endl;
          //     std::cout << std::endl;
          //   }
          //   lostScintillation++;
          //   continue;
          // }
          Int_t binx = module_ID_map->GetXaxis()->FindBin(showerX);
          Int_t biny = module_ID_map->GetYaxis()->FindBin(showerY);
          Float_t xMod = module_ID_map->GetXaxis()->GetBinCenter(binx);
          Float_t yMod = module_ID_map->GetYaxis()->GetBinCenter(biny);

          // module number
          int module_number = module_ID_map->GetBinContent(binx,biny);

          int cellMapIndex;
          std::map<int, int> ::const_iterator iterCellMap = mapOfCellMaps.find(showerModuleType);
          if(iterCellMap != mapOfCellMaps.end())
          {
            // iter is item pair in the map. The value will be accessible as `iter->second`.
            cellMapIndex = (iterCellMap->second);
          }
          else
          {
            if(verbose > 0)
            {
              std::cout << "WARNING: Map of cell maps not found for showerModuleType " << showerModuleType << std::endl;
              std::cout << "         This event will be ignored! Check the input files provided." << std::endl;
              std::cout << std::endl;
            }
            lostScintillation++;
            continue;
          }



          // cell number is done later...


          // int cell_global_number = moduleNumberOfBeforeMap[showerModuleType] + cell_number;
          // calculate detector from module x y index and cell number

          // std::cout << "iii" << std::endl;
          float deltaT = calibrations[calIndex].calibration_point[xIndex][yIndex][zIndex][eIndex][sIndex].delay->GetRandom();
          // assign a time of arrival
          // sum of energy deposition moment (showerT) + photon generation (scintillator time profile) + propagation time (time pdf)
          // TH1F *tHisto = tHistoMap.at(showerMaterialNumber);
          float generationTime = tHisto->GetRandom();
          float totalT =  showerT + generationTime + deltaT;

          float fakeZ ;


          // distinguish the two possibilities
          // 1) no longitudinal segmentation
          // 2) longitudinal segmentation
          int moduleSections;
          std::map<int, int> ::const_iterator iterModSec = moduleSectionsMap.find(showerModuleType);
          if(iterModSec != moduleSectionsMap.end())
          {
            // iter is item pair in the map. The value will be accessible as `iter->second`.
            moduleSections = (iterModSec->second);
          }
          else
          {
            if(verbose > 0)
            {
              std::cout << "WARNING: Module sections map not found for showerModuleType " << showerModuleType << std::endl;
              std::cout << "         This event will be ignored! Check the input files provided." << std::endl;
              std::cout << std::endl;
            }
            lostScintillation++;
            continue;
          }
          // fetch also module min and max and separation
          float moduleMin;
          std::map<int, float> ::const_iterator iterModMin = moduleMinMap.find(showerModuleType);
          if(iterModMin != moduleMinMap.end())
          {
            // iter is item pair in the map. The value will be accessible as `iter->second`.
            moduleMin = (iterModMin->second);
          }
          else
          {
            if(verbose > 0)
            {
              std::cout << "WARNING: Module min map not found for showerModuleType " << showerModuleType << std::endl;
              std::cout << "         This event will be ignored! Check the input files provided." << std::endl;
              std::cout << std::endl;
            }
            lostScintillation++;
            continue;
          }

          float moduleMax;
          std::map<int, float> ::const_iterator iterModMax = moduleMaxMap.find(showerModuleType);
          if(iterModMax != moduleMaxMap.end())
          {
            // iter is item pair in the map. The value will be accessible as `iter->second`.
            moduleMax = (iterModMax->second);
          }
          else
          {
            if(verbose > 0)
            {
              std::cout << "WARNING: Module max map not found for showerModuleType " << showerModuleType << std::endl;
              std::cout << "         This event will be ignored! Check the input files provided." << std::endl;
              std::cout << std::endl;
            }
            lostScintillation++;
            continue;
          }

          float moduleSeparation;
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
            lostScintillation++;
            continue;
          }



          if(moduleSections == 1)
          {
            // you cannot use the mean of cell Z, it's always 0!
            // so use min and max of the cell, give the sIndex
            // sIndex = 0 means negative z
            // sIndex = 1 means positive z
            if(sIndex == 0)
            {
              fakeZ = moduleMin;

            }
            else
            {
              fakeZ = moduleMax;
            }



          }
          else // then it has to be == 2, given check at the beginning of program
          {

            if(showerZ > moduleSeparation) // "second" part of the module, so exit is for max z
            {
              fakeZ = moduleMax;
            }
            else
            {
              fakeZ = moduleMin;
            }
          }

          Int_t c_binx    = cell_map_per_module[cellMapIndex][sIndex]->GetXaxis()->FindBin(showerX - xMod);
          Int_t c_biny    = cell_map_per_module[cellMapIndex][sIndex]->GetYaxis()->FindBin(showerY - yMod);
          // cell number
          int cell_number = cell_map_per_module[cellMapIndex][sIndex]->GetBinContent(c_binx,c_biny);

          if(verbose > 1)
          {
            // std::cout << "showerCellID " << showerCellID << std::endl;
            std::cout << "cell_number " << cell_number << std::endl;
            std::cout << "module_number " << module_number << std::endl;
            std::cout << "calibrations[calIndex].calibration_point[xIndex][yIndex][zIndex][eIndex][sIndex].delay " << calibrations[calIndex].calibration_point[xIndex][yIndex][zIndex][eIndex][sIndex].delay->GetName() << std::endl;

            std::cout << "showerT " << showerT << std::endl;
            std::cout << "generationTime " << generationTime << std::endl;
            std::cout << "deltaT " << deltaT << std::endl;
            std::cout << "totalT " << totalT << std::endl;
            std::cout << "moduleSections" << moduleSections << std::endl;
            std::cout << "moduleMin " << moduleMin << std::endl;
            std::cout << "moduleMax" << moduleMax << std::endl;
            std::cout << "moduleSeparation" << moduleSeparation << std::endl;
            std::cout << "fakeZ " << fakeZ << std::endl;
          }




          totalNumberOfOpticalPhotonsOut++;
          // same for hybrid scintillation part
          hb_event = showerEvent;
          hb_x = showerX;
          hb_y = showerY;
          hb_z = fakeZ;
          hb_detector = cell_number;
          hb_module   = module_number;
          hb_moduleType = showerModuleType;
          hb_timestamp = totalT;
          hb_t_deposition  = showerT;
          hb_t_generation  = generationTime;
          hb_t_propagation = deltaT;
          hb_photon_energy = photon_energy;
          hb_processName = "Scintillation";
          hb_processNumber = 0;
          hb_primaryPositionOnAbsorberX = shower_primaryPositionOnAbsorberX;
          hb_primaryPositionOnAbsorberY = shower_primaryPositionOnAbsorberY;
          hb_primaryPositionOnAbsorberZ = shower_primaryPositionOnAbsorberZ;
          hb_primaryMomentumOnAbsorberX = shower_primaryMomentumOnAbsorberX;
          hb_primaryMomentumOnAbsorberY = shower_primaryMomentumOnAbsorberY;
          hb_primaryMomentumOnAbsorberZ = shower_primaryMomentumOnAbsorberZ;
          hb_primaryEnergyOnAbsorber    = shower_primaryEnergyOnAbsorber   ;
          outTreeHybrid->Fill();



        }

      } // end of running on all produced
    } // end of showerIsInCrystal


    // counter++;
    // if((counter % (nEntriesShower/10))==0)
    // {
    //   std::cout << ((100*counter)/nEntriesShower) << "% done... " << std::endl;
    // }
    // counter++;
    // if((nEntriesShower/10) != 0)
    // {
    //   if((counter % (nEntriesShower/10))==0)
    //   {
    //     if(nEntriesShower != 0)
    //     {
    //       std::cout << ((100*counter)/nEntriesShower) << "% done... " << std::endl;
    //     }
    //   }
    // }

  }

  std::cout << "Done." << std::endl;
  std::cout << "lostShower                     = " << lostShower << std::endl;
  std::cout << "totalNumberOfOpticalPhotons    = " << totalNumberOfOpticalPhotons << std::endl;
  std::cout << "totalNumberOfOpticalPhotonsOut = " << totalNumberOfOpticalPhotonsOut << std::endl;
  std::cout << "lostScintillation              = " << lostScintillation << std::endl;

 // complete hybrid with cerenkov propagated by simulation
 // -------------------------- //
 // TChain of photons for hybrid //
 // -------------------------- //
 TTree* cer_photons = (TTree*) inputFile->Get("photons");
 std::cout << std::endl;
 int         cer_event   ;
 // int         cer_front_back    ;
 // int         cer_pmt_number    ;
 float       cer_vertX         ;
 float       cer_vertY         ;
 float       cer_vertZ         ;
 float       cer_vertMomentumX ;
 float       cer_vertMomentumY ;
 float       cer_vertMomentumZ ;
 float       cer_PositionX     ;
 float       cer_PositionY     ;
 float       cer_PositionZ     ;
 float       cer_PreMomentumX  ;
 float       cer_PreMomentumY  ;
 float       cer_PreMomentumZ  ;
 float       cer_PostMomentumX ;
 float       cer_PostMomentumY ;
 float       cer_PostMomentumZ ;
 float       cer_globalTime    ;
 float       cer_localTime     ;
 float       cer_PhotonEnergy  ;
 // int         cer_cellID        ;
 // int         cer_moduleID      ;
 int         cer_moduleType    ;
 // float       cer_module_x      ;
 // float       cer_module_y      ;
 // float       cer_module_z      ;
 std::string *cer_processName = 0;
 float cer_primaryPositionOnAbsorberX    ;
 float cer_primaryPositionOnAbsorberY    ;
 float cer_primaryPositionOnAbsorberZ    ;
 float cer_primaryMomentumOnAbsorberX    ;
 float cer_primaryMomentumOnAbsorberY    ;
 float cer_primaryMomentumOnAbsorberZ    ;
 float cer_primaryEnergyOnAbsorber       ;

 TBranch *b_cer_event         ;
 // TBranch *b_cer_front_back    ;
 // TBranch *b_cer_pmt_number    ;
 TBranch *b_cer_vertX         ;
 TBranch *b_cer_vertY         ;
 TBranch *b_cer_vertZ         ;
 TBranch *b_cer_vertMomentumX ;
 TBranch *b_cer_vertMomentumY ;
 TBranch *b_cer_vertMomentumZ ;
 TBranch *b_cer_PositionX     ;
 TBranch *b_cer_PositionY     ;
 TBranch *b_cer_PositionZ     ;
 TBranch *b_cer_PreMomentumX  ;
 TBranch *b_cer_PreMomentumY  ;
 TBranch *b_cer_PreMomentumZ  ;
 TBranch *b_cer_PostMomentumX ;
 TBranch *b_cer_PostMomentumY ;
 TBranch *b_cer_PostMomentumZ ;
 TBranch *b_cer_globalTime    ;
 TBranch *b_cer_localTime     ;
 TBranch *b_cer_PhotonEnergy  ;
 TBranch *b_cer_moduleType  ;
 TBranch *b_cer_primaryPositionOnAbsorberX    ;
 TBranch *b_cer_primaryPositionOnAbsorberY    ;
 TBranch *b_cer_primaryPositionOnAbsorberZ    ;
 TBranch *b_cer_primaryMomentumOnAbsorberX    ;
 TBranch *b_cer_primaryMomentumOnAbsorberY    ;
 TBranch *b_cer_primaryMomentumOnAbsorberZ    ;
 TBranch *b_cer_primaryEnergyOnAbsorber       ;

 cer_photons->SetBranchAddress("event"         , &cer_event         , &b_cer_event);
 // cer_photons->SetBranchAddress("front_back"    , &cer_front_back    , &b_cer_front_back);
 // cer_photons->SetBranchAddress("pmt_number"    , &cer_pmt_number    , &b_cer_pmt_number);
 cer_photons->SetBranchAddress("vertX"         , &cer_vertX         , &b_cer_vertX);
 cer_photons->SetBranchAddress("vertY"         , &cer_vertY         , &b_cer_vertY);
 cer_photons->SetBranchAddress("vertZ"         , &cer_vertZ         , &b_cer_vertZ);
 cer_photons->SetBranchAddress("vertMomentumX" , &cer_vertMomentumX , &b_cer_vertMomentumX);
 cer_photons->SetBranchAddress("vertMomentumY" , &cer_vertMomentumY , &b_cer_vertMomentumY);
 cer_photons->SetBranchAddress("vertMomentumZ" , &cer_vertMomentumZ , &b_cer_vertMomentumZ);
 cer_photons->SetBranchAddress("PositionX"     , &cer_PositionX     , &b_cer_PositionX);
 cer_photons->SetBranchAddress("PositionY"     , &cer_PositionY     , &b_cer_PositionY);
 cer_photons->SetBranchAddress("PositionZ"     , &cer_PositionZ     , &b_cer_PositionZ);
 cer_photons->SetBranchAddress("PreMomentumX"  , &cer_PreMomentumX  , &b_cer_PreMomentumX);
 cer_photons->SetBranchAddress("PreMomentumY"  , &cer_PreMomentumY  , &b_cer_PreMomentumY);
 cer_photons->SetBranchAddress("PreMomentumZ"  , &cer_PreMomentumZ  , &b_cer_PreMomentumZ);
 cer_photons->SetBranchAddress("PostMomentumX" , &cer_PostMomentumX , &b_cer_PostMomentumX);
 cer_photons->SetBranchAddress("PostMomentumY" , &cer_PostMomentumY , &b_cer_PostMomentumY);
 cer_photons->SetBranchAddress("PostMomentumZ" , &cer_PostMomentumZ , &b_cer_PostMomentumZ);
 cer_photons->SetBranchAddress("globalTime"    , &cer_globalTime    , &b_cer_globalTime);
 cer_photons->SetBranchAddress("localTime"     , &cer_localTime     , &b_cer_localTime);
 cer_photons->SetBranchAddress("PhotonEnergy"  , &cer_PhotonEnergy  , &b_cer_PhotonEnergy);
 // cer_photons->SetBranchAddress("cellID"        , &cer_cellID );
 // cer_photons->SetBranchAddress("moduleID"      , &cer_moduleID );
 cer_photons->SetBranchAddress("moduleType"    , &cer_moduleType ,&b_cer_moduleType);
 // cer_photons->SetBranchAddress("module_x"      , &cer_module_x );
 // cer_photons->SetBranchAddress("module_y"      , &cer_module_y );
 // cer_photons->SetBranchAddress("module_z"      , &cer_module_z );
 cer_photons->SetBranchAddress("processName"   , &cer_processName );
 cer_photons->SetBranchAddress("primary_PositionOnAbsorberX" ,&cer_primaryPositionOnAbsorberX    ,&b_cer_primaryPositionOnAbsorberX);
 cer_photons->SetBranchAddress("primary_PositionOnAbsorberY" ,&cer_primaryPositionOnAbsorberY    ,&b_cer_primaryPositionOnAbsorberY);
 cer_photons->SetBranchAddress("primary_PositionOnAbsorberZ" ,&cer_primaryPositionOnAbsorberZ    ,&b_cer_primaryPositionOnAbsorberZ);
 cer_photons->SetBranchAddress("primary_MomentumOnAbsorberX" ,&cer_primaryMomentumOnAbsorberX    ,&b_cer_primaryMomentumOnAbsorberX);
 cer_photons->SetBranchAddress("primary_MomentumOnAbsorberY" ,&cer_primaryMomentumOnAbsorberY    ,&b_cer_primaryMomentumOnAbsorberY);
 cer_photons->SetBranchAddress("primary_MomentumOnAbsorberZ" ,&cer_primaryMomentumOnAbsorberZ    ,&b_cer_primaryMomentumOnAbsorberZ);
 cer_photons->SetBranchAddress("primary_EnergyOnAbsorber" ,&cer_primaryEnergyOnAbsorber       ,&b_cer_primaryEnergyOnAbsorber   );



  long int nEntriesCer = cer_photons->GetEntries();
  long long int lostCerenkov = 0;
  // long int nEntriesOptPhotGen = opticalPhotonGen->GetEntries();
  // std::cout << nEntriesOptPhotGen <<  std::endl;
  std::cout << "Cerenkov photons = " << nEntriesCer <<  std::endl;
  std::cout << "Running on cerenkov photons..." << std::endl;
  long long int totalNumberOfCherekovPhotons = 0;
  counter = 0;
  feedbackDivision = (int) (nEntriesCer/20);
  if(feedbackDivision == 0) feedbackDivision = 1;
  for(int iEntry = 0 ; iEntry < nEntriesCer ; iEntry++)
  {
   if (  iEntry % feedbackDivision == 0 ) { std::cout << (int) ((double) iEntry/nEntriesCer * 100) << "% done..." << std::endl;}
   cer_photons->GetEvent(iEntry);
   // std::cout << cer_processName->data() << std::endl;
   std::string cerString = "Cerenkov";
   // get only cherenkov
   if(cer_processName->data() == cerString)
   {

     // find cell
     // int cell_number = -1;
     // cell_number = cer_cellID;
     //
     // if(cell_number == -1 )
     // {
     //   std::cout << "WARNING: cell not found for Cerenkov event located at " << cer_vertX  << " " << cer_vertY << " " << cer_vertZ << std::endl;
     //   std::cout << std::endl;
     //   lostCerenkov++;
     //   // return 1;
     // }
     totalNumberOfCherekovPhotons++;

     // int module_number      = moduleNumberOfBeforeMap[cer_moduleType] + cer_moduleID;

     Int_t binx = module_ID_map->GetXaxis()->FindBin(cer_PositionX);
     Int_t biny = module_ID_map->GetYaxis()->FindBin(cer_PositionY);
     Float_t xMod = module_ID_map->GetXaxis()->GetBinCenter(binx);
     Float_t yMod = module_ID_map->GetYaxis()->GetBinCenter(biny);

     // module number
     int module_number = module_ID_map->GetBinContent(binx,biny);

     int cellMapIndex;
     std::map<int, int> ::const_iterator iterCellMap = mapOfCellMaps.find(cer_moduleType);
     if(iterCellMap != mapOfCellMaps.end())
     {
       // iter is item pair in the map. The value will be accessible as `iter->second`.
       cellMapIndex = (iterCellMap->second);
     }
     else
     {
       if(verbose > 0)
       {
         std::cout << "WARNING: Map of cell maps not found for cer_moduleType " << cer_moduleType << std::endl;
         std::cout << "         This event will be ignored! Check the input files provided." << std::endl;
         std::cout << std::endl;
       }
       lostCerenkov++;
       continue;
     }


     float moduleSeparation;
     std::map<int, float> ::const_iterator iterModSep = moduleSeparationMap.find(cer_moduleType);
     if(iterModSep != moduleSeparationMap.end())
     {
       // iter is item pair in the map. The value will be accessible as `iter->second`.
       moduleSeparation = (iterModSep->second);
     }
     else
     {
       if(verbose > 0)
       {
         std::cout << "WARNING: Module separation not found for cer_moduleType " << cer_moduleType << std::endl;
         std::cout << "         This event will be ignored! Check the input files provided." << std::endl;
         std::cout << std::endl;
       }
       lostCerenkov++;
       continue;
     }



     int sIndex = 0;
     if(cer_PositionZ < moduleSeparation)
     {
       sIndex = 0;
     }
     else
     {
       sIndex = 1;
     }

     Int_t c_binx    = cell_map_per_module[cellMapIndex][sIndex]->GetXaxis()->FindBin(cer_PositionX - xMod);
     Int_t c_biny    = cell_map_per_module[cellMapIndex][sIndex]->GetYaxis()->FindBin(cer_PositionY - yMod);
     // cell number
     int cell_number = cell_map_per_module[cellMapIndex][sIndex]->GetBinContent(c_binx,c_biny);



     // save cerenkov
     hb_x = cer_PositionX;
     hb_y = cer_PositionY;
     hb_z = cer_PositionZ;

     hb_event = cer_event;
     hb_detector = cell_number;
     hb_module   = module_number;
     hb_moduleType = cer_moduleType;
     hb_timestamp = cer_globalTime;
     hb_t_deposition  = cer_globalTime - cer_localTime; // deposition is global time of arrival minus propagation time (since generation time is 0)
     hb_t_generation  = 0.0; // Cerenkov generation is instantaneous
     hb_t_propagation = cer_localTime; // propagation time is local time
     hb_photon_energy = cer_PhotonEnergy;
     hb_processName = "Cerenkov";
     hb_processNumber = 1;
     hb_primaryPositionOnAbsorberX = cer_primaryPositionOnAbsorberX;
     hb_primaryPositionOnAbsorberY = cer_primaryPositionOnAbsorberY;
     hb_primaryPositionOnAbsorberZ = cer_primaryPositionOnAbsorberZ;
     hb_primaryMomentumOnAbsorberX = cer_primaryMomentumOnAbsorberX;
     hb_primaryMomentumOnAbsorberY = cer_primaryMomentumOnAbsorberY;
     hb_primaryMomentumOnAbsorberZ = cer_primaryMomentumOnAbsorberZ;
     hb_primaryEnergyOnAbsorber    = cer_primaryEnergyOnAbsorber   ;

     outTreeHybrid->Fill();

   }


   // counter++;
   // if((counter % (nEntriesCer/10))==0)
   // {
   //   std::cout << ((100*counter)/nEntriesCer) << "% done... " << std::endl;
   // }
   // counter++;
   // if((nEntriesCer/10) != 0)
   // {
   //   if((counter % (nEntriesCer/10))==0)
   //   {
   //     if(nEntriesCer != 0)
   //     {
   //       std::cout << ((100*counter)/nEntriesCer) << "% done... " << std::endl;
   //     }
   //   }
   // }

 }

 std::cout << "done" << std::endl;
 std::cout << "totalNumberOfCherekovPhotons = " << totalNumberOfCherekovPhotons << std::endl;
 std::cout << "lostCerenkov                 = " << lostCerenkov << std::endl;





  std::cout << "Writing output to file " << outputFileName << std::endl;
  // std::map<int,TH1F*>::iterator it_e = eHistoMap.begin();
  // std::map<int,TH1F*>::iterator it_t = tHistoMap.begin();

  std::stringstream sPhotonSeed;
  sPhotonSeed << randomSeed;
  TNamed PhotonSeedNameD("PhotonSeed",sPhotonSeed.str().c_str());


  outputFile->cd();
  // outTreePara->Write();
  // outTreeFromEnergyDep->Write();
  outTreeHybrid->Write();
  out_shower->Write();
  for (std::map<int,TH1F*>::iterator it=eHistoMap.begin(); it!=eHistoMap.end(); ++it)
    it->second->Write();
  for (std::map<int,TH1F*>::iterator it=tHistoMap.begin(); it!=tHistoMap.end(); ++it)
    it->second->Write();
  // while (it_e != eHistoMap.end())
	// {
	// 	it_e->second->Write();
  // }
  // while (it_t != tHistoMap.end())
	// {
	// 	it_t->second->Write();
  // }

  // enHisto->Write();
  // tHisto->Write();
  out_modules->Write();
  out_absorbers->Write();
  out_cells->Write();
  out_holes->Write();
  out_fibres->Write();
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


  TDirectory *histoConf = outputFile->mkdir("Configuration");
  histoConf->cd();
  SeedNameD.Write();
  HostNameD.Write();
  PWDNameD.Write();
  ConfigNameD.Write();
  GpsNameD.Write();
  PrimariesNameD.Write();
  PhotonSeedNameD.Write();
  TNamed FeedbackNameD("Parameters",feedbackString.str().c_str());
  TNamed CommandLineNameD("CommandLine",CommandLineString.c_str());
  FeedbackNameD.Write();
  CommandLineNameD.Write();

  // FIXME - extend to multiple cali input
  // gDirectory->cd();
  // TDirectory *histoDir = outputFile->mkdir("Hybrid_Calibration");
  // histoDir->cd();
  // for(int ix = 0 ; ix < x.size() ; ix++)
  // {
  //   for(int iy = 0 ; iy < y.size() ; iy++)
  //   {
  //     for(int iz = 0 ; iz < z.size() ; iz++)
  //     {
  //       for(int ie = 0 ; ie < energy.size() ; ie++)
  //       {
  //         for(int iSide = 0; iSide < nSides ; iSide++)
  //         {
  //           calibration_point[ix][iy][iz][ie][iSide].delay->Write();
  //         }
  //       }
  //     }
  //   }
  // }
  outputFile->Close();
  std::cout << "Done. Goodbye!" << std::endl;
  // need to kill ROOT badly, but this needs a fix
  // too much is loaded in memory, most likely
  gROOT->ProcessLine(".qqqqqqqqqqqqqqqqqqqq");
  // explicitly delete
  // for(int i = 0 ; i < calibrations.size();i++)
  // {
  //   for(int ix = 0 ; ix < calibrations[i].xn ; ix++)
  //   {
  //     for(int iy = 0 ; iy < calibrations[i].yn ; iy++)
  //     {
  //       for(int iz = 0 ; iz < calibrations[i].zn ; iz++)
  //       {
  //         for(int ie = 0 ; ie < calibrations[i].en ; ie++)
  //         {
  //           for(int iSide = 0; iSide < 2 ; iSide++)
  //           {
  //             delete calibrations[i].calibration_point[ix][iy][iz][ie][iSide].delay;
  //           }
  //         }
  //       }
  //     }
  //   }
  // }
  //
  // std::cout << "Done deleting!" << std::endl;

  return 0;
}


void importCalibrationFile(std::string calibrationFileName,calibration_t &calibration)
{
  //----------------------------------//
  // OPTICAL CALIBRATION FILE         //
  //----------------------------------//
  std::ifstream input_file(calibrationFileName.c_str(), std::ios::binary);
  calibration.fileName = calibrationFileName.c_str();
  // read header
  int bins;
  float pdfStart,pdfEnd;
  int xn;
  int yn;
  int zn;
  int en;
  float x_v,y_v,z_v,e_v;
  std::vector <float > x,y,z,energy;

  std::cout << "Reading calibration file " << calibrationFileName.c_str() << std::endl;
  // store ints
  input_file.read( (char*)&bins , sizeof(int));
  input_file.read( (char*)&pdfStart , sizeof(float));
  input_file.read( (char*)&pdfEnd , sizeof(float));
  input_file.read( (char*)&xn , sizeof(int));
  input_file.read( (char*)&yn , sizeof(int));
  input_file.read( (char*)&zn , sizeof(int));
  input_file.read( (char*)&en , sizeof(int));
  std::cout << bins << " "
            << pdfStart << " "
            << pdfEnd << " "
            << xn << " "
            << yn << " "
            << zn << " "
            << en << " "
            << std::endl;
  // read vectors

  for(int i = 0 ; i < xn; i++)
  {
    input_file.read( (char*)&x_v , sizeof(float));
    x.push_back(x_v);
  }
  for(int i = 0 ; i < yn; i++)
  {
    input_file.read( (char*)&y_v , sizeof(float));
    y.push_back(y_v);
  }
  for(int i = 0 ; i < zn; i++)
  {
    input_file.read( (char*)&z_v , sizeof(float));
    z.push_back(z_v);
  }
  for(int i = 0 ; i < en; i++)
  {
    input_file.read( (char*)&e_v , sizeof(float));
    energy.push_back(e_v);
  }

  float xmin   = 0;
  float xstep  = 1;
  float ymin   = 0;
  float ystep  = 1;
  float zmin   = 0;
  float zstep  = 1;
  float emin   = 1;
  float estep  = 1;
  int   nSides = 2;

  xmin = x[0];
  ymin = y[0];
  zmin = z[0];
  emin = energy[0];
  if(x.size() == 1)
  {
    xstep = 1;
  }
  else
  {
    xstep = x[1] - x[0];
  }
  if(y.size() == 1)
  {
    ystep = 1;
  }
  else
  {
    ystep = y[1] - y[0];
  }
  if(z.size() == 1)
  {
    zstep = 1;
  }
  else
  {
    zstep = z[1] - z[0];
  }
  if(energy.size() == 1)
  {
    estep = 1;
  }
  else
  {
    estep = energy[1] - energy[0];
  }



  // FEEDBACK
  std::cout << std::endl;
  std::cout << "----------------------------- " << std::endl;
  std::cout << "Calibration points " << std::endl;
  std::cout << "----------------------------- " << std::endl;
  std::cout << std::endl;
  std::cout << "---> X: " << std::endl;
  std::cout << "Number of points in X           = " << xn << std::endl;
  std::cout << "Step length           [mm]      = " << xstep << std::endl;
  std::cout << "Starting X coordinate [mm]      = " << xmin  << std::endl;
  std::cout << "List of X points      [mm]      = " ;
  for(unsigned int ix = 0; ix < x.size(); ix++)
  {
    std::cout << x[ix] << " ";
  }
  std::cout << std::endl;

  std::cout << "---> Y: " << std::endl;
  std::cout << "Number of points in Y           = " << yn << std::endl;
  std::cout << "Step length           [mm]      = " << ystep << std::endl;
  std::cout << "Starting Y coordinate [mm]      = " << ymin  << std::endl;
  std::cout << "List of Y points      [mm]      = " ;
  for(unsigned int iy = 0; iy < y.size(); iy++)
  {
    std::cout << y[iy] << " ";
  }
  std::cout << std::endl;

  std::cout << "---> Z: " << std::endl;
  std::cout << "Number of points in Z           = " << zn << std::endl;
  std::cout << "Step length           [mm]      = " << zstep << std::endl;
  std::cout << "Starting Z coordinate [mm]      = " << zmin  << std::endl;
  std::cout << "List of Z points      [mm]      = " ;
  for(unsigned int iz = 0; iz < z.size(); iz++)
  {
    std::cout << z[iz] << " ";
  }
  std::cout << std::endl;
  std::cout << "---> ENERGY: " << std::endl;
  std::cout << "Number of points in ENERGY      = " << en << std::endl;
  std::cout << "Step length                [eV] = " << estep << std::endl;
  std::cout << "Starting ENERGY coordinate [eV] = " << emin  << std::endl;
  std::cout << "List of ENERGY points      [eV] = " ;
  for(unsigned int ie = 0; ie < energy.size(); ie++)
  {
    std::cout << energy[ie] << " ";
  }
  std::cout << std::endl;
  std::cout << "----------------------------- " << std::endl;
  std::cout << std::endl;
  // ------------

  // copy vectors
  calibration.x = x;
  calibration.y = y;
  calibration.z = z;
  calibration.energy = energy;
  // std::cout << "Waiting for ENTER.." << std::endl;
  // getchar();

  // prepare calibration structs
  // calibration_point_t *****calibration_point;
  calibration.calibration_point = new calibration_point_t****[x.size()];
  for(int ix = 0 ; ix < x.size() ; ix++)
  {
    calibration.calibration_point[ix] = new calibration_point_t***[y.size()];
    for(int iy = 0 ; iy < y.size() ; iy++)
    {
      calibration.calibration_point[ix][iy] = new calibration_point_t**[z.size()];
      for(int iz = 0 ; iz < z.size() ; iz++)
      {
        calibration.calibration_point[ix][iy][iz] = new calibration_point_t*[energy.size()];
        for(int ie = 0 ; ie < energy.size() ; ie++)
        {
          calibration.calibration_point[ix][iy][iz][ie] = new calibration_point_t[nSides];
          for(int iSide = 0; iSide < nSides ; iSide++)
          {
            std::stringstream ssname;
            ssname << "delayHisto_moduleType_" << calibration.type << " " << x[ix] << " " << y[iy] << " " << z[iz] << " " << energy[ie] << " " << iSide;
            // calibration.calibration_point[ix][iy][iz][ie][iSide].x = x[ix];
            // calibration.calibration_point[ix][iy][iz][ie][iSide].y = y[iy];
            // calibration.calibration_point[ix][iy][iz][ie][iSide].z = z[iz];
            // calibration.calibration_point[ix][iy][iz][ie][iSide].energy = energy[ie];
            calibration.calibration_point[ix][iy][iz][ie][iSide].eff = 0;
            calibration.calibration_point[ix][iy][iz][ie][iSide].delay = new TH1I(ssname.str().c_str(),ssname.str().c_str(),bins,pdfStart,pdfEnd);
          }
        }
      }
    }
  }

  // std::cout << "Waiting for ENTER 2.." << std::endl;
  // getchar();

  int CaliPoints = xn*yn*zn*en*nSides;
  int feedbackDivision = (int) (CaliPoints/20);
  if(feedbackDivision == 0) feedbackDivision = 1;
  std::cout << "Reading calibration points..." << std::endl;
  for(int iCali = 0; iCali < CaliPoints; iCali++)
  {
    // if (  iCali % feedbackDivision == 0 ) { std::cout << (int) ((double) iCali/CaliPoints * 100) << "% done..." << std::endl;}
    int ix   ;
    int iy   ;
    int iz   ;
    int ie   ;
    int iSide;
    float x_i  ;
    float y_i  ;
    float z_i  ;
    int e_i    ;
    int s_i    ;
    float eff;
    float value ;

    input_file.read( (char*)&ix , sizeof(int));
    input_file.read( (char*)&iy , sizeof(int));
    input_file.read( (char*)&iz , sizeof(int));
    input_file.read( (char*)&ie , sizeof(int));
    input_file.read( (char*)&iSide , sizeof(int));

    input_file.read( (char*)&x_i , sizeof(float));
    input_file.read( (char*)&y_i , sizeof(float));
    input_file.read( (char*)&z_i , sizeof(float));
    input_file.read( (char*)&s_i , sizeof(int));
    input_file.read( (char*)&e_i , sizeof(float));
    input_file.read( (char*)&eff , sizeof(float));

    // std:: cout << ix << " " << " " << iy << " " << iz << " " << ie << " " << iSide << std::endl;

    // calibration.calibration_point[ix][iy][iz][ie][iSide].x = x_i;
    // calibration.calibration_point[ix][iy][iz][ie][iSide].y = y_i;
    // calibration.calibration_point[ix][iy][iz][ie][iSide].z = z_i;
    // calibration.calibration_point[ix][iy][iz][ie][iSide].s = s_i;
    // calibration.calibration_point[ix][iy][iz][ie][iSide].energy = e_i;
    // calibration.calibration_point[ix][iy][iz][ie][iSide].primaries = 0; // not needed now
    calibration.calibration_point[ix][iy][iz][ie][iSide].eff = eff;
    // std::cout << eff << std::endl;

    for(int b = 1 ; b < bins+1; b++)
    {
      input_file.read( (char*)&value , sizeof(float));
      calibration.calibration_point[ix][iy][iz][ie][iSide].delay->SetBinContent(b,(Int_t) value);
    }

     // std::cout << calibration_point[ix][iy][iz][ie][iSide].x << std::endl;
    // std::cout << calibration_point_out.binContents[2] << std::endl;
    // for(int i = 0 ; i < bins ; i++)
    // {
      // calibration_point[ix][iy][iz][ie][iSide].delay->SetBinContent(i+1,calibration_point_out.binContents[i]);
    // }
    // calibration_point[ix][iy][iz][ie][iSide].delay;

  }
  // std::cout << "Waiting for ENTER 3.." << std::endl;
  // getchar();

  // save values in cali struct

  calibration.xmin   = xmin   ;
  calibration.xstep  = xstep  ;
  calibration.ymin   = ymin   ;
  calibration.ystep  = ystep  ;
  calibration.zmin   = zmin   ;
  calibration.zstep  = zstep  ;
  calibration.emin   = emin   ;
  calibration.estep  = estep  ;
  calibration.nSides = nSides ;
  calibration.xn     = xn     ;
  calibration.yn     = yn     ;
  calibration.zn     = zn     ;
  calibration.en     = en     ;
  

  // try to read what comes after 
  float moduleCenter;
  input_file.read( (char*)&moduleCenter, sizeof(float));

  if(input_file)
  {
    std::cout << "New type of calibration data file, I will read the footer too " << std::endl;
    
    // flag footer to true
    calibration.footer = true;
    
    // declare values

    // read values 
        
    // write in the structure
    calibration.moduleCenter = moduleCenter;
    relevant_data_t para;
    input_file.read( (char*)&para, sizeof(relevant_data_t));
    calibration.parameters = para;
  }
  else 
  {
    std::cout << "Old type of calibration data file"<< std::endl;
  }
  input_file.close();
  std::cout << "Finished importing calibration points." << std::endl;
  std::cout << std::endl;
  // ---------------------------------//


}

void trim( std::string& st )
{
	// Remove leading and trailing whitespace
	static const char whitespace[] = " \n\t\v\r\f";
	st.erase( 0, st.find_first_not_of(whitespace) );
	st.erase( st.find_last_not_of(whitespace) + 1U );
}

void usage()
{
  std::cout << "\t\t" << "[-i | --input]          <input file>   " << std::endl
            << "\t\t" << "[-c | --calibration]    <CSV list of optical calibration files>    " << std::endl
            << "\t\t" << "[--type]                <CSV list of module types corresponding to calibration files, in the same order of calibration list!!> " << std::endl
            << "\t\t" << "[-o | --output]         <output file name>    " << std::endl
            << "\t\t" << "[-p | --photonSeed]     <seed for TRandom3>   - default = 0 , i.e. seed automatically computed via a TUUID object" << std::endl
            << "\t\t" << "[--external]            <external file name>  - external ROOT file with data for another crystal material (use with caution...) - default = ignored" << std::endl
            << "\t\t" << "[--verbose]             <verbosity level (0, 1 or 2)> - default = 0  " << std::endl
            // << "\t\t" << "[--repo]                <folder where with library of calibrations - default = \"\" " << std::endl
            << "\t\t" << std::endl;
}


//*-*-*
//------------------------------------------------
// G4Scintillation implementation:
// G4double MeanNumberOfPhotons;
//
//   // Birk's correction via fEmSaturation and specifying scintillation by
//   // by particle type are physically mutually exclusive
//
//   if (fScintillationByParticleType)
//      MeanNumberOfPhotons = ScintillationYield;
//   else if (fEmSaturation)
//      MeanNumberOfPhotons = ScintillationYield*
//        (fEmSaturation->VisibleEnergyDepositionAtAStep(&aStep));
//   else
//      MeanNumberOfPhotons = ScintillationYield*TotalEnergyDeposit;
//
//   G4int NumPhotons;
//
//   if (MeanNumberOfPhotons > 10.)
//   {
//     G4double sigma = ResolutionScale * std::sqrt(MeanNumberOfPhotons);
//     NumPhotons = G4int(G4RandGauss::shoot(MeanNumberOfPhotons,sigma)+0.5);
//   }
//   else
//   {
//     NumPhotons = G4int(G4Poisson(MeanNumberOfPhotons));
//   }
//------------------------------------------------

// ===============> from literature, it's not clear if Birks' law applies to inorganic scintillators.
// so use standard one
