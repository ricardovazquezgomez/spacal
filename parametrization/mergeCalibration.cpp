// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch

// compile with
// g++ -o ../build/mergeCalibration mergeCalibration.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore


// program to merge optical calibration files


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
#include "TEnv.h"



#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <getopt.h>

#include <sys/stat.h>
#include <dirent.h>

#include "regions.h"
#include "Data.cpp"

void usage();

struct calibration_point_t
{
  float x;
  float y;
  float z;
  int   s;
  float energy;
  int primaries;
  float eff;
  // std::vector<float> globalTimes;
  TH1F *delay;
};

struct calibration_point_output_t
{
  int ix;
  int iy;
  int iz;
  int is;
  int ie;
  float x;
  float y;
  float z;
  int   s;
  float e;
  float eff;
  float *binCenters;
  float *binContents;
};




int main(int argc, char** argv)
{

  if(argc < 2)
  {
    std::cout << argv[0] << std::endl;
    usage();
    return 1;
  }

  gEnv->GetValue("TFile.Recover", 0);

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


  std::string inputFilePrefix = "";
  std::string inputFolderName = "./";
  std::string outputFileName = "calibration.data";
  int NofPrimaries = 100000;
  float xmin  = 0;
  float xstep = 1;
  int   xn    = 1;
  float ymin  = 0;
  float ystep = 1;
  int   yn    = 1;
  float zmin  = 0;
  float zstep = 1;
  int   zn    = 1;
  float emin  = 1;
  float estep = 1;
  int   en    = 1;
  int nSides  = 2;
  int pdfBins  = 5000;
  float pdfStart = 0; // no reason for not starting in 0!
  float pdfEnd = 10; // in ns

  float ScintillationYield = 0.0 ; // GaGG_Ce_Mg
  float ResolutionScale = 1.0;   // for now, but needs to be changed!

  UInt_t photonSeed = 0;

  int strategy = 0;


  static struct option longOptions[] =
  {
    { "folder", required_argument, 0, 0 },
    { "prefix", required_argument, 0, 0 },
    { "output", required_argument, 0, 0 },
    { "pdfStart", required_argument, 0, 0 },
    { "pdfEnd", required_argument, 0, 0 },
    { NULL, 0, 0, 0 }
  };

  while(1) {
    int optionIndex = 0;
    int c = getopt_long(argc, argv, "f:p:o:", longOptions, &optionIndex);
    if (c == -1) {
      break;
    }
    if (c == 'f'){
      inputFolderName = (char *)optarg;
    }
    else if (c == 'p'){
      inputFilePrefix = (char *)optarg;
    }
    else if (c == 'o'){
      outputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 0){
      inputFolderName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 1){
      inputFilePrefix = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 2){
      outputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 3){
      pdfBins = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 4){
      pdfEnd = atof((char *)optarg);
    }
    else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
      usage();
      return 1;
    }
  }

  if(inputFilePrefix == "")
  {
    std::cout << std::endl;
    std::cout << "ERROR! You need to provide an input file prefix!" << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }


  std::stringstream feedbackString;
  feedbackString << "|***************************************************************|" << std::endl;
  feedbackString << "| HYBRID MONTE CARLO FOR SPACAL-RD                              |" << std::endl;
  feedbackString << "| MERGING OPTICAL CALIBRATION FILES                             |" << std::endl;
  feedbackString << "|***************************************************************|" << std::endl;
  feedbackString << "" << std::endl;
  feedbackString << "COMMAND LINE PARAMETERS:" << std::endl;
  feedbackString << "Input folder                = " << inputFolderName <<std::endl;
  feedbackString << "Input file prefix           = " << inputFilePrefix <<std::endl;
  feedbackString << "Output file                 = " << outputFileName <<std::endl;
  feedbackString << "PDF bins                    = " << pdfBins <<std::endl;
  feedbackString << "PDF end                     = " << pdfEnd <<std::endl;
  feedbackString << "" << std::endl;
  // trandom
  // TRandom3 *rand = new TRandom3(photonSeed);
  // UInt_t randomSeed = rand->GetSeed();

  

  //----------------------------------//
  // OPTICAL CALIBRATION FILEs        //
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
    std::cout << "ERROR! Input files do not exists! Aborting." << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }

  // assuming calibration has been done on the same gps (and if not, the user is to blame...)
  // we take the nxyze vectors from first file
  TFile *firstFile = TFile::Open(listInputFiles[0].c_str());
  std::vector<int>   *pN;
  std::vector<float> *pX;
  std::vector<float> *pY;
  std::vector<float> *pZ;
  std::vector<float> *pE;
  gDirectory->GetObject("optN",pN);
  gDirectory->GetObject("optX",pX);
  gDirectory->GetObject("optY",pY);
  gDirectory->GetObject("optZ",pZ);
  gDirectory->GetObject("optE",pE);
  std::vector<int> primaries = pN[0];
  std::vector<float> AllX = pX[0];
  std::vector<float> AllY = pY[0];
  std::vector<float> AllZ = pZ[0];
  std::vector<float> AllE = pE[0];
  
  // also get the modules and absorbers values 
  TTree *modules = (TTree*)   firstFile->Get("modules");
  TTree *absorbers = (TTree*) firstFile->Get("absorbers");
  TTree *cells = (TTree*)     firstFile->Get("cells");
  TTree *holes = (TTree*)     firstFile->Get("holes");
  TTree *fibres = (TTree*)    firstFile->Get("fibres");
  TTree* shower = (TTree*)    firstFile->Get("shower");
  // and extract absorber and module z position

  // load also the parameters ttree, if it's there
  TTree* parametersTree =NULL;
  parametersTree = (TTree*)    firstFile->Get("parameters");
  Data_t parameters; // this is the footer for optical calibration files 
  bool parametersFound = false;
  // and read...
  if(parametersTree)
  {
    parametersFound = true;
    feedbackString << ">>>>> parameters found " << std::endl;
    // read parameters and fill struct     
    readParametersFillStruc(parametersTree,parameters);
    // in optical calibrations, there MUST be only 1 entry
    parametersTree->GetEntry(0);
  }
  


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
      // moduleCenterZ.push_back(modules_regions[iMod].position_z);
    }
  }

  // calc absorbers 
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
  std::map<int, int>   absoberNumberMap; // map from module type to number of modules
  std::map<int, float> absoberPositionMap; // map from module type to z center

  for(int i = 0 ; i < absorbersTypes.size(); i++)
  {
    int modCounter = 0;
    for(int iAbs = 0 ; iAbs < absorbers_regions.size(); iAbs++)
    {
      if(absorbers_regions[iAbs].type == absorbersTypes[i])
      {
        modCounter++;
      }
    }
    absoberNumberMap.insert(std::make_pair(absorbersTypes[i],modCounter));
    absoberPositionMap.insert(std::make_pair(absorbersTypes[i],absorbersPositionZ[i]));
  }

  // lots of maps
  std::map<int, int>   moduleNumberMap; // map from module type to number of modules
  std::map<int, int>   moduleNumberOfBeforeMap; // map from module type to number of other modules
  std::map<int, int>   moduleSectionsMap; // map from module type to number of segments
  std::map<int, float> moduleMinMap; // map from module type to min z
  std::map<int, float> moduleMaxMap; // map from module type to max z
  std::map<int, float> moduleSeparationMap; // map from module type to z of separation
  std::map<int, float> modulePositionMap; // map from module type to z center

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

  // sum all modules
  // std::map<int,int>::iterator it_moduleNumberMap = moduleNumberMap.begin();

  int totModules = 0;
  for (std::map<int,int>::iterator it=moduleNumberMap.begin(); it!=moduleNumberMap.end(); ++it)
  {
    moduleNumberOfBeforeMap.insert(std::make_pair(it->first,totModules));
    totModules += it->second;
  }

  feedbackString << "Total number of modules                 = "<< totModules << std::endl;
  feedbackString << "Total number of abs                     = "<< absorbersTypes.size() << std::endl;




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
    feedbackString << "Type " << moduleTypes[i] << " position at " << modulePositionMap.at(moduleTypes[i]) << " mm " << std::endl;
    feedbackString << "Type " << moduleTypes[i] << " separation at " << moduleSeparationMap.at(moduleTypes[i]) << " mm " << std::endl;
  }
  


  for(int i = 0 ; i < absorbersTypes.size(); i++)
  {
    feedbackString << "Abs type " << absorbersTypes[i] << " position in module at " << absoberPositionMap.at(absorbersTypes[i]) << " mm " << std::endl;
  }

  // calculate abs center (here assuming that in optical calibrations, there is only 1 module and 1 absorber...)
  float moduleCenter = modulePositionMap.at(0) + absoberPositionMap.at(0);
  feedbackString << "Abs center in world " << moduleCenter << " mm " << std::endl;
  
  std::cout << feedbackString.str() << std::endl;

  // fill x y z e vectors, avoiding repetitions
  std::vector<float> x,y,z,energy;
  // x
  for(unsigned int i = 0; i < AllX.size(); i++)
  {
    bool isAlreadyThere = false;
    for(unsigned int j = 0; j < x.size(); j++)
    {
      if(x[j] == AllX[i])
      {
        isAlreadyThere = true;
      }
    }
    if(!isAlreadyThere)
    {
      x.push_back(AllX[i]);
    }
  }
  // y
  for(unsigned int i = 0; i < AllY.size(); i++)
  {
    bool isAlreadyThere = false;
    for(unsigned int j = 0; j < y.size(); j++)
    {
      if(y[j] == AllY[i])
      {
        isAlreadyThere = true;
      }
    }
    if(!isAlreadyThere)
    {
      y.push_back(AllY[i]);
    }
  }
  // z
  for(unsigned int i = 0; i < AllZ.size(); i++)
  {
    bool isAlreadyThere = false;
    for(unsigned int j = 0; j < z.size(); j++)
    {
      if(z[j] == AllZ[i])
      {
        isAlreadyThere = true;
      }
    }
    if(!isAlreadyThere)
    {
      z.push_back(AllZ[i]);
    }
  }
  // energy
  for(unsigned int i = 0; i < AllE.size(); i++)
  {
    bool isAlreadyThere = false;
    for(unsigned int j = 0; j < energy.size(); j++)
    {
      if(energy[j] == AllE[i])
      {
        isAlreadyThere = true;
      }
    }
    if(!isAlreadyThere)
    {
      energy.push_back(AllE[i]);
    }
  }

  // sort vectors, calc num, min and step
  sort(x     .begin(),x     .end());
  sort(y     .begin(),y     .end());
  sort(z     .begin(),z     .end());
  sort(energy.begin(),energy.end());
  xn = x.size();
  yn = y.size();
  zn = z.size();
  en = energy.size();
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

  TFile *outputFile = new TFile("calibration_histograms.root","RECREATE");
  outputFile->cd();
  // prepare calibration structs
  calibration_point_t *****calibration_point;
  calibration_point = new calibration_point_t****[x.size()];
  for(int ix = 0 ; ix < x.size() ; ix++)
  {
    calibration_point[ix] = new calibration_point_t***[y.size()];
    for(int iy = 0 ; iy < y.size() ; iy++)
    {
      calibration_point[ix][iy] = new calibration_point_t**[z.size()];
      for(int iz = 0 ; iz < z.size() ; iz++)
      {
        calibration_point[ix][iy][iz] = new calibration_point_t*[energy.size()];
        for(int ie = 0 ; ie < energy.size() ; ie++)
        {
          calibration_point[ix][iy][iz][ie] = new calibration_point_t[nSides];
          for(int iSide = 0; iSide < nSides ; iSide++)
          {
            std::stringstream ssname;
            ssname << "delayHisto " << x[ix] << " " << y[iy] << " " << z[iz] << " " << energy[ie] << " " << iSide;
            calibration_point[ix][iy][iz][ie][iSide].x = x[ix];
            calibration_point[ix][iy][iz][ie][iSide].y = y[iy];
            calibration_point[ix][iy][iz][ie][iSide].z = z[iz];
            calibration_point[ix][iy][iz][ie][iSide].s = iSide;
            calibration_point[ix][iy][iz][ie][iSide].energy = energy[ie];
            calibration_point[ix][iy][iz][ie][iSide].primaries = 0;
            calibration_point[ix][iy][iz][ie][iSide].eff = 0;
            // calibration_point[ix][iy][iz][ie][iSide].globalTimes.clear();
            calibration_point[ix][iy][iz][ie][iSide].delay = new TH1F(ssname.str().c_str(),ssname.str().c_str(),pdfBins,pdfStart,pdfEnd);
          }
        }
      }
    }
  }

  TFile *inputFile = ((TFile *)0);
  TTree *photons = ((TTree*) 0);
  // now run on all input files
  for(unsigned int i = 0 ; i < listInputFiles.size(); i++)
  {
    // if (!(listInputFiles && listInputFiles[i] && (*(listInputFiles[i])))) continue; // just a precaution
    if (inputFile) delete inputFile; // just a precaution
    if(i == 0)
    {
      // keep the first file...
      inputFile = firstFile;
    }
    else
    {
      inputFile = new TFile(listInputFiles[i].c_str());
    }
    if (!inputFile || inputFile->IsZombie() || inputFile->TestBit(TFile::kRecovered))
    {
      std::cout << "Skipping file " << listInputFiles[i] << std::endl;
      continue;
    }
    std::cout << "Opening file " << listInputFiles[i] << std::endl;
    inputFile->cd();
    // run on points vector, increment the primaries per point
    // again taking the vectors...
    std::vector<int>   *input_pN;
    std::vector<float> *input_pX;
    std::vector<float> *input_pY;
    std::vector<float> *input_pZ;
    std::vector<float> *input_pE;
    gDirectory->GetObject("optN",input_pN);
    gDirectory->GetObject("optX",input_pX);
    gDirectory->GetObject("optY",input_pY);
    gDirectory->GetObject("optZ",input_pZ);
    gDirectory->GetObject("optE",input_pE);
    std::vector<int>   input_primaries = input_pN[0];
    std::vector<float> input_AllX      = input_pX[0];
    std::vector<float> input_AllY      = input_pY[0];
    std::vector<float> input_AllZ      = input_pZ[0];
    std::vector<float> input_AllE      = input_pE[0];
    for(unsigned int i = 0; i < input_primaries.size(); i++)
    {
      int xIndex = (int) roundf((input_AllX[i] - xmin)/xstep);
      int yIndex = (int) roundf((input_AllY[i] - ymin)/ystep);
      int zIndex = (int) roundf((input_AllZ[i] - zmin)/zstep);
      int eIndex = (int) roundf((input_AllE[i] - emin)/estep);
      calibration_point[xIndex][yIndex][zIndex][eIndex][0].primaries += input_primaries[i];
      calibration_point[xIndex][yIndex][zIndex][eIndex][1].primaries += input_primaries[i];
    }
    // TFile *inputFile = TFile::Open(listInputFiles[i].c_str());
    // if(photons) delete photons;
    photons = (TTree*) inputFile->Get("photons");
    if (!photons) continue; // requested ROOT file does not exist or is unreadable
    // std::cout << photons << std::endl;
    // set variables and branches
    int   run           ;
    int   event         ;
    float vertX         ;
    float vertY         ;
    float vertZ         ;
    float PositionX     ;
    float PositionY     ;
    float PositionZ     ;
    float globalTime    ;
    float PhotonEnergy  ;
    TBranch *b_run           ;
    TBranch *b_event         ;
    TBranch *b_vertX         ;
    TBranch *b_vertY         ;
    TBranch *b_vertZ         ;
    TBranch *b_PositionX     ;
    TBranch *b_PositionY     ;
    TBranch *b_PositionZ     ;
    TBranch *b_globalTime    ;
    TBranch *b_PhotonEnergy  ;
    // set addresses
    photons->SetBranchAddress("run"           , &run           , &b_run);
    photons->SetBranchAddress("event"         , &event         , &b_event);
    photons->SetBranchAddress("vertX"         , &vertX         , &b_vertX);
    photons->SetBranchAddress("vertY"         , &vertY         , &b_vertY);
    photons->SetBranchAddress("vertZ"         , &vertZ         , &b_vertZ);
    photons->SetBranchAddress("PositionX"     , &PositionX     , &b_PositionX);
    photons->SetBranchAddress("PositionY"     , &PositionY     , &b_PositionY);
    photons->SetBranchAddress("PositionZ"     , &PositionZ     , &b_PositionZ);
    photons->SetBranchAddress("globalTime"    , &globalTime    , &b_globalTime);
    photons->SetBranchAddress("PhotonEnergy"  , &PhotonEnergy  , &b_PhotonEnergy);

    int counterCalib = 0;
    int nEventsCalib =  photons->GetEntries();
    // int feedbackDivision = (int) (nEventsCalib/20);
    // if(feedbackDivision == 0) feedbackDivision = 1;
    std::cout << "Filling calibration points array..." << std::endl;
    for(int iEvent = 0 ; iEvent < nEventsCalib ; iEvent++)
    {
      // if (  iEvent % feedbackDivision == 0 ) { std::cout << (int) ((double) iEvent/nEventsCalib * 100) << "% done..." << std::endl;}
      // std::cout << "in "<< iEvent << std::endl;
      photons->GetEvent(iEvent);
      // find index of x, y, z and energy, fast
      // knowing the starting paramters, it should be easy...
      // what is not easy is to cope with the rounding errors!
      int xIndex = (int) roundf((vertX - xmin)/xstep);
      int yIndex = (int) roundf((vertY - ymin)/ystep);
      int zIndex = (int) roundf((vertZ - zmin)/zstep);
      int eIndex = (int) roundf((PhotonEnergy - emin)/estep);
      int sIndex = 0;
      if(PositionZ >0)
      {
        sIndex = 1;
      }

      // if(i ==1) std::cout << &calibration_point[xIndex][yIndex][zIndex][eIndex][sIndex] << std::endl;
      // fill histo of approp calib point
      outputFile->cd();
      calibration_point[xIndex][yIndex][zIndex][eIndex][sIndex].delay->Fill(globalTime);
      // calibration_point[xIndex][yIndex][zIndex][eIndex][sIndex].globalTimes.push_back(globalTime);
    }
    // std::cout << "after" << std::endl;
    inputFile->Close();
  }
  // std::cout << "pre " << std::endl;
  // run on calibration_point to set the eff values
  for(int ix = 0 ; ix < x.size() ; ix++)
  {
    for(int iy = 0 ; iy < y.size() ; iy++)
    {
      for(int iz = 0 ; iz < z.size() ; iz++)
      {
        for(int ie = 0 ; ie < energy.size() ; ie++)
        {
          for(int iSide = 0; iSide < nSides ; iSide++)
          {
            // std::cout << &calibration_point[ix][iy][iz][ie][iSide] << std::endl;
            // create histo
            // std::stringstream ssname;
            // ssname << "delayHisto " << x[ix] << " " << y[iy] << " " << z[iz] << " " << energy[ie] << " " << iSide;
            // calibration_point[ix][iy][iz][ie][iSide].delay = new TH1F(ssname.str().c_str(),ssname.str().c_str(),pdfBins,pdfStart,pdfEnd);
            // fill histo
            // for(int i = 0; i < calibration_point[ix][iy][iz][ie][iSide].globalTimes.size();i++)
            // {
            //   calibration_point[ix][iy][iz][ie][iSide].delay->Fill(calibration_point[ix][iy][iz][ie][iSide].globalTimes[i]);
            // }
            // std::cout << "post " << std::endl;
            calibration_point[ix][iy][iz][ie][iSide].eff = ( (float) calibration_point[ix][iy][iz][ie][iSide].delay->GetEntries()) / ( (float) calibration_point[ix][iy][iz][ie][iSide].primaries);
          }
        }
      }
    }
  }
  std::cout << "\nFinished reading optical calibration files" << std::endl;
  //----------------------------------//


  // FILE *calibrationFile;
  // calibrationFile = fopen("calibration.data","wb");
  std::ofstream calibrationFile(outputFileName.c_str(), std::ios::binary);
  // calibrationFile.write((char*)&calibration_point, sizeof(calibration_point));

  // write an header
  // bins
  // pdfStart
  // pdfEnd
  // xn
  // yn
  // zn
  // en
  // x
  // y
  // z
  // energy
  // bin info, take from first calibration point
  int bins = calibration_point[0][0][0][0][0].delay->GetNbinsX();
  calibrationFile.write((char*)&bins, sizeof(int));
  calibrationFile.write((char*)&pdfStart, sizeof(float));
  calibrationFile.write((char*)&pdfEnd, sizeof(float));
  calibrationFile.write((char*)&xn, sizeof(int));
  calibrationFile.write((char*)&yn, sizeof(int));
  calibrationFile.write((char*)&zn, sizeof(int));
  calibrationFile.write((char*)&en, sizeof(int));

  for(int i = 0 ; i < x.size(); i++)
  {
    calibrationFile.write((char*)&x[i], sizeof(float));
  }
  for(int i = 0 ; i < y.size(); i++)
  {
    calibrationFile.write((char*)&y[i], sizeof(float));
  }
  for(int i = 0 ; i < z.size(); i++)
  {
    calibrationFile.write((char*)&z[i], sizeof(float));
  }
  for(int i = 0 ; i < energy.size(); i++)
  {
    calibrationFile.write((char*)&energy[i], sizeof(float));
  }

  // write in a binary file (in this case ROOT files can go to...)
  for(int ix = 0 ; ix < x.size() ; ix++)
  {
    for(int iy = 0 ; iy < y.size() ; iy++)
    {
      for(int iz = 0 ; iz < z.size() ; iz++)
      {
        for(int ie = 0 ; ie < energy.size() ; ie++)
        {
          for(int iSide = 0; iSide < nSides ; iSide++)
          {
            // calibration_point_output_t calibration_point_out;
            calibrationFile.write((char*)&ix, sizeof(int));
            calibrationFile.write((char*)&iy, sizeof(int));
            calibrationFile.write((char*)&iz, sizeof(int));
            calibrationFile.write((char*)&ie, sizeof(int));
            calibrationFile.write((char*)&iSide, sizeof(int));

            calibrationFile.write((char*)&calibration_point[ix][iy][iz][ie][iSide].x, sizeof(float));
            calibrationFile.write((char*)&calibration_point[ix][iy][iz][ie][iSide].y, sizeof(float));
            calibrationFile.write((char*)&calibration_point[ix][iy][iz][ie][iSide].z, sizeof(float));
            calibrationFile.write((char*)&calibration_point[ix][iy][iz][ie][iSide].s, sizeof(int));
            calibrationFile.write((char*)&calibration_point[ix][iy][iz][ie][iSide].energy, sizeof(float));
            calibrationFile.write((char*)&calibration_point[ix][iy][iz][ie][iSide].eff, sizeof(float));
            // std::cout << calibration_point[ix][iy][iz][ie][iSide].eff << std::endl;

            for(int b = 1 ; b < bins+1; b++)
            {
              float value = calibration_point[ix][iy][iz][ie][iSide].delay->GetBinContent(b);
              calibrationFile.write( (char*)&value, sizeof(float));
            }
          }
        }
      }
    }
  }
  // write a footer
  // it will consist only of the moduleCenter followed by the Data_t struct, if it's there
  // moduleCenter
  calibrationFile.write((char*)&moduleCenter, sizeof(float));
  // calibrationFile.write((char*)&moduleCenter, sizeof(float));
  if(parametersFound)
  {
    relevant_data_t relevant;
    fillRelevant(relevant,parameters);
    // std::cout << "cut " << relevant.defaultCut << std::endl;
    calibrationFile.write((char*)&relevant, sizeof(relevant));
  }
  
  calibrationFile.close();



  gDirectory->cd();
  TDirectory *histoDir = outputFile->mkdir("Hybrid_Calibration");
  histoDir->cd();
  for(int ix = 0 ; ix < x.size() ; ix++)
  {
    for(int iy = 0 ; iy < y.size() ; iy++)
    {
      for(int iz = 0 ; iz < z.size() ; iz++)
      {
        for(int ie = 0 ; ie < energy.size() ; ie++)
        {
          for(int iSide = 0; iSide < nSides ; iSide++)
          {
            calibration_point[ix][iy][iz][ie][iSide].delay->Write();
          }
        }
      }
    }
  }
  outputFile->Close();
  std::cout << "Done. Goodbye!" << std::endl;

  return 0;
}



void usage()
{
  std::cout << "\t\t" << "[-f | --folder]         <folder path>   " << std::endl
            << "\t\t" << "[-p | --prefix]         <input files prefix>   " << std::endl
            << "\t\t" << "[-o | --output]         <output file name>    " << std::endl
            << "\t\t" << "[--pdfBins]             <number of divisions for sampling of propagation time PDFs>   - default = 5000" << std::endl
            << "\t\t" << "[--pdfEnd]              <end in ns of propagation time PDFs (begin is always 0ns)>    - default = 10" << std::endl
            << "\t\t" << "                        N.B. the combination of pdfBins and pdfEnd determins the LSB of PDFs (default = 10ns/5000 = 2ps)" << std::endl
            << "\t\t" << std::endl;
}
