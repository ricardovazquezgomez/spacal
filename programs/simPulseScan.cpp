// compile with
// g++ -o ../build/simPulseScan simPulseScan.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore

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
#include "TEnv.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <getopt.h>
#include <numeric>

#include <iomanip>
#include <sys/stat.h>
#include <dirent.h>


// info about the command line usage
void usage()
{
  std::cout << std::endl;
  std::cout << "\t\t" << "[-f | --folder]  <input folder>   " << std::endl
            << "\t\t" << "[-i | --input]   <input file prefix>   " << std::endl
            << "\t\t" << "[-o | --output]  <output file name>   " << std::endl
            << "\t\t" << "[--energy]       <energy of incident beam [GeV] - mandatory >  " << std::endl
            << "\t\t" << "[--type]         <csv list of module types - default = 0 >  " << std::endl
            << "\t\t" << "[--module]       <csv list of module numbers in sim output world, for each module type - default = 0 >  " << std::endl
            << "\t\t" << "[--cells]        <csv list of number of cells per module type - default = 2 >  " << std::endl
            << "\t\t" << "[--lastFront]    <csv list of last front cell for each module type - default = 0 >  " << std::endl
            << "\t\t" << "[--fraction]     <csv list of sampling fractions per module type - default = 1>  " << std::endl
            << "\t\t" << "[--caliFront]    <csv list of front calibration Ph.El./MeV factors - default = 1>  " << std::endl
            << "\t\t" << "[--caliBack]     <csv list of back calibration Ph.El./MeV factors - default = 1>  " << std::endl
            << "\t\t" << "[--verbose]      <verbosity level - default = 0>  " << std::endl
            << "\t\t" << "[--fl]           <FL factor in pulse formation - default = 1>  " << std::endl
            << "\t\t" << "[--combine_front]           <front cell number for combined t histogram - default = -1>  " << std::endl
            << "\t\t" << "[--combine_back]            <back cell number for combined t histogram - default = -1>  " << std::endl

            << "\t\t" << std::endl;
}


void read_directory(const std::string& name, std::vector<std::string> &v)
{
  // #######################################
  // ### list files in directory
  // ### taken from
  // ### http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html
    DIR* dirp = opendir(name.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        v.push_back(dp->d_name);
    }
    closedir(dirp);
}


// struct directory_t
// {
//   int energy;
//   std::string name;
//   std::vector<std::string> listInputFiles;
//   TH1F *histoLightFront;
//   TH1F *enHistoCali;
//   TH1F *tHistoFront;
//   TH1F *tHistoBack;
//   std::vector<float> counts;
// };

struct module_t
{
  // parameters
  int number;
  int type;
  int last_front;
  int n_cells;
  float fb_calibration;
  float front_calibration;
  float back_calibration;
  float sampling_fraction;
  
  
  // data vectors
  std::vector<float>* v_l_cell;
  std::vector<float>* v_t_cell;
  std::vector<float>  v_combined_t;
  std::vector<float>  v_quick_combined_t;

  float *sigma;
  float *offset;

  
  
  // histograms
  TH1F** h_l_cell;
  TH1F** h_phel_mev_cell;
  TH1F** h_t_cell;
  TH1F* h_t_combined;
  
  TH1F* h_t_quick_combined;
  int combine_front;
  int combine_back;
  // TH1F* h_energy;
  // int sections;
  // int *nCellsPerSections;
  // TH1F *h_EnDepInCrystal   ;
  // TH1F *h_EnDepNotInCrystal;
  // TH1F *h_EnDepTotal       ;
  // bool **saveCell;
  // TH1F ***h_Cell_EnDepInCrystal   ;
  // cell_t **cell;
};

void compute_histogram_parameters(std::vector<float> values,float &max,float &min,int &bins_total)
{
  // fill histograms
  // first calc mean and stdevs
  float avg = 0. ;
  float std = 0. ;
  unsigned int n = values.size();
  if ( n != 0) {
     avg = accumulate( values.begin(), values.end(), 0.0) / n; 
     double sq_sum = std::inner_product(values.begin(), values.end(), values.begin(), 0.0);
     std = std::sqrt(sq_sum / values.size() - avg * avg);
  }
  // calculate appropriate bin size 
  // by assuming 100 bins in the +/- 4 stdevs range
  int bins = 100;
  float range = (avg + 2.0*std) - (avg - 2.0*std);
  float bin_size = range / bins;
  max = (avg + 10.0*std);
  min = (avg - 10.0*std);
  if(min < 0 ) min = 0;
  bins_total = (int) std::round((max-min) / bin_size);
}

int main (int argc, char **argv)
{

  // --- FAST DICTIONARIES COMPILATION ON THE FLY
  gInterpreter -> GenerateDictionary("std::vector<std::vector<float>>", "vector");
  gInterpreter -> GenerateDictionary("std::vector<std::vector<std::vector<float>>>", "vector");
  // gStyle->SetOptStat(0000);
  

  //------------------------------------------------//
  // PARSE INPUT                                    //
  //------------------------------------------------//
  // prepare variable, giving defaults
  std::string inputFilePrefix = "";
  std::string inputFolderName = "./";
  std::string outputFileName = "";
  int verbose = 0;
  bool multiType = false;
  bool multiFrac = false;
  bool multiCaliFront = false;
  bool multiCaliBack = false;
  bool multiMod  = false;
  bool multiLast = false;
  bool multiCell = false;
  std::vector<int> typeNumber;
  std::string typeString = "";
  std::vector<float> fracValue;
  std::string fracString = "";
  std::vector<float> caliFrontValue;
  std::string caliFrontString = "";
  std::vector<float> caliBackValue;
  std::string caliBackString = "";
  std::vector<float> modValue;
  std::string modString = "";
  std::vector<float> lastValue;
  std::string lastString = "";
  std::vector<float> cellValue;
  std::string cellString = "";
  float energyGeV = -1;
  float fl = 1.;
  int combine_front = -1;
  int combine_back = -1;
  bool combine = false;
  static struct option longOptions[] =
  {
			{ "input", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { "folder", required_argument, 0, 0 },
      { "verbose", required_argument, 0, 0 },
      { "type", required_argument, 0, 0 },
      { "fraction", required_argument, 0, 0 },
      { "caliFront", required_argument, 0, 0 },
      { "caliBack", required_argument, 0, 0 },
      { "energy", required_argument, 0, 0 },
      { "module", required_argument, 0, 0 },
      { "lastFront", required_argument, 0, 0 },
      { "cells", required_argument, 0, 0 },
      { "fl", required_argument, 0, 0 },
      { "combine_front", required_argument, 0, 0 },
      { "combine_back", required_argument, 0, 0 },
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
      caliFrontString = (char *)optarg;
      multiCaliFront = true;
    }
    else if (c == 0 && optionIndex == 7){
      caliBackString = (char *)optarg;
      multiCaliBack = true;
    }
    else if (c == 0 && optionIndex == 8){
      energyGeV = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 9){
      modString = (char *)optarg;
      multiMod = true;
    }
    else if (c == 0 && optionIndex == 10){
      lastString = (char *)optarg;
      multiLast = true;
    }
    else if (c == 0 && optionIndex == 11){
      cellString = (char *)optarg;
      multiCell = true;
    }
    else if (c == 0 && optionIndex == 12){
      fl = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 13){
      combine_front = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 14){
      combine_back  = atoi((char *)optarg);
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
  if(energyGeV == -1)
  {
    std::cout << std::endl;
    std::cout << "ERROR! You NEEd to the energy of incident beam!!!" << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }
  //--------------------------------------------------
  

  //------------------- 
  if(combine_back > -1)
  {
    if(combine_front > -1)
    {
      combine = true;
    }
  }
  //------------------- 

  //------------------------------------------------//
  // SAMPLING FRACTIONS AND F/B CALIBRATION         //
  //------------------------------------------------//
  // do some mandatory check...
  // the "multi" boolean variables need to be either all true or all false, so
  if(multiFrac && multiCaliBack && multiCaliFront && multiType && multiMod && multiLast && multiCell) // all true, fine
  {
    // nothing to do
  }
  else if ((!multiFrac) && (!multiCaliFront) && (!multiCaliBack) && (!multiType) && (!multiMod) && (!multiLast) && (!multiCell)) // all false, fine
  {
    // nothing to do


  }
  else 
  {
    std::cout << std::endl;
    std::cout << "ERROR! If you specify type, fractions, and calibrations at the same time. Otherwise, do not specify anything" << std::endl;
    std::cout << "multiType       = " << multiType << std::endl;
    std::cout << "multiFrac       = " << multiFrac << std::endl;
    std::cout << "multiCaliFront  = " << multiCaliFront << std::endl;
    std::cout << "multiCaliBack   = " << multiCaliBack << std::endl;
    std::cout << "multiMod        = " << multiMod  << std::endl;
    std::cout << "multiLast       = " << multiLast << std::endl;
    std::cout << "multiCell       = " << multiCell << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }

  // parse sampling types, fractions, and calibrations 
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
  if(multiCaliFront)
  {
    std::stringstream cali_stream(caliFrontString); //create string stream from the string
    while(cali_stream.good()) {
      std::string substr;
      getline(cali_stream, substr, ','); //get first string delimited by comma
      float t_number;
      std::istringstream ( substr ) >> t_number;
      caliFrontValue.push_back(t_number);
    }
  }
  if(multiCaliBack)
  {
    std::stringstream cali_stream(caliBackString); //create string stream from the string
    while(cali_stream.good()) {
      std::string substr;
      getline(cali_stream, substr, ','); //get first string delimited by comma
      float t_number;
      std::istringstream ( substr ) >> t_number;
      caliBackValue.push_back(t_number);
    }
  }
  if(multiMod)
  {
    std::stringstream mod_stream(modString); //create string stream from the string
    while(mod_stream.good()) {
      std::string substr;
      getline(mod_stream, substr, ','); //get first string delimited by comma
      int t_number;
      std::istringstream ( substr ) >> t_number;
      modValue.push_back(t_number);
    }
  }
  if(multiLast)
  {
    std::stringstream last_stream(lastString); //create string stream from the string
    while(last_stream.good()) {
      std::string substr;
      getline(last_stream, substr, ','); //get first string delimited by comma
      int t_number;
      std::istringstream ( substr ) >> t_number;
      lastValue.push_back(t_number);
    }
  }
  if(multiCell)
  {
    std::stringstream cell_stream(cellString); //create string stream from the string
    while(cell_stream.good()) {
      std::string substr;
      getline(cell_stream, substr, ','); //get first string delimited by comma
      int t_number;
      std::istringstream ( substr ) >> t_number;
      cellValue.push_back(t_number);
    }
  }
  // check that they are equal in number
  if(
     (typeNumber.size() != fracValue.size()) ||
     (typeNumber.size() != caliFrontValue.size()) ||
     (typeNumber.size() != caliBackValue.size()) ||
     (typeNumber.size() != modValue.size()) ||
     (typeNumber.size() != lastValue.size()) ||
     (typeNumber.size() != cellValue.size()) 
    )
  {
    std::cout << std::endl;
      std::cout << "ERROR! You need to specify an equal number of types, fractions, and calibrations!" << std::endl;
      std::cout << "Types        = " << typeNumber.size() << std::endl;
      std::cout << "Fractions    = " << fracValue.size()  << std::endl;
      std::cout << "Cali Front   = " << caliFrontValue.size()  << std::endl;
      std::cout << "Cali Back    = " << caliBackValue.size()  << std::endl;
      std::cout << "Modules      = " << modValue.size()  << std::endl;
      std::cout << "LastFront    = " << lastValue.size()  << std::endl;
      std::cout << "Cells        = " << cellValue.size()  << std::endl;
      std::cout << "See program usage below..." << std::endl;
      std::cout << std::endl;
      std::cout << argv[0];
      usage();
      return 1;
  }
  // put in a map, from crystal type to fraction
  std::map<int, float> fractionMap;
  std::map<int, float> calibrationMapFront;
  std::map<int, float> calibrationMapBack;
  std::map<int, int> moduleMap;
  std::map<int, int> lastFrontMap;
  std::map<int, int> cellMap;
  bool usingFractionMap = true;
  bool usingCalibrationMapFront = true;
  bool usingCalibrationMapBack = true;
  bool usingModuleMap = true;
  bool usingLastMap = true;
  bool usingCellMap = true;
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
  //
  if((typeNumber.size() == 0) && (caliFrontValue.size() == 0))
  {
    usingCalibrationMapFront = false;
  }
  else 
  {
    for(unsigned int i = 0 ; i < typeNumber.size(); i++)
    {
      calibrationMapFront.insert(std::make_pair(typeNumber[i],caliFrontValue[i]));
    }
  }
  //
  if((typeNumber.size() == 0) && (caliBackValue.size() == 0))
  {
    usingCalibrationMapBack = false;
  }
  else 
  {
    for(unsigned int i = 0 ; i < typeNumber.size(); i++)
    {
      calibrationMapBack.insert(std::make_pair(typeNumber[i],caliBackValue[i]));
    }
  }
  //
  if((typeNumber.size() == 0) && (modValue.size() == 0))
  {
    usingModuleMap = false;
  }
  else 
  {
    for(unsigned int i = 0 ; i < typeNumber.size(); i++)
    {
      moduleMap.insert(std::make_pair(typeNumber[i],modValue[i]));
    }
  }
  //
  if((typeNumber.size() == 0) && (lastValue.size() == 0))
  {
    usingLastMap = false;
  }
  else 
  {
    for(unsigned int i = 0 ; i < typeNumber.size(); i++)
    {
      lastFrontMap.insert(std::make_pair(typeNumber[i],lastValue[i]));
    }
  }
  //
  if((typeNumber.size() == 0) && (cellValue.size() == 0))
  {
    usingCellMap = false;
  }
  else 
  {
    for(unsigned int i = 0 ; i < typeNumber.size(); i++)
    {
      cellMap.insert(std::make_pair(typeNumber[i],cellValue[i]));
    }
  }
  //--------------------------------------------------


  
  //------------------------------------------------//
  // INPUT FEEDBACK                                 //
  //------------------------------------------------//
  // compose a feedback string for info that we will also save directly in the output file
  std::stringstream feedbackString;
  feedbackString << "|***************************************************************|" << std::endl;
  feedbackString << "| INPUT FEEDBACK                                                |" << std::endl;
  feedbackString << "|***************************************************************|" << std::endl;
  feedbackString << "Folder                                 = " << inputFolderName << std::endl;
  feedbackString << "File prefix                            = " << inputFilePrefix << std::endl;
  feedbackString << "Output file name                       = " << outputFileName << std::endl;
  feedbackString << "Beam energy [GeV]                      = " << energyGeV << std::endl;
  feedbackString << "verbose                                = " << verbose << std::endl;
  feedbackString << "multiType                              = " << multiType << std::endl;
  feedbackString << "multiFrac                              = " << multiFrac << std::endl;
  feedbackString << "multiCaliFront                         = " << multiCaliFront << std::endl;
  feedbackString << "multiCaliBack                          = " << multiCaliBack << std::endl;
  feedbackString << "multiMod                               = " << multiMod  << std::endl;
  feedbackString << "multiLast                              = " << multiLast << std::endl;
  feedbackString << "multiCell                              = " << multiCell << std::endl;
  feedbackString << "Types                                  = " << typeNumber.size() << std::endl;
  feedbackString << "Fractions                              = " << fracValue.size()  << std::endl;
  feedbackString << "Calibrations Front                     = " << caliFrontValue.size()  << std::endl;
  feedbackString << "Calibrations Back                      = " << caliBackValue.size()  << std::endl;
  feedbackString << "LastFront                              = " << lastValue.size()  << std::endl;
  feedbackString << "Cells                                  = " << cellValue.size()  << std::endl;
  feedbackString << "FL                                     = " << fl << std::endl;
  feedbackString << "Combine                                = " << combine << std::endl;
  feedbackString << "Combine front                          = " << combine_front << std::endl;
  feedbackString << "Combine back                           = " << combine_back << std::endl;
  feedbackString << std::setw(4) << "Type"     
                 << std::setw(20) << "Sampling Fraction" 
                 << std::setw(20) <<"f Calibration" 
                 << std::setw(20) <<"b Calibration" 
                 << std::setw(20) <<"Module" 
                 << std::setw(20) <<"Cells" 
                 << std::setw(20) <<"Last Front" 
                 << std::endl;
  for(unsigned int i = 0 ; i < typeNumber.size(); i++)
  {
    feedbackString << std::setw(4) << typeNumber[i]  
                   << std::setw(20) << fractionMap.at(typeNumber[i])
                   << std::setw(20) << calibrationMapFront.at(typeNumber[i])
                   << std::setw(20) << calibrationMapBack.at(typeNumber[i])
                   << std::setw(20) << moduleMap.at(typeNumber[i])
                   << std::setw(20) << cellMap.at(typeNumber[i])
                   << std::setw(20) << lastFrontMap.at(typeNumber[i])
                   << std::endl;
  }



  feedbackString << "|***************************************************************|" << std::endl;

  // print the feedback string to screen
  feedbackString << std::endl;
  std::cout << feedbackString.str() << std::endl;
  //--------------------------------------------------
  
  TFile *fOut = new TFile(outputFileName.c_str(),"RECREATE");


  //----------------------------------//
  // INPUT FILES                      //
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
  else 
  {
    std::cout << "Input files found." << std::endl;
  }
  //----------------------------------//

  

  //----------------------------------//
  // MODULES VECTOR                   //
  //----------------------------------//
  std::vector<module_t> module;
  for(unsigned int i = 0 ; i < typeNumber.size(); i++)
  {
    module_t t_module;
    t_module.number            = moduleMap.at(typeNumber[i]);
    t_module.type              = typeNumber[i];
    t_module.last_front        = lastFrontMap.at(typeNumber[i]);
    t_module.n_cells           = cellMap.at(typeNumber[i]);
    t_module.front_calibration = calibrationMapFront.at(typeNumber[i]);
    t_module.back_calibration  = calibrationMapBack.at(typeNumber[i]);
    t_module.sampling_fraction = fractionMap.at(typeNumber[i]);
    t_module.v_l_cell          = new std::vector<float>[cellMap.at(typeNumber[i])];
    t_module.v_t_cell          = new std::vector<float>[cellMap.at(typeNumber[i])];
    t_module.h_l_cell          = new TH1F*[cellMap.at(typeNumber[i])];
    t_module.h_phel_mev_cell   = new TH1F*[cellMap.at(typeNumber[i])];
    t_module.h_t_cell          = new TH1F*[cellMap.at(typeNumber[i])];
    // if(combine) t_module.h_t_combined      = new TH1F*;
    t_module.sigma             = new float[cellMap.at(typeNumber[i])];
    t_module.offset            = new float[cellMap.at(typeNumber[i])];

    t_module.combine_front = combine_front;
    t_module.combine_back  = combine_back;
    module.push_back(t_module);
  }

  //----------------------------------//
  // VECTORS ETC                      //
  //----------------------------------//
  
  int n_modules = module.size();
  

  std::vector<float> v_l_front;
  std::vector<float> v_l_back;
  std::vector<float> v_t_front;
  std::vector<float> v_t_back;

  std::vector<float> v_calibrated_light;
  std::vector<float> v_total_phel_mev;
  TH1F* h_energy;
  TH1F* h_total_phel_mev;
  
  std::cout << "Running on input files..." << std::endl;
  for(int iFile = 0 ; iFile < listInputFiles.size(); iFile++)
  {
    std::cout << "Opening " << listInputFiles[iFile].c_str() << std::endl;
    TFile f (listInputFiles[iFile].c_str(), "READ");
    auto TT = (TTree*) f.Get("tree");
    // std::string moduleName = "mod4";
    // std::vector<std::vector<float>> *VVec=0;
    
    std::vector<int>   * VMod = 0;
    std::vector<int>   ** VPh = 0;
    std::vector<float> ** VTime = 0;
    VPh = new std::vector<int>*[n_modules];
    VTime = new std::vector<float>*[n_modules];
    for (int iMod = 0; iMod < n_modules; iMod++)
    {
      VPh[iMod] = 0;
      VTime[iMod] = 0;
    }
    int entry = 0;
    int totalLight = 0;

    // --- Assign branches' addresses
    TT->SetBranchAddress("Total_Light", &totalLight);
    TT->SetBranchAddress("modulesHit", &VMod);
    for (int iMod = 0; iMod < n_modules; iMod++)
    {
      std::string moduleName = "mod";
      std::stringstream smod;
      smod << moduleName << module[iMod].number;
      moduleName = smod.str();
      TT->SetBranchAddress( (moduleName + "_ph").c_str(), &VPh[iMod]);
      TT->SetBranchAddress( (moduleName + "_t").c_str(), &VTime[iMod]);
    }
    
    auto entries = TT->GetEntries();
    for(int iEvent = 0 ; iEvent < entries; iEvent++)
    {
      TT->GetEntry(iEvent);
      // a counter of calibrated ph.el. values 
      float calibrated_light = 0;
      float total_phel_mev = 0;

      // loop on modules
      
      for (int iMod = 0; iMod < n_modules; iMod++)
      {
        double combine_t_stamp = 0;
        //loop on cells
        for (size_t iCell = 0; iCell < VPh[iMod]->size(); ++iCell)
        {
          total_phel_mev +=( VPh[iMod]->at(iCell))/(1000.*energyGeV);
          module[iMod].v_l_cell[iCell].push_back(VPh[iMod]->at(iCell));
          module[iMod].v_t_cell[iCell].push_back(VTime[iMod]->at(iCell));
          if(iCell > module[iMod].last_front)
          {
            calibrated_light += (fl * VPh[iMod]->at(iCell)) / (module[iMod].sampling_fraction * module[iMod].back_calibration);
          }
          else 
          {
            calibrated_light += (fl * VPh[iMod]->at(iCell)) / (module[iMod].sampling_fraction * module[iMod].front_calibration);
          }

          // quick combination
          if(iCell == module[iMod].combine_front)
          {
            combine_t_stamp += VTime[iMod]->at(iCell);
          } 
          if(iCell == module[iMod].combine_back)
          {
            combine_t_stamp += VTime[iMod]->at(iCell);
          }

        }
        combine_t_stamp = combine_t_stamp/2.0;
        module[iMod].v_quick_combined_t.push_back(combine_t_stamp);
      }
      
      // fill energy vector
      v_calibrated_light.push_back(calibrated_light);
      v_total_phel_mev.push_back(total_phel_mev);
    }
    f.Close();
  }
  //----------------------------------//



  
  //----------------------------------//
  // HISTOS ETC                      //
  //----------------------------------//
  float max,min;
  int bins_total;
  for (int iMod = 0; iMod < n_modules; iMod++)
  {
    for (int iCell = 0; iCell < module[iMod].n_cells; iCell++)
    {
      std::stringstream sname;
      
      // ph.el. absolute
      compute_histogram_parameters(module[iMod].v_l_cell[iCell],max,min,bins_total);
      // and fill
      sname << "Phel_mod_" << module[iMod].number << "_cell_" << iCell;
      module[iMod].h_l_cell[iCell] = new TH1F(sname.str().c_str(),"",bins_total,min,max);
      sname.str("");
      sname << "Ph.El. in module " << module[iMod].number << " cell " << iCell;
      module[iMod].h_l_cell[iCell]->SetTitle(sname.str().c_str());
      sname.str("");
      module[iMod].h_l_cell[iCell]->GetXaxis()->SetTitle("Ph.El.");
      module[iMod].h_l_cell[iCell]->GetYaxis()->SetTitle("Counts");
      for(unsigned int i = 0; i < module[iMod].v_l_cell[iCell].size(); i++)
      {
        module[iMod].h_l_cell[iCell]->Fill(module[iMod].v_l_cell[iCell].at(i));
      }

      // ph.el. / MeV
      // compute_histogram_parameters(module[iMod].v_l_cell[iCell],max,min,bins_total);
      // and fill
      sname << "Phel_mev_mod_" << module[iMod].number << "_cell_" << iCell;
      module[iMod].h_phel_mev_cell[iCell] = new TH1F(sname.str().c_str(),"",bins_total,min/(1000*energyGeV),max/(1000.*energyGeV));
      sname.str("");
      sname << "Ph.El. per MeV in module " << module[iMod].number << " cell " << iCell;
      module[iMod].h_phel_mev_cell[iCell]->SetTitle(sname.str().c_str());
      sname.str("");
      module[iMod].h_phel_mev_cell[iCell]->GetXaxis()->SetTitle("Ph.El. / MeV");
      module[iMod].h_phel_mev_cell[iCell]->GetYaxis()->SetTitle("Counts");
      for(unsigned int i = 0; i < module[iMod].v_l_cell[iCell].size(); i++)
      {
        module[iMod].h_phel_mev_cell[iCell]->Fill((module[iMod].v_l_cell[iCell].at(i))/(1000.*energyGeV));
      }

      // t
      // TH1F *histo_t_front = new TH1F("histo_t_front","Timing front section",1000,41,48);
      sname << "t_mod_" << module[iMod].number << "_cell_" << iCell;
      module[iMod].h_t_cell[iCell] = new TH1F(sname.str().c_str(),"",1000,41,48);
      sname.str("");
      sname << "Timestamps in module " << module[iMod].number << " cell " << iCell;
      module[iMod].h_t_cell[iCell]->SetTitle(sname.str().c_str());
      sname.str("");

      module[iMod].h_t_cell[iCell]->GetXaxis()->SetTitle("t [ns]");
      module[iMod].h_t_cell[iCell]->GetYaxis()->SetTitle("Counts");
      for(unsigned int i = 0; i < module[iMod].v_t_cell[iCell].size(); i++)
      {
        module[iMod].h_t_cell[iCell]->Fill(module[iMod].v_t_cell[iCell].at(i));
      }

      // extract sigmas and offset 
      // sigma
      module[iMod].sigma[iCell] = module[iMod].h_t_cell[iCell]->GetStdDev();
      // offset, temporary
      module[iMod].offset[iCell] = module[iMod].h_t_cell[iCell]->GetMean();
    }

    // find mean offset
    int moduleOffsets = 0;
    double meanOffset = 0.;
    for (int iCell = 0; iCell < module[iMod].n_cells; iCell++)
    {
      meanOffset += module[iMod].offset[iCell];
      moduleOffsets++;
    }
    meanOffset = meanOffset / moduleOffsets;
    // rescale all offsets to mean
    for (int iCell = 0; iCell < module[iMod].n_cells; iCell++)
    {
      module[iMod].offset[iCell] = module[iMod].offset[iCell] - meanOffset;
    } 
  }


  //----------------------------------//
  // FRONT BACK COMBINATION           //
  //----------------------------------//
  for(int iFile = 0 ; iFile < listInputFiles.size(); iFile++)
  {
    TFile f (listInputFiles[iFile].c_str(), "READ");
    auto TT = (TTree*) f.Get("tree");
    // std::string moduleName = "mod4";
    // std::vector<std::vector<float>> *VVec=0;
    
    std::vector<int>   * VMod = 0;
    std::vector<int>   ** VPh = 0;
    std::vector<float> ** VTime = 0;
    VPh = new std::vector<int>*[n_modules];
    VTime = new std::vector<float>*[n_modules];
    for (int iMod = 0; iMod < n_modules; iMod++)
    {
      VPh[iMod] = 0;
      VTime[iMod] = 0;
    }
    int entry = 0;
    int totalLight = 0;

    // --- Assign branches' addresses
    TT->SetBranchAddress("Total_Light", &totalLight);
    TT->SetBranchAddress("modulesHit", &VMod);
    for (int iMod = 0; iMod < n_modules; iMod++)
    {
      std::string moduleName = "mod";
      std::stringstream smod;
      smod << moduleName << module[iMod].number;
      moduleName = smod.str();
      TT->SetBranchAddress( (moduleName + "_ph").c_str(), &VPh[iMod]);
      TT->SetBranchAddress( (moduleName + "_t").c_str(), &VTime[iMod]);
    }
    
    auto entries = TT->GetEntries();
    for(int iEvent = 0 ; iEvent < entries; iEvent++)
    {
      TT->GetEntry(iEvent);
      // a counter of calibrated ph.el. values 
      float calibrated_light = 0;
      float total_phel_mev = 0;

      // loop on modules
      for (int iMod = 0; iMod < n_modules; iMod++)
      {
        // calculate combined time stamp per module 
        float sumWeights = 0.;
        float timestamp = 0.;
        //loop on cells
        for (size_t iCell = 0; iCell < VPh[iMod]->size(); ++iCell)
        {
          float weigth = 1.0/(pow(module[iMod].sigma[iCell],2));
          sumWeights += weigth;
          timestamp += (VTime[iMod]->at(iCell) - module[iMod].offset[iCell]) * weigth;
        }
        timestamp = timestamp / sumWeights;
        module[iMod].v_combined_t.push_back(timestamp);
      }
    }
    f.Close();
  }
  //----------------------------------//


  // histograms of combined t 
  for (int iMod = 0; iMod < n_modules; iMod++)
  {
    
    std::stringstream sname;
    // t
    // TH1F *histo_t_front = new TH1F("histo_t_front","Timing front section",1000,41,48);
    sname << "t_mod_" << module[iMod].number << "_combined";
    module[iMod].h_t_combined = new TH1F(sname.str().c_str(),"",1000,41,48);
    sname.str("");
    sname << "Combined timestamps in module " << module[iMod].number;
    module[iMod].h_t_combined->SetTitle(sname.str().c_str());
    sname.str("");
    module[iMod].h_t_combined->GetXaxis()->SetTitle("t [ns]");
    module[iMod].h_t_combined->GetYaxis()->SetTitle("Counts");
    for(unsigned int i = 0; i < module[iMod].v_combined_t.size(); i++)
    {
      module[iMod].h_t_combined->Fill(module[iMod].v_combined_t.at(i));
    }
  }

  // quick histograms of combined t for just one cell
  for (int iMod = 0; iMod < n_modules; iMod++)
  {
    
    std::stringstream sname;
    // t
    // TH1F *histo_t_front = new TH1F("histo_t_front","Timing front section",1000,41,48);
    sname << "t_mod_" << module[iMod].number << "_quick_combined";
    module[iMod].h_t_quick_combined = new TH1F(sname.str().c_str(),"",1000,41,48);
    sname.str("");
    sname << "Quick combined timestamps in module " << module[iMod].number;
    module[iMod].h_t_quick_combined->SetTitle(sname.str().c_str());
    sname.str("");
    module[iMod].h_t_quick_combined->GetXaxis()->SetTitle("t [ns]");
    module[iMod].h_t_quick_combined->GetYaxis()->SetTitle("Counts");
    for(unsigned int i = 0; i < module[iMod].v_quick_combined_t.size(); i++)
    {
      module[iMod].h_t_quick_combined->Fill(module[iMod].v_quick_combined_t.at(i));
    }
  }


  // energy 
  // TH1F *histo_t_front = new TH1F("histo_t_front","Timing front section",1000,41,48);
  compute_histogram_parameters(v_calibrated_light,max,min,bins_total);
  // and fill
  std::stringstream sname;
  sname << "Calibrated_Energy";
  h_energy = new TH1F(sname.str().c_str(),"",bins_total,min,max);
  sname.str("");
  sname << "Calibrated Energy histogram";
  h_energy->SetTitle(sname.str().c_str());
  sname.str("");
  h_energy->GetXaxis()->SetTitle("A.U.");
  h_energy->GetYaxis()->SetTitle("Counts");
  for(unsigned int i = 0; i < v_calibrated_light.size(); i++)
  {
    h_energy->Fill(v_calibrated_light.at(i));
  }

  
  // photoelectrons per mev total 
  // TH1F *histo_t_front = new TH1F("histo_t_front","Timing front section",1000,41,48);
  compute_histogram_parameters(v_total_phel_mev,max,min,bins_total);
  // and fill
  // std::stringstream sname;
  sname << "Phel_mev";
  h_total_phel_mev = new TH1F(sname.str().c_str(),"",bins_total,min,max);
  sname.str("");
  sname << "Ph.el. per MeV entire prototype";
  h_total_phel_mev->SetTitle(sname.str().c_str());
  sname.str("");
  h_total_phel_mev->GetXaxis()->SetTitle("Ph.el. / MeV");
  h_total_phel_mev->GetYaxis()->SetTitle("Counts");
  for(unsigned int i = 0; i < v_total_phel_mev.size(); i++)
  {
    h_total_phel_mev->Fill(v_total_phel_mev.at(i));
  }

  
  fOut->cd();
  h_energy->Write();
  h_total_phel_mev->Write();
  for (int iMod = 0; iMod < n_modules; iMod++)
  {
    for (int iCell = 0; iCell < module[iMod].n_cells; iCell++)
    {
      module[iMod].h_l_cell[iCell]->Write();
    }
  }
  for (int iMod = 0; iMod < n_modules; iMod++)
  {
    for (int iCell = 0; iCell < module[iMod].n_cells; iCell++)
    {
      module[iMod].h_t_cell[iCell]->Write();
    }
    module[iMod].h_t_combined->Write();
    module[iMod].h_t_quick_combined->Write();
  }
  for (int iMod = 0; iMod < n_modules; iMod++)
  {
    for (int iCell = 0; iCell < module[iMod].n_cells; iCell++)
    {
      module[iMod].h_phel_mev_cell[iCell]->Write();
    }
  }
  // histo_light_back->Write();
  // histo_t_front->Write();
  // histo_t_back->Write();
  fOut->Close();
  
  return 0;

}

  