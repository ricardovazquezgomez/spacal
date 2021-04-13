// compile with
// g++ -o ../build/hybridScan hybridScan.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore

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

#include "ROOT/RDataFrame.hxx"

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <getopt.h>
#include <numeric>


#include <sys/stat.h>
#include <dirent.h>



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

void usage()
{
  std::cout << "\t\t" << "[-f | --folder]  <input folder>   " << std::endl
            << "\t\t" << "[-i | --input]   <input file prefix>   " << std::endl
            << "\t\t" << "[--calibration]  <f/b calibration factor - default = 1. >  " << std::endl
            << "\t\t" << "[--type]         <module type - default = 0 >  " << std::endl
            << "\t\t" << "[--verbose]      <verbosity level - default = 0>  " << std::endl
            << "\t\t" << std::endl;
}

int main (int argc, char **argv)
{
  ROOT::EnableImplicitMT(); // Enable ROOT's implicit multi-threading

  //------------------------------------------------//
  // PARSE INPUT                                    //
  //------------------------------------------------//
  // prepare variable, giving defaults
  std::string inputFilePrefix = "";
  std::string inputFolderName = "./";
  // std::string outputFileName = "";
  int verbose = 0;
  int moduleTypeNumber = 0;
  float calibration_factor = 1.;
  // int eventsNumber = 1;
  static struct option longOptions[] =
      {
          {"input", required_argument, 0, 0},
          {"folder", required_argument, 0, 0},
          {"verbose", required_argument, 0, 0},
          {"type", required_argument, 0, 0},
          {"calibration", required_argument, 0, 0},
          {NULL, 0, 0, 0}};
  // read input
  while (1)
  {
    int optionIndex = 0;
    int c = getopt_long(argc, argv, "i:f:", longOptions, &optionIndex);
    if (c == -1)
    {
      break;
    }
    if (c == 'i')
    {
      inputFilePrefix = (char *)optarg;
    }
    else if (c == 'f')
    {
      inputFolderName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 0)
    {
      inputFilePrefix = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 1)
    {
      inputFolderName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 2)
    {
      verbose = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 3)
    {
      moduleTypeNumber = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 4)
    {
      calibration_factor = atof((char *)optarg);
    }
    else
    {
      std::cout << "Usage: " << argv[0] << std::endl;
      usage();
      return 1;
    }
  }
  // check mandatory input
  if (inputFilePrefix == "")
  {
    std::cout << std::endl;
    std::cout << "ERROR! You need to provide an input file name!" << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }
  // if (outputFileName == "")
  // {
  //   std::cout << std::endl;
  //   std::cout << "ERROR! You need to provide an output file name!" << std::endl;
  //   std::cout << "See program usage below..." << std::endl;
  //   std::cout << std::endl;
  //   std::cout << argv[0];
  //   usage();
  //   return 1;
  // }

  //----------------------------------//
  // INPUT FILES                   //
  //----------------------------------//
  std::vector<std::string> v;
  read_directory(inputFolderName, v);
  // extract files with correct prefix
  std::vector<std::string> listInputFiles;
  for (unsigned int i = 0; i < v.size(); i++)
  {
    if (!v[i].compare(0, inputFilePrefix.size(), inputFilePrefix))
    {
      listInputFiles.push_back(inputFolderName + "/" + v[i]);
    }
  }
  // check if it's empty
  if (listInputFiles.size() == 0)
  {
    std::cout << std::endl;
    std::cout << "ERROR! Some input files do not exists! Aborting." << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }

  TFile *inputFile = ((TFile *)0);
  // TTree *shower = ((TTree *)0);

  for (unsigned int i = 0; i < listInputFiles.size(); i++)
  {
    if (inputFile)
      delete inputFile; // just a precaution
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

    
    ROOT::RDataFrame d("hybrid", inputFile);

    // LAMBDAS, filters etc 
    auto filterModuleType = [&moduleTypeNumber](int moduleType) { return moduleType == moduleTypeNumber; };
    auto filterFront = [](float z){return z < 0;};     // here we assume front z will always be < 0, and back > 0. Could be problematic for more general code
    auto filterBack = [](float z) {return z > 0;};     // here we assume front z will always be < 0, and back > 0. Could be problematic for more general code
    int eventNumber = 0;
    auto filterEvent = [&eventNumber](int event){return event == eventNumber;};

    // header
    std::cout << "event\t"
              << "front\t"
              << "back(corr)\t"
              << "f/(f+b)\t" 
              << std::endl;

    // loop on events
    auto maxEvent = d.Max("event");
    std::cout << "Total number of primary events = " << (*maxEvent)+1 << std::endl;
    for (int iEvent = 0; iEvent < (*maxEvent) + 1; iEvent++)
    {
      eventNumber = iEvent;
      auto Nfront = d.Filter(filterEvent, {"event"})
                        .Filter(filterModuleType, {"moduleType"}) // take only this module type
                        .Filter(filterFront, {"z"})               // count front photons
                        .Count();
      auto Nback = d.Filter(filterEvent, {"event"})
                       .Filter(filterModuleType, {"moduleType"}) // take only this module type
                       .Filter(filterBack, {"z"})                // count back photons
                       .Count();
      float f = (float)*Nfront;
      float b = (float)*Nback;
      float bCorr = b * calibration_factor;

      // output to terminal 
      std::cout << iEvent << "\t" << f << "\t" << bCorr << "\t" << f / (f + bCorr) << std::endl;
    }
    
  }

  return 0;

    
}