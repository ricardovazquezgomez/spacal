// compile with
// g++ -o ../build/readoutScan readoutScan.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore

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


// struct directory_t
// {
//   int energy;
//   std::string name;
//   std::vector<std::string> listInputFiles;
  
//   TH1F *enHistoCali;
//   TH1F *enHistoCali_front;
//   TH1F *enHistoCali_back;
//   TH1F *enHistoNoCali;
//   TH1F *enHistoNoCali_front;
//   TH1F *enHistoNoCali_back;

//   TH1F *enHistoAllCali;
//   TH1F *enHistoAllCali_front;
//   TH1F *enHistoAllCali_back;
//   TH1F *enHistoAllNoCali;
//   TH1F *enHistoAllNoCali_front;
//   TH1F *enHistoAllNoCali_back;

//   TH1F *tHistoFront;
//   TH1F *tHistoBack;
//   std::vector<float> counts;
// };


// void read_directory(const std::string& name, std::vector<std::string> &v)
// {
//   // #######################################
//   // ### list files in directory
//   // ### taken from
//   // ### http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html
//   DIR* dirp = opendir(name.c_str());
//   struct dirent * dp;
//   while ((dp = readdir(dirp)) != NULL) {
//     v.push_back(dp->d_name);
//   }
//   closedir(dirp);
// }

void usage()
{
  std::cout << "\t\t" << "[-f | --folder]  <input folder>   " << std::endl
            << "\t\t" << "[-i | --input]   <input file prefix>   " << std::endl
            << "\t\t" << "[--calibration]  <f/b calibration factor - default = 1. >  " << std::endl
            // << "\t\t" << "[--type]         <module type - default = 0 >  " << std::endl
            << "\t\t" << "[--modules]      <number of modules - default = 1 >  " << std::endl
            << "\t\t" << "[--last]         <last front cell number - default = 1 >  " << std::endl
            << "\t\t" << "[--verbose]      <verbosity level - default = 0>  " << std::endl
            << "\t\t" << std::endl;
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
  // std::string outputFileName = "";
  bool singleEvent = false;
  int verbose = 0;
  int moduleTypeNumber = 0;
  float calibration_factor = 1.;
  int last_front_cell_select = 1;
  int modules = 1;
  // int eventsNumber = 1;
  static struct option longOptions[] =
      {
          {"input", required_argument, 0, 0},
          {"folder", required_argument, 0, 0},
          {"verbose", required_argument, 0, 0},
          {"type", required_argument, 0, 0},
          {"calibration", required_argument, 0, 0},
          {"last", required_argument, 0, 0},
          {"modules", required_argument, 0, 0},
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
    else if (c == 0 && optionIndex == 5)
    {
      last_front_cell_select = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 6)
    {
      modules = atoi((char *)optarg);
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

  // float c_factor_pb = 1.53811;
  // float c_factor_w = 1.47525;
  // int last_front_cell_pb =  15;
  // int last_front_cell_w  =  63;

  // int wmods  = 32;
  // int pmods  = 144;
  // int totmods = wmods + pmods;
  TFile *inputFile = ((TFile *)0);

  for (unsigned int i = 0; i < listInputFiles.size(); i++)
  {
    if (inputFile)
      delete inputFile; // just a precaution
    // if(i == 0)
    // {
    //   // keep the first file...
    //   inputFile = firstFile;
    // }

    inputFile = new TFile(listInputFiles[i].c_str(), "READ");

    if (!inputFile || inputFile->IsZombie() || inputFile->TestBit(TFile::kRecovered))
    {
      std::cout << "Skipping file " << listInputFiles[i] << std::endl;
      continue;
    }
    std::cout << "Opening file " << listInputFiles[i] << std::endl;

    auto TT = (TTree*) inputFile->Get("tree");

    std::vector<int>   * VPh = 0;
    std::vector<int>   * VMod = 0;
    std::vector<float> * VTime = 0;
    int entry = 0;
    int totalLight = 0;
    // int chosenModule = 121;
    auto entries = TT->GetEntries();

    float counts = 0;
    float no_cali_counts_front = 0;
    float no_cali_counts_back = 0;
    float cali_counts_front = 0;
    float cali_counts_back = 0;
    float all_counts_cali = 0;
    float all_counts_nocali = 0;
    float all_no_cali_counts_front = 0;
    float all_no_cali_counts_back = 0;
    float all_cali_counts_front = 0;
    float all_cali_counts_back = 0;
    int last_front_cell;
    float c_factor;

    
    for (int iEvent = 0; iEvent < entries; iEvent++)
    {
      counts = 0;
      no_cali_counts_front = 0;
      no_cali_counts_back = 0;
      cali_counts_front = 0;
      cali_counts_back = 0;
      all_counts_cali = 0;
      all_counts_nocali = 0;
      all_no_cali_counts_front = 0;
      all_no_cali_counts_back = 0;
      all_cali_counts_front = 0;
      all_cali_counts_back = 0;
      for (int iMod = 0; iMod < modules; iMod++)
      {
        last_front_cell = last_front_cell_select;
        c_factor = calibration_factor;
        // if(iMod > pmods -1) 
        // {
        //   last_front_cell = last_front_cell_w;
        //   c_factor = c_factor_w;
        // }
        
        std::string moduleName = "mod";
        std::stringstream smod;
        smod << moduleName << iMod;
        moduleName = smod.str();

        // --- Assign branches' addresses
        TT->SetBranchAddress("Total_Light", &totalLight);
        TT->SetBranchAddress("modulesHit", &VMod);
        TT->SetBranchAddress((moduleName + "_ph").c_str(), &VPh);
        TT->SetBranchAddress((moduleName + "_t").c_str(), &VTime);

        TT->GetEntry(iEvent);

          
        for (size_t iCell = 0; iCell < VPh->size(); ++iCell)
        {
            
          // std::cout << VPh->at(iCell) << std::endl;
          if (iCell > last_front_cell)
          {
            // if(iMod == chosenModule)
            // {
            //   counts += c_factor * VPh->at(iCell);
            //   cali_counts_back += c_factor * VPh->at(iCell);
            //   no_cali_counts_back += VPh->at(iCell);
            // }
            
            all_counts_nocali       += VPh->at(iCell);
            all_no_cali_counts_back += VPh->at(iCell);
            all_counts_cali         += c_factor * VPh->at(iCell);
            all_cali_counts_back    += c_factor * VPh->at(iCell);
          }
          else
          {
            // if (iMod == chosenModule)
            // {
            //   counts += VPh->at(iCell);
            //   cali_counts_front += VPh->at(iCell);
            //   no_cali_counts_front += VPh->at(iCell);
            // }
            all_counts_cali += VPh->at(iCell);
            all_counts_nocali += VPh->at(iCell);
            all_no_cali_counts_front += VPh->at(iCell);
            all_cali_counts_front += VPh->at(iCell);
          }

          // if(iCell == 4) directories[iDir].tHisto4->Fill(V->at(iCell);)
        }
          
        // std::cout << VTime->at(front_cell) << " " << VTime->at(back_cell) << std::endl;
        // directories[iDir].tHistoFront->Fill(VTime->at(front_cell));
        // directories[iDir].tHistoBack->Fill(VTime->at(back_cell));
          
      }

      std::cout << iEvent << "\t" 
                << all_cali_counts_front << "\t" 
                << all_cali_counts_back << "\t" 
                << ((float) all_cali_counts_front) / (((float) all_cali_counts_front) + ((float) all_cali_counts_back)) << std::endl;
        
    }
      
    // f.Close();




  }

  // int front_cell = 10; //irrelevant in this version, we are just interested in energy
  // int back_cell = 26;  //irrelevant in this version, we are just interested in energy

  // std::ifstream fFile("calibration_factor.txt");
  // std::ifstream lFile("last_front_cell.txt");
  // std::ifstream cellFile("timing_cells.txt");

  // fFile >> c_factor;
  // lFile >> last_front_cell;
  // cellFile >> front_cell >> back_cell;

  // std::vector<directory_t> directories;
  // std::stringstream sname;
  // // create directory list
  // for(int i = 0 ; i < energies.size() ; i++)
  // {
  //   directory_t directory;

  //   sname << energies[i] << "GeV";

  //   directory.energy = energies[i];
  //   directory.name = sname.str();

  //   sname.str("");
  //   sname << "enHistoNoCali_" << energies[i] << "GeV";
  //   directory.enHistoNoCali = new TH1F(sname.str().c_str(),sname.str().c_str(),8e3,0,8e6);
    
  //   sname.str("");
  //   sname << "enHistoNoCali_front_" << energies[i] << "GeV";
  //   directory.enHistoNoCali_front = new TH1F(sname.str().c_str(),sname.str().c_str(),8e3,0,8e6);

  //   sname.str("");
  //   sname << "enHistoNoCali_back_" << energies[i] << "GeV";
  //   directory.enHistoNoCali_back = new TH1F(sname.str().c_str(), sname.str().c_str(), 8e3, 0, 8e6);


  //   sname.str("");
  //   sname << "enHistoCali_" << energies[i] << "GeV";
  //   directory.enHistoCali = new TH1F(sname.str().c_str(), sname.str().c_str(), 8e3, 0, 8e6);

  //   sname.str("");
  //   sname << "enHistoCali_front_" << energies[i] << "GeV";
  //   directory.enHistoCali_front = new TH1F(sname.str().c_str(),sname.str().c_str(),8e3,0,8e6);

  //   sname.str("");
  //   sname << "enHistoCali_back_" << energies[i] << "GeV";
  //   directory.enHistoCali_back = new TH1F(sname.str().c_str(), sname.str().c_str(), 8e3, 0, 8e6);



  //   sname.str("");
  //   sname << "enHistoAllCali_" << energies[i] << "GeV";
  //   directory.enHistoAllCali = new TH1F(sname.str().c_str(), sname.str().c_str(), 8e3, 0, 8e6);

  //   sname.str("");
  //   sname << "enHistoAllCali_front_" << energies[i] << "GeV";
  //   directory.enHistoAllCali_front = new TH1F(sname.str().c_str(), sname.str().c_str(), 8e3, 0, 8e6);

  //   sname.str("");
  //   sname << "enHistoAllCali_back_" << energies[i] << "GeV";
  //   directory.enHistoAllCali_back = new TH1F(sname.str().c_str(), sname.str().c_str(), 8e3, 0, 8e6);

  //   sname.str("");
  //   sname << "enHistoAllNoCali_" << energies[i] << "GeV";
  //   directory.enHistoAllNoCali = new TH1F(sname.str().c_str(), sname.str().c_str(), 8e3, 0, 8e6);

  //   sname.str("");
  //   sname << "enHistoAllNoCali_front_" << energies[i] << "GeV";
  //   directory.enHistoAllNoCali_front = new TH1F(sname.str().c_str(), sname.str().c_str(), 8e3, 0, 8e6);

  //   sname.str("");
  //   sname << "enHistoAllNoCali_back_" << energies[i] << "GeV";
  //   directory.enHistoAllNoCali_back = new TH1F(sname.str().c_str(), sname.str().c_str(), 8e3, 0, 8e6);



  //   sname.str("");
  //   sname << "tHistoFront_" << energies[i] << "GeV";
  //   directory.tHistoFront = new TH1F(sname.str().c_str(),sname.str().c_str(),2000,81,88);
    
  //   sname.str("");
  //   sname << "tHistoBack_" << energies[i] << "GeV";
  //   directory.tHistoBack = new TH1F(sname.str().c_str(),sname.str().c_str(),2000,81,88);

  //   // directory.counts = 0;

  //   directories.push_back(directory);
  //   sname.str("");
  // }

  // loop on dir, find files, put in list
  // for(int iDir = 0 ; iDir < directories.size(); iDir++)
  // {
  //   std::vector<std::string> v;
  //   std::string inputFolderName;
  //   inputFolderName = std::string("./") + directories[iDir].name + std::string("/");
  //   // inputFolderName = std::string("/eos/experiment/spacal/Simulations/byYanxi/spacal-simulation/PulseFormation2/") + directories[iDir].name + std::string("/");
  //   read_directory(inputFolderName, v);
  //   // extract files with correct prefix
  //   // std::vector<std::string> listInputFiles;
  //   for(unsigned int i = 0 ; i < v.size() ; i++)
  //   {
  //     if(!v[i].compare(0,inputFilePrefix.size(),inputFilePrefix))
  //     {
  //       directories[iDir].listInputFiles.push_back(inputFolderName + "/" + v[i]);
  //     }
  //   }

  // }

  // // check output
  // for(int iDir = 0 ; iDir < directories.size(); iDir++)
  // {
  //   std::cout << "Directory " << directories[iDir].name << std::endl;
  //   std::cout << "Energy    " << directories[iDir].energy << std::endl;
  //   std::cout << "Files:    " << std::endl;
  //   for(int i = 0 ; i < directories[iDir].listInputFiles.size(); i++)
  //   {
  //     std::cout << directories[iDir].listInputFiles[i] << std::endl;
  //   }
  // }

  // // loop on dir, loop on files, fill histos
  // for(int iDir = 0 ; iDir < directories.size(); iDir++)
  // {
  //   std::cout << "Directory " << directories[iDir].name << std::endl;
  //   for(int iFile = 0 ; iFile < directories[iDir].listInputFiles.size(); iFile++)
  //   {
      
  //     TFile f(directories[iDir].listInputFiles[iFile].c_str(), "READ");
  //     auto TT = (TTree*) f.Get("tree");
  //     // std::string moduleName = "mod0";
      

  //     // std::vector<std::vector<float>> *VVec=0;
  //     std::vector<int>   * VPh = 0;
  //     std::vector<int>   * VMod = 0;
  //     std::vector<float> * VTime = 0;
  //     int entry = 0;
  //     int totalLight = 0;
  //     int chosenModule = 121;

  //     auto entries = TT->GetEntries();

  //     float counts = 0;
  //     float no_cali_counts_front = 0;
  //     float no_cali_counts_back = 0;
  //     float cali_counts_front = 0;
  //     float cali_counts_back = 0;
  //     float all_counts_cali = 0;
  //     float all_counts_nocali = 0;
  //     float all_no_cali_counts_front = 0;
  //     float all_no_cali_counts_back = 0;
  //     float all_cali_counts_front = 0;
  //     float all_cali_counts_back = 0;
  //     int last_front_cell;
  //     float c_factor;

  //     for (int iEvent = 0; iEvent < 1; iEvent++)
  //     {
  //       for (int iMod = 0; iMod < totmods; iMod++)
  //       {
  //         last_front_cell = last_front_cell_pb;
  //         c_factor = c_factor_pb;
  //         if(iMod > pmods -1) 
  //         {
  //           last_front_cell = last_front_cell_w;
  //           c_factor = c_factor_w;
  //         }
          

  //         std::string moduleName = "mod";
  //         std::stringstream smod;
  //         smod << moduleName << iMod;
  //         moduleName = smod.str();

  //         // --- Assign branches' addresses
  //         TT->SetBranchAddress("Total_Light", &totalLight);
  //         TT->SetBranchAddress("modulesHit", &VMod);
  //         TT->SetBranchAddress((moduleName + "_ph").c_str(), &VPh);
  //         TT->SetBranchAddress((moduleName + "_t").c_str(), &VTime);

  //         TT->GetEntry(iEvent);

          
  //         for (size_t iCell = 0; iCell < VPh->size(); ++iCell)
  //         {
            
  //           // std::cout << VPh->at(iCell) << std::endl;
  //           if (iCell > last_front_cell)
  //           {
  //             if(iMod == chosenModule)
  //             {
  //               counts += c_factor * VPh->at(iCell);
  //               cali_counts_back += c_factor * VPh->at(iCell);
  //               no_cali_counts_back += VPh->at(iCell);
  //             }
  //             all_counts_cali += c_factor * VPh->at(iCell);
  //             all_counts_nocali += VPh->at(iCell);
  //             all_no_cali_counts_back += VPh->at(iCell);
  //             all_cali_counts_back   += c_factor * VPh->at(iCell);
  //           }
  //           else
  //           {
  //             if (iMod == chosenModule)
  //             {
  //               counts += VPh->at(iCell);
  //               cali_counts_front += VPh->at(iCell);
  //               no_cali_counts_front += VPh->at(iCell);
  //             }
  //             all_counts_cali += VPh->at(iCell);
  //             all_counts_nocali += VPh->at(iCell);
  //             all_no_cali_counts_front += VPh->at(iCell);
  //             all_cali_counts_front += VPh->at(iCell);
  //           }

  //           // if(iCell == 4) directories[iDir].tHisto4->Fill(V->at(iCell);)
  //         }
          
  //         // std::cout << VTime->at(front_cell) << " " << VTime->at(back_cell) << std::endl;
  //         // directories[iDir].tHistoFront->Fill(VTime->at(front_cell));
  //         // directories[iDir].tHistoBack->Fill(VTime->at(back_cell));
          
  //       }
  //       directories[iDir].enHistoNoCali->Fill(totalLight);
  //       directories[iDir].enHistoNoCali_front->Fill(no_cali_counts_front);
  //       directories[iDir].enHistoNoCali_back->Fill(no_cali_counts_back);
  //       directories[iDir].enHistoCali->Fill(counts);
  //       directories[iDir].enHistoCali_front->Fill(cali_counts_front);
  //       directories[iDir].enHistoCali_back->Fill(cali_counts_back);

  //       directories[iDir].enHistoAllCali->Fill(all_counts_cali);
  //       directories[iDir].enHistoAllCali_front->Fill(all_cali_counts_front);
  //       directories[iDir].enHistoAllCali_back->Fill(all_cali_counts_back);
  //       directories[iDir].enHistoAllNoCali->Fill(all_counts_nocali);
  //       directories[iDir].enHistoAllNoCali_front->Fill(all_no_cali_counts_front);
  //       directories[iDir].enHistoAllNoCali_back->Fill(all_no_cali_counts_back);

  //       directories[iDir].counts.push_back(counts);
  //     }
      
  //     f.Close();
  //   }
  //   std::cout << "Directory " << directories[iDir].name << std::endl;
  // }

  // // graph
  // for(int iDir = 0 ; iDir < directories.size(); iDir++)
  // {

  //   std::cout << "Gr Directory " << directories[iDir].name << std::endl;

  //   TF1  *gauss1 = new TF1("gauss1","gaus");
  //   directories[iDir].enHistoNoCali->Fit(gauss1);
  //   float min1 = gauss1->GetParameter(1) - 10.0*gauss1->GetParameter(2);
  //   float max1 = gauss1->GetParameter(1) + 10.0*gauss1->GetParameter(2);
  //   directories[iDir].enHistoNoCali->GetXaxis()->SetRangeUser(min1,max1);

  //   xEnNoCali.push_back( directories[iDir].energy );
  //   exEnNoCali.push_back(0.);
  //   yEnNoCali.push_back((directories[iDir].enHistoNoCali)->GetStdDev()/(directories[iDir].enHistoNoCali)->GetMean());
  //   eyEnNoCali.push_back( TMath::Sqrt(TMath::Power( (directories[iDir].enHistoNoCali)->GetStdDevError() / (directories[iDir].enHistoNoCali)->GetStdDev() ,2) + TMath::Power( (directories[iDir].enHistoNoCali)->GetMeanError() / (directories[iDir].enHistoNoCali)->GetMean(),2)  ) * (directories[iDir].enHistoNoCali)->GetStdDev()/(directories[iDir].enHistoNoCali)->GetMean());

  //   TF1  *gauss2 = new TF1("gauss2","gaus");
  //   directories[iDir].enHistoCali->Fit(gauss2);
  //   float min2 = gauss2->GetParameter(1) - 10.0*gauss2->GetParameter(2);
  //   float max2 = gauss2->GetParameter(1) + 10.0*gauss2->GetParameter(2);
  //   directories[iDir].enHistoCali->GetXaxis()->SetRangeUser(min2,max2);

  //   xEnCali.push_back( directories[iDir].energy );
  //   exEnCali.push_back(0.);
  //   yEnCali.push_back((directories[iDir].enHistoCali)->GetStdDev()/(directories[iDir].enHistoCali)->GetMean());
  //   eyEnCali.push_back( TMath::Sqrt(TMath::Power( (directories[iDir].enHistoCali)->GetStdDevError() / (directories[iDir].enHistoCali)->GetStdDev() ,2) + TMath::Power( (directories[iDir].enHistoCali)->GetMeanError() / (directories[iDir].enHistoCali)->GetMean(),2)  ) * (directories[iDir].enHistoCali)->GetStdDev()/(directories[iDir].enHistoCali)->GetMean());



  //   TF1  *gaussTfront = new TF1("gaussTfront","gaus");
  //   directories[iDir].tHistoFront->Fit(gaussTfront);
  //   float minTfront = gaussTfront->GetParameter(1) - 10.0*gaussTfront->GetParameter(2);
  //   float maxTfront = gaussTfront->GetParameter(1) + 10.0*gaussTfront->GetParameter(2);
  //   directories[iDir].tHistoFront->GetXaxis()->SetRangeUser(minTfront,maxTfront);

  //   xtFront.push_back( directories[iDir].energy );
  //   extFront.push_back( 0. );
  //   ytFront.push_back ( (directories[iDir].tHistoFront)->GetStdDev()      );
  //   eytFront.push_back( (directories[iDir].tHistoFront)->GetStdDevError() );

  //   TF1  *gaussTback = new TF1("gaussTback","gaus");
  //   directories[iDir].tHistoBack->Fit(gaussTback);
  //   float minTback = gaussTback->GetParameter(1) - 10.0*gaussTback->GetParameter(2);
  //   float maxTback = gaussTback->GetParameter(1) + 10.0*gaussTback->GetParameter(2);
  //   directories[iDir].tHistoBack->GetXaxis()->SetRangeUser(minTback,maxTback);

  //    xtBack.push_back( directories[iDir].energy );
  //   extBack.push_back( 0. );
  //    ytBack.push_back ( (directories[iDir].tHistoBack)->GetStdDev()      );
  //   eytBack.push_back( (directories[iDir].tHistoBack)->GetStdDevError() );

  //   // yt4.push_back(directories[iDir].tHisto4->GetRMS());
  //   // xt13.push_back( directories[iDir].energy );
  //   // yt13.push_back(directories[iDir].tHisto13->GetRMS());


  //   // MANUALLY CALC EN RES
  //   // calc mean
  //   double sum = std::accumulate(std::begin(directories[iDir].counts), std::end(directories[iDir].counts), 0.0);
  //   double m =  sum / directories[iDir].counts.size();
  //   //calc stdev
  //   double accum = 0.0;
  //   std::for_each (std::begin(directories[iDir].counts), std::end(directories[iDir].counts), [&](const double d) {
  //     accum += (d - m) * (d - m);
  //   });
  //   double stdev = sqrt(accum / (directories[iDir].counts.size()-1));
  //   // save enres
  //   xEn_manual.push_back( directories[iDir].energy );
  //   yEn_manual.push_back(stdev/m);


  // }

  // TGraph *grEnNoCali = new TGraphErrors(xEnNoCali.size(),&xEnNoCali[0],&yEnNoCali[0],&exEnNoCali[0],&eyEnNoCali[0]);
  // grEnNoCali-> SetName("enResNoCali");
  // grEnNoCali->SetTitle("enResNoCali");
  // grEnNoCali->GetXaxis()->SetTitle("Energy [GeV]");
  // grEnNoCali->GetYaxis()->SetTitle("Energy Resolution [sigma E / E]");

  // TGraph *grEnCali = new TGraphErrors(xEnCali.size(),&xEnCali[0],&yEnCali[0],&exEnCali[0],&eyEnCali[0]);
  // grEnCali-> SetName("enResCali");
  // grEnCali->SetTitle("enResCali");
  // grEnCali->GetXaxis()->SetTitle("Energy [GeV]");
  // grEnCali->GetYaxis()->SetTitle("Energy Resolution [sigma E / E]");

  // TGraph *grEn_manual = new TGraph(xEn_manual.size(),&xEn_manual[0],&yEn_manual[0]);
  // grEn_manual->SetName("enRes_manual");
  // grEn_manual->SetTitle("enRes_manual");
  // grEn_manual->GetXaxis()->SetTitle("Energy [GeV]");
  // grEn_manual->GetYaxis()->SetTitle("Energy Resolution [sigma E / E]");


  // TGraph *grtFront = new TGraph(xtFront.size(),&xtFront[0],&ytFront[0]);
  // grtFront->SetName("tResFront");
  // grtFront->SetTitle("tResFront");
  // grtFront->GetXaxis()->SetTitle("Energy [GeV]");
  // grtFront->GetYaxis()->SetTitle("Time Resolution sigma [ns]");
  // grtFront->SetLineColor(kBlue);
  // grtFront->SetMarkerColor(kBlue);
  // grtFront->SetMarkerStyle(20);
  // grtFront->SetMarkerSize(1.5);

  // TGraph *grtBack = new TGraph(xtBack.size(),&xtBack[0],&ytBack[0]);
  // grtBack->SetName("tResBack");
  // grtBack->SetTitle("tResBack");
  // grtBack->GetXaxis()->SetTitle("Energy [GeV]");
  // grtBack->GetYaxis()->SetTitle("Time Resolution sigma [ns]");
  // grtBack->SetLineColor(kRed);
  // grtBack->SetMarkerColor(kRed);
  // grtBack->SetMarkerStyle(21);
  // grtBack->SetMarkerSize(1.5);

  // TMultiGraph *mg = new TMultiGraph();
  // mg->Add(grtFront,"lp");
  // mg->Add(grtBack,"lp");
  // mg->GetXaxis()->SetTitle("Energy [GeV]");
  // mg->GetYaxis()->SetTitle("Time Resolution sigma [ns]");

  // TCanvas *c3 = new TCanvas("c3","c3",1200,800);
  // c3->SetLogx();
  // c3->SetGrid();
  // mg->Draw("a");
  // TLegend *legend = new TLegend(0.5,0.7,0.893,0.89,"");
  // // legend->SetFillStyle(1);
  // legend->AddEntry(grtFront,"Front cell","l");
  // legend->AddEntry(grtBack,"Back cell","l");
  // legend->Draw();
  // gPad->Modified();
  // mg->GetXaxis()->SetLimits(0.9,110);
  // mg->SetMinimum(0.);
  // mg->SetMaximum(0.3);

  // std::cout << "Saving..." << std::endl;
  // fOut->cd();
  // grEnNoCali->Write();
  // grEnCali->Write();
  // grEn_manual->Write();
  // mg->Write();
  // c3->Write();
  // grtFront->Write();
  // grtBack->Write();
  // for(int iDir = 0 ; iDir < directories.size(); iDir++)
  // {
  //   directories[iDir].enHistoNoCali->Write();
  //   directories[iDir].enHistoNoCali_front->Write();
  //   directories[iDir].enHistoNoCali_back->Write();
  //   directories[iDir].enHistoCali->Write();
  //   directories[iDir].enHistoCali_front->Write();
  //   directories[iDir].enHistoCali_back->Write();

  //   directories[iDir].enHistoAllNoCali->Write();
  //   directories[iDir].enHistoAllNoCali_front->Write();
  //   directories[iDir].enHistoAllNoCali_back->Write();
  //   directories[iDir].enHistoAllCali->Write();
  //   directories[iDir].enHistoAllCali_front->Write();
  //   directories[iDir].enHistoAllCali_back->Write();

  //   directories[iDir].tHistoFront->Write();
  //   directories[iDir].tHistoBack->Write();
  // }
  // fOut->Close();






  return 0;

}
