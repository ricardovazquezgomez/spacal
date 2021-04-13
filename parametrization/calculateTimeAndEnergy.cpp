// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch

// compile with
// g++ -o ../build/calculateTimeAndEnergy calculateTimeAndEnergy.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore

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
#include <algorithm>    // std::sort

#include <sys/stat.h>
#include <dirent.h>

void usage();

// list files in directory
// taken from
// http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html
void read_directory(const std::string& name, std::vector<std::string> &v)
{
    DIR* dirp = opendir(name.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        v.push_back(dp->d_name);
    }
    closedir(dirp);
}

int main(int argc, char** argv)
{

  // check args
  if(argc < 2)
  {
    std::cout << argv[0] << std::endl;
    usage();
    return 1;
  }

  gEnv->GetValue("TFile.Recover", 0);

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

  // variables
  std::string inputFolderName = "./";
  std::string inputFilePrefix = "";
  std::string outputFilePrefix = "plots_res";
  float en_value = 1.; // GeV
  float angle_x  = 0.; // Deg
  float angle_y  = 0.; // Deg

  // parse command line arguments
  static struct option longOptions[] =
  {
    { "input", required_argument, 0, 0 },
    { "output", required_argument, 0, 0 },
    { "folder", required_argument, 0, 0 },
    { "energy", required_argument, 0, 0 },
    { "angle_x", required_argument, 0, 0 },
    { "angle_y", required_argument, 0, 0 },
    { NULL, 0, 0, 0 }
  };
  while(1) {
    int optionIndex = 0;
    int c = getopt_long(argc, argv, "i:o:f:", longOptions, &optionIndex);
    if (c == -1) {
      break;
    }
    else if (c == 'i'){
      inputFilePrefix = (char *)optarg;
    }
    else if (c == 'o'){
      outputFilePrefix = (char *)optarg;
    }
    else if (c == 'f'){
      inputFolderName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 0){
      inputFilePrefix = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 1){
      outputFilePrefix = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 2){
      inputFolderName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 3){
      en_value = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 4){
      angle_x = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 5){
      angle_y = atof((char *)optarg);
    }
    else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
      usage();
      return 1;
    }
  }

  // check if input file prefix is given
  if(inputFilePrefix == "")
  {
    std::cout << std::endl;
    std::cout << "ERROR! You need to an input file prefix!" << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }

  // feedback
  std::stringstream feedbackString;
  feedbackString << std::endl;
  feedbackString << PWDstring << std::endl;
  feedbackString << CommandLineString << std::endl;
  feedbackString << std::endl;
  feedbackString << "|***************************************************************|" << std::endl;
  feedbackString << "| HYBRID MONTE CARLO FOR SPACAL-RD                              |" << std::endl;
  feedbackString << "| QUICK PLOTS                                                   |" << std::endl;
  feedbackString << "|***************************************************************|" << std::endl;
  feedbackString << "" << std::endl;
  feedbackString << "COMMAND LINE PARAMETERS:" << std::endl;
  feedbackString << "Input folder            = " << inputFolderName <<std::endl;
  feedbackString << "Input file prefix       = " << inputFilePrefix <<std::endl;
  feedbackString << "Output file prefix      = " << outputFilePrefix <<std::endl;
  feedbackString << "Energy [GeV]            = " << en_value <<std::endl;
  feedbackString << "Angle x [deg]           = " << angle_x <<std::endl;
  feedbackString << "Angle y [deg]           = " << angle_y <<std::endl;
  feedbackString << std::endl;

  std::cout << feedbackString.str() << std::endl;

  // output files
  std::string outputFileRootName = outputFilePrefix + ".root";
  std::string outputFileTextName = outputFilePrefix + ".txt";
  TFile *outputFileRoot = new TFile(outputFileRootName.c_str(),"RECREATE");
  std::ofstream outputFileText;
  outputFileText.open(outputFileTextName.c_str(),std::ifstream::out);
  outputFileText << feedbackString.str() << std::endl;

  // open files
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

  // -------------------------- //
  // TChain of input files      //
  // -------------------------- //
  TChain *tree =  new TChain("tree"); // read input files
  for(unsigned int i = 0 ; i < listInputFiles.size(); i++)
  {
    // add file only if ok...
    TFile *inputFile = ((TFile *)0);
    inputFile = new TFile(listInputFiles[i].c_str());
    if (!inputFile || inputFile->IsZombie() || inputFile->TestBit(TFile::kRecovered) || inputFile->GetNkeys() == 0)
    {
      std::cout << "Skipping file " << listInputFiles[i] << std::endl;
    }
    else
    {
      std::cout << "Adding file " << listInputFiles[i] << std::endl;
      tree->Add(listInputFiles[i].c_str());
    }
  }
  std::cout << "|----------------------------------------|" << std::endl;
  std::cout << std::endl;
  // -------------------------- //

  // do histograms
  ULong64_t max = 0;;
  for(int i = 0; i < 18 ; i++)
  {
    std::stringstream sname;
    sname << "ch" << i;
    max += tree->GetMaximum(sname.str().c_str());
  }
  TH1F *energyH  = new TH1F("energyH","Sum of all 18 cells",100,0,max);
  energyH->GetXaxis()->SetTitle("Total charge collected [A.U.]");
  energyH->GetYaxis()->SetTitle("Counts");
  TH1F *time4   = new TH1F("time4" ,"Time histogram Cell 4",1000,-5,5);
  time4->GetXaxis()->SetTitle("Time [ns]");
  time4->GetYaxis()->SetTitle("Counts");
  TH1F *time13  = new TH1F("time13","Time histogram Cell 13",1000,-5,5);
  time13->GetXaxis()->SetTitle("Time [ns]");
  time13->GetYaxis()->SetTitle("Counts");
  tree->Draw("ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15+ch16+ch17 >> energyH"," (ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15+ch16+ch17 ) > 0");
  tree->Draw("t4  >> time4","t4 > -1000");
  tree->Draw("t13 >> time13","t13 > -1000");
  //
  TH1F *u_time4   = new TH1F("u_time4" ,"Time histogram Cell 4",1000,-5,5);
  u_time4->GetXaxis()->SetTitle("Time [ns]");
  u_time4->GetYaxis()->SetTitle("Counts");
  TH1F *u_time13  = new TH1F("u_time13","Time histogram Cell 13",1000,-5,5);
  u_time13->GetXaxis()->SetTitle("Time [ns]");
  u_time13->GetYaxis()->SetTitle("Counts");
  tree->Draw("u_t4  >> u_time4","u_t4 > -1000");
  tree->Draw("u_t13 >> u_time13","u_t13 > -1000");
  //
  TH1F *gen_time4   = new TH1F("gen_time4" ,"Time histogram Cell 4",1000,-5,5);
  gen_time4->GetXaxis()->SetTitle("Time [ns]");
  gen_time4->GetYaxis()->SetTitle("Counts");
  TH1F *gen_time13  = new TH1F("gen_time13","Time histogram Cell 13",1000,-5,5);
  gen_time13->GetXaxis()->SetTitle("Time [ns]");
  gen_time13->GetYaxis()->SetTitle("Counts");
  tree->Draw("t_gen4  >> gen_time4","t_gen4 > -1000");
  tree->Draw("t_gen13 >> gen_time13","t_gen13 > -1000");

  // need also to fit to avoid weird events and the 0s issue (later removed from propagateHybrid)
  // energy
  float minFitEnergy = energyH->GetMean() - 2.0*energyH->GetRMS();
  float maxFitEnergy = energyH->GetMean() + 2.0*energyH->GetRMS();
  float maxValue     = energyH->GetXaxis()->GetBinCenter(energyH->GetMaximumBin());
  std::cout << minFitEnergy << " "<< maxFitEnergy << " " <<  maxValue << std::endl;
  TF1* energyGauss = new TF1("energyGauss","gaus",minFitEnergy,maxFitEnergy);
  energyGauss->SetParameter(0,maxValue); // norm
  energyGauss->SetParameter(1,energyH->GetMean()); // mean
  energyGauss->SetParameter(2,energyH->GetRMS()); // rms
  // t4
  float minFitt4       = time4->GetMean() - 2.0*time4->GetRMS();
  float maxFitt4       = time4->GetMean() + 2.0*time4->GetRMS();
  float maxValuet4     = time4->GetXaxis()->GetBinCenter(time4->GetMaximumBin());
  std::cout << minFitt4 << " "<< maxFitt4 << " " <<  maxValuet4 << std::endl;
  TF1* t4Gauss = new TF1("t4Gauss","gaus",minFitt4,maxFitt4);
  t4Gauss->SetParameter(0,maxValuet4); // norm
  t4Gauss->SetParameter(1,time4->GetMean()); // mean
  t4Gauss->SetParameter(2,time4->GetRMS()); // rms
  // t13
  float minFitt13       = time13->GetMean() - 2.0*time13->GetRMS();
  float maxFitt13       = time13->GetMean() + 2.0*time13->GetRMS();
  float maxValuet13     = time13->GetXaxis()->GetBinCenter(time13->GetMaximumBin());
  std::cout << minFitt13 << " "<< maxFitt13 << " " <<  maxValuet13 << std::endl;
  TF1* t13Gauss = new TF1("t13Gauss","gaus",minFitt13,maxFitt13);
  t13Gauss->SetParameter(0,maxValuet13); // norm
  t13Gauss->SetParameter(1,time13->GetMean()); // mean
  t13Gauss->SetParameter(2,time13->GetRMS()); // rms

  // fit
  energyH->Fit(energyGauss,"QR");
  time4->Fit(t4Gauss,"QR");
  time13->Fit(t13Gauss,"QR");
  // save to text file

  // then restric range og histograms according to fits
  energyH->GetXaxis()->SetRangeUser(energyGauss->GetParameter(1)-10.0*energyGauss->GetParameter(2),energyGauss->GetParameter(1)+10.0*energyGauss->GetParameter(2));
  // time4->GetXaxis()->SetRangeUser(t4Gauss->GetParameter(1)-5.0*t4Gauss->GetParameter(2),t4Gauss->GetParameter(1)+5.0*t4Gauss->GetParameter(2));
  // time13->GetXaxis()->SetRangeUser(t13Gauss->GetParameter(1)-5.0*t13Gauss->GetParameter(2),t13Gauss->GetParameter(1)+5.0*t13Gauss->GetParameter(2));

  outputFileText << en_value << " " << angle_x << " " << angle_y << " "
                 << energyH->GetMean() << " " << energyH->GetRMS() << " "
                 << time4->GetMean()   << " " << time4->GetRMS() << " "
                 << time13->GetMean()  << " " << time13->GetRMS() << " "
                 << energyGauss->GetParameter(1) <<  " " << energyGauss->GetParameter(2) << " "
                 << t4Gauss->GetParameter(1) <<  " " << t4Gauss->GetParameter(2) << " "
                 << t13Gauss->GetParameter(1) <<  " " << t13Gauss->GetParameter(2) << " "
                 << u_time4->GetMean()   << " " << u_time4->GetRMS() << " "
                 << u_time13->GetMean()  << " " << u_time13->GetRMS() << " "
                 << gen_time4->GetMean()   << " " << gen_time4->GetRMS() << " "
                 << gen_time13->GetMean()  << " " << gen_time13->GetRMS() << " "
                 << std::endl;

  // save to root file
  std::cout << std::endl;
  std::cout << "Saving data to file " << outputFileRootName << std::endl;
  outputFileRoot->cd();
  energyH->Write();
  time4->Write();
  time13->Write();
  outputFileRoot->Close();

  std::cout << "Done." << std::endl;
  return 0;
}





void usage()
{
  std::cout << "\t\t" << "[ [-f|--folder] <input folder>            path to input folder] " << std::endl
            << "\t\t" << "[ [-i|--prefix] <input prefix>            prefix of input files] " << std::endl
            << "\t\t" << "[ [-o|--output] <output files prefx>      output file prefix] " << std::endl
            << "\t\t" << std::endl;
}
