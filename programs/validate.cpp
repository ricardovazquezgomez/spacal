// compile with
// g++ -o ../build/validate validate.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore

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

#include <sys/stat.h>
#include <dirent.h>


Double_t doComparison(TCanvas* canvas, TH1F *ray, TH1F* com, TString Title, TString Title_ray, TString Title_comp,bool rise);
void usage();

// from https://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
inline bool exists_test (const std::string& name)
{
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

int main(int argc, char** argv)
{

  if(argc < 3)
  {
    std::cout << argv[0];
    usage();
    return 1;
  }

  gStyle->SetOptStat(0);
  gEnv->GetValue("TFile.Recover", 0);
  std::string out_prefix    = argv[1];
  std::string hybrid_prefix = argv[2];
  int maxEvents     = atoi(argv[3]);

  // feedback
  std::cout << "Prefix of full sim files    = " << out_prefix << std::endl;
  std::cout << "Prefix of hybrid files      = " << hybrid_prefix << std::endl;
  std::cout << "Max number of events        = " << maxEvents << std::endl;

  int binShort = 50;
  int binLong = 2000;

  TH1F* ksHistoShort0 = new TH1F("ksHistoShort0","Distribution of KS test results - Detector 0, first 2ns",100,0,1);
  TH1F* ksHistoShort1 = new TH1F("ksHistoShort1","Distribution of KS test results - Detector 1, first 2ns",100,0,1);
  TH1F* ksHistoLong0  = new TH1F("ksHistoLong0" ,"Distribution of KS test results - Detector 0, full pulse [0 - 2microsec]" ,100,0,1);
  TH1F* ksHistoLong1  = new TH1F("ksHistoLong1" ,"Distribution of KS test results - Detector 1, full pulse [0 - 2microsec]" ,100,0,1);

  // try opening all files with input prefix and numbers from 0 to maxEvents
  TFile *outputFile = new TFile("testOut.root","RECREATE");
  for(int i = 0 ; i < maxEvents ; i++)
  {
    std::stringstream sname;
    std::stringstream svar;
    std::stringstream sdet;
    std::stringstream scanvas;

    sname << out_prefix  << i << ".root";
    std::string outFileName = sname.str().c_str();
    sname.str("");
    sname << hybrid_prefix  << i << ".root";
    std::string hybFileName = sname.str().c_str();
    sname.str("");


    if(exists_test(outFileName) && exists_test(hybFileName))
    {
      std::cout << "Running on "<< outFileName << " and " << hybFileName << std::endl;
      // then open each file
      // ROOT::RDataFrame outDataFrame("photons", outFileName);
      // ROOT::RDataFrame outDataFrame("hybrid" , hybFileName);
      TFile *outFile = new TFile(outFileName.c_str());
      TFile *hybFile = new TFile(hybFileName.c_str());
      if (outFile->IsZombie() || hybFile->IsZombie())
      {
        //something very wrong, cannot use this file, exit
        std::cout << "Corrupted file(s), skipping "<< outFileName << " and " << hybFileName << std::endl;
      }
      else
      {
        if(outFile->TestBit(TFile::kRecovered) || hybFile->TestBit(TFile::kRecovered))
        {
          //the Recover procedure has been run when opening the file
          //and the Recover was successful
          std::cout << "WARNING, at least one file was recovered out of "<< outFileName << " and " << hybFileName << std::endl;
        }

        TTree *photons = (TTree*) outFile->Get("photons");
        TTree *hybrid = (TTree*)  hybFile->Get("hybrid");

        outputFile->cd();

        // ------> full ray tracing histograms
        // short histograms, det 0
        sname << "full_short_Det0_" << i;
        svar  << "globalTime >> " << sname.str().c_str();
        TH1F* fullShortDet0 = new TH1F(sname.str().c_str(),sname.str().c_str(),binShort,0,2);
        photons->Draw(svar.str().c_str(),"PositionZ < 0");
        sname.str("");
        svar.str("");
        // short histograms, det 1
        sname << "full_short_Det1_" << i;
        svar  << "globalTime >> " << sname.str().c_str();
        TH1F* fullShortDet1 = new TH1F(sname.str().c_str(),sname.str().c_str(),binShort,0,2);
        photons->Draw(svar.str().c_str(),"PositionZ > 0");
        sname.str("");
        svar.str("");
        // long histograms, det 0
        sname << "full_long_Det0_" << i;
        svar  << "globalTime >> " << sname.str().c_str();
        TH1F* fullLongDet0 = new TH1F(sname.str().c_str(),sname.str().c_str(),binLong,0,2000);
        photons->Draw(svar.str().c_str(),"PositionZ < 0");
        sname.str("");
        svar.str("");
        // long histograms, det 1
        sname << "full_long_Det1_" << i;
        svar  << "globalTime >> " << sname.str().c_str();
        TH1F* fullLongDet1 = new TH1F(sname.str().c_str(),sname.str().c_str(),binLong,0,2000);
        photons->Draw(svar.str().c_str(),"PositionZ > 0");
        sname.str("");
        svar.str("");

        // ------> hybrid histograms
        // short histograms, det 0
        sname << "hybr_short_Det0_" << i;
        svar  << "timestamp >> " << sname.str().c_str();
        TH1F* hybrShortDet0 = new TH1F(sname.str().c_str(),sname.str().c_str(),binShort,0,2);
        hybrid->Draw(svar.str().c_str(),"z < 0");
        sname.str("");
        svar.str("");
        // short histograms, det 1
        sname << "hybr_short_Det1_" << i;
        svar  << "timestamp >> " << sname.str().c_str();
        TH1F* hybrShortDet1 = new TH1F(sname.str().c_str(),sname.str().c_str(),binShort,0,2);
        hybrid->Draw(svar.str().c_str(),"z > 0");
        sname.str("");
        svar.str("");
        // long histograms, det 0
        sname << "hybr_long_Det0_" << i;
        svar  << "timestamp >> " << sname.str().c_str();
        TH1F* hybrLongDet0 = new TH1F(sname.str().c_str(),sname.str().c_str(),binLong,0,2000);
        hybrid->Draw(svar.str().c_str(),"z < 0");
        sname.str("");
        svar.str("");
        // long histograms, det 1
        sname << "hybr_long_Det1_" << i;
        svar  << "timestamp >> " << sname.str().c_str();
        TH1F* hybrLongDet1 = new TH1F(sname.str().c_str(),sname.str().c_str(),binLong,0,2000);
        hybrid->Draw(svar.str().c_str(),"z > 0");
        sname.str("");
        svar.str("");

        outputFile->cd();
        // do the comparison and the test
        // ---> short pulses
        sdet  << "Detector 0 - [0 - 2ns] - Event " << i;
        scanvas << "comparison_0 - [0 - 2ns] - Event " << i;
        TCanvas *comparison_s_0 = new TCanvas(scanvas.str().c_str(),scanvas.str().c_str(),1200,800);
        Double_t s_0 = doComparison(comparison_s_0,fullShortDet0,hybrShortDet0,sdet.str().c_str(), "Full Raytracing", "Hybrid MC",true);
        ksHistoShort0->Fill(s_0);
        comparison_s_0->Write();
        sdet.str("");
        scanvas.str("");

        sdet  << "Detector 1 - [0 - 2ns] - Event " << i;
        scanvas << "comparison_1 - [0 - 2ns] - Event " << i;
        TCanvas *comparison_s_1 = new TCanvas(scanvas.str().c_str(),scanvas.str().c_str(),1200,800);
        Double_t s_1 = doComparison(comparison_s_1,fullShortDet1,hybrShortDet1,sdet.str().c_str(), "Full Raytracing", "Hybrid MC",true);
        ksHistoShort1->Fill(s_1);
        comparison_s_1->Write();
        sdet.str("");
        scanvas.str("");

        // ---> long pulses
        sdet  << "Detector 0 - [0 - 2 microsec] - Event " << i;
        scanvas << "comparison_0 - [0 - 2 microsec] - Event " << i;
        TCanvas *comparison_l_0 = new TCanvas(scanvas.str().c_str(),scanvas.str().c_str(),1200,800);
        Double_t l_0 = doComparison(comparison_l_0,fullLongDet0,hybrLongDet0,sdet.str().c_str(), "Full Raytracing", "Hybrid MC",false);
        ksHistoLong0->Fill(l_0);
        comparison_l_0->Write();
        sdet.str("");
        scanvas.str("");

        sdet  << "Detector 1 - [0 - 2 microsec] - Event " << i;
        scanvas << "comparison_1 - [0 - 2 microsec] - Event " << i;
        TCanvas *comparison_l_1 = new TCanvas(scanvas.str().c_str(),scanvas.str().c_str(),1200,800);
        Double_t l_1 = doComparison(comparison_l_1,fullLongDet1,hybrLongDet1,sdet.str().c_str(), "Full Raytracing", "Hybrid MC",false);
        ksHistoLong1->Fill(l_1);
        comparison_l_1->Write();
        sdet.str("");
        scanvas.str("");




        outputFile->cd();
        fullShortDet0->Write();
        hybrShortDet0->Write();
        fullShortDet1->Write();
        hybrShortDet1->Write();
        fullLongDet0->Write();
        hybrLongDet0->Write();
        fullLongDet1->Write();
        hybrLongDet1->Write();

        outFile->Close();
        hybFile->Close();

      }



    }



  }


  outputFile->cd();
  TCanvas *c_ksHistoShort0 = new TCanvas("c_ksHistoShort0","c_ksHistoShort0",1200,800);
  ksHistoShort0->Draw();

  TCanvas *c_ksHistoShort1 = new TCanvas("c_ksHistoShort1","c_ksHistoShort1",1200,800);
  ksHistoShort1->Draw();

  TCanvas *c_ksHistoLong0 = new TCanvas("c_ksHistoLong0","c_ksHistoLong0",1200,800);
  ksHistoLong0->Draw();

  TCanvas *c_ksHistoLong1 = new TCanvas("c_ksHistoLong1","c_ksHistoLong1",1200,800);
  ksHistoLong1->Draw();

  c_ksHistoShort0->Write();
  c_ksHistoShort1->Write();
  c_ksHistoLong0->Write();
  c_ksHistoLong1->Write();
  ksHistoShort0->Write();
  ksHistoShort1->Write();
  ksHistoLong0->Write();
  ksHistoLong1->Write();
  outputFile->Close();
  return 0;
}

Double_t doComparison(TCanvas* canvas, TH1F *ray, TH1F* com, TString Title, TString Title_ray, TString Title_comp, bool rise)
{

  // clone histos
  TH1F* ray_0 = (TH1F*) ray->Clone();
  TH1F* com_0 = (TH1F*) com->Clone();
  ray_0->GetXaxis()->SetTitle("time from primary emission [ns]");
  ray_0->GetYaxis()->SetTitle("Counts");
  com_0->GetXaxis()->SetTitle("time from primary emission [ns]");
  com_0->GetYaxis()->SetTitle("Counts");
  canvas->cd();
  canvas->SetGrid();
  if(!rise) canvas->SetLogy();
  ray_0->Draw("HIST");
  ray_0->SetLineWidth(3);
  ray_0->SetLineColor(kRed);
  ray_0->SetTitle(Title);
  com_0->SetMarkerStyle(kFullCircle);
  com_0->SetMarkerColor(kBlue);
  com_0->Draw("PE same");

  // TRatioPlot *rp = new TRatioPlot(ray_0, com_0);
  // rp->Draw();
  // gPad->Modified(); gPad->Update(); // make sure it’s really (re)drawn
  // TPad *p = rp->GetUpperPad();

  TLegend *leg;
  if(rise)
  {
    leg = gPad->BuildLegend(0.55,0.11,0.89,0.35);
  }
  else
  {
    leg = gPad->BuildLegend(0.55,0.65,0.89,0.89);
  }


  // TLegend *leg = gPad->BuildLegend(0.55,0.11,0.89,0.35);
  TList *prim = leg->GetListOfPrimitives();

  TIter next(prim);
  TObject *obj;
  TLegendEntry *le;
  int i = 0;

  while ((obj = next()))
  {
    le = (TLegendEntry*)obj;
    i++;
    if (i==1) le->SetLabel(Title_ray);
    if (i==2) le->SetLabel(Title_comp);
  }
  // p->Modified(); p->Update(); // make sure it’s really (re)drawn

  Double_t ks = ray_0->KolmogorovTest(com_0);
  // std::cout << ks << std::endl;
  return ks;


}

void usage()
{
  std::cout << "\t\t" << " [out_prefix] [hybrid_prefix] [max events] " << std::endl
            << "\t\t" << std::endl;
}
