// compile with
// g++ -o ../build/splitFluxFile splitFluxFile.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore

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

void usage();



int main(int argc, char** argv)
{

  if(argc < 3)
  {
    std::cout << argv[0];
    usage();
    return 1;
  }

  TString input_filename = argv[1];
  Int_t maxEvents = -1;
  std::string subfolder = "fluxFiles";

  if(argc > 2)
  {
    maxEvents = atoi(argv[2]);
  }
  if(argc > 3)
  {
    subfolder = argv[3];
  }


  struct stat sb;

  if (stat(subfolder.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode))
  {
    std::cout << "Subfolder exists" << std::endl;
  }
  else
  {
    // std::cout << "Creating directory " << subfolder.c_str() << std::endl;
    // const int dir_err = mkdir(subfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    // if (-1 == dir_err)
    // {
    std::cout << "ERROR! Directory " << subfolder.c_str() << " does not exists! Aborting..." << std::endl;
    exit(1);
    // }
  }


  TTree* source_tree;
  TFile* fInput;
  //input file
  Int_t    runNumber;
  Int_t    evtNumber;
  Int_t    evtIndex;
  Double_t prod_vertex_x;
  Double_t prod_vertex_y;
  Double_t prod_vertex_z;
  Double_t entry_x;
  Double_t entry_y;
  Double_t entry_z;
  Double_t px;
  Double_t py;
  Double_t pz;
  Int_t    G4index;
  Int_t    pdgID;
  Int_t    mother_index;
  Double_t charge;
  Double_t timing;
  Double_t eKinetic;
  Double_t eTot;
  Int_t    treeIndex;

  fInput = new TFile(input_filename.Data());
  source_tree = (TTree *) fInput->Get("tree");

  source_tree->SetBranchAddress("runNumber",&runNumber);
  source_tree->SetBranchAddress("evtNumber",&evtNumber);
  source_tree->SetBranchAddress("evtIndex",&evtIndex);
  source_tree->SetBranchAddress("prod_vertex_x",&prod_vertex_x);
  source_tree->SetBranchAddress("prod_vertex_y",&prod_vertex_y);
  source_tree->SetBranchAddress("prod_vertex_z",&prod_vertex_z);
  source_tree->SetBranchAddress("entry_x",&entry_x);
  source_tree->SetBranchAddress("entry_y",&entry_y);
  source_tree->SetBranchAddress("entry_z",&entry_z);
  source_tree->SetBranchAddress("px",&px);
  source_tree->SetBranchAddress("py",&py);
  source_tree->SetBranchAddress("pz",&pz);
  source_tree->SetBranchAddress("G4index",&G4index);
  source_tree->SetBranchAddress("pdgID",&pdgID);
  source_tree->SetBranchAddress("mother_index",&mother_index);
  source_tree->SetBranchAddress("charge",&charge);
  source_tree->SetBranchAddress("timing",&timing);
  source_tree->SetBranchAddress("eKinetic",&eKinetic);
  source_tree->SetBranchAddress("eTot",&eTot);
  source_tree->SetBranchAddress("treeIndex",&treeIndex);

  //output file
  TTree *out_tree = NULL;
  Int_t    out_runNumber;
  Int_t    out_evtNumber;
  Int_t    out_evtIndex;
  Double_t out_prod_vertex_x;
  Double_t out_prod_vertex_y;
  Double_t out_prod_vertex_z;
  Double_t out_entry_x;
  Double_t out_entry_y;
  Double_t out_entry_z;
  Double_t out_px;
  Double_t out_py;
  Double_t out_pz;
  Int_t    out_G4index;
  Int_t    out_pdgID;
  Int_t    out_mother_index;
  Double_t out_charge;
  Double_t out_timing;
  Double_t out_eKinetic;
  Double_t out_eTot;
  Int_t    out_treeIndex;

  source_tree->GetEntry(0);
  Int_t firstEvtID = std::max(evtIndex,0);

  Int_t maxEvtIndex = source_tree->GetMaximum("evtIndex");
  std::cout << "Input bunch crossings " << maxEvtIndex << std::endl;
  std::cout << "First event ID        " << firstEvtID << std::endl;




  TFile *fOut = NULL;

  long int nEvents = source_tree->GetEntries();
  Int_t bunchIndex = -1;
  std::stringstream sname;
  int fileCounter = 0;

  std::cout << "Creating files..." << std::endl;
  for(int i = 0 ; i < nEvents; i++)
  {
    source_tree->GetEntry(i);


    // std::cout << i << " " << evtIndex << " " <<  bunchIndex << std::endl;
    if(evtIndex < firstEvtID) continue;

    // check if index changes
    if(evtIndex != bunchIndex)
    {
      bunchIndex = evtIndex;
      if(evtIndex == firstEvtID)
      {
        // bunchIndex++;
        // std::cout << "new bunch " << bunchIndex << std::endl;
        sname << subfolder.c_str() << "/flux_" << fileCounter << ".root";
        fileCounter++;
        fOut = new TFile(sname.str().c_str(),"recreate");
        sname.str("");

        out_tree = new TTree("tree","tree");
        out_tree->Branch("runNumber"    ,&out_runNumber     ,"runNumber/I" );
        out_tree->Branch("evtNumber"    ,&out_evtNumber     ,"evtNumber/I" );
        out_tree->Branch("evtIndex"     ,&out_evtIndex      ,"evtIndex /I" );
        out_tree->Branch("prod_vertex_x",&out_prod_vertex_x ,"prod_vertex_x/D" );
        out_tree->Branch("prod_vertex_y",&out_prod_vertex_y ,"prod_vertex_y/D" );
        out_tree->Branch("prod_vertex_z",&out_prod_vertex_z ,"prod_vertex_z/D" );
        out_tree->Branch("entry_x"      ,&out_entry_x       ,"entry_x/D" );
        out_tree->Branch("entry_y"      ,&out_entry_y       ,"entry_y/D" );
        out_tree->Branch("entry_z"      ,&out_entry_z       ,"entry_z/D" );
        out_tree->Branch("px"           ,&out_px            ,"px/D" );
        out_tree->Branch("py"           ,&out_py            ,"py/D" );
        out_tree->Branch("pz"           ,&out_pz            ,"pz/D" );
        out_tree->Branch("G4index"      ,&out_G4index       ,"G4index/I" );
        out_tree->Branch("pdgID"        ,&out_pdgID         ,"pdgID/I" );
        out_tree->Branch("mother_index" ,&out_mother_index  ,"mother_index/I" );
        out_tree->Branch("charge"       ,&out_charge        ,"charge/D" );
        out_tree->Branch("timing"       ,&out_timing        ,"timing/D" );
        out_tree->Branch("eKinetic"     ,&out_eKinetic      ,"eKinetic/D" );
        out_tree->Branch("eTot"         ,&out_eTot          ,"eTot/D" );
        out_tree->Branch("treeIndex"    ,&out_treeIndex     ,"treeIndex/I" );
      }
      else
      {


        fOut->cd();
        out_tree->Write();
        fOut->Close();


        if (maxEvents != -1)
        {
          if(evtIndex > (maxEvents-1)) break;
        }


        // bunchIndex++;

        sname << subfolder.c_str() << "/flux_" << fileCounter << ".root";
        fileCounter++;
        fOut = new TFile(sname.str().c_str(),"recreate");
        sname.str("");
        // std::cout << "new bunch " << bunchIndex << std::endl;

        out_tree = new TTree("tree","tree");
        out_tree->Branch("runNumber"    ,&out_runNumber     ,"runNumber/I" );
        out_tree->Branch("evtNumber"    ,&out_evtNumber     ,"evtNumber/I" );
        out_tree->Branch("evtIndex"     ,&out_evtIndex      ,"evtIndex /I" );
        out_tree->Branch("prod_vertex_x",&out_prod_vertex_x ,"prod_vertex_x/D" );
        out_tree->Branch("prod_vertex_y",&out_prod_vertex_y ,"prod_vertex_y/D" );
        out_tree->Branch("prod_vertex_z",&out_prod_vertex_z ,"prod_vertex_z/D" );
        out_tree->Branch("entry_x"      ,&out_entry_x       ,"entry_x/D" );
        out_tree->Branch("entry_y"      ,&out_entry_y       ,"entry_y/D" );
        out_tree->Branch("entry_z"      ,&out_entry_z       ,"entry_z/D" );
        out_tree->Branch("px"           ,&out_px            ,"px/D" );
        out_tree->Branch("py"           ,&out_py            ,"py/D" );
        out_tree->Branch("pz"           ,&out_pz            ,"pz/D" );
        out_tree->Branch("G4index"      ,&out_G4index       ,"G4index/I" );
        out_tree->Branch("pdgID"        ,&out_pdgID         ,"pdgID/I" );
        out_tree->Branch("mother_index" ,&out_mother_index  ,"mother_index/I" );
        out_tree->Branch("charge"       ,&out_charge        ,"charge/D" );
        out_tree->Branch("timing"       ,&out_timing        ,"timing/D" );
        out_tree->Branch("eKinetic"     ,&out_eKinetic      ,"eKinetic/D" );
        out_tree->Branch("eTot"         ,&out_eTot          ,"eTot/D" );
        out_tree->Branch("treeIndex"    ,&out_treeIndex     ,"treeIndex/I" );
      }
    }


    out_runNumber     = runNumber;
    out_evtNumber     = evtNumber;
    out_evtIndex      = evtIndex;
    out_prod_vertex_x = prod_vertex_x;
    out_prod_vertex_y = prod_vertex_y;
    out_prod_vertex_z = prod_vertex_z;
    out_entry_x       = entry_x;
    out_entry_y       = entry_y;
    out_entry_z       = entry_z;
    out_px            = px;
    out_py            = py;
    out_pz            = pz;
    out_G4index       = G4index;
    out_pdgID         = pdgID;
    out_mother_index  = mother_index;
    out_charge        = charge;
    out_timing        = timing;
    out_eKinetic      = eKinetic;
    out_eTot          = eTot;
    out_treeIndex     = treeIndex;

    out_tree->Fill();

  }
  std::cout << "Created " << fileCounter << " files." << std::endl;



  return 0;
}


void usage()
{
  std::cout << "\t\t" << " <input_file_name> [max events] [subfolder]" << std::endl
            << "\t\t" << std::endl;
}
