#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TH1F.h"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RResultPtr.hxx"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <dirent.h>


// ### D E F I N E S ###
#define FILE_N    "pruned"       
#define TREE_N    "tree"          
#define MOD_TS_N  "mod0_t"
#define MOD_LO_N  "mod0_ph"
// #define LCh4_N    "ch4"           
// #define LCh13_N   "ch13"              
// #define TCh4_N    "t_pulse4"      
// #define TCh13_N   "t_pulse13"        
// // #define TCh4_N    "t_pulse_deposition4"      
// // #define TCh13_N   "t_pulse_deposition13"        
// // #define TCh4_N    "t_pulse_depoAndScint4"      
// // #define TCh13_N   "t_pulse_depoAndScint13"        
// // #define TCh4_N    "t_pulse_depoAndTrans4"      
// // #define TCh13_N   "t_pulse_depoAndTrans13"        
#define TL_N      "Total_Light"     
#define CFD_N     "CFD_Thresh"      
#define ENER_N    "folderEnergy"  
// #define CELL_1    4
// #define CELL_2    13


#define T_MIN     25
#define T_MAX     50
#define NBINS     75

#define NCELLSMAX 1000

struct dataCell_t
{
  // std::vector<float> energies;
  std::vector<float> LY;
  std::vector<float> err_LY;
  std::vector<float> sigma_t;
  std::vector<float> err_sigma_t;
};

// struct dataEnergy_t
// {
//   float energy;
//   float LYtot;
//   float sigma_LYtot;
//   float err_LYtot;
//   float err_sigma_LYtot;
  
//   TH1F h_LYtot;

//   std::vector<TH1F> h_LY;
//   std::vector<TH1F> h_t;

//   std::vector<float> LY;
//   std::vector<float> sigma_t;
//   std::vector<float> err_sigma_t;

// }




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




void analyseSignals (std::string path="-1", std::string outFileName="pulseAnalysed", std::string filesKeyWord=FILE_N, std::string folderKeyWord="GeV")
{
  time_t _start = clock();
  // ROOT::EnableImplicitMT();
  // --- Check input values
  if (path == "-1") { std::cout << "No path given in input. Usage:\n root \'analyseSignals.cxx(\"path\",  outFileName=\"pulseAnalysed\", filesKeyWord=" << filesKeyWord << ", folderKeyword=\"GeV\")\'\n"; return;}
  else std::cout << "Looking for files matching " << filesKeyWord << " into folders with " << folderKeyWord << " in " << path <<"\n";


  std::cout << "\n";
  std::cout << "----------------------------------\n";
  std::cout << "FILE_N   : " << filesKeyWord.c_str() << "\n";
  std::cout << "TREE_N   : " << TREE_N << "\n";
  std::cout << "----------------------------------\n";
  // std::cout << "CELL_1   : " << CELL_1 << "\n";
  // std::cout << "CELL_2   : " << CELL_2 << "\n";
  std::cout << "TL_N     : " << TL_N   << "\n";
  std::cout << "CFD_N    : " << CFD_N  << "\n";
  std::cout << "ENER_N   : " << ENER_N << "\n";
  std::cout << "----------------------------------\n";
  std::cout << std::endl;


  gStyle -> SetPalette(kRainBow);
  gStyle -> SetOptStat(1111);
  gStyle -> SetOptFit (1111);



  // --- Add '/' at the end of the input path if needed
  if (path.back() != '/') path += '/';



  // --- Find the folders matching the keyword
  std::vector <std::string> dataFolders, folders;
  std::vector <std::string> dataFiles, files;
  // std::vector <float> energies;



  read_directory( path, folders);
  for (auto && folder : folders) {
    if (folder.find(folderKeyWord) != string::npos) {
      dataFolders.emplace_back(path + folder + "/");                                                            
      read_directory((path + folder + "/"), files);
      for (auto && file : files) {
        if (file.find(filesKeyWord) != string::npos) {
          dataFiles.emplace_back(path + folder + "/" + file);
        }
      }
      files.clear();
    }    
  }

  // --- Feedback
  std::cout << "Found " << dataFolders.size() << " folders and " << dataFiles.size() << " files.\n";
  for (auto && elem : dataFiles) std::cout << elem << "\n";
  std::cout << std::endl;





  // --- Open the Output File
  TFile outFile ( (path + outFileName + ".root").c_str(), "RECREATE");



  // --- Prepare vectors of energies
  std::vector<float> energies;  
  std::vector<float> err_energies;  
  std::vector<float> v_enRes;
  std::vector<float> v_err_enRes;
  std::vector<dataCell_t> v_cells (NCELLSMAX);















  std::cout << "START LOOPING OVER THE FILES.\nLooking for branches " << MOD_LO_N << " and " << MOD_TS_N /*<< " and cells " << CELL_1 << " & " << CELL_2 */ << std::endl;
  for ( auto && file : dataFiles)   // LOOPING OVER THE ENERGIES!
  {
    std::cout << "\n\n... Opening " << file << "..." << std::endl;
    TFile inFile (file.c_str(), "READ");








    auto TT = (TTree*) inFile.Get(TREE_N);
    int light = 0;
    std::vector<float> *ts = 0;
    std::vector<int> *lo = 0;

    TT -> SetBranchAddress(TL_N, &light);
    TT -> SetBranchAddress(MOD_TS_N, &ts);
    TT -> SetBranchAddress(MOD_LO_N, &lo);
    
    
    







    
    auto energy = TT -> GetMaximum (ENER_N);
    energies.emplace_back(energy);
    err_energies.emplace_back(energy * 0.01);
    std::cout << "ENERGY MAX: " << energy << std::endl;

    auto max_LY = TT -> GetMaximum (TL_N);
    auto min_LY = TT -> GetMinimum (TL_N);
    auto h_LYtot = TH1F( ("Energy_H_"      + std::to_string(int(energy)) + "GeV").c_str(), ("Energy_H_"      + std::to_string(int(energy)) + "GeV").c_str(), NBINS,  min_LY, max_LY);


    int nCells = 0;
    std::vector<TH1F> h_LY;
    std::vector<TH1F> h_t;
    















    auto nEv = TT -> GetEntries();
    for ( int iEv = 0; iEv < nEv; ++iEv)
    {
      TT -> GetEvent(iEv);

      if (light <= 0) continue;



      // --- Initialise the histos
      if (nCells < ts->size()) {
        nCells = ts->size();

        std::cout << "- Initialising the vectors of TH1F\n";
        std::cout << "nCells: " << nCells << std::endl;

        for (size_t iCell = 0; iCell < ts->size(); ++iCell) {
          h_LY.emplace_back( TH1F ( ("H_LY_Cell"+std::to_string(iCell) + "_"      + std::to_string(int(energy)) + "GeV").c_str(), ("H_LY_Cell"+std::to_string(iCell) + "_"      + std::to_string(int(energy)) + "GeV").c_str(), int( NBINS*(max_LY)/(max_LY - min_LY) ),  0, max_LY) );
          h_t .emplace_back( TH1F ( ("H_TS_Cell"+std::to_string(iCell) + "_"      + std::to_string(int(energy)) + "GeV").c_str(), ("H_TS_Cell"+std::to_string(iCell) + "_"      + std::to_string(int(energy)) + "GeV").c_str(), int( NBINS* (T_MAX - T_MIN) ),  T_MIN, T_MAX) );
        }

        std::cout << "Histo LY: " << h_LY.size() << "\tHisto T: " << h_t.size() << std::endl;
      }





      // --- Update the Histos
      h_LYtot.Fill(light);
      for (size_t iCell = 0; iCell < ts->size(); ++iCell) {
        if (ts->at(iCell) > 0) h_t .at(iCell).Fill(ts->at(iCell));
        if (lo->at(iCell) > 0) h_LY.at(iCell).Fill(lo->at(iCell));
      }

    } // <- Event loop








    if (v_cells.size() != h_t.size()) {
      std::cout << std::endl << std::endl;
      std::cout << "######################################################################################################" << std::endl;
      std::cout << "### Reshaping the cells vectors. If this message appears more than once, it's a (serious) problem. ###" << std::endl;
      std::cout << "######################################################################################################" << std::endl;
      std::cout << std::endl << std::endl;
      v_cells.resize(h_t.size());
    }  








    // --- Get the values
    //    at this stage I have histos of LYtot and the vectors of histos, one for each cell.

    // TF1 fgaus ("fgaus", "gaus");

    // - LYtot
    // h_LYtot.Fit("fgaus");
    v_enRes.push_back       ( h_LYtot.GetStdDev() / h_LYtot.GetMean() );
    v_err_enRes.push_back ( h_LYtot.GetStdDev() / h_LYtot.GetMean() * sqrt( pow(h_LYtot.GetStdDevError() / h_LYtot.GetStdDev(),2) + pow(h_LYtot.GetMeanError() / h_LYtot.GetMean(),2)));

    outFile.cd();
    h_LYtot.Write();



    for (int iCell = 0; iCell < nCells; ++iCell) {
      // - LY
      v_cells.at(iCell).LY.push_back (h_LY.at(iCell).GetMean());
      v_cells.at(iCell).err_LY.push_back (h_LY.at(iCell).GetMeanError());
    
      // - TRes
      v_cells.at(iCell).sigma_t.push_back (h_t.at(iCell).GetStdDev());
      v_cells.at(iCell).err_sigma_t.push_back (h_t.at(iCell).GetStdDevError());

      // --- Write to disk
      std::string dirName = "Cell_" + std::to_string(iCell);
      outFile.cd();
      if (! outFile.GetDirectory(dirName.c_str())) outFile.mkdir(dirName.c_str());
      outFile.cd(dirName.c_str());

      h_LY.at(iCell).Write();
      h_t .at(iCell).Write();
    }


    std::cout << std::endl;
  } //  <- END Loop over the files




  outFile.cd();

  // --- GRAPH: EnRes vs En 
  TCanvas c_enRes ("c_enRes", "c_enRes", 1000, 600);
  c_enRes.SetGrid();
  TGraphErrors g_energy (energies.size(), energies.data(), v_enRes.data(), err_energies.data(), v_err_enRes.data());
  g_energy.SetMarkerStyle(20);
  g_energy.SetMarkerSize(1.3);
  g_energy.SetMarkerColor(kRed+1);
  g_energy.GetXaxis()->SetTitle("Energy [GeV]");
  g_energy.GetYaxis()->SetTitle("#sigma_{E}/E");
  g_energy.Draw("AP");
  TF1 fitEnRes ("fitEnRes", "sqrt( pow([0],2)/x + pow([1],2) )");
  fitEnRes.SetParameter(0, 0.1);
  fitEnRes.SetParameter(1, 0.01);
  g_energy.Fit("fitEnRes");
  g_energy.Write("EnergyRes_vs_Energy");
  c_enRes.Write();





  // --- GRAPHs: CELLS
  for (size_t iCell = 0; iCell < v_cells.size(); ++iCell) {
    std::string dirName = "Cell_" + std::to_string(iCell);
    outFile.cd(dirName.c_str());


    // - Time   
    TCanvas c_TRes (("c_TRes_"+dirName).c_str(), ("c_TRes_"+dirName).c_str(), 1000, 600);
    c_TRes.SetGrid();
    TGraphErrors g_time (energies.size(), energies.data(), v_cells.at(iCell).sigma_t.data(), err_energies.data(), v_cells.at(iCell).err_sigma_t.data());
    g_time.SetMarkerStyle(20);
    g_time.SetMarkerSize(1.3);
    g_time.SetMarkerColor(kBlue+1);
    g_time.GetXaxis()->SetTitle("Energy [GeV]");
    g_time.GetYaxis()->SetTitle("#sigma_{t} [ns]");
    g_time.Draw("AP");
    g_time.Write(("TRes_"+dirName).c_str());
    c_TRes.Write();

    // - LY
    TCanvas c_LY (("c_LY_"+dirName).c_str(), ("c_LY_"+dirName).c_str(), 1000, 600);
    c_LY.SetGrid();
    TGraphErrors g_LY (energies.size(), energies.data(), v_cells.at(iCell).LY.data(), err_energies.data(), v_cells.at(iCell).err_LY.data());
    g_LY.SetMarkerStyle(20);
    g_LY.SetMarkerSize(1.3);
    g_LY.GetXaxis()->SetTitle("Energy [GeV]");
    g_LY.GetYaxis()->SetTitle("Photons Detected");
    g_LY.Draw("AP");
    g_LY.Write(("LY_"+dirName).c_str());
    c_LY.Write();

  }








  outFile.Close();
  std::cout << "\n";
  std::cout << "Saved to: " << (path + outFileName + ".root").c_str() << "\n";
  std::cout << "### Duration: " << (std::clock() - _start) / (double) CLOCKS_PER_SEC << "s\n";
  std::cout << "### Standard end of the program. ###" << std::endl;
  gApplication->Terminate();
}
