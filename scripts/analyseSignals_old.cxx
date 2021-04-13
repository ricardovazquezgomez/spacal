#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
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
#define LCh4_N    "ch4"           
#define LCh13_N   "ch13"              
#define TCh4_N    "t_pulse4"      
#define TCh13_N   "t_pulse13"        
// #define TCh4_N    "t_pulse_deposition4"      
// #define TCh13_N   "t_pulse_deposition13"        
// #define TCh4_N    "t_pulse_depoAndScint4"      
// #define TCh13_N   "t_pulse_depoAndScint13"        
// #define TCh4_N    "t_pulse_depoAndTrans4"      
// #define TCh13_N   "t_pulse_depoAndTrans13"        
#define TL_N      "Total_Light"     
#define CFD_N     "CFD_Thresh"      
#define ENER_N    "folderEnergy"  


#define NBINS     75

// struct dataEnergy_t
// {
//   float zero = -1;

//   // - Energy
//   float trueEn      = -1;        
//   float measEn      = -1;            
//   float sigmaEn     = -1;       
//   float enRes       = -1;         
//   float sigmaEnRes  = -1;         

//   // - Time
//   float tRes_Ch4              = -1;      
//   float tRes_Ch13             = -1;     
//   float sigmatRes_Ch4         = -1;
//   float sigmatRes_Ch13        = -1;
//   float sigma_sigmatRes_Ch4   = -1;
//   float sigma_sigmatRes_Ch13  = -1;

// };

struct dataThreshold_t
{
  float CFD_Threshold;
  std::vector<double> Energy_Col;

  std::vector <ROOT::RDF::RResultPtr<double>> v_lTotMax;    
  std::vector <ROOT::RDF::RResultPtr<double>> v_lTotMin;    
  std::vector <ROOT::RDF::RResultPtr<double>> v_lTot_Mean; 
  std::vector <ROOT::RDF::RResultPtr<double>> v_lTot_Sigma; 


  std::vector <ROOT::RDF::RResultPtr<double>> v_lMax_Ch4;  
  std::vector <ROOT::RDF::RResultPtr<double>> v_lMax_Ch13; 
  std::vector <ROOT::RDF::RResultPtr<double>> v_lMin_Ch4;  
  std::vector <ROOT::RDF::RResultPtr<double>> v_lMin_Ch13; 
  std::vector <ROOT::RDF::RResultPtr<double>> v_lCh4_Mean;  
  std::vector <ROOT::RDF::RResultPtr<double>> v_lCh13_Mean; 
  std::vector <ROOT::RDF::RResultPtr<double>> v_lCh4_Sigma;  
  std::vector <ROOT::RDF::RResultPtr<double>> v_lCh13_Sigma; 


  std::vector <ROOT::RDF::RResultPtr<double>> v_tMax_Ch4;   
  std::vector <ROOT::RDF::RResultPtr<double>> v_tMax_Ch13;  
  std::vector <ROOT::RDF::RResultPtr<double>> v_tMin_Ch4;   
  std::vector <ROOT::RDF::RResultPtr<double>> v_tMin_Ch13;  
  std::vector <ROOT::RDF::RResultPtr<double>> v_tCh4_Mean;  
  std::vector <ROOT::RDF::RResultPtr<double>> v_tCh13_Mean; 
  std::vector <ROOT::RDF::RResultPtr<double>> v_tCh4_Sigma; 
  std::vector <ROOT::RDF::RResultPtr<double>> v_tCh13_Sigma;


  std::vector <ROOT::RDF::RResultPtr<TH1D>> h_lTot ;
  std::vector <ROOT::RDF::RResultPtr<TH1D>> h_lCh4 ;
  std::vector <ROOT::RDF::RResultPtr<TH1D>> h_lCh13;
  std::vector <ROOT::RDF::RResultPtr<TH1D>> h_tCh4 ;
  std::vector <ROOT::RDF::RResultPtr<TH1D>> h_tCh13;


};





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
  std::cout << "FILE_N   : " << FILE_N << "\n";
  std::cout << "TREE_N   : " << TREE_N << "\n";
  std::cout << "----------------------------------\n";
  std::cout << "LCh4_N   : " << LCh4_N << "\n"; 
  std::cout << "LCh13_N  : " << LCh13_N<< "\n";
  std::cout << "TCh4_N   : " << TCh4_N << "\n";
  std::cout << "TCh13_N  : " << TCh13_N<< "\n";
  std::cout << "TL_N     : " << TL_N   << "\n";
  std::cout << "CFD_N    : " << CFD_N  << "\n";
  std::cout << "ENER_N   : " << ENER_N << "\n";
  std::cout << "----------------------------------\n";
  std::cout << "\n";


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



  // // --- Prepare the map for the groupby
  // //    This map will contain all the information at a given threshold.
  // //    The key will be the threshold itself, the struct will contain a series of vectors.
  // std::map<float, dataThreshold_t> map;





  // --- Create the DataFrame and filter out events <= 0
  ROOT::RDataFrame dPreFilt ( TREE_N, dataFiles);
  auto dFilt  = dPreFilt.Filter ( (std::string(TL_N)   + " > 0").c_str() )
                        .Filter ( (std::string(TCh4_N) + " > 0").c_str() )
                        .Filter ( (std::string(TCh13_N)+ " > 0").c_str() );
  // auto dEn  = dPreFilt.Filter ( (std::string(TL_N)   + " > 0").c_str() ); 
  // auto dT   = dPreFilt.Filter ( (std::string(TL_N)   + " > 0").c_str() ) 
  //                     .Filter ( (std::string(TCh4_N) + " > 0").c_str() )
  //                     .Filter ( (std::string(TCh13_N)+ " > 0").c_str() );



  // --- Get all the different CFD threshold values and energies
  //      Take the column and remove the non unique elements.
  //      Then we will perform a for loop on all the different CFD Thresholds.
  auto CFD_Col    = *dPreFilt.Take<float>(CFD_N);
  auto Energy_Col = *dPreFilt.Take<float>(ENER_N);

  std::sort(CFD_Col.begin(), CFD_Col.end());
  CFD_Col.erase( std::unique( CFD_Col.begin(), CFD_Col.end()), CFD_Col.end() );

  std::cout << "Thresholds found: ";
  for (auto & threshold : CFD_Col) std::cout << threshold << " ";
  std::cout << std::endl;


  // -- Now all the energies
  std::sort(Energy_Col.begin(), Energy_Col.end());
  Energy_Col.erase( std::unique( Energy_Col.begin(), Energy_Col.end()), Energy_Col.end() );

  std::cout << "Energies found: ";
  for (auto & energy : Energy_Col) std::cout << energy << " ";
  std::cout << std::endl;


  // --- Build the struct
  std::vector<dataThreshold_t> v_data;

  // --- Prepare the lambdas
  std::vector<std::function<bool(float)>> Energy_Lambdas;
  for ( auto & Energy : Energy_Col) Energy_Lambdas.emplace_back ( [&Energy] (const float & energy) ->bool{ /*std::cout << "En: " << Energy << std::endl; */return Energy == energy; });

  std::vector<std::function<bool(float)>> CFD_Lambdas;
  for ( auto & CFD : CFD_Col) CFD_Lambdas.emplace_back ( [&CFD] (const float & thresh) ->bool{ /*std::cout << "CFD: " <<  CFD << std::endl; */return CFD == thresh; });

  // double CFD_thresh, trueEn;
  // auto FilterThreshold  = [&CFD_thresh] (const float & thresh) ->bool{ return CFD_thresh == thresh;};
  // auto FilterEnergy     = [&trueEn]     (const float & energy) ->bool{ return trueEn     == energy;};







  // --- FIRST LOOP ---
  std::cout << "LOOP 1.\nCFD\tEnergies"  << "\n";
  // --> THRESHOLD -->
  for (size_t iCFD = 0; iCFD < CFD_Col.size(); ++iCFD) {

    std::cout << CFD_Col.at(iCFD) << "\t";
    dataThreshold_t data;

    data.CFD_Threshold = CFD_Col.at(iCFD);
    for (auto & energy : Energy_Col) data.Energy_Col.emplace_back(energy);

    // --> ENERGY -->
    for (size_t iEn = 0; iEn < Energy_Col.size(); ++ iEn) {
      std::cout << Energy_Col.at(iEn) << " ";      
      // --- Get statistical quantities
      data.v_lTotMax    .push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Max(TL_N));
      data.v_lTotMin    .push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Min(TL_N));
      data.v_lTot_Mean  .push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Mean(TL_N));
      data.v_lTot_Sigma .push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).StdDev(TL_N));


      data.v_lMax_Ch4   .push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Max(LCh4_N));
      data.v_lMax_Ch13  .push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Max(LCh13_N));
      data.v_lMin_Ch4   .push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Min(LCh4_N));
      data.v_lMin_Ch13  .push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Min(LCh13_N));
      data.v_lCh4_Mean  .push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Mean(LCh4_N));
      data.v_lCh13_Mean .push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Mean(LCh13_N));
      data.v_lCh4_Sigma .push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).StdDev(LCh4_N));
      data.v_lCh13_Sigma.push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).StdDev(LCh13_N));

      data.v_tMax_Ch4   .push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Max(TCh4_N));
      data.v_tMax_Ch13  .push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Max(TCh13_N));
      data.v_tMin_Ch4   .push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Min(TCh4_N));
      data.v_tMin_Ch13  .push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Min(TCh13_N));
      data.v_tCh4_Mean  .push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Mean(TCh4_N));
      data.v_tCh13_Mean .push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Mean(TCh13_N));
      data.v_tCh4_Sigma .push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).StdDev(TCh4_N));
      data.v_tCh13_Sigma.push_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).StdDev(TCh13_N));



      // RDFs.emplace_back(dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Mean(TL_N));
    }
    std::cout << "\n";
    // <-- ENERGY <--

    v_data.emplace_back(data);
    // std::cout << "Values: ";
    // for ( auto & RDF : RDFs) std::cout << *RDF << " ";
  }
  // <-- THRESHOLD <--

  std::cout << "\nPRINT.\nCFD\tEnergies\n";
  for ( auto & data: v_data) {
    std::cout << data.CFD_Threshold << " ";
    for (auto & en: data.v_lTot_Mean) std::cout << *en << " ";
    std::cout << "\n";
  }











  // --- LOOP ANEW TO MAKE THE HISTOS ---
  // --- SECOND LOOP ---
  std::cout << "\n.LOOP 2.\nCFD\tEnergies"  << "\n";
  // --> THRESHOLD -->
  for (size_t iCFD = 0; iCFD < v_data.size(); ++iCFD) {

    float Thresh = v_data.at(iCFD).CFD_Threshold;
    float Energy = -1;

    std::string CFD_suffix = "CFD_" + std::to_string(int(Thresh*100)) + "_perc";
    std::cout << Thresh << "\t";

    // --> ENERGY -->
    for (size_t iEn = 0; iEn < v_data.at(iCFD).Energy_Col.size(); ++ iEn) {
      Energy = v_data.at(iCFD).Energy_Col.at(iEn);
      std::cout << Energy << " ";      

      // --- Get statistical quantities
      v_data.at(iCFD).h_lTot .push_back( dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Histo1D({("Energy_H_"      + std::to_string(int(Energy)) + "GeV_" + CFD_suffix).c_str(), ("Energy_H_"      + std::to_string(int(Energy)) + "GeV_" + CFD_suffix).c_str(), NBINS, /*(*v_data[iCFD].v_lTotMin[iEn])   , (*v_data[iCFD].v_lTotMax[iEn])  */ (*v_data[iCFD].v_lTot_Mean[iEn])   - 5*(*v_data[iCFD].v_lTot_Sigma[iEn])  , (*v_data[iCFD].v_lTot_Mean[iEn])   + 5*(*v_data[iCFD].v_lTot_Sigma[iEn])    }, TL_N )     );
      v_data.at(iCFD).h_lCh4 .push_back( dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Histo1D({("Photons_H_Ch4_" + std::to_string(int(Energy)) + "GeV_" + CFD_suffix).c_str(), ("Photons_H_Ch4_" + std::to_string(int(Energy)) + "GeV_" + CFD_suffix).c_str(), NBINS, (*v_data[iCFD].v_lMin_Ch4[iEn])  , (*v_data[iCFD].v_lMax_Ch4[iEn])  /*(*v_data[iCFD].v_lCh4_Mean[iEn])  - 5*(*v_data[iCFD].v_lCh4_Sigma[iEn]) , (*v_data[iCFD].v_lCh4_Mean[iEn])  + 5*(*v_data[iCFD].v_lCh4_Sigma[iEn]) */ }, LCh4_N )   );
      v_data.at(iCFD).h_lCh13.push_back( dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Histo1D({("Photons_H_Ch13_"+ std::to_string(int(Energy)) + "GeV_" + CFD_suffix).c_str(), ("Photons_H_Ch13_"+ std::to_string(int(Energy)) + "GeV_" + CFD_suffix).c_str(), NBINS, (*v_data[iCFD].v_lMin_Ch13[iEn]) , (*v_data[iCFD].v_lMax_Ch13[iEn]) /*(*v_data[iCFD].v_lCh13_Mean[iEn]) - 5*(*v_data[iCFD].v_lCh13_Sigma[iEn]), (*v_data[iCFD].v_lCh13_Mean[iEn]) + 5*(*v_data[iCFD].v_lCh13_Sigma[iEn])*/ }, LCh13_N )  );
      v_data.at(iCFD).h_tCh4 .push_back( dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Histo1D({("Time_Res_Ch4_"  + std::to_string(int(Energy)) + "GeV_" + CFD_suffix).c_str(), ("Time_Res_Ch4_"  + std::to_string(int(Energy)) + "GeV_" + CFD_suffix).c_str(), NBINS, /*(*v_data[iCFD].v_tMin_Ch4[iEn])  , (*v_data[iCFD].v_tMax_Ch4[iEn]) */ (*v_data[iCFD].v_tCh4_Mean[iEn])  - 5*(*v_data[iCFD].v_tCh4_Sigma[iEn]) , (*v_data[iCFD].v_tCh4_Mean[iEn])  + 5*(*v_data[iCFD].v_tCh4_Sigma[iEn])   }, TCh4_N )   );
      v_data.at(iCFD).h_tCh13.push_back( dFilt.Filter(CFD_Lambdas.at(iCFD), {CFD_N}).Filter(Energy_Lambdas.at(iEn), {ENER_N}).Histo1D({("Time_Res_Ch13_" + std::to_string(int(Energy)) + "GeV_" + CFD_suffix).c_str(), ("Time_Res_Ch13_" + std::to_string(int(Energy)) + "GeV_" + CFD_suffix).c_str(), NBINS, /*(*v_data[iCFD].v_tMin_Ch13[iEn]) , (*v_data[iCFD].v_tMax_Ch13[iEn])*/ (*v_data[iCFD].v_tCh13_Mean[iEn]) - 5*(*v_data[iCFD].v_tCh13_Sigma[iEn]), (*v_data[iCFD].v_tCh13_Mean[iEn]) + 5*(*v_data[iCFD].v_tCh13_Sigma[iEn])  }, TCh13_N )  );
    }
    std::cout << "\n";
    // <-- ENERGY <--
  }
  // <-- THRESHOLD <--





  // --- THIRD LOOP ---
  // --- LOOP ANEW TO PLOT AND SAVE DATA ---
  std::cout << "\n.LOOP 3.\nCFD\tEnergies"  << "\n";
  TMultiGraph mg_tResCh4, mg_tResCh13;
  int iColor = 0;  
  // --> THRESHOLD -->
  for (auto & data : v_data)
  {

  // for (size_t iCFD = 0; iCFD < v_data.size(); ++iCFD) {
    // float Thresh = v_data.at(iCFD).CFD_Threshold;
    float Thresh = data.CFD_Threshold;
    std::string CFD_suffix = "CFD_" + std::to_string(int(Thresh*100)) + "_perc";
    std::cout << Thresh << "\t";

    gStyle -> SetPalette(kRainBow);
    gStyle -> SetOptStat(1111);
    gStyle -> SetOptFit (1111);
    // --- For every threshold prepare vectors for the plots ---
    std::vector<double> energy_resolution, sigma_energy_resolution, time_resolution_Ch4, sigma_time_resolution_Ch4, time_resolution_Ch13, sigma_time_resolution_Ch13;

    // --- MAKE DIRECTORIES IN THE OUTFILE
    std::string dirName = "Threshold_" + CFD_suffix + "/";
    outFile.cd();
    outFile.mkdir(dirName.c_str());
    outFile.cd(dirName.c_str());

    // --> ENERGY -->
    size_t energies = data.Energy_Col.size();
    for (size_t iEn = 0; iEn < energies; ++iEn) {
      std::cout << data.Energy_Col.at(iEn) << " ";      
      gStyle -> SetPalette(kRainBow);


      // --- GAUSSIAN FITS TO TOTAL ENERGY AND TIME RESOLUTIONS ---
      TF1 fgaus ("fgaus", "gaus");

      // - Energy resolution
      // // - FIT 
      // data.h_lTot .at(iEn) -> Fit("fgaus", "Q");
      // energy_resolution.emplace_back      ( fgaus.GetParameter(2)/fgaus.GetParameter(1) );
      // sigma_energy_resolution.emplace_back( fgaus.GetParameter(2)/fgaus.GetParameter(1) * sqrt( pow(fgaus.GetParError(2)/fgaus.GetParameter(2),2) + pow(fgaus.GetParError(1)/fgaus.GetParameter(1),2)));
      // - Sigma
      energy_resolution.emplace_back      ( data.h_lTot.at(iEn)->GetStdDev() / data.h_lTot.at(iEn)->GetMean() );
      sigma_energy_resolution.emplace_back( data.h_lTot.at(iEn)->GetStdDev() / data.h_lTot.at(iEn)->GetMean() * sqrt( pow(data.h_lTot.at(iEn)->GetStdDevError()/data.h_lTot.at(iEn)->GetStdDev(),2) + pow(data.h_lTot.at(iEn)->GetMeanError()/data.h_lTot.at(iEn)->GetMean(),2)));

      // - Time resolution ch4
      data.h_tCh4 .at(iEn) -> Fit("fgaus", "Q");
      time_resolution_Ch4.emplace_back(fgaus.GetParameter(2));
      sigma_time_resolution_Ch4.emplace_back(fgaus.GetParError(2));
    

      // - Time resolution ch13
      data.h_tCh13.at(iEn) -> Fit("fgaus", "Q");
      time_resolution_Ch13.emplace_back(fgaus.GetParameter(2));
      sigma_time_resolution_Ch13.emplace_back(fgaus.GetParError(2));


      // ---> Save Stuff
      data.h_lTot .at(iEn)->Write();
      data.h_lCh4 .at(iEn)->Write();
      data.h_lCh13.at(iEn)->Write();
      data.h_tCh4 .at(iEn)->Write();
      data.h_tCh13.at(iEn)->Write();
    }
    std::cout << "\n";
    // <-- ENERGY <--



    // --- PLOT VS ENERGY
    std::vector<double> zeros (energies, 0.);

    // --- Energy resolution vs. Energy
    TCanvas c_EnRes ("c_EnRes", "c_EnRes", 800, 600);
    TGraphErrors g_EnRes (energies, &data.Energy_Col[0], &energy_resolution[0], &zeros[0], &sigma_energy_resolution[0]);
    g_EnRes.SetMarkerStyle(20);
    g_EnRes.SetMarkerColor(kRed+1);
    g_EnRes.SetMarkerSize(1.3);
    g_EnRes.GetXaxis() -> SetTitle("Energy [GeV]");
    g_EnRes.GetYaxis() -> SetTitle("#sigma_{E} / E");

    std::cout << "ENERGY RESOLUTION VS ENERGY FIT:" << std::endl;
    TF1 fitEnRes ("fitEnRes", "sqrt( ([0]/sqrt(x))**2 + ([1])**2 )");
    fitEnRes.SetParameter(0, 0.1);
    fitEnRes.SetParameter(1, 0.01);
    fitEnRes.SetParName(0, "Sampling");
    fitEnRes.SetParName(1, "Constant");
    g_EnRes.Fit("fitEnRes");
    g_EnRes.Draw("AP");
    g_EnRes.Write(("Energy Resolution - " + CFD_suffix).c_str());
    c_EnRes.SetGrid();
    c_EnRes.Write();




    // --- Time resolution vs. Energy
    TCanvas c_tResCh4 ("c_tResCh4", "c_tResCh4", 800, 600);
    TGraphErrors *g_tResCh4 = new TGraphErrors (energies, &data.Energy_Col[0], &time_resolution_Ch4[0], &zeros[0], &sigma_time_resolution_Ch4[0]);
    g_tResCh4 -> SetMarkerStyle(20);
    g_tResCh4 -> SetMarkerColor(kBlue+iColor);
    g_tResCh4 -> SetMarkerSize(1.3);
    g_tResCh4 -> GetXaxis() -> SetTitle("Energy [GeV]");
    g_tResCh4 -> GetYaxis() -> SetTitle("#sigma_{T} [ns]");
    g_tResCh4 -> SetTitle(CFD_suffix.c_str());

    g_tResCh4 -> Draw("AP");
    g_tResCh4 -> Write(("Time Resolution Ch4 - " + CFD_suffix).c_str());
    c_tResCh4.SetGrid();
    c_tResCh4.Write();

    mg_tResCh4.Add(g_tResCh4);

    TCanvas c_tResCh13 ("c_tResCh13", "c_tResCh13", 800, 600);
    TGraphErrors *g_tResCh13 = new TGraphErrors (energies, &data.Energy_Col[0], &time_resolution_Ch13[0], &zeros[0], &sigma_time_resolution_Ch13[0]);
    g_tResCh13 -> SetMarkerStyle(20);
    g_tResCh13 -> SetMarkerColor(kBlue+iColor);
    g_tResCh13 -> SetMarkerSize(1.3);
    g_tResCh13 -> GetXaxis() -> SetTitle("Energy [GeV]");
    g_tResCh13 -> GetYaxis() -> SetTitle("#sigma_{T} [ns]");
    g_tResCh13 -> SetTitle(CFD_suffix.c_str());

    g_tResCh13 -> Draw("AP");
    g_tResCh13 -> Write(("Time Resolution Ch13 - " + CFD_suffix).c_str());
    c_tResCh13.SetGrid();
    c_tResCh13.Write();

    mg_tResCh13.Add(g_tResCh13);

    iColor += 3;
  }
  // <-- THRESHOLD <--


  outFile.cd();
  TCanvas c_tResCh4  ("c_tResCh4", "c_tResCh4", 800, 600);
  mg_tResCh4.GetXaxis() -> SetTitle("Energy [GeV]");
  mg_tResCh4.GetYaxis() -> SetTitle("#sigma_{T} [ns]");
  mg_tResCh4.Draw("AP");
  c_tResCh4.BuildLegend();
  mg_tResCh4.Write ("Time Resolution Ch4 for different thresholds");
  c_tResCh4.SetGrid();
  c_tResCh4.Write();

  TCanvas c_tResCh13 ("c_tResCh13", "c_tResCh13", 800, 600);
  mg_tResCh13.GetXaxis() -> SetTitle("Energy [GeV]");
  mg_tResCh13.GetYaxis() -> SetTitle("#sigma_{T} [ns]");
  mg_tResCh13.Draw("AP");
  c_tResCh13.BuildLegend();
  mg_tResCh13.Write("Time Resolution Ch13 for different thresholds");
  c_tResCh13.SetGrid();
  c_tResCh13.Write();







  // --- NOW DRAW STUFF


  // for (auto & threshold : CFD_Col) {
  //   CFD_thresh = threshold;
    
  //   for (auto & energy : Energy_Col) {
  //     data.v_lTot_Mean  .push_back(dEn.Mean(TL_N));
  //     data.v_lTot_Sigma .push_back(dEn.StdDev(TL_N));

  //   }

  // }


// #####################################################################
// #####################################################################
// #####################################################################
// #####################################################################
// #####################################################################



//     // --- Book all the operations ---
//     // --> Loop over the thresholds -->
//     std::cout << "Book the operations.\nThreshold\tEnergy\n";
//     for ( auto && threshold : CFD_Col)
//     {
//       dataThreshold_t data;
//       data.CFD_Threshold = threshold;
//       for (auto && energy : Energy_Col) data.Energy_Col.emplace_back(energy);
//       // data.Energy_Col = Energy_Col;


//       CFD_thresh = threshold;
//       std::cout << CFD_thresh << "\t";
//       auto dThresh = dFilt.Filter(FilterThreshold, {CFD_N});


//       // --> Loop over the energies -->
//       for (auto & energy : data.Energy_Col)
//       {
//         trueEn  = energy;
//         std::cout << trueEn << " ";
//         auto dEn = dThresh.Filter(FilterEnergy, {ENER_N});

//         // --- Get statistical quantities
//         data.v_lTotMax    .push_back(dEn.Max(TL_N));
//         data.v_lTotMin    .push_back(dEn.Min(TL_N));
//         data.v_lTot_Mean  .push_back(dEn.Mean(TL_N));
//         data.v_lTot_Sigma .push_back(dEn.StdDev(TL_N));


//         data.v_lMax_Ch4   .push_back(dEn.Max(LCh4_N));
//         data.v_lMax_Ch13  .push_back(dEn.Max(LCh13_N));
//         data.v_lMin_Ch4   .push_back(dEn.Min(LCh4_N));
//         data.v_lMin_Ch13  .push_back(dEn.Min(LCh13_N));


//         data.v_tMax_Ch4   .push_back(dEn.Max(TCh4_N));
//         data.v_tMax_Ch13  .push_back(dEn.Max(TCh13_N));
//         data.v_tMin_Ch4   .push_back(dEn.Min(TCh4_N));
//         data.v_tMin_Ch13  .push_back(dEn.Min(TCh13_N));
//         data.v_tCh4_Mean  .push_back(dEn.Mean(TCh4_N));
//         data.v_tCh13_Mean .push_back(dEn.Mean(TCh13_N));
//         data.v_tCh4_Sigma .push_back(dEn.StdDev(TCh4_N));
//         data.v_tCh13_Sigma.push_back(dEn.StdDev(TCh13_N));

//         // --- Trigger the loop. Not ideal, but I have some problems with RDF
//         //      so I must use this trick.
//         *dEn.Max(TL_N);

//       }
//       // <-- Loop over the energies <--

//       v_data.push_back(data);
//       std::cout << std::endl;

//     }
//     // <-- Loop over the thresholds <--
//     std::cout << "\nDone.\n\n";









// for (auto && elem : v_data[0].v_lTot_Mean) std::cout << "Energy: " << *elem << std::endl;









//     // --- Once more to book the histograms, and perform the first loop.
//     // --> Loop over the thresholds -->
//     std::cout << "Book the Histos. Bin number: " << NBINS /*<< "for energy histos, " << 2*NBINS << " for time histos."*/ << "\nThreshold\tEnergy\n";
//     for ( size_t iThresh = 0; iThresh < CFD_Col.size(); ++iThresh)
//     {
//       CFD_thresh = v_data[iThresh].CFD_Threshold;
//       std::cout << CFD_thresh << "\t";
//       auto dThresh = dFilt.Filter(FilterThreshold, {CFD_N});

//       std::string CFD_suffix = "CFD_" + std::to_string(int(std::round(CFD_thresh*100))) + "perc";


//       // --> Loop over the energies -->
//       for ( size_t iEn = 0; iEn < v_data[iThresh].Energy_Col.size(); ++iEn)
//       {
//         trueEn  = v_data[iThresh].Energy_Col[iEn];
//         std::cout << trueEn << " ";
//         auto dEn = dThresh.Filter(FilterEnergy, {ENER_N});

//         // --- Book the histos
//         v_data[iThresh].h_lTot  .push_back(dEn.Histo1D({("Energy_H_"      + std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str(), ("Energy_H_"      + std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str(), NBINS, (*v_data[iThresh].v_lTotMin[iEn])   * 0.9, (*v_data[iThresh].v_lTotMax[iEn])  * 1.1 }, TL_N )     );
//         v_data[iThresh].h_lCh4  .push_back(dEn.Histo1D({("Photons_H_Ch4_" + std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str(), ("Photons_H_Ch4_" + std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str(), NBINS, (*v_data[iThresh].v_lMin_Ch4[iEn])  * 0.9, (*v_data[iThresh].v_lMax_Ch4[iEn]) * 1.1 }, LCh4_N )   );
//         v_data[iThresh].h_lCh13 .push_back(dEn.Histo1D({("Photons_H_Ch13_"+ std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str(), ("Photons_H_Ch13_"+ std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str(), NBINS, (*v_data[iThresh].v_lMin_Ch13[iEn]) * 0.9, (*v_data[iThresh].v_lMax_Ch13[iEn])* 1.1 }, LCh13_N )  );
//         v_data[iThresh].h_tCh4  .push_back(dEn.Histo1D({("Time_Res_Ch4_"  + std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str(), ("Time_Res_Ch4_"  + std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str(), NBINS, (*v_data[iThresh].v_tMin_Ch4[iEn])  * 0.999, (*v_data[iThresh].v_tMax_Ch4[iEn]) * 1.001 }, TCh4_N )   );
//         v_data[iThresh].h_tCh13 .push_back(dEn.Histo1D({("Time_Res_Ch13_" + std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str(), ("Time_Res_Ch13_" + std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str(), NBINS, (*v_data[iThresh].v_tMin_Ch13[iEn]) * 0.999, (*v_data[iThresh].v_tMax_Ch13[iEn])* 1.001 }, TCh13_N )  );

//         // --- Trigger the loop
//         *v_data[iThresh].h_tCh13[iEn];
        
//       }
//       // <-- Loop over the energies <--
//       std::cout << std::endl;
//     }
//     // <-- Loop over the thresholds <--
//     std::cout << "\nDone.\n\n";













//     // --- Do the calculations and save the data.
//     std::cout << "Do the calculation and save the trees.\nThreshold\tEnergy\n";
//     for ( size_t iThresh = 0; iThresh < CFD_Col.size(); ++iThresh)
//     {
//       CFD_thresh = v_data[iThresh].CFD_Threshold;
//       std::cout << CFD_thresh << " ";
//       std::string CFD_suffix = "CFD_" + std::to_string(int(std::round(CFD_thresh*100))) + "perc";

//       // --- Define a Tree
//       TTree tree (("tree_"+CFD_suffix).c_str(), ("tree_"+CFD_suffix).c_str());
//       double measEn, sigmaEn, enRes, sigmaEnRes, tRes_Ch4, tRes_Ch13, sigmatRes_Ch4, sigmatRes_Ch13, sigma_sigmatRes_Ch4, sigma_sigmatRes_Ch13;
//       tree.Branch(  "CFD_thresh",         &CFD_thresh,          "CFD_thresh/D" );
//       tree.Branch(  "trueEn",             &trueEn,              "trueEn/D"   );
//       tree.Branch(  "measEnergy",         &measEn,              "measEnergy/D" );
//       tree.Branch(  "sigmaEn",            &sigmaEn,             "sigmaEn/D" );
//       tree.Branch(  "enRes",              &enRes,               "enRes/D"   );
//       tree.Branch(  "sigmaEnRes",         &sigmaEnRes,          "sigmaEnRes/D");
//       tree.Branch(  "tRes_4",             &tRes_Ch4,            "tRes_Ch4/D" );
//       tree.Branch(  "tRes_13",            &tRes_Ch13,           "tRes_Ch13/D" );
//       tree.Branch(  "sigmatRes_4",        &sigmatRes_Ch4,       "sigmatRes_Ch4/D" );
//       tree.Branch(  "sigmatRes_13",       &sigmatRes_Ch13,      "sigmatRes_Ch13/D" );
//       tree.Branch(  "sigma_sigmatRes_4",  &sigma_sigmatRes_Ch4, "sigma_sigmatRes_Ch4/D" );
//       tree.Branch(  "sigma_sigmatRes_13", &sigma_sigmatRes_Ch13,"sigma_sigmatRes_Ch13/D" );

//       // --- Change folder in the output .root file
//       std::string dirName = "Threshold_" + CFD_suffix + "/";
//       outFile.cd();
//       outFile.mkdir(dirName.c_str());
//       outFile.cd(dirName.c_str());

//       std::vector<double> v_measEn, v_sigmaEn, v_enRes, v_sigmaEnRes, v_tRes_Ch4, v_tRes_Ch13, v_zero,
//                           v_sigmatRes_Ch4, v_sigmatRes_Ch13, v_sigma_sigmatRes_Ch4, v_sigma_sigmatRes_Ch13;


//       // --> Loop over the energies -->
//       for ( size_t iEn = 0; iEn < v_data[iThresh].Energy_Col.size(); ++iEn)
//       {
//         // -- Save parameter in a tree --
//         trueEn              = v_data[iThresh].Energy_Col[iEn];
//         std::cout << trueEn << " ";

//         // measEn              = *v_data[iThresh].v_lTot_Mean [iEn];
//         // sigmaEn             = *v_data[iThresh].v_lTot_Sigma[iEn];
//         // enRes               = sigmaEn / measEn;
//         // sigmaEnRes          = enRes * sqrt(  pow( v_data[iThresh].h_lTot[iEn]->GetMeanError() / measEn ,2)  + pow( v_data[iThresh].h_lTot[iEn]->GetStdDevError() / sigmaEn ,2) ) ;
//         tRes_Ch4            = *v_data[iThresh].v_tCh4_Mean[iEn];
//         tRes_Ch13           = *v_data[iThresh].v_tCh13_Mean[iEn];
//         sigmatRes_Ch4       = *v_data[iThresh].v_tCh4_Sigma[iEn];
//         sigmatRes_Ch13      = *v_data[iThresh].v_tCh13_Sigma[iEn];
//         sigma_sigmatRes_Ch4 = v_data[iThresh].h_tCh4 [iEn] -> GetStdDevError();
//         sigma_sigmatRes_Ch13= v_data[iThresh].h_tCh13[iEn] -> GetStdDevError();



//         TF1 gaus      ("gaus", "gaus");
//         TF1 gaus_Ch4  ("gausCh4", "gaus");
//         TF1 gaus_Ch13 ("gausCh13", "gaus");
//         v_data[iThresh].h_lTot[iEn]  -> Fit("gaus"    , "Q");
//         v_data[iThresh].h_tCh4[iEn]  -> Fit("gausCh4" , "Q");
//         v_data[iThresh].h_tCh13[iEn] -> Fit("gausCh13", "Q");

//         measEn              = gaus.GetParameter(1);
//         sigmaEn             = gaus.GetParameter(2);
//         enRes               = sigmaEn / measEn;
//         sigmaEnRes          = enRes * sqrt(  pow( gaus.GetParError(1) / measEn ,2)  + pow( gaus.GetParError(2) / sigmaEn ,2) ) ;
//         // tRes_Ch4            = gaus_Ch4.GetParameter(1);
//         // tRes_Ch13           = gaus_Ch13.GetParameter(1);
//         // sigmatRes_Ch4       = gaus_Ch4.GetParameter(2);
//         // sigmatRes_Ch13      = gaus_Ch13.GetParameter(2);
//         // sigma_sigmatRes_Ch4 = gaus_Ch4.GetParError(2);
//         // sigma_sigmatRes_Ch13= gaus_Ch13.GetParError(2);

//         v_data[iThresh].h_lTot [iEn] ->  Write();
//         v_data[iThresh].h_lCh4 [iEn] ->  Write();
//         v_data[iThresh].h_lCh13[iEn] ->  Write();
//         v_data[iThresh].h_tCh4 [iEn] ->  Write();
//         v_data[iThresh].h_tCh13[iEn] ->  Write();

//         tree.Fill();


//         // --- Vectors for the plots. Seriously, who thought of GetVal()?!
//         v_zero                .emplace_back(0.                  );
//         v_measEn              .emplace_back(measEn              );
//         v_sigmaEn             .emplace_back(sigmaEn             );
//         v_enRes               .emplace_back(enRes               );
//         v_sigmaEnRes          .emplace_back(sigmaEnRes          );
//         v_tRes_Ch4            .emplace_back(tRes_Ch4            );
//         v_tRes_Ch13           .emplace_back(tRes_Ch13           );
//         v_sigmatRes_Ch4       .emplace_back(sigmatRes_Ch4       );
//         v_sigmatRes_Ch13      .emplace_back(sigmatRes_Ch13      );
//         v_sigma_sigmatRes_Ch4 .emplace_back(sigma_sigmatRes_Ch4 );
//         v_sigma_sigmatRes_Ch13.emplace_back(sigma_sigmatRes_Ch13);

//       }
//       // <-- Loop over the energies <--

//     tree.Write();
//     std::cout << std::endl;




//     TF1 fitEnRes ("fitEnRes", "sqrt( ([0]/sqrt(x))**2 + ([1])**2  )");
//     fitEnRes.SetParameter(0, 0.1);
//     fitEnRes.SetParameter(1, 0.01);

//     // --- Energy resolution
//     TCanvas c_EnRes ("EnRes", "EnRes", 800, 600);
//     TGraphErrors g_EnRes (v_data[iThresh].Energy_Col.size(), &v_data[iThresh].Energy_Col[0], &v_enRes[0], &v_zero[0], &v_sigmaEnRes[0] );
//     g_EnRes.SetMarkerStyle(20);
//     g_EnRes.SetMarkerColor(kRed+2);
//     g_EnRes.SetMarkerSize (1.3);
  
//     g_EnRes.Fit("fitEnRes");
//     g_EnRes.Draw("AP");
//     g_EnRes.Write("g_EnRes");
//     c_EnRes.SetGrid();
//     c_EnRes.Write();



//     TCanvas c_tRes_Ch4 ("tRes_Ch4", "tRes_Ch4", 800, 600);
//     TGraphErrors g_tRes_Ch4 (v_data[iThresh].Energy_Col.size(), &v_data[iThresh].Energy_Col[0], &v_sigmatRes_Ch4[0], &v_zero[0], &v_sigma_sigmatRes_Ch4[0] );
//     g_tRes_Ch4.SetMarkerStyle(20);
//     g_tRes_Ch4.SetMarkerColor(kBlue+2);
//     g_tRes_Ch4.SetMarkerSize (1.3);

//     g_tRes_Ch4.Draw("AP");
//     g_tRes_Ch4.Write("g_tRes4");
//     c_tRes_Ch4.SetGrid();
//     c_tRes_Ch4.Write();


//     TCanvas c_tRes_Ch13 ("tRes_Ch13", "tRes_Ch13", 800, 600);
//     TGraphErrors g_tRes_Ch13 (v_data[iThresh].Energy_Col.size(), &v_data[iThresh].Energy_Col[0], &v_sigmatRes_Ch13[0], &v_zero[0], &v_sigma_sigmatRes_Ch13[0] );
//     g_tRes_Ch13.SetMarkerStyle(20);
//     g_tRes_Ch13.SetMarkerColor(kGreen+2);
//     g_tRes_Ch13.SetMarkerSize (1.3);

//     g_tRes_Ch13.Draw("AP");
//     g_tRes_Ch13.Write("g_tRes13");
//     c_tRes_Ch13.SetGrid();
//     c_tRes_Ch13.Write();



//     TCanvas c_tRes_mg ("tRes_mg", "tRes_mg", 800, 600);
//     TMultiGraph mg_tRes;
//     mg_tRes.Add(&g_tRes_Ch4);
//     mg_tRes.Add(&g_tRes_Ch13);

//     mg_tRes.Write("g_tRes_mg");
//     mg_tRes.Draw("AP");

//     auto legend = new TLegend(0.65,0.7,0.85,0.8);
//     legend->AddEntry(&g_tRes_Ch4 ,"Front","ep");
//     legend->AddEntry(&g_tRes_Ch13,"Back","ep");
//     legend->Draw();

//     c_tRes_mg.SetGrid();
//     c_tRes_mg.Write();


//     }
//     // <-- Loop over the thresholds <--
//     std::cout << "\nDone.\n\n";




  

    































    
    // for ( auto && threshold : CFD_Col)
    // {
    //   double trueEn, measEn, sigmaEn, enRes, sigmaEnRes, tRes_Ch4, tRes_Ch13, sigmatRes_Ch4, sigmatRes_Ch13, sigma_sigmatRes_Ch4, sigma_sigmatRes_Ch13, CFD_thresh;

    //   auto FilterThreshold  = [&CFD_thresh] (const float & thresh) ->float{ return CFD_thresh == thresh;};
    //   auto FilterEnergy     = [&trueEn]     (const float & energy) ->float{ return trueEn     == energy;};

    //   CFD_thresh = threshold;
    //   std::cout << CFD_thresh << " ";

    //   std::string CFD_suffix = "CFD_" + std::to_string(int(std::round(CFD_thresh*100))) + "perc";

    //   // --- Define a Tree
    //   TTree tree (("tree_"+CFD_suffix).c_str(), ("tree_"+CFD_suffix).c_str());
    //   // std::vector<double> trueEn_V, en_V, sigmaEn_V, enRes_V, tRes_Ch4_V, tRes_Ch13_V, sigmatRes_Ch4_V, sigmatRes_Ch13_V, CFD_thresh_V;
    //   // std::vector<double> zero_V, sigmaEnRes_V, sigma_sigmatRes_Ch4_V, sigma_sigmatRes_Ch13_V;
    //   tree.Branch(  "CFD_thresh",         &CFD_thresh,          "CFD_thresh/D" );
    //   tree.Branch(  "trueEn",             &trueEn,              "trueEn/D"   );
    //   tree.Branch(  "measEnergy",         &measEn,              "measEnergy/D" );
    //   tree.Branch(  "sigmaEn",            &sigmaEn,             "sigmaEn/D" );
    //   tree.Branch(  "enRes",              &enRes,               "enRes/D"   );
    //   tree.Branch(  "sigmaEnRes",         &sigmaEnRes,          "sigmaEnRes/D");
    //   tree.Branch(  "tRes_4",             &tRes_Ch4,            "tRes_Ch4/D" );
    //   tree.Branch(  "tRes_13",            &tRes_Ch13,           "tRes_Ch13/D" );
    //   tree.Branch(  "sigmatRes_4",        &sigmatRes_Ch4,       "sigmatRes_Ch4/D" );
    //   tree.Branch(  "sigmatRes_13",       &sigmatRes_Ch13,      "sigmatRes_Ch13/D" );
    //   tree.Branch(  "sigma_sigmatRes_4",  &sigma_sigmatRes_Ch4, "sigma_sigmatRes_Ch4/D" );
    //   tree.Branch(  "sigma_sigmatRes_13", &sigma_sigmatRes_Ch13,"sigma_sigmatRes_Ch13/D" );



    //   // --- Change folder in the output .root file
    //   std::string dirName = "Threshold_" + CFD_suffix;
    //   outFile.cd();
    //   outFile.mkdir(dirName.c_str());
    //   outFile.cd(dirName.c_str());





    //   // -- Filter on the threshold value
    //   auto dThresh = dFilt.Filter(FilterThreshold, {CFD_N});


    //   // --> Loop over the energies -->
    //   for (auto & energy : Energy_Col)
    //   {
    //     trueEn = energy;
    //     auto dEn = dThresh.Filter(FilterEnergy, {ENER_N});

    //     // --- Get statistical quantities, also to draw the histograms 
    //     auto lTotMax    = dEn.Max(TL_N);
    //     auto lTotMin    = dEn.Min(TL_N);
    //     auto lTot_Mean  = dEn.Mean(TL_N);
    //     auto lTot_Sigma = dEn.StdDev(TL_N);


    //     auto lMax_Ch4   = dEn.Max(LCh4_N);
    //     auto lMax_Ch13  = dEn.Max(LCh13_N);
    //     auto lMin_Ch4   = dEn.Min(LCh4_N);
    //     auto lMin_Ch13  = dEn.Min(LCh13_N);


    //     auto tMax_Ch4   = dEn.Max(TCh4_N);
    //     auto tMax_Ch13  = dEn.Max(TCh13_N);
    //     auto tMin_Ch4   = dEn.Min(TCh4_N);
    //     auto tMin_Ch13  = dEn.Min(TCh13_N);
    //     auto tCh4_Mean  = dEn.Mean(TCh4_N);
    //     auto tCh13_Mean = dEn.Mean(TCh13_N);
    //     auto tCh4_Sigma = dEn.StdDev(TCh4_N);
    //     auto tCh13_Sigma= dEn.StdDev(TCh13_N);



    //     // --- Build the histos
    //     auto h_lTot  = dEn.Histo1D({("Energy_H_"      + std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str(), ("Energy_H_"      + std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str() , 100, (*lTotMin) * 0.9, (*lTotMax) * 1.1 }, TL_N );
    //     auto h_lCh4  = dEn.Histo1D({("Photons_H_Ch4_" + std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str(), ("Photons_H_Ch4_" + std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str() , 100, (*lMin_Ch4) * 0.9, (*lMax_Ch4) * 1.1 }, LCh4_N );
    //     auto h_lCh13 = dEn.Histo1D({("Photons_H_Ch13_"+ std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str(), ("Photons_H_Ch13_"+ std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str() , 100, (*lMin_Ch13) * 0.9, (*lMax_Ch13) * 1.1 }, LCh13_N );
    //     auto h_tCh4  = dEn.Histo1D({("Time_Res_Ch4_"  + std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str(), ("Time_Res_Ch4_"  + std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str() , 100, (*tMin_Ch4) * 0.9, (*tMax_Ch4) * 1.1 }, TCh4_N );
    //     auto h_tCh13 = dEn.Histo1D({("Time_Res_Ch13_" + std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str(), ("Time_Res_Ch13_" + std::to_string(int(trueEn)) + "GeV_" + CFD_suffix).c_str() , 100, (*tMin_Ch13) * 0.9, (*tMax_Ch13) * 1.1 }, TCh13_N );

    //     TF1 gaus      ("gaus", "gaus");
    //     TF1 gaus_Ch4  ("gausCh4", "gaus");
    //     TF1 gaus_Ch13 ("gausCh13", "gaus");
    //     h_lTot  -> Fit("gaus"    , "Q");
    //     h_lCh4  -> Fit("gausCh4" , "Q");
    //     h_lCh13 -> Fit("gausCh13", "Q");


    //     // -> Save parameter in a tree and in vectors ->
    //     measEn              = gaus.GetParameter(1);
    //     sigmaEn             = gaus.GetParameter(2);
    //     enRes               = sigmaEn / measEn;
    //     sigmaEnRes          = enRes * sqrt(  pow( gaus.GetParError(1) / measEn ,2)  + pow( gaus.GetParError(2) / sigmaEn ,2) ) ;
    //     tRes_Ch4            = *tCh4_Mean;
    //     tRes_Ch13           = *tCh13_Mean;
    //     sigmatRes_Ch4       = * tCh4_Sigma;
    //     sigmatRes_Ch13      = * tCh13_Sigma;
    //     sigma_sigmatRes_Ch4 = h_tCh4  -> GetStdDevError();
    //     sigma_sigmatRes_Ch13= h_tCh13 -> GetStdDevError();


    //     h_lTot  ->  Write();
    //     h_lCh4  ->  Write();
    //     h_lCh13 ->  Write();
    //     h_tCh4  ->  Write();
    //     h_tCh13 ->  Write();

    //     tree.Fill();

    //   }
    //   // <-- Loop over the energies <--

    // tree.Write();




      // --- Update the map entry
      // mapPos -> second.zero_V          .emplace_back(0.);
      // mapPos -> second.trueEn_V        .emplace_back(trueEn);
      // mapPos -> second.measEn_V        .emplace_back(measEn);
      // mapPos -> second.sigmaEn_V       .emplace_back(sigmaEn);
      // mapPos -> second.enRes_V         .emplace_back(enRes);
      // mapPos -> second.sigmaEnRes_V    .emplace_back(sigmaEnRes);
      // mapPos -> second.tRes_Ch4_V      .emplace_back(tRes_Ch4);
      // mapPos -> second.tRes_Ch13_V     .emplace_back(tRes_Ch13);
      // mapPos -> second.sigmatRes_Ch4_V .emplace_back(sigmatRes_Ch4);
      // mapPos -> second.sigmatRes_Ch13_V.emplace_back(sigmatRes_Ch13);
      // mapPos -> second.sigma_sigmatRes_Ch4_V .emplace_back(sigma_sigmatRes_Ch4);
      // mapPos -> second.sigma_sigmatRes_Ch13_V.emplace_back(sigma_sigmatRes_Ch13);


 







  // // --> Now loop over the map keys to produce the plots -->
  // for (auto&& threshold : map)
  // {
  //   dirName = "Plots_CFD_Threshold_" + std::to_string(threshold.first) + "/";
  //   outFile.cd();
  //   outFile.mkdir(dirName.c_str());
  //   outFile.cd(dirName.c_str());

  //   TGraphErrors enResGraph ( threshold.second.trueEn_V.size(), &threshold.second.trueEn_V[0], &threshold.second.enRes_V[0],          &threshold.second.zero_V[0], &threshold.second.sigmaEnRes_V[0] );
  //   TGraphErrors tResG_Ch4  ( threshold.second.trueEn_V.size(), &threshold.second.trueEn_V[0], &threshold.second.sigmatRes_Ch4_V[0],  &threshold.second.zero_V[0], &threshold.second.sigma_sigmatRes_Ch4_V[0] );
  //   TGraphErrors tResG_Ch13 ( threshold.second.trueEn_V.size(), &threshold.second.trueEn_V[0], &threshold.second.sigmatRes_Ch13_V[0], &threshold.second.zero_V[0], &threshold.second.sigma_sigmatRes_Ch13_V[0] );

  //   TMultiGraph  tRes_mg; 

  //   TCanvas cEn (("Energy_Resolution_threshold_" + std::to_string(threshold.first)).c_str(), "Energy_Resolution", 800, 600);
  //   enResGraph.GetXaxis() -> SetTitle("Energy [GeV]");
  //   enResGraph.GetYaxis() -> SetTitle("#frac{#sigma}{Energy}");
  //   enResGraph.SetMarkerStyle(20);
  //   enResGraph.SetMarkerSize(1.3);
  //   enResGraph.SetMarkerColor(kBlue);
  //   enResGraph.Draw("AP");
  //   cEn.SetGrid();

  //   TCanvas ct_Ch4 (("Time_Resolution_Front_threshold_" + std::to_string(threshold.first)).c_str(), "Time_Resolution_Front", 800, 600);
  //   tResG_Ch4.GetXaxis() -> SetTitle("Energy [GeV]");
  //   tResG_Ch4.GetYaxis() -> SetTitle("#sigma_{t} [ns]");
  //   tResG_Ch4.SetMarkerStyle(20);
  //   tResG_Ch4.SetMarkerSize(1.5);
  //   tResG_Ch4.SetMarkerColor(kRed+1);
  //   tResG_Ch4.Draw("AP");
  //   ct_Ch4.SetGrid();

  //   TCanvas ct_Ch13 (("Time_Resolution_Back_threshold_" + std::to_string(threshold.first)).c_str(), "Time_Resolution_Back", 800, 600);
  //   tResG_Ch13.GetXaxis() -> SetTitle("Energy [GeV]");
  //   tResG_Ch13.GetYaxis() -> SetTitle("#sigma_{t} [ns]");
  //   tResG_Ch13.SetMarkerStyle(20);
  //   tResG_Ch13.SetMarkerSize(1.5);
  //   tResG_Ch13.SetMarkerColor(kGreen+3);
  //   tResG_Ch13.Draw("AP");
  //   ct_Ch13.SetGrid();


  //   tRes_mg.Add(&tResG_Ch4);
  //   tRes_mg.Add(&tResG_Ch13);
  //   TCanvas ct_Chs (("Time_Resolutions_threshold_" + std::to_string(threshold.first)).c_str(), "Time_Resolutions", 800, 600);
  //   tRes_mg.SetTitle("Time Resolution - Central Cells");
  //   tRes_mg.GetXaxis() -> SetTitle("Energy [GeV]");
  //   tRes_mg.GetYaxis() -> SetTitle("#sigma_{t} [ns]");
  //   tRes_mg.GetYaxis() -> SetRangeUser (0., 1);
  //   tRes_mg.GetXaxis() -> SetRangeUser (0., 120.);
  //   tRes_mg.Draw("AP");
  //   ct_Chs.SetGrid();

  //   auto legend = new TLegend(0.65,0.7,0.85,0.8);
  //   legend->AddEntry(&tResG_Ch4 ,"Front","ep");
  //   legend->AddEntry(&tResG_Ch13,"Back","ep");
  //   legend->Draw();






  // TF1 fitEnRes ("fitEnRes", "sqrt( ([0]/sqrt(x))**2 + ([1])**2  )");
  // fitEnRes.SetParameter(0, 0.1);
  // fitEnRes.SetParameter(1, 0.01);
  // enResGraph.Fit("fitEnRes");



  // // enResGraph.SetDirectory(gDirectory);
  // // tResG_Ch4 .SetDirectory(gDirectory);
  // // tResG_Ch13.SetDirectory(gDirectory);
  // // tRes_mg   .SetDirectory(gDirectory);
  // enResGraph.Write(("Energy_Resolution_Graph_threshold_" + std::to_string(threshold.first)).c_str());
  // tResG_Ch4 .Write(("Time_Resolution_Ch4_Graph_threshold_"  + std::to_string(threshold.first)).c_str());
  // tResG_Ch13.Write(("Time_Resolution_Ch13_Graph_threshold_" + std::to_string(threshold.first)).c_str());
  // tRes_mg.Write();

  // // cEn    .SetDirectory(gDirectory);
  // // ct_Ch4 .SetDirectory(gDirectory);
  // // ct_Ch13.SetDirectory(gDirectory);
  // // ct_Chs .SetDirectory(gDirectory);
  // cEn    .Write();
  // ct_Ch4 .Write();
  // ct_Ch13.Write();
  // ct_Chs .Write();
  // }

  outFile.Close();
  std::cout << "\n";
  std::cout << "Saved to: " << (path + outFileName + ".root").c_str() << "\n";
  std::cout << "### Duration: " << (std::clock() - _start) / (double) CLOCKS_PER_SEC << "s\n";
  std::cout << "### Standard end of the program. ###" << std::endl;
  gApplication->Terminate();
}
