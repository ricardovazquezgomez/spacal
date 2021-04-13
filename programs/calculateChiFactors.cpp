// compile with
// g++ -o ../build/calculateChiFactors calculateChiFactors.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore

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
#include "../parametrization/regions.h"


// info about the command line usage
void usage()
{
  std::cout << std::endl;
  std::cout << "\t\t" << "[-e | --efile]   <input energy deposition file name>   " << std::endl
            << "\t\t" << "[-p | --pfile]   <input pulse file name>   " << std::endl
            << "\t\t" << "[-o | --output]  <output file name>   " << std::endl
            << "\t\t" << std::endl
            << "\t\t" << "Restricted analysis: you can specify either a module, or a module plus a front and back cell numbers" << std::endl
            << "\t\t" << "[--module]       <restric analysis to a specific module number>   " << std::endl
            << "\t\t" << "[--front]        <front cell for restricted analysis>   " << std::endl
            << "\t\t" << "[--back]         <back cell for restricted analysis>   " << std::endl
            // << "\t\t" << "[--energy]       <energy of incident beam [GeV] - mandatory >  " << std::endl
            // << "\t\t" << "[--pdgid]        <PDG ID of incident beam - mandatory >  " << std::endl
            // << "\t\t" << "[--type]         <csv list of module types - default = 0 >  " << std::endl
            // << "\t\t" << "[--module]       <csv list of module numbers in sim output world, for each module type - default = 0 >  " << std::endl
            // << "\t\t" << "[--cells]        <csv list of number of cells per module type - default = 2 >  " << std::endl
            // << "\t\t" << "[--lastFront]    <csv list of last front cell for each module type - default = 0 >  " << std::endl
            // << "\t\t" << "[--fraction]     <csv list of sampling fractions per module type - default = 1>  " << std::endl
            // << "\t\t" << "[--caliFront]    <csv list of front calibration Ph.El./MeV factors - default = 1>  " << std::endl
            // << "\t\t" << "[--caliBack]     <csv list of back calibration Ph.El./MeV factors - default = 1>  " << std::endl
            << "\t\t" << "[--verbose]      <verbosity level - default = 0>  " << std::endl
            // << "\t\t" << "[--fl]           <FL factor in pulse formation - default = 1>  " << std::endl

            << "\t\t" << std::endl;
}


// void read_directory(const std::string& name, std::vector<std::string> &v)
// {
//   // #######################################
//   // ### list files in directory
//   // ### taken from
//   // ### http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html
//     DIR* dirp = opendir(name.c_str());
//     struct dirent * dp;
//     while ((dp = readdir(dirp)) != NULL) {
//         v.push_back(dp->d_name);
//     }
//     closedir(dirp);
// }

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
  float *sigma;
  float *offset;
  // histograms
  TH1F** h_l_cell;
  TH1F** h_phel_mev_cell;
  TH1F** h_t_cell;
  TH1F* h_t_combined;
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


struct event_t
{
  float primaryEnergy;
  int   primaryPDGID;
  float endep_total;     // total energy deposited in module 
  float endep_front;     // energy in the front section
  float endep_back;      // energy in back section
  float endep_cry_total; // total energy deposited in crystals
  float endep_cry_front; // energy in the front crystals
  float endep_cry_back;  // energy in back crystals
  int   light_front;
  int   light_back;
};


int main (int argc, char **argv)
{
  

  if(argc < 2)
  {
    std::cout << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }
  //-------------------------------------------//
  // DICTIONARIES                              //
  //-------------------------------------------//
  // --- FAST DICTIONARIES COMPILATION ON THE FLY
  gInterpreter -> GenerateDictionary("std::vector<std::vector<float>>", "vector");
  gInterpreter -> GenerateDictionary("std::vector<std::vector<std::vector<float>>>", "vector");
  // gStyle->SetOptStat(0000);
  //-------------------------------------------//





  //-------------------------------------------//
  // INPUT PARAMETERS                          //
  //-------------------------------------------//
  std::string energyFileName = "";
  std::string pulseFileName = "";
  std::string outputFileName = "";
  int verbose = 0;
  int module = -1;
  int front = -1;
  int back = -1;

  static struct option longOptions[] =
  {
			{ "efile", required_argument, 0, 0 },
      { "pfile", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { "verbose", required_argument, 0, 0 },
      { "module", required_argument, 0, 0 },
      { "front", required_argument, 0, 0 },
      { "back", required_argument, 0, 0 },
			{ NULL, 0, 0, 0 }
	};
  // read input
  while(1) {
		int optionIndex = 0;
		int c = getopt_long(argc, argv, "e:p:o:", longOptions, &optionIndex);
		if (c == -1) {
			break;
		}
		if (c == 'e'){
			energyFileName= (char *)optarg;
    }
		else if (c == 'p'){
      pulseFileName= (char *)optarg;
    }
    else if (c == 'o'){
      outputFileName = (char *)optarg;
    }
		else if (c == 0 && optionIndex == 0){
      energyFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 1){
      pulseFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 2){
      outputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 3){
      verbose = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 4){
      module = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 5){
      front = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 6){
      back = atoi((char *)optarg);
    }
		else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
	}
  // check mandatory input
  if(energyFileName == "")
  {
    std::cout << std::endl;
    std::cout << "ERROR! You need to provide an energy file name!" << std::endl;
    std::cout << "See program usage below..." << std::endl;
    std::cout << std::endl;
    std::cout << argv[0];
    usage();
    return 1;
  }
  if(pulseFileName == "")
  {
    std::cout << std::endl;
    std::cout << "ERROR! You need to provide a pulse file name!" << std::endl;
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
  //-------------------------------------------//


  //-------------------------------------------//
  // RESTRICTED ANALYSIS                       //
  //-------------------------------------------//
  // understand if analysis is restricted
  bool restricted = false;
  if(module > -1)
  {
    restricted = true;
  }
  if(front > -1)
  {
    if(!restricted)
    {
      std::cout << std::endl;
      std::cout << "ERROR! If you specify a front cell number, you also need to set a module number!" << std::endl;
      std::cout << "See program usage below..." << std::endl;
      std::cout << std::endl;
      std::cout << argv[0];
      usage();
      return 1;
    }

    if(back == -1)
    {
      std::cout << std::endl;
      std::cout << "ERROR! If you specify a front cell number, you also need to provide a back cell number!" << std::endl;
      std::cout << "See program usage below..." << std::endl;
      std::cout << std::endl;
      std::cout << argv[0];
      usage();
      return 1;
    }
  }
  if(back > -1)
  {
    if(!restricted)
    {
      std::cout << std::endl;
      std::cout << "ERROR! If you specify a front cell number, you also need to set a module number!" << std::endl;
      std::cout << "See program usage below..." << std::endl;
      std::cout << std::endl;
      std::cout << argv[0];
      usage();
      return 1;
    }

    if(front == -1)
    {
      std::cout << std::endl;
      std::cout << "ERROR! If you specify a back cell number, you also need to provide a front cell number!" << std::endl;
      std::cout << "See program usage below..." << std::endl;
      std::cout << std::endl;
      std::cout << argv[0];
      usage();
      return 1;
    }
  }



  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");


  //-------------------------------------------//
  // MAPS                                      //
  //-------------------------------------------//
  maps_t maps;
  TFile *inputFile = TFile::Open(energyFileName.c_str());
  int last_front = readStructure(inputFile,maps);
  //-------------------------------------------//

  





  
  //-------------------------------------------//
  // ENERGY DEPOSITIONS                        //
  //-------------------------------------------//
  
  // prepare vector 
  std::vector<event_t> event;
  // get ttree 
  TTree *shower = ((TTree*) 0);
  shower = (TTree*) inputFile->Get("shower");
  // read input shower
  //----------------------------------//
  // INPUT SHOWER                     //
  //----------------------------------//
  int              showerRun              ;
  int              showerEvent            ;
  int              showerIsInCrystal      ;
  int              showerModuleType       ;
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
  TBranch *b_showerModuleType       ;
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
  shower->SetBranchAddress("moduleType"         ,&showerModuleType      ,&b_showerModuleType      );
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
  // //------------------------------------------------//


  int eventID = -1;
  long int nEntriesShower = shower->GetEntries();
  std::cout << "Energy deposition events = " << nEntriesShower <<  std::endl;
  std::cout << "Running on energy depositions..." << std::endl;
  long int counter = 0;
  int lostScintillation = 0;
  int lostShower = 0;
  float totalEnergyPerEvent = 0;
  float totalEnergyPerEventFront = 0;
  float totalEnergyPerEventBack = 0;
  float totalCrystalEnergyPerEvent = 0;
  float totalCrystalEnergyPerEventFront = 0;
  float totalCrystalEnergyPerEventBack = 0;
  int loop_primaryPDGID = 0;
  float loop_primaryEnergy = 0.;
  for(int iEntry = 0 ; iEntry < nEntriesShower ; iEntry++)
  {
    shower->GetEvent(iEntry);
    // summation of energy needs to be performed only until a new event arrives in the list
    // an event corresponds to 1 primary shot. the info is stored in event
    // so we start be setting a event indicator, eventID, to -1. At the first loop step,
    // the first value of variable event will be 0, so what comes in the next if statement
    // (event saving and clean up) will not be executed, and instead the event accumulation starts
    // and the end of event accumulation the eventID becomes = event, then a new step of the cycle
    // starts. Now eventID is = 0, so "event saving and clean up" will not be performed until
    // the value of event becomes != 0. When that happens, the variables were the enegies were
    // summed are used to filled the relevant histograms, than they are cleared and a new accumulation
    // starts. Of course this means that after the end of the loop, the last event still needs to be saved
    //
    //
    //------------------------------//
    // event saving and clean up    //
    //------------------------------//
    if((eventID != showerEvent) && (eventID != -1)) // new event, dump the previous
    {
      // prepare event 
      event_t temp_event;
      temp_event.primaryPDGID = loop_primaryPDGID;
      temp_event.primaryEnergy = loop_primaryEnergy;
      temp_event.endep_total = totalEnergyPerEvent;
      temp_event.endep_front = totalEnergyPerEventFront;
      temp_event.endep_back  = totalEnergyPerEventBack;
      temp_event.endep_cry_total = totalCrystalEnergyPerEvent;
      temp_event.endep_cry_front = totalCrystalEnergyPerEventFront;
      temp_event.endep_cry_back  = totalCrystalEnergyPerEventBack;

      // fill
      event.push_back(temp_event);

      // cleanup
      totalEnergyPerEvent = 0;
      totalEnergyPerEventFront = 0;
      totalEnergyPerEventBack = 0;
      totalCrystalEnergyPerEvent = 0;
      totalCrystalEnergyPerEventFront = 0;
      totalCrystalEnergyPerEventBack = 0;
    }

    // skip accumulation, if restricted and conditions not met
    if(restricted)
    {
      // check module number 
      Int_t binx   = maps.module_ID_map->GetXaxis()->FindBin(showerX);
      Int_t biny   = maps.module_ID_map->GetYaxis()->FindBin(showerY);
      Float_t xMod = maps.module_ID_map->GetXaxis()->GetBinCenter(binx);
      Float_t yMod = maps.module_ID_map->GetYaxis()->GetBinCenter(biny);
      // module number
      int module_number = maps.module_ID_map->GetBinContent(binx,biny);
      if(module_number != module) 
      {
        if(verbose > 0)
        {
          std::cout << "WARNING: Skipping accumulation of event because of --module choice " << std::endl;
          std::cout << "         Energy deposition module number " << module_number << " does not correspond to input module " << module << std::endl;
          std::cout << "         This event will be ignored! " << std::endl;
          std::cout << std::endl;
        }
        continue;
      }
      // check front and back cell numbers
      if( (front > -1) && (back > -1) )
      {
        // fetch cell map index
        int cellMapIndex;
        std::map<int, int> ::const_iterator iterCellMap = maps.mapOfCellMaps.find(showerModuleType);
        if(iterCellMap != maps.mapOfCellMaps.end())
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
          continue;
        }

        float moduleSeparation = 0.;
        //fetch separation point in z
        std::map<int, float> ::const_iterator iterModSep = maps.moduleSeparationMap.find(showerModuleType);
        if(iterModSep != maps.moduleSeparationMap.end())
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
          //lostScintillation++;
          continue;
        }
        // find sIndex 
        int sIndex = -1;
        if(showerZ < moduleSeparation) // then it's front
        {
          sIndex = 0;
        }
        else
        {
          sIndex = 1;
        }

        // check if cell is front or back 
        Int_t c_binx    = maps.cell_map_per_module[cellMapIndex][sIndex]->GetXaxis()->FindBin(showerX - xMod);
        Int_t c_biny    = maps.cell_map_per_module[cellMapIndex][sIndex]->GetYaxis()->FindBin(showerY - yMod);
          // cell number
        int cell_number = maps.cell_map_per_module[cellMapIndex][sIndex]->GetBinContent(c_binx,c_biny);

        if((cell_number == front) || (cell_number == back))
        {
          // keep event, so do nothing here
        }
        else 
        {
          if(verbose > 0)
          {
            std::cout << "WARNING: Skipping accumulation of event because of --front or --back choice " << std::endl;
            std::cout << "         Energy deposition cell number " << cell_number << " does not corrispond to either front (" << front << ") or back (" << back << ")" << std::endl;
            std::cout << "         This event will be ignored! " << std::endl;
            std::cout << std::endl;
          }
          continue;
        }

      }
    }
    
    float moduleSeparation = 0.;
    //fetch separation point in z
    std::map<int, float> ::const_iterator iterModSep = maps.moduleSeparationMap.find(showerModuleType);
    if(iterModSep != maps.moduleSeparationMap.end())
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
      //lostScintillation++;
      continue;
    }

    // accumulate total energy
    totalEnergyPerEvent += showerTotalEnDep;
    if(showerZ < moduleSeparation) // then it's front
    {
      totalEnergyPerEventFront += showerTotalEnDep;
    }
    else
    {
      totalEnergyPerEventBack += showerTotalEnDep;
    }
    if(showerIsInCrystal)
    {
      totalCrystalEnergyPerEvent += showerTotalEnDep;
      if(showerZ < moduleSeparation) // then it's front
      {
        totalCrystalEnergyPerEventFront += showerTotalEnDep;
      }
      else
      {
        totalCrystalEnergyPerEventBack += showerTotalEnDep;
      }
    }
    
    // update event ID
    eventID = showerEvent;
    loop_primaryPDGID = shower_primaryPDGID;
    loop_primaryEnergy = shower_primaryEnergy;


  }

// prepare last event 
  event_t temp_event;
  temp_event.primaryPDGID = loop_primaryPDGID;
  temp_event.primaryEnergy = loop_primaryEnergy;
  temp_event.endep_total = totalEnergyPerEvent;
  temp_event.endep_front = totalEnergyPerEventFront;
  temp_event.endep_back  = totalEnergyPerEventBack;
  temp_event.endep_cry_total = totalCrystalEnergyPerEvent;
  temp_event.endep_cry_front = totalCrystalEnergyPerEventFront;
  temp_event.endep_cry_back  = totalCrystalEnergyPerEventBack;
  // fill
  event.push_back(temp_event);
  // cleanup
  totalEnergyPerEvent = 0;
  totalEnergyPerEventFront = 0;
  totalEnergyPerEventBack = 0;
  totalCrystalEnergyPerEvent = 0;
  totalCrystalEnergyPerEventFront = 0;
  totalCrystalEnergyPerEventBack = 0;

  inputFile->Close();
  std::cout << "...done." << std::endl; 


  //-------------------------------------------//
  // LIGHT OUTPUT                              //
  //-------------------------------------------//
  

  TFile f (pulseFileName.c_str(), "READ");
  auto TT = (TTree*) f.Get("tree");
  std::vector<int>   * VPh = 0;
  std::vector<int>   * VMod = 0;
  std::vector<float> * VTime = 0;
  int totalLight = 0;
  
  std::string moduleName = "mod0";
  if(restricted)
  {
    std::stringstream smodule;
    smodule << "mod" << module;
    moduleName = smodule.str();
  }

  TT->SetBranchAddress("Total_Light", &totalLight);
  TT->SetBranchAddress("modulesHit", &VMod);
  TT->SetBranchAddress((moduleName + "_ph").c_str(), &VPh);
  TT->SetBranchAddress((moduleName + "_t").c_str(), &VTime);

  auto entries = TT->GetEntries();
  for(int iEvent = 0 ; iEvent < entries; iEvent++)
  {
    TT->GetEntry(iEvent);
    //loop on cells
    float light_front = 0.;
    float light_back  = 0.;
    for (size_t iCell = 0; iCell < VPh->size(); ++iCell)
    {
      // skip if restricted and conditions not met 

      if(restricted)
      {
        if((iCell == front) || (iCell == back))
        {
          // keep so do nothing
        }
        else 
        {
          if(verbose > 0)
          {
            std::cout << "WARNING: Skipping ph.el of event because of --front or --back choice " << std::endl;
            std::cout << "         ph.el. cell number " << iCell << " does not corrispond to either front (" << front << ") or back (" << back << ")" << std::endl;
            std::cout << "         This event will be ignored! " << std::endl;
            std::cout << std::endl;
          }
          continue;
        }
      }
      if(iCell > last_front)
      {
        light_back += VPh->at(iCell);
      }
      else 
      {
        light_front += VPh->at(iCell);
      } 
    }
    event[iEvent].light_front = light_front;
    event[iEvent].light_back  = light_back ;
  }
  f.Close();






  // fill a TTree 
  TTree* chi = new TTree ("chi","chi") ;
  Double_t primaryEnergy  ;
  Int_t    primaryPDGID   ;
  Double_t endep_total    ;
  Double_t endep_front    ; 
  Double_t endep_back     ;
  Double_t endep_cry_total; 
  Double_t endep_cry_front; 
  Double_t endep_cry_back ;
  Int_t    light_front    ; 
  Int_t    light_back     ;
  chi->Branch("primaryEnergy_mev"   , &primaryEnergy   , "primaryEnergy_mev/D");
  chi->Branch("primaryPDGID"    , &primaryPDGID    , "primaryPDGID/I");
  chi->Branch("endep_total_mev"     , &endep_total     , "endep_total_mev/D");
  chi->Branch("endep_front_mev"     , &endep_front     , "endep_front_mev/D");
  chi->Branch("endep_back_mev"      , &endep_back      , "endep_back_mev/D");
  chi->Branch("endep_cry_total_mev" , &endep_cry_total , "endep_cry_total_mev/D");
  chi->Branch("endep_cry_front_mev" , &endep_cry_front , "endep_cry_front_mev/D");
  chi->Branch("endep_cry_back_mev" , &endep_cry_back  , "endep_cry_back_mev/D");
  chi->Branch("photoElectrons_front"     , &light_front     , "photoElectrons_front/I");
  chi->Branch("photoElectrons_back"      , &light_back      , "photoElectrons_back/I");

  for(unsigned int i = 0 ; i < event.size(); i++)
  {
    if(verbose > 0)
    {
      std::cout << i << "\t" 
                << event[i].primaryEnergy   << "\t" 
                << event[i].primaryPDGID    << "\t" 
                << event[i].endep_total     << "\t" 
                << event[i].endep_front     << "\t" 
                << event[i].endep_back      << "\t"
                << event[i].endep_cry_total << "\t" 
                << event[i].endep_cry_front << "\t" 
                << event[i].endep_cry_back  << "\t"
                << event[i].light_front     << "\t" 
                << event[i].light_back      << "\t"
                << std::endl; 
    }

    primaryEnergy   = event[i].primaryEnergy   ;
    primaryPDGID    = event[i].primaryPDGID    ;
    endep_total     = event[i].endep_total     ;
    endep_front     = event[i].endep_front     ;
    endep_back      = event[i].endep_back      ;
    endep_cry_total = event[i].endep_cry_total ;
    endep_cry_front = event[i].endep_cry_front ;
    endep_cry_back  = event[i].endep_cry_back  ;
    light_front     = event[i].light_front     ;
    light_back      = event[i].light_back      ;

    chi->Fill();
  }


  std::cout << "Writing output to file " << outputFileName << std::endl;
  outputFile->cd();
  chi->Write();
  outputFile->Close();
  std::cout << "Done. Goobye. " << std::endl;
  
  return 0;

}
 