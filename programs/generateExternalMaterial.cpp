// Marco Pizzichemi 04.06.2020 marco.pizzichemi@cern.ch

// compile with
// g++ -o ../build/generateExternalMaterial generateExternalMaterial.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore


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
#include "TGraphDelaunay.h"
#include "TVector.h"

#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <getopt.h>
#include <limits>
#include <dirent.h>
#include "TEnv.h"

// #include "../parametrization/regions.h"
//

// struct cell_t
// {
//   bool save;
//   TH1F *h_EnDepInCrystal   ;
//   // TH1F *h_EnDepNotInCrystal;
//   // TH1F *h_EnDepTotal       ;
// };

// struct module_t
// {
//   bool save;
//   int type;
//   int sections;
//   int *nCellsPerSections;
//   TH1F *h_EnDepInCrystal   ;
//   TH1F *h_EnDepNotInCrystal;
//   TH1F *h_EnDepTotal       ;
//   bool **saveCell;
//   TH1F ***h_Cell_EnDepInCrystal   ;
//   // cell_t **cell;
// };

void usage()
{
  std::cout << "\t\t" << "[-f | --folder]  <input folder>   " << std::endl
            << "\t\t" << "[-i | --input]   <input file prefix>   " << std::endl
            << "\t\t" << "[-o | --output]  <output file name>  " << std::endl
            << "\t\t" << "[--verbose]      <verbosity level - default = 0>  " << std::endl
            // << "\t\t" << "[--max]          <max of energy histograms [MeV] - default = 100000>  " << std::endl
            // << "\t\t" << "[--bins]         <bins of energy histograms      - default = 1000>  " << std::endl
            // << "\t\t" << "Choice of level of energy summing (at least 1, can be all 3):" << std::endl
            // << "\t\t" << "[--abs]          sum energy on each absorber   " << std::endl
            // << "\t\t" << "[--cell]         sum energy on each cell       " << std::endl
            // << "\t\t" << "[--cry]          sum energy on each crystal    " << std::endl
            // << "\t\t" << "[--debug]        output some debug    " << std::endl
            // << "\t\t" << "[-c | --cell]    <cell number>       " << std::endl
            << "\t\t" << std::endl;
}



int main(int argc, char** argv)
{

  //------------------------------------------------//
  // CHECK CMD LINE INPUT                           //
  //------------------------------------------------//
  if(argc < 2)
  {
    std::cout << argv[0] << std::endl;
    usage();
    return 1;
  }
  // end of CHECK CMD LINE INPUT
  //------------------------------------------------//

  // gEnv->GetValue("TFile.Recover", 0);



  //------------------------------------------------//
  // PARSE INPUT                                    //
  //------------------------------------------------//
  // prepare variable, giving defaults
  std::string inputFilePrefix = "";
  std::string inputFolderName = "./";
  std::string outputFileName = "";
  int verbose = 0;
  static struct option longOptions[] =
  {
			{ "rise", required_argument, 0, 0 },
      { "decay", required_argument, 0, 0 },
      { "rise", required_argument, 0, 0 },
      { "rise", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      // { "folder", required_argument, 0, 0 },
      // { "verbose", required_argument, 0, 0 },
      // { "min", required_argument, 0, 0 },
      // { "folder", no_argument, 0, 0 },
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


  //----------------------------------//
  // INPUT FILES                   //
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
  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");


  //------------------------------------------------//
  // CALO STRUCTURE                                 //
  //------------------------------------------------//
  // calo structures from first file
  // Need to take info on
  // how many modules
  // their phyiscal position and dimensions
  // how many cells per modules
  // their phyiscal position and dimensions

  // open first file (we suppose all the input file corresponds to the same simulation, otherwise crazy user...)
  TFile *_file0 = TFile::Open(listInputFiles[0].c_str());
  // get structure TTrees
  TTree *cells = (TTree*) _file0->Get("cells");
  TTree *modules = (TTree*) _file0->Get("modules");
  // calculate regions (volumes of space corresponding to...)
  std::vector<region_t> cells_regions   = CalculateRegions(cells);
  std::vector<region_t> modules_regions = CalculateRegions(modules);

  // produce the module ID TH2I map
  // this TH2I has 1 pixel per module, and the bin limits are the module limits (in x and y)
  int det_counter = 0;
  outputFile->cd();
  TH2I *module_ID_map = ComputeElementMap(modules,"module_ID_map",-1,-1,det_counter);


  // std::cout << module_ID_map << std::endl;

  // get info on the module types, sections, separations
  std::vector<int> moduleTypes;
  std::vector<int> moduleSections;
  std::vector<float> moduleSeparationZ;
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
    if(!moduleTypeIsAlreadyThere) // new module type found, add type, sections, separation Z to lists
    {
      moduleTypes.push_back(modules_regions[iMod].type);
      // also get sections
      moduleSections.push_back(modules_regions[iMod].sections);
      moduleSeparationZ.push_back(modules_regions[iMod].separation_z);
    }
  }


  // lots of maps
  std::map<int, int>   moduleNumberMap;         // map from module type to number of modules
  std::map<int, int>   moduleSectionsMap;       // map from module type to number of segments
  std::map<int, float> moduleMinMap;            // map from module type to min z
  std::map<int, float> moduleMaxMap;            // map from module type to max z
  std::map<int, float> moduleSeparationMap;     // map from module type to z of separation
  std::map<int, int>   mapOfCellMaps;           // map from showerModuleType to i-index of cell_map_per_module (then j is 0 or 1)
  std::map<int, int>   mapOfFrontCellPerModuleType;           // map from module type to number of front cells
  std::map<int, int>   mapOfBackCellPerModuleType;           // map from module type to number of back cells

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
  // int totModules = 0;
  // for (std::map<int,int>::iterator it=moduleNumberMap.begin(); it!=moduleNumberMap.end(); ++it)
  // {
  //   moduleNumberOfBeforeMap.insert(std::make_pair(it->first,totModules));
  //   totModules += it->second;
  // }
  int Nmodules = modules_regions.size();
  std::cout << "Total number of modules                 = "<< Nmodules << std::endl;
  int NofModuleTypes    = moduleTypes.size();
  std::cout << "Number of module types                 = "<< NofModuleTypes << std::endl;

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
      if(j == 0)
      {
        mapOfFrontCellPerModuleType.insert(std::make_pair(moduleTypes[i],cell_map_per_module[i][j]->GetEntries()));
      }
      else
      {
        mapOfBackCellPerModuleType.insert(std::make_pair(moduleTypes[i],cell_map_per_module[i][j]->GetEntries()));
      }
      mapOfCellMaps.insert(std::make_pair(moduleTypes[i],i));
    }
  }

  // vector of found
  std::vector< std::vector<bool> > foundLimit ;
  for(int i = 0 ; i < moduleTypes.size(); i++)
  {
    moduleSectionsMap.insert(std::make_pair(moduleTypes[i],moduleSections[i]));
    moduleSeparationMap.insert(std::make_pair(moduleTypes[i],moduleSeparationZ[i]));
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
          moduleMaxMap.insert(std::make_pair(moduleTypes[i],cells_regions[iCell].max[2]));
          foundLimit[i][0] = true;
          moduleMinMap.insert(std::make_pair(moduleTypes[i],cells_regions[iCell].min[2]));
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
              moduleMaxMap.insert(std::make_pair(moduleTypes[i],cells_regions[iCell].max[2]));
            }
          }
          else
          {
            if(foundLimit[i][0] == false)
            {
              foundLimit[i][0] = true;
              moduleMinMap.insert(std::make_pair(moduleTypes[i],cells_regions[iCell].min[2]));
            }
          }
        }

      }
    }
  }

  for(int i = 0 ; i < moduleTypes.size(); i++)
  {
    std::cout << "Type " << moduleTypes[i] << " has " << moduleSectionsMap.at(moduleTypes[i]) << " sections " << std::endl;
    std::cout << "Type " << moduleTypes[i] << " has " << moduleNumberMap.at(moduleTypes[i]) << " modules " << std::endl;
    // std::cout << "Type " << moduleTypes[i] << " has " << moduleNumberOfBeforeMap.at(moduleTypes[i]) << " modules before" << std::endl;
    std::cout << "Type " << moduleTypes[i] << " min " << moduleMinMap.at(moduleTypes[i]) << " mm " << std::endl;
    std::cout << "Type " << moduleTypes[i] << " max " << moduleMaxMap.at(moduleTypes[i]) << " mm " << std::endl;
    std::cout << "Type " << moduleTypes[i] << " separation at " << moduleSeparationMap.at(moduleTypes[i]) << " mm " << std::endl;
    std::cout << "Type " << moduleTypes[i] << " cells front " << mapOfFrontCellPerModuleType.at(moduleTypes[i]) <<  std::endl;
    std::cout << "Type " << moduleTypes[i] << " cells back " << mapOfBackCellPerModuleType.at(moduleTypes[i]) <<  std::endl;
  }

  // so now for each shower event, we can read showerModuleType, and get moduleSections, moduleNumber, moduleNumberOfBefore, moduleMin, moduleMax, moduleSeparation
  _file0->Close();
  //------------------------------------------------//


  outputFile->cd();
  // ---------------------------//
  // ACCUMULATORS               //
  // ---------------------------//
  // for each event, need to accumulate shower per module, section, cell
  // but modules can have different number of cells
  float *moduleEnDepInCrystal    = new float[modules_regions.size()];
  float *moduleEnDepNotInCrystal = new float[modules_regions.size()];
  float *moduleEnDepTotal        = new float[modules_regions.size()];

  float ***cellEnDepInCrystal    = new float**[modules_regions.size()];
  // Initialize

  // histo bin, min, man
  int bins  = 100;
  float min = 0.;
  float max = 1000.;

  module_t* module = new module_t[modules_regions.size()];
  for(int iMod = 0 ; iMod < modules_regions.size(); iMod++)
  {
    // std::cout << iMod << " ";
    // inizialize not to save
    module[iMod].save = false;

    // Initialize counters
    moduleEnDepInCrystal   [iMod] = 0.;
    moduleEnDepNotInCrystal[iMod] = 0.;
    moduleEnDepTotal       [iMod] = 0.;

    // Initialize histograms
    std::stringstream sname;
    sname << "Module_" << iMod << "_EnergyInCrystals";
    module[iMod].h_EnDepInCrystal    = new TH1F(sname.str().c_str(),sname.str().c_str(),bins,min,max);
    module[iMod].h_EnDepInCrystal->GetXaxis()->SetTitle("Energy [MeV]");
    module[iMod].h_EnDepInCrystal->GetYaxis()->SetTitle("Counts");
    sname.str("");

    sname << "Module_" << iMod << "_EnergyNotInCrystals";
    module[iMod].h_EnDepNotInCrystal = new TH1F(sname.str().c_str(),sname.str().c_str(),bins,min,max);
    module[iMod].h_EnDepNotInCrystal->GetXaxis()->SetTitle("Energy [MeV]");
    module[iMod].h_EnDepNotInCrystal->GetYaxis()->SetTitle("Counts");
    sname.str("");

    sname << "Module_" << iMod << "_EnergyTotal";
    module[iMod].h_EnDepTotal        = new TH1F(sname.str().c_str(),sname.str().c_str(),bins,min,max);
    module[iMod].h_EnDepTotal->GetXaxis()->SetTitle("Energy [MeV]");
    module[iMod].h_EnDepTotal->GetYaxis()->SetTitle("Counts");
    sname.str("");
    // Initialize cells
    int type = modules_regions[iMod].type;
    int sections = moduleSectionsMap.at(type);
    module[iMod].type = type;
    module[iMod].sections = sections;
    // module[iMod].cell = new cell_t*[sections];
    module[iMod].saveCell = new bool*[sections];
    module[iMod].h_Cell_EnDepInCrystal = new TH1F**[sections];
    module[iMod].nCellsPerSections = new int[sections];
    cellEnDepInCrystal[iMod] = new float*[sections];

    for(int iSec = 0; iSec < sections; iSec++)
    {
      int nCells;
      if(iSec == 0)
      {
        nCells = mapOfFrontCellPerModuleType.at(type);
      }
      else
      {
        nCells = mapOfBackCellPerModuleType.at(type);
      }
      module[iMod].nCellsPerSections[iSec] = nCells;
      cellEnDepInCrystal[iMod][iSec] = new float[nCells];

      // module[iMod].cell[iSec] = new cell_t[nCells];
      module[iMod].saveCell[iSec] = new bool[nCells];
      module[iMod].h_Cell_EnDepInCrystal[iSec] = new TH1F*[nCells];

      for(int iCel = 0; iCel < nCells; iCel++)
      {
        module[iMod].saveCell[iSec][iCel] = false;


        cellEnDepInCrystal[iMod][iSec][iCel] = 0; // <--------------
        std::stringstream bname;

        bname << "Cell_" << iCel ;
        if(iSec == 0)
        {
          bname << "_front_";
        }
        else
        {
          bname << "_back_";
        }

        sname << bname.str() << "Module_" << iMod << "_EnergyInCrystals";
        module[iMod].h_Cell_EnDepInCrystal[iSec][iCel] = new TH1F(sname.str().c_str(),sname.str().c_str(),bins,min,max);
        module[iMod].h_Cell_EnDepInCrystal[iSec][iCel]->GetXaxis()->SetTitle("Energy [MeV]");
        module[iMod].h_Cell_EnDepInCrystal[iSec][iCel]->GetYaxis()->SetTitle("Counts");
        sname.str("");

        // sname << bname.str() << "Module_" << iMod << "_EnergyNotInCrystals";
        // module[iMod].cell[iSec][iCel].h_EnDepNotInCrystal = new TH1F(sname.str().c_str(),sname.str().c_str(),bins,min,max);
        // module[iMod].cell[iSec][iCel].h_EnDepNotInCrystal->GetXaxis()->SetTitle("Energy [MeV]");
        // module[iMod].cell[iSec][iCel].h_EnDepNotInCrystal->GetYaxis()->SetTitle("Counts");
        // sname.str("");
        //
        // sname << bname.str() << "Module_" << iMod << "_EnergyTotal";
        // module[iMod].cell[iSec][iCel].h_EnDepTotal = new TH1F(sname.str().c_str(),sname.str().c_str(),bins,min,max);
        // module[iMod].cell[iSec][iCel].h_EnDepTotal->GetXaxis()->SetTitle("Energy [MeV]");
        // module[iMod].cell[iSec][iCel].h_EnDepTotal->GetYaxis()->SetTitle("Counts");
        // sname.str("");
      }
    }
  }





  //----------------------------------//
  // LOOP ON INPUT                    //
  //----------------------------------//


  TFile *inputFile = ((TFile *)0);
  TTree *shower = ((TTree*) 0);
  // now run on all input files
  for(unsigned int i = 0 ; i < listInputFiles.size(); i++)
  {
    if (inputFile) delete inputFile; // just a precaution
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
    inputFile->cd();
    shower = (TTree*) inputFile->Get("shower");
    if (!shower) continue; // requested ROOT file does not exist or is unreadable



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
    int feedbackDivision;
    long int nEntriesShower = shower->GetEntries();
    std::cout << "Energy deposition events = " << nEntriesShower <<  std::endl;
    std::cout << "Running on energy depositions..." << std::endl;
    long int counter = 0;
    int lostScintillation = 0;
    int lostShower = 0;
    feedbackDivision = (int) (nEntriesShower/20);
    if(feedbackDivision == 0) feedbackDivision = 1;
    for(int iEntry = 0 ; iEntry < nEntriesShower ; iEntry++)
    {
      if (  iEntry % feedbackDivision == 0 ) { std::cout << (int) ((double) iEntry/nEntriesShower * 100) << "% done..." << std::endl;}
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
      // std::cout << eventID << " " << showerEvent;
      if((eventID != showerEvent) && (eventID != -1)) // new event, dump the previous
      {
        // std::cout << " in";
        for(int iMod = 0 ; iMod < modules_regions.size(); iMod++)
        {
          if(moduleEnDepTotal[iMod] != 0)
          {
            module[iMod].save = true;
            module[iMod].h_EnDepInCrystal   ->Fill(moduleEnDepInCrystal   [iMod]);
            module[iMod].h_EnDepNotInCrystal->Fill(moduleEnDepNotInCrystal[iMod]);
            module[iMod].h_EnDepTotal       ->Fill(moduleEnDepTotal       [iMod]);
            for(int iSec = 0; iSec < module[iMod].sections; iSec++)
            {
              for(int iCel = 0; iCel <  module[iMod].nCellsPerSections[iSec]; iCel++)
              {
                if(cellEnDepInCrystal[iMod][iSec][iCel] != 0)
                {
                  module[iMod].saveCell[iSec][iCel] = true;
                  // outputFile->cd();
                  // module[iMod].cell[iSec][iCel].h_EnDepInCrystal->Fill(cellEnDepInCrystal[iMod][iSec][iCel]);
                  module[iMod].h_Cell_EnDepInCrystal[iSec][iCel]->Fill(cellEnDepInCrystal[iMod][iSec][iCel]);
                }

              }
            }
          }


        }
        // cleanup
        for(int iMod = 0 ; iMod < modules_regions.size(); iMod++)
        {
          moduleEnDepInCrystal   [iMod] = 0.;
          moduleEnDepNotInCrystal[iMod] = 0.;
          moduleEnDepTotal       [iMod] = 0.;
          for(int iSec = 0; iSec < module[iMod].sections; iSec++)
          {
            for(int iCel = 0; iCel <  module[iMod].nCellsPerSections[iSec]; iCel++)
            {
              cellEnDepInCrystal[iMod][iSec][iCel] = 0.;
            }
          }
        }
      }
      // std::cout << std::endl;

      // get module number
      Int_t binx = module_ID_map->GetXaxis()->FindBin(showerX);
      Int_t biny = module_ID_map->GetYaxis()->FindBin(showerY);
      Float_t xMod = module_ID_map->GetXaxis()->GetBinCenter(binx);
      Float_t yMod = module_ID_map->GetYaxis()->GetBinCenter(biny);
      int module_number = module_ID_map->GetBinContent(binx,biny);

      // accumulate on modules

      if(showerIsInCrystal)  moduleEnDepInCrystal   [module_number] += showerTotalEnDep;
      if(!showerIsInCrystal) moduleEnDepNotInCrystal[module_number] += showerTotalEnDep;
      moduleEnDepTotal[module_number] += showerTotalEnDep;

      //update event id
      eventID = showerEvent;

      // get cell number
      // cell will be skipped if module type not found in showerModuleType, so only en dep in crystals can be saved!

      if(showerIsInCrystal)
      {
        // first, cell map index
        int cellMapIndex;
        std::map<int, int> ::const_iterator iterCellMap = mapOfCellMaps.find(showerModuleType);
        if(iterCellMap != mapOfCellMaps.end())
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
          lostScintillation++;
          continue;
        }
        // then, get how many module sections are there
        int moduleSections;
        std::map<int, int> ::const_iterator iterModSec = moduleSectionsMap.find(showerModuleType);
        if(iterModSec != moduleSectionsMap.end())
        {
          // iter is item pair in the map. The value will be accessible as `iter->second`.
          moduleSections = (iterModSec->second);
        }
        else
        {
          if(verbose > 0)
          {
            std::cout << "WARNING: Module sections map not found for showerModuleType " << showerModuleType << std::endl;
            std::cout << "         This event will be ignored! Check the input files provided." << std::endl;
            std::cout << std::endl;
          }
          lostScintillation++;
          continue;
        }
        // get separation point
        float moduleSeparation;
        std::map<int, float> ::const_iterator iterModSep = moduleSeparationMap.find(showerModuleType);
        if(iterModSep != moduleSeparationMap.end())
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
          lostScintillation++;
          continue;
        }

        int fb = -1; // front_back, front = 0 , back = 1
        if(moduleSections == 1)
        {
          fb = 0; // only one section per module, so cell is only one
        }
        else // then it has to be == 2, given check at the beginning of program
        {
          if(showerZ > moduleSeparation) // "second" part of the module, so exit is for max z
          {
            fb = 1;
          }
          else
          {
            fb = 0; // front, i.e. z < 0
          }
        }

        Int_t c_binx    = cell_map_per_module[cellMapIndex][fb]->GetXaxis()->FindBin(showerX - xMod);
        Int_t c_biny    = cell_map_per_module[cellMapIndex][fb]->GetYaxis()->FindBin(showerY - yMod);
        int cell_number = cell_map_per_module[cellMapIndex][fb]->GetBinContent(c_binx,c_biny);


        cellEnDepInCrystal[module_number][fb][cell_number] += showerTotalEnDep;
      }


    }
    // fill last event!
    // ...
    // std::cout << " in";
    for(int iMod = 0 ; iMod < modules_regions.size(); iMod++)
    {
      if(moduleEnDepTotal[iMod] != 0)
      {
        module[iMod].save = true;
        module[iMod].h_EnDepInCrystal   ->Fill(moduleEnDepInCrystal   [iMod]);
        module[iMod].h_EnDepNotInCrystal->Fill(moduleEnDepNotInCrystal[iMod]);
        module[iMod].h_EnDepTotal       ->Fill(moduleEnDepTotal       [iMod]);
        for(int iSec = 0; iSec < module[iMod].sections; iSec++)
        {
          for(int iCel = 0; iCel <  module[iMod].nCellsPerSections[iSec]; iCel++)
          {
            if(cellEnDepInCrystal[iMod][iSec][iCel] != 0)
            {
              module[iMod].saveCell[iSec][iCel] = true;
              // outputFile->cd();
              // module[iMod].cell[iSec][iCel].h_EnDepInCrystal->Fill(cellEnDepInCrystal[iMod][iSec][iCel]);
              module[iMod].h_Cell_EnDepInCrystal[iSec][iCel]->Fill(cellEnDepInCrystal[iMod][iSec][iCel]);
            }
          }
        }
      }


    }
    // cleanup
    for(int iMod = 0 ; iMod < modules_regions.size(); iMod++)
    {
      moduleEnDepInCrystal   [iMod] = 0.;
      moduleEnDepNotInCrystal[iMod] = 0.;
      moduleEnDepTotal       [iMod] = 0.;
      for(int iSec = 0; iSec < module[iMod].sections; iSec++)
      {
        for(int iCel = 0; iCel <  module[iMod].nCellsPerSections[iSec]; iCel++)
        {
          cellEnDepInCrystal[iMod][iSec][iCel] = 0.;
        }
      }
    }



  }





  std::cout << "Writing output to file " << outputFileName << std::endl;
  outputFile->cd();
  module_ID_map->Write();
  // detector_ID_map_front->Write();
  // detector_ID_map_back->Write();
  for(int i = 0; i < NofModuleTypes; i++)
  {
    for(int j = 0; j < 2 ; j++)
    {
      cell_map_per_module[i][j]->Write();
    }
  }

  for(int iMod = 0 ; iMod < modules_regions.size(); iMod++)
  {
    if(module[iMod].save)
    {
      module[iMod].h_EnDepInCrystal   ->Write();
      module[iMod].h_EnDepNotInCrystal->Write();
      module[iMod].h_EnDepTotal       ->Write();
      for(int iSec = 0; iSec < module[iMod].sections; iSec++)
      {
        for(int iCel = 0; iCel <  module[iMod].nCellsPerSections[iSec]; iCel++)
        {
          if(module[iMod].saveCell[iSec][iCel])
          {
            module[iMod].h_Cell_EnDepInCrystal[iSec][iCel]->Write();
          }
        }
      }
    }

  }

  // manual delete
  for(int iMod = 0 ; iMod < modules_regions.size(); iMod++)
  {


      for(int iSec = 0; iSec < module[iMod].sections; iSec++)
      {
        for(int iCel = 0; iCel <  module[iMod].nCellsPerSections[iSec]; iCel++)
        {
          if(module[iMod].saveCell[iSec][iCel])
          {
            delete module[iMod].h_Cell_EnDepInCrystal[iSec][iCel];
          }
        }
      }

  }
  std::cout << "Done." << std::endl;
  // gROOT->ProcessLine(".qqqqqqqqqqqqqqqqqqqq");
  outputFile->Close();
  return 0;
}

  //
  //
  //
  // //------------------------------------------------//
  // // CALCULATE EN RES.                              //
  // //------------------------------------------------//
  // // prepare vars to sum energies
  // // one variable for each abs, cell and fiber , so dinamic arrays
  // float *absorberEnergyInCrystal;
  // float *absorberEnergyNotInCrystal;
  // float *cellEnergyInCrystal;
  // float *cellEnergyNotInCristal;
  // float *fibreEnergy;
  // // prepare debug variables
  // float *debug_absorberEnergyInCrystal;
  // float *debug_absorberEnergyNotInCrystal;
  // float *debug_cellEnergyInCrystal;
  // float *debug_cellEnergyNotInCristal;
  // float *debug_fibreEnergy;
  // // prepare histograms to plot energies
  // // one histo for each abs, cell and fiber, std::vectors of histos
  // std::stringstream sname;
  // std::stringstream stitle;
  // std::vector<TH1F*> abs_histoIn;
  // std::vector<TH1F*> abs_histoTot;
  // std::vector<TH1F*> cell_histoIn;
  // std::vector<TH1F*> cell_histoTot;
  // std::vector<TH1F*> fiber_histoIn;
  //
  //
  // std::vector<TH1F*> histo_pair_cell_Tot;
  // std::vector<TH1F*> histo_pair_cell_In;
  // std::vector<TH1F*> histo_ratio_pair_cell_Tot;
  // std::vector<TH1F*> histo_ratio_pair_cell_In;
  // float *pair_cellEnergyInCrystal;
  // float *pair_cellEnergyNotInCrystal;
  // float *ratio_pair_cellEnergyInCrystal;
  // float *ratio_pair_cellEnergyNotInCrystal;
  //
  // //hardcoded
  //
  // float energyInFrontInCry = 0.;
  // float energyInBackInCry  = 0.;
  // float energyInFrontNotInCry = 0.;
  // float energyInBackNotInCry  = 0.;
  // TH1F *h_energyInFrontInCry = new TH1F("EnergyInFrontInCry","EnergyInFrontInCry",bins,min,max);
  // h_energyInFrontInCry->GetXaxis()->SetTitle("MeV");
  // TH1F *h_energyInBackInCry  = new TH1F("EnergyInBackInCry" ,"EnergyInBackInCry" ,bins,min,max);
  // h_energyInBackInCry->GetXaxis()->SetTitle("MeV");
  // TH1F *h_energyInFrontInAbs = new TH1F("EnergyInFrontInAbs","EnergyInFrontInAbs",bins,min,max);
  // h_energyInFrontInAbs->GetXaxis()->SetTitle("MeV");
  // TH1F *h_energyInBackInAbs  = new TH1F("EnergyInBackInAbs" ,"EnergyInBackInAbs" ,bins,min,max);
  // h_energyInBackInAbs->GetXaxis()->SetTitle("MeV");
  // pair_cellEnergyInCrystal = new float[cells_regions.size()/2];
  // pair_cellEnergyNotInCrystal = new float[cells_regions.size()/2];
  // ratio_pair_cellEnergyInCrystal = new float[cells_regions.size()/2];
  // ratio_pair_cellEnergyNotInCrystal = new float[cells_regions.size()/2];
  // for(int i = 0 ; i < cells_regions.size()/2; i++)
  // {
  //   pair_cellEnergyInCrystal[i]          = 0.0;
  //   pair_cellEnergyNotInCrystal[i]       = 0.0;
  //   ratio_pair_cellEnergyInCrystal[i]    = 0.0;
  //   ratio_pair_cellEnergyNotInCrystal[i] = 0.0;
  // }
  // for(int i = 0 ; i < cells_regions.size()/2; i++)
  // {
  //   sname << "EnInCrystals_Cells_" << i << "_" << i+cells_regions.size()/2;
  //   stitle << "Energy deposited in crystals of Cells " << i << " and " << i+cells_regions.size()/2;
  //   TH1F *t_histoIn = new TH1F(sname.str().c_str(),stitle.str().c_str(),bins,min,max);
  //   t_histoIn->GetXaxis()->SetTitle("MeV");
  //   histo_pair_cell_In.push_back(t_histoIn);
  //   sname.str("");
  //   stitle.str("");
  // }
  // for(int i = 0 ; i < cells_regions.size()/2; i++)
  // {
  //   sname << "TotEn_Cells_" << i << "_" << i+cells_regions.size()/2;
  //   stitle << "Total energy deposited in Cells " << i << " and " << i+cells_regions.size()/2 << " (crystals + absorber + air)";
  //   TH1F *t_histoIn = new TH1F(sname.str().c_str(),stitle.str().c_str(),bins,min,max);
  //   t_histoIn->GetXaxis()->SetTitle("MeV");
  //   histo_pair_cell_Tot.push_back(t_histoIn);
  //   sname.str("");
  //   stitle.str("");
  // }
  // for(int i = 0 ; i < cells_regions.size()/2; i++)
  // {
  //   sname << "RatioBackFront_inCrystals_Cells" << i << "_" << i+cells_regions.size()/2;;
  //   stitle << "Ratio Back/Front, energy deposited in crystals of Cells " << i << " and " << i+cells_regions.size()/2;
  //   TH1F *t_histoIn = new TH1F(sname.str().c_str(),stitle.str().c_str(),bins,0,100);
  //   t_histoIn->GetXaxis()->SetTitle("Ratio back/front");
  //   histo_ratio_pair_cell_In.push_back(t_histoIn);
  //   sname.str("");
  //   stitle.str("");
  // }
  // for(int i = 0 ; i < cells_regions.size()/2; i++)
  // {
  //   sname << "RatioBackFront_Tot_Cells" << i << "_" << i+cells_regions.size()/2;;
  //   stitle << "Ratio Back/Front, total energy deposited in Cells " << i << " and " << i+cells_regions.size()/2 << " (crystals + absorber + air)";
  //   TH1F *t_histoIn = new TH1F(sname.str().c_str(),stitle.str().c_str(),bins,0,100);
  //   t_histoIn->GetXaxis()->SetTitle("Ratio back/front");
  //   histo_ratio_pair_cell_Tot.push_back(t_histoIn);
  //   sname.str("");
  //   stitle.str("");
  // }
  //
  //
  // if(sum_abs)
  // {
  //   absorberEnergyInCrystal = new float[absobers_regions.size()];
  //   absorberEnergyNotInCrystal = new float[absobers_regions.size()];
  //   if(debug) // some debug
  //   {
  //     debug_absorberEnergyInCrystal = new float[absobers_regions.size()];
  //     debug_absorberEnergyNotInCrystal = new float[absobers_regions.size()];
  //   }
  //   for(int i = 0 ; i < absobers_regions.size(); i++)
  //   {
  //     absorberEnergyInCrystal[i] = 0;
  //     absorberEnergyNotInCrystal[i] = 0;
  //     if(debug) // some debug
  //     {
  //       debug_absorberEnergyInCrystal[i] = 0;
  //       debug_absorberEnergyNotInCrystal[i] = 0;
  //     }
  //   }
  //   for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //   {
  //     sname << "EnInCrystals_Absorber" << iAbs;
  //     stitle << "Energy deposited in crystals of Absorber " << iAbs;
  //     TH1F *t_histoIn = new TH1F(sname.str().c_str(),stitle.str().c_str(),bins,min,max);
  //     t_histoIn->GetXaxis()->SetTitle("MeV");
  //     abs_histoIn.push_back(t_histoIn);
  //     sname.str("");
  //     stitle.str("");
  //   }
  //   for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //   {
  //     sname << "EnTot_Absorber" << iAbs;
  //     stitle << "Total Energy in Absorber " << iAbs;
  //     TH1F *t_histoIn = new TH1F(sname.str().c_str(),stitle.str().c_str(),bins,min,max);
  //     t_histoIn->GetXaxis()->SetTitle("MeV");
  //     abs_histoTot.push_back(t_histoIn);
  //     sname.str("");
  //     stitle.str("");
  //   }
  // }
  // if(sum_cell)
  // {
  //   cellEnergyInCrystal     = new float[cells_regions.size()];
  //   cellEnergyNotInCristal     = new float[cells_regions.size()];
  //   if(debug) // some debug
  //   {
  //     debug_cellEnergyInCrystal     = new float[cells_regions.size()];
  //     debug_cellEnergyNotInCristal     = new float[cells_regions.size()];
  //   }
  //   for(int i = 0 ; i < cells_regions.size(); i++)
  //   {
  //     cellEnergyInCrystal[i] = 0;
  //     cellEnergyNotInCristal[i] = 0;
  //     if(debug) // some debug
  //     {
  //       debug_cellEnergyInCrystal[i] = 0;
  //       debug_cellEnergyNotInCristal[i] = 0;
  //     }
  //   }
  //   for(int i = 0 ; i < cells_regions.size(); i++)
  //   {
  //     sname << "EnInCrystals_Cell" << i;
  //     stitle << "Energy deposited in crystals of Cell " << i;
  //     TH1F *t_histoIn = new TH1F(sname.str().c_str(),stitle.str().c_str(),bins,min,max);
  //     t_histoIn->GetXaxis()->SetTitle("MeV");
  //     cell_histoIn.push_back(t_histoIn);
  //     sname.str("");
  //     stitle.str("");
  //   }
  //   for(int i = 0 ; i < cells_regions.size(); i++)
  //   {
  //     sname << "EnTot_Cell" << i;
  //     stitle << "Total energy deposited in Cell " << i << " (crystals + absorber + air)";
  //     TH1F *t_histoTot = new TH1F(sname.str().c_str(),stitle.str().c_str(),bins,min,max);
  //     t_histoTot->GetXaxis()->SetTitle("MeV");
  //     cell_histoTot.push_back(t_histoTot);
  //     sname.str("");
  //     stitle.str("");
  //   }
  // }
  // if(sum_cry)
  // {
  //   fibreEnergy    = new float[fibres_regions.size()];
  //   if(debug)  // some debug
  //   {
  //     debug_fibreEnergy    = new float[fibres_regions.size()];
  //   }
  //   for(int i = 0 ; i < fibres_regions.size(); i++)
  //   {
  //     fibreEnergy[i] = 0;
  //     if(debug)  // some debug
  //     {
  //       debug_fibreEnergy = 0;
  //     }
  //   }
  //   for(int i = 0 ; i < fibres_regions.size(); i++)
  //   {
  //     sname << "EnIn_Crystal" << i;
  //     stitle << "Total energy deposited in Crystal " << i ;
  //     TH1F *t_histoTot = new TH1F(sname.str().c_str(),stitle.str().c_str(),bins,min,max);
  //     t_histoTot->GetXaxis()->SetTitle("MeV");
  //     fiber_histoIn.push_back(t_histoTot);
  //     sname.str("");
  //     stitle.str("");
  //   }
  // }
  //
  // // loop on events
  // long int nEvents = shower->GetEntries();
  // int eventID = -1;
  // long long int counter = 0;
  // float dumb_sum = 0.;
  // for(int i = 0 ; i < nEvents ; i++) // loop on ttree entries
  // {
  //   bool inAbs = false;
  //   bool inCry = false;
  //   bool inDumb = false;
  //   // get the entry
  //   shower->GetEvent(i);
  //
  //   // now summation of energy needs to be performed only until a new event arrives in the list
  //   // an event corresponds to 1 primary shot. the info is stored in event
  //   // so we start be setting a event indicator, eventID, to -1. At the first loop step,
  //   // the first value of variable event will be 0, so what comes in the next if statement
  //   // (event saving and clean up) will not be executed, and instead the event accumulation starts
  //   // and the end of event accumulation the eventID becomes = event, then a new step of the cycle
  //   // starts. Now eventID is = 0, so "event saving and clean up" will not be performed until
  //   // the value of event becomes != 0. When that happens, the variables were the enegies were
  //   // summed are used to filled the relevant histograms, than they are cleared and a new accumulation
  //   // starts. Of course this means that after the end of the loop, the last event still needs to be saved
  //
  //
  //   //------------------------------//
  //   // event saving and clean up    //
  //   //------------------------------//
  //   if((eventID != event) && (eventID != -1)) // new event, dump the previous
  //   {
  //     //fill histograms
  //     if(sum_abs)
  //     {
  //       for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //       {
  //         abs_histoIn[iAbs]->Fill(absorberEnergyInCrystal[iAbs]);
  //         abs_histoTot[iAbs]->Fill(absorberEnergyInCrystal[iAbs] + absorberEnergyNotInCrystal[iAbs]);
  //       }
  //     }
  //     if(sum_cell)
  //     {
  //       h_energyInFrontInCry->Fill(energyInFrontInCry);
  //       h_energyInFrontInAbs->Fill(energyInFrontInCry + energyInFrontNotInCry);
  //       h_energyInBackInCry->Fill(energyInBackInCry);
  //       h_energyInBackInAbs->Fill(energyInBackInCry + energyInBackNotInCry);
  //       for(int iCell = 0 ; iCell < cells_regions.size(); iCell++)
  //       {
  //         cell_histoIn[iCell]->Fill(cellEnergyInCrystal[iCell]);
  //         cell_histoTot[iCell]->Fill(cellEnergyInCrystal[iCell] + cellEnergyNotInCristal[iCell]);
  //       }
  //       for(int iCell = 0 ; iCell < cells_regions.size()/2; iCell++)
  //       {
  //         histo_pair_cell_In[iCell]->Fill(pair_cellEnergyInCrystal[iCell]);
  //         histo_pair_cell_Tot[iCell]->Fill(pair_cellEnergyInCrystal[iCell] + pair_cellEnergyNotInCrystal[iCell]);
  //         histo_ratio_pair_cell_In[iCell]->Fill(cellEnergyInCrystal[iCell+9]/cellEnergyInCrystal[iCell]);
  //         histo_ratio_pair_cell_Tot[iCell]->Fill( (cellEnergyInCrystal[iCell+9] + cellEnergyNotInCristal[iCell+9])/(cellEnergyInCrystal[iCell] + cellEnergyNotInCristal[iCell]) );
  //       }
  //
  //     }
  //     if(sum_cry)
  //     {
  //       for(int iCry = 0 ; iCry < fibres_regions.size(); iCry++)
  //       {
  //         fiber_histoIn[iCry]->Fill(fibreEnergy[iCry]);
  //       }
  //     }
  //
  //     if(debug)// DEBUG output sums
  //     {
  //       std::cout << std::endl;
  //       std::cout << "EVENT " << eventID << std::endl;
  //       if(sum_abs)
  //       {
  //         for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //         {
  //           std::cout << "Absorber " << iAbs << std::endl;
  //           std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
  //                     << absorberEnergyInCrystal[iAbs] << std::endl;
  //           std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
  //                     << absorberEnergyNotInCrystal[iAbs] << std::endl;
  //         }
  //       }
  //       if(sum_cell)
  //       {
  //
  //         for(int iCell = 0 ; iCell < cells_regions.size(); iCell++)
  //         {
  //           std::cout << "Cell " << iCell << std::endl;
  //           std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
  //                     << cellEnergyInCrystal[iCell] << std::endl;
  //           std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1)
  //                     << cellEnergyNotInCristal[iCell] << std::endl;
  //
  //         }
  //
  //       }
  //       if(sum_cry)
  //       {
  //         float sum_in_cry = 0;
  //         for(int iCry = 0 ; iCry < fibres_regions.size(); iCry++)
  //         {
  //           // if(fibreEnergy[iCry] != 0)
  //           // {
  //           //   std::cout << fibreEnergy[iCry] << std::endl;
  //           // }
  //           sum_in_cry += fibreEnergy[iCry];
  //         }
  //         std::cout << "SUM IN CRYSTALS = " << sum_in_cry << std::endl;
  //       }
  //       std::cout << "dumb sum = " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << dumb_sum << std::endl;
  //       std::cout << std::endl;
  //     } // end of debug
  //
  //
  //
  //     //clean variables
  //     dumb_sum = 0;
  //     energyInFrontInCry = 0.0;
  //     energyInFrontNotInCry = 0.0;
  //     energyInBackInCry = 0.0;
  //     energyInBackNotInCry = 0.0;
  //
  //     if(sum_abs)
  //     {
  //       for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //       {
  //         absorberEnergyInCrystal[iAbs] = 0;
  //         absorberEnergyNotInCrystal[iAbs] = 0;
  //       }
  //     }
  //     if(sum_cell)
  //     {
  //       for(int iCell = 0 ; iCell < cells_regions.size(); iCell++)
  //       {
  //         cellEnergyInCrystal[iCell] = 0;
  //         cellEnergyNotInCristal[iCell] = 0;
  //       }
  //       for(int iCell = 0 ; iCell < cells_regions.size()/2; iCell++)
  //       {
  //         pair_cellEnergyInCrystal[iCell] = 0;
  //         pair_cellEnergyNotInCrystal[iCell] = 0;
  //       }
  //     }
  //     if(sum_cry)
  //     {
  //       for(int iCry = 0 ; iCry < fibres_regions.size(); iCry++)
  //       {
  //         fibreEnergy[iCry] = 0;
  //       }
  //     }
  //   }
  //   // end of event saving and clean up //
  //   //------------------------------//
  //
  //
  //   //------------------------------//
  //   // event accumulation           //
  //   //------------------------------//
  //
  //   if(crystalID != -1)
  //   {
  //     inDumb = true;
  //     dumb_sum += showerIonizingEnDep;
  //   }
  //
  //   if(sum_abs)
  //   {
  //     //check if event inside absorber
  //     for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //     {
  //       if( (showerX >= absobers_regions[iAbs].min[0]) &&
  //           (showerX <  absobers_regions[iAbs].max[0]) &&
  //           (showerY >= absobers_regions[iAbs].min[1]) &&
  //           (showerY <  absobers_regions[iAbs].max[1]) &&
  //           (showerZ >= absobers_regions[iAbs].min[2]) &&
  //           (showerZ <  absobers_regions[iAbs].max[2])   ) // inside the abs region
  //       {
  //         if(showerIsInCrystal) /// if in crystal
  //         {
  //           // for(int iCell = 0 ; iCell < cells_regions.size(); iCell++) // find cell
  //           // {
  //           //   if( (showerX >= cells_regions[iCell].min[0]) &&
  //           //       (showerX <  cells_regions[iCell].max[0]) &&
  //           //       (showerY >= cells_regions[iCell].min[1]) &&
  //           //       (showerY <  cells_regions[iCell].max[1]) &&
  //           //       (showerZ >= cells_regions[iCell].min[2]) &&
  //           //       (showerZ <  cells_regions[iCell].max[2])   ) // inside this cell region
  //           //   {
  //           //
  //           //     // if(showerIsInCrystal)
  //           //     // {
  //           //     //   cellEnergyInCrystal[iCell] += showerIonizingEnDep;
  //           //     // }
  //           //     // else
  //           //     // {
  //           //     //   cellEnergyNotInCristal[iCell] += showerIonizingEnDep;
  //           //     // }
  //           //   }
  //           // }
  //
  //
  //
  //           absorberEnergyInCrystal[iAbs] += showerIonizingEnDep;
  //           inAbs = true;
  //         }
  //         else
  //         {
  //           absorberEnergyNotInCrystal[iAbs] += showerIonizingEnDep;
  //         }
  //       }
  //     }
  //   }
  //   if(sum_cell)
  //   {
  //     //check if event inside cell
  //     for(int iCell = 0 ; iCell < cells_regions.size(); iCell++)
  //     {
  //       if( (showerX >= cells_regions[iCell].min[0]) &&
  //           (showerX <  cells_regions[iCell].max[0]) &&
  //           (showerY >= cells_regions[iCell].min[1]) &&
  //           (showerY <  cells_regions[iCell].max[1]) &&
  //           (showerZ >= cells_regions[iCell].min[2]) &&
  //           (showerZ <  cells_regions[iCell].max[2])   ) // inside this cell region
  //       {
  //
  //         if(showerIsInCrystal)
  //         {
  //           cellEnergyInCrystal[iCell] += showerIonizingEnDep;
  //           if(iCell > 8)
  //           {
  //             pair_cellEnergyInCrystal[iCell - 9] += showerIonizingEnDep;
  //             energyInBackInCry  += showerIonizingEnDep;
  //           }
  //           else
  //           {
  //             pair_cellEnergyInCrystal[iCell] += showerIonizingEnDep;
  //             energyInFrontInCry  += showerIonizingEnDep;
  //           }
  //         }
  //         else
  //         {
  //           cellEnergyNotInCristal[iCell] += showerIonizingEnDep;
  //           if(iCell > 8)
  //           {
  //             pair_cellEnergyNotInCrystal[iCell - 9] += showerIonizingEnDep;
  //             energyInBackNotInCry  += showerIonizingEnDep;
  //           }
  //           else
  //           {
  //             pair_cellEnergyNotInCrystal[iCell] += showerIonizingEnDep;
  //             energyInFrontNotInCry  += showerIonizingEnDep;
  //           }
  //         }
  //       }
  //     }
  //   }
  //   if(sum_cry)
  //   {
  //     //check if event inside fiber
  //     for(int iCry = 0 ; iCry < fibres_regions.size(); iCry++)
  //     {
  //       if( (showerX >= fibres_regions[iCry].min[0]) &&
  //           (showerX <= fibres_regions[iCry].max[0]) &&
  //           (showerY >= fibres_regions[iCry].min[1]) &&
  //           (showerY <= fibres_regions[iCry].max[1]) &&
  //           (showerZ >= fibres_regions[iCry].min[2]) &&
  //           (showerZ <= fibres_regions[iCry].max[2])   ) // inside this fiber region
  //       {
  //         // no need to check if it is in a crystal...
  //         inCry = true;
  //         fibreEnergy[iCry] += showerIonizingEnDep;
  //       }
  //
  //     }
  //   }
  //   // end of event accumulation
  //   //------------------------------//
  //
  //   //DEBUG
  //   if(debug)  // some debug
  //   {
  //     if(inAbs == true)
  //     {
  //       if(inCry == false)
  //       {
  //         std::cout << "found " << event << " " << i << std::endl;
  //         std::cout << showerIsInCrystal << " "
  //                   << std::setprecision(std::numeric_limits<double>::digits10 + 1) << showerX << " "
  //                   << std::setprecision(std::numeric_limits<double>::digits10 + 1) << showerY << " "
  //                   << std::setprecision(std::numeric_limits<double>::digits10 + 1) << showerZ << " "
  //                   << std::endl;
  //       }
  //     }
  //
  //     if(inAbs == false)
  //     {
  //       if(inCry == true)
  //       {
  //         std::cout << "weird " << event << " " << i << std::endl;
  //         std::cout << showerIsInCrystal << " "
  //                   << std::setprecision(std::numeric_limits<double>::digits10 + 1) << showerX << " "
  //                   << std::setprecision(std::numeric_limits<double>::digits10 + 1) << showerY << " "
  //                   << std::setprecision(std::numeric_limits<double>::digits10 + 1) << showerZ << " "
  //                   << std::endl;
  //       }
  //     }
  //   }
  //
  //
  //
  //   // update event id
  //   eventID = event;
  //   // give feedback to user
  //   counter++;
  //   int perc = ((100*counter)/nEvents);
  //   if(!debug)
  //   {
  //     if( (perc % 10) == 0 )
  //     {
  //       std::cout << "\r";
  //       std::cout << perc << "% done... ";
  //     }
  //   }
  // } // end of loop on ttree entries
  //
  // //fill last event
  // if(sum_abs)
  // {
  //   for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //   {
  //     abs_histoIn[iAbs]->Fill(absorberEnergyInCrystal[iAbs]);
  //     abs_histoTot[iAbs]->Fill(absorberEnergyInCrystal[iAbs] + absorberEnergyNotInCrystal[iAbs]);
  //   }
  // }
  // if(sum_cell)
  // {
  //   for(int iCell = 0 ; iCell < cells_regions.size(); iCell++)
  //   {
  //     cell_histoIn[iCell]->Fill(cellEnergyInCrystal[iCell]);
  //     cell_histoTot[iCell]->Fill(cellEnergyInCrystal[iCell] + cellEnergyNotInCristal[iCell]);
  //   }
  // }
  // if(sum_cry)
  // {
  //   for(int iCry = 0 ; iCry < fibres_regions.size(); iCry++)
  //   {
  //     fiber_histoIn[iCry]->Fill(fibreEnergy[iCry]);
  //   }
  // }
  //
  // if(debug)// DEBUG output sums
  // {
  //   std::cout << std::endl;
  //   std::cout << "EVENT " << eventID << std::endl;
  //   if(sum_abs)
  //   {
  //     for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //     {
  //       std::cout << "Absorber " << iAbs << std::endl;
  //       std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << absorberEnergyInCrystal[iAbs] << std::endl;
  //       std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << absorberEnergyNotInCrystal[iAbs] << std::endl;
  //     }
  //   }
  //   if(sum_cell)
  //   {
  //     for(int iCell = 0 ; iCell < cells_regions.size(); iCell++)
  //     {
  //       std::cout << "Cell " << iCell << std::endl;
  //       std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << cellEnergyInCrystal[iCell] << std::endl;
  //       std::cout << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << cellEnergyNotInCristal[iCell] << std::endl;
  //     }
  //   }
  //   if(sum_cry)
  //   {
  //     float sum_in_cry = 0;
  //     for(int iCry = 0 ; iCry < fibres_regions.size(); iCry++)
  //     {
  //       // if(fibreEnergy[iCry] != 0)
  //       // {
  //       //   std::cout << fibreEnergy[iCry] << std::endl;
  //       // }
  //       sum_in_cry += fibreEnergy[iCry];
  //     }
  //     std::cout << "SUM IN CRYSTALS = " << sum_in_cry << std::endl;
  //   }
  //
  //   std::cout << "dumb sum = " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) <<  dumb_sum << std::endl;
  //   std::cout << std::endl;
  // } // end of debug
  // // end of CALCULATE EN RES
  // //------------------------------------------------//
  //
  //
  // //------------------------------------------------//
  // // SAVE TO OUTPUT                                 //
  // //------------------------------------------------//
  // std::string outputRootFileName = outputFileName + ".root";
  // std::string outputTextFileName = outputFileName + ".txt";
  // // save to text file
  // std::ofstream outputTextFile;
  // std::cout <<std::endl;
  // std::cout << "Saving results in "<< outputTextFileName.c_str() << std::endl ;
  // outputTextFile.open (outputTextFileName.c_str(),std::ofstream::out);
  // outputTextFile << "Abs Number\t FIT Fraction Tot En[MeV]\t FIT Energy Res. In Cry (sigma/E) \tFraction Tot En[MeV]\t Energy Res. In Cry (sigma/E)" << std::endl;
  // if(sum_abs)
  // {
  //   // fit total and crystal abs spectra
  //   for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //   {
  //     TF1 *gaussIn = new TF1("gaussIn","gaus",min,max);
  //     abs_histoIn[iAbs]->Fit(gaussIn,"Q");
  //
  //     TF1 *gaussTot = new TF1("gaussTot","gaus",min,max);
  //     abs_histoTot[iAbs]->Fit(gaussTot,"Q");
  //
  //     outputTextFile << iAbs << "\t"
  //                    << gaussTot->GetParameter(1) / max << "\t"
  //                    << gaussIn->GetParameter(2) / gaussIn->GetParameter(1) << "\t"
  //                    << abs_histoTot[iAbs]->GetMean() / max << "\t"
  //                    << abs_histoIn[iAbs]->GetRMS() / abs_histoIn[iAbs]->GetMean() << "\t"
  //                    << std::endl;
  //   }
  // }
  // outputTextFile.close();
  // std::cout <<std::endl;
  // // save to ROOT file
  // std::cout << "Saving histograms to "<< outputRootFileName.c_str() << " ..." ;
  // TFile *outputFile = new TFile(outputRootFileName.c_str(),"RECREATE");
  // outputFile->cd();
  // if(sum_abs)
  // {
  //   for(int iAbs = 0; iAbs < absobers_regions.size(); iAbs++)
  //   {
  //     abs_histoIn[iAbs]->Write();
  //     abs_histoTot[iAbs]->Write();
  //   }
  // }
  // if(sum_cell)
  // {
  //   h_energyInFrontInCry->Write();
  //   h_energyInFrontInAbs->Write();
  //   h_energyInBackInCry->Write();
  //   h_energyInBackInAbs->Write();
  //   for(int iCell = 0 ; iCell < cells_regions.size(); iCell++)
  //   {
  //     cell_histoIn[iCell]->Write();
  //     cell_histoTot[iCell]->Write();
  //   }
  //   for(int iCell = 0 ; iCell < cells_regions.size()/2; iCell++)
  //   {
  //     histo_pair_cell_In[iCell]->Write();
  //     histo_pair_cell_Tot[iCell]->Write();
  //     histo_ratio_pair_cell_In[iCell]->Write();
  //     histo_ratio_pair_cell_Tot[iCell]->Write();
  //   }
  //
  // }
  // if(sum_cry)
  // {
  //   for(int iCry = 0 ; iCry < fibres_regions.size(); iCry++)
  //   {
  //     fiber_histoIn[iCry]->Write();
  //   }
  // }
  // outputFile->Close();
  // // end of SAVE TO OUTPUT
  // //------------------------------------------------//
  //
  //
  // std::cout <<" done!" << std::endl;
