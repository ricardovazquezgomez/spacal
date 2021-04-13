// c++ -O2 -std=c++17 -o pruneSigFormOutput pruneSigFormOutput.cpp `root-config --cflags --glibs` 

#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "ROOT/RDataFrame.hxx"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <map>
#include <dirent.h>


// ### D E F I N E S ###
// #define TREE_N    "tree"
// #define TL_N      "Total_Light"
// #define TCh4_N    "t_pulse4"
// #define TCh13_N   "t_pulse13"
// #define LCh4_N    "ch4"
// #define LCh13_N   "ch13"
// #define CFD_N     "CFD_Thresh"

// #define FILE_SELECTION ( file.find(excludeS) == std::string::npos  && !file.compare(0,std::string(filePref).size(),std::string(filePref)) )


void PrintUsage () { std::cout << " Usage: " << std::endl; std::cout << " ./pruneSigFormOutput [-i DataInputFolder] [-k files' key word] [-b pruned branches' key word] [-f folders' key word] [-t input/output tree name] [-o output file name]" << std::endl; }
std::vector<std::string>  Prune (std::vector<std::string> columnNames, const std::string& keyWord)
{
  // ################################################################
  // ### Takes a vector and remove the entries starting with keyword.
  // ### rfind finds the last occurence of keyWord starting no later than 0,
  // ### therefore starting at 0. Returns the position, which must be 0, or npos.
  auto it = columnNames.begin();

  while (it != columnNames.end()) {
    if ( it->rfind(keyWord, 0) == 0) {          
      it = columnNames.erase(it);
    } 

    else ++it;    
  }

  return columnNames;
  
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




int main (int argc, char **argv)
{
  // --> Read the Command Line  -->
  std::string path = "-1", outFileName = "prunedOutput", fileKeyWord = "OutTrigd", folderKeyWord = "GeV", branchKeyWord = "pulse", treeName = "tree";
  if  ( argc < 2 ) {  
    PrintUsage();
    return -1;
  }

  else {
    for (int i = 0; i < argc; i++) {
      if      ( std::string(argv[i])  == "-i") path         = argv [i+1]; 
      else if ( std::string(argv[i])  == "-k") fileKeyWord  = argv [i+1];
      else if ( std::string(argv[i])  == "-b") branchKeyWord= argv [i+1];
      else if ( std::string(argv[i])  == "-f") folderKeyWord= argv [i+1];
      else if ( std::string(argv[i])  == "-t") treeName     = argv [i+1];
      else if ( std::string(argv[i])  == "-o") outFileName  = argv [i+1];
    }
  }

  if (path == "-1") {std::cout << "No input folder path provided, the program will be stopped.\n"; return -1; }

  std::cout << "Looking for files matching " << fileKeyWord << " into folders with " << folderKeyWord << " in " << path <<"\n";
  std::cout << "Pruning branches containing the keyword: " << branchKeyWord << std::endl;
  std::cout << "\n### Warning: energy values will be taken from the folders names.\n";




  // --- Add '/' at the end of the input path if needed
  if (path.back() != '/') path += '/';



  // --- Find the folders matching the keyword
  std::vector <std::string> dataFolders, folders;
  std::vector <float> energies;



  read_directory( path, folders);
  for (auto && elem : folders) {
    if (elem.find(folderKeyWord) != std::string::npos) {
      dataFolders.emplace_back(path + elem + "/");                                                                        // Save the path
      energies.emplace_back   (std::atof( elem.erase (elem.find(folderKeyWord), folderKeyWord.size()).c_str() ));         // Save only the energy
    }    
  }

  // --- Feedback
  std::cout << "Found " << dataFolders.size() << " folders.\n";
  std::cout << "Energies: ";
  for (auto && elem : energies) std::cout << elem << " ";
  std::cout << std::endl;


 













  // #################################
  // --> LOOP OVER ALL THE FOLDERS -->
  for (size_t iFol = 0; iFol < dataFolders.size(); ++iFol)
  {
    time_t _start = clock();

    std::string dataFolder = dataFolders.at(iFol);
    float trueEn = energies.at(iFol);

    // -- Extract the data files from the folders
    std::cout << "\n###\nOpening folder " << dataFolder << "\n";
    std::cout << "Energy: " << trueEn << "\n";
    std::vector<std::string> files, dataFiles;
    read_directory (dataFolder, files);

    for (auto && file : files)  {
      if (file.find(fileKeyWord) != std::string::npos) {
        dataFiles.emplace_back(dataFolder + file);      
      }    
    }

    if (dataFiles.size()) {
      std::cout << "Found: " << dataFiles.size() << " files.\n";
    }
    else {
      std::cout << "No files found. Skipping the folder.\n";
      continue;
    }


    ROOT::RDataFrame d (treeName.c_str(), dataFiles);
    // for (auto && dataFile : dataFiles) chain -> Add(dataFile.c_str());

    // --- Add the energy branch
    // chains[iFol] -> Branch("trueEn", &trueEn, "trueEn/F");
    // for (Long64_t i = 0; i < chains[iFol]->GetEntries(); ++i) chains[iFol] -> Fill(); 
    auto dEn = d.Define ("folderEnergy", [&trueEn](){return trueEn;} );


    dEn.Snapshot(treeName.c_str(), (dataFolder + outFileName + ".root").c_str(), Prune(dEn.GetColumnNames(), branchKeyWord));
    std::cout << "### Data of folder " << dataFolder << " grouped, pruned, and saved at " << (dataFolder + outFileName + ".root").c_str() << "\n";
    std::cout << "### Duration: " << (std::clock() - _start) / (double) CLOCKS_PER_SEC << "s\n";
  } // <-- For loop - END <--




  // --> Now Chain all the Chains together. Is it even possible?



  // outFile.Close();
  std::cout << "\n";
  std::cout << "Saved to: " << (path).c_str() << "\n";
  std::cout << "### Standard end of the program. ###" << std::endl;
}
