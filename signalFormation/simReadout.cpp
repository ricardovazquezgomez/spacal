#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <functional>


#include "ConfigFile.hh"
#include "ApplySpatialFilters.hh"
// #include "ApplyQEFilters.hh"
#include "AssignDetector.hh"



#include "LHamamatsuR7899.hh"
#include "LHamamatsuR7600U_20.hh"
#include "LHamamatsuH13543.hh"



#include "TROOT.h"
#include "ROOT/RDataFrame.hxx"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TRandom3.h"
#include "TMath.h"



#include <dirent.h>
#include "../parametrization/regions.h"

/* - ### - ### - ### - ### - ### - */
struct event_t
{
    // Each event
    // Module --> Detector --> Datum
    int nEvent      = -1;
    int nPhotonsTot =  0;
    std::vector <int> modulesHit;
    std::vector <std::vector<int>> nPhotonsDet;
    std::vector <std::vector<std::vector<float>>>  detectorTimeStamps;   
    std::vector <std::vector<std::vector<float>>>  digiPulse;
    std::vector <std::vector<std::vector<float>>>  digiPulseDepo;
    std::vector <std::vector<std::vector<float>>>  digiPulseScint;
    std::vector <std::vector<std::vector<float>>>  digiPulseDepoTrans;
};

// --| Objects
class data_t 
{
    std::map <int, event_t> _M;

    const int   _nModules;
    const int   _sampling;
    const float _gate;
    const float _delay;
    const float _baseline;
    const float _timeBin;
    const float _sigmaNoise;

    std::vector<int> _moduleTypes;
    std::vector<int> _nFrontDetPerType;
    std::vector<int> _nBackDetPerType;
    std::map<int, std::vector<int>> _moduleMap;

    std::string _nBWave;
    std::string _nBWaveDepo;
    std::string _nBWaveScint;
    std::string _nBWaveDepoTrans;
    std::string _nBPhot;
    std::string _nBMod;


 public:
    data_t (const ConfigFile & C, const int & nModules, std::map<int, std::vector<int>> & moduleMap);
    ~data_t(){};

    // --- Members
    void UpdateMap          (const int &nEv, const float &nTime, const int & nMod, const int & modType, const int & nDet, const LPhD & photoDet, TRandom3 *randGen);
    // void UpdateMapAllPulses (const int &nEv, const float &nTime, const float & t_dep, const float & t_gen, const float & t_prop, const int & nDet, const LPhD & photoDet, TRandom3 *randGen);
    void FillPulse          (std::vector<float> & Pulse, const float &nTime, const LPhD & photoDet, const float & randSPTR, const float & randAmpl);
    void PrintTheMap() { for (auto& iEv : _M) std::cout << "Key:" << iEv.first << "\tEvent: " << iEv.second.nEvent << "\tTotal number of photons: " << iEv.second.nPhotonsTot << std::endl;};
    void SaveToTTree          (TTree* TT);
    // void SaveToTTreeAllPulses (TTree* TT);
    void SaveTimeStamps       (TTree* TT, int nTS_Saved);

    // --- Get members
    inline const int   Get_nModules()  {return _nModules;};
    inline const int   Get_sampling()  {return _sampling;};
    inline const float Get_gate()      {return _gate;};
    inline const float Get_delay()     {return _delay;};
    inline const float Get_baseline()  {return _baseline;};
    inline const float Get_sigmaNoise(){return _sigmaNoise;};
    inline const int   Get_nCellFront(const int & moduleType) { return _moduleMap[moduleType].at(0); };
    inline const int   Get_nCellBack (const int & moduleType) { return _moduleMap[moduleType].at(1); };

};
/* - ### - ### - ### - ### - ### - */




// --| Functions
void PrintUsage () { std::cout << " Usage: " << std::endl; std::cout << " ./simReadout [-c ConfigFile] [-i DataInputFile] [-f FiltersFile] [-o OutputFileName] [-g PMT Gain]" << std::endl; }
bool Detection  (int event, float x, float y, float z, float timestamp, float PhotonEnergy) {return true;}
void GetGeometry (const std::string & inputFileName, int & nModules, std::map<int, std::vector<int>> & moduleMap,  std::vector<TH2I> & geoMaps);







int main (int argc, char** argv)
{

    // --> Read the Command Line
    std::string ConfigFileName = "ConfigFile.cfg", InputFile = "-1", OutputFile = "-1", FiltersFile="-1";
    float input_Gain = -1.;
    if  ( argc < 2 ) {  
        PrintUsage();
        return -1; }

    else {
        for (int i = 0; i < argc; i++)
        {
            if      ( std::string(argv[i])  == "-c") ConfigFileName = argv [i+1];
            else if ( std::string(argv[i])  == "-i") InputFile      = argv [i+1]; 
            else if ( std::string(argv[i])  == "-o") OutputFile     = argv [i+1]; 
            else if ( std::string(argv[i])  == "-f") FiltersFile    = argv [i+1]; 
            else if ( std::string(argv[i])  == "-g") input_Gain     = stof(argv [i+1]); 
           
        }
    }
    if (OutputFile == "-1") {std::cout << "No output file path and name provided, the program will be stopped.\n"; return -1; }

    // --> Create the configuration Object and initialise some useful variables
    ConfigFile C (ConfigFileName);
    int verbosity       = C.read<int> ("verbosity", 0);
    int doAllPulses     = C.read<int> ("doAllPulses", 0);
    float timeBin       = C.read<float>("gate", 204.800) / C.read<float>("sampling", 1024);

    float usr_Gain      = C.read<float>("Gain", -1);
    float usr_SPTR      = C.read<float>("SPTR", -1);
    float usr_QE_ScaleF = C.read<float>("QE_ScaleFactor", -1);

    std::string nEvent  = C.read<std::string> ("nameEv", "event");
    std::string nX      = C.read<std::string> ("nameX", "x");
    std::string nY      = C.read<std::string> ("nameY", "y");
    std::string nZ      = C.read<std::string> ("nameZ", "z");
    std::string nDetSim = C.read<std::string> ("nameDetSim", "detector");
    std::string nModSim = C.read<std::string> ("nameModSim", "module");
    std::string modType = C.read<std::string> ("nameModType", "moduleType");
    std::string nTime   = C.read<std::string> ("nameTime", "timestamp");
    std::string nTDepo  = C.read<std::string> ("nameTimeDepo", "t_deposition");
    std::string nTGen   = C.read<std::string> ("nameTimeGen ", "t_generation");
    std::string nTProp  = C.read<std::string> ("nameTimeProp", "t_propagation");
    std::string nEner   = C.read<std::string> ("nameEnergy", "PhotonEnergy");
    std::string nWave   = "PhotonWavelength";
    std::string nDetHit = "DetectorHit";
    std::string nModHit = "ModuleHit";
    std::string modTypeHit = "ModuleTypeHit";


    //
    // --> START
    int nThreads = C.read<int> ("multiThread", 0);
    if      (nThreads == 0) ROOT::EnableImplicitMT();                       // Enable Multithreading. Implicit evaluation of the number of threads
    else if (nThreads >  0) ROOT::EnableImplicitMT(nThreads);               // Enable Multithreading. Suggesting the number of threads.

    if (verbosity)      std::cout << "\n" << R"(
                                  
    ,o888888o.     8 888888888o   
 . 8888     `88.   8 8888    `88. 
,8 8888       `8b  8 8888     `88 
88 8888        `8b 8 8888     ,88 
88 8888         88 8 8888.   ,88' 
88 8888         88 8 888888888P'  
88 8888        ,8P 8 8888         
`8 8888       ,8P  8 8888         
 ` 8888     ,88'   8 8888         
    `8888888P'     8 8888         
        
)" <<
    "L.M. Apr 2020\n"
    << std::endl;



    // --> INITIALIZE THE RANDOM
    TRandom3 *randGen = new TRandom3 (C.read<int>("randomSeed", 0));         // Initialise the random object. It will be passed to the functions


    // --> UNDERSTAND THE GEOMETRY OF THE INPUT FILE
    //      27 May 2020   As a step forward for the simulation of the full LHCb ECAL, variable geometries are required.
    //      In this part the programs get from the input file the division into regions, modules and cells.
    int nModules = 0;
    std::map<int, std::vector<int>> moduleMap;
    std::vector<TH2I> geoMaps;
    GetGeometry (InputFile, nModules, moduleMap, geoMaps);
    // std::vector<int> moduleTypes, nFrontDetPerType, nBackDetPerType;
    // GetGeometry (InputFile, moduleTypes, nFrontDetPerType, nBackDetPerType);


    // --> CREATE THE DATAFRAME
    //      Let's assume that the .root file coming in is already as good as I need it.
    //      I can always prepare some other function to chain the different root files to place just before this step.
    if (verbosity) std::cout << "\n### Loading the Dataframe..." << std::endl;
    std::string TTreeName = C.read<std::string> ("TTreeName", "hybrid"); 
    ROOT::RDataFrame d (TTreeName.c_str(), InputFile.c_str());                      // Open the Dataframe
    std::vector<std::string> colNames   = d.GetColumnNames();                       // Get names of the columns



    // --> CREATE THE DATABASE
    //      This maps will be employed as a data structure to gather total charge and the waveform of each event,
    //      using as a key the event number.
    data_t data (C, nModules, moduleMap);



    // --> CREATE THE READOUT
    if (verbosity) std::cout << "### Initializing the readout..." << std::endl;
    // -- Check if the input gain is valid and in case replace the one from configuration file.
    if ( input_Gain >= 0. && !isnan( input_Gain )){
        std::cout << "\n### WARNING: replacing the gain from the configuration file with the command line input! ###\n\n";
        usr_Gain = input_Gain;
    }

    std::unique_ptr<LPhD>  photoDet;
    std::string PhD_Type = C.read<std::string>("PhD_Type", "Hama_R7899");
    if      (PhD_Type == "Hama_maPMT") {
        if ( usr_SPTR >= 0 || usr_QE_ScaleF >= 0)   photoDet = std::make_unique<LHamamatsuH13543>(timeBin, usr_Gain, usr_SPTR, usr_QE_ScaleF);
        else                                        photoDet = std::make_unique<LHamamatsuH13543>(timeBin, usr_Gain);
    }
    // else if (PhD_Type == "FBK_RGB")     photoDet = std::make_unique<LHamamatsuR7899>();/*(LHamamatsuR7899*) new LHamamatsuR7899();*/
    else if (PhD_Type == "Hama_R7899") {
        if ( usr_SPTR >= 0 || usr_QE_ScaleF >= 0)   photoDet = std::make_unique<LHamamatsuR7899> (timeBin, usr_Gain, usr_SPTR, usr_QE_ScaleF);
        else                                        photoDet = std::make_unique<LHamamatsuR7899> (timeBin, usr_Gain);
    }
    else if (PhD_Type == "Hama_R7600U_20") {
        if ( usr_SPTR >= 0 || usr_QE_ScaleF >= 0)   photoDet = std::make_unique<LHamamatsuR7600U_20> (timeBin, usr_Gain, usr_SPTR, usr_QE_ScaleF);
        else                                        photoDet = std::make_unique<LHamamatsuR7600U_20> (timeBin, usr_Gain);
    }
    else if (PhD_Type == "Ideal")       { std::cout << "### Ideal readout selected! ###\n"; photoDet = std::make_unique<LPhD>(timeBin); }      
    else {
        if ( usr_SPTR >= 0 || usr_QE_ScaleF >= 0)   photoDet = std::make_unique<LHamamatsuR7899> (timeBin, usr_Gain, usr_SPTR, usr_QE_ScaleF);
        else                                        photoDet = std::make_unique<LHamamatsuR7899> (timeBin, usr_Gain);
    }

    if (verbosity) photoDet -> PrintInfo();


    // --> EXTRACT THE FILTERS
    //      Open the .root file with the filters and retrieve the filters.
    //      Let's extract them.
    // TFile * FiltersRootFile = new TFile(FiltersFile.c_str(),"READ");
    // if (!FiltersRootFile->IsOpen()) std::cout << "### WARNING: Filters file " << FiltersFile << " not open. No spatial filters!" << std::endl;

    TGraph2D *SEGraph   = 0;                                                
    std::vector<TCutG> *Borders = NULL;                                     // Detector selection --- NOT USED!

    // ******************************** TEMPORARY SOLUTION ****************************
    float UnifLossFac   = C.read<float>("UnifLossFac", 1.);
    float UnifSpatEff   = 1. / UnifLossFac;
    std::cout << "### Using a uniform loss factor (photons / photons detected) = " << UnifLossFac << " that is a spatial efficiency = " << UnifSpatEff << "\n\n\n";
    // ******************************** TEMPORARY SOLUTION ****************************

    // FiltersRootFile -> GetObject(C.read<std::string>("SEGraphName","SpatialEfficiencyGraph").c_str(), SEGraph); // Extract the spatial efficiencies

    // if (!SEGraph)               std::cout << "### WARNING: No spatial efficiency given. Efficiency = 1 ###" << std::endl; 
    // if (!Borders)               std::cout << "### WARNING: No detectors borders given. The label present in the input file will be used. ###" << std::endl;    









    // ###################################
    // ---> B O O K   T H E   L O O P --->
    //      But actually do it too, 'cause I am using Foreach.
    if (verbosity) std::cout << "\t.Starting the analysis." << std::endl;

    auto evTonm = [] (float & energy) ->float{ return float(1240.) / float(energy); };            // Lambda to convert from energy to nanometers


    auto dFilt = d.Define ( nWave.c_str(), evTonm, {nEner.c_str()} )                // Do all the filters!
                    .Filter ( [&UnifSpatEff, &SEGraph, &randGen] (float& x, float& y, float& z) ->bool{ return ApplySpatialFilters(UnifSpatEff, SEGraph, x, y, z, randGen);},                   {nX.c_str(), nY.c_str(), nZ.c_str()} )                    // Apply spatial efficiency
                    // .Define ( nModHit.c_str(), [&Borders] (float& x, float& y, float& z, int& nDetSim) ->int{return AssignDetector(x, y, z, nDetSim, Borders);},     {nX.c_str(), nY.c_str(), nZ.c_str(), nDetSim.c_str()} )   // Assign each photon to a detector
                    // .Define ( nDetHit.c_str(), [&Borders] (float& x, float& y, float& z, int& nDetSim) ->int{return AssignDetector(x, y, z, nDetSim, Borders);},     {nX.c_str(), nY.c_str(), nZ.c_str(), nDetSim.c_str()} )   // Assign each photon to a detector
                    // .Filter ( [&data] (int & nDet) ->bool{return ( nDet >= 0 && nDet < data.Get_nDetectors());},    {nDetHit.c_str()})                                        // Filter the event not assigned to any of the detector we are looking for
                    .Define ( nModHit.c_str(), []   (int& nModSim) ->int{return nModSim;},     {nModSim.c_str()} )   // Assign each photon to a module
                    .Define ( modTypeHit.c_str(),[] (int& modType) ->int{return modType;},     {modType.c_str()} )   // Get the module type
                    .Define ( nDetHit.c_str(), []   (int& nDetSim) ->int{return nDetSim;},     {nDetSim.c_str()} )   // Assign each photon to a detector
                    .Filter ( [&data] (int & nMod) ->bool{return ( nMod >= 0 && nMod < data.Get_nModules());},      {nModHit.c_str()})                                        // Filter the event not assigned to any of the detector we are looking for
                    // .Filter ( [&data] (int & nDet) ->bool{return ( nDet >= 0 && nDet < data.Get_nDetectors());},    {nDetSim.c_str()})                                        // Filter the event not assigned to any of the detector we are looking for
                    .Filter ( [&photoDet, &randGen] (float& wave, int & nDet) ->bool{ if (randGen->Rndm() < photoDet->Get_QE(wave)) return true; else return false;},{nWave.c_str(), nDetHit.c_str()});                        // Apply QE of each detector - nDetHit not used but left for future upgrades, if necessary





    // --- DO IT
    // if (doAllPulses) dFilt.Foreach( [&data, &photoDet, &randGen] (int & eventNum, float & timeStamp, int & nDet, float & t_dep, float & t_gen, float & t_prop)
    //                     { data.UpdateMapAllPulses(eventNum, timeStamp, t_dep, t_gen, t_prop, nDet, *photoDet, randGen); },
    //                     { nEvent.c_str(), nTime.c_str(), nDetHit.c_str(), nTDepo.c_str(), nTGen.c_str(), nTProp.c_str() } );


/*    else  */           dFilt.Foreach([&data, &photoDet, &randGen] (int & eventNum, float & timeStamp, int & nMod, int & modType, int & nDet)
                        { data.UpdateMap(eventNum, timeStamp, nMod, modType, nDet, *photoDet, randGen); },
                        { nEvent.c_str(), nTime.c_str(), nModHit.c_str(), modTypeHit.c_str(), nDetHit.c_str() } );
    
    // <-- ANALYSIS - DONE <--

    if (verbosity) std::cout << "\t.Done." << std::endl;
    // #######################









    if (verbosity) data.PrintTheMap();

    if (C.read <int> ("saveFiltered", 0)) { std::cout << "SNAPSHOT!" << std::endl; dFilt.Snapshot(TTreeName.c_str(), (OutputFile + "_FilteredSnapshot.root").c_str(), dFilt.GetColumnNames()); }











    // #######################
    // --> SAVE IN A TFILE -->
    std::cout << "### ...Exporting data --> TTree at " << (OutputFile + ".root").c_str() << "\n"; 
    TFile *outFile = new TFile ((OutputFile + ".root").c_str(), "RECREATE");
    TTree *outTree = new TTree (C.read<std::string>("OutTreeName", "DATA").c_str(), C.read<std::string>("OutTreeName", "DATA").c_str());

    // if (doAllPulses) data.SaveToTTreeAllPulses(outTree);    // Save to TTree the map 
/*    else  */           data.SaveToTTree(outTree);             // Save to TTree the map with all the additional pulses

    // --- Save facultative material
    if (C.read<int> ("saveInfo", 1)) {
        if (verbosity) std::cout << "### Saving QE and SPR graphs...\n";
        photoDet -> Get_QEGraph() -> Write("QE");
        photoDet -> Get_SPRGraph() -> Write("SPR");
        photoDet -> Get_fSPR_AmplitudeDistr() -> Write("SPR_AmplitudeDistribution");
        photoDet -> Get_hSPTR() -> Write("_hSPTR");

        for (auto && map : geoMaps) map.Write();
    }
    // --> Close Everything
    outTree -> Write();
    outFile -> Close();

    
    // <-- SAVE - DONE <--
    // #######################


    // --- Just a recap of the photodetector characteristics, might be useful for debugging
    if (verbosity) {
        std::cout << "Printing again the photodetector's information...\n";
        photoDet -> PrintInfo();
    }

    if (verbosity) std::cout << "\nStandard end of the program.\n\t# Farewell #" << std::endl;
    return 0;
}








//
// ############### I M P L E M E N T A T I O N ###############
//
data_t::data_t (const ConfigFile & C, const int & nModules, std::map<int, std::vector<int>> & moduleMap) : _nModules (nModules),
                                                            // _moduleTypes      (infoGeometry.at(0)),
                                                            // _nFrontDetPerType (infoGeometry.at(1)),
                                                            // _nBackDetPerType  (infoGeometry.at(2)),
                                                            _sampling   (C.read<int>("sampling", 1024)),
                                                            _gate       (C.read<float>("gate", 204.800)),  
                                                            _delay      (C.read<float>("delay", 0)),
                                                            _baseline   (C.read<float>("baseline", 0.)),
                                                            _sigmaNoise (C.read<float>("sigmaNoise", 0.)),
                                                            _timeBin    ( float(_gate) / float(_sampling) )                                        
{
    std::cout << "Modules found: " << _nModules << std::endl;

    // --- Initialise the map
    for (auto& iMod : moduleMap) {
        _moduleMap [iMod.first] = iMod.second;
        std::cout << "Module Type: " << _moduleMap.find(iMod.first)->first << "\t # Cells front: " << _moduleMap[iMod.first].at(1) << " back: " << _moduleMap[iMod.first].at(1) << std::endl;
    }

    // --- Initialise the strings
    _nBPhot = C.read<std::string> ("OutPhotName", "pmt");
    _nBWave = C.read<std::string> ("OutWaveName", "pulse");
    _nBMod  = C.read<std::string> ("OutModName", "mod");
    _nBWaveDepo     = C.read<std::string> ("OutWaveNDepo",     "pulse_deposition");
    _nBWaveScint    = C.read<std::string> ("OutWaveNScint",    "pulse_depoAndScint");
    _nBWaveDepoTrans= C.read<std::string> ("OutWaveNDepoTran", "pulse_depoAndTrans");
};








void data_t::UpdateMap (const int &nEv, const float &nTime, const int & nMod, const int & modType, const int & nDet, const LPhD & photoDet, TRandom3 *randGen)
{
    // --- Look for the event in the map
    auto nEvPos = _M.find(nEv);


    // --- If the event is new, add it to the map and return an iterator to it;
    if ( nEvPos == _M.end() ) {
        // std::cout << "### NEW EVENT FOUND! iev: " << nEv << std::endl;
        // std::cout << "Initialising " << _nModules << " modules" << std::endl;
        event_t newEntry;

        newEntry.nEvent = nEv;
        newEntry.nPhotonsDet        = std::vector <std::vector<int>> (_nModules);
        newEntry.detectorTimeStamps = std::vector <std::vector<std::vector<float>>> (_nModules);
        newEntry.digiPulse          = std::vector <std::vector<std::vector<float>>> (_nModules);//, std::vector<float> (_sampling, _baseline));


        nEvPos = _M.insert ( std::make_pair (nEv, newEntry) ).first;           // FIND returns a pair, where the first is the iterator pointing at the position of insertion
        // std::cout << "### NEW EVENT Inserted!" << std::endl;
    }


    // --- If the module has never been hit before, initialise all the cells in it.
    //      The check is performed looking at the size of the vectors of photons detected 
    //      in each cell of the nMod module.
    if ( nEvPos->second.nPhotonsDet.at(nMod).size() <= 0 ) {
        // std::cout << "### ### NEW MODULE HIT! " << " nMod: " << nMod  << " moduleType: " << modType << std::endl;
        // std::cout << "                   " << " # modules: " << _nModules << std::endl;
        // -- Add it to the modules hit
        // std::cout << "Adding to the modules hit" << std::endl;
        nEvPos->second.modulesHit.emplace_back(nMod);

        // -- Initialise the module's cells
        // std::cout << "Getting number of cells" << std::endl;
        auto nOfCells = this->Get_nCellFront(modType) + this->Get_nCellBack(modType) ;
        // std::cout << "modType: " << modType << " nOfCells: " << nOfCells << std::endl;
        // std::cout << "Initialise photons per cell" << std::endl;
        nEvPos->second.nPhotonsDet.at(nMod)           = std::vector <int> (nOfCells, 0);
        // std::cout << "Initialise detector timestamps" << std::endl;
        nEvPos->second.detectorTimeStamps.at(nMod)    = std::vector <std::vector<float>> (nOfCells);
        // std::cout << "Initialise pulses" << std::endl;
        nEvPos->second.digiPulse.at(nMod)             = std::vector <std::vector<float>> (nOfCells, std::vector<float> (_sampling, _baseline));

        // -- Add random electronic noise to the pulses of the module
        for ( int iDet=0; iDet < nEvPos->second.digiPulse.at(nMod).size(); ++iDet) {
            for (int iBin=0; iBin < _sampling; ++iBin) {
                nEvPos->second.digiPulse.at(nMod).at(iDet).at(iBin) += randGen -> Gaus (0, _sigmaNoise);
            }
        }
        // std::cout << "### ### NEW MODULE Created." << std::endl;
    }


    // std::cout << "Updating the event: nMod " << nMod << "   nDet " << nDet << std::endl;
    // --- UPDATE THE EVENT
    // -> Tot Photons ->
    nEvPos->second.nPhotonsTot ++;                            // Tot photons

    // -> Photons per Cell ->
    nEvPos->second.nPhotonsDet.at(nMod).at(nDet) ++;

    // -> SinglePhel response ->
    float randSPTR = randGen -> Gaus (0, photoDet.Get_SPTR());  // SPTR
    float randAmpl = photoDet.GetRand_Amplitude();              // Simulate Single Photoelectron Pulse's Amplitude
    FillPulse( nEvPos -> second.digiPulse.at(nMod).at(nDet), nTime, photoDet, randSPTR, randAmpl);
}






/*
void data_t::UpdateMapAllPulses (const int &nEv, const float &nTime, const float & t_dep, const float & t_gen, const float & t_prop, const int & nDet, const LPhD & photoDet, TRandom3 *randGen)
{
    // ###################################################################################
    // ### The idea is to apply the same time smearing to all the different timestamps
    // ######

    auto nEvPos = _M.find(nEv);     // Look for the event in the map

    // --> If it does not exist, create it! -->
    if ( nEvPos == _M.end() )       // If it is not present, create it!
    {
        // --- Initialise the entry
        event_t newEntry;

        // std::cout << "New event: " << nEv << "\n";   // DEBUG
        newEntry.nEvent = nEv;
        newEntry.nPhotonsDet        = std::vector <int> (_nDetectors, 0);
        newEntry.detectorTimeStamps = std::vector <std::vector<float>> (_nDetectors);
        newEntry.digiPulse          = std::vector <std::vector<float>> (_nDetectors, std::vector<float> (_sampling, _baseline));
        newEntry.digiPulseDepo      = std::vector <std::vector<float>> (_nDetectors, std::vector<float> (_sampling, _baseline));
        newEntry.digiPulseScint     = std::vector <std::vector<float>> (_nDetectors, std::vector<float> (_sampling, _baseline));
        newEntry.digiPulseDepoTrans = std::vector <std::vector<float>> (_nDetectors, std::vector<float> (_sampling, _baseline));

        // -- Add random electronic noise
        float noise = 0;
        for ( int iDet=0; iDet < _nDetectors; iDet++)
        {
            for (int iBin=0; iBin < _sampling; iBin++)
            {
                noise = randGen -> Gaus (0, _sigmaNoise);
                newEntry.digiPulse.at(iDet).at(iBin)          += noise;
                newEntry.digiPulseDepo.at(iDet).at(iBin)      += noise;
                newEntry.digiPulseScint.at(iDet).at(iBin)     += noise;
                newEntry.digiPulseDepoTrans.at(iDet).at(iBin) += noise;
            }
        }

        // ---> Update the event -->
        // -> Tot Photons ->
        newEntry.nPhotonsTot ++;                                    // Tot photons

        // -> Timestamps ->
        newEntry.detectorTimeStamps.at(nDet).emplace_back(nTime);   // timestamps

        // -> SinglePhel response ->
        float randSPTR = randGen -> Gaus (0, photoDet.Get_SPTR());  // SPTR
        float randAmpl = photoDet.GetRand_Amplitude();              // Simulate Single Photoelectron Pulse's Amplitude
        FillPulse( newEntry.digiPulse.at(nDet)         , nTime       , photoDet, randSPTR, randAmpl );
        FillPulse( newEntry.digiPulseDepo.at(nDet)     , t_dep       , photoDet, randSPTR, randAmpl );
        FillPulse( newEntry.digiPulseScint.at(nDet)    , t_dep+t_gen , photoDet, randSPTR, randAmpl );
        FillPulse( newEntry.digiPulseDepoTrans.at(nDet), t_dep+t_prop, photoDet, randSPTR, randAmpl );




        _M.insert ( std::make_pair (nEv, newEntry) );

    }
    // <-- If it does not exist, create it! <--

    // --> Else, update it! -->
    else
    {
        // Update the event
        nEvPos -> second.nPhotonsTot ++;                                // Tot photons 
        nEvPos -> second.detectorTimeStamps.at(nDet).emplace_back(nTime);   // timestamps

        // -> SinglePhel response ->
        float randSPTR = randGen -> Gaus (0, photoDet.Get_SPTR());  // SPTR
        float randAmpl = photoDet.GetRand_Amplitude();              // Simulate Single Photoelectron Pulse's Amplitude
        FillPulse( nEvPos -> second.digiPulse.at(nDet)         , nTime       , photoDet, randSPTR, randAmpl );
        FillPulse( nEvPos -> second.digiPulseDepo.at(nDet)     , t_dep       , photoDet, randSPTR, randAmpl );
        FillPulse( nEvPos -> second.digiPulseScint.at(nDet)    , t_dep+t_gen , photoDet, randSPTR, randAmpl );
        FillPulse( nEvPos -> second.digiPulseDepoTrans.at(nDet), t_dep+t_prop, photoDet, randSPTR, randAmpl );
    }
    // <-- Else, update it! <--
    
}





*/





void data_t::FillPulse(std::vector<float> & Pulse, const float &nTime, const LPhD & photoDet, const float & randSPTR, const float & randAmpl)
{
    float timeStamp = nTime + _delay;                           // Overall constant delay
    timeStamp += randSPTR;                                      // Detector SPTR

    int binsOn = photoDet.Get_pulseLen()/_timeBin;              // bins turned on by the pulse
    int nBin = timeStamp / _timeBin + 1;                        // Start always from the following bin to be > than timeStamp
    for (int iBin = 0; iBin < binsOn; iBin ++)                  // sample the pulse. iBin < binsOn 'cause you start from the timestamp bin + 1
    {
        if (nBin + iBin > _sampling-1 || nBin + iBin < 0) continue;
        Pulse.at(nBin + iBin) += randAmpl * photoDet.sPhelPulseInterp( (nBin + iBin) * _timeBin - timeStamp );
    }

}




void data_t::SaveToTTree(TTree* TT)
{
    int LYTot = 0;
    int nModulesHit = 0;
    std::vector<int> modulesHit {0};
    std::vector<std::vector<int>> nPhotonsDet (_nModules);
    std::vector<std::vector<std::vector<float>>> digiPulse (_nModules);

    // --- Create the branches
    TT -> Branch ("Total_Light", &LYTot, "Total_Light/I");
    TT -> Branch ("modulesHit", &modulesHit);

    for (int iMod = 0; iMod < _nModules; ++iMod)
    {
        TT -> Branch ( (_nBMod + std::to_string(iMod) + "_" + _nBPhot).c_str(), &(nPhotonsDet.at(iMod)) );
        TT -> Branch ( (_nBMod + std::to_string(iMod) + "_" + _nBWave).c_str(), &(digiPulse.at(iMod)) );
    }


    // ---> Fill the TTree
    for (auto&& iEv : _M)
    {
        LYTot = iEv.second.nPhotonsTot;
        modulesHit  = iEv.second.modulesHit;
        for (int iMod = 0; iMod < _nModules; ++iMod)
        {
            nPhotonsDet.at(iMod) = iEv.second.nPhotonsDet.at(iMod);
            digiPulse.at(iMod) = iEv.second.digiPulse.at(iMod);
        }

        TT -> Fill();
    }


}








/*
void data_t::SaveToTTreeAllPulses(TTree* TT)
{
    // Initialize the TTree
    int LYTot = 0;
    std::vector<std::vector<float>> pulseRow            (_nDetectors, std::vector<float> (_sampling, 0)); 
    std::vector<std::vector<float>> pulseRowDepo        (_nDetectors, std::vector<float> (_sampling, 0)); 
    std::vector<std::vector<float>> pulseRowScint       (_nDetectors, std::vector<float> (_sampling, 0)); 
    std::vector<std::vector<float>> pulseRowDepoTrans   (_nDetectors, std::vector<float> (_sampling, 0)); 
    std::vector<int> LYRow                    (_nDetectors, 0);                  

    TT -> Branch ("Total_Light", &LYTot, "Total_Light/I");
    for (int iDet = 0; iDet < _nDetectors; iDet++)
    {
        TT -> Branch ( (_nBPhot + std::to_string(iDet)).c_str()   , &LYRow.at(iDet), (_nBPhot + std::to_string(iDet) + "/I").c_str() );
        TT -> Branch ( (_nBWave + std::to_string(iDet)).c_str()           , pulseRow.at(iDet).data()          , (_nBWave + std::to_string(iDet) + "[" + std::to_string(_sampling) + "]/F").c_str());
        TT -> Branch ( (_nBWaveDepo     + std::to_string(iDet)).c_str()   , pulseRowDepo.at(iDet).data()      , (_nBWaveDepo     + std::to_string(iDet) + "[" + std::to_string(_sampling) + "]/F").c_str());
        TT -> Branch ( (_nBWaveScint    + std::to_string(iDet)).c_str()   , pulseRowScint.at(iDet).data()     , (_nBWaveScint    + std::to_string(iDet) + "[" + std::to_string(_sampling) + "]/F").c_str());
        TT -> Branch ( (_nBWaveDepoTrans+ std::to_string(iDet)).c_str()   , pulseRowDepoTrans.at(iDet).data() , (_nBWaveDepoTrans+ std::to_string(iDet) + "[" + std::to_string(_sampling) + "]/F").c_str());
    }

    for (auto&& iEv : _M)
    {
        LYTot = iEv.second.nPhotonsTot;
        for ( int iDet=0; iDet < iEv.second.detectorTimeStamps.size(); iDet++) LYRow.at(iDet) = iEv.second.detectorTimeStamps.at(iDet).size();
        pulseRow          = iEv.second.digiPulse;
        pulseRowDepo      = iEv.second.digiPulseDepo;
        pulseRowScint     = iEv.second.digiPulseScint;
        pulseRowDepoTrans = iEv.second.digiPulseDepoTrans;

        TT -> Fill();
    }

}
*/





void GetGeometry (const std::string & inputFileName, int & nModules, std::map<int, std::vector<int>> & moduleMap,  std::vector<TH2I> & geoMaps)
{
    std::cout << "\n### Retrieving geometrical information from the input file...\n";
    TFile fileIn (inputFileName.c_str(), "READ");
    TTree *modules = (TTree*) fileIn.Get("modules");


    std::vector<int> moduleTypes;
    std::vector<int> nFrontDetPerType;
    std::vector<int> nBackDetPerType;

    auto modules_regions = CalculateRegions(modules);

    // --- Find how many module types are involved
    for (int iMod = 0; iMod < modules_regions.size(); ++iMod) {
        bool moduleTypeIsAlreadyThere = false;

        for (size_t im = 0; im < moduleTypes.size(); ++im) {
            if (modules_regions[iMod].type == moduleTypes[im]) moduleTypeIsAlreadyThere = true;
        }

        if (!moduleTypeIsAlreadyThere) moduleTypes.push_back(modules_regions[iMod].type);
    }



    // --- Count the module types
    int Nmodules  = modules_regions.size();
    nModules = Nmodules;
    int NmodTypes = moduleTypes.size();
    
    std::vector <int> nModulePerType(NmodTypes, 0);

    nFrontDetPerType.resize(NmodTypes);
    nBackDetPerType .resize(NmodTypes);

    // --- Save the moduleType map
    auto modMap = (TH2I*) fileIn.Get("module_ID_map");
    if (modMap) {
        geoMaps.emplace_back(TH2I (*modMap));
        delete modMap;
    }

    // --- Get info on the cells    
    for (int i = 0; i < NmodTypes; ++i)
    {
        std::stringstream snameFront;
        std::stringstream snameBack;
        snameFront  << "cell_map_front_" << moduleTypes[i];
        snameBack   << "cell_map_back_"  << moduleTypes[i];
        nFrontDetPerType[i] = ((TH2I*) fileIn.Get(snameFront.str().c_str()))->GetEntries();
        nBackDetPerType [i] = ((TH2I*) fileIn.Get(snameBack .str().c_str()))->GetEntries();

        // // --- Save the info
        geoMaps.emplace_back(TH2I (*(TH2I*) fileIn.Get(snameFront.str().c_str())) );
        geoMaps.emplace_back(TH2I (*(TH2I*) fileIn.Get(snameBack.str().c_str()))  );
    }



    // --- Pass the info to the map and print it 
    for(size_t im = 0; im < moduleTypes.size(); ++im)
    {
        std::cout << "Module Type n " << moduleTypes[im] << " has " << nFrontDetPerType[im] << " modules front" << std::endl;
        std::cout << "Module Type n " << moduleTypes[im] << " has " << nBackDetPerType[im]  << " modules back"  << std::endl;

        moduleMap[ moduleTypes[im] ] = std::vector<int> {nFrontDetPerType[im], nBackDetPerType[im]};
    }


    delete modules;
    fileIn.Close();


}
