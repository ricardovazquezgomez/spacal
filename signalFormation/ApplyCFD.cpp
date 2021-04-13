// #################################################################################################
// ### Program that applies CFD triggering to the pulses coming out of OP/simReadout
// ### Actually one can apply almost any trigger, they just need to modify the ApplyDefines function
// ### adding some other functions instead of TriggerCFD!!!
// ######
//
#include <iostream>
#include <string>
#include <vector>

#include "ConfigFile.hh"

#include "TROOT.h"
#include "TCanvas.h"
#include "TLine.h"
#include "ROOT/RDataFrame.hxx"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TStyle.h"
#include "TKey.h"

struct triggerInfo_t
{
    // ################################################################
    // ### A nice struct to help me recursively create the RDF columns.
    // ######
    int sampling;
    int baseLineSamples;
    int fitSamplesNum;
    float timesBase;
    int polarity;
    int zeroTrigger;
    float peakFraction;
    float timeBin;

    int convert_to_ns;
};



void PrintUsage () { std::cout << " Usage: " << std::endl; std::cout << " ./ApplyCFD [-c ConfigFile] [-i DataInputFile] [-o OutputFile without extension] [-t peak fraction for threshold]" << std::endl; }
//double interpolatePulse(float* data, unsigned int length, int polarity, int Tstart,int baseLineSamples, int samplesNum, /*float* res, float &area,*/ float timesBase /*,double cutoffBaseline,int ignore*/);
double TriggerCFD (float* data, unsigned int length, int polarity, int tStart, int baseLineSamples, int samplesNum, float timesBase, float peakFraction);
std::vector<std::string> Whitelist(std::vector<std::string> names, std::vector<std::string> blackList);
void GetGeometry (const std::string & inputFileName, std::vector<TH2I> & geoMaps);

















int main (int argc, char** argv)
{
    // --> Read the Command Line
    std::string ConfigFileName = "ConfigFile.cfg", InputFile = "-1", OutputFile = "-1";
    float input_peakFraction =  -1.;
    if  ( argc < 2 ) {
        PrintUsage();
        return -1; }

    else {
        for (int i = 0; i < argc; i++)
        {
            if      ( std::string(argv[i])  == "-c") ConfigFileName      = argv [i+1];
            else if ( std::string(argv[i])  == "-i") InputFile           = argv [i+1];
            else if ( std::string(argv[i])  == "-o") OutputFile          = argv [i+1];
            else if ( std::string(argv[i])  == "-t") input_peakFraction  = stof(argv [i+1]);
        }
    }
    if (InputFile  == "-1") {std::cout << "### ERROR!!! No input  file specified!\n"; return -1;}
    if (OutputFile == "-1") {std::cout << "### ERROR!!! No output file specified!\n"; return -1;}





    // --> Create the configuration Object and initialise some useful variables
    ConfigFile C (ConfigFileName);

    int verbosity           = C.read<int> ("verbosity", 0);

    // --> Let's do this.
    if (verbosity)      std::cout << "\n" << R"(
 _______  ______ _____  ______  ______ _______  ______
    |    |_____/   |   |  ____ |  ____ |______ |_____/
    |    |    \_ __|__ |_____| |_____| |______ |    \_

)" <<
    "L.M. Apr, Jun 2020\n\n";








    // --- Initialise variables
    int nThreads            = C.read<int> ("multiThread", -1);
    int savePulse           = C.read<int> ("savePulse", 0);
    int skimPulses          = C.read<int> ("skimPulses", 1);

    float timeBin           = C.read<float>("gate", 204.800) / C.read<float>("sampling", 1024);
    float peakFraction      = C.read<float> ("peakFraction", 0.5);

    // -- Check if the input peak fraction is valid and in case replace the one from configuration file.
    if ( input_peakFraction >= 0. && input_peakFraction <= 1. && !isnan( input_peakFraction )){
        std::cout << "\n### Replacing the peak fraction from the configuration file with the command line input! ###\n\n";
        peakFraction = input_peakFraction;
    }

    triggerInfo_t trig          = {C.read<int> ("sampling", 1024), C.read<int> ("baseLineSamples", 20),
                                    C.read<int> ("fitSamplesNum", 3), C.read<float> ("timesBase", 10),
                                    C.read<int> ("PulseSign", -1), C.read<int> ("zeroTrigger", 0.),
                                    peakFraction, timeBin,
                                    C.read<int> ("convert_to_ns", 0)};

    std::string nBMod           = C.read<std::string> ("OutModName", "mod");
    std::string nBTimeS         = C.read<std::string> ("OutTimeStamp", "t");
    std::string nBWave          = C.read<std::string> ("OutWaveName", "pulse");
    // std::string nBWaveDepo      = C.read<std::string> ("OutWaveName", "pulse_deposition");
    // std::string nBWaveScint     = C.read<std::string> ("OutWaveName", "pulse_depoAndScint");
    // std::string nBWaveDepoTrans = C.read<std::string> ("OutWaveName", "pulse_depoAndTrans");
    std::string TTreeName       = C.read<std::string> ("inTreeName", "tree");
    std::string outTTreeN       = C.read<std::string> ("OutTreeNameCFD", "tree");




    // --> FEEDBACK
    if (verbosity)
    {
        std::cout << "#### F E E D B A C K ###\n";
        std::cout << "Input waveform branches prefix      : " << nBWave << "\n";
        std::cout << "Output timestamp branches prefix    : " << nBTimeS << "\n";
        std::cout << "Input TTree name                    : " << TTreeName << "\n";
        std::cout << "Output TTree name                   : " << TTreeName << "\n";
        std::cout << "Verbosity                           : " << TTreeName << "\n";
        std::cout << "Number of threads                   : " << nThreads << "\n";
        std::cout << "Save pulses of event                : " << savePulse << "\n";
        std::cout << "Skim pulse branches [Y/N]           : " << skimPulses << "\n";
        std::cout << "------------------------------------\n";
        std::cout << "Sampling bins                       : " << trig.sampling << "\n";
        std::cout << "Samples to get the baseline         : " << trig.baseLineSamples << "\n";
        std::cout << "Samples for the linear fit (+/-)    : " << trig.fitSamplesNum << "\n";
        std::cout << "Minimum threshold in sigma units    : " << trig.timesBase << "\n";
        std::cout << "Pulse polarity                      : " << trig.polarity << "\n";
        std::cout << "Starting point in samples units     : " << trig.zeroTrigger << "\n";
        std::cout << "Constant fraction of the pulse peak : " << trig.peakFraction << "\n";
        std::cout << "Convert timestamps to time          : " << trig.convert_to_ns << "\n";
        std::cout << "Samples to time conversion factor   : " << trig.timeBin << "\n";
        std::cout << "\n";
    }

    if (trig.convert_to_ns) std::cout << "\n### WARNING: Timestamps will be calculated in time units (conversion factor: " << trig.timeBin << ")! ###\n";
    else                    std::cout << "\n### WARNING: Timestamps will be calculated in digitizer's clock units! ###\n";



    // --> Get the maps
    std::vector<TH2I> geoMaps;
    GetGeometry (InputFile, geoMaps);

    // --> Open the TFile
    TFile fin (InputFile.c_str(), "READ");
    TFile fout((OutputFile+".root").c_str(), "RECREATE");
    TTree *TTin = (TTree*) fin.Get(TTreeName.c_str());







    // --> Prepare the vectors of strings with the name of the branches
    if (verbosity) std::cout << "... Looking for all the branches starting with \"" << nBMod << "\" containing \"" << nBWave << "\"" << std::endl;
    std::vector<std::string> colNames;

    auto listNames = TTin -> GetListOfBranches();
    for (int i = 0; i < listNames -> GetEntries(); ++i) colNames.push_back(listNames -> At(i) -> GetName());

    std::vector <std::string> newColNames, pulseColNames;
	for (auto &&colName : colNames)	{
		if ( (!colName.compare(0, nBMod.size(), nBMod)) && (colName.find(nBWave) != std::string::npos)  ) {
			pulseColNames.push_back(colName);
            newColNames.push_back(colName.erase(colName.find(nBWave), nBWave.size()) + nBTimeS);
        }
    }





    // --- Assign branches
    std::vector<TBranch*> bTimestamps;
    std::vector<std::vector<float>> mTimestamps (newColNames.size());
    std::vector< std::vector<std::vector<float>>* > modulesPulses (pulseColNames.size(), 0);

    for (size_t iB = 0; iB < pulseColNames.size(); ++iB) {
        TTin -> SetBranchAddress(pulseColNames.at(iB).c_str(), &(modulesPulses.at(iB)));                // <-- Old branches
        bTimestamps.emplace_back( TTin -> Branch(newColNames.at(iB).c_str(), &mTimestamps.at(iB)) );    // <-- New branhces
    }




    // --- START the reading
    std::cout << "\n\t... Starting the analysis ..." << std::endl;
    int nEv = TTin -> GetEntries();
    for (int iEv = 0; iEv < nEv; ++iEv)                                                     // -> Event
    {
        // --- Start working on a new event
        std::cout << "Event: " << iEv << std::endl;
        TTin -> GetEvent(iEv);









        // --- Get the timestamps
        for (size_t iM = 0; iM < pulseColNames.size(); ++iM)                                    // -> Module
        {
            mTimestamps.at(iM).clear();

            for (size_t iCell = 0; iCell < modulesPulses.at(iM)->size(); ++iCell) {             // -> Cell



                if (trig.convert_to_ns) mTimestamps.at(iM).emplace_back( trig.timeBin * TriggerCFD (&modulesPulses.at(iM)->at(iCell)[0],
                                                                        trig.sampling, trig.polarity, trig.zeroTrigger, trig.baseLineSamples,
                                                                        trig.fitSamplesNum, trig.timesBase, trig.peakFraction) );



                else                    mTimestamps.at(iM).emplace_back(                TriggerCFD (&modulesPulses.at(iM)->at(iCell)[0],
                                                                        trig.sampling, trig.polarity, trig.zeroTrigger, trig.baseLineSamples,
                                                                        trig.fitSamplesNum, trig.timesBase, trig.peakFraction) );



            } // <- Cell

            bTimestamps.at(iM) -> Fill();
        } // <- Module


    } // <- Event
    std::cout << " ... done." << std::endl;






    // --- Now Write() to a new file the TTree
    //      If desired (configfile) skim the pulse branches.
    //
    //      NOTE: This is quite the weak spot of this little program. First cloning and then writing.
    //          Does this include 2 loops over the data? But on the other hand, the only alternative would be
    //          to store all the other branches in memory whilst looping over the pulses, which implies knowing
    //          all their types (does it?) and a consequent loss of generality.
    //          So I'll leave it as it is for the time being, I 'spose.
    //
    std::cout << "\n... Writing to file..." << std::endl;
    if (skimPulses) for ( auto && column : pulseColNames) TTin -> SetBranchStatus (column.c_str(), 0);

    auto TTout = TTin->CloneTree(0);
    TTout->CopyEntries(TTin);

    fout.cd();
    for (auto && map : geoMaps) map.Write();
    TTout -> Write();
    fout.Close();
    fin.Close();




    if (verbosity) std::cout << "Data saved to " << (OutputFile+".root").c_str() << "\n";
    // <-- DONE. <--

    if (verbosity) std::cout << ".Standard end of the program.\n";
    if (verbosity) std::cout << "# # #   F A R E W E L L   # # #\n";
    return 0;

}




















double TriggerCFD (float* data, unsigned int length, int polarity, int tStart, int baseLineSamples, int samplesNum, float timesBase, float peakFraction)
{
    // ###############################################################################
    // ### Function to do CFD triggering. It returns a double which is the timestamp.
    // ### Modified version of the one by M. Pizzichemi, CERN.
    // ### data                 --> array of data samples
    // ### length               --> length of the array "data"
    // ### polarity             --> polarity of the pulse, either > 0 or < 0
    // ### tStart               --> bin from which the analysis start, everything below is ignored
    // ### baseLineSamples      --> How many samples employ to calculate the baseline
    // ### samplesNum           --> How many samples before and after the threshold employ to do the linear regression of the crosspoint
    // ### timesBase            --> How many times the sigmaBaseline to detect the start of the pulse
    // ### peakFraction         --> Fraction of the pulse height at which set the threshold to do CFD.
    // ######

    unsigned int i,j;
    double crosspoint    = -1.0;
    double baseline      = 0.0;
    double sigmaBaseline = 0;
    double minBaseline   = FLT_MAX;
    double maxBaseline   = -FLT_MAX;

    // -- Mean Calculation
    for ( i = int(tStart); i < baseLineSamples; i++) baseline += float(data[i]) / float(baseLineSamples);

    // -- RMS
    for ( i = int(tStart); i < baseLineSamples; i++) sigmaBaseline += pow(float(data[i]) - baseline, 2 );
    sigmaBaseline = fabs( sqrt( sigmaBaseline / float(baseLineSamples)) );


    int pulseStartPoint = 0;
    int foundStart = 0;

    // -- Look for the start of the pulse
    for ( i = int(tStart); i < length - 1; i++)
    {
        if (foundStart) break;
        if ( polarity < 0)
        {
            if (data[i] < baseline - fabs( timesBase * sigmaBaseline))
            {
                pulseStartPoint = i;
                foundStart = 1;
            }
        }
        else if ( polarity > 0)
        {
            if (data[i] > baseline + fabs( timesBase * sigmaBaseline))
            {
                pulseStartPoint = i;
                foundStart = 1;
            }
        }
    }
    if (!foundStart) { /*printf("Start not found!\n"); */return -FLT_MAX; }// If the pulse is not higher than the noise, return -FLT_MAX

    // -- Look for the max/min
    double thresholdReal;
    double secondBaseline = 0;
    double minPulse = FLT_MAX;
    double maxPulse = -FLT_MAX;
    if (polarity < 0)
    {
        for (i=0; i < length; i++)
        {
            if (data[i] < minPulse) minPulse = data[i];
        }
        secondBaseline = double(minPulse);
    }
    else if (polarity > 0)
    {
        for (i=0; i < length; i++)
        {
            if (data[i] > maxPulse) maxPulse = data[i];
        }
        secondBaseline = double(maxPulse);
    }

    thresholdReal = baseline + (secondBaseline - baseline) * peakFraction;    // Threshold at some fraction of the pulse height
    double amplitude = fabs (secondBaseline - baseline);                     // Get the amplitude, might come in handy in future

    // Now take +/- samplesNum bins and do the linear fit
    double slope     = 0.;
    double intercept = 0.;
    double sumx      = 0.;
    double sumy      = 0.;
    double sumxy     = 0.;
    double sumx2     = 0.;
    int n = samplesNum * 2;

    if (polarity < 0)
    {
        for ( i = int(tStart); i < length-1 - samplesNum; i++)
        {
            if ((data[i] >= thresholdReal) && (data[i+1] < thresholdReal) )
            {
                for ( j = i-samplesNum; j < i+samplesNum; j++)
                {
                    sumx    += float(j);
                    sumy    += data[j];
                    sumxy   += j*data[j];
                    sumx2   += j*j;
                }
                slope       = (n*sumxy - sumx*sumy) / (n*sumx2 - sumx*sumx);
                intercept   = (sumy - slope * sumx) / double(n);
                crosspoint  = (thresholdReal - intercept) / slope;
                break;
            }
        }
    }
    else if (polarity > 0)
    {
        for ( i = int(tStart); i < length-1 - samplesNum; i++)
        {
            if ((data[i] <= thresholdReal) && (data[i+1] > thresholdReal) )
            {
                for ( j = i-samplesNum; j < i+samplesNum; j++)
                {
                    sumx    += float(j);
                    sumy    += data[j];
                    sumxy   += j*data[j];
                    sumx2   += j*j;
                }
                slope       = (n*sumxy - sumx*sumy) / (n*sumx2 - sumx*sumx);
                intercept   = (sumy - slope * sumx) / double(n);
                crosspoint  = (thresholdReal - intercept) / slope;
                break;
            }
        }

    }

// printf("----------------------------\n");                   // DEBUG
// printf("baseline      : %f\n", baseline);                   // DEBUG
// printf("secondBaseline: %f\n", secondBaseline);             // DEBUG
// printf("thresholdRead : %f\n", thresholdReal);              // DEBUG
// printf("slope         : %f\n", slope);                      // DEBUG
// printf("intercept     : %f\n", intercept);                  // DEBUG
// printf("crosspoint    : %f\n", crosspoint);                 // DEBUG
// printf("----------------------------\n");                   // DEBUG

return crosspoint;
}   // <-- Trigger CFD - END <--
















std::vector<std::string> Whitelist(std::vector<std::string> names, std::vector<std::string> blackList)
{
    // --- Not needed anymore, but so pretty I'll leave it.
    std::vector<std::string> whiteList;

    for ( auto && name : names) {
        if ( std::find( blackList.begin(), blackList.end(), name) == blackList.end()) {
            // std::cout << "Whitelisting: " << name << std::endl;  // DEBUG
            whiteList.push_back(name);
        }
    }

    return whiteList;
}











void GetGeometry (const std::string & inputFileName, std::vector<TH2I> & geoMaps)
{
    std::cout << "\n... Retrieving geometrical information from the input file...\n";
    TFile fileIn (inputFileName.c_str(), "READ");



    // --- CopyPasta from ROOT Tutorials "loopdir.C"
    //      Copy anything that inherits from TH2I
   TIter next(fileIn.GetListOfKeys());
   TKey *key;
   while ((key = (TKey*)next())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
        if (!cl->InheritsFrom("TH2I")) continue;
        TH2I *h = (TH2I*)key->ReadObj();
        // std::cout << h->GetName() << std::endl;
        geoMaps.emplace_back(TH2I (*(TH2I*) fileIn.Get(h->GetName())) );
        delete h;
   }

    delete key;
    fileIn.Close();


}
