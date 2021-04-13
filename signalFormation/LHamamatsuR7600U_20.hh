#ifndef LHAMAMATSUR7600U_20_HH
#define LHAMAMATSUR7600U_20_HH

// ##################################################################
// ### Class to simulate the Hamamatsu R7600U-20 (MCD-PMT) of 2020
// ### - SPTR (TTS) from Barcelona measurements
// ### - QE from datasheet
// ### - The rest is assumed identical to the Hamamatsu R7899-20, unfortunately info is missing.
// ###

#include <iostream>
#include <vector>

#include "LPhD.hh"

class LHamamatsuR7600U_20 : public LPhD
{
 protected:
    const float _R = 50;                // Resistance
    const float _q = 1.602176634e-19;   // Elementary electric charge
    const float _tau_d = 1;             // [ns]

    float _G = 1e6;                     // Gain        
    // float _DefaultAmp = 5./6.;                 // Mode of the gamma distribution
    float _DefaultAmp;                  // Fixed amplitude, then multiblied by the fluctuation factor


 public:
    LHamamatsuR7600U_20(const float & timeBin, const float & G=-1);
    LHamamatsuR7600U_20(const float & timeBin, const float& SPTR, const float& QE_ScaleFactor, const float& G);           // Another constructor which will be used to change SPTR and QE and G from configfile
    ~LHamamatsuR7600U_20(){};

    float const sPhelPulse(const float & t) const;

};

LHamamatsuR7600U_20::LHamamatsuR7600U_20(const float & timeBin, const float & G) : LPhD(timeBin)
{
    _detName = "Hamamatsu_R7600U-20";

    _pulseSign  = -1.;    
    _SPTR       = 0.030;            
    _pulseLen   = 10 * _tau_d;

    if (G >= 0) _G = G;
    _DefaultAmp        = _pulseSign * _R * _G * _q / pow(_tau_d,3) * 1e9 ;

    std::cout << "PMT Gain: " << _G << "\n";
    std::cout << "PMT DefaultAmp: " << _DefaultAmp << "\n";

    // _Amplitude  = _pulseSign * _DefaultAmp * _R * _G * _q / pow(_tau_d,3) * 1e9 ; // the last part is the conversion factor from [s] to [ns]

    // std::cout << "AMPLITUDE: " << this->_Amplitude << std::endl;
    // std::cout << "SPTR: " << this->_SPTR << std::endl;
    // this->PrintInfo();


    // -> Initialize QE ->
    std::vector<double> wavelengths = {
                                        283.9286,
                                        287.5000,
                                        289.2857,
                                        292.8571,
                                        294.6429,
                                        296.4286,
                                        298.2143,
                                        300.0000,
                                        301.7857,
                                        303.5714,
                                        305.3571,
                                        310.7143,
                                        316.0714,
                                        321.4286,
                                        335.7143,
                                        353.5714,
                                        378.5714,
                                        398.2143,
                                        416.0714,
                                        435.7143,
                                        455.3571,
                                        478.5714,
                                        505.3571,
                                        530.3571,
                                        555.3571,
                                        578.5714,
                                        601.7857,
                                        626.7857,
                                        651.7857,
                                        675.0000,
                                        700.0000,
                                        723.2143,
                                        748.2143,
                                        771.4286,
                                        796.4286,
                                        817.8571,
                                        835.7143,
                                        844.6429,
                                        851.7857,
                                        857.1429,
                                        862.5000,
                                        866.0714,
                                        867.8571,
                                        871.4286,
                                        875.0000,
                                        876.7857,
                                        880.3571,
                                        882.1429,
                                        885.7143,
                                        889.2857,
                                        891.0714,
                                        894.6429,
                                        898.2143,
                                        900.0000,
                                        903.5714,
                                        907.1429,
                                        908.9286,
                                        914.2857,
                                        917.8571
                                        };

    std::vector<double> QEs = {
                                    0.0014,
                                    0.0019,
                                    0.0025,
                                    0.0031,
                                    0.0042,
                                    0.0056,
                                    0.0075,
                                    0.0097,
                                    0.0128,
                                    0.0169,
                                    0.0223,
                                    0.0294,
                                    0.0388,
                                    0.0512,
                                    0.0637,
                                    0.0762,
                                    0.0858,
                                    0.1005,
                                    0.1225,
                                    0.1464,
                                    0.1682,
                                    0.1857,
                                    0.1894,
                                    0.1785,
                                    0.1649,
                                    0.1523,
                                    0.1380,
                                    0.1275,
                                    0.1178,
                                    0.1110,
                                    0.1025,
                                    0.0947,
                                    0.0875,
                                    0.0808,
                                    0.0732,
                                    0.0650,
                                    0.0523,
                                    0.0412,
                                    0.0306,
                                    0.0232,
                                    0.0176,
                                    0.0139,
                                    0.0103,
                                    0.0078,
                                    0.0059,
                                    0.0045,
                                    0.0034,
                                    0.0026,
                                    0.0019,
                                    0.0015,
                                    0.0011,
                                    0.0008,
                                    0.0006,
                                    0.0005,
                                    0.0004,
                                    0.0003,
                                    0.0002,
                                    0.0002,
                                    0.0001
                                    };


    if (_gQE) delete _gQE;
    _gQE = new TGraph(wavelengths.size(), wavelengths.data(), QEs.data() );
    _gQE -> GetXaxis() -> SetTitle("Wavelength [nm]");
    _gQE -> GetYaxis() -> SetTitle("Quantum Efficiency");

    Build_singlePhelPulse();
    if (_gSPR) delete _gSPR;
    _gSPR = new TGraph(_singlePhelPulse.size(), _times.data(), _singlePhelPulse.data() );
    _gSPR -> GetXaxis() -> SetTitle("Time [ns]");
    _gSPR -> GetYaxis() -> SetTitle("Amplitude [V]");





    if (_hSPTR) delete _hSPTR;
    TRandom3 random;
    _hSPTR = new TH1F ("_hSPTR", "_hSPTR", 100, -5*this->Get_SPTR(), +5*this->Get_SPTR());
    for (int i = 0; i < 1000; ++i) _hSPTR -> Fill ( random.Gaus(0, this->Get_SPTR()) );        




    // --- Distribution of single photoelectron fluctuations
    //    Beta distribution, based on Yuri Gouz' experience with this PMT.
    if (_fSPR_AmplitudeDistr) delete _fSPR_AmplitudeDistr;
    _fSPR_AmplitudeDistr = new TF1 ( "_fSPR", "ROOT::Math::gamma_pdf(x, [0], [1], 0.)", 0., 5.);
    _fSPR_AmplitudeDistr -> SetParameter(0, 6.);
    _fSPR_AmplitudeDistr -> SetParameter(1, 1./6.);

}






LHamamatsuR7600U_20::LHamamatsuR7600U_20(const float & timeBin, const float & G, const float & SPTR, const float & QE_ScaleFactor) : LHamamatsuR7600U_20 (timeBin, G)
{
    // --- Different constructor which takes in input an arbitrary SPTR value and a scale factor and
    //      changes the SPTR and rescales the QE. And the GAIN, too!!!!
    std::cout << "\n !!! WARNING: custom PhD constructor with user-defined SPTR = " << SPTR << " and scaling QE by " << QE_ScaleFactor << "\n";

    if (SPTR >= 0) {
        _SPTR   = SPTR;
        _usr_def_parameters = true;

        if (_hSPTR) delete _hSPTR;
        TRandom3 random;
        _hSPTR = new TH1F ("_hSPTR", "_hSPTR", 100, -5*this->Get_SPTR(), +5*this->Get_SPTR());
        for (int i = 0; i < 10000; ++i) _hSPTR -> Fill ( random.Gaus(0, this->Get_SPTR()) );        

    }
    else std::cout << "SPTR < 0. Not changed.\n";

    if (QE_ScaleFactor >= 0) {
        for (int i=0; i < _gQE->GetN(); i++) _gQE->GetY()[i] *= QE_ScaleFactor;
        _usr_def_parameters = true;
    }
    else std::cout << "QE_ScaleFactor < 0. Not changed.\n";

}










float const LHamamatsuR7600U_20::sPhelPulse(const float & t) const
{
    // ##########################################################
    // ### Single Photoelectron pulse!
    // ######
    return _DefaultAmp * pow( (t), 2) * exp(- (t) / _tau_d);
    
}




#endif