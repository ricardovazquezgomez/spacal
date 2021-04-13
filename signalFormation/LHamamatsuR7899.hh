#ifndef LHAMAMATSUR7899_HH
#define LHAMAMATSUR7899_HH

#include <iostream>
#include <vector>

#include "LPhD.hh"

class LHamamatsuR7899 : public LPhD
{
 protected:
    const float _R = 50;                // Resistance
    const float _q = 1.602176634e-19;   // Elementary electric charge
    const float _tau_d = 1;             // [ns]

    float _G = 1e6;                     // Gain        
    // float _DefaultAmp = 5./6.;                 // Mode of the gamma distribution
    float _DefaultAmp;                  // Fixed amplitude, then multiblied by the fluctuation factor


 public:
    LHamamatsuR7899(const float & timeBin, const float & G=-1);
    LHamamatsuR7899(const float & timeBin, const float& SPTR, const float& QE_ScaleFactor, const float& G);           // Another constructor which will be used to change SPTR and QE and G from configfile
    ~LHamamatsuR7899(){};

    float const sPhelPulse(const float & t) const;

};

LHamamatsuR7899::LHamamatsuR7899(const float & timeBin, const float & G) : LPhD(timeBin)
{
    _detName = "Hamamatsu_R7899-20";

    _pulseSign  = -1.;    
    _SPTR       = 0.276;
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
    std::vector<double> wavelengths = {350.182,
                                        360.335,
                                        370.127,
                                        379.922,
                                        389.726,
                                        399.532,
                                        409.704,
                                        419.877,
                                        429.689,
                                        439.865,
                                        450.043,
                                        459.861,
                                        469.684,
                                        479.871,
                                        489.334,
                                        499.893,
                                        510.087,
                                        519.915,
                                        530.104,
                                        539.927,
                                        550.116,
                                        559.953,
                                        570.162,
                                        579.996,
                                        589.820,
                                        600.004,
                                        610.190,
                                        620.011,
                                        630.194,
                                        640.377,
                                        650.194,
                                        660.010,
                                        670.191,
                                        680.005,
                                        690.182,
                                        699.998};

    std::vector<double> QEs = {0.150,
                                    0.169,
                                    0.187,
                                    0.201,
                                    0.209,
                                    0.215,
                                    0.218,
                                    0.222,
                                    0.223,
                                    0.223,
                                    0.222,
                                    0.217,
                                    0.210,
                                    0.201,
                                    0.189,
                                    0.174,
                                    0.160,
                                    0.147,
                                    0.137,
                                    0.128,
                                    0.118,
                                    0.099,
                                    0.072,
                                    0.055,
                                    0.046,
                                    0.039,
                                    0.032,
                                    0.025,
                                    0.020,
                                    0.015,
                                    0.011,
                                    0.009,
                                    0.005,
                                    0.004,
                                    0.004,
                                    0.002};
    
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






LHamamatsuR7899::LHamamatsuR7899(const float & timeBin, const float & G, const float & SPTR, const float & QE_ScaleFactor) : LHamamatsuR7899 (timeBin, G)
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










float const LHamamatsuR7899::sPhelPulse(const float & t) const
{
    // ##########################################################
    // ### Single Photoelectron pulse!
    // ######
    return _DefaultAmp * pow( (t), 2) * exp(- (t) / _tau_d);
    
}




#endif