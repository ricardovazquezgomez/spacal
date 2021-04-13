#ifndef LHAMAMATSUH13543_HH
#define LHAMAMATSUH13543_HH

#include <iostream>
#include <vector>

#include "LPhD.hh"

class LHamamatsuH13543 : public LPhD
{
 protected:
 // --- Parameters assumed identical to LHamamatsuR7899
    const float _R = 50;                // Resistance
    const float _q = 1.602176634e-19;   // Elementary electric charge
    const float _tau_d = 7;             // [ns]

    float _G = 1e6;                     // Gain
    float _DefaultAmp;                  // Fixed amplitude, then multiplied by the fluctuating factor


 public:
    LHamamatsuH13543(const float & timeBin, const float & G = -1);
    LHamamatsuH13543(const float & timeBin, const float & G, const float & SPTR, const float & QE_ScaleFactor);           // Another constructor which will be used to change SPTR and QE from configfile
    ~LHamamatsuH13543(){};

    float const sPhelPulse(const float & t) const;

};

LHamamatsuH13543::LHamamatsuH13543(const float & timeBin, const float & G) : LPhD(timeBin)
{
    _detName = "Hamamatsu_H13543_ma-PMT";

    _pulseSign  = -1.;    
    _SPTR       = 0.1486;
    _pulseLen   = 10 * _tau_d;

    if ( G >= 0.) _G = G;

    _DefaultAmp        = _pulseSign * _R * _G * _q / pow(_tau_d,3) * 1e9 ;
    
    std::cout << "PMT Gain: " << _G << "\n";
    std::cout << "PMT DefaultAmp: " << _DefaultAmp << "\n";

    // _Amplitude  = _pulseSign * _DefaultAmp * _R * _G * _q / pow(_tau_d,3) * 1e9 ; // the last part is the conversion factor from [s] to [ns]

    // std::cout << "AMPLITUDE: " << this->_Amplitude << std::endl;
    // std::cout << "SPTR: " << this->_SPTR << std::endl;
    // std::cout << "PulseSign: " << _pulseSign << "\n";
    // std::cout << "SPTR: " << _SPTR << "\n";
    // std::cout << "pulseLen: " << _pulseLen << "\n";
    // std::cout << "Amplitude: " << _Amplitude << "\n";


    // this->PrintInfo();

    // -> Initialize QE ->
    std::vector<double> wavelengths = {269.9360341151,
273.3475479744,
280.1705756930,
288.6993603412,
295.5223880597,
298.9339019190,
304.0511727079,
312.5799573561,
321.1087420043,
327.9317697228,
338.1663113006,
350.1066098081,
362.0469083156,
375.6929637527,
387.6332622601,
399.5735607676,
409.8081023454,
418.3368869936,
428.5714285714,
440.5117270789,
449.0405117271,
457.5692963753,
464.3923240938,
471.2153518124,
479.7441364606,
488.2729211087,
498.5074626866,
508.7420042644,
515.5650319829,
520.6823027719,
527.5053304904,
532.6226012793,
537.7398720682,
542.8571428571,
549.6801705757,
556.5031982942,
563.3262260128,
570.1492537313,
575.2665245203,
582.0895522388,
592.3240938166,
600.8528784648,
609.3816631130,
617.9104477612,
624.7334754797,
631.5565031983,
638.3795309168,
645.2025586354,
650.3198294243,
655.4371002132,
660.5543710021,
663.9658848614,
669.0831556503,
672.4946695096,
677.6119402985,
682.7292110874,
686.1407249467,
689.5522388060,
694.6695095949,
698.0810234542,
701.4925373134,
704.9040511727,
708.3155650320,
711.7270788913,
715.1385927505,
720.2558635394,
723.6673773987,
727.0788912580,
730.4904051173,
735.6076759062,
739.0191897655,
742.4307036247,
745.8422174840,
750.9594882729,
756.0767590618,
759.4882729211,
762.8997867804,
768.0170575693,
771.4285714286,
776.5458422175,
779.9573560768,
785.0746268657,
790.1918976546};
        
        
        
        
        
        
        
      
    std::vector<double> QEs = {.15458442409,
.18430812176,
.22123482942,
.24962992683,
.26795916348,
.29107195349,
.31444549612,
.33608371124,
.35276840934,
.35869468002,
.36806795333,
.37619061283,
.37775483658,
.38520483823,
.38767729884,
.38324629873,
.37367511125,
.37301977761,
.36411184263,
.34762999217,
.33466159985,
.32842373220,
.31163036040,
.29575953302,
.27975760104,
.26974016867,
.25443097869,
.24009186412,
.22815236471,
.20950885685,
.19178521649,
.17286451755,
.15878884354,
.14587250612,
.13360284514,
.12009960580,
.11003125127,
.10082175161,
.094432998029,
.086550553883,
.075962704833,
.068151456630,
.058893960978,
.052858894654,
.045821119545,
.040480738415,
.035099306940,
.030436749596,
.027482157602,
.024352691434,
.021178064744,
.018816950504,
.016365593782,
.014541967911,
.012648756703,
.010797269840,
.009595097718,
.008367819218,
.007279777083,
.006349029427,
.005642743222,
.004921530660,
.004292598523,
.003744125791,
.003145023769,
.002635588429,
.002299017372,
.001968037072,
.001653295637,
.001359833597,
.001164152412,
.000978047298,
.000821711029,
.000675953033,
.000545695985,
.000449948172,
.000378064874,
.000311049836,
.000256490946,
.000203236058,
.000164467807,
.000130328986,
.000103280882};
    
    _gQE = new TGraph(wavelengths.size(), wavelengths.data(), QEs.data() );
    _gQE -> GetXaxis() -> SetTitle("Wavelength [nm]");
    _gQE -> GetYaxis() -> SetTitle("Quantum Efficiency");


    Build_singlePhelPulse();
    _gSPR = new TGraph(_singlePhelPulse.size(), _times.data(), _singlePhelPulse.data() );
    _gSPR -> GetXaxis() -> SetTitle("Time [ns]");
    _gSPR -> GetYaxis() -> SetTitle("Amplitude [V]");


    if (_hSPTR) delete _hSPTR;
    TRandom3 random;
    _hSPTR = new TH1F ("_hSPTR", "_hSPTR", 100, -5*this->Get_SPTR(), +5*this->Get_SPTR());
    for (int i = 0; i < 10000; ++i) _hSPTR -> Fill ( random.Gaus(0, this->Get_SPTR()) );        


}



LHamamatsuH13543::LHamamatsuH13543(const float & timeBin, const float & G, const float & SPTR, const float & QE_ScaleFactor) : LHamamatsuH13543 (timeBin, G)
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





float const LHamamatsuH13543::sPhelPulse(const float & t) const
{
    // ##########################################################
    // ### Single Photoelectron pulse!
    // ######
    // std::cout << _Amplitude * pow( (t), 2) * exp(- (t) / _tau_d) << std::endl;
    return _DefaultAmp * pow( (t), 2) * exp(- (t) / _tau_d);
    
}




#endif