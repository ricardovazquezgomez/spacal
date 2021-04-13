#ifndef LPHD_HH
#define LPHD_HH

#include <vector>
#include <iostream>
#include <string>

// --- To Draw SPR and QE
#include "TGraph.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1F.h"
#include "TRandom3.h"

class LPhD
{
 protected:

  bool _usr_def_parameters = false;

  float _Amplitude;
  float _pulseSign;
  float _pulseLen;
  float _SPTR;


  float _timeBin;          


  std::string _detName;

  std::vector<float> _singlePhelPulse;
  std::vector<float> _singlePhelSlopes;
  std::vector<float> _times;  

  TGraph* _gQE;
  TGraph* _gSPR;                                      
  TF1*  _fSPR_AmplitudeDistr;                           // TF1 function representing the distribution of single photoelectron pulses' amplitude. Used to simulate fluctuations in it
  TH1F* _hSPTR;

 public:
  LPhD( const float & timeBin );
  ~LPhD(){ delete _gQE; delete _gSPR; delete _hSPTR; delete _fSPR_AmplitudeDistr; };

  // --- Setters
  inline void SetRand_Amplitude () { _Amplitude = _fSPR_AmplitudeDistr -> GetRandom(); std::cout << "_Amplitude changed: " << _Amplitude << std::endl;};

  // --- Getters
  inline float const & Get_SPTR() const { return _SPTR; };
  inline float const & Get_pulseLen() const { return _pulseLen; };
  inline float const Get_QE (const float & wavelength) const { return _gQE->Eval(wavelength); };
  inline float const Get_Amplitude () const { return _Amplitude; };
  inline float const GetRand_Amplitude () const { return _fSPR_AmplitudeDistr -> GetRandom(); };

  inline TGraph* Get_QEGraph () const { return _gQE; };
  inline TGraph* Get_SPRGraph () const { return _gSPR; };
  inline TF1* Get_fSPR_AmplitudeDistr () const { return _fSPR_AmplitudeDistr; }
  inline TH1F* Get_hSPTR() const { return _hSPTR; };

  inline std::vector<float> const Get_singlePhelPulse () const { return _singlePhelPulse; };

  // --- Members
  virtual void const PrintInfo() const;

  void const Build_singlePhelPulse();

  inline float const sPhelPulseInterp(const float & t) const { int i = t/_timeBin; return _singlePhelPulse[i] + (t - _times[i]) * (_singlePhelPulse[i+1] - _singlePhelPulse[i]) / _timeBin ;};
  virtual float const sPhelPulse(const float & t) const { return _Amplitude;};
};

LPhD::LPhD (const float & timeBin) : _Amplitude (1.),
                                      _pulseSign (-1.),
                                      _pulseLen (0.),
                                      _SPTR (0.),
                                      _timeBin(timeBin),
                                      _detName ("Base")
{
  std::vector<double> wavelengths = {1,1200};
  std::vector<double> QEs = {1,1};
  _gQE = new TGraph(wavelengths.size(), wavelengths.data(), QEs.data() );
  _gQE -> GetXaxis() -> SetTitle("Wavelength [nm]");
  _gQE -> GetYaxis() -> SetTitle("Quantum Efficiency");

  Build_singlePhelPulse();
  _gSPR = new TGraph(_singlePhelPulse.size(), _times.data(), _singlePhelPulse.data() );
  _gSPR -> GetXaxis() -> SetTitle("Time [ns]");
  _gSPR -> GetYaxis() -> SetTitle("Amplitude [V]");



  // if (_hSPTR) delete _hSPTR;     // Provokes crashes on lxplus sometimes, not always. Why? What have I misunderstood?
  TRandom3 random;
  _hSPTR = new TH1F ("_hSPTR", "_hSPTR", 100, -5*this->Get_SPTR(), +5*this->Get_SPTR());
  for (int i = 0; i < 10000; ++i) _hSPTR -> Fill ( random.Gaus(0, this->Get_SPTR()) );        


  



  // --- Distribution of single photoelectron fluctuations
  //    For the LPhD will be a very narrow normalized gaussian centered around 1.
  _fSPR_AmplitudeDistr = new TF1 ( "_fSPR", "gausn");
  _fSPR_AmplitudeDistr -> SetParameter(0, 1.);
  _fSPR_AmplitudeDistr -> SetParameter(1, 1.);
  _fSPR_AmplitudeDistr -> SetParameter(2, 0.001);
}


void const LPhD::PrintInfo() const
{
  std::cout << "\n####################\n";
  std::cout << "### Detector name  : " << _detName << "\n";
  std::cout << "### SPTR [ns]      : " << _SPTR << "\n";
  std::cout << "### sPulse len [ns]: " << _pulseLen << "\n";
  std::cout << "### sPulse sign[ns]: " << _pulseSign << "\n";
  std::cout << "### Amplitude      : " << _Amplitude << "\n";
  if (_usr_def_parameters) std::cout << "### WARNING: Employed user-defined parameters!\n";
  std::cout << "####################\n\n";
}

void const LPhD::Build_singlePhelPulse() 
{
  // ###################################
  // ### Here I prepare the sampled single photoelectron pulse and a linear interpolation to speed up the program
  // ######

  // --- Make sure the vectors are empty
  _singlePhelPulse.clear();       
  _singlePhelSlopes.clear();
  _times.clear();

  // --- Sampled pulse
  for (int i=0; i <= int(_pulseLen/_timeBin); i++)  _singlePhelPulse.emplace_back( sPhelPulse(i*_timeBin) );                        
  _singlePhelPulse.emplace_back(0);

  // --- Sampled times                   
  for (int i=0; i <= int(_pulseLen/_timeBin); i++)  _times.emplace_back( i*_timeBin );                                               
  _times.emplace_back( (int(_pulseLen/_timeBin)+1) * _timeBin );

  // --- Slopes for the interpolation
  for (int i=0; i < _singlePhelPulse.size()-1; i++) _singlePhelSlopes.emplace_back( (_singlePhelPulse.at(i+1) - _singlePhelPulse.at(i)) / _timeBin );      
  _singlePhelSlopes.emplace_back( (0 - _singlePhelPulse.at(_singlePhelPulse.size()-1)) / _timeBin );
}


#endif