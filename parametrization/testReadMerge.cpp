// Marco Pizzichemi 11.03.2020 marco.pizzichemi@cern.ch

// compile with
// g++ -o ../build/testReadMerge testReadMerge.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore


// program to merge optical calibration files


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


#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <getopt.h>

#include <sys/stat.h>
#include <dirent.h>

#include "regions.h"

std::ifstream::pos_type filesize(const char* filename)
{
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}


struct calibration_point_t
{
  float x;
  float y;
  float z;
  int   s;
  float energy;
  int primaries;
  float eff;
  TH1F *delay;
};

struct calibration_point_output_t
{
  int ix;
  int iy;
  int iz;
  int is;
  int ie;
  float x;
  float y;
  float z;
  int   s;
  float e;
  float eff;
  float *binCenters;
  float *binContents;
};




// struct header_t
// {
//   int bins;
//   int xn;
//   int yn;
//   int zn;
//   int en;
//   std::vector<float> x;
//   std::vector<float> y;
//   std::vector<float> z;
//   std::vector<float> energy;
// };

int main(int argc, char** argv)
{


  // IMPORTED
  std::string inputFile = argv[1];
  std::ifstream input_file(inputFile.c_str(), std::ios::binary);
  // read header
  int bins;
  float pdfStart,pdfEnd;
  int xn;
  int yn;
  int zn;
  int en;
  float x_v,y_v,z_v,e_v;
  std::vector <float > x,y,z,energy;

  // store ints
  input_file.read( (char*)&bins , sizeof(int));
  input_file.read( (char*)&pdfStart , sizeof(float));
  input_file.read( (char*)&pdfEnd , sizeof(float));
  input_file.read( (char*)&xn , sizeof(int));
  input_file.read( (char*)&yn , sizeof(int));
  input_file.read( (char*)&zn , sizeof(int));
  input_file.read( (char*)&en , sizeof(int));
  std::cout << bins << " "
            << pdfStart << " "
            << pdfEnd << " "
            << xn << " "
            << yn << " "
            << zn << " "
            << en << " "
            << std::endl;
  // read vectors
  for(int i = 0 ; i < xn; i++)
  {
    input_file.read( (char*)&x_v , sizeof(float));
    x.push_back(x_v);
  }
  for(int i = 0 ; i < yn; i++)
  {
    input_file.read( (char*)&y_v , sizeof(float));
    y.push_back(y_v);
  }
  for(int i = 0 ; i < zn; i++)
  {
    input_file.read( (char*)&z_v , sizeof(float));
    z.push_back(z_v);
  }
  for(int i = 0 ; i < en; i++)
  {
    input_file.read( (char*)&e_v , sizeof(float));
    energy.push_back(e_v);
  }

  float xmin  = 0;
  float xstep = 1;
  float ymin  = 0;
  float ystep = 1;
  float zmin  = 0;
  float zstep = 1;
  float emin  = 1;
  float estep = 1;
  int   nSides  = 2;

  xmin = x[0];
  ymin = y[0];
  zmin = z[0];
  emin = energy[0];
  if(x.size() == 1)
  {
    xstep = 1;
  }
  else
  {
    xstep = x[1] - x[0];
  }
  if(y.size() == 1)
  {
    ystep = 1;
  }
  else
  {
    ystep = y[1] - y[0];
  }
  if(z.size() == 1)
  {
    zstep = 1;
  }
  else
  {
    zstep = z[1] - z[0];
  }
  if(energy.size() == 1)
  {
    estep = 1;
  }
  else
  {
    estep = energy[1] - energy[0];
  }


  // FEEDBACK
  std::cout << std::endl;
  std::cout << "----------------------------- " << std::endl;
  std::cout << "Calibration points " << std::endl;
  std::cout << "----------------------------- " << std::endl;
  std::cout << std::endl;
  std::cout << "---> X: " << std::endl;
  std::cout << "Number of points in X           = " << xn << std::endl;
  std::cout << "Step length           [mm]      = " << xstep << std::endl;
  std::cout << "Starting X coordinate [mm]      = " << xmin  << std::endl;
  std::cout << "List of X points      [mm]      = " ;
  for(unsigned int ix = 0; ix < x.size(); ix++)
  {
    std::cout << x[ix] << " ";
  }
  std::cout << std::endl;

  std::cout << "---> Y: " << std::endl;
  std::cout << "Number of points in Y           = " << yn << std::endl;
  std::cout << "Step length           [mm]      = " << ystep << std::endl;
  std::cout << "Starting Y coordinate [mm]      = " << ymin  << std::endl;
  std::cout << "List of Y points      [mm]      = " ;
  for(unsigned int iy = 0; iy < y.size(); iy++)
  {
    std::cout << y[iy] << " ";
  }
  std::cout << std::endl;

  std::cout << "---> Z: " << std::endl;
  std::cout << "Number of points in Z           = " << zn << std::endl;
  std::cout << "Step length           [mm]      = " << zstep << std::endl;
  std::cout << "Starting Z coordinate [mm]      = " << zmin  << std::endl;
  std::cout << "List of Z points      [mm]      = " ;
  for(unsigned int iz = 0; iz < z.size(); iz++)
  {
    std::cout << z[iz] << " ";
  }
  std::cout << std::endl;
  std::cout << "---> ENERGY: " << std::endl;
  std::cout << "Number of points in ENERGY      = " << en << std::endl;
  std::cout << "Step length                [eV] = " << estep << std::endl;
  std::cout << "Starting ENERGY coordinate [eV] = " << emin  << std::endl;
  std::cout << "List of ENERGY points      [eV] = " ;
  for(unsigned int ie = 0; ie < energy.size(); ie++)
  {
    std::cout << energy[ie] << " ";
  }
  std::cout << std::endl;
  std::cout << "----------------------------- " << std::endl;
  std::cout << std::endl;
  // ------------


  // long long int file0N = filesize(inputFile.c_str()) /  sizeof(calibration_point_output_t);
  // std::cout << "Calibration points in file " << inputFile.c_str() << " = " << file0N << std::endl;
  // input_file.read((char*)&master, sizeof(master));


  // prepare the array

  // prepare calibration structs
  calibration_point_t *****calibration_point;
  calibration_point = new calibration_point_t****[x.size()];
  for(int ix = 0 ; ix < x.size() ; ix++)
  {
    calibration_point[ix] = new calibration_point_t***[y.size()];
    for(int iy = 0 ; iy < y.size() ; iy++)
    {
      calibration_point[ix][iy] = new calibration_point_t**[z.size()];
      for(int iz = 0 ; iz < z.size() ; iz++)
      {
        calibration_point[ix][iy][iz] = new calibration_point_t*[energy.size()];
        for(int ie = 0 ; ie < energy.size() ; ie++)
        {
          calibration_point[ix][iy][iz][ie] = new calibration_point_t[nSides];
          for(int iSide = 0; iSide < nSides ; iSide++)
          {
            std::stringstream ssname;
            ssname << "delayHisto " << x[ix] << " " << y[iy] << " " << z[iz] << " " << energy[ie] << " " << iSide;
            calibration_point[ix][iy][iz][ie][iSide].x = x[ix];
            calibration_point[ix][iy][iz][ie][iSide].y = y[iy];
            calibration_point[ix][iy][iz][ie][iSide].z = z[iz];
            calibration_point[ix][iy][iz][ie][iSide].energy = energy[ie];
            calibration_point[ix][iy][iz][ie][iSide].eff = 0;
            calibration_point[ix][iy][iz][ie][iSide].delay = new TH1F(ssname.str().c_str(),ssname.str().c_str(),bins,pdfStart,pdfEnd);
          }
        }
      }
    }
  }


  calibration_point_output_t calibration_point_out;
  calibration_point_out.binCenters = new float [bins];
  calibration_point_out.binContents = new float [bins];
  //
  // // input_file.read((char*)&calibration_point, sizeof(calibration_point_output_t));
  // // std::cout << calibration_point.x << std::endl;
  // while( input_file.read((char*)&calibration_point_out, sizeof(calibration_point_out)) )
  int CaliPoints = xn*yn*zn*en*nSides;
  std::cout << CaliPoints << std::endl;
  for(int iCali = 0; iCali < CaliPoints; iCali++)
  {
    int ix   ;
    int iy   ;
    int iz   ;
    int ie   ;
    int iSide;
    float x_i  ;
    float y_i  ;
    float z_i  ;
    int e_i    ;
    int s_i    ;
    float eff;
    float value ;

    input_file.read( (char*)&ix , sizeof(int));
    input_file.read( (char*)&iy , sizeof(int));
    input_file.read( (char*)&iz , sizeof(int));
    input_file.read( (char*)&ie , sizeof(int));
    input_file.read( (char*)&iSide , sizeof(int));

    input_file.read( (char*)&x_i , sizeof(float));
    input_file.read( (char*)&y_i , sizeof(float));
    input_file.read( (char*)&z_i , sizeof(float));
    input_file.read( (char*)&s_i , sizeof(int));
    input_file.read( (char*)&e_i , sizeof(float));
    input_file.read( (char*)&eff , sizeof(float));

    // std:: cout << ix << " " << " " << iy << " " << iz << " " << ie << " " << iSide << std::endl;

    calibration_point[ix][iy][iz][ie][iSide].x = x_i;
    calibration_point[ix][iy][iz][ie][iSide].y = y_i;
    calibration_point[ix][iy][iz][ie][iSide].z = z_i;
    calibration_point[ix][iy][iz][ie][iSide].s = s_i;
    calibration_point[ix][iy][iz][ie][iSide].energy = e_i;
    calibration_point[ix][iy][iz][ie][iSide].primaries = 0; // not needed now
    calibration_point[ix][iy][iz][ie][iSide].eff = eff;

    for(int b = 1 ; b < bins+1; b++)
    {
      input_file.read( (char*)&value , sizeof(float));
      calibration_point[ix][iy][iz][ie][iSide].delay->SetBinContent(b,value);
    }

     // std::cout << calibration_point[ix][iy][iz][ie][iSide].x << std::endl;
    // std::cout << calibration_point_out.binContents[2] << std::endl;
    // for(int i = 0 ; i < bins ; i++)
    // {
      // calibration_point[ix][iy][iz][ie][iSide].delay->SetBinContent(i+1,calibration_point_out.binContents[i]);
    // }
    // calibration_point[ix][iy][iz][ie][iSide].delay;

  }
  // IMPORTED -- end

  TFile *outputFile = new TFile("histograms.root","RECREATE");
  outputFile->cd();

  TDirectory *histoDir = outputFile->mkdir("Hybrid_Calibration");
  histoDir->cd();
  for(int ix = 0 ; ix < x.size() ; ix++)
  {
    for(int iy = 0 ; iy < y.size() ; iy++)
    {
      for(int iz = 0 ; iz < z.size() ; iz++)
      {
        for(int ie = 0 ; ie < energy.size() ; ie++)
        {
          for(int iSide = 0; iSide < nSides ; iSide++)
          {
            calibration_point[ix][iy][iz][ie][iSide].delay->Write();
          }
        }
      }
    }
  }
  outputFile->Close();
  std::cout << "Done. Goodbye!" << std::endl;




  return 0;
}
