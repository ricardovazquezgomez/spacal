{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0000);


  TFile *_fileDepo = TFile::Open("energyResolutionDepo.root");
  TGraphErrors *grEnDep = (TGraphErrors*) _fileDepo->Get("Energy resolution (energy deposition alone) - Mono Pb+Poly");
  grEnDep->SetLineColor(kRed);
  grEnDep->SetMarkerColor(kRed);
  grEnDep->SetMarkerStyle(20);
  grEnDep->SetMarkerSize(1.5);

  TFile *_filePulses = TFile::Open("output.root");
  // TGraphErrors *grPulseNoCali   = (TGraphErrors*) _filePulses->Get("enResNoCali");
  // grPulseNoCali->SetLineColor(kGreen);
  // grPulseNoCali->SetMarkerColor(kGreen);
  // grPulseNoCali->SetMarkerStyle(21);
  // grPulseNoCali->SetMarkerSize(1.5);
  TGraphErrors *grPulseWithCali = (TGraphErrors*) _filePulses->Get("enResCali");
  grPulseWithCali->SetLineColor(kBlue);
  grPulseWithCali->SetMarkerColor(kBlue);
  grPulseWithCali->SetMarkerStyle(22);
  grPulseWithCali->SetMarkerSize(1.5);
  TGraphErrors *grPulseBack = (TGraphErrors*) _filePulses->Get("enResBack");
  grPulseBack->SetLineColor(kGreen);
  grPulseBack->SetMarkerColor(kGreen);
  grPulseBack->SetMarkerStyle(21);
  grPulseBack->SetMarkerSize(1.5);

  TGraphErrors *grPulseHalf = (TGraphErrors*) _filePulses->Get("enResHalf");
  grPulseHalf->SetLineColor(kRed);
  grPulseHalf->SetMarkerColor(kRed);
  grPulseHalf->SetMarkerStyle(20);
  grPulseHalf->SetMarkerSize(1.5);



  // gStyle->SetPalette(55);
  TMultiGraph *mg = new TMultiGraph();
  mg->Add(grEnDep,"p");
  mg->Add(grPulseBack,"p");
  mg->Add(grPulseWithCali,"p");

  mg->SetTitle("Mono Pb + Polystyrene 12x12 cm SPACAL - 3+3 deg electrons");
  mg->SetTitle("Energy Resolution - Single/Double Readout Comparison - 3#circ + 3#circ");

  mg->GetXaxis()->SetTitle("Energy [GeV]");
  mg->GetYaxis()->SetTitle("Energy Resolution [sigma E / E]");
  mg->GetYaxis()->SetTitle("#sigma_{E}/E");

  TCanvas *c1 = new TCanvas("c1","c1",1200,700);
  c1->SetLogx();
  c1->SetGrid();
  mg->Draw("aP");
  TLegend *legend = new TLegend(0.4,0.7,0.893,0.89,"");
  // legend->SetFillStyle(1);
  legend->AddEntry(grEnDep,"Light transport OFF - Only energy deposition","pl");
  legend->AddEntry(grPulseBack,"Light transport ON - Only back photodetectors","pl");
  legend->AddEntry(grPulseWithCali,"Light transport ON - Front and back photodetectors","pl");

  legend->Draw();
  gPad->Modified();
  mg->GetXaxis()->SetLimits(0.9,110);
  mg->SetMinimum(0.);
  mg->SetMaximum(0.12);
  c1->Draw();











  TMultiGraph *mg2 = new TMultiGraph();
  // mg2->Add(grEnDep,"p");
  // mg2->Add(grPulseBack,"p");
  mg2->Add(grPulseWithCali,"p");
  mg2->Add(grPulseHalf,"p");

  mg2->SetTitle("Mono Pb + Polystyrene 12x12 cm SPACAL - 3+3 deg electrons");
  mg2->GetXaxis()->SetTitle("Energy [GeV]");
  mg2->GetYaxis()->SetTitle("#frac{#sigma_E}{E}");
  TCanvas *c3 = new TCanvas("c3","c3",1000,800);
  c3->SetLogx();
  c3->SetGrid();
  mg2->Draw("aP PMC");
  TLegend *legend2 = new TLegend(0.5,0.7,0.893,0.89,"");
  // legend->SetFillStyle(1);
  // legend2->AddEntry(grEnDep,"Only Energy Deposition","p");
  // legend2->AddEntry(grPulseBack,"With optical photons - only back detector","p");
  legend2->AddEntry(grPulseWithCali,"Both detectors - All optical photons","p");
  legend2->AddEntry(grPulseHalf,"Both detectors - 50% optical photons","p");

  legend2->Draw();
  gPad->Modified();
  mg2->GetXaxis()->SetLimits(0.9,110);
  mg2->SetMinimum(0.);
  mg2->SetMaximum(0.20);


  TCanvas *c2 = new TCanvas("c2","c2",1000,800);
  c2->SetLogx();
  c2->SetGrid();
  grPulseWithCali->Draw("ap");
  grPulseWithCali->SetTitle("Mono Pb + Poly - With optical photons, with calibration");
  TF1* fEn = new TF1 ("fEn","TMath::Sqrt(TMath::Power([0]/TMath::Sqrt(x),2)+TMath::Power([1],2))",0,100);
  fEn->SetParameter(0,0.08);
  fEn->SetParameter(1,0.04);
  fEn->SetLineColor(kBlue);
  grPulseWithCali->Fit(fEn);


}
