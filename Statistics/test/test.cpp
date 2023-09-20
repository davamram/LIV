#include <TLorentzVector.h>
#include <TVector3.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TFormula.h>
#include <TF1.h>
#include <TSystem.h>
#include <TClonesArray.h>
#include <TLeaf.h>
#include <TChain.h>
#include <TObject.h>
#include <TLimit.h>
#include <TVectorD.h>
#include <TRandom3.h>
#include <TError.h>
#include <TMath.h>

//#include <iostream>
#include <string.h>
//#include <stdio.h>
//#include <stdlib.h>

Double_t ContinuousPoisson(Double_t *x, Double_t *params){
  Double_t A = params[0]; // Le paramètre A
  Double_t arg = x[0];
  return TMath::Exp(-A) * TMath::Power(A, arg) / TMath::Gamma(arg+1);
}

Double_t PoissonCdf(Double_t *x, Double_t *params){
  Double_t A = params[0];
  Double_t arg = x[0];
  return TMath::Gamma(arg+1, A)/TMath::Gamma(arg+1);
}

double GetContinuousPoissonCdf(double mean, double rand){
  TF1 *f = new TF1("f", PoissonCdf, floor(mean/5), 5*ceil(mean), 1); // plutot faire en fonction des sigma
  f->SetParameter(0, mean);

  cout<<"Rand : "<<rand<<" / min : "<<f->Eval(floor(mean/5))<<" / max : "<<f->Eval(ceil(mean)*5)<<endl;

  if(rand <= f->Eval(floor(mean/5))) return floor(mean/5);
  if(rand >= f->Eval(ceil(mean)*5)) return 5*ceil(mean);
  cout<<"X : "<<f->GetX(rand)<<endl;
  return f->GetX(rand);
}

double GetContinuousPoisson(double mean){
  TF1 *f = new TF1("f", ContinuousPoisson, mean/2, 2*ceil(mean), 1);
  f->SetParameter(0, mean);

  // TCanvas *canvas = new TCanvas("Nom_du_canvas", "Titre_du_canvas");
  // f->Draw();
  // TLine *line = new TLine(10, 0, 10, 10);
  // line->Draw("SAME");
  // canvas->Update();
  // TString nomFichier = Form("TF1_%.2f.png", mean);
  // canvas->SaveAs(nomFichier);

  double rand = f->GetRandom();
  delete f;
  //delete canvas;
  return rand;
}

Double_t ContinuousPoissonInt(Double_t *x, Double_t *params){
  Double_t A = params[0]; // Le paramètre A
  Double_t arg = x[0];
  arg= round(arg);
  return TMath::Exp(-A) * TMath::Power(A, arg) / TMath::Gamma(arg+1);
}

double GetContinuousPoissonInt(double mean){
  TF1 *f = new TF1("f", ContinuousPoissonInt, mean/2, 2*mean, 1);
  f->SetParameter(0, mean);
  double rand = f->GetRandom();
  delete f;
  return rand;
}

void test(){
  
  double mean=10;
  TH1F *histo1 = new TH1F("Distribution", "Distribution", 50, mean/2, 2*mean);
  TRandom3* rnd = new TRandom3();
  rnd->SetSeed(10111999);
  for(int i = 0; i<10000; i++){
    //histo1->Fill(GetContinuousPoisson(mean));
    histo1->Fill(GetContinuousPoissonCdf(mean, rnd->Uniform(0,1)));
  }
  TF1 *f1 = new TF1("f", PoissonCdf, mean/2, 2*mean, 1);
  f1->SetParameter(0, mean);

  TCanvas *canvas = new TCanvas("Nom_du_canvas", "Titre_du_canvas");
  //histo1->Scale(1/histo1->Integral());
  histo1->SetLineColor(kBlue);
  histo1->SetLineWidth(2);
  histo1->Draw("HIST");
  f1->SetLineColor(kGreen);
  f1->SetLineWidth(2);
  f1->Draw("SAME");
  cout<<f1->Integral(0, 3*mean)<<endl;
  cout<<histo1->Integral()<<endl;

  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(histo1, "Random Poisson", "l");
  legend->AddEntry(f1, "Poisson function", "l");
  legend->Draw("SAME");

  canvas->Update();
  canvas->SaveAs("Prout2.png");
}