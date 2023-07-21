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



double* GetPoissonDistribution(double mean, int ntrials){

  TRandom3* rnd = new TRandom3();
  rnd->SetSeed(10111999);
  double* P = new double[ntrials];

  for (int i=0; i<ntrials; i++){
    P[i] = rnd->PoissonD(mean);
  }

  return P;
}

double* GetDataPoissonDistribution_Fluctuate(double N, double N_err, int ntrials){

  double mean;

  double N_fluctuated;

  TRandom3* rnd = new TRandom3();
  double* P = new double[ntrials];

  for (int i=0; i<ntrials; i++){

    N_fluctuated=0;
    while (N_fluctuated<=0){
      N_fluctuated = rnd->Gaus(N, N_err);
    }

    P[i] = rnd->PoissonD(N_fluctuated);

  }

  return P;
}

double CalcLikelihood(double n, double x){

  double l = exp(-x)*TMath::Power(x,n)/TMath::Gamma(n+1);
  return l;
}

double CalcLogLikelihoodRatio(double Ndata, double Ns, double Nb){

  double lr;
  if (Nb>0) lr = CalcLikelihood(Ndata, Ns)/CalcLikelihood(Ndata, Nb);

  return log(lr);
}

double GetSignificanceStandardDeviation(double alpha){

  //cout << "alpha="<<alpha<<endl;

  double SdDeviation;
  double alpha_test = 1;
  double signif_test = 0;

  while (alpha_test>1-alpha){
    signif_test += 0.01;
    alpha_test = 1 - ROOT::Math::erf(signif_test/sqrt(2));
    //cout << "signif_test="<<signif_test<<" alpha_test="<<alpha_test<<endl;
  }

  SdDeviation = signif_test;

  return SdDeviation;
}

TH1F* PlotGedankenExpMinus2LLR_MultiChannels(string Name, double** GedankenExp, double* Ns, double* Nb, int ntrials, int nChannels){

  TH1F* GedankenExpMinus2LLR = new TH1F(Name.c_str(), Name.c_str(), 100000, -100, 100);

  double LLR=0;
  for (int i=0; i<ntrials; i++){
    LLR=0;
    for (int j=0; j<nChannels; j++){
      LLR += CalcLogLikelihoodRatio(GedankenExp[j][i], Ns[j], Nb[j]);
    }
    GedankenExpMinus2LLR->Fill(-2*LLR);
  }
  double a = GedankenExpMinus2LLR->Integral();
  GedankenExpMinus2LLR->Scale(1/a);

  return GedankenExpMinus2LLR;

}

double ComputeCL_MultiChannels_hypSB(TH1F* HistoGedankenExp, double* Ns, double* Nb, double* Ndata, int nChannels, double &LLR_data, string name){

  LLR_data = 0;
  double LLR_formula = 0;
  for (int i=0; i<nChannels; i++){
  
    LLR_data += CalcLogLikelihoodRatio(Ndata[i], Ns[i], Nb[i]);
    
  }

  double xmin = HistoGedankenExp->GetXaxis()->GetXmin();
  double xmax = HistoGedankenExp->GetXaxis()->GetXmax();
  double binwidth = (xmax-xmin)/((double)HistoGedankenExp->GetNbinsX());

  int bin1 = (-2*LLR_data - xmin)/binwidth;
  int bin2 = (xmax - xmin)/binwidth;

  double CL = HistoGedankenExp->Integral(bin1, bin2);

  return CL;
}

double ComputeCL_MultiChannels_hypB(TH1F* HistoGedankenExp, double* Ns, double* Nb, double* Ndata, int nChannels, double &LLR_data, string name){

  LLR_data = 0;
  double LLR_formula = 0;
  for (int i=0; i<nChannels; i++){
  
    LLR_data += CalcLogLikelihoodRatio(Ndata[i], Ns[i], Nb[i]);
    
  }

  double xmin = HistoGedankenExp->GetXaxis()->GetXmin();
  double xmax = HistoGedankenExp->GetXaxis()->GetXmax();
  double binwidth = (xmax-xmin)/((double)HistoGedankenExp->GetNbinsX());

  int bin1 = (-2*LLR_data - xmin)/binwidth;

  double CL = HistoGedankenExp->Integral(1, bin1);

  return CL;
}

void ColorBinsUp(TH1F* h1, TH1F*& h2, Double_t LLR, int color) {

    h2 = (TH1F*)h1->Clone("h2");

    // Color only what I want
    int binLLR = h2->GetXaxis()->FindBin(LLR);
    h2->GetXaxis()->SetRangeUser(LLR, h2->GetXaxis()->GetXmax());
    h2->SetFillColor(kRed);
    h2->SetLineColor(kRed);
}
void ColorBinsDown(TH1F* h1, TH1F*& h2, Double_t LLR, int color) {

    h2 = (TH1F*)h1->Clone("h2");

    // Color only what I want
    h2->GetXaxis()->SetRangeUser(h2->GetXaxis()->GetXmin(), LLR);
    h2->SetFillColor(kGreen);
    h2->SetLineColor(kGreen);
}

void SetXrange(TH1F*& h){
  int min=1;
  int max=h->GetNbinsX();
  double maxBin = h->GetMaximum();
  for(int i=0; i<h->GetNbinsX();i++){
    if(h->GetBinContent(i)>0.0001*maxBin){
      min=i;
      cout<<"min : "<<min<<endl;
      break;
    }
  }
  for(int i=h->GetNbinsX(); i>0;i--){
    if(h->GetBinContent(i)>0.001*maxBin){
      max=i;
      cout<<"max : "<<max<<endl;
      break;
    }
  }
  h->GetXaxis()->SetRange(min,max);
}

void DrawHist(TH1F* hH0, TH1F* hH1, double LLR_data, string name){
  // hH0->Rebin(100000/100);
  // hH1->Rebin(100000/100);
  hH0->SetTitle("Test Statistic");
  TCanvas* c1 = new TCanvas("canvas", "Canvas", 800, 600);
  double maxHist = TMath::Max(hH0->GetMaximum(), hH1->GetMaximum());
  hH0->SetMaximum(1.3 * maxHist);
  // I only do it for this one and not for the others. I could do a more sophisticated thing to take into account others but it is not very usefull
  SetXrange(hH0);
  hH0->SetLineColor(kBlue);
  hH0->SetLineWidth(2);
  hH1->SetLineColor(kViolet);
  hH1->SetLineWidth(2);
  hH0->Draw("HIST C");
  hH1->Draw("HIST C SAME");

  //color for LLR

  TH1F* hH1_color;
  TH1F* hH0_color;
  ColorBinsUp(hH1, hH1_color, LLR_data, kRed);
  ColorBinsDown(hH0, hH0_color, LLR_data, kGreen);
  hH1_color->Draw("SAME C HIST");
  hH0_color->Draw("SAME C HIST");

  TLine* line1 = new TLine(LLR_data, 0, LLR_data, maxHist);
  line1->SetLineColor(kBlack);
  line1->SetLineWidth(2);
  line1->Draw("SAME");

  TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Coordonnées de la légende dans le canvas (x1, y1, x2, y2)
  legend->AddEntry(hH0, "Expected for SM", "l"); // "l" signifie ligne, vous pouvez également utiliser "f" pour le remplissage
  legend->AddEntry(hH1, "Expected for LIV", "l");
  legend->AddEntry(hH0_color, "P(SM)", "f");
  legend->AddEntry(hH1_color, "P(LIV)", "f");
  legend->AddEntry(line1, "Observed (in SM-MC sample)", "l");
  legend->Draw();



  c1->SaveAs(name.c_str());
}

double* GenerateToyExperiment_MultiChannels(double* Ns, double* Nb, double* Ndata, int nChannels, int ntrials){

  double* result = new double[3];

  double** GedankenExp_SBhyp = new double*[nChannels];
  double** GedankenExp_Bhyp = new double*[nChannels];

  for (int i=0; i<nChannels; i++){
    GedankenExp_SBhyp[i] = GetPoissonDistribution(Ns[i], ntrials);
    GedankenExp_Bhyp[i] = GetPoissonDistribution(Nb[i], ntrials);
  }

  TH1F* HistoGedankenExp_SBhyp = PlotGedankenExpMinus2LLR_MultiChannels("GedankenExp_SBhyp", GedankenExp_SBhyp, Ns, Nb, ntrials, nChannels);
  TH1F* HistoGedankenExp_Bhyp = PlotGedankenExpMinus2LLR_MultiChannels("GedankenExp_Bhyp", GedankenExp_Bhyp, Ns, Nb, ntrials, nChannels);

  cout << endl;

  //double Ndata = Ns+Nb;
  cout << "******************************"<<endl;
  cout << "Experience avec Signal modifié"<<endl;
  cout << "******************************"<<endl;
  double LLR_data;
  double CLsb = ComputeCL_MultiChannels_hypSB(HistoGedankenExp_SBhyp, Ns, Nb, Ndata, nChannels, LLR_data, "hSB.pdf");
  cout << "CLsb="<<CLsb<<endl;
  double CLb = ComputeCL_MultiChannels_hypB(HistoGedankenExp_Bhyp, Ns, Nb, Ndata, nChannels, LLR_data, "hB.pdf");
  cout << "CLb="<<CLb<<" S="<<GetSignificanceStandardDeviation(CLb)<<endl; //C'est le CLb observe dans les donnees S+B => a utiliser pour la significance
  cout << "CLs=CLsb/CLb="<<CLsb/CLb<<endl;

  DrawHist(HistoGedankenExp_Bhyp, HistoGedankenExp_SBhyp, LLR_data, "hCLs.pdf");

  result[0] = CLsb;
  result[1] = CLb;
  result[2] = CLsb/CLb;

  delete [] GedankenExp_SBhyp;
  delete [] GedankenExp_Bhyp;

  return result;
}

double* GenerateToyExperiment_MultiChannels_withSyst(double* Ns, double* Nb, double* Ns_err, double* Nb_err, double* Ndata, int nChannels, int ntrials){

  double* result = new double[3];

  double** GedankenExp_SBhyp = new double*[nChannels];
  double** GedankenExp_Bhyp = new double*[nChannels];

  for (int i=0; i<nChannels; i++){

    GedankenExp_Bhyp[i] = GetDataPoissonDistribution_Fluctuate(Nb[i], Nb_err[i], ntrials);
    GedankenExp_SBhyp[i] = GetDataPoissonDistribution_Fluctuate(Ns[i], Ns_err[i], ntrials);

  }

  TH1F* HistoGedankenExp_SBhyp = PlotGedankenExpMinus2LLR_MultiChannels("GedankenExp_SBhyp", GedankenExp_SBhyp, Ns, Nb, ntrials, nChannels);
  TH1F* HistoGedankenExp_Bhyp = PlotGedankenExpMinus2LLR_MultiChannels("GedankenExp_Bhyp", GedankenExp_Bhyp, Ns, Nb, ntrials, nChannels);

  cout << endl;
  cout << "******************************"<<endl;
  cout << "Experience avec Signal modifié (Fluctuation)"<<endl;
  cout << "******************************"<<endl;
  double LLR_data;
  double CLsb = ComputeCL_MultiChannels_hypSB(HistoGedankenExp_SBhyp, Ns, Nb, Ndata, nChannels, LLR_data, "hSB_sys.pdf");
  cout<<"LLR_data is "<<LLR_data<<endl;
  double CLb = ComputeCL_MultiChannels_hypSB(HistoGedankenExp_Bhyp, Ns, Nb, Ndata, nChannels, LLR_data, "hB_sys.pdf");
  cout<<"LLR_data is "<<LLR_data<<endl;
  cout << "CLb="<<CLb<<" S="<<GetSignificanceStandardDeviation(CLb)<<endl; //C'est le CLb observe dans les donnees S+B => a utiliser pour la significance
  cout << "CLsb="<<CLsb<<endl;
  cout << "CLs=CLsb/CLb="<<CLsb/CLb<<endl;

  DrawHist(HistoGedankenExp_Bhyp, HistoGedankenExp_SBhyp, LLR_data, "hCLs.pdf");

  result[0] = CLsb;
  result[1] = CLb;
  result[2] = CLsb/CLb;

  delete [] GedankenExp_SBhyp;
  delete [] GedankenExp_Bhyp;
  return result;
}

void cls_multiBin(){
    bool test = true;
    int ntry=1000000;
    gErrorIgnoreLevel = kError;
    // Fill the pseudo-histograms
    
    if(test){
      int n = 2;
      double Ns[] = {9, 6};
      double Ns_err[] = {1, 2};
      double Nb[] = {10, 10};
      double Nb_err[] = {3, 4};
      double Ndata[] = {9, 6};
      double Ndata_err[] = {3, 2};

      GenerateToyExperiment_MultiChannels(Ns, Nb, Nb, n, ntry);
      //GenerateToyExperiment_MultiChannels_withSyst(Ns, Nb, Ns_err, Nb_err, Nb, n, ntry);
    }
    else{
      // double Ndata[] = {2.284340e+00, 1.067280e+00, 5.424560e-01, 2.317930e-01, 8.713560e-02, 3.699980e-02, 1.789660e-02, 8.427990e-03, 3.592930e-03, 1.498080e-03, 6.316090e-04, 2.499870e-04, 7.636430e-05, 1.399750e-05, 1.098470e-06, 2.185840e-07};
      // double Nb[] = {3.937259e-09, 1.476472e-08, 1.968630e-09, 4.921574e-09, 2.952944e-09, 0.000000e+00, 5.905888e-09, 9.843147e-09, 1.814830e-08, 2.829905e-08, 7.382361e-08, 2.167133e-07, 1.407718e-05, 8.503187e-06, 8.722506e-07, 7.657968e-08};
      // double Nb_err[] = {3.682969e-09, 3.812235e-09, 1.392031e-09, 1.841484e-09, 1.392031e-09, 1.704883e-09, 2.088047e-09, 2.982924e-09, 4.489251e-09, 6.698586e-09, 1.218028e-08, 2.075701e-08, 5.633294e-08, 2.760972e-08, 7.937478e-09, 2.338654e-09};
      // double Ns[] = {3.937259e-09, 1.476472e-08, 1.968630e-09, 4.921574e-09, 2.952944e-09, 0.000000e+00, 5.905888e-09, 9.843147e-09, 1.814830e-08, 2.829905e-08, 7.382361e-08, 2.167133e-07, 1.407718e-05, 8.503187e-06, 8.052680e-07, 0.000000e+00};
      // double Ns_err[] = {3.682969e-09, 3.812235e-09, 1.392031e-09, 1.841484e-09, 1.392031e-09, 1.704883e-09, 2.088047e-09, 2.982924e-09, 4.489251e-09, 6.698586e-09, 1.218028e-08, 2.075701e-08, 5.633294e-08, 2.760972e-08, 7.629868e-09, 0.000000e+00};
      int n = 5;
      double Ndata[] = {2.499870e-04, 7.636430e-05, 1.399750e-05, 1.098470e-06, 2.185840e-07};
      double Nb[] = {2.167133e-07, 1.407718e-05, 8.503187e-06, 8.722506e-07, 7.657968e-08};
      double Nb_err[] = {2.075701e-08, 5.633294e-08, 2.760972e-08, 7.937478e-09, 2.338654e-09};
      double Ns[] = {2.167133e-07, 1.407718e-05, 8.503187e-06, 8.052680e-07, 2.000000e-09};
      double Ns_err[] = {2.075701e-08, 5.633294e-08, 2.760972e-08, 7.629868e-09, 0.000000e+00};


      //Convert into a number of event
      double Energy_bin[] = {1.5e+02, 2e+02, 4e+02, 5e+02, 5e+02};
      //double Energy_bin[] = {0.25e+02, 0.25e+02, 0.25e+02, 0.5e+02, 0.5e+02, 0.5e+02, 0.5e+02, 0.7e+02, 0.8e+02, 1e+02, 1e+02, 1.5e+02, 2e+02, 4e+02, 5e+02, 5e+02};
      double luminosity = 36e03;
      double factor;

      for(int i = 0; i<n; i++){

        factor = Energy_bin[i]*luminosity;
        Ndata[i] *= factor;
        cout<<"Ndata : "<<Ndata[i]<<endl;
        Nb[i] *=  factor;
        cout<<"Nb : "<<Nb[i]<<endl;
        Ns[i] *= factor;
        cout<<"Ns : "<<Ns[i]<<endl;
        Nb_err[i] *= factor;
        Ns_err[i] *= factor;

      }
      GenerateToyExperiment_MultiChannels(Ns, Nb, Nb, n, ntry);
      //GenerateToyExperiment_MultiChannels_withSyst(Ns, Nb, Ns_err, Nb_err, Nb, n, ntry);
    }

    // GenerateToyExperiment_MultiChannels(Ns, Nb, Ndata, n, 100000);
    // GenerateToyExperiment_MultiChannels_withSyst(Ns, Nb, Ns_err, Nb_err, Ndata, n, 100000);
    // GenerateToyExperiment_MultiChannels(Ns, Nb, Nb, n, 100000);
    // GenerateToyExperiment_MultiChannels_withSyst(Ns, Nb, Ns_err, Nb_err, Nb, n, 100000);
}