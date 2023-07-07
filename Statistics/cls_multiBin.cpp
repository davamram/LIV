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

  TH1F* GedankenExpMinus2LLR = new TH1F(Name.c_str(), Name.c_str(), 10000, -100, 100);

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

double ComputeCL_MultiChannels_hypSB(TH1F* HistoGedankenExp, double* Ns, double* Nb, double* Ndata, int nChannels){

  double xmin = HistoGedankenExp->GetXaxis()->GetXmin();
  double xmax = HistoGedankenExp->GetXaxis()->GetXmax();
  double binwidth = (xmax-xmin)/((double)HistoGedankenExp->GetNbinsX());

  double LLR_data = 0;
  double LLR_formula = 0;
  for (int i=0; i<nChannels; i++){
  
    LLR_data += CalcLogLikelihoodRatio(Ndata[i], Ns[i], Nb[i]);
    
  }

  int bin1 = (-2*LLR_data - xmin)/binwidth;
  int bin2 = (xmax - xmin)/binwidth;

  double CL = HistoGedankenExp->Integral(bin1, bin2);

  return CL;
}

double* GenerateToyExperiment_MultiChannels(double* Ns, double* Nb, double* Ndata, int nChannels, int ntrials){

  double result[6];

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

  cout << "Experience avec Signal modifié"<<endl;
  double CLsb = ComputeCL_MultiChannels_hypSB(HistoGedankenExp_SBhyp, Ns, Nb, Ndata, nChannels);
  double CLb = ComputeCL_MultiChannels_hypSB(HistoGedankenExp_Bhyp, Ns, Nb, Ndata, nChannels);
  cout << "CLb="<<CLb<<" S="<<GetSignificanceStandardDeviation(CLb)<<endl; //C'est le CLb observe dans les donnees S+B => a utiliser pour la significance
  cout << "CLsb="<<CLsb<<endl;
  cout << "CLs=CLsb/CLb="<<CLsb/CLb<<endl;

  result[0] = CLsb;
  result[1] = CLb;
  result[2] = CLsb/CLb;

  return result;
}

void cls_multiBin(){
    int n = 2;
    double Ns[] = {8, 7};
    double Nb[] = {10, 10};
    double Ndata[] = {9, 6};
    GenerateToyExperiment_MultiChannels(Ns, Nb, Ndata, n, 100000);
}