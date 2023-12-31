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
#include <ctime>
#include "ExtractValues.h"
#include <Math/SpecFuncMathCore.h>



Double_t PoissonCdf(Double_t *x, Double_t *params){
  Double_t A = params[0];
  Double_t arg = x[0];
  return ROOT::Math::inc_gamma_c(arg+1, A);
}

double GetContinuousPoissonCdf(double mean, double rand){
  double a = mean - 5*TMath::Sqrt(mean) < 0 ? 0 : mean - 5*TMath::Sqrt(mean);
  double b = mean + 5*TMath::Sqrt(mean);
  TF1 *f = new TF1("f", PoissonCdf, a, b, 1);

  f->SetParameter(0, mean);

  if(rand <= f->Eval(a)) return a;
  if(rand >= f->Eval(b)) return b;
  return f->GetX(rand);
}

Double_t ContinuousPoisson(Double_t *x, Double_t *params){
  Double_t A = params[0];
  if(A<=0) return 0;
  Double_t arg = x[0];
  return TMath::Exp(-A) * TMath::Power(A, arg) / TMath::Gamma(arg + 1);
}

double GetContinuousPoisson(double mean){
  TF1 *f = new TF1("f", ContinuousPoisson, 0, 3*ceil(mean), 1);
  f->SetParameter(0, mean);
  double rand = f->GetRandom();
  delete f;
  return rand;
}


double* GetPoissonDistribution(double mean, int ntrials){

  TRandom3* rnd = new TRandom3();
  rnd->SetSeed(10111999);
  double* P = new double[ntrials];

  for (int i=0; i<ntrials; i++){
    //P[i] = rnd->PoissonD(mean);
    P[i] = GetContinuousPoissonCdf(mean, rnd->Uniform(0,1));
  }

  return P;
}

double* GetDataPoissonDistribution_Fluctuate(double N, double N_err_u, double N_err_d, int ntrials, double Ndata){

  double mean;

  double N_fluctuated_u; double N_fluctuated_d; double N_fluctuated;
  
  // TH1F *histo = new TH1F("Distribution", "Distribution", 100, N-5*TMath::Sqrt(N), N+5*TMath::Sqrt(N));
  TH1F *histo = new TH1F("Distribution", "Distribution", 100, N/5, 3*N);

  TH1F *histoG = new TH1F("DistributionGaus", "DistributionGaus", 100, N-3*N_err_d, N+3*N_err_u);

  TRandom3* rnd = new TRandom3(10111999);
  double* P = new double[ntrials];

  double int_d = ROOT::Math::normal_cdf_c(0, N_err_d, N)-0.5;
  double int_tot = ROOT::Math::normal_cdf_c(0, N_err_d, N);
  double normFactor = int_d/int_tot;

  for (int i=0; i<ntrials; i++){

    N_fluctuated_u=0;
    N_fluctuated_d=0;
    N_fluctuated=0;

    double draw = rnd->Uniform(0,1);
    while (N_fluctuated<=0){
      do{
        N_fluctuated_d = rnd->Gaus(N, N_err_d);
      } while((N_fluctuated_d >= N || N_fluctuated_d <= 0) && N!=0);
      do{
        N_fluctuated_u = rnd->Gaus(N, N_err_u);
      } while(N_fluctuated_u <= N);
      // See detail in presentation -> This ensure that i have a continuous distribution
      N_fluctuated = draw < normFactor ? N_fluctuated_d : N_fluctuated_u;
    }

    //P[i] = rnd->PoissonD(N_fluctuated);

    P[i] = GetContinuousPoissonCdf(N_fluctuated, rnd->Uniform(0,1));

    histo->Fill(P[i]);
    histoG->Fill(N_fluctuated);

  }
  TCanvas *canvas = new TCanvas("Nom_du_canvas", "Titre_du_canvas");
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  TLine *mcLine = new TLine(N, 0, N, histo->GetMaximum());
  TLine *dataLine = new TLine(Ndata, 0, Ndata, histo->GetMaximum());
  mcLine->SetLineColor(kRed);
  dataLine->SetLineColor(kGreen);
  TString nameData = Form("N_data = %.2f", Ndata);
  TString nameMC = Form("N_mc = %.2f", N);
  legend->AddEntry(mcLine, nameMC, "l");
  legend->AddEntry(dataLine, nameData, "l");
  histo->Draw();
  mcLine->Draw("SAME");
  dataLine->Draw("SAME");
  legend->Draw("SAME");
  canvas->Update();
  TString nomFichier = Form("Poisson/Poisson_Distribution_%.2f.png", N);
  canvas->SaveAs(nomFichier);
  histoG->Draw();
  mcLine->Draw("SAME");
  dataLine->Draw("SAME");
  legend->Draw("SAME");
  canvas->Update();
  nomFichier = Form("Gauss/Gauss_Distribution_%.2f.png", N);
  canvas->SaveAs(nomFichier);
  canvas->Close();
  delete histo;
  delete histoG;
  delete canvas;

  return P;
}

double CalcLikelihood(double n, double x){

  double l = exp(-x)*TMath::Power(x,n)/TMath::Gamma(n+1);
  return l;
}

double CalcLogLikelihoodRatioTrue(double Ndata, double Ns, double Nb){

  double lr;
  if (Nb>0) lr = CalcLikelihood(Ndata, Ns)/CalcLikelihood(Ndata, Nb);
  return log(lr);
}

double CalcLogLikelihoodRatio(double Ndata, double Ns, double Nb){
  // If the bin value is 0, i should merge it xith the previous one. To discuss with Nicolas
  // if(Nb==0 || Ns==0) cout<<"Bin value is 0. LLR impossible to compute. Setting bin value to 1e-15"<<endl;
  if(Nb==0) Nb=1e-15;
  if(Ns==0) Ns=1e-15;
  double llr = Nb-Ns+Ndata*log(Ns/Nb);
  return llr;
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
      // if(LLR<-50) cout<<"Pseudo LLR : "<<LLR<<" with Ndata, Ns, Nb = "<<GedankenExp[j][i]<<" ; "<<Ns[j]<<" ; "<<Nb[j]<<endl;
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
    // cout<<"LLR step is in SB : "<<LLR_data<<" with Ndata, Ns, Nb = "<<Ndata[i]<<" ; "<<Ns[i]<<" ; "<<Nb[i]<<endl;
    
  }

  double xmin = HistoGedankenExp->GetXaxis()->GetXmin();
  double xmax = HistoGedankenExp->GetXaxis()->GetXmax();
  double binwidth = (xmax-xmin)/((double)HistoGedankenExp->GetNbinsX());

  int bin1 = (-2*LLR_data - xmin)/binwidth;
  int bin2 = (xmax - xmin)/binwidth;
  cout<<"LLR is : "<<LLR_data<<endl;
  double CL = HistoGedankenExp->Integral(bin1, bin2);

  return CL;
}

double ComputeCL_MultiChannels_hypB(TH1F* HistoGedankenExp, double* Ns, double* Nb, double* Ndata, int nChannels, double &LLR_data, string name){

  LLR_data = 0;
  double LLR_formula = 0;
  for (int i=0; i<nChannels; i++){
  
    LLR_data += CalcLogLikelihoodRatio(Ndata[i], Ns[i], Nb[i]);
    // cout<<"LLR step is in B : "<<LLR_data<<" with Ndata, Ns, Nb = "<<Ndata[i]<<" ; "<<Ns[i]<<" ; "<<Nb[i]<<endl;

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
  hH0->SetLineColor(kGreen);
  hH0->SetLineWidth(2);
  hH1->SetLineColor(kRed);
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

  TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(hH0, "Expected for SM", "l");
  legend->AddEntry(hH1, "Expected for LIV", "l");
  // legend->AddEntry(hH0_color, "P(SM)", "f");
  // legend->AddEntry(hH1_color, "P(LIV)", "f");
  legend->AddEntry(line1, "Observed in data", "l");
  legend->Draw();



  c1->SaveAs(name.c_str());
  delete c1;
}

double* GenerateToyExperiment_MultiChannels(double* Ns, double* Nb, double* Ndata, int nChannels, int ntrials){

  double* result = new double[3];

  double** GedankenExp_SBhyp = new double*[nChannels];
  double** GedankenExp_Bhyp = new double*[nChannels];

  for (int i=0; i<nChannels; i++){
    cout<<"Channel "<<i<<" ...";
    GedankenExp_SBhyp[i] = GetPoissonDistribution(Ns[i], ntrials);
    GedankenExp_Bhyp[i] = GetPoissonDistribution(Nb[i], ntrials);
    cout<<" DONE"<<endl;
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
  double CLb = 1 - ComputeCL_MultiChannels_hypB(HistoGedankenExp_Bhyp, Ns, Nb, Ndata, nChannels, LLR_data, "hB.pdf");
  cout << "CLb="<<CLb<<" S="<<GetSignificanceStandardDeviation(CLb)<<endl; //C'est le CLb observe dans les donnees S+B => a utiliser pour la significance
  CLb = ComputeCL_MultiChannels_hypSB(HistoGedankenExp_Bhyp, Ns, Nb, Ndata, nChannels, LLR_data, "hB.pdf");
  cout << "CLb2="<<CLb<<" S="<<GetSignificanceStandardDeviation(CLb)<<endl; //C'est le CLb observe dans les donnees S+B => a utiliser pour la significance
  cout << "CLs=CLsb/CLb="<<CLsb/CLb<<endl;

  DrawHist(HistoGedankenExp_Bhyp, HistoGedankenExp_SBhyp, LLR_data, "hCLs.pdf");

  result[0] = CLsb;
  result[1] = CLb;
  result[2] = CLsb/CLb;

  delete [] GedankenExp_SBhyp;
  delete [] GedankenExp_Bhyp;

  return result;
}

double* GenerateToyExperiment_MultiChannels_withSyst(double* Ns, double* Nb, double* Ns_err_u, double* Ns_err_d, double* Nb_err_u, double* Nb_err_d, double* Ndata, int nChannels, int ntrials, string Energie){

  clock_t start = clock();

  double* result = new double[3];

  double** GedankenExp_SBhyp = new double*[nChannels];
  double** GedankenExp_Bhyp = new double*[nChannels];

  for (int i=0; i<nChannels; i++){
    GedankenExp_Bhyp[i] = GetDataPoissonDistribution_Fluctuate(Nb[i], Nb_err_u[i], Nb_err_d[i], ntrials, Ndata[i]);
    GedankenExp_SBhyp[i] = GetDataPoissonDistribution_Fluctuate(Ns[i], Ns_err_u[i], Ns_err_d[i], ntrials, Ndata[i]);

  }

  TH1F* HistoGedankenExp_SBhyp = PlotGedankenExpMinus2LLR_MultiChannels("GedankenExp_SBhyp", GedankenExp_SBhyp, Ns, Nb, ntrials, nChannels);
  TH1F* HistoGedankenExp_Bhyp = PlotGedankenExpMinus2LLR_MultiChannels("GedankenExp_Bhyp", GedankenExp_Bhyp, Ns, Nb, ntrials, nChannels);

  cout << endl;
  cout << "******************************"<<endl;
  cout << "Experience avec Signal modifié (Fluctuation)"<<endl;
  cout << "******************************"<<endl;
  double LLR_data;
  double CLsb = ComputeCL_MultiChannels_hypSB(HistoGedankenExp_SBhyp, Ns, Nb, Ndata, nChannels, LLR_data, "hSB_sys.pdf");
  double CLb = 1 - ComputeCL_MultiChannels_hypB(HistoGedankenExp_Bhyp, Ns, Nb, Ndata, nChannels, LLR_data, "hB_sys.pdf");
  cout << "CLb="<<CLb<<" S="<<GetSignificanceStandardDeviation(CLb)<<endl; //C'est le CLb observe dans les donnees S+B => a utiliser pour la significance
  cout << "CLsb="<<CLsb<<endl;
  cout << "CLs=CLsb/CLb="<<CLsb/CLb<<endl;

  string Name = "hCLs_"+Energie+".pdf";
  DrawHist(HistoGedankenExp_Bhyp, HistoGedankenExp_SBhyp, LLR_data, Name);

  result[0] = CLsb;
  result[1] = CLb;
  result[2] = CLsb/CLb;

  delete [] GedankenExp_SBhyp;
  delete [] GedankenExp_Bhyp;
  
  double time;

  cout << "Done in : " << (double)(clock() - start) / CLOCKS_PER_SEC << " s" << endl;

  return result;
}

void cls_multiBin(int eValue){

    gROOT->SetBatch(kTRUE);
    gSystem->Load("biblio.so");

    //int eValue = 2050;

    bool test = false;
    int ntry=100000;
    gErrorIgnoreLevel = kError;
    // Fill the pseudo-histograms

    //HepData Numbers
    //Je pourrais faire des dictionnaires mais pas envie de perdre du temps pour ça
    int n = 16;

    double Atlas_down[] = {-7.827176e-02, -3.551370e-02, -1.697738e-02, -7.084197e-03, -2.575804e-03, -1.157590e-03, -5.549091e-04, -2.714763e-04, -1.185141e-04, -5.092944e-05, -2.129382e-05, -8.701203e-06, -2.891808e-06, -5.895342e-07, -7.166077e-08, -8.281072e-09};
    double Atlas_up[] = {5.233613e-02, 2.433757e-02, 1.214272e-02, 5.235209e-03, 1.963962e-03, 8.571483e-04, 4.382003e-04, 2.079510e-04, 9.104776e-05, 3.976086e-05, 1.729810e-05, 7.207799e-06, 2.486231e-06, 5.904051e-07, 7.280856e-08, 9.336763e-09};
    double Atlas_stat[] = {0.00295715, 0.00177455, 0.000834419, 0.000268998, 0.000112464, 6.68152e-05, 3.31831e-05, 1.32362e-05, 6.20673e-06, 2.69945e-06, 1.42323e-06, 5.45151e-07, 1.94116e-07, 3.17007e-08, 4.10542e-09, 5.15698e-10};

    double Nb[] = {2.224240e+00, 1.017430e+00, 5.141330e-01, 2.227860e-01, 8.262710e-02, 3.580690e-02, 1.729010e-02, 8.139830e-03, 3.521510e-03, 1.464270e-03, 6.108830e-04, 2.379370e-04, 7.349860e-05, 1.409570e-05, 1.395380e-06, 1.220970e-07};
    double Nb_err[] = {7.074313e-02, 4.497576e-02, 2.768121e-02, 1.581030e-02, 7.730480e-03, 7.547453e-03, 3.877605e-03, 2.543693e-03, 1.195816e-03, 8.785437e-04, 1.057928e-03, 3.362813e-04, 5.268441e-05, 1.021463e-07, 2.762140e-08, 1.643521e-09};
    double Nb_err_d[n]; double Nb_err_u[n];

    double Ndata[] = {2.284340e+00, 1.067280e+00, 5.424560e-01, 2.317930e-01, 8.713560e-02, 3.699980e-02, 1.789660e-02, 8.427990e-03, 3.592930e-03, 1.498080e-03, 6.316090e-04, 2.499870e-04, 7.636430e-05, 1.399750e-05, 1.098470e-06, 2.185840e-07};
    double Ndata_err[] = {8.230996e-02, 3.598052e-02, 1.722778e-02, 7.057156e-03, 2.743396e-03, 1.133563e-03, 5.585615e-04, 2.722863e-04, 1.227772e-04, 5.432952e-05, 2.607210e-05, 1.193473e-05, 4.651671e-06, 1.256509e-06, 2.793959e-07, 1.267303e-07};
    double Ndata_stat[] = {0.0060483, 0.00124098, 0.000875161, 0.000396898, 0.000243451, 0.000156815, 0.000108274, 6.29188e-05, 3.85399e-05, 2.22023e-05, 1.45279e-05, 7.53505e-06, 3.58288e-06, 1.11595e-06, 2.74616e-07, 1.262e-07};

    std::vector<double> values = getValues(eValue, 2);
    std::vector<double> errs = getValues(eValue, 3);
    if(values.size()!=n || errs.size()!=n) cout<<"Big problem"<<endl;
    double Ns[n];
    double Ns_err[n];
    for (int i=0; i<n; i++){
        Ns[i]=values[i];
        Ns_err[i]=errs[i];
    }
    double Ns_err_d[n]; double Ns_err_u[n];

    double Energy_bin[] = {0.25e+02, 0.25e+02, 0.25e+02, 0.5e+02, 0.5e+02, 0.5e+02, 0.5e+02, 0.7e+02, 0.8e+02, 1e+02, 1e+02, 1.5e+02, 2e+02, 4e+02, 5e+02, 5e+02};
    double luminosity = 36e03;
    double factor;

    for(int i = 0; i<n; i++){
    // Convert in number of events and add systematics in quadrature
    factor = Energy_bin[i]*luminosity;

    Atlas_up[i] *= factor;
    Atlas_down[i] *= factor;
    Atlas_stat[i] *= factor;

    Ndata[i] *= factor;
    Ndata_err[i] *= factor;
    Ndata_stat[i] *= factor;


    Nb[i] *=  factor;
    Nb_err[i] *= factor;
    Nb_err_u[i] = sqrt( Nb_err[i]*Nb_err[i] + Atlas_up[i]*Atlas_up[i] + Ndata_err[i]*Ndata_err[i] - Atlas_stat[i]*Atlas_stat[i] - Ndata_stat[i]*Ndata_stat[i] );
    Nb_err_d[i] = sqrt( Nb_err[i]*Nb_err[i] + Atlas_down[i]*Atlas_down[i] + Ndata_err[i]*Ndata_err[i] - Atlas_stat[i]*Atlas_stat[i] - Ndata_stat[i]*Ndata_stat[i] );

    Ns[i] *= factor;
    Ns_err[i] *=factor;
    Ns_err_d[i] = sqrt( Ns_err[i]*Ns_err[i] + Atlas_down[i]*Atlas_down[i] + Ndata_err[i]*Ndata_err[i] - Atlas_stat[i]*Atlas_stat[i] - Ndata_stat[i]*Ndata_stat[i] );
    Ns_err_u[i] = sqrt( Ns_err[i]*Ns_err[i] + Atlas_up[i]*Atlas_up[i] + Ndata_err[i]*Ndata_err[i] - Atlas_stat[i]*Atlas_stat[i] - Ndata_stat[i]*Ndata_stat[i] );

    }
    cout<<"\n#### Results for E="+to_string(eValue)+" ####"<<endl;
    GenerateToyExperiment_MultiChannels_withSyst(Ns, Nb, Ns_err_u, Ns_err_d, Nb_err_u, Nb_err_d, Ndata, n, ntry, to_string(eValue));
}
