#ifndef SELECTPHOTON_H
#define SELECTPHOTON_H

#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/LeadingParticlesFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

#include <iostream>
#include <functional>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctime>

struct Data{
	double mFermion;
	double KTr;
	double c;
	double hbar;
	double alpha;
	double CrossSection;
	double Luminosity;
	std::vector <double> Distance;
	double Masses[9]={511e-6,105.7e-3,1.777,2.16e-3,4.67e-3,93e-3,1.27,4.18,172.5}; //Masses des fermions du SM (hors neutrinos) dans l'ordre e,mu,tau,u,d,s,c,b,t
	std::string Names[9]={"electron","muon","tau","up quark","down quark","strange quark","charm quark","bottom quark","top quark"};
	int Identifier;
	double ThresholdEnergy;
	double TriggerEfficiency;
	double IdentificationEfficiency;
	double ReconstructionEfficiency;
}; //Structure contenant les constantes du probleme

double integrate(std::function<double(double)> func, double a, double b, int numSteps);
double LIVParameter(int particule, Data d);
void initialize(Data &d);
double CalcGamma(double E, double eb, double eh, Data d);
double integrate(Data d, double E, double a, double b);
double partialwidth(double x, Data d, double E);
double Threshold(Data d);
double SurvivalProb(double gamma, double x);
double EBar(double egamma, Data d);
double DispersionRelationPhoton(double egamma, Data d);
double DispersionRelationFermion(double efermion, Data d);
double costheta(double epart, double egamma, Data d);
void normalize(double (&vec)[3]);
void generatebase(double px, double py, double pz, double (&pun)[3], double (&pdeux)[3], double (&ptrois)[3]);
double angle(double normphoton, double normfermion, double normantifermion, double cosangle);
int HasSurvived(Rivet::FourMomentum &momFermion, Rivet::FourMomentum &Fermion, Rivet::FourMomentum &AntiFermion);
double calcDeltaPhi(Rivet::FourMomentum Fermion, Rivet::FourMomentum AntiFermion);
double Reweight(double Et);

#endif