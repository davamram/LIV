#include <iostream>
#include <functional>
#include <cmath>
#include <vector>
#include <string>
#include <cstdlib>
#include <ctime>


double LIVParameter(double E) //Revoie le parametre d'invariance de Lorentz pour avoir le seuil souhaite, compte-tenu de la particule consideree
{
	return 1/(1-((E*E)/(2*0.000511*0.000511)));
}


double partialwidth(double x, double E, double etr) //Formule de la largeur de desintegration partielle implementee pour etre utilisee avec TF1
{
    double par[7]; //masse fermion, KTr, hbar, c, alpha, Egamma, 1
    par[0] = 0.000511; par[1] = LIVParameter(etr); par[2] = 1; par[3] = 1; par[4]=1.0/137; par[5] = E; par[6] = 1;
    double me2c4 = par[0]*par[0]*par[3]*par[3]*par[3]*par[3]; //Les parametres sont dans l'ordre de Data (donc 0 pour mFermion, 1 pour KTr, ...), l'avant-dernier est EGamma
    double num = par[4]*((1-par[1])*(2*par[1]*x*(par[5]-x)+(1+par[1])*me2c4)-par[1]*par[5]*par[5]);
    double den = par[2]*(1+par[1])*(1+par[1])*(sqrt(1-par[1]*par[1]))*par[5]*par[5];
    // std::cout<<num<<" / "<<den<<std::endl;
    return num/(den*par[6]); //Le dernier est un facteur de normalisation
}

double integrate(double E, double a, double b, double etr)
{
    int numSteps = 1000;
    double h = (b - a) / numSteps;  // Largeur de chaque subdivision
    double integral = 0.0;

    for (int i = 0; i < numSteps; i++) {
        double x1 = a + i * h;      // Point de départ du segment
        double x2 = a + (i + 1) * h;  // Point de fin du segment
        double fx1 = partialwidth(x1, E, etr);      // Valeur de la fonction au point de départ
        double fx2 = partialwidth(x2, E, etr);      // Valeur de la fonction au point de fin
        integral += (fx1 + fx2) * h / 2.0;  // Approximation de l'intégrale du segment par un trapèze
            }

    return integral;
}

double EBar(double egamma, double etr) //Variable dont depend le domaine de validite de cos(theta)
{
	return sqrt((1+LIVParameter(etr))/(1-LIVParameter(etr))*(egamma*egamma+2*(1/LIVParameter(etr)-1)*0.000511*0.000511));
}

double Gamma(double E, double eb, double eh, double etr){
    double alpha = 1.0/137;
    double k = LIVParameter(etr);
    double m = 0.000511;
    double num = alpha*((1-k)*(1+k)*m*m*(eh-eb)-k*E*E*(eh-eb)+(1-k)*(k*(E*(eh*eh-eb*eb)-2*(eh*eh*eh-eb*eb*eb)/3)));
    double den = (1+k)*(1+k)*sqrt(1-k)*E*E;
    return num/den;
}

int main(){

    double E=2000;
    double etr=1900;
    double ebar=EBar(E,etr);
    std::cout<<"Ebar : "<<ebar<<std::endl;
    std::cout<<"Kappa : "<<LIVParameter(E)<<std::endl;
    std::cout<<"Gamma is : "<<integrate(E, 1.0/2*(E-ebar), 1.0/2*(E+ebar), etr)<<std::endl;
    std::cout<<"Gamma is also (analytique): "<<Gamma(E, 1.0/2*(E-ebar), 1.0/2*(E+ebar), etr)<<std::endl;
    return 0;
}