#include "lib/SelectPhoton.hh"
using namespace std;

double LIVParameter(int particule, Data d) //Revoie le parametre d'invariance de Lorentz pour avoir le seuil souhaite, compte-tenu de la particule consideree
{
	return 1/(1-((d.ThresholdEnergy*d.ThresholdEnergy)/(2*d.Masses[particule]*d.Masses[particule])));
}

void initialize(Data &d) //Initialise les valeurs
{
	d.Identifier=8; //Indice de la particule dans Masses (ex : 0 pour un electron, 3 pour un quark up, etc.) 8 quark top ?
	d.ThresholdEnergy=2200; //Energie seuil souhaitee pour la desintegration en GeV
    d.mFermion=d.Masses[d.Identifier]; //Masse du fermion en GeV/c2
    d.KTr=LIVParameter(d.Identifier,d); //Parametre de violation de l'invariance de Lorentz //Ordre de grandeur : -1e-13 pour l'electron et -1e-2 pour le top
    d.hbar=1;
    d.c=1;
    d.alpha=1.0/137;
	d.CrossSection=2.163e-02*1000; // Cross-section en fb
	d.Luminosity=36; //Luminosite en fb-1
	d.Distance.push_back(10); //Distance(s) de vol en um
	d.Distance.push_back(100);
	d.Distance.push_back(200);
	d.Distance.push_back(1000);
	d.TriggerEfficiency=1.00; //Efficacite de trigger pour un photon
	d.IdentificationEfficiency=0.95; //Efficacite d'identification pour un photon
	d.ReconstructionEfficiency=0.70; //Efficacite de reconstruction pour un photon
}

double integrate(Data d, double E, double a, double b)
{
    int numSteps = 1000;
    double h = (b - a) / numSteps;  // Largeur de chaque subdivision
    double integral = 0.0;

    for (int i = 0; i < numSteps; i++) {
        double x1 = a + i * h;      // Point de départ du segment
        double x2 = a + (i + 1) * h;  // Point de fin du segment
        double fx1 = partialwidth(x1, d, E);      // Valeur de la fonction au point de départ
        double fx2 = partialwidth(x2, d, E);      // Valeur de la fonction au point de fin
        integral += (fx1 + fx2) * h / 2.0;  // Approximation de l'intégrale du segment par un trapèze
    }

    return integral;
}

double partialwidth(double x, Data d, double E) //Formule de la largeur de desintegration partielle implementee pour etre utilisee avec TF1
{
    double par[7]; //masse fermion, KTr, hbar, c, alpha, Egamma, 1
    par[0] = d.mFermion; par[1] = d.KTr; par[2] = d.hbar; par[3] = d.c; par[4]=d.alpha; par[5] = E; par[6] = 1;
    double me2c4 = par[0]*par[0]*par[3]*par[3]*par[3]*par[3]; //Les parametres sont dans l'ordre de Data (donc 0 pour mFermion, 1 pour KTr, ...), l'avant-dernier est EGamma
    double num = par[4]*((1-par[1])*(2*par[1]*x*(par[5]-x)+(1+par[1])*me2c4)-par[1]*par[5]*par[5]);
    double den = par[2]*(1+par[1])*(1+par[1])*(sqrt(1-par[1]*par[1]))*par[5]*par[5];
    return num/(den*par[6]); //Le dernier est un facteur de normalisation
}

double Threshold(Data d) //Renvoie l'energie seuil a partir des donnees du probleme
{
	return (2*d.mFermion*d.c*d.c*sqrt((1-d.KTr)/(-2*d.KTr)));
}

double SurvivalProb(double gamma, double x) //Renvoie la probabilite qu'un photon survive
{
	return exp(-gamma*x*5.07e9); //Avec x en um
}

double EBar(double egamma, Data d) //Variable dont depend le domaine de validite de cos(theta)
{
	return sqrt((1+d.KTr)/(1-d.KTr)*(egamma*egamma+2*(1/d.KTr-1)*d.mFermion*d.mFermion));
}

double DispersionRelationPhoton(double egamma, Data d) //Renvoie la norme de PGamma a partir de l'energie et de la relation de dispersion
{
        return (egamma*sqrt(1+d.KTr))/(d.c*sqrt(1-d.KTr));
}

double DispersionRelationFermion(double efermion, Data d) //Renvoie la norme de PFermion a partir de l'energie et de la relation de dispersion
{
        return sqrt(efermion*efermion-d.mFermion*d.mFermion);
}

double costheta(double epart, double egamma, Data d) //Renvoie le cosinus de  l'angle entre PFermion et PAntiFermion (angle theta)
{
        double me2c4=d.mFermion*d.mFermion*d.c*d.c*d.c*d.c;
        double num=epart*(egamma-epart)+((d.KTr)/(1-d.KTr))*egamma*egamma+me2c4;
        double den=sqrt((epart*epart-me2c4)*((egamma-epart)*(egamma-epart)-me2c4));
        return num/den;
}

void normalize(double (&vec)[3]) //Normalise un 3-vecteur
{
        double norme=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
        for(int i=0;i<3;i++)
        {
                vec[i]/=norme;
        }
}

void generatebase(double px, double py, double pz, double (&pun)[3], double (&pdeux)[3], double (&ptrois)[3]) //Genere la base ORTHONORMEE autour de PGamma
{
        pun[0]=px;
        pun[1]=py;
        pun[2]=pz;
        pdeux[0]=py;
        pdeux[1]=-px;
        pdeux[2]=0;
        ptrois[0]=px*pz;
        ptrois[1]=py*pz;
        ptrois[2]=-(px*px+py*py);
        normalize(pun);
        normalize(pdeux);
        normalize(ptrois);
}

double angle(double normphoton, double normfermion, double normantifermion, double cosangle) //Renvoie l'angle entre PGamma et PFermion (angle alpha)
{
        if((normfermion+normantifermion*cosangle)/normphoton > 1) return 0; //This did not happened in ROOT, but here sometime I have cos(alpha)>1 (1.000001)
        return acos((normfermion+normantifermion*cosangle)/normphoton);
}

// Check if the photon has survived and if not, if the fermions created are reconstructed as photon. If so, the photon momentum become the fermion momentum
int HasSurvived(Rivet::FourMomentum &momPhoton, Rivet::FourMomentum &Fermion, Rivet::FourMomentum &AntiFermion){
    Data d;
    initialize(d);
    // Only photon above Etresh would desintegrated in LIV
    double E = momPhoton.E();
    
    if(E<d.ThresholdEnergy) return 1;

    double EBasse = 0.5*(E-EBar(E,d)); //On determine les bornes d'intégration
    double EHaute = 0.5*(E+EBar(E,d));
    double Gamma = integrate(d, E, EBasse, EHaute);
    // This is a problem : every number generated during the same second are the same. Need to find an other way
    srand(static_cast<unsigned int>(std::time(0)));
    // Génération d'un nombre aléatoire entre 0 et 1
    double Draw = static_cast<double>(std::rand()) / RAND_MAX;
    if(Draw<1-SurvivalProb(Gamma, d.Distance[d.Distance.size()-1])){
        
        double EFermion = rand()/RAND_MAX*(EHaute-EBasse)+EBasse;
        double EAntiFermion = E-EFermion;
        
        if(EFermion < 0 || EAntiFermion <0) cout<<"Error : Negative Energy. EBasse = "<<EBasse<<" | EHaute = "<<EHaute<<endl;

        double PGamma[3]; double PPerpendicular[3]; double P3[3]; double PFermion[3]; double PAntiFermion[3];
        double NormPGamma=DispersionRelationPhoton(E,d);
        double NormPFermion=DispersionRelationFermion(EFermion,d);
        double NormPAntiFermion=DispersionRelationFermion(EAntiFermion,d);
        generatebase(momPhoton.px(),momPhoton.py(),momPhoton.pz(),PGamma,PPerpendicular,P3);
        double cosangle=costheta(EFermion,E,d);
        double alpha=angle(NormPGamma,NormPFermion,NormPAntiFermion,cosangle);
        double phi=rand()*2*M_PI;
        for(int n=0;n<3;n++) //Reconstruction de PFermion
        {
            PFermion[n]=NormPFermion*(PGamma[n]*cos(alpha)+sin(alpha)*(cos(phi)*PPerpendicular[n]+sin(phi)*P3[n]));
            // cout<<"PFermion "<<n<<" is : "<<PFermion[n]<<endl;
            // cout<<NormPFermion<<" / "<<PGamma[n]<<" / "<<alpha<<" / "<<phi<<" / "<<PPerpendicular[n]<<" / "<<P3[n]<<endl;
            // cout<<"***"<<endl;
            PAntiFermion[n]=(PGamma[n]*NormPGamma)-PFermion[n];
        }
        Fermion.setXYZE(PFermion[0], PFermion[1], PFermion[2], EFermion);
        AntiFermion.setXYZE(PAntiFermion[0], PAntiFermion[1], PAntiFermion[2], EAntiFermion);
        // Disintegration probability
        double desE = std::rand() / RAND_MAX;
        double desP = std::rand() / RAND_MAX;

        if(desE < 0.015 && desP<0.015 && Fermion.pT()>1000 && AntiFermion.pT()>1000){
          momPhoton = Fermion.pT() > AntiFermion.pT() ? Fermion : AntiFermion;
          return 10;
        }
        else if(desE < 0.015 && Fermion.pT()>1000){
          momPhoton = Fermion;
          return 10;
        }
        else if(desP < 0.015 && AntiFermion.pT()>1000){
          momPhoton = AntiFermion;
          return 10;
        }
        return 0;
    }
    return 1;
}