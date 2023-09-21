import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
from math import exp
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

def calcEbar(Eg, mf, K):
    return sqrt( (1+K)/(1-K) * (Eg**2 +2*(1/K-1)*mf**2))

def calcKappa(Etr, m):
    return 1/(1-(Etr**2/(2*m**2)))

def calcCosTheta(Ef, Eg, K, mt):
    num = Ef*(Eg-Ef) + K/(1-K)*Eg**2 + mt**2
    den = sqrt( (Ef**2 - mt**2) * ((Eg-Ef)**2 -mt**2))
    return num/den

def fillTab(Eg, Ef, K, mt):
    Result=[]
    result=[]

    for eg in Eg:
        result=[]
        Ebar = calcEbar(eg, mt, K)
        for ef in Ef:
            if( ef < 1/2*(eg+Ebar) and ef > 1/2*(eg-Ebar)):
                result.append(calcCosTheta(ef, eg, K, mt))
            else:
                result.append(0)
        Result.append(result)
    return Result

def makeHeatMap(Etr, mt, v_max, name):

    nbpoint = 1000
    Ktr = calcKappa(Etr, mt)

    Egamma = np.linspace(Etr + 10, v_max, nbpoint)
    Ef = np.linspace(10, v_max, nbpoint)

    CosTheta = fillTab(Egamma, Ef, Ktr, mt)
    CosTheta = np.array(CosTheta)
    CosTheta = np.transpose(CosTheta) # It's quicker than reverse because we'd have to calculate Ebar more often

    # Heatmap
    plt.imshow(CosTheta, cmap='viridis', vmin=0.5, vmax = 1, extent=[Egamma.min(), Egamma.max(), Ef.min(), Ef.max()], aspect='auto', origin='lower')

    plt.colorbar(label='Cos($\\theta$)')

    plt.xlabel('$E_\\gamma$')
    plt.ylabel('$E_{\mathrm{fermion}}$')
    plt.title('Cos($\\theta$) for Etr={}GeV'.format(Etr))

    filename = "{}/heatMapCosTheta_{}_{}".format(name, Etr, v_max)
    path="Plots/post/theta/"
    plt.savefig(path+filename+".png")
    plt.savefig(path+filename+".pdf")
    plt.clf()

def plotTheta(Etr_min, Etr_max, Egamma, m, name):
    ETR = np.linspace(Etr_min, Etr_max, 11)

    cmap = plt.get_cmap('rainbow')

    for Etr in ETR:
        k=calcKappa(Etr, m)
        Ebar= calcEbar(Egamma, m, k)
        Ef = np.linspace(1/2*(Egamma - Ebar), 1/2*(Egamma + Ebar), 1000)
        CosTheta = np.array([calcCosTheta(ef, Egamma, k, m) for ef in Ef])
        CosTheta = np.clip(CosTheta, -1, 1) # Precision error. Sometime we have 1.000000000002
        Theta = np.arccos(CosTheta)*180/np.pi
        color = cmap((Etr - Etr_min) / (Etr_max - Etr_min))
        plt.plot(Ef, Theta, label="$Etr$={}GeV".format(Etr), color=color)
    #plt.legend(loc='upper right', bbox_to_anchor=(1.1, 1.1))
    plt.xlabel('$E_{\mathrm{top}}$ (GeV)')
    plt.ylabel('$\\theta$ (°)')
    plt.title('Evolution of cos($\\theta$) for $E_\\gamma$={}'.format(Egamma))

    norm = Normalize(vmin=Etr_min, vmax=Etr_max)
    scalar_map = ScalarMappable(cmap=cmap, norm=norm)
    scalar_map.set_array([])

    ax = plt.gca()
    plt.colorbar(scalar_map, label='Etr', ax=ax)
    plt.savefig("Plots/post/theta/{}/theta_{}.png".format(name, Egamma))
    plt.savefig("Plots/post/theta/{}/theta_{}.pdf".format(name, Egamma))
    plt.clf()

def pGamma(K,E):
    return E*sqrt((1+K)/(1-K))

def PlotInvariantMass(Etr, m, v_max, name):
    Eg = np.linspace(Etr+10, v_max, 100)
    mass = []
    for eg in Eg:
        pgamma = pGamma(calcKappa(Etr,m), eg)
        mass.append(sqrt(eg**2 - pgamma**2))
    plt.plot(Eg,mass, label='Invariant mass')
    plt.xlabel('$E_\\gamma$ (GeV)')
    plt.ylabel('$m_{\mathrm{inv}}$')
    plt.title('Evolution of m for Etr = {}'.format(Etr))
    plt.savefig('Plots/post/mass/{}/mass_{}.png'.format(name, Etr))
    plt.clf()

def partialwidth(x, m, k, e):
    par = [m, k, 1, 1, 1/137, e, 1]
    me2c4 = par[0]*par[0]*par[3]*par[3]*par[3]*par[3]
    num = par[4]*((1-par[1])*(2*par[1]*x*(par[5]-x)+(1+par[1])*me2c4)-par[1]*par[5]*par[5])
    den = par[2]*(1+par[1])*(1+par[1])*(sqrt(1-par[1]*par[1]))*par[5]*par[5]
    return num/(den*par[6])

def Gamma(m, k, e):
    Ebar = calcEbar(e, m, k)
    a = 1/2*(e-Ebar)
    b = 1/2*(e+Ebar)
    numSteps = 1000
    h = (b - a) / numSteps
    integral = 0.0

    for i in range(numSteps):
        x1 = a + i * h
        x2 = a + (i + 1) * h
        fx1 = partialwidth(x1, m, k, e)
        fx2 = partialwidth(x2, m, k, e)
        integral += (fx1 + fx2) * h / 2.0

    return integral

def PlotProba(m, Etr, name):
    Eg = np.linspace(Etr+10, Etr+110, 3)
    if m < 1:
        X = np.linspace(1, 6.4*Etr-7800, 51)
    else:
        X = np.linspace(1e-9, 11e-9, 51)
    for eg in Eg:
        k = calcKappa(Etr, m)
        gamma = Gamma(m, k, eg)
        P=[exp(-gamma*x*5.07e9) for x in X]
        plt.plot(X, P, label="Probability for $E_\\gamma$ = {}".format(eg))

    plt.xlabel('Distance ($\\mu$ m)')
    plt.ylabel('Probability')
    plt.title('Probability of disintegration for Etr = {}'.format(Etr))
    plt.legend()
    plt.savefig('Plots/post/proba/{}/proba_{}.png'.format(name, Etr))
    plt.clf()
    


# Constant for top
mt = 172.5
me = 0.000511
mass = {"el" : me, "top" : mt}

Etr = np.linspace(1500, 9000, 31)
Egamma = np.linspace(1610, 2610, 11)
v_max=10000
for name, m in mass.items():
    for etr in Etr:
        #makeHeatMap(etr, m, v_max, name)
        PlotProba(m, etr, name)

        
    #for eg in Egamma:
        #plotTheta(1500, eg-10, eg, m, name)