import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
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

def makeHeatMap(Etr, mt, v_max):

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

    filename = "heatMapCosTheta_{}_{}".format(Etr, v_max)
    path="Plots/post/"
    plt.savefig(path+filename+".png")
    plt.savefig(path+filename+".pdf")
    plt.clf()

def plotTheta(Etr_min, Etr_max, Egamma, m):
    ETR = np.linspace(Etr_min, Etr_max, 11)

    cmap = plt.get_cmap('rainbow')

    for Etr in ETR:
        k=calcKappa(Etr, m)
        Ebar= calcEbar(Egamma, m, k)
        Ef = np.linspace(1/2*(Egamma - Ebar), 1/2*(Egamma + Ebar), 1000)
        CosTheta = np.array([calcCosTheta(ef, Egamma, k, m) for ef in Ef])
        CosTheta2=[]
        for ef in Ef:
            CosTheta2.append(calcCosTheta(ef, Egamma, k, m))
        Theta = np.arccos(CosTheta)*180/np.pi
        color = cmap((Etr - Etr_min) / (Etr_max - Etr_min))
        plt.plot(Ef, Theta, label="$Etr$={}GeV".format(Etr), color=color)
    #plt.legend(loc='upper right', bbox_to_anchor=(1.1, 1.1))
    plt.xlabel('$E_{\mathrm{top}}$ (GeV)')
    plt.ylabel('$\\theta$ (°)')
    plt.title('Evolution of cos($\\theta$) for $E_\\gamma$={}'.format(Egamma))

    norm = Normalize(vmin=Etr_min, vmax=Etr_max)
    scalar_map = ScalarMappable(cmap=cmap, norm=norm)
    scalar_map.set_array([])  # Vous devez définir un tableau vide ici

    # Ajoutez la barre de couleur à votre graphique
    plt.colorbar(scalar_map, label='Etr')
    plt.savefig("Plots/post/theta_{}.png".format(Egamma))
    plt.savefig("Plots/post/theta_{}.pdf".format(Egamma))
    plt.clf()



# Constant for top
mt = 172.5

Etr = np.linspace(1500, 9000, 16)
v_max=10000
for etr in Etr:
    makeHeatMap(etr, mt, v_max)

Egamma = np.linspace(1610, 2610, 11)
for eg in Egamma:
    plotTheta(1500, eg-10, eg, mt)
