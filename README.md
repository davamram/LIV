# RivetLIV
# Projet d'analyse des photons isolés avec les données du LHC en 2012

Ce projet contient du code pour effectuer une analyse des photons isolés en utilisant les données Atlas. 

## Prérequis

Avant d'exécuter ce code, assurez-vous d'avoir les éléments suivants installés :

- Rivet (version 3.1.7) : [lien vers Rivet](https://rivet.hepforge.org/)
- Les bibliothèques Rivet (Rivet Analysis, Rivet Projections, Rivet Tools, etc.)
- ROOT 6.26
- Python 3.9

## Instructions d'utilisation

### Analyse des évènements
1. Clonez ce dépôt GitHub sur votre machine locale en utilisant la commande suivante :
``` git clone https://github.com/davamram/LIV ```
2. Accédez au répertoire du dépôt cloné :
``` cd LIV ```
3. Compiler le code à l'aide des commandes rivet :
``` rivet-build TEST_ANALYSIS.cc lib/SelectPhoton.cc ```
4. Exécutez le code en fournissant les fichiers de données d'entrée nécessaires :
``` rivet --pwd -a TEST_ANALYSIS <my_events>.hepmc.gz ```
5. Tracer les histogrammes à partir des fichiers yoda :
``` rivet-mkhtml --errs --no-weight Rivet.yoda -o <PlotPath>/<Etr>GeV```

### Traitement Statistique

1. Aller dans le repertoire Statistic :
```cd Statistic```

2. Changer le filePath dans le fichier Tools/ExtractValues.cpp
```cpp
    std::string filePath = "/home/amram/Documents/LorentzPhotons/Rivet/LIV/Plots/Sherpa/Reweight/" + std::to_string(energy) + "GeV/TEST_ANALYSIS/d01-x01-y01.dat";
```

3. Recompiler la librairie :
```g++ -shared -o biblio.so Tools/ExtractValues.cpp -fPIC `root-config --cflags --glibs` ```

3. Lancer le calcul CLs avec root : 
```root -l -q "cls_multiBin.cpp($value)" >> results/cls_etr_$value.txt```