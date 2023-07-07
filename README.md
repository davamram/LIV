# RivetLIV
# Projet d'analyse des photons isolés avec les données du LHC en 2012

Ce projet contient du code pour effectuer une analyse des photons isolés en utilisant les données Atlas. 

## Prérequis

Avant d'exécuter ce code, assurez-vous d'avoir les éléments suivants installés :

- Rivet (version 3.1.7) : [lien vers Rivet](https://rivet.hepforge.org/)
- Les bibliothèques Rivet (Rivet Analysis, Rivet Projections, Rivet Tools, etc.)

## Instructions d'utilisation

1. Clonez ce dépôt GitHub sur votre machine locale en utilisant la commande suivante :
``` git clone https://github.com/davamram/LIV ```
2. Accédez au répertoire du dépôt cloné :
``` cd LIV ```
3. Compiler le code à l'aide des commandes rivet :
``` rivet-build TEST_ANALYSIS.cc lib/SelectPhoton.cc ```
4. Exécutez le code en fournissant les fichiers de données d'entrée nécessaires :
``` rivet --pwd -a TEST_ANALYSIS <my_events>.hepmc.gz ```
5. Tracer les histogrammes à partir des fichiers yoda :
``` rivet-mkhtml --errs --no-weight Rivet.yoda ```
