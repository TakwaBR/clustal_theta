# Projet court - ALIGNEMENT MULTIPLE HEURISTIQUE PAR LA METHODE CLUSTAL

## Introduction
Ce projet vise à effectuer un alignement multiple de séquences en utilisant l'algorithme de Needlman et Wunsch pour un premier alignement global par pair. Ensuite l'algortihme UPGMA est utilisé pour regrouper les séquences les plus proches. Finalement l'alignement multiple est fait (le code n'est pas complet pour cette étape).

## Installation
Pour utiliser ce projet il faut tout d'abord installer [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

Ensuite il faut cloner ce répertoire
```
git clone https://github.com/TakwaBR/clustal_theta.git
cd clustal_theta
```
Créer l'environnement nécessaire
```
conda env create -f projet_clustal.yml
```
Activer l'environnement
```
conda activate projet_clustal
```
## Utilisation
Pour faire tourner le code, il faut lancer le script pair alignment.py
```
cd src/
python3 pair_alignment.py
```
Deux autres scripts se trouvent dans ce dossier, le script upgma.py et le script multiple_alignment qui n'est pas complet.\
Dans le dossier data/ se trouvent les fichiers contenant 10 séquences de la protéine p53 au format fasta et la matrice de substitution BLOSUM62 au format csv
