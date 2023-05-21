Ce projet est composé de trois fichiers principaux :

1. **ProjectPreProcessor** : Ce fichier est responsable de la création du maillage et de la modification des fichiers `problem.txt` et `mesh.txt`. Il génère ces fichiers en fonction des paramètres spécifiés.

2. **Projet** : Ce fichier est chargé de résoudre le problème en utilisant les fichiers `problem.txt` et `mesh.txt` générés par le préprocesseur. Il calcule les solutions et génère les fichiers `U` et `V` contenant les résultats.

3. **ProjetPostProcessor** : Ce fichier lit les fichiers `U` et `V` et crée une visualisation du maillage à partir des données de solution.


Tous sa ce trouve à  https://github.com/Renkinios/Projet_element_finis_elasticity_lineaire.git
