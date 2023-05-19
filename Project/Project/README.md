# Instructions d'utilisation - Cholesky et Renumérotation

Ce fichier README fournit des instructions sur la manière d'utiliser les options de renumérotation et de résolution de systèmes linéaires à l'aide de l'algorithme de Cholesky ou d'un solveur plein dans notre application.

## Renumérotation

Pour changer le type de renumérotation, vous pouvez modifier la variable `renumType` dans votre code. Les options disponibles pour `renumType` sont :

- `FEM_NO` : Aucune renumérotation n'est effectuée.
- `FEM_XNUM` : Renumérotation basée sur la numérotation en X.
- `FEM_YNUM` : Renumérotation basée sur la numérotation en Y.

Assurez-vous de choisir la valeur appropriée pour `renumType` en fonction de vos besoins.

## Solveur Cholesky

Pour utiliser le solveur Cholesky, vous pouvez modifier la variable `solverType` dans votre code. Les options disponibles pour `solverType` sont :

- `FEM_FULL` : Résolution du système linéaire en utilisant la méthode du solveur plein.
- `FEM_Cholesky` : Résolution du système linéaire en utilisant l'algorithme de Cholesky.

Veillez à sélectionner la valeur appropriée pour `solverType` en fonction de vos exigences.

## Exemple d'utilisation

Voici un exemple de code montrant comment utiliser les options de renumérotation et de résolution de systèmes linéaires :

```python
# Définition des variables de configuration
renumType = FEM_XNUM
solverType = FEM_Cholesky

# Utilisation des variables dans votre code
# ...
