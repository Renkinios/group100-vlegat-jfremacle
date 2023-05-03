/*
 *  main.c
 *  Library for EPL1110 : Finite Elements for dummies
 *  Elasticite lineaire plane
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "glfem.h"


int main(void)
{  
    
    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();
    geoMeshRead("../data/mesh.txt"); //ici pour la renumératasion des élements et des noeuds with     femMeshRenumber()
    femProblem* theProblem = femElasticityCreate(theGeometry,"../data/problem.txt");
    femElasticityPrint(theProblem);
    return 0;  
    // pas oublié de free le number 
}

 
