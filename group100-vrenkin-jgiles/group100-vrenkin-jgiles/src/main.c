/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Code de calcul
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */
 
#include "fem.h"
#include <time.h>

int main(void)
{  
    femGeo* theGeometry = geoGetGeometry();   
    geoMeshRead("../data/mesh.txt");
    femRenumType renumType = FEM_NO  ;
    femSolverType solverType = FEM_Cholesky;
    femProblem* theProblem = femElasticityRead(theGeometry,"../data/problem.txt",solverType,renumType);
    femElasticityPrint(theProblem);
    clock_t tic = clock();

    double *theSoluce = femElasticitySolve(theProblem); 
    printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
    
    femNodes *theNodes = theGeometry->theNodes;
    femFieldWrite(theNodes->nNodes,2,&theSoluce[0],"../data/U.txt");
    femFieldWrite(theNodes->nNodes,2,&theSoluce[1],"../data/V.txt");
    femElasticityFree(theProblem); 
    geoFree();
    return 0;  
}
 
