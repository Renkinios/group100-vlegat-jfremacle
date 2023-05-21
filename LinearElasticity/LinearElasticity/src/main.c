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
#include <time.h>

int main(void)
{  
    printf("\n\n    V : Mesh and displacement norm \n");
    printf("    D : Domains \n");
    printf("    N : Next domain highlighted\n\n\n");


    double Lx = 1.0;
    double Ly = 1.0;
      
    geoInitialize();
    femGeo* theGeometry = geoGetGeometry();
    
    theGeometry->LxPlate     =  Lx;
    theGeometry->LyPlate     =  Ly;     
    theGeometry->h           =  Lx * 0.05;    
    theGeometry->elementType = FEM_QUAD;
  
    // geoMeshGenerate(); 
    // geoMeshImport();
    // geoSetDomainName(0,"Symmetry");
    // geoSetDomainName(7,"Bottom");
    // geoSetDomainName(6,"Top");
    geoMeshGenerateGeo();
    geoMeshImport();
    geoSetDomainName(0,"Contre poid gauche");
    geoSetDomainName(1,"Contre poid haut");
    geoSetDomainName(3,"Saut");
    geoSetDomainName(7,"Contre poid bas");
    geoSetDomainName(6,"Contre poid droit");
//
//  -2- Creation probleme 
//
    double E   = 2.e9;
    double nu  = 0.37;
    double rho = 1.22e3; 
    double g   = 9.81;
    double m = -700; 
    femProblem* theProblem = femElasticityCreate(theGeometry,E,nu,rho,g,PLANAR_STRESS,FEM_FULL,FEM_NO);
    femElasticityAddBoundaryCondition(theProblem,"Contre poid droit",DIRICHLET_X,0.0) ; 
    femElasticityAddBoundaryCondition(theProblem,"Saut",NEUMANN_T, -3.2600000e+6) ;     
    femElasticityAddBoundaryCondition(theProblem,"Contre poid bas",DIRICHLET_Y, 0.0 ) ;
    clock_t tic = clock();
    double *theSoluce = femElasticitySolve(theProblem) ; 
    printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
    
   
//  
//  -3- Deformation du maillage pour le plot final
//      Creation du champ de la norme du deplacement
//
    


    femNodes *theNodes = theGeometry->theNodes;
    double deformationFactor = 1;
    double *normDisplacement = malloc(theNodes->nNodes * sizeof(double));
    
    for (int i=0; i<theNodes->nNodes; i++){
        theNodes->X[i] += theSoluce[2*i+0]*deformationFactor;
        theNodes->Y[i] += theSoluce[2*i+1]*deformationFactor;
        normDisplacement[i] = sqrt(theSoluce[2*i+0]*theSoluce[2*i+0] + 
                                   theSoluce[2*i+1]*theSoluce[2*i+1]); }
    int n = theGeometry->theNodes->nNodes;
    printf("---------------------------------------------------\n") ;
    for (int i=0; i<n; i++){
        if(theNodes->Y[i] >0.8 && theNodes->X[i] > 0.8)
            printf(" Y[i] , X[i]: %f, %f \n",theNodes->Y[i],theNodes->X[i]);
    }
    // double hMin = femMin(normDisplacement,theNodes->nNodes);  
    // double hMax = femMax(normDisplacement,theNodes->nNodes);  
    // printf(" ==== Minimum displacement          : %14.7e \n",hMin);
    // printf(" ==== Maximum displacement          : %14.7e \n",hMax);
 
//
//  -4- Visualisation du maillage
//  
    
    int mode = 1; 
    int domain = 0;
    int freezingButton = FALSE;
    double t, told = 0;
    char theMessage[MAXNAME];
   
 
    GLFWwindow* window = glfemInit("EPL1110 : Linear elasticity ");
    glfwMakeContextCurrent(window);

    do {
        int w,h;
        glfwGetFramebufferSize(window,&w,&h);
        glfemReshapeWindows(theGeometry->theNodes,w,h);

        t = glfwGetTime();  
        if (glfwGetKey(window,'D') == GLFW_PRESS) { mode = 0;}
        if (glfwGetKey(window,'V') == GLFW_PRESS) { mode = 1;}
        if (glfwGetKey(window,'N') == GLFW_PRESS && freezingButton == FALSE) { domain++; freezingButton = TRUE; told = t;}
        
        if (t-told > 0.5) {freezingButton = FALSE; }
        if (mode == 1) {
            glfemPlotField(theGeometry->theElements,normDisplacement);
            glfemPlotMesh(theGeometry->theElements); 
            sprintf(theMessage, "Number of elements : %d ",theGeometry->theElements->nElem);
            glColor3f(1.0,0.0,0.0); glfemMessage(theMessage); }
        if (mode == 0) {
            domain = domain % theGeometry->nDomains;
            glfemPlotDomain( theGeometry->theDomains[domain]); 
            sprintf(theMessage, "%s : %d ",theGeometry->theDomains[domain]->name,domain);
             glColor3f(1.0,0.0,0.0); glfemMessage(theMessage);  }
            
         glfwSwapBuffers(window);
         glfwPollEvents();
    } while( glfwGetKey(window,GLFW_KEY_ESCAPE) != GLFW_PRESS &&
             glfwWindowShouldClose(window) != 1 );
            
    // Check if the ESC key was pressed or the window was closed

    free(normDisplacement);  
    geoFinalize();
    glfwTerminate(); 
    
    exit(EXIT_SUCCESS);
    return 0;  
}

 