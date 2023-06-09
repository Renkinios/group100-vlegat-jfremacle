
/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#ifndef _FEM_H_
#define _FEM_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "gmshc.h"


#define ErrorScan(a)   femErrorScan(a,__LINE__,__FILE__)
#define ErrorGmsh(a)   femErrorGmsh(a,__LINE__,__FILE__)
#define Error(a)       femError(a,__LINE__,__FILE__)
#define Warning(a)     femWarning(a,  __LINE__, __FILE__)
#define FALSE 0 
#define TRUE  1
#define MAXNAME 256

typedef enum {FEM_TRIANGLE,FEM_QUAD} femElementType;
typedef enum {DIRICHLET_X,DIRICHLET_Y,DIRICHLET_N,DIRICHLET_T,
              NEUMANN_X,NEUMANN_Y,NEUMANN_N,NEUMANN_T} femBoundaryType;
typedef enum {PLANAR_STRESS,PLANAR_STRAIN,AXISYM} femElasticCase;
typedef enum {FEM_FULL,FEM_BAND,FEM_ITER,FEM_Cholesky,FEM_Cholesky_band} femSolverType;
typedef enum {FEM_NO,FEM_XNUM,FEM_YNUM} femRenumType;
static const double _gaussDos2Eta[2]     = { -0.577350269189626, 0.577350269189626};
static const double _gaussDos2Weight[2]  = { 1, 1};
typedef struct {
    int nNodes;
    double *X;
    double *Y;
    int*number;
} femNodes;

typedef struct {
    int nLocalNode;
    int nElem;
    int *elem;
    femNodes *nodes;
} femMesh;

typedef struct {
    femMesh *mesh;
    int nElem;
    int *elem;
    char name[MAXNAME];
} femDomain;

typedef struct {
    double LxPlate, LyPlate;
    double h;
    femElementType elementType;
    double (*geoSize)(double x, double y);
    femNodes *theNodes;
    femMesh  *theElements;
    femMesh  *theEdges;
    int nDomains;
    femDomain **theDomains;
} femGeo;

typedef struct {
    int n;
    void (*x2)(double *xsi, double *eta);
    void (*phi2)(double xsi, double eta, double *phi);
    void (*dphi2dx)(double xsi, double eta, double *dphidxsi, double *dphideta);
} femDiscrete;
    
typedef struct {
    int n;
    const double *xsi;
    const double *eta;
    const double *weight;
} femIntegration;

typedef struct {
    double *B;
    double **A;
    int size;
} femFullSystem;

typedef struct {
    double *B;
    double **A;        
    int size;
    int band;        
} femBandSystem;

typedef struct {
    double *R;
    double *D;
    double *S;
    double *X; 
    double error;      
    int size;
    int iter;        
} femIterativeSolver;

typedef struct {
    femDomain* domain;
    femBoundaryType type; 
    double value;
} femBoundaryCondition;

typedef struct {
    femSolverType type;
    femFullSystem *local;
    void *solver;
} femSolver;

typedef struct {
    double E,nu,rho,g;
    double A,B,C;
    int planarStrainStress;
    int nBoundaryConditions;
    femBoundaryCondition **conditions;  
    int *constrainedNodes; 
    int *contrainteEdges ;
    femGeo *geometry;
    femDiscrete *space;
    femIntegration *rule;
    femSolver *solver;
} femProblem;


void                geoMeshGenerateGeo() ; 
void                geoInitialize();
femGeo*             geoGetGeometry();
double              geoSize(double x, double y);
double              geoSizeDefault(double x, double y);
void                geoSetSizeCallback(double (*geoSize)(double x, double y));
void                geoMeshGenerate();
void                geoMeshImport();
void                geoMeshPrint();
void                geoMeshWrite(const char *filename);
void                geoMeshRead(const char *filename);
void                geoSetDomainName(int iDomain, char *name);
int                 geoGetDomain(char *name);
void                geoFinalize();

femProblem *          femElasticityCreate(femGeo* theGeometry, 
                            double E, double nu, double rho, double g, femElasticCase iCase, femSolverType type,femRenumType renumType);
void                 femElasticityFree(femProblem *theProblem);
void                 femElasticityPrint(femProblem *theProblem);
void                 femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value);
double*              femElasticitySolve(femProblem *theProblem);

femIntegration*      femIntegrationCreate(int n, femElementType type);
void                 femIntegrationFree(femIntegration *theRule);

femDiscrete*         femDiscreteCreate(int n, femElementType type);
void                 femDiscreteFree(femDiscrete* mySpace);
void                 femDiscretePrint(femDiscrete* mySpace);
void                 femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta);
void                 femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi);
void                 femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta);

femSolver*          femSolverCreate(int sizeLoc);
void                femSolverAssemble(femSolver* mySolver, double *Aloc, double *Bloc, double *Uloc,int *map, int nLoc,int*renumber);
femSolver *         femSolverFullCreate(int size);
femSolver*          femSolverBandCreate(int size, int sizeLoc, int band);
void                femSolverFree(femSolver *mySolver);
void                femSolverInit(femSolver *mySolver);
double              femSolverGet(femSolver *mySolver,int i,int j);
void                femSolverPrint(femSolver *mySolver);
double *            femSolverEliminate(femSolver *mySolver);

femFullSystem*       femFullSystemCreate(int size);
void                 femFullSystemFree(femFullSystem* mySystem);
void                 femFullSystemInit(femFullSystem* mySystem);
void                 femFullSystemPrint(femFullSystem* mySystem);
void                 femFullSystemPrintInfos(femFullSystem* mySystem);
double*              femFullSystemEliminate(femFullSystem* mySystem);
void                 femFullSystemAssemble(femFullSystem* mySystem, double *Aloc, double *Bloc, int *map, int nLoc);
double               femFullSystemGet(femFullSystem* mySystem, int i, int j);

femBandSystem*       femBandSystemCreate(int size, int band);
void                 femBandSystemFree(femBandSystem* myBandSystem);
void                 femBandSystemInit(femBandSystem *myBand);
void                 femBandSystemPrint(femBandSystem *myBand);
void                 femBandSystemPrintInfos(femBandSystem *myBand);
double*              femBandSystemEliminate(femBandSystem *myBand);
void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc,int*renumber);
double               femBandSystemGet(femBandSystem* myBandSystem, int i, int j);

femIterativeSolver*  femIterativeSolverCreate(int size);
void                 femIterativeSolverFree(femIterativeSolver* mySolver);
void                 femIterativeSolverInit(femIterativeSolver* mySolver);
void                 femIterativeSolverPrint(femIterativeSolver* mySolver);
void                 femIterativeSolverPrintInfos(femIterativeSolver* mySolver);
double*              femIterativeSolverEliminate(femIterativeSolver* mySolver);
void                 femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc);
double               femIterativeSolverGet(femIterativeSolver* mySolver, int i, int j);
int                  femIterativeSolverConverged(femIterativeSolver *mySolver);

double              femMin(double *x, int n);
double              femMax(double *x, int n);
void                femError(char *text, int line, char *file);
void                femErrorScan(int test, int line, char *file);
void                femErrorGmsh(int test, int line, char *file);
void                femWarning(char *text, int line, char *file);
void                getEdge(femProblem *problem,int iEdge,double *jac,double *nx,double *ny,int *map) ;
void                femMeshRenumber(femMesh *theMesh, femRenumType renumType);

double          **choleskyDecomposition(double **A, int n);
double          *solvecholesky(double **L, double *b, int n);

void femSystemConstrainNEUMANN(femProblem *theProblem,femFullSystem *mySystem, 
                             int iEdge, double value,femBoundaryType type,femElasticCase iCase) ;
void  femSystemConstrainDIRICHLETXY(femFullSystem *mySystem, 
                             int myNode, double myValue,femBoundaryType myType) ;
double *            reorder(double *B,femNodes*Nodes,int* renumber);
#endif