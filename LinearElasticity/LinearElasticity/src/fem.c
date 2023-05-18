/*
 *  fem.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

femGeo theGeometry;

femGeo *geoGetGeometry()                        { return &theGeometry; }

double geoSizeDefault(double x, double y)       { return theGeometry.h; }

double geoGmshSize(int dim, int tag, double x, double y, double z, double lc, void *data)
                                                { return theGeometry.geoSize(x,y);    }
void geoInitialize() 
{
    int ierr;
    theGeometry.geoSize = geoSizeDefault;
    gmshInitialize(0,NULL,1,0,&ierr);                         ErrorGmsh(ierr);
    gmshModelAdd("MyGeometry",&ierr);                         ErrorGmsh(ierr);
    gmshModelMeshSetSizeCallback(geoGmshSize,NULL,&ierr);     ErrorGmsh(ierr);
    theGeometry.theNodes = NULL;
    theGeometry.theElements = NULL;
    theGeometry.theEdges = NULL;
    theGeometry.nDomains = 0;
    theGeometry.theDomains = NULL;
}

void geoFinalize() 
{
    int ierr;
    
    if (theGeometry.theNodes) {
        free(theGeometry.theNodes->X);
        free(theGeometry.theNodes->Y);
        free(theGeometry.theNodes); }
    if (theGeometry.theElements) {
        free(theGeometry.theElements->elem);
        free(theGeometry.theElements); }
    if (theGeometry.theEdges) {
        free(theGeometry.theEdges->elem);
        free(theGeometry.theEdges); }
    for (int i=0; i < theGeometry.nDomains; i++) {
        free(theGeometry.theDomains[i]->elem);
        free(theGeometry.theDomains[i]);  }
    free(theGeometry.theDomains);
    gmshFinalize(&ierr); ErrorGmsh(ierr);
}

void geoSetSizeCallback(double (*geoSize)(double x, double y)) 
{
    theGeometry.geoSize = geoSize; }


void geoMeshImport() 
{
    int ierr;
    
    /* Importing nodes */
    
    size_t nNode,n,m,*node;
    double *xyz,*trash;
    gmshModelMeshRenumberNodes(&ierr);                        ErrorGmsh(ierr);
    gmshModelMeshGetNodes(&node,&nNode,&xyz,&n,
                         &trash,&m,-1,-1,0,0,&ierr);          ErrorGmsh(ierr);                         
    femNodes *theNodes = malloc(sizeof(femNodes));
    theNodes->nNodes = nNode;
    theNodes->X = malloc(sizeof(double)*(theNodes->nNodes));
    theNodes->Y = malloc(sizeof(double)*(theNodes->nNodes));
    for (int i = 0; i < theNodes->nNodes; i++){
        theNodes->X[i] = xyz[3*node[i]-3];
        theNodes->Y[i] = xyz[3*node[i]-2]; }
    theGeometry.theNodes = theNodes;
    theNodes->number = malloc(sizeof(int)*theNodes->nNodes); 
    for (int i = 0; i < theNodes->nNodes; i++) 
          theNodes->number[i] = i; 
    gmshFree(node);
    gmshFree(xyz);
    gmshFree(trash);
    printf("Geo     : Importing %d nodes \n",theGeometry.theNodes->nNodes);
       
    /* Importing elements */
    /* Pas super joli : a ameliorer pour eviter la triple copie */
        
    size_t nElem, *elem;
    gmshModelMeshGetElementsByType(1,&elem,&nElem,
                               &node,&nNode,-1,0,1,&ierr);    ErrorGmsh(ierr);
    femMesh *theEdges = malloc(sizeof(femMesh));
    theEdges->nLocalNode = 2;
    theEdges->nodes = theNodes;
    theEdges->nElem = nElem;  
    theEdges->elem = malloc(sizeof(int)*2*theEdges->nElem);
    for (int i = 0; i < theEdges->nElem; i++)
        for (int j = 0; j < theEdges->nLocalNode; j++)
            theEdges->elem[2*i+j] = node[2*i+j]-1;  
    theGeometry.theEdges = theEdges;
    int shiftEdges = elem[0];
    gmshFree(node);
    gmshFree(elem);
    printf("Geo     : Importing %d edges \n",theEdges->nElem);
  
    gmshModelMeshGetElementsByType(2,&elem,&nElem,
                               &node,&nNode,-1,0,1,&ierr);    ErrorGmsh(ierr);
    if (nElem != 0) {
      femMesh *theElements = malloc(sizeof(femMesh));
      theElements->nLocalNode = 3;
      theElements->nodes = theNodes;
      theElements->nElem = nElem;  
      theElements->elem = malloc(sizeof(int)*3*theElements->nElem);
      for (int i = 0; i < theElements->nElem; i++)
          for (int j = 0; j < theElements->nLocalNode; j++)
              theElements->elem[3*i+j] = node[3*i+j]-1;  
      theGeometry.theElements = theElements;
      gmshFree(node);
      gmshFree(elem);
      printf("Geo     : Importing %d triangles \n",theElements->nElem); }
    
    int nElemTriangles = nElem;
    gmshModelMeshGetElementsByType(3,&elem,&nElem,
                               &node,&nNode,-1,0,1,&ierr);    ErrorGmsh(ierr);
    if (nElem != 0 && nElemTriangles != 0)  
      Error("Cannot consider hybrid geometry with triangles and quads :-(");                       
                               
    if (nElem != 0) {
      femMesh *theElements = malloc(sizeof(femMesh));
      theElements->nLocalNode = 4;
      theElements->nodes = theNodes;
      theElements->nElem = nElem;  
      theElements->elem = malloc(sizeof(int)*4*theElements->nElem);
      for (int i = 0; i < theElements->nElem; i++)
          for (int j = 0; j < theElements->nLocalNode; j++)
              theElements->elem[4*i+j] = node[4*i+j]-1;  
      theGeometry.theElements = theElements;
      gmshFree(node);
      gmshFree(elem);
      printf("Geo     : Importing %d quads \n",theElements->nElem); }

    
    /* Importing 1D entities */
  
    int *dimTags;
    gmshModelGetEntities(&dimTags,&n,1,&ierr);        ErrorGmsh(ierr);
    theGeometry.nDomains = n/2;
    theGeometry.theDomains = malloc(sizeof(femDomain*)*n/2);
    printf("Geo     : Importing %d entities \n",theGeometry.nDomains);

    for (int i=0; i < n/2; i++) {
        int dim = dimTags[2*i+0];
        int tag = dimTags[2*i+1];
        femDomain *theDomain = malloc(sizeof(femDomain)); 
        theGeometry.theDomains[i] = theDomain;
        theDomain->mesh = theEdges;
        sprintf(theDomain->name, "Entity %d ",tag-1);
        int *elementType;
        size_t nElementType, **elementTags, *nElementTags, nnElementTags, **nodesTags, *nNodesTags, nnNodesTags; 
        gmshModelMeshGetElements(&elementType, &nElementType, &elementTags, &nElementTags, &nnElementTags, &nodesTags, &nNodesTags, &nnNodesTags, dim, tag, &ierr);
        theDomain->nElem = nElementTags[0];
        theDomain->elem = malloc(sizeof(int)*2*theDomain->nElem); 
        for (int j = 0; j < theDomain->nElem; j++) {
            theDomain->elem[j] = elementTags[0][j] - shiftEdges; }
        printf("Geo     : Entity %d : %d elements \n",i,theDomain->nElem);
        gmshFree(nElementTags);
        gmshFree(nNodesTags);
        gmshFree(elementTags);
        gmshFree(nodesTags);
        gmshFree(elementType); }
    gmshFree(dimTags);
 
    return;

}

void geoMeshPrint() 
{
   femNodes *theNodes = theGeometry.theNodes;
   if (theNodes != NULL) {
      printf("Number of nodes %d \n", theNodes->nNodes);
      for (int i = 0; i < theNodes->nNodes; i++) {
        printf("%6d : %14.7e %14.7e \n",i,theNodes->X[i],theNodes->Y[i]); }}
   femMesh *theEdges = theGeometry.theEdges;
   if (theEdges != NULL) {
     printf("Number of edges %d \n", theEdges->nElem);
     int *elem = theEdges->elem;
     for (int i = 0; i < theEdges->nElem; i++) {
        printf("%6d : %6d %6d \n",i,elem[2*i],elem[2*i+1]); }}
   femMesh *theElements = theGeometry.theElements;
   if (theElements != NULL) {
     if (theElements->nLocalNode == 3) {
        printf("Number of triangles %d \n", theElements->nElem);
        int *elem = theElements->elem;
        for (int i = 0; i < theElements->nElem; i++) {
            printf("%6d : %6d %6d %6d\n",i,elem[3*i],elem[3*i+1],elem[3*i+2]); }}
     if (theElements->nLocalNode == 4) {
        printf("Number of quads %d \n", theElements->nElem);
        int *elem = theElements->elem;
        for (int i = 0; i < theElements->nElem; i++) {
            printf("%6d : %6d %6d %6d %6d\n",i,elem[4*i],elem[4*i+1],elem[4*i+2],elem[4*i+3]); }}}
   int nDomains = theGeometry.nDomains;
   printf("Number of domains %d\n", nDomains);
   for (int iDomain = 0; iDomain < nDomains; iDomain++) {
      femDomain *theDomain = theGeometry.theDomains[iDomain];
      printf("  Domain : %6d \n", iDomain);
      printf("  Name : %s\n", theDomain->name);
      printf("  Number of elements : %6d\n", theDomain->nElem);
      for (int i=0; i < theDomain->nElem; i++){
 //         if (i != theDomain->nElem  && (i % 10) != 0)  printf(" - ");
          printf("%6d",theDomain->elem[i]);
          if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) printf("\n"); }
      printf("\n"); }
  
  
}


void geoMeshWrite(const char *filename) 
{
   FILE* file = fopen(filename,"w");
 
   femNodes *theNodes = theGeometry.theNodes;
   fprintf(file, "Number of nodes %d \n", theNodes->nNodes);
   for (int i = 0; i < theNodes->nNodes; i++) {
      fprintf(file,"%6d : %14.7e %14.7e \n",i,theNodes->X[i],theNodes->Y[i]); }
      
   femMesh *theEdges = theGeometry.theEdges;
   fprintf(file,"Number of edges %d \n", theEdges->nElem);
   int *elem = theEdges->elem;
   for (int i = 0; i < theEdges->nElem; i++) {
      fprintf(file,"%6d : %6d %6d \n",i,elem[2*i],elem[2*i+1]); }
      
   femMesh *theElements = theGeometry.theElements;
   if (theElements->nLocalNode == 3) {
      fprintf(file,"Number of triangles %d \n", theElements->nElem);
      elem = theElements->elem;
      for (int i = 0; i < theElements->nElem; i++) {
          fprintf(file,"%6d : %6d %6d %6d\n",i,elem[3*i],elem[3*i+1],elem[3*i+2]); }}
   if (theElements->nLocalNode == 4) {
      fprintf(file,"Number of quads %d \n", theElements->nElem);
      elem = theElements->elem;
      for (int i = 0; i < theElements->nElem; i++) {
          fprintf(file,"%6d : %6d %6d %6d %6d\n",i,elem[4*i],elem[4*i+1],elem[4*i+2],elem[4*i+3]); }}
     
   int nDomains = theGeometry.nDomains;
   fprintf(file,"Number of domains %d\n", nDomains);
   for (int iDomain = 0; iDomain < nDomains; iDomain++) {
      femDomain *theDomain = theGeometry.theDomains[iDomain];
      fprintf(file,"  Domain : %6d \n", iDomain);
      fprintf(file,"  Name : %s\n", theDomain->name);
      fprintf(file,"  Number of elements : %6d\n", theDomain->nElem);
      for (int i=0; i < theDomain->nElem; i++){
          fprintf(file,"%6d",theDomain->elem[i]);
          if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) fprintf(file,"\n"); }
      fprintf(file,"\n"); }
    
   fclose(file);
}

void geoMeshRead(const char *filename) 
{
   FILE* file = fopen(filename,"r");
   
   int trash, *elem;
   
   femNodes *theNodes = malloc(sizeof(femNodes));
   theGeometry.theNodes = theNodes;
   ErrorScan(fscanf(file, "Number of nodes %d \n", &theNodes->nNodes));
   theNodes->X = malloc(sizeof(double)*(theNodes->nNodes));
   theNodes->Y = malloc(sizeof(double)*(theNodes->nNodes));
   for (int i = 0; i < theNodes->nNodes; i++) {
       ErrorScan(fscanf(file,"%d : %le %le \n",&trash,&theNodes->X[i],&theNodes->Y[i]));} 

   femMesh *theEdges = malloc(sizeof(femMesh));
   theGeometry.theEdges = theEdges;
   theEdges->nLocalNode = 2;
   theEdges->nodes = theNodes;
   ErrorScan(fscanf(file, "Number of edges %d \n", &theEdges->nElem));
   theEdges->elem = malloc(sizeof(int)*theEdges->nLocalNode*theEdges->nElem);
   for(int i=0; i < theEdges->nElem; ++i) {
        elem = theEdges->elem;
        ErrorScan(fscanf(file, "%6d : %6d %6d \n", &trash,&elem[2*i],&elem[2*i+1])); }
  
   femMesh *theElements = malloc(sizeof(femMesh));
   theGeometry.theElements = theElements;
   theElements->nLocalNode = 0;
   theElements->nodes = theNodes;
   char elementType[MAXNAME];  
   ErrorScan(fscanf(file, "Number of %s %d \n",elementType,&theElements->nElem));  
   if (strncasecmp(elementType,"triangles",MAXNAME) == 0) {
      theElements->nLocalNode = 3;
      theElements->elem = malloc(sizeof(int)*theElements->nLocalNode*theElements->nElem);
      for(int i=0; i < theElements->nElem; ++i) {
          elem = theElements->elem;
          ErrorScan(fscanf(file, "%6d : %6d %6d %6d \n", 
                    &trash,&elem[3*i],&elem[3*i+1],&elem[3*i+2])); }}
   if (strncasecmp(elementType,"quads",MAXNAME) == 0) {
      theElements->nLocalNode = 4;
      theElements->elem = malloc(sizeof(int)*theElements->nLocalNode*theElements->nElem);
      for(int i=0; i < theElements->nElem; ++i) {
          elem = theElements->elem;
          ErrorScan(fscanf(file, "%6d : %6d %6d %6d %6d \n", 
                    &trash,&elem[4*i],&elem[4*i+1],&elem[4*i+2],&elem[4*i+3])); }}
           
   ErrorScan(fscanf(file, "Number of domains %d\n", &theGeometry.nDomains));
   int nDomains = theGeometry.nDomains;
   theGeometry.theDomains = malloc(sizeof(femDomain*)*nDomains);
   for (int iDomain = 0; iDomain < nDomains; iDomain++) {
      femDomain *theDomain = malloc(sizeof(femDomain)); 
      theGeometry.theDomains[iDomain] = theDomain;
      theDomain->mesh = theEdges; 
      ErrorScan(fscanf(file,"  Domain : %6d \n", &trash));
      ErrorScan(fscanf(file,"  Name : %[^\n]s \n", (char*)&theDomain->name));
      ErrorScan(fscanf(file,"  Number of elements : %6d\n", &theDomain->nElem));
      theDomain->elem = malloc(sizeof(int)*2*theDomain->nElem); 
      for (int i=0; i < theDomain->nElem; i++){
          ErrorScan(fscanf(file,"%6d",&theDomain->elem[i]));
          if ((i+1) != theDomain->nElem  && (i+1) % 10 == 0) ErrorScan(fscanf(file,"\n")); }}
    
   fclose(file);
}

void geoSetDomainName(int iDomain, char *name) 
{
    if (iDomain >= theGeometry.nDomains)  Error("Illegal domain number");
    if (geoGetDomain(name) != -1)         Error("Cannot use the same name for two domains");
    sprintf(theGeometry.theDomains[iDomain]->name,"%s",name);
} 

int geoGetDomain(char *name)
{
    int theIndex = -1;
    int nDomains = theGeometry.nDomains;
    for (int iDomain = 0; iDomain < nDomains; iDomain++) {
        femDomain *theDomain = theGeometry.theDomains[iDomain];
        if (strncasecmp(name,theDomain->name,MAXNAME) == 0)
            theIndex = iDomain;  }
    return theIndex;
            
}

static const double _gaussQuad4Xsi[4]    = {-0.577350269189626,-0.577350269189626, 0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Eta[4]    = { 0.577350269189626,-0.577350269189626,-0.577350269189626, 0.577350269189626};
static const double _gaussQuad4Weight[4] = { 1.000000000000000, 1.000000000000000, 1.000000000000000, 1.000000000000000};
static const double _gaussTri3Xsi[3]     = { 0.166666666666667, 0.666666666666667, 0.166666666666667};
static const double _gaussTri3Eta[3]     = { 0.166666666666667, 0.166666666666667, 0.666666666666667};
static const double _gaussTri3Weight[3]  = { 0.166666666666667, 0.166666666666667, 0.166666666666667};

femIntegration *femIntegrationCreate(int n, femElementType type)
{
    femIntegration *theRule = malloc(sizeof(femIntegration));
    if (type == FEM_QUAD && n == 4) {
        theRule->n      = 4;
        theRule->xsi    = _gaussQuad4Xsi;
        theRule->eta    = _gaussQuad4Eta;
        theRule->weight = _gaussQuad4Weight; }
    else if (type == FEM_TRIANGLE && n == 3) {
        theRule->n      = 3;
        theRule->xsi    = _gaussTri3Xsi;
        theRule->eta    = _gaussTri3Eta;
        theRule->weight = _gaussTri3Weight; }
    else Error("Cannot create such an integration rule !");
    return theRule; 
}

void femIntegrationFree(femIntegration *theRule)
{
    free(theRule);
}

void _q1c0_x(double *xsi, double *eta) 
{
    xsi[0] =  1.0;  eta[0] =  1.0;
    xsi[1] = -1.0;  eta[1] =  1.0;
    xsi[2] = -1.0;  eta[2] = -1.0;
    xsi[3] =  1.0;  eta[3] = -1.0;
}

void _q1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = (1.0 + xsi) * (1.0 + eta) / 4.0;  
    phi[1] = (1.0 - xsi) * (1.0 + eta) / 4.0;
    phi[2] = (1.0 - xsi) * (1.0 - eta) / 4.0;
    phi[3] = (1.0 + xsi) * (1.0 - eta) / 4.0;
}

void _q1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] =   (1.0 + eta) / 4.0;  
    dphidxsi[1] = - (1.0 + eta) / 4.0;
    dphidxsi[2] = - (1.0 - eta) / 4.0;
    dphidxsi[3] =   (1.0 - eta) / 4.0;
    dphideta[0] =   (1.0 + xsi) / 4.0;  
    dphideta[1] =   (1.0 - xsi) / 4.0;
    dphideta[2] = - (1.0 - xsi) / 4.0;
    dphideta[3] = - (1.0 + xsi) / 4.0;

}

void _p1c0_x(double *xsi, double *eta) 
{
    xsi[0] =  0.0;  eta[0] =  0.0;
    xsi[1] =  1.0;  eta[1] =  0.0;
    xsi[2] =  0.0;  eta[2] =  1.0;
}

void _p1c0_phi(double xsi, double eta, double *phi)
{
    phi[0] = 1 - xsi - eta;  
    phi[1] = xsi;
    phi[2] = eta;
}

void _p1c0_dphidx(double xsi, double eta, double *dphidxsi, double *dphideta)
{
    dphidxsi[0] = -1.0;  
    dphidxsi[1] =  1.0;
    dphidxsi[2] =  0.0;
    dphideta[0] = -1.0;  
    dphideta[1] =  0.0;
    dphideta[2] =  1.0;

}


femDiscrete *femDiscreteCreate(int n, femElementType type)
{
    femDiscrete *theSpace = malloc(sizeof(femDiscrete));
    if (type == FEM_QUAD && n == 4) {
        theSpace->n       = 4;
        theSpace->x2      = _q1c0_x;
        theSpace->phi2    = _q1c0_phi;
        theSpace->dphi2dx = _q1c0_dphidx; }
    else if (type == FEM_TRIANGLE && n == 3) {
        theSpace->n       = 3;
        theSpace->x2      = _p1c0_x;
        theSpace->phi2    = _p1c0_phi;
        theSpace->dphi2dx = _p1c0_dphidx; }
    else Error("Cannot create such a discrete space !");
    return theSpace; 
}

void femDiscreteFree(femDiscrete *theSpace)
{
    free(theSpace);
}

void femDiscreteXsi2(femDiscrete* mySpace, double *xsi, double *eta)
{
    mySpace->x2(xsi,eta);
}

void femDiscretePhi2(femDiscrete* mySpace, double xsi, double eta, double *phi)
{
    mySpace->phi2(xsi,eta,phi);
}

void femDiscreteDphi2(femDiscrete* mySpace, double xsi, double eta, double *dphidxsi, double *dphideta)
{
    mySpace->dphi2dx(xsi,eta,dphidxsi,dphideta);
}

void femDiscretePrint(femDiscrete *mySpace)
{
    int i,j;
    int n = mySpace->n;
    double xsi[4], eta[4], phi[4], dphidxsi[4], dphideta[4];
    
    femDiscreteXsi2(mySpace,xsi,eta);
    for (i=0; i < n; i++) {
        
        femDiscretePhi2(mySpace,xsi[i],eta[i],phi);
        femDiscreteDphi2(mySpace,xsi[i],eta[i],dphidxsi,dphideta);

        for (j=0; j < n; j++)  {
            printf("(xsi=%+.1f,eta=%+.1f) : ",xsi[i],eta[i]);
            printf(" phi(%d)=%+.1f",j,phi[j]);  
            printf("   dphidxsi(%d)=%+.1f",j,dphidxsi[j]);  
            printf("   dphideta(%d)=%+.1f \n",j,dphideta[j]);  }
        printf(" \n"); }
}
femSolver *femSolverCreate(int sizeLoc)
{
    femSolver *mySolver = malloc(sizeof(femSolver));
    mySolver->local = femFullSystemCreate(sizeLoc);
    return (mySolver);
}
femSolver *femSolverFullCreate(int size) {
    femSolver *mySolver = malloc(sizeof(femSolver)); 
    mySolver ->type = FEM_FULL;
    mySolver->solver = (void *)femFullSystemCreate(size); 
    //printf("system%p\n",mySolver->solver);
    return(mySolver);
}

femSolver *femSolverBandCreate(int size, int sizeLoc, int band)
{
    femSolver *mySolver = femSolverCreate(sizeLoc);
    mySolver->type = FEM_BAND;
    mySolver->solver = (femSolver *)femBandSystemCreate(size,band);
    return(mySolver);
}

femSolver *femSolverIterativeCreate(int size, int sizeLoc)
{
    femSolver *mySolver = femSolverCreate(sizeLoc);
    mySolver->type = FEM_ITER;
    mySolver->solver = (femSolver *)femIterativeSolverCreate(size);
    return(mySolver);
}
void femSolverFree(femSolver *mySolver)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemFree((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemFree((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverFree((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
    femFullSystemFree(mySolver->local);
    free(mySolver);
}
void femSolverInit(femSolver *mySolver)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemInit((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemInit((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverInit((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
}
double femSolverGet(femSolver *mySolver,int i,int j)
{
    double value = 0;
    switch (mySolver->type) {
        case FEM_FULL : value = femFullSystemGet((femFullSystem *)mySolver->solver,i,j); break;
        case FEM_BAND : value = femBandSystemGet((femBandSystem *)mySolver->solver,i,j); break;
        case FEM_ITER : value = (i==j); break;
        default : Error("Unexpected solver type"); }
    return(value);
}

void femSolverPrint(femSolver *mySolver)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemPrint((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemPrint((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverPrint((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
}

void femSolverPrintInfos(femSolver *mySolver)
{
    switch (mySolver->type) {
        case FEM_FULL : femFullSystemPrintInfos((femFullSystem *)mySolver->solver); break;
        case FEM_BAND : femBandSystemPrintInfos((femBandSystem *)mySolver->solver); break;
        case FEM_ITER : femIterativeSolverPrintInfos((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
}

  
void femSolverAssemble(femSolver* mySolver, double *Aloc, double *Bloc, double *Uloc,int *map, int nLoc,int*renumber)
{

    switch (mySolver->type) {
        case FEM_FULL : femFullSystemAssemble((femFullSystem *)mySolver->solver,Aloc,Bloc,map,nLoc); break;
        // case FEM_BAND : femBandSystemAssemble((femBandSystem *)mySolver->solver,Aloc,Bloc,map,nLoc,renumber); break; 
        default : Error("Unexpected solver type"); }
}
    

  
// double *femSolverEliminate(femSolver *mySolver)
// {
//     double *soluce;
//     switch (mySolver->type) {
//         case FEM_FULL : soluce = femFullSystemEliminate((femFullSystem *)mySolver->solver); break;
//         case FEM_BAND : soluce = femBandSystemEliminate((femBandSystem *)mySolver->solver); break;
//         case FEM_ITER : soluce = femIterativeSolverEliminate((femIterativeSolver *)mySolver->solver); break;
//         default : Error("Unexpected solver type"); }
//     return(soluce);
// }



int femSolverConverged(femSolver *mySolver)
{
    int  testConvergence;
    switch (mySolver->type) {
        case FEM_FULL : testConvergence = 1; break;
        case FEM_BAND : testConvergence = 1; break;
        case FEM_ITER : testConvergence = femIterativeSolverConverged((femIterativeSolver *)mySolver->solver); break;
        default : Error("Unexpected solver type"); }
    return(testConvergence);
}


femFullSystem *femFullSystemCreate(int size)
{
    femFullSystem *theSystem = malloc(sizeof(femFullSystem));
    theSystem->A = malloc(sizeof(double*) * size); 
    theSystem->B = malloc(sizeof(double) * size * (size+1));
    theSystem->A[0] = theSystem->B + size;  
    theSystem->size = size;
    int i;
    for (i=1 ; i < size ; i++) 
        theSystem->A[i] = theSystem->A[i-1] + size;
    femFullSystemInit(theSystem);

    return theSystem; 
}

void femFullSystemFree(femFullSystem *theSystem)
{
    free(theSystem->A);
    free(theSystem->B);
    free(theSystem);
}
void femFullSystemInit(femFullSystem *mySystem)
{
    int i,size = mySystem->size;
    for (i=0 ; i < size*(size+1) ; i++) 
        mySystem->B[i] = 0;}
void femFullSystemPrint(femFullSystem *mySystem)
{
    double  **A, *B;
    int     i, j, size;
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    for (i=0; i < size; i++) {
        for (j=0; j < size; j++)
            if (A[i][j] == 0)  printf("         ");   
            else               printf(" %+.1e",A[i][j]);
        printf(" :  %+.1e \n",B[i]); }
}

void femFullSystemPrintInfos(femFullSystem *mySystem)
{
    int  size = mySystem->size;
    printf(" \n");
    printf("    Full Gaussian elimination \n");
    printf("    Storage informations \n");
    printf("    Matrix size      : %8d\n",size);
    printf("    Bytes required   : %8d\n",(int)sizeof(double)*size*(size+1));     
}

double femFullSystemGet(femFullSystem* myFullSystem, int myRow, int myCol)
{
    return(myFullSystem->A[myRow][myCol]); 
}

void femFullSystemAssemble(femFullSystem* mySystem, double *Aloc, double *Bloc, int *map, int nLoc)
{
    int i,j;
    for (i = 0; i < nLoc; i++) { 
        for(j = 0; j < nLoc; j++) {
            mySystem->A[map[i]][map[j]] += Aloc[i*nLoc+j]; }
    mySystem->B[map[i]] += Bloc[i]; }
}
double* femFullSystemEliminate(femFullSystem *mySystem)
{
    double  **A, *B, factor;
    int     i, j, k, size;
    
    A    = mySystem->A;
    B    = mySystem->B;
    size = mySystem->size;
    
    /* Gauss elimination */
    
    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-8 ) {
            printf("Pivot index %d  ",k);
            printf("Pivot value %e  ",A[k][k]);
            Error("Cannot eliminate with such a pivot"); }
        for (i = k+1 ; i <  size; i++) {
            factor = A[i][k] / A[k][k];
            for (j = k+1 ; j < size; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
    
    /* Back-substitution */
    
    for (i = size-1; i >= 0 ; i--) {
        factor = 0;
        for (j = i+1 ; j < size; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }
    
    return(mySystem->B);    
}
femBandSystem *femBandSystemCreate(int size, int band)
{
    femBandSystem *myBandSystem = malloc(sizeof(femBandSystem));
    myBandSystem->B = malloc(sizeof(double)*size*(band+1));
    myBandSystem->A = malloc(sizeof(double*)*size);        
    myBandSystem->size = size;
    myBandSystem->band = band;
    myBandSystem->A[0] = myBandSystem->B + size;
    int i;
    for (i=1 ; i < size ; i++) 
        myBandSystem->A[i] = myBandSystem->A[i-1] + band - 1;
    femBandSystemInit(myBandSystem);
    return(myBandSystem);
}
 
void femBandSystemFree(femBandSystem *myBandSystem)
{
    free(myBandSystem->B);
    free(myBandSystem->A); 
    free(myBandSystem);
}
 
void femBandSystemInit(femBandSystem *myBandSystem)
{
    int i;
    int size = myBandSystem->size;
    int band = myBandSystem->band;
    for (i=0 ; i < size*(band+1) ; i++) 
        myBandSystem->B[i] = 0;        
}
 
void femBandSystemPrint(femBandSystem *myBand)
{
    double  **A, *B;
    int     i, j, band, size;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    for (i=0; i < size; i++) {
        for (j=i; j < i+band; j++)
            if (A[i][j] == 0) printf("         ");   
            else              printf(" %+.1e",A[i][j]);
        printf(" :  %+.1e \n",B[i]); }
}
  
void femBandSystemPrintInfos(femBandSystem *myBand)
{
    int size = myBand->size;
    int band = myBand->band;
    printf(" \n");
    printf("    Banded Gaussian elimination \n");
    printf("    Storage informations \n");
    printf("    Matrix size      : %8d\n",size);
    printf("    Matrix band      : %8d\n",band);
    printf("    Bytes required   : %8d\n",(int)sizeof(double)*size*(band+1));     
}


double femBandSystemGet(femBandSystem* myBandSystem, int myRow, int myCol)
{
    double value = 0;
    if (myCol >= myRow && myCol < myRow+myBandSystem->band)  value = myBandSystem->A[myRow][myCol]; 
    return(value);
}

void femSystemConstrainNEUMANN(femProblem *theProblem,femFullSystem *mySystem, 
                             int iEdge, double value,femBoundaryType type,femElasticCase iCase){
    
    double jac,nx,ny,x[2],y[2],dx,dy,len,r;
    double  **A, *B;
    int     i, size,Nmap[2];
    A    = mySystem->A ; 
    B    = mySystem->B ; 
    int nodeL= theProblem->geometry->theEdges->elem[iEdge*2];
    int nodeR= theProblem->geometry->theEdges->elem[iEdge*2+1];
    double *X=theProblem->geometry->theNodes->X;
    double *Y=theProblem->geometry->theNodes->Y;
    x[0]=X[nodeL];
    x[1]=X[nodeR];
    y[0]=Y[nodeL];
    y[1]=Y[nodeR];
    // printf("x[0], x[1], y[0], y[1] %f %f %f %f\n",x[0],x[1],y[0],y[1]);  

    dx=x[1]-x[0];
    dy=y[1]-y[0];
    len=sqrt(dx*dx+dy*dy);
    nx=dy/(len);
    ny=-dx/(len);
    jac=len/2;
    int*renumber = theProblem->geometry->theNodes->number;
    Nmap[0]=renumber[nodeL];
    Nmap[1]=renumber[nodeR];
    printf("Nmap[0], Nmap[1] %d %d\n",Nmap[0],Nmap[1]);
    printf("value*jac*ny : %f\n",value*jac*ny);
    if(iCase == AXISYM){
        for (int iInteg=0;iInteg<2;iInteg++)
        {
            double phi[2] = {(1. - _gaussDos2Eta[iInteg]) / 2., (1. + _gaussDos2Eta[iInteg]) / 2.};
            for (size_t i = 0; i < 2; i++)
            {
                r += phi[i] * x[i];
            }
            
            for(int i=0;i<2;i++)
            {
                if (type==NEUMANN_X) 
                    B[2*Nmap[i]]+= phi[i]*value*jac*_gaussDos2Weight[iInteg] *r;
                if (type==NEUMANN_Y) 
                    B[2*Nmap[i]+1]+= phi[i]*value*jac*_gaussDos2Weight[iInteg]*r;
                if (type==NEUMANN_N)
                {
                    B[2*Nmap[i]]+= phi[i]*value*jac*nx*_gaussDos2Weight[iInteg]*r;
                    B[2*Nmap[i]+1]+= phi[i]*value*jac*ny*_gaussDos2Weight[iInteg]*r;
                }
                if (type==NEUMANN_T)
                {
                    B[2*Nmap[i]]+= phi[i]*value*jac*ny*_gaussDos2Weight[iInteg]*r;
                    B[2*Nmap[i]+1]-= phi[i]*value*jac*nx*_gaussDos2Weight[iInteg]*r;
                }
            }
        }
    }
    else{
        for (int i=0;i<2;i++)
        {
            if (type==NEUMANN_X) 
                B[2*Nmap[i]]+= value*jac;
            if (type==NEUMANN_Y) 
                B[2*Nmap[i]+1]+= value*jac;
            if (type==NEUMANN_N)
            {
                B[2*Nmap[i]]+= value*jac*nx;
                B[2*Nmap[i]+1]+= value*jac*ny;
            }
            if (type==NEUMANN_T)
            {

                B[2*Nmap[i]]+= value*jac*ny ;
                B[2*Nmap[i]+1]-= value*jac*nx;

                

            }
        }
    }
}
femIterativeSolver *femIterativeSolverCreate(int size)
{
    femIterativeSolver *mySolver = malloc(sizeof(femIterativeSolver));
    mySolver->R = malloc(sizeof(double)*size*4);      
    mySolver->D = mySolver->R + size;       
    mySolver->S = mySolver->R + size*2;       
    mySolver->X = mySolver->R + size*3;       
    mySolver->size = size;
    femIterativeSolverInit(mySolver);
    return(mySolver);
}

void femIterativeSolverFree(femIterativeSolver *mySolver)
{
    free(mySolver->R);
    free(mySolver);
}

void femIterativeSolverInit(femIterativeSolver *mySolver)
{
    int i;
    mySolver->iter = 0;
    mySolver->error = 10.0e+12;
    for (i=0 ; i < mySolver->size*4 ; i++) 
        mySolver->R[i] = 0;        
}
 
void femIterativeSolverPrint(femIterativeSolver *mySolver)
{
    double  *R;
    int     i, size;
    R    = mySolver->R;
    size = mySolver->size;

    for (i=0; i < size; i++) {
        printf("%d :  %+.1e \n",i,R[i]); }
}

void femIterativeSolverPrintInfos(femIterativeSolver *mySolver)
{
    if (mySolver->iter == 1)     printf("\n    Iterative solver \n");
    printf("    Iteration %4d : %14.7e\n",mySolver->iter,mySolver->error);
}

int femIterativeSolverConverged(femIterativeSolver *mySolver)
{
    int  testConvergence = 0;
    if (mySolver->iter  > 3000)     testConvergence = -1;
    if (mySolver->error < 10.0e-6)  testConvergence = 1;
    return(testConvergence);
}

void femIterativeSolverAssemble(femIterativeSolver* mySolver, double *Aloc, double *Bloc, double *Uloc, int *map, int nLoc)
{
    int i,j;
    for (i = 0; i < nLoc; i++) { 
        int myRow = map[i];
        mySolver->R[myRow] -= Bloc[i];
        for(j = 0; j < nLoc; j++) {
            mySolver->R[myRow] += Aloc[i*nLoc+j]*Uloc[j]; }}
}



double *femIterativeSolverEliminate(femIterativeSolver *mySolver)
{
    mySolver->iter++;
    double error = 0.0; int i;
    for (i=0; i < mySolver->size; i++) {
        error += (mySolver->R[i])*(mySolver->R[i]);
        mySolver->X[i] = -mySolver->R[i]/5.0; 
        mySolver->R[i] = 0.0; }
        
    mySolver->error = sqrt(error);
    return(mySolver->X);
}
double *theGlobalCoord;

int compare(const void *nodeOne, const void *nodeTwo) 
{
    int *iOne = (int *)nodeOne;
    int *iTwo = (int *)nodeTwo;
    double diff = theGlobalCoord[*iOne] - theGlobalCoord[*iTwo];
    if (diff < 0)    return  1;
    if (diff > 0)    return -1;
    return  0;  
}

void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    int i, *inverse;
    
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < theMesh->nodes->nNodes; i++) 
                theMesh->nodes->number[i] = i;
            break;
        case FEM_XNUM : 
            inverse = malloc(sizeof(int)*theMesh->nodes->nNodes);
            for (i = 0; i < theMesh->nodes->nNodes; i++) 
                inverse[i] = i; 
            theGlobalCoord = theMesh->nodes->X;
            qsort(inverse, theMesh->nodes->nNodes, sizeof(int), compare);
            for (i = 0; i < theMesh->nodes->nNodes; i++)
                theMesh->nodes->number[inverse[i]] = i;
            free(inverse);  
            break;
        case FEM_YNUM : 
            inverse = malloc(sizeof(int)*theMesh->nodes->nNodes);
            for (i = 0; i < theMesh->nodes->nNodes; i++) 
                inverse[i] = i; 
            theGlobalCoord = theMesh->nodes->Y;
            qsort(inverse, theMesh->nodes->nNodes, sizeof(int), compare);
            for (i = 0; i < theMesh->nodes->nNodes; i++)
                theMesh->nodes->number[inverse[i]] = i;
            free(inverse);  
            break;
        default : Error("Unexpected renumbering option"); }
}
int femMeshComputeBand(femMesh *theMesh)
{
    int iElem,j,myMax,myMin,myBand,map[4];
    int nLocal = theMesh->nLocalNode;
    myBand = 0;
    for(iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; ++j) 
            map[j] = theMesh->nodes->number[theMesh->elem[iElem*nLocal+j]];
        myMin = map[0];
        myMax = map[0];
        for (j=1; j < nLocal; j++) {
            myMax = fmax(map[j],myMax);
            myMin = fmin(map[j],myMin); }
        if (myBand < (myMax - myMin)) myBand = myMax - myMin; }         
    return(++myBand);
}
void femBandSystemAssemble(femBandSystem* myBandSystem, double *Aloc, double *Bloc, int *map, int nLoc,int*renumber)
{
    int i,j,mapX_row[4],mapY_row[4],mapX_col[4],mapY_col[4];
    int colX,colY,rowX,rowY;

    for (i = 0; i < nLoc; i++) { 
           mapX_row[i] = 2*renumber[map[i]];
           mapY_row[i] = 2*renumber[map[i]] + 1;
           rowX=2*i;
           rowY=2*i+1;
        for(j = 0; j < nLoc; j++) {
            //printf("control assemblage i: %d map:%d j:%d map:%d\n",i,map[i],j,map[j]);   

           mapX_col[j] = 2*renumber[map[j]];
           mapY_col[j] = 2*renumber[map[j]] + 1; 
           colX=2*j;
           colY=2*j+1;  

         
        /*   if ( (i==1) & (j==2))
           {
            printf("i: %d j:%d, nloc:%d\n",i,j,nLoc);
            printf("map[i]:%d\n",map[i]);
            printf("map[j]:%d\n",map[j]);
            printf("renumber map[i]:%d\n",renumber[map[i]]);
            printf("renumber map[j]:%d\n",renumber[map[j]]);   
            printf("mapX_row[i]:%d\n",mapX_row[i]);
            printf("mapY_row[i]:%d\n",mapY_row[i]);
            printf("mapY_col[j]:%d\n",mapY_col[j]);
            printf("mapX_col[j]:%d\n",mapX_col[j]);
  
            printf("rowX[i]:%d\n",rowX);
            printf("rowY[i]:%d\n",rowY);
            printf("colX[j]:%d\n",colX);
            printf("colY[j]:%d\n",colY);

            //printf("(2*i)*(2*nLoc+(2*j): %d\n",(2*i+1)*(2*nLoc)+(2*j+1));
            //printf("Aloc: %f\n",Aloc[(2*i+1)*(2*nLoc)+(2*j+1)]);
            //printf("Aloc[3][3] pointer: %p\n",&(Aloc[(2*1+1)*(2*nLoc)+(2*1+1)]));
           }*/
           


            if (mapX_col[j] >= mapX_row[i])  
            myBandSystem->A[mapX_row[i]][mapX_col[j]] += Aloc[rowX*2*nLoc+colX]; 
            
            if (mapY_col[j] >= mapY_row[i])  
            myBandSystem->A[mapY_row[i]][mapY_col[j]] += Aloc[rowY*2*nLoc+colY]; 
            
            if (mapX_col[j] >= mapY_row[i])  
            myBandSystem->A[mapY_row[i]][mapX_col[j]] += Aloc[rowY*2*nLoc+colX]; 
            
            if (mapY_col[j] >= mapX_row[i])  
            myBandSystem->A[mapX_row[i]][mapY_col[j]] += Aloc[rowX*2*nLoc+colY]; 
         /*if ( (i==1) & (j==1))
           {
            printf("A %d %d : %f\n",mapX_row[i],mapX_col[j],myBandSystem->A[mapX_row[i]][mapX_col[j]]);
            printf("A %d %d : %f\n",mapY_row[i],mapY_col[j],myBandSystem->A[mapY_row[i]][mapY_col[j]]);
            //printf("A %d %d %d  \n",mapY_row[i],myBandSystem->band,mapY_col[j]);
    
            //printf("A %d  : %f\n",mapY_row[i]*myBandSystem->band+mapY_col[j],myBandSystem->B[mapY_row[i]*myBandSystem->band+mapY_col[j]+myBandSystem->size]);
    
            printf("A %d %d : %f\n",mapY_row[i],mapX_col[j],myBandSystem->A[mapY_row[i]][mapX_col[j]]);
            printf("A %d %d : %f\n",mapX_row[i],mapY_col[j],myBandSystem->A[mapX_row[i]][mapY_col[j]]);
           }   */     
        } /*for j*/
       if ( mapX_row[i] < myBandSystem->size)
        {
        myBandSystem->B[mapX_row[i]] += Bloc[i*2]; 
        }
        else 
        {
            printf("Error: assemblage B X i: %d map:%d\n",i,map[i]);
        }
       if ( mapY_row[i] < myBandSystem->size)
        {
        myBandSystem->B[mapY_row[i]] += Bloc[i*2+1]; 
        }
        else 
        {
            printf("Error: assemblage B Y i: %d map:%d\n",i,map[i]);
        }


        }/*for i*/
}


double  *femBandSystemEliminate(femBandSystem *myBand)
{
    double  **A, *B, factor;
    int     i, j, k, jend, size, band;
    A    = myBand->A;
    B    = myBand->B;
    size = myBand->size;
    band = myBand->band;

    /* Incomplete Cholesky factorization */ 

    for (k=0; k < size; k++) {
        if ( fabs(A[k][k]) <= 1e-4 ) {
            Error("Cannot eleminate with such a pivot"); }
        jend = fmin(k + band,size);
        for (i = k+1 ; i <  jend; i++) {
            factor = A[k][i] / A[k][k];
            for (j = i ; j < jend; j++) 
                A[i][j] = A[i][j] - A[k][j] * factor;
            B[i] = B[i] - B[k] * factor; }}
        
    /* Back-substitution */

    for (i = (size-1); i >= 0 ; i--) {
        factor = 0;
        jend = fmin(i + band,size);
        for (j = i+1 ; j < jend; j++)
            factor += A[i][j] * B[j];
        B[i] = ( B[i] - factor)/A[i][i]; }

    return(myBand->B);
}


void  femSystemConstrainDIRICHLETXY(femFullSystem *mySystem, 
                             int myNode, double myValue,femBoundaryType myType) 
{ 
    double  **A, *B;
    int     i, size;
    A    = mySystem->A ; 
    B    = mySystem->B ; 
    size = mySystem->size;
    if(myType == DIRICHLET_X || myType == DIRICHLET_Y){
     for (i=0; i < size; i++) {
        B[i] -= myValue * A[i][myNode];
        A[i][myNode] = 0; }
    
        for (i=0; i < size; i++) 
            A[myNode][i] = 0; 
        
        A[myNode][myNode] = 1;
        B[myNode] = myValue;
    }
}


femProblem *femElasticityCreate(femGeo* theGeometry, 
                  double E, double nu, double rho, double g, femElasticCase iCase, femSolverType type,femRenumType renumType)
{
    femProblem *theProblem = malloc(sizeof(femProblem));
    theProblem->E   = E;
    theProblem->nu  = nu;
    theProblem->g   = g;
    theProblem->rho = rho;
    
    if (iCase == PLANAR_STRESS) {
        theProblem->A = E/(1-nu*nu);
        theProblem->B = E*nu/(1-nu*nu);
        theProblem->C = E/(2*(1+nu)); }
    else if (iCase == PLANAR_STRAIN || iCase == AXISYM) {
        theProblem->A = E*(1-nu)/((1+nu)*(1-2*nu));
        theProblem->B = E*nu/((1+nu)*(1-2*nu));
        theProblem->C = E/(2*(1+nu)); }

    theProblem->planarStrainStress = iCase;
    theProblem->nBoundaryConditions = 0;
    theProblem->conditions = NULL;
    
    int size = 2*theGeometry->theNodes->nNodes;
    int number_Edges= theGeometry->theEdges->nElem;
    theProblem->contrainteEdges = malloc(number_Edges*sizeof(int)) ; 
    theProblem->constrainedNodes = malloc(size*sizeof(int));
    for (int i=0; i < size; i++) 
        theProblem->constrainedNodes[i] = -1;
    for(int i=0; i < number_Edges; i++)
        theProblem->contrainteEdges[i] = -1;
    theProblem->geometry = theGeometry;  
    if (theGeometry->theElements->nLocalNode == 3) {
        theProblem->space    = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule     = femIntegrationCreate(3,FEM_TRIANGLE); }
    if (theGeometry->theElements->nLocalNode == 4) {
        theProblem->space    = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule     = femIntegrationCreate(4,FEM_QUAD); }
    femMeshRenumber(theProblem->geometry->theElements,renumType);   
    int band ; 
    switch (type) {
        case FEM_FULL : 
                theProblem->solver = femSolverFullCreate(size); 
                break;
        case FEM_BAND : 
                band = femMeshComputeBand(theGeometry->theElements); //theProblem->mesh
                theProblem->solver = femSolverBandCreate(size,
                                                         2*theGeometry->theElements->nLocalNode,band); break;
        case FEM_ITER : 
               theProblem->solver = femSolverIterativeCreate(size,
                                                             2*theGeometry->theElements->nLocalNode); break;
        case FEM_Cholesky : 
                theProblem->solver = femSolverFullCreate(size); break;
        case FEM_Cholesky_band : 
                band = femMeshComputeBand(theGeometry->theElements); //theProblem->mesh
                theProblem->solver = femSolverBandCreate(size,
                                                         2*theGeometry->theElements->nLocalNode,band); break;
        default : Error("Unexpected solver option"); }

    return theProblem;
}

void femElasticityFree(femProblem *theProblem)
{
    femSolverFree(theProblem->solver);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    free(theProblem->conditions);
    free(theProblem->constrainedNodes);
    free(theProblem);
}
    
void femElasticityAddBoundaryCondition(femProblem *theProblem, char *nameDomain, femBoundaryType type, double value)
{
    int iDomain = geoGetDomain(nameDomain);
    if (iDomain == -1)  Error("Undefined domain :-(");

    femBoundaryCondition* theBoundary = malloc(sizeof(femBoundaryCondition));
    theBoundary->domain = theProblem->geometry->theDomains[iDomain];
    theBoundary->value = value;
    theBoundary->type = type;
    theProblem->nBoundaryConditions++;
    int size = theProblem->nBoundaryConditions;
    
    if (theProblem->conditions == NULL)
        theProblem->conditions = malloc(size*sizeof(femBoundaryCondition*));
    else 
        theProblem->conditions = realloc(theProblem->conditions, size*sizeof(femBoundaryCondition*));
    theProblem->conditions[size-1] = theBoundary;
    
    
    int shift;
    if (type == DIRICHLET_X )  shift = 0;      
    if (type == DIRICHLET_Y ) shift = 1;    
    int *elem = theBoundary->domain->elem;
    int nElem = theBoundary->domain->nElem;

    printf("nElem %d \n",nElem);
    for (int e=0; e<nElem; e++) {
        for (int i=0; i<2; i++) {
            int node = theBoundary->domain->mesh->elem[2*elem[e]+i];
            theProblem->constrainedNodes[2*node+shift] = size-1; 
        }
        if ((type == NEUMANN_X)|| (type == NEUMANN_Y)||(type == NEUMANN_N)|| (type == NEUMANN_T))
            theProblem ->contrainteEdges[elem[e]] =  size-1 ; 
    }

}

void femElasticityPrint(femProblem *theProblem)  
{    
    printf("\n\n ======================================================================================= \n\n");
    printf(" Linear elasticity problem \n");
    printf("   Young modulus   E   = %14.7e [N/m2]\n",theProblem->E);
    printf("   Poisson's ratio nu  = %14.7e [-]\n",theProblem->nu);
    printf("   Density         rho = %14.7e [kg/m3]\n",theProblem->rho);
    printf("   Gravity         g   = %14.7e [m/s2]\n",theProblem->g);
    
    if (theProblem->planarStrainStress == PLANAR_STRAIN)  printf("   Planar strains formulation \n");
    if (theProblem->planarStrainStress == PLANAR_STRESS)  printf("   Planar stresses formulation \n");
    if (theProblem->planarStrainStress == AXISYM)         printf("   Axisymmetric formulation \n");

    
    printf("   Boundary conditions : \n");
    for(int i=0; i < theProblem->nBoundaryConditions; i++) {
          femBoundaryCondition *theCondition = theProblem->conditions[i];
          double value = theCondition->value;
          printf("  %20s :",theCondition->domain->name);
          if (theCondition->type==DIRICHLET_X)  printf(" imposing %9.2e as the horizontal displacement  \n",value);
          if (theCondition->type==DIRICHLET_Y)  printf(" imposing %9.2e as the vertical displacement  \n",value); }
    printf(" ======================================================================================= \n\n");
}




double femMin(double *x, int n) 
{
    double myMin = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMin = fmin(myMin,x[i]);
    return myMin;
}

double femMax(double *x, int n) 
{
    double myMax = x[0];
    int i;
    for (i=1 ;i < n; i++) 
        myMax = fmax(myMax,x[i]);
    return myMax;
}


void femError(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    exit(69);                                                 
}

void femErrorGmsh(int ierr, int line, char *file)                                  
{ 
    if (ierr == 0)  return;
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in %s at line %d : \n  error code returned by gmsh %d\n", file, line, ierr);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");
    gmshFinalize(NULL);                                        
    exit(69);                                                 
}

void femErrorScan(int test, int line, char *file)                                  
{ 
    if (test >= 0)  return;
    
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Error in fscanf or fgets in %s at line %d : \n", file, line);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");   
    exit(69);                                       
}

void femWarning(char *text, int line, char *file)                                  
{ 
    printf("\n-------------------------------------------------------------------------------- ");
    printf("\n  Warning in %s at line %d : \n  %s\n", file, line, text);
    printf("--------------------------------------------------------------------- Yek Yek !! \n\n");                                              
}
femProblem* femElasticityRead(femGeo* theGeometry, const char *filename, femSolverType type,femRenumType renumType)
{
    FILE* file = fopen(filename,"r");
    femProblem *theProblem = malloc(sizeof(femProblem));
    theProblem->nBoundaryConditions = 0;
    theProblem->conditions = NULL;
    
    int size = 2*theGeometry->theNodes->nNodes;
    theProblem->constrainedNodes = malloc(size*sizeof(int));
    int nbr_Edges= theGeometry->theEdges->nElem;
    theProblem->contrainteEdges = malloc(nbr_Edges*sizeof(int)) ;

    for (int i=0; i < size; i++) 
        theProblem->constrainedNodes[i] = -1;
    for (int i=0; i < nbr_Edges; i++) 
    
         theProblem->contrainteEdges[i] = -1;
    
    theProblem->geometry = theGeometry;  
    if (theGeometry->theElements->nLocalNode == 3) {
        theProblem->space    = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule     = femIntegrationCreate(3,FEM_TRIANGLE); }
    if (theGeometry->theElements->nLocalNode == 4) {
        theProblem->space    = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule     = femIntegrationCreate(4,FEM_QUAD); }
    int band ;
    switch (type) {
        case FEM_FULL : 
                theProblem->solver = femSolverFullCreate(size); break;
        case FEM_BAND : 
                band = femMeshComputeBand(theGeometry->theElements); //theProblem->mesh
                theProblem->solver = femSolverBandCreate(size,
                                                         2*theGeometry->theElements->nLocalNode,band); break;
        case FEM_ITER : 
               theProblem->solver = femSolverIterativeCreate(size,
                                                             2*theGeometry->theElements->nLocalNode); break;
        default : Error("Unexpected solver option"); }

    return theProblem;


    char theLine[MAXNAME];
    char theDomain[MAXNAME];
    char theArgument[MAXNAME];
    double value;
    double typeCondition;

    while (!feof(file)) {
        ErrorScan(fscanf(file,"%19[^\r\n]s \r\n",(char *)&theLine));
        if (strncasecmp(theLine,"Type of problem     ",19) == 0) {
            ErrorScan(fscanf(file,":  %[^\n]s \n",(char *)&theArgument));
            if (strncasecmp(theArgument,"Planar stresses",13) == 0)
               theProblem->planarStrainStress = PLANAR_STRESS; 
            if (strncasecmp(theArgument,"Planar strains",13) == 0)
               theProblem->planarStrainStress = PLANAR_STRAIN; 
            if (strncasecmp(theArgument,"Axi-symetric problem",13) == 0)
               theProblem->planarStrainStress = AXISYM; }
        if (strncasecmp(theLine,"Young modulus       ",19) == 0) {
            ErrorScan(fscanf(file,":  %le\n",&theProblem->E)); }
        if (strncasecmp(theLine,"Poisson ratio       ",19) == 0) {
            ErrorScan(fscanf(file,":  %le\n",&theProblem->nu)); }
        if (strncasecmp(theLine,"Mass density        ",19) == 0) {
            ErrorScan(fscanf(file,":  %le\n",&theProblem->rho)); }
        if (strncasecmp(theLine,"Gravity             ",19) == 0) {
            ErrorScan(fscanf(file,":  %le\n",&theProblem->g)); }
        if (strncasecmp(theLine,"Boundary condition  ",19) == 0) {
            ErrorScan(fscanf(file,":  %19s = %le : %[^\n]s\n",(char *)&theArgument,&value,(char *)&theDomain));
            if (strncasecmp(theArgument,"Dirichlet-X",19) == 0)
                typeCondition = DIRICHLET_X;
            if (strncasecmp(theArgument,"Dirichlet-Y",19) == 0)
                typeCondition = DIRICHLET_Y;                
            if (strncasecmp(theArgument,"Neumann-X",19) == 0)
                typeCondition = NEUMANN_X;
            if (strncasecmp(theArgument,"Neumann-Y",19) == 0)
                typeCondition = NEUMANN_Y;                
            femElasticityAddBoundaryCondition(theProblem,theDomain,typeCondition,value); }
        ErrorScan(fscanf(file,"\r\n")); }
 
    int iCase = theProblem->planarStrainStress;
    double E = theProblem->E;
    double nu = theProblem->nu;
    
    if (iCase == PLANAR_STRESS) {
        theProblem->A = E/(1-nu*nu);
        theProblem->B = E*nu/(1-nu*nu);
        theProblem->C = E/(2*(1+nu)); }
    else if (iCase == PLANAR_STRAIN || iCase == AXISYM) {
        theProblem->A = E*(1-nu)/((1+nu)*(1-2*nu));
        theProblem->B = E*nu/((1+nu)*(1-2*nu));
        theProblem->C = E/(2*(1+nu)); }
    femMeshRenumber(theProblem->geometry->theElements,renumType);

    fclose(file);
    return theProblem;
}
double **choleskyDecomposition(double **A, int n) {
    int i, j, k;
    double **L = malloc(n * sizeof(double *));
    for (i = 0; i < n; i++) {
        L[i] = malloc(n * sizeof(double));
    }
    L[0][0] = sqrt(A[0][0]);
    for (i = 0; i < n; i++) {
        for (j = 0; j < (i + 1); j++) {
            double sum = 0.0;

            if (i == j) {
                for (k = 0; k < j; k++) {
                    sum += pow(L[j][k], 2);
                }
                L[j][j] = sqrt(A[j][j] - sum);
            } else {
                for (k = 0; k < j; k++) {
                    sum += L[i][k] * L[j][k];
                }
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }
        }
    }
    return L;
}
double *solvecholesky(double **L, double *b, int n) {
    int i, j;
    double y[n];
    double *x = malloc(n * sizeof(double));
    y[0] = b[0] / L[0][0];
    for (i = 1; i < n; i++) {
        double sum = 0.0;
        for (j = 0; j < i; j++) {
            sum += L[i][j] * y[j];
        }
        y[i] = (b[i] - sum) / L[i][i];
    }
    x[n - 1] = y[n - 1] / L[n - 1][n - 1];
    for (i = n - 2; i >= 0; i--) {
        double sum = 0.0;
        for (j = i + 1; j < n; j++) {
            sum += L[j][i] * x[j];
        }
        x[i] = (y[i] - sum) / L[i][i];
    }
    return x;
}
double *            reorder(double *B,femNodes*Nodes,int* renumber)
{
  double *C;
  C=malloc(sizeof(double)*Nodes->nNodes*2);
 
  for (int i = 0; i < Nodes->nNodes; ++i)
  {
  C[2*i]=B[2*renumber[i]];
  C[2*i+1]=B[2*renumber[i]+1];
  }  
   for (int i = 0; i < 2*(Nodes->nNodes); ++i)
  {
  B[i]=C[i];
  }  
  free(C);
 return B;
}