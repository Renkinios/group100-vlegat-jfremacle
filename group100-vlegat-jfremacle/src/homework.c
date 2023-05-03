#include "fem.h"

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
    femNodes *nodes = theMesh->nodes;
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < nodes->nNodes; i++) 
                nodes->number[i] = i;
            break;
        case FEM_XNUM : 
            inverse = malloc(sizeof(int)*nodes->nNodes);
            for (i = 0; i < nodes->nNodes; i++) 
                inverse[i] = i; 
            theGlobalCoord = nodes->X;
            qsort(inverse, nodes->nNodes, sizeof(int), compare);
            for (i = 0; i < nodes->nNodes; i++)
                nodes->number[inverse[i]] = i;
            free(inverse);  
            break;
        case FEM_YNUM : 
            inverse = malloc(sizeof(int)*nodes->nNodes);
            for (i = 0; i < nodes->nNodes; i++) 
                inverse[i] = i; 
            theGlobalCoord = nodes->Y;
            qsort(inverse, nodes->nNodes, sizeof(int), compare);
            for (i = 0; i < nodes->nNodes; i++)
                nodes->number[inverse[i]] = i;
            free(inverse);  
            break;
        default : Error("Unexpected renumbering option"); }
}


double *femElasticitySolve(femProblem *theProblem)
{

    femFullSystem  *theSystem = theProblem->system;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theMesh = theGeometry->theElements;
    
    
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4];
    int iElem,iInteg,iEdge,i,j,d,map[4],mapX[4],mapY[4];
    
    int nLocal = theMesh->nLocalNode;

    double a   = theProblem->A;
    double b   = theProblem->B;
    double c   = theProblem->C;      
    double rho = theProblem->rho;
    double g   = theProblem->g;
    double **A = theSystem->A;
    double *B  = theSystem->B;
    
    
  //
  //  A faire :-)
  //                
                
                
  
    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value); }}
                            
    return femFullSystemEliminate(theSystem);
}
