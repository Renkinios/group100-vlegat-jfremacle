#include "fem.h"

double *X_or_Y ; 
int compare_node(const void *a,const void *b ){
    int node_a = *((int*)a) ; 
    int node_b = *((int*)b) ;
    return X_or_Y[node_b] - X_or_Y[node_a] ;  

}


void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    
    int ierr;
    double r = w/4;
    int idRect = gmshModelOccAddRectangle(0.0,0.0,0.0,w,h,-1,0.0,&ierr); 
    int idDisk = gmshModelOccAddDisk(w/2.0,h/2.0,0.0,r,r,-1,NULL,0,NULL,0,&ierr); 
    int idSlit = gmshModelOccAddRectangle(w/2.0,h/2.0-r,0.0,w,2.0*r,-1,0.0,&ierr); 
    int rect[] = {2,idRect};
    int disk[] = {2,idDisk};
    int slit[] = {2,idSlit};

    gmshModelOccCut(rect,2,disk,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    gmshModelOccCut(rect,2,slit,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr); 
    gmshModelOccSynchronize(&ierr); 

    if (theGeometry->elementType == FEM_QUAD) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshOptionSetNumber("Mesh.RecombineAll",1,&ierr);
        gmshOptionSetNumber("Mesh.Algorithm",8,&ierr);  
        gmshOptionSetNumber("Mesh.RecombinationAlgorithm",1.0,&ierr); 
        gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  
        gmshModelMeshGenerate(2,&ierr);  }
  
    if (theGeometry->elementType == FEM_TRIANGLE) {
        gmshOptionSetNumber("Mesh.SaveAll",1,&ierr);
        gmshModelMeshGenerate(2,&ierr);  }
 
    return;
}
void femMeshRenumber(femMesh *theMesh, femRenumType renumType)
{
    int i, *inverse;
    femNodes *nodes = theMesh->nodes;
    switch (renumType) {
        case FEM_NO :
            for (i = 0; i < nodes->nNode; i++) 
                nodes->number[i] = i;
            break;
        case FEM_XNUM : 
            inverse = malloc(sizeof(int)*nodes->nNode);
            for (i = 0; i < nodes->nNode; i++) 
                inverse[i] = i; 
            theGlobalCoord = nodes->X;
            qsort(inverse, nodes->nNode, sizeof(int), compare);
            for (i = 0; i < nodes->nNode; i++)
                nodes->number[inverse[i]] = i;
            free(inverse);  
            break;
        case FEM_YNUM : 
            inverse = malloc(sizeof(int)*nodes->nNode);
            for (i = 0; i < nodes->nNode; i++) 
                inverse[i] = i; 
            theGlobalCoord = theMesh->Y;
            qsort(inverse, nodes->nNode, sizeof(int), compare);
            for (i = 0; i < nodes->nNode; i++)
                theMesh->number[inverse[i]] = i;
            free(inverse);  
            break;
        default : Error("Unexpected renumbering option"); }
}
    int node[nodes->nNodes] ; 
    for (size_t i = 0; i < nodes->nNodes; i++)
        node[i] = i ; 
    switch (renumType) {
        case FEM_NO :
            break;
        case FEM_XNUM : 
            X_or_Y = nodes ->X ;  
            qsort(node,nodes->nNodes,sizeof(int),compare_node) ; 
            break;
        case FEM_YNUM : 
            X_or_Y = nodes->Y ;  
            qsort(node,nodes->nNodes,sizeof(int),compare_node) ;  // regarde ou sont les plus grand_noeud pour après pouvoir les remplacé 
            break;            

        default : Error("Unexpected renumbering option"); 
        }
    for (size_t i = 0; i < nodes->nNodes; i++)
        nodes -> number[node[i]] = i ; // remplace par les premier noeud  
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
