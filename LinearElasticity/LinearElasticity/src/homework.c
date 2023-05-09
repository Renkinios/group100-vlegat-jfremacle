#include "fem.h"




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
    
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];} 
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            // x like r y like z
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; } 
                for (i = 0; i < theSpace->n; i++) { 
                    for(j = 0; j < theSpace->n; j++) {
                        if(theProblem->planarStrainStress == PLANAR_STRESS || theProblem->planarStrainStress == PLANAR_STRAIN){
                            A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + 
                                                    dphidy[i] * c * dphidy[j]) * jac * weight;                                                                                            
                            A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + 
                                                    dphidy[i] * c * dphidx[j]) * jac * weight;                                                                                           
                            A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + 
                                                    dphidx[i] * c * dphidy[j]) * jac * weight;                                                                                            
                            A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + 
                                                    dphidx[i] * c * dphidx[j]) * jac * weight;
                        }
                        if(theProblem->planarStrainStress == AXISYM ) {
                             A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] *x[i] + 
                                                    dphidy[i] * c * dphidy[j] * x[i] + phi[i]*(b*dphidx[j]) + a*phi[j]/x[i]) * jac * weight;                                                                                            
                            A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] *x[i]+ x[i] *
                                                    dphidy[i] * c * dphidx[j] + phi[i] *b* dphidy[j]) * jac * weight;                                                                                           
                            A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] * x[i]+ 
                                                    dphidx[i] * c * dphidy[j] * x[i] + dphidy[i] * b * phi[j]) * jac * weight;                                                                                            
                            A[mapY[i]][mapY[j]] += (dphidy[i] * a * x[i] * dphidy[j] + 
                                                    dphidx[i] * c * x[i] *dphidx[j]) * jac * weight;
                        }
                    }
                }
            
    
             for (i = 0; i < theSpace->n; i++) {
                if(theProblem->planarStrainStress == PLANAR_STRESS || theProblem->planarStrainStress == PLANAR_STRAIN)
                    B[mapY[i]] -= phi[i] * g * rho * jac * weight; 
                if(theProblem->planarStrainStress == AXISYM)
                    B[mapY[i]] -= phi[i] * x[i] * g * rho * jac * weight * x[i];

  
    int *theConstrainedNodes = theProblem->constrainedNodes; 

    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femFullSystemConstrain(theSystem,i,value,theProblem->conditions[theConstrainedNodes[i]]->type); 
           }
        }
            
            
    double *sol ;
    double **L ;
    switch (theProblem->system->type)
    {
    case FEM_Cholesky:
        L  = choleskyDecomposition(A,theSystem->size);  
        sol = solvecholesky(L,theSystem->B,theSystem->size) ;
        break;
    case FEM_FULL : 
        sol = femFullSystemEliminate(theSystem);
        break;
    default:
        Error("Unexpected solver option");
        break;
    }   
    return sol; 
}
