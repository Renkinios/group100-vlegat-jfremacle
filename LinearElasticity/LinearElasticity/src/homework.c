#include "fem.h"




void geoMeshGenerateGeo() {

    femGeo* theGeometry = geoGetGeometry();
    // geoSetSizeCallback(geoSize);   

       int ierr;
    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
    double r = w/4;
    double lc = theGeometry->h;

    double e = 0.07;
    double HAUT = 1 ; 
    double L = 1.8 ;
    double jump = 0.2 ;
    int p1 = gmshModelGeoAddPoint(0, 0, 0., lc, 1, &ierr);
    int p2 = gmshModelGeoAddPoint(0, HAUT+e, 0., lc, 2, &ierr);
    int p3 = gmshModelGeoAddPoint(HAUT/4,    HAUT+e, 0., lc, 3, &ierr);
    int p4 = gmshModelGeoAddPoint(L-jump,  HAUT+e, 0., lc, 4, &ierr);
    int p5 = gmshModelGeoAddPoint(L,    HAUT+e, 0., lc, 5, &ierr);
    int p6 = gmshModelGeoAddPoint(L,      HAUT, 0., lc, 6, &ierr);
    int p7 = gmshModelGeoAddPoint(HAUT/4,     HAUT, 0., lc, 7, &ierr);
    int p8 = gmshModelGeoAddPoint(HAUT/4,   0, 0., lc, 8, &ierr);


    int l1 = gmshModelGeoAddLine(p1, p2, 1, &ierr);
    int l2 = gmshModelGeoAddLine(p2, p3, 2, &ierr);
    int l3 = gmshModelGeoAddLine(p3, p4, 3, &ierr);
    int l4 = gmshModelGeoAddLine(p4, p5, 4, &ierr);
    int l5 = gmshModelGeoAddLine(p5, p6, 5, &ierr);
    int l6 = gmshModelGeoAddLine(p6, p7, 6, &ierr);
    int l7 = gmshModelGeoAddLine(p7, p8, 7, &ierr);
    int l8 = gmshModelGeoAddLine(p8, p1, 8, &ierr); 

    int lTags[] = {l1, l2, l3, l4, l5, l6,l7,l8};
    int c1[] = {1};
    c1[0] = gmshModelGeoAddCurveLoop(lTags, 8, 1, 0, &ierr);  
    int s1 = gmshModelGeoAddPlaneSurface(c1, 1, 1, &ierr);
    gmshModelGeoSynchronize(&ierr);


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


 //   gmshFltkRun(&ierr);
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
    femMesh        *theEdges = theGeometry->theEdges;
    
    
    double x[4],y[4],phi[4],dphidxsi[4],dphideta[4],dphidx[4],dphidy[4],r;
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
                dydeta += y[i]*dphideta[i]; 
                r += x[i]*phi[i];
            }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            
            for (i = 0; i < theSpace->n; i++) {    
                dphidx[i] = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                dphidy[i] = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; 
            }            
            for (i = 0; i < theSpace->n; i++) { 
                for(j = 0; j < theSpace->n; j++) {
                    if (theProblem->planarStrainStress == PLANAR_STRESS || theProblem->planarStrainStress == PLANAR_STRAIN)
                    {
                        A[mapX[i]][mapX[j]] += (dphidx[i] * a * dphidx[j] + 
                                                dphidy[i] * c * dphidy[j]) * jac * weight;                                                                                            
                        A[mapX[i]][mapY[j]] += (dphidx[i] * b * dphidy[j] + 
                                                dphidy[i] * c * dphidx[j]) * jac * weight;                                                                                           
                        A[mapY[i]][mapX[j]] += (dphidy[i] * b * dphidx[j] + 
                                                dphidx[i] * c * dphidy[j]) * jac * weight;                                                                                            
                        A[mapY[i]][mapY[j]] += (dphidy[i] * a * dphidy[j] + 
                                                dphidx[i] * c * dphidx[j]) * jac * weight; 
                    }
                    if (theProblem->planarStrainStress == AXISYM)
                    {
                        A[mapX[i]][mapX[j]] += (dphidx[i] * a * r * dphidx[j] + 
                                            dphidy[i] * c * r * dphidy[j] + dphidx[i] * b * phi[j] + phi[i]*(b*dphidx[j] + a*phi[j]/r)) * jac * weight;                                                                                            
                        A[mapX[i]][mapY[j]] += (dphidx[i] * b * r * dphidy[j] + 
                                            dphidy[i] * c * r * dphidx[j] + b*phi[i]*dphidy[j]) * jac * weight;                                                                                           
                        A[mapY[i]][mapX[j]] += (dphidy[i] * b * r * dphidx[j] + 
                                            dphidx[i] * c * r * dphidy[j] + b*dphidy[i]*phi[j]) * jac * weight;                                                                                            
                        A[mapY[i]][mapY[j]] += (dphidy[i] * a * r * dphidy[j] + 
                                            dphidx[i] * c * r * dphidx[j]) * jac * weight; 
                    }
                    
                }
            }
            for (i = 0; i < theSpace->n; i++) {
                if (theProblem->planarStrainStress == AXISYM)
                {
                    B[mapY[i]] -= phi[i] * g * rho * jac * weight * r ;
                }       
                if (theProblem->planarStrainStress == PLANAR_STRESS || theProblem->planarStrainStress == PLANAR_STRAIN)
                {
                    B[mapY[i]] -= phi[i] * g * rho * jac * weight ;
                }
            }
        }
    }  
    int *theConstrainedEdges = theProblem->contrainteEdges; 
    for  (int iEdge=0; iEdge < theEdges->nElem; iEdge++)
    {
        if (theConstrainedEdges[iEdge] != -1) {
            double value = theProblem->conditions[theConstrainedEdges[iEdge]]->value;
            int type = theProblem->conditions[theConstrainedEdges[iEdge]]->type;
            double jac,nx,ny;
            for (int iInteg=0;iInteg<2;iInteg++)
            {
                double phi[2] = {(1. - _gaussDos2Eta[iInteg]) / 2., (1. + _gaussDos2Eta[iInteg]) / 2.};
                int Nmap[2]; 
                getEdge(theProblem,iEdge,&jac,&nx,&ny,Nmap); 
                for(int i=0;i<2;i++)
                {
                    if (type==NEUMANN_X) B[2*Nmap[i]]+= phi[i]*value*jac*_gaussDos2Weight[iInteg];
                    if (type==NEUMANN_Y) B[2*Nmap[i]+1]+= phi[i]*value*jac*_gaussDos2Weight[iInteg];
                    if (type==NEUMANN_N)
                    {
                        if (nx!=0) B[2*Nmap[i]]+= phi[i]*value*jac*nx*_gaussDos2Weight[iInteg];
                        if (ny!=0) B[2*Nmap[i]+1]+= phi[i]*value*jac*ny*_gaussDos2Weight[iInteg];
                    }
                    if (type==NEUMANN_T)
                    {
                        if (ny!=0) B[2*Nmap[i]]+= phi[i]*value*jac*ny*_gaussDos2Weight[iInteg];
                        if (nx!=0) B[2*Nmap[i]+1]-= phi[i]*value*jac*nx*_gaussDos2Weight[iInteg];
                    }
                }
            }
        }
    }
    //     for(i = 0; i < theProblem->nBoundaryConditions; i++){
    //     femBoundaryCondition *theCondition = theProblem->conditions[i];
    //     if((theCondition->type == NEUMANN_N)||(theCondition->type == NEUMANN_T)||(theCondition->type == NEUMANN_X)||(theCondition->type == NEUMANN_Y)){
    //         double val = theCondition->value;
    //         femDomain *Dom = theCondition->domain;
    //         for(j = 0; j < theEdges->nElem; j++){

    //             for(int k = 0; k < 2; k++){
    //                 double phi[2] = {(1.0 - _gaussDos2Eta[k])/2,(1 - _gaussDos2Eta[k])/2};
    //                 int nodeL= theProblem->geometry->theEdges->elem[j*2];
    //                 printf("nodeL = %d\n", nodeL);
    //                 int nodeR= theProblem->geometry->theEdges->elem[j*2+1];
    //                 printf("nodeR = %d\n", nodeR);
    //                 double *X=theProblem->geometry->theNodes->X;
    //                 double *Y=theProblem->geometry->theNodes->Y;
    //                 double X1=X[nodeL];
    //                 double X2=X[nodeR];
    //                 double Y1=Y[nodeL];
    //                 double Y2=Y[nodeR];
    //                 printf(" X2 = %lf , Y2 = %lf",X2,Y2) ;  
    //                 printf(" X1 = %lf , Y1 = %lf",X1,Y1) ;
    //                 double length = sqrt((X1 - X2)*(X1 - X2) + (Y1-Y2)*(Y1-Y2));
    //                 double jac = length/2;
    //                 if(theCondition->type == NEUMANN_T){
    //                     double g[2] = {val*(X2 - X1)/length, val*(Y2 - Y1)/length}; // Calcul du vecteur tangent
    //                     B[2*(2*(Dom->elem[j]))] += jac*g[0] *phi[k];
    //                     B[2*(2*(Dom->elem[j])) + 1] += jac*g[1] *phi[k];
    //                     B[2*(2*(Dom->elem[j]) + 1)] += jac*g[0]* phi[k];
    //                     B[2*(2*(Dom->elem[j]) + 1) + 1] += jac*g[1] *phi[k];
    //                 }
    //                 if(theCondition->type == NEUMANN_N){
    //                     double g[2] = {val*(Y1 - Y2)/length, val*(X2 - X1)/length}; // Calcul du vecteur normal
    //                     B[2*(2*(Dom->elem[j]))] += jac*g[0] *phi[k];
    //                     B[2*(2*(Dom->elem[j])) + 1] += jac*g[1]*phi[k];
    //                     B[2*(2*(Dom->elem[j]) + 1)] += jac*g[0] *phi[k];
    //                     B[2*(2*(Dom->elem[j]) + 1) + 1] += jac*g[1]*phi[k];
    //                 }
    //                 if(theCondition->type == NEUMANN_X){
    //                     B[2*(2*(Dom->elem[j]))] += jac*val*phi[k];
    //                     B[2*(2*(Dom->elem[j]) + 1)] += jac*val*phi[k];
    //                 }
    //                 if(theCondition->type == NEUMANN_Y){
    //                     B[2*(2*(Dom->elem[j])) + 1] += jac*val*phi[k];
    //                     B[2*(2*(Dom->elem[j]) + 1) + 1] += jac*val*phi[k];
    //                 }
    //             }
    //         }
    //     }
    // }
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