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

double *femElasticitySolve(femProblem *theProblem)
{
    femFullSystem  *theSystem = (femFullSystem  *)theProblem->solver->solver;
    femIntegration *theRule = theProblem->rule;
    femDiscrete    *theSpace = theProblem->space;
    femGeo         *theGeometry = theProblem->geometry;
    femNodes       *theNodes = theGeometry->theNodes;
    femMesh        *theEdges = theGeometry->theEdges;
    femMesh        *theMesh = theGeometry->theElements;
    int            *renumber = theGeometry->theNodes->number;
    
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
    // printf("\n=====================================\n");
    // printf("    Constrained number_1 : \n");
    // for (size_t i = 0; i < theSystem->size; i++)
    // {
    //     printf(" B : %f \t ", B[i]);
    // }
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        for (j=0; j < nLocal; j++) {
            map[j]  = theMesh->elem[iElem*nLocal+j];
            mapX[j] = 2*map[j];
            mapY[j] = 2*map[j] + 1;
            x[j]    = theNodes->X[map[j]];
            y[j]    = theNodes->Y[map[j]];} 
        
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            // printf("xsi = %f\n",xsi) ;  
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            // femDiscretePrint(theSpace);  
            // printf("phi = %f\n",phi[0]) ;
            double dxdxsi = 0.0;
            double dxdeta = 0.0;
            double dydxsi = 0.0; 
            double dydeta = 0.0;
            for (i = 0; i < theSpace->n; i++) {  
                dxdxsi += x[i]*dphidxsi[i];       
                dxdeta += x[i]*dphideta[i];   
                dydxsi += y[i]*dphidxsi[i];   
                dydeta += y[i]*dphideta[i]; 
                r      += x[i]*phi[i];
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
        // printf("\n=====================================\n");
        // printf("    Constrained number_2 : \n");
        // for (size_t i = 0; i < theSystem->size; i++)
        // {
        //     printf(" B : %f \t", B[i]);
        // }
    
    for  (int iEdge=0; iEdge < theEdges->nElem; iEdge++)
    {
        
        if (theConstrainedEdges[iEdge] != -1) {
            // printf("iEdges : %d  \n", iEdge);  
            double value = theProblem->conditions[theConstrainedEdges[iEdge]]->value;
            int type = theProblem->conditions[theConstrainedEdges[iEdge]]->type;
            femSystemConstrainNEUMANN(theProblem, theSystem, iEdge, value, type,theProblem->planarStrainStress);
        }
    }
    int *theConstrainedNodes = theProblem->constrainedNodes;     
    for (int i=0; i < theSystem->size; i++) {
        if (theConstrainedNodes[i] != -1) {
            double value = theProblem->conditions[theConstrainedNodes[i]]->value;
            femSystemConstrainDIRICHLETXY(theSystem,i,value,theProblem->conditions[theConstrainedNodes[i]]->type); 
        }
    }
        printf("\n=====================================\n");
        printf("    Constrained number_3     : \n");
        for (size_t i = 0; i < theSystem->size; i++)
        {
            printf(" B : %f \t", B[i]);
        }
    double *sol ;
    double **L ;
    // switch (theSolver->type)
    // {
    // case FEM_Cholesky:
    //     L  = choleskyDecomposition(A,theSystem->size);  
    //     sol = solvecholesky(L,theSystem->B,theSystem->size) ;
    //     break;
    // case FEM_FULL : 
        sol = femFullSystemEliminate(theSystem);
        // break;
    // default:
    //     Error("Unexpected solver option");
    //     break;
    // }   
    return sol; 
}