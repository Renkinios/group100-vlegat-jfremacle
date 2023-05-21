#include "fem.h"


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
            y[j]    = theNodes->Y[map[j]];
        } 
        
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
            r = 0.0;
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
            femSystemConstrainNEUMANN(theProblem, theSystem, iEdge, value, type,theProblem->planarStrainStress);
        }
    }
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
