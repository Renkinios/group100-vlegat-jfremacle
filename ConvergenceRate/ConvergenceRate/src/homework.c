//Inspir� de devoirs de l'ann�e pass�e
#include"fem.h"

double convergenceSource(double x, double y){

//Trouve �quation � l'aide de SymPy en python
    
    double racine = 5*M_SQRT2*(x+y) - 8;
    double powerRac = 4*pow(racine, 2);
    double powPowerRac = pow(powerRac+1,2);
    double racPowerRac = 10*M_SQRT2*(powerRac + 1);
    double arctg = atan(10*M_SQRT2*(x+y)-16);
    double mail = 2*(-800*x*y*(x-1)*(y-1)*racine + (x*(x-1) + y*(y-1))* powPowerRac*arctg + racPowerRac*(x*y*(x-1) + x*y*(y-1) + x*(x-1)*(y-1) + y*(x-1)*(y-1))) / powPowerRac;
    return -mail;
}

double convergenceSoluce(double x, double y, double *u){
    double a = sqrt(2.0);
    double arctg = atan(10*a*(x+y)-16);
    double power = (4*pow(5*a*(x+y)-8, 2) + 1);
    u[0] = x*y*(1-x)*(1-y)*atan(20*(x+y)/a - 16);
    u[1] = y*(y - 1)*(10*a*x*(x - 1) + (2*x - 1)*power*arctg)/power;
    u[2] = x*(x - 1)* ( 10*a*y*(y - 1) + (2*y - 1)*power*arctg )/power;
    return u[0];
}


double convergenceEstimateRate(double *errors, int n, double ratio){
    double rate = 0;
    for (int i = 0; i < (n-1); i++){
        rate += log(errors[i]/errors[i+1]);
    }
    double err = rate/ ((n-1)*log(ratio)); //Moyenne sur n-1 mesures
    return err;
}


void femDiffusionComputeError(femDiffusionProblem *theProblem, 
                                    double(*soluce)(double,double,double*)){
    femMesh *theMesh = theProblem -> mesh;
    femIntegration *theRule = theProblem -> rule;
    femDiscrete *theSpace = theProblem -> space;
    
    int node = theMesh->nLocalNode;
    
    // Initialise normes erreur
    double *SL2 = &theProblem->errorSoluceL2;
    *SL2 = 0.0;
    double *SH1 = &theProblem->errorSoluceH1;
    *SH1 = 0.0;
    double *IL2 = &theProblem->errorInterpolationL2;
    *IL2 = 0.0;
    double *IH1 = &theProblem->errorInterpolationH1;
    *IH1 = 0.0;
    
    //Initialisation variable utile pour la suite
    double *solbuf = malloc(3*sizeof(double));
    double *dphi_dksi=malloc(node*sizeof(double));
    double *dphi_deta=malloc(node*sizeof(double));
    double *phi= malloc(node*sizeof(double));
    double x,y,uh,ut,duhdx,duhdy,dutdx,dutdy;
    double dx_dksi,dx_deta,dy_deta,dy_dksi;
    double dphi_dxLoc, dphi_dyLoc;
    double jacobian;
    double xLoc[4], yLoc[4], uLoc[4], map [4];
    
    if(theSpace ->n > 4) Error("Unexpected discrete space size");
    
    for(int i = 0; i < theMesh->nElem; i++){
        femDiffusionMeshLocal(theProblem, i, (int*)map, (int*)map, xLoc, yLoc, uLoc);
        for(int j = 0; j < theSpace->n; j++){
            theSpace->phi2(theRule->xsi[j], theRule->eta[j], phi);
            theSpace->dphi2dx(theRule->xsi[j], theRule->eta[j], dphi_dksi, dphi_deta);
            x = 0.0; y = 0.0; uh = 0.0; ut = 0.0; duhdx = 0.0; duhdy = 0.0; dutdx = 0.0; dutdy = 0.0; dx_dksi = 0.0; dx_deta = 0.0; dy_dksi = 0.0; dy_deta = 0.0;
            for(int k = 0; k < node; k++){
                dx_dksi += xLoc[k]*dphi_dksi[k];
                dx_deta += xLoc[k]*dphi_deta[k];
                dy_dksi += yLoc[k]*dphi_dksi[k];
                dy_deta += yLoc[k]*dphi_deta[k];
                uh += uLoc[k]*phi[k];
                ut += soluce(xLoc[k], yLoc[k],map)*phi[k];
                
                jacobian = dx_dksi*dy_deta - dx_deta*dy_dksi;
            }
            for(int l = 0; l < node; l++){
                
                jacobian = (jacobian >= 0)? (jacobian): (-jacobian);
                
                dphi_dxLoc = (dphi_dksi[l]*dy_deta - dphi_deta[l]*dy_dksi)/jacobian;
                dphi_dyLoc = (dphi_deta[l]*dx_dksi - dphi_dksi[l]*dx_deta)/jacobian;
                dutdx += soluce(xLoc[l], yLoc[l], map)*dphi_dxLoc;
                dutdy += soluce(xLoc[l], yLoc[l], map)*dphi_dyLoc;
                duhdx += uLoc[l]*dphi_dxLoc;
                duhdy += uLoc[l]*dphi_dyLoc;
                x += xLoc[l]*phi[l];
                y += yLoc[l]*phi[l];
            }
            soluce(x,y,solbuf);
            *SL2 += pow(solbuf[0] - uh, 2) * jacobian *theRule->weight[j];
            *SH1 += (pow(solbuf[0]-uh,2)+pow(solbuf[1]-duhdx,2)+pow(solbuf[2]-duhdy,2)) * jacobian * theRule->weight[j];
            *IL2 += pow(solbuf[0] - ut, 2) * jacobian * theRule->weight[j];
            *IH1 += (pow(solbuf[0]-ut,2)+pow(solbuf[1]-dutdx,2)+pow(solbuf[2]-dutdy,2)) * jacobian * theRule->weight[j];
        }
    }
    *SL2 = sqrt(*SL2);
    *SH1 = sqrt(*SH1);
    *IL2 = sqrt(*IL2);
    *IH1 = sqrt(*IH1);
    free(solbuf);
    free(dphi_dksi);
    free(dphi_deta);
    free(phi);
}




femMesh *femMeshCreateBasicSquare(int n){
    femMesh *theMesh = malloc(sizeof(femMesh));

    // Toujours le maillage �l�mentaire :-)
    // A modifier [begin]
    
    theMesh->nNode = (n+1)*(n+1);
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
    double *X = theMesh->X;
    double *Y = theMesh->Y;
    double Xh = (double) -1/n;
    double Yh;
    
    
    theMesh->nElem = 2*n*n;
    theMesh->nLocalNode = 3;
    theMesh->elem = malloc(sizeof(int)*3*theMesh->nElem);
    int *elem = theMesh->elem;
    int i = 0;
    int k = 0;
    
    for(int colonne = 0; colonne < n; colonne++){
        Yh = 0;
        Xh += (double) 1/n;
        for(int ligne = 0; ligne < n; ligne++){
            elem[i] = k;
            elem[i+1] = k + n + 2;
            elem[i+2] = k + 1;
            elem[i+3] = k + n + 1;
            elem[i+4] = k + n + 2;
            elem[i+5] = k;
            
            X[k] = Xh;
            X[k+1] = Xh;
            X[k+n+1] = Xh + (double) 1/n;
            X[k+n+2] = Xh + (double) 1/n;
            
            Y[k] = Yh;
            Y[k+1] = Yh + (double) 1/n;
            Y[k+n+1] = Yh;
            Y[k+n+2] = Yh + (double) 1/n;
            
            i += 6;
            k += 1;
            Yh += (double) 1/n;
        }
        k++;
    }

    theMesh->number = malloc(sizeof(int)*theMesh->nNode); 
    for (int i = 0; i < theMesh->nNode; i++) 
          theMesh->number[i] = i;     
    return theMesh;
}
