#include"fem.h"


double convergenceSource(double x, double y)
{
    double f = 2*(-800*x*y*(x - 1)*(y - 1)*(5*M_SQRT2*(x + y) - 8) + (x*(x - 1) 
      + y*(y - 1))*pow(4*pow(5*M_SQRT2*(x + y) - 8, 2) + 1, 2)*atan(10*M_SQRT2*(x + y) - 16) 
      + 10*M_SQRT2*(4*pow(5*M_SQRT2*(x + y) - 8,2) + 1)*(x*y*(x - 1) 
      + x*y*(y - 1) + x*(x - 1)*(y - 1) 
      + y*(x - 1)*(y - 1)))/pow(4*pow(5*M_SQRT2*(x + y) - 8, 2) + 1,2);
    return -f; 
}

double convergenceSoluce(double x, double y, double *u)
{
    u[0] = x*y*(1-x)*(1-y)*atan(20*(x+y)/M_SQRT2 - 16);
    u[1] = y*(y - 1)*(10*M_SQRT2*x*(x - 1)
      + (2*x - 1)*(4*pow(5*M_SQRT2*(x + y) - 8, 2) + 1)*atan(10*M_SQRT2*(x + y) - 16))
      /(4*pow(5*M_SQRT2*(x + y) - 8, 2) + 1);
    u[2] = x*(x - 1)*(10*M_SQRT2*y*(y - 1) 
      + (2*y - 1)*(4*pow(5*M_SQRT2*(x + y) - 8, 2) + 1)*atan(10*M_SQRT2*(x + y) - 16))
      /(4*pow(5*M_SQRT2*(x + y) - 8, 2) + 1);
    return u[0];
}

double convergenceEstimateRate(double *errors, int n, double ratio)
{

     
     double rate = 0;
     for (int i=0; i < n-1; i++)  
         rate += log(errors[i]/errors[i+1])/log(2);
     return rate/(n-1);
}


void femDiffusionComputeError(femDiffusionProblem *theProblem, 
                                    double(*soluce)(double,double,double*))
{
    femMesh *theMesh = theProblem->mesh;
    femIntegration *theRule = theProblem->rule;
    femDiscrete *theSpace = theProblem->space;

    if (theSpace->n > 4) Error("Unexpected discrete space size !");    
    double Xloc[4],Yloc[4],phi[4],dphidxsi[4],dphideta[4];
    double Uloc[4],Utilt[4],uRef[3];
    int iEdge,iElem,iInteg,i,j,map[4],ctr[4];
    
    
    theProblem->errorSoluceL2 = 0.0;
    theProblem->errorSoluceH1 = 0.0;
    theProblem->errorInterpolationL2 = 0.0;
    theProblem->errorInterpolationH1 = 0.0;
    
    for (iElem = 0; iElem < theMesh->nElem; iElem++) {
        femDiffusionMeshLocal(theProblem,iElem,map,ctr,Xloc,Yloc,Uloc); 
        for (i = 0; i < theSpace->n; i++)
            Utilt[i] = soluce(Xloc[i],Yloc[i],uRef);
        for (iInteg=0; iInteg < theRule->n; iInteg++) {    
            double xsi    = theRule->xsi[iInteg];
            double eta    = theRule->eta[iInteg];
            double weight = theRule->weight[iInteg];  
            femDiscretePhi2(theSpace,xsi,eta,phi);
            femDiscreteDphi2(theSpace,xsi,eta,dphidxsi,dphideta);
            double x = 0;
            double y = 0;
            double u = 0;
            double utilt = 0;
            double dudx = 0;
            double dudy = 0;
            double dutiltdx = 0;
            double dutiltdy = 0;
            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0; 
            double dydeta = 0;
            for (i = 0; i < theSpace->n; i++) {  
                x      += Xloc[i]*phi[i];
                y      += Yloc[i]*phi[i]; 
                u      += Uloc[i]*phi[i]; 
                utilt  += Utilt[i]*phi[i]; 
                dxdxsi += Xloc[i]*dphidxsi[i];       
                dxdeta += Xloc[i]*dphideta[i];   
                dydxsi += Yloc[i]*dphidxsi[i];   
                dydeta += Yloc[i]*dphideta[i]; }
            double jac = fabs(dxdxsi * dydeta - dxdeta * dydxsi);
            for (i = 0; i < theSpace->n; i++) {    
                double dphidx = (dphidxsi[i] * dydeta - dphideta[i] * dydxsi) / jac;       
                double dphidy = (dphideta[i] * dxdxsi - dphidxsi[i] * dxdeta) / jac; 
                dudx += Uloc[i]*dphidx;
                dudy += Uloc[i]*dphidy;
                dutiltdx += Utilt[i]*dphidx;
                dutiltdy += Utilt[i]*dphidy;}            
            
            soluce(x,y,uRef);
            double e     = (u - uRef[0]);          
            double dedx  = (dudx - uRef[1]);
            double dedy  = (dudy - uRef[2]);           
            theProblem->errorSoluceL2 += jac * e*e * weight; 
            theProblem->errorSoluceH1 += jac * (e*e + dedx*dedx + dedy*dedy) * weight;
            
            e    = (utilt - uRef[0]);
            dedx = (dutiltdx - uRef[1]);
            dedy = (dutiltdy - uRef[2]);           
            theProblem->errorInterpolationL2 += jac * e*e * weight;              
            theProblem->errorInterpolationH1 += jac * (e*e + dedx*dedx + dedy*dedy) * weight;}}
    
     theProblem->errorSoluceL2 = sqrt(theProblem->errorSoluceL2);
     theProblem->errorSoluceH1 = sqrt(theProblem->errorSoluceH1);
     theProblem->errorInterpolationL2 = sqrt(theProblem->errorInterpolationL2);
     theProblem->errorInterpolationH1 = sqrt(theProblem->errorInterpolationH1);
}



femMesh *femMeshCreateBasicSquare(int n) 
{
    femMesh *theMesh = malloc(sizeof(femMesh));
    
    theMesh->nNode = (n+1)*(n+1);
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
    
    double h = 1.0 / n;
    int iNode = 0;    
    for (int i = 0; i < n+1; i++) {
        double x = i*h;
        for (int j = 0; j < n+1; j++) {
            double y = j*h;
            theMesh->X[iNode] = x;
            theMesh->Y[iNode] = y; 
            iNode++; }}
            
    theMesh->nElem = n*n*2;
    theMesh->nLocalNode = 3;
    theMesh->elem = malloc(sizeof(int)*3*theMesh->nElem);
    int iElem = 0, shift;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int shift = i*(n+1)+j;
            int *elem = &(theMesh->elem[iElem*3]);
            elem[0] = shift;
            elem[1] = shift+n+2;
            elem[2] = shift+1;
            elem[3] = shift+n+1;
            elem[4] = shift+n+2;
            elem[5] = shift;
            iElem = iElem+2; }}
    
            
    theMesh->number = malloc(sizeof(int)*theMesh->nNode); 
    for (int i = 0; i < theMesh->nNode; i++) 
          theMesh->number[i] = i;     
    return theMesh;
}
