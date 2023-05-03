#include "fem.h"

double femTriangleInterpolate(double u[3], double xsi, double eta)
{
    return u[0] * xsi + u[1] * eta + u[2] * (1.0 - xsi - eta);
}

#ifndef NORHOSTEEL
double inertiaSteelRho()
{
    double rho = 7700;
    return rho;
}
#endif

#ifndef NOINTEGRATE
double inertiaIntegrate(femMesh *theMesh, femIntegration *theRule, double rho)
{
    double jacobian;
    double I = 0;
    double xLoc, x[3];
    double yLoc, y[3];
    int i, j, k;
    
    for(i = 0; i < theMesh->nElem; i++) {
        for(j = 0; j < 3; j++) {
            x[j] = theMesh->X[theMesh->elem[3*i+j]] ;
            y[j] = theMesh->Y[theMesh->elem[3*i+j]] ; }
        jacobian = fabs((x[1] - x[0]) * (y[2] - y[0]) - (y[1] - y[0]) * (x[2] - x[0]));
        for(k = 0; k < theRule->n; k++) {
            xLoc = femTriangleInterpolate(x, theRule->xsi[k], theRule->eta[k]);
            yLoc = femTriangleInterpolate(y, theRule->xsi[k], theRule->eta[k]);
            I += (xLoc*xLoc + yLoc*yLoc) * jacobian * theRule->weight[k]; }}
     return I * rho * 1.0e-13;
}
#endif

#ifndef NOREAD
femMesh *inertiaMeshRead(const char *filename)
{
    femMesh *theMesh = malloc(sizeof(femMesh));

    int i,trash,*elem;
    
    FILE* file = fopen(filename,"r");
    if (file == NULL) Error("No mesh file !");

    ErrorScan(fscanf(file, "Number of nodes %d \n", &theMesh->nNode));
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
    for (i = 0; i < theMesh->nNode; ++i) {
        ErrorScan(fscanf(file,"%d : %le %le \n",&trash,&theMesh->X[i],&theMesh->Y[i])); }

  
    ErrorScan(fscanf(file, "Number of triangles %d \n", &theMesh->nElem));
    theMesh->elem = malloc(sizeof(int)*3*theMesh->nElem);
    theMesh->nLocalNode = 3;
    for (i = 0; i < theMesh->nElem; ++i) {
        elem = &(theMesh->elem[i*3]);
        ErrorScan(fscanf(file,"%d : %d %d %d\n", &trash,&elem[0],&elem[1],&elem[2])); }
 
    fclose(file);
    return theMesh;
}
#endif

#ifndef NOFREE
void inertiaMeshFree(femMesh *theMesh)
{
    free(theMesh->X);
    free(theMesh->Y);
    free(theMesh->elem);
    free(theMesh);
}
#endif