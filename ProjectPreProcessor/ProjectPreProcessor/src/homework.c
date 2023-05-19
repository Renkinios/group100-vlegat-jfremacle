#include "fem.h"

//
// Ici, vous pouvez définir votre géométrie :-)
//  (1) Raffiner intelligemment.... (yes )
//  (2) Construire la geometrie avec OpenCascade 
//  (3) Construire la geometrie avec les outils de GMSH
//  (4) Obtenir la geometrie en lisant un fichier .geo de GMSH


double geoSize(double x, double y) {

    femGeo* theGeometry = geoGetGeometry();
    return theGeometry->h * (1.0 - 0.5 * x);
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
