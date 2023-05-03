#include "fem.h"




double geoSize(double x, double y){

    femGeo* theGeometry = geoGetGeometry();
    
    double h = theGeometry->h;
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    double h0 = theGeometry->hNotch;
    double d0 = theGeometry->dNotch;
  
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
    double h1 = theGeometry->hHole;
    double d1 = theGeometry->dHole;

    double d_h= sqrt(pow((x1-x),2)+pow((y1-y),2));
    double d_n= sqrt(pow((x0-x),2)+pow((y0-y),2));
    double h_hole,h_notch,h_final;
    if (d_h<=r1+d1)
        h_hole=(h1+(d_h-r1)*(2/d1)*h1)*pow((d_h-(r1+d1))/(d1),2)
        +(h+(d_h-r1-d1)*(-2/d1)*h)*pow(((d_h-r1)/d1),2);
    else 
        h_hole=h;

    if (d_n<=r0+d0) 
        h_notch=(h0+(d_n-r0)*(2/d0)*h0)*pow((d_n-(r0+d0))/(d0),2) 
        +(h+(d_n-r0-d0)*(-2/d0)*h)*pow(((d_n-r0)/d0),2);           
    else 
        h_notch=h;
    return h_notch<h_hole?h_notch:h_hole;
    
//   
// Your contribution ends here :-)
//

}


#define ___ 0

void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();

    double w = theGeometry->LxPlate;
    double h = theGeometry->LyPlate;
     
    double x0 = theGeometry->xNotch;
    double y0 = theGeometry->yNotch;
    double r0 = theGeometry->rNotch;
    
    
    double x1 = theGeometry->xHole;
    double y1 = theGeometry->yHole;
    double r1 = theGeometry->rHole;
 
//
//  -1- Construction de la g�om�trie avec OpenCascade
//      On cr�e le rectangle
//      On cr�e les deux cercles
//      On soustrait les cercles du rectangle :-)
//
 
    int ierr;
    int idPlate = gmshModelOccAddRectangle(-w/2,- h/2,0, w, h, -1   ,0 ,&ierr);   
    ErrorGmsh(ierr);
    int idNotch = gmshModelOccAddDisk(x0, y0, 0,r0, r0 ,-1 ,NULL,0,NULL,0,&ierr); 
    ErrorGmsh(ierr);
    int idHole  =   gmshModelOccAddDisk(x1, y1, 0, r1, r1, -1,NULL,0,NULL,0,&ierr); 
    ErrorGmsh(ierr);
    
    int plate[] = {2,idPlate};
    int notch[] = {2,idNotch};
    int hole[]  = {2,idHole};
    gmshModelOccCut(plate,2,notch,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    ErrorGmsh(ierr);
    gmshModelOccCut(plate,2,hole,2,NULL,NULL,NULL,NULL,NULL,-1,1,1,&ierr);
    ErrorGmsh(ierr);
 
//
//  -2- D�finition de la fonction callback pour la taille de r�f�rence
//      Synchronisation de OpenCascade avec gmsh
//      G�n�ration du maillage (avec l'option Mesh.SaveAll :-)
                  
   
    geoSetSizeCallback(geoSize);
                                  
    gmshModelOccSynchronize(&ierr);       
    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
    gmshModelMeshGenerate(2, &ierr);  
       
//
//  Generation de quads :-)
//
//    gmshOptionSetNumber("Mesh.SaveAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.RecombineAll", 1, &ierr);
//    gmshOptionSetNumber("Mesh.Algorithm", 8, &ierr);  chk(ierr);
//    gmshOptionSetNumber("Mesh.RecombinationAlgorithm", 1.0, &ierr);  chk(ierr);
//    gmshModelGeoMeshSetRecombine(2,1,45,&ierr);  chk(ierr);
//    gmshModelMeshGenerate(2, &ierr);  
   
 
//
//  Plot of Fltk
//
//   gmshFltkInitialize(&ierr);
//   gmshFltkRun(&ierr);  chk(ierr);
//
    
}