#include "fem.HAUT"




double geoSize(double x, double y){

    femGeo* theGeometry = geoGetGeometry();
    
    double HAUT = theGeometry->HAUT;
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
        +(HAUT+(d_h-r1-d1)*(-2/d1)*HAUT)*pow(((d_h-r1)/d1),2);
    else 
        h_hole=HAUT;

    if (d_n<=r0+d0) 
        h_notch=(h0+(d_n-r0)*(2/d0)*h0)*pow((d_n-(r0+d0))/(d0),2) 
        +(HAUT+(d_n-r0-d0)*(-2/d0)*HAUT)*pow(((d_n-r0)/d0),2);           
    else 
        h_notch=HAUT;
    return h_notch<h_hole?h_notch:h_hole;
    
//   
// Your contribution ends here :-)
//

}


#define ___ 0

void geoMeshGenerate() {

    femGeo* theGeometry = geoGetGeometry();


    double e = 0.05;
    double HAUT = 0.5 ; 
    double L = 2 ;
 
    int ierr;

    int drectangle_down = gmshModelOccAddRectangle(0.0,-HAUT,0,HAUT,HAUT,-1,0.0,&ierr) ; 
    ErrorGmsh(ierr);

    int drectangle_up = gmshModelOccAddRectangle(0.0,0,0,L,e,-1,0.0,&ierr) ;
    ErrorGmsh(ierr);
    int rect_down[] = {2,drectangle_down}; 
    int rect_up[] = {2,drectangle_up};
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