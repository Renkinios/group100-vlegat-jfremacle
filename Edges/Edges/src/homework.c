#include "fem.h"



# ifndef NOEXPAND
void edgesExpand(femEdges *theEdges)
{
    int *x = theEdges ->mesh-> elem ;
    theEdges ->nBoundary = theEdges ->nEdge ;
    // creation arrete 
    for (size_t i = 0; i < theEdges->mesh->nElem; i++)
    {
        theEdges->edges[3*i].node[0] = x[3*i] ;
        theEdges->edges[3*i].node[1] = x[3*i+1] ;
        theEdges->edges[3*i].elem[0] = i ; 
        theEdges->edges[3*i].elem[1] = -1 ; 
        theEdges->edges[3*i+1].node[0] = x[3*i+1] ;
        theEdges->edges[3*i+1].node[1] = x[3*i+2] ;
        theEdges->edges[3*i+1].elem[0] = i ; 
        theEdges->edges[3*i+1].elem[1] = -1 ; 
        theEdges->edges[3*i+2].node[0] = x[3*i+2] ;
        theEdges->edges[3*i+2].node[1] = x[3*i] ;
        theEdges->edges[3*i+2].elem[0] = i ; 
        theEdges->edges[3*i+2].elem[1] = -1 ; 
    }
}

# endif
# ifndef NOSORTù

void edgesSort(femEdges *theEdges)
{
	qsort(theEdges->edges, theEdges->nEdge, sizeof(femEdge), edgesCompare);
}

# endif
# ifndef NOCOMPARE

int edgesCompare(const void* e0, const void *e1)
{   
    femEdge* edge0 = (femEdge*) e0;
    femEdge* edge1 = (femEdge*) e1;

    int minNode0 = fmin(edge0->node[0], edge0->node[1]),minNode1 = fmin(edge1->node[0], edge1->node[1]);
    int maxNode0 = fmax(edge0->node[0], edge0->node[1]),maxNode1 = fmax(edge1->node[0], edge1->node[1]);

    if (minNode0 < minNode1) 
        return 1 ;
    if (minNode0 > minNode1)
        return -1;
    if (maxNode0 < maxNode1)
        return 1;
    if (maxNode0 > maxNode1)
        return -1;

    return 0;
}

# endif
# ifndef NOSHRINK

void edgesShrink(femEdges *theEdges)
{
    // Ici commence votre contribution a faire :-)

    int n = 0;          // Nouveau nombre total de segments : A MODIFIER
    int nBoundary = 0;  // Nombre de segments frontieres : A MODIFIER

    femEdge *edges = theEdges->edges;
    for (size_t i = 0; i < theEdges->nEdge; i++){
        if (i==theEdges->nEdge -1 || (edgesCompare((const void *) &edges[i], (const void *) &edges[i+1]) != 0)){
            edges[n] = edges[i];
            nBoundary++;
        }
        else{
            edges[i].elem[1] = edges[i+1].elem[0] ;
            edges[n] = edges[i];
            i++;
        }
        n++;
    }
    theEdges->edges = realloc(theEdges->edges, n * sizeof(femEdge));
    theEdges->nEdge = n;
    theEdges->nBoundary = nBoundary;
}

# endif
# ifndef NOBOUNDARYLENGTH

double edgesBoundaryLength(femEdges *theEdges)
{
    double length_boudary = 0.0; 
    for (size_t i = 0; i < theEdges->nEdge; i++) 
    {
        if (theEdges->edges[i].elem[1] == -1)
        {
            int node_1 = theEdges->edges[i].node[0] , node_2 = theEdges->edges[i].node[1] ; 
            double delta_X = theEdges->mesh ->X[node_1] - theEdges->mesh ->X[node_2], delta_Y = theEdges->mesh ->Y[node_1] - theEdges->mesh ->Y[node_2] ; 
            length_boudary += sqrt(delta_X*delta_X+delta_Y*delta_Y) ; 
        }
    }
    return length_boudary ; 
}

# endif
