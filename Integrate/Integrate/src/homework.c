#include <stdio.h>
#include <math.h>
#include "glfem.h"

static double integrate(double x[3], double y[3], double (*f)(double, double))
{
    double xis[3] =     { 1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0 };
    double es[3] =      { 1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0 };
    double weights[3] = { 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0 };
    
    double jac = (x[0] - x[1]) * (y[0] - y[2]) - (x[0] - x[2]) * (y[0] - y[1]);
    if (jac < 0.0)
        jac = -jac;

    double xy[2][3];
    double I = 0;
    for (int i = 0; i < 3; i++)
    {
        xy[0][i] = x[0] * xis[i] + x[1] * es[i] + x[2] * (1.0 - xis[i] - es[i]);
        xy[1][i] = y[0] * xis[i] + y[1] * es[i] + y[2] * (1.0 - xis[i] - es[i]);
        I += weights[i] * f(xy[0][i], xy[1][i]) * jac;
    }

    glfemSetColor(GLFEM_BLACK);
    glfemDrawElement(x, y, 3);

    glfemSetColor(GLFEM_BLUE);
    glfemDrawNodes(x, y, 3);

    glfemSetColor(GLFEM_RED);
    glfemDrawNodes(xy[0], xy[1], 3);

    return I;
};


double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{
    if(n<=0){
        return integrate(x,y,f);
    }else{
        double x_divide[4][3] = {{x[0],(x[1]+x[0])/2.0,(x[2]+x[0])/2.0},{(x[1]+x[0])/2.0,x[1],(x[2]+x[1])/2.0},{(x[2]+x[0])/2.0,(x[2]+x[1])/2.0,x[2]},{(x[1]+x[0])/2.0,(x[2]+x[1])/2.0,(x[0]+x[2])/2.0}};
        double y_divide[4][3] = {{y[0],(y[1]+y[0])/2.0,(y[2]+y[0])/2.0},{(y[1]+y[0])/2.0,y[1],(y[2]+y[1])/2.0},{(y[2]+y[0])/2.0,(y[2]+y[1])/2.0,y[2]},{(y[1]+y[0])/2.0,(y[2]+y[1])/2.0,(y[0]+y[2])/2.0}} ;
        double I = 0.0 ; 
        for (size_t i = 0; i < 4; i++)
        {
            I += integrateRecursive(x_divide[i],y_divide[i],f,n-1);
            
        }
        return I;
    }
}