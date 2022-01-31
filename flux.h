#ifndef FLUX_H
#define FLUX_H

#include <math.h>
#include <algorithm>
#include <iostream>

void flux(double* UL, double* UR, double* n, double* f, double *ss)
{
    double g = 9.8;

    //Left state
    double hL = UL[0];
    double uL = UL[1]/hL;
    double vL = UL[2]/hL;
    double unL = uL*n[0] + vL*n[1];
    
    //Left flux
    double FL[3];
    FL[0] = hL*unL;
    FL[1] = hL*uL*unL + 0.5*g*hL*hL*n[0];
    FL[2] = hL*vL*unL + 0.5*g*hL*hL*n[1];
    
    //Right state
    double hR = UR[0];
    double uR = UR[1]/hR;
    double vR = UR[2]/hR;
    double unR = uR*n[0] + vR*n[1];
    
    //Right flux
    double FR[3];
    FR[0] = hR*unR;
    FR[1] = hR*uR*unR + 0.5*g*hR*hR*n[0];
    FR[2] = hR*vR*unR + 0.5*g*hR*hR*n[1];
    
    //Difference in states
    double dh = hR - hL;
    double dhu = UR[1] - UL[1];
    double dhv = UR[2] - UL[2];

    /*
    std::cout << "FL = (" << FL[0] << ", " << FL[1] << ", " << FL[2] << ")" << std::endl;
    std::cout << "FR = (" << FR[0] << ", " << FR[1] << ", " << FR[2] << ")" << std::endl;
    std::cout << "d = (" << dh << ", " << dhu << ", " << dhv << ")" << std::endl;
    */

    //Roe average
    double ha = 0.5*(hL + hR);
    double ua = 0.5*(uL + uR);
    double va = 0.5*(vL + vR);
    double c = sqrt(g*ha);
    double un = ua*n[0] + va*n[1];

    //Eigenvalues
    double lambda[3];
    lambda[0] = un;
    lambda[1] = un - c;
    lambda[2] = un + c;

    /*
    std::cout << "average = (" << ha << ", " << ua << ", " << va << ")" << std::endl;
    std::cout << "c = " << c << std::endl;
    std::cout << "un = " << un << std::endl;
    std::cout << "lambda = (" << lambda[0] << ", " << lambda[1] << ", " << lambda[2] << ")" << std::endl;
    */

    //Entropy fix
    double epsilon = 0.01*c;
    for (int i = 0; i < 3; i++)
    {
        if (fabs(lambda[i]) < epsilon)
        {
            lambda[i] = 0.5*(epsilon + lambda[i]*lambda[i]/epsilon);
        }
        lambda[i] = fabs(lambda[i]);
    }

    /*
    std::cout << "lambda = (" << lambda[0] << ", " << lambda[1] << ", " << lambda[2] << ")" << std::endl;
    */

    //Coefficients in the stabilization term
    double s1 = 0.5*(lambda[1] + lambda[2]);
    double s2 = 0.5*(lambda[1] - lambda[2]);
    double A1 = s2/c*un + s1;
    double B1 = -s2/c*(dhu*n[0] + dhv*n[1]);
    double A2 = lambda[0];
    double B2 = (A1 - A2)*dh + B1;
    double C2 = ((lambda[0] - s1)*un - s2*c)*dh + (s1 - lambda[0])*(dhu*n[0] + dhv*n[1]);

    //Flux assembly
    f[0] = 0.5*(FL[0] + FR[0]) - 0.5*(A1*dh + B1);
    f[1] = 0.5*(FL[1] + FR[1]) - 0.5*(A2*dhu + B2*ua + C2*n[0]);
    f[2] = 0.5*(FL[2] + FR[2]) - 0.5*(A2*dhv + B2*va + C2*n[1]);
    *ss = std::max(std::max(lambda[0], lambda[1]), lambda[2]);

}

#endif