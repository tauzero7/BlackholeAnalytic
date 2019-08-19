
#include <iostream>
#include "bhAnalytic.h"


int main(int argc, char* argv[]) {
    double xi, ksi, ksiCrit, xf, phi1, phi2, xmin, phimin, phiinf;
    
    gsl_set_error_handler_off();
    
    // A geodesic starting at xi,phi=0 and ending at xi,phi=pi must start with
    // and initial direction ksi = ...
    
    ksi = radians(30.0);
    printf("%.12f\n", degrees(phiInfty(1e-7, &ksi)));
    return 0;
    
    
    double rs = 2.0;
    xi = rs / 4.0;
    calcKsiCrit(xi, ksiCrit);
    printf("ksi_crit = %.12f° ; %.12f°\n", degrees(ksiCrit), 180 - degrees(ksiCrit));
    printf("ksi_crit = %.12f ; %.12f\n", ksiCrit, PI - ksiCrit);


    xf = rs / atof(argv[1]);
    phi2 = atof(argv[2]);
    
    double ksi_lo = radians(atof(argv[3]));
    double ksi_hi = radians(atof(argv[4]));
    
    ksi = getKsiFromEndpoints(xi, xf, phi2, ksi_lo, ksi_hi);
    printf("ksi = %.10f° (%.10f)\n\n", 
        ksi * RAD_TO_DEG, ksi);


    //xf = calcXf(xi, (atof(argv[1])), atof(argv[2]));
    //printf("rf = %.12f\n", rs/xf);
    
    /*
    ksi = atof(argv[1]);
    phi2 = atof(argv[2]);
    
    xf = calcXf(xi, radians(ksi), phi2);
    std::cout << rs / xf << std::endl;
    */
    return 0;
}
