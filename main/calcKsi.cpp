
#include <iostream>
#include "bhAnalytic.h"


int main(int argc, char* argv[]) {
    double xi, ksi, ksiCrit, xf, phi1, phi2, xmin, phimin, phiinf;
    
    gsl_set_error_handler_off();
    
    // A geodesic starting at xi,phi=0 and ending at xf,phif must start with
    // and initial direction ksi = ...
       
    double sn,cn,dn;
    double u = atof(argv[1]);
    double m = atof(argv[2]);
    gsl_sf_elljac_e(u, m, &sn, &cn, &dn);
    printf("%f %f %f\n", sn, cn, dn);
    return 0;
       
    double rs = 2.0;
    xi = rs / 10.0;    
    xf = rs / atof(argv[1]);    
    phi2 = radians(atof(argv[2]));

    double ksi_lo = radians(atof(argv[3]));
    double ksi_hi = radians(atof(argv[4]));
    ksi = getKsiFromEndpoints(xi, xf, phi2, ksi_lo, ksi_hi);
    printf("ksi = %.10f° (%.10f)\n\n",  ksi * RAD_TO_DEG, ksi);

    return 0;
}
