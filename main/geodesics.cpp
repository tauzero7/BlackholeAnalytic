/**
 * @file    geodesics.cpp
 * @author  Thomas Mueller
 * 
 * @brief   Null geodesics in the Schwarzschild spacetime.
 * 
 * 
 * @copyright (c) 2017 Thomas Mueller
 *
 * This file is part of the BlackholeAnayltic.
 * 
 * BlackholeAnayltic is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * BlackholeAnayltic is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with BlackholeAnayltic.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include "bhAnalytic.h"


int main(int argc, char* argv[]) {
    double xi, ksi, ksiCrit, xf, phi1, phi2, xmin, phimin, phiinf;
    
    gsl_set_error_handler_off();
    
    printf("===============================================================\n");
    
    // A geodesic with starting position xi, final position xf, and initial
    // direction ksi covers an azimuthal angle phi1 (and phi2)
    xi = 0.1;
    xf = 0.2;
    ksi = 15.0 * DEG_TO_RAD;
    getAzimuthFromEndpoints(xi, xf, ksi, phi1, phi2);
    calcKsiCrit(xi, ksiCrit);
    
    printf("Start position at xi = %f with ksi = %f° has final position xf = %f for\n",
        xi, ksi * RAD_TO_DEG, xf);
    printf("\t phi1 = %12.8f° ... phi1 = %f\n", phi1 * RAD_TO_DEG, phi1);
    printf("\t phi2 = %12.8f° ... phi2 = %f\n\n", phi2 * RAD_TO_DEG, phi2);
    printf("Critical angle reads ksiCrit = %.12f°\n\n", ksiCrit * RAD_TO_DEG);
    
    // The point of closes approach for this geodesic reads...
    xmin = getXmin(xi, ksi);
    phimin = getPhiMin(xi, xmin, ksi);
    printf("The point of closest approach xmin = %f is reached at phimin = %f\n\n",
        xmin, phimin);
        
    // And in the limit, the geodesic will reach x=0 for phi = ...
    phiinf = phiInfty(xi, &ksi);
    printf("In the limit r->inf (x=0), the geodesic's azimuth reads phiInfty = %f° ...  phiInfty = %f\n\n",
        phiinf * RAD_TO_DEG, phiinf);
        
        
    // If a geodesic starts at xi with initial angle ksi, the radials position xf
    // for given phi reads...
    xi = 0.1;
    ksi = 26.0 * DEG_TO_RAD;
    phi1 = PI;
    xf = calcXf(xi, ksi, phi1);
    printf("A geodesic starting at xi = %f with ksi = %f° will end at xf = %f, phif = %f°\n", 
        xi, ksi * RAD_TO_DEG, xf, phi1 * RAD_TO_DEG);
    printf("\t rf = rs/xf = 1/xf = %f\n\n", 1.0/xf);
    
    //getAzimuthFromEndpoints(xi, xf, ksi, phi1, phi2);
    //printf("\t phi1 = %.12f  phi2 = %.12f\n", phi1, phi2);
    
    // A geodesic starting at xi,phi=0 and ending at xi,phi=pi must start with
    // and initial direction ksi = ...
    xi = 0.2;
    ksi = getKsiFromEndpoints(xi, xi, PI);
    printf("A geodesic starting at xi = %f,phi=0 and ending at xi,phi=pi must start with\n\t ksi = %.10f° (%.10f)\n\n", 
        xi, ksi * RAD_TO_DEG, ksi);
    
    printf("===============================================================\n");
    return 0;
}
