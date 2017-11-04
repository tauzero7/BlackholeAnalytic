/**
 * @file    ring.cpp
 * @author  Thomas Mueller
 * 
 * @brief   Einstein rings in the Schwarzschild spacetime.
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
    double xi, ksi, ksiCrit, ksi1, ksi2, ksi3, xmin1, xmin2, xmin3;
    double phimin1, phimin2, phimin3, T1, T2, T3;
    
    gsl_set_error_handler_off();
    
    printf("===============================================================\n");
    
    // scaled observer position xi = rs/robs with rs=Schwarzschlid radius
    // and robs=distance of the observer to black hole.    
    // The critical angle corresponds to the black hole shadow.
    xi = 0.1;
    calcKsiCrit(xi, ksiCrit);    
    printf("An observer at xi = %f has a critical angle ksi_crit = %.12f° = %.12f (rad).\n\n", 
        xi, ksiCrit * RAD_TO_DEG, ksiCrit);
    
    
    // If an observer emits a flash of light and sees a first Einstein ring
    // with an opening angle of ksi has a distance xi to the black hole.
    ksi = 30.0 * DEG_TO_RAD;
    xi = getDistanceFromRing(ksi);
    printf("Observer distance for opening angle ksi = %f° reads xi = %.10f.\n\n", 
        ksi * RAD_TO_DEG, xi);
    
    
    // If an observer at distance xi to the black hole emits a flash of light,
    // s/he will see Einstein rings of order n=1,2,3 with opening angles ksi(n).
    // The corresponding null geodesic reaches a point of closest approach xmin 
    // at phi = n*180 and needs dT to return to the observer.
    xi = 0.1;
    printf("For an observer distance of xi = %f, the opening angle\n", xi);
    
    ksi1 = getKsiFromDistance(xi, 1);
    xmin1 = getXmin(xi, ksi1);
    phimin1 = getPhiMin(xi, xmin1, ksi1);    
    printf("\t ksi(n=1) = %14.12f°   xmin(n=1) = %.8f   phimin(n=1) = %f°\n", ksi1 * RAD_TO_DEG, xmin1, phimin1 * RAD_TO_DEG);
    
    ksi2 = getKsiFromDistance(xi, 2);
    xmin2 = getXmin(xi, ksi2);
    phimin2 = getPhiMin(xi, xmin2, ksi2);
    printf("\t ksi(n=2) = %14.12f°   xmin(n=2) = %.8f   phimin(n=2) = %f°\n", ksi2 * RAD_TO_DEG, xmin2, phimin2 * RAD_TO_DEG);
    
    ksi3 = getKsiFromDistance(xi, 3);
    xmin3 = getXmin(xi, ksi3);
    phimin3 = getPhiMin(xi, xmin3, ksi3);
    printf("\t ksi(n=3) = %14.12f°   xmin(n=3) = %.8f   phimin(n=3) = %f°\n\n", ksi3 * RAD_TO_DEG, xmin3, phimin3 * RAD_TO_DEG);
        
    printf("\t ksi(n=1) = %14.12f (rad)\n", ksi1);
    printf("\t ksi(n=2) = %14.12f (rad)\n", ksi2);
    printf("\t ksi(n=3) = %14.12f (rad)\n\n", ksi3);
    
    getKsiAndTFromDistance(xi, 1, ksi1, T1);
    printf("\t dT(n=1)/rs = %.8f\n", T1);
    
    getKsiAndTFromDistance(xi, 2, ksi2, T2);
    printf("\t dT(n=2)/rs = %.8f\n", T2);
    
    getKsiAndTFromDistance(xi, 3, ksi3, T3);
    printf("\t dT(n=3)/rs = %.8f\n", T3);
         
    printf("===============================================================\n");
    return 0;
}
