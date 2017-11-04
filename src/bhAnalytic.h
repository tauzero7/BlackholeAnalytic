/**
 * @file    bhAnalytic.h
 * @author  Thomas Mueller
 * 
 * @brief   Analytic methods for null geodesics in the Schwarzschild spacetime.
 * 
 *  Note that any initial position "ri" are given as inverse scaled values
 *  "xi = rs/ri" with Schwarzschild radius rs.
 * 
 *  For a detailed discussion of the analytic solution, have a look into
 * 
 *  [TM2008]
 *    T.Mueller, "Einstein rings as a tool for estimating distances and the mass 
 *       of a Schwarzschild black hole", Phys. Rev. D 77, 124042 (2008).
 *    DOI: 10.1103/PhysRevD.77.124042
 * 
 * 
 * Dependencies:
 *   GNU Scientific Library
 * 
 * 
 * @copyright (c) 2017 Thomas Mueller
 *
 * This file is part of the BlackholeAnalytic.
 * 
 * BlackholeAnalytic is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * BlackholeAnalytic is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with BlackholeAnalytic.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef BH_ANALYTIC_H
#define BH_ANALYTIC_H

#include <iostream>
#include <cmath>
#include <cassert>

#include "defs.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_ellint.h>
#include <gsl/gsl_sf_elljac.h>
#include <gsl/gsl_integration.h>

/**
 * @brief Use Brent method to find root.
 * @param function
 * @param params
 * @param x_lo
 * @param x_hi
 * @param x
 * @return status
 */
int gslBrent(double function(double, void*), void *params, 
    double x_lo, double x_hi, double &x);


/**
 * 
 */
double gslIntegral(double function(double, void*), void *params,
    double x_lo, double x_hi);


/**
 * @brief Calculate a² parameter depending on xi and ksi.
 * @param xi
 * @param ksi
 * @return a^2
 */
double calcAqua(const double xi, const double ksi);


/**
 * @brief Calculate critical angle depending on observer position xi.
 *   The critical angle is given with respect to the direction to the black hole.
 * 
 * @param xi    Observer position.
 * @param ksi   Reference to critical angle (in radians)
 * @return true if critical angle exists
 */
bool calcKsiCrit(const double xi, double &ksi);


/**
 * @brief Calculate initial angle ksi for a null geodesic starting at xi
 *   and ending at azimuth phif at infinity.
 * @param xi    Starting position.
 * @param phif  Azimuth angle for end point at infinity.
 * @return ksi (in radians).
 */
double calcKsiPhiInfty(const double xi, const double phif);


/**
 * @brief  Calculate final position ...
 * @param xi
 * @param ksi
 * @param phi
 * @return
 */
double calcXf(const double xi, const double ksi, const double phi);


/**
 * 
 */
double calcXf_impl(double ksi, void *params);


/**
 * @brief Calculate limiting observer position xi.
 *    If an Einstein ring occurs for an observation angle ksi, the observer
 *    must have a minimum distance xi < xiLimit to the black hole. Otherwise,
 *    the light ray would come out of the black hole.
 * 
 * @param ksi  Angle of Einstein ring (radians).
 * @return
 */
double calcXiLimit(const double ksi);


/**
 * @brief
 * @param xi
 * @param xf
 * @param ksi
 * @param phi1
 * @param phi2
 * @return number of phis (0,1,2)
 */
int getAzimuthFromEndpoints(const double xi, const double xf, const double ksi,
    double &phi1, double &phi2);


/**
 * @brief Calculate maximum distance of a null geodesic to the black hole for
 *   the starting position xi closer than the photon orbit and direction ksi.
 * @param xi
 * @param ksi
 * @return
 */
double getXmax(const double xi, const double ksi);


/**
 * @brief Calculate minimum distance of a null geodesic to the black hole for
 *   the starting position xi and initial direction ksi.
 * @param xi
 * @param ksi
 * @return Scaled inimum distance xmin = rs/rmin; equals 0 if there is no 
 *   minimum distance.
 */
double getXmin(const double xi, const double ksi);


/**
 * @brief
 * @param aqua
 * @param x1
 * @param x2
 * @param x3
 */
void calcXs(const double aqua, double &x1, double &x2, double &x3);


/**
 * @brief Calculate observer distance from black hole depending on the opening 
 * angle ksi of the first Einstein ring produced by a flash of light.
 * @param ksi  Opening angle of first Einstein ring (in radians).
 * @return scaled distance xi
 */
double getDistanceFromRing(const double ksi);


/**
 * @brief Calculate opening angle ksi of Einstein ring due to a flash of light
 *   by the observer depending on the observer's scaled distance xi and the order
 *   of the flash n.
 * @param xi  Scaled observer distance to black hole.
 * @param n   Order of Einstein ring (n=1,2,3,...)
 * @return opening angle ksi (in radians).
 */
double getKsiFromDistance(const double xi, const int n);

void getKsiAndTFromDistance(const double xi, const int n, double &ksi, double &T);


/**
 * @brief
 */
double getKsiFromEndpoints(const double xi, const double xf, const double phi);


/**
 * @brief Calculate azimuth angle phi for a null geodesic at the point of
 *   closest approach to the black hole.
 * @param xi
 * @param xmin
 * @param ksi
 * @return azimuth angle (in radians)
 */
double getPhiMin(const double xi, const double xmin, const double ksi);


/**
 * @brief Differential equation for coordinate time t wrt. x.
 * @param x
 * @param params
 * @return
 */
double dtdx(double x, void *params);


/**
 * @brief Implicit function [TM2008, Eq(21)] for ksi depending on observer 
 *   position xi and order of Einstein ring n.
 * @param ksi     Opening angle of Einstein ring.
 * @param params  Observer position xi and order of ring n.
 * @return 
 */
double ksi_impl(double ksi, void *params);


// phi_infinity fuer Startort x_i und Startwinkel ksi
/**
 * @brief Calculate the azimuth angle for a null geodesic starting at xi with an
 *  initial direction ksi.
 * @param xi       Starting position.
 * @param params   Pointer to ksi
 * @return phi_infinity (in radians)
 */
double phiInfty(double xi, void* params);


/**
 * @brief Implicit version of phiInfty.
 * @param ksi
 * @param params
 */
double phiInfty_impl(double ksi, void *params);


#endif // BH_ANALYTIC_H
