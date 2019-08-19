/**
 * @file    defs.h
 * @author  Thomas Mueller
 * 
 * @brief   ...
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
#ifndef BH_DEFS_H
#define BH_DEFS_H

#include <gsl/gsl_sf_ellint.h>

const double brentAcc = 1.0e-10;
const gsl_mode_t tMode = 8;
const int maxIter = 250;

constexpr double one_third = 1.0 / 3.0;
constexpr double two_third = 2.0 / 3.0;
constexpr double one_over_ts = 1.0 / 27.0;
constexpr double four_over_ts = 4.0 / 27.0;

constexpr double PI = 3.141592653589793238463;
constexpr double RAD_TO_DEG = 180.0 / PI;
constexpr double DEG_TO_RAD = PI / 180.0;

#define sign(x) ((x >= 0.0) ? (1.0) : (-1.0))


typedef struct {
    double xi;
    double ksi;
} dt_params;

typedef struct {
    int n;
    double ksi;
} xf_params;


typedef struct {
    int n;
    double xi;
} nx_params;


typedef struct {
    double ksi0;
    double ksi1;
} xfs_params;


typedef struct {
    double xi;
    double xf;
    double phif;
} geod_params;

#endif // BH_DEFS_H
