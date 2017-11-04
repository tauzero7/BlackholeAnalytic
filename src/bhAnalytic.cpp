/**
 * @file    bhAnalytic.cpp
 * @author  Thomas Mueller
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

#include <cassert>
#include "bhAnalytic.h"


int gslBrent(double function(double, void*), void *params,  
    double x_lo, double x_hi, double &x)
{
    int status = 0;
    int iter = 0, max_iter = 200;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;

    gsl_function F;
    F.function = function;
    F.params   = params;

    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);

    status = gsl_root_fsolver_set(s, &F, x_lo, x_hi);
    if (status) return status;

    do {
      iter++;
      status = gsl_root_fsolver_iterate(s);
      if (status) {break;}
      
      x = gsl_root_fsolver_root(s);
      x_lo = gsl_root_fsolver_x_lower(s);
      x_hi = gsl_root_fsolver_x_upper(s);
      status = gsl_root_test_interval(x_lo, x_hi, 0, brentAcc);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    
    return status;
}


double gslIntegral(double function(double, void*), void *params,
    double x_lo, double x_hi)
{
    gsl_integration_workspace* w = gsl_integration_workspace_alloc(1000);
 
    double result, error;
    gsl_function F;
    F.function = function;
    F.params = params;
    
    gsl_integration_qags(&F, x_lo, x_hi, 1e-6, 0.0, 1000, w, &result, &error);   
    gsl_integration_workspace_free(w);
    return result;
}


double calcAqua(const double xi, const double ksi) {
    double sinksi = sin(ksi);
    return xi * xi * (1.0 - xi) / (sinksi * sinksi);
}


int getAzimuthFromEndpoints(const double xi, const double xf, const double ksi,
    double &phi1, double &phi2)
{
    assert(ksi > 0.0 && ksi < PI);
    
    double sk   = sin(ksi);
    double aqua = calcAqua(xi, ksi);

    double x1,x2,x3;
    calcXs(aqua,x1,x2,x3);
    
    double xmin = getXmin(xi, ksi);
    double phimin = getPhiMin(xi, xmin, ksi);
  
    phi1 = phi2 = 0.0;
        
    // geodesic cannot read xf:
    if (xf > xmin) {
        return 0;
    }

    double m = sqrt((x1-x3)/(x2-x3));
    double p1 = 2.0/sqrt(x2-x3)*gsl_sf_ellint_F(asin(1.0/sqrt((xf-x2)/(xf-x1))/m),m,tMode);
    double p2 = 2.0/sqrt(x2-x3)*gsl_sf_ellint_F(asin(1.0/sqrt((xi-x2)/(xi-x1))/m),m,tMode);

    if (ksi > 0.5 * PI) {
        if (xf > xi) {
            return 0;
        }
        else {
            phi1 = phi2 = p1 - p2;
            return 1;
        }
    }
    else {
        if (xf < xi) {
            phi1 = phi2 = p1 + p2;
            return 1;
        }
        else {
            phi1 = p2 - p1;
            phi2 = p1 + p2;
            return 2;
        }
    }
    return 0;
}


double getDistanceFromRing(const double ksi) {
    if (ksi < PI * 0.25) {
        return -two_third * cos((PI + 2.0*ksi) * one_third) + one_third;
    }
    else {
        return two_third * cos((2.0 * PI - 2.0 * ksi) * one_third) + one_third;
    }
}


double getKsiFromDistance(const double xi, const int n) {
    nx_params par = {n, xi};
    double ksiCrit;
    calcKsiCrit(xi, ksiCrit);
    
    double ksi;
    int status = gslBrent(ksi_impl, &par, ksiCrit + 1e-9, PI * 0.5 - 1e-5, ksi);
    if (status) {
        fprintf(stderr, "failed, gsl_errno=%d ... %s:%d\n", status, __FILE__, __LINE__);
    }
    return ksi;
}


void getKsiAndTFromDistance(const double xi, const int n, double &ksi, double &T) {
    nx_params par = {n, xi};
    double ksiCrit;
    calcKsiCrit(xi, ksiCrit);    
    int status = gslBrent(ksi_impl, &par, ksiCrit + 1e-9, PI * 0.5 - 1e-5, ksi);
    if (status) {
        fprintf (stderr, "failed, gsl_errno=%d ... %s:%d\n", status, __FILE__, __LINE__);
    }
    
    // default
    T = 0.0;
    
    dt_params tpar = {xi, ksi};
    if (ksi < 0.5 * PI) {
        double xmin = getXmin(xi, ksi);
        T = 2.0 * gslIntegral(dtdx, &tpar, xi, xmin);
        return;
    }
    
    status = gslBrent(ksi_impl, &par, 0.5*PI + 1e-5, PI - ksiCrit - 1e-9, ksi);
    if (status) {
        fprintf (stderr, "failed, gsl_errno=%d ... %s:%d\n", status, __FILE__, __LINE__);
    }
    
    if (ksi > 0.5 * PI) {
        double xmax = getXmax(xi, ksi);
        T = -2.0 * gslIntegral(dtdx, &tpar, xi, xmax);
    }    
}


double getKsiFromEndpoints(const double xi, const double xf, const double phi) {
    geod_params par = {xi, xf, phi};
    int n = fmod(phi, PI) + 1;
    double ksi_low = getKsiFromDistance(xi, n) + 1e-9;
    double ksi_hi  = calcKsiPhiInfty(xi, phi);
    //std::cerr << ksi_low << " " << ksi_hi << std::endl;    
    //std::cerr << calcXf_impl(ksi_low, &par) << std::endl;
    //std::cerr << calcXf_impl(ksi_hi, &par) << std::endl;

    double ksi = 0.0;
    int status = gslBrent(calcXf_impl, &par, ksi_low, ksi_hi, ksi);
    if (status) {
        fprintf (stderr, "failed, gsl_errno=%d ... %s:%d\n", status, __FILE__, __LINE__);
    }
    return ksi;
}


double getPhiMin(const double xi, const double xmin, const double ksi) {
    double aqua = xmin * xmin * (1.0 - xmin);
    
    double x1,x2,x3;
    calcXs(aqua,x1,x2,x3);

    if (aqua < 4.0 * one_over_ts) {
        double val1 = sqrt((xi - x2) / (xi - x1));
        double val2 = sqrt(fabs((xmin - x2) / (xmin - x1)));
    
        double k = sqrt((x1 - x3) / (x2 - x3));
        double ellipf = 2.0 / sqrt(x2-x3)
            * (gsl_sf_ellint_F(asin(1.0/(val1*k)), k, tMode) 
                - gsl_sf_ellint_F(asin(1.0/val2*k), k, tMode));
        return ellipf;
    }
  
    return 0.0;
}


bool calcKsiCrit(const double xi, double &ksi) {
    if (xi < 0.0 || xi > two_third) {
        return false;
    }
    
    double sk = sqrt(6.75 * xi * xi * (1.0 - xi));
    if (fabs(sk) <= 1.0) {
        ksi = asin(sk);
        return true;
    }
        
    return false;
}


double calcKsiPhiInfty(const double xi, const double phif) {
    double ksi, ksiCrit;
    calcKsiCrit(xi, ksiCrit);
    ksiCrit += 1.0e-8;

    geod_params p;
    p.xi   = xi;
    p.phif = phif;

    int status = gslBrent(&phiInfty_impl, &p, ksiCrit, PI - ksiCrit, ksi);
    if (status) {
        fprintf (stderr, "failed, gsl_errno=%d ... %s:%d\n", status, __FILE__, __LINE__);
    }
    return ksi;
}



double calcXf(const double xi, const double ksi, const double phi) {
    assert(xi > 0.0);
  
    double x1, x2, x3, m2;
    calcXs(calcAqua(xi,ksi), x1, x2, x3);
    m2 = (x1 - x3) / (x2 - x3);
    
    double u = 0.5 * sqrt(x2 - x3) * phi;

    double sn,cn,dn;
    gsl_sf_elljac_e(u, m2, &sn, &cn, &dn);

    double f1 = sn * sqrt((x2-x1)*(x2-x1)*(xi-x3)/((x2-x3)*(xi-x1)*(xi-x1)));
    double f2 = sqrt((xi-x2)/(xi-x1)) * cn * dn;
    double f3 = 1.0-(x1-x3)*(xi-x2)/(x2-x3)/(xi-x1) * sn * sn;

    double SN = (f1+f2)/f3;
    double SN2 = SN * SN;
    return (x2 - x1 * SN2) / (1.0 - SN2);
}


double calcXf_impl(double ksi, void *params) {
    geod_params *p = static_cast<geod_params*>(params);
    double xi = p->xi;
    double xf = p->xf;
    double phi = p->phif;
    return calcXf(xi, ksi, phi) - xf;
}


double calcXiGrenz(const double ksi) {
    double psi;
    if (ksi > 0.0 && ksi <= 0.25 * PI) {
        return -two_third * cos(one_third * (PI + 2.0 * ksi)) + one_third;
    }
    else if (ksi > 0.25 * PI && ksi < PI) {
        return two_third * cos(two_third * (PI - ksi)) + one_third;
    }
    return 0.0;
}


double getXmax(const double xi, const double ksi) {
    double a = xi * sqrt(1.0 - xi) / sin(ksi);
    double w = 0.5 * sqrt(3.0) * a;
    return w / cos(one_third * acos(-3.0 * w));
}


double getXmin(const double xi, const double ksi) {
    if (ksi==0.0) {
        return 0.0;
    }

    double b = 0.5 * sqrt(3.0) * xi * sqrt(1.0 - xi) / sin(ksi);
    double val = -3.0*b;
  
    if (fabs(val) > 1.0) {
        return 0.0;
    }

    double psidr = acos(val) * one_third;
    return b / cos(psidr);
}


void calcXs(const double aqua, double &x1, double &x2, double &x3) {
    double q   = aqua * 0.5 - one_over_ts;
    double rho = sign(q) * one_third;
  
    double psi  = acos(q / pow(rho,3.0));
    double psidr = psi * one_third;

    x1 =  2.0 * rho * cos(PI * one_third + psidr) + one_third;
    x2 =  2.0 * rho * cos(PI * one_third - psidr) + one_third;
    x3 = -2.0 * rho * cos(psidr) + one_third;
  
    if (sign(x1 * x2) < 0.0) {
        //x1 =  2.0 * rho * cos(PI * one_third + psidr) + one_third;
        //x3 =  2.0 * rho * cos(PI * one_third - psidr) + one_third;
        //x2 = -2.0 * rho * cos(psidr) + one_third;
        std::swap(x2, x3);
    }
}


double dtdx(double x, void *params) {
    dt_params* par = static_cast<dt_params*>(params);
    double xi = par->xi;
    double ksi = par->ksi;
    double aqua = calcAqua(xi, ksi);

    return 1.0 / (x * x * (1.0 - x) * sqrt(1.0 - (1.0 - x) * x * x / aqua));
}



double ksi_impl(double ksi, void *params) {
    nx_params* par = static_cast<nx_params*>(params);
    int n = par->n;
    double xi = par->xi;
    
    double a2,x1,x2,x3;
    
    a2 = calcAqua(xi, ksi);
    if (a2 >= 4.0 * one_over_ts) {
        return 0.0;
    }
    
    calcXs(a2, x1, x2, x3);
    double m = sqrt((x1 - x3) / (x2 - x3));
    double x = sqrt((xi - x2) / (xi - x1));
    
    double ellf;
    if (m <= 1.0) {
        if (m * x < 1.0) {
            ellf = 2.0/sqrt(x2-x3)*gsl_sf_ellint_F(asin(m*x), m, tMode);
        }
        else {
            ellf = 2.0/sqrt(x2-x3)*gsl_sf_ellint_F(asin(1.0/(m*x)),m, tMode);
        }
    }
    else if (m > 1.0 && m * x < 1.0) {
        ellf = 2.0/sqrt(x2-x3)*gsl_sf_ellint_F(asin(m*x),1.0/m, tMode) / m;
    }
    else {
        // uuups
    }

// TODO ???
    if (ksi < 0.5 * PI) {
        return ellf - n * PI;
    }
    else {
        return ellf - n * PI;
    }
}


double phiInfty(double xi, void* params) {
    double ksi = *(static_cast<double*>(params));
      
    double param_q = 0.5 * calcAqua(xi,ksi) - one_over_ts;
    double param_rho = sign(param_q) * one_third;
    double psi       = acos(param_q / pow(param_rho, 3.0));
    //std::cerr << ksi << " " << param_q << " " << param_rho << " " << param_q / pow(param_rho, 3.0) << " " << psi << std::endl;
    
    double x1 = 2.0 * param_rho * cos((PI + psi) * one_third) + one_third;
    double x2, x3;
    if (param_rho > 0.0) {
        x2 =  2.0 * param_rho * cos((PI - psi) * one_third) + one_third;
        x3 = -2.0 * param_rho * cos(psi * one_third) + one_third;
    }
    else {
        x3 =  2.0 * param_rho * cos((PI - psi) * one_third) + one_third;
        x2 = -2.0 * param_rho * cos(psi * one_third) + one_third;
    }
    //std::cerr << x1 << " " << x2 << " " << x3 << std::endl;
    
    double m = sqrt((x1 - x3) / (x2 - x3));

    double F1 = gsl_sf_ellint_F(asin(1.0 / (m * sqrt(x2 / x1))), m, 8);
    double F2 = gsl_sf_ellint_F(asin(1.0 / (m * sqrt((xi - x2) / (xi - x1)))), m, 8);
  
    double pm = 1.0;
    if (ksi > 0.5 * PI) pm = -1.0;

    return 2.0 / sqrt(x2 - x3) * (F1 + pm * F2);
}


double phiInfty_impl(double ksi, void *params) {
    geod_params* p = static_cast<geod_params*>(params);
    double xi = p->xi;
    double phif = p->phif;

    return phiInfty(xi, &ksi) - phif;
}

