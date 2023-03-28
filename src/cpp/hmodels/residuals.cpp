#define _USE_MATH_DEFINES

#include <cmath>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include "../CTvalues_optimization.hpp"
#include "simulate.hpp"
#include "residuals.hpp"

namespace residuals {


/* ------------------------------------------------------------------------------
 |  This file provides the definitions for calculating the residuals for
 |  the different model forms
 |
 |  These models are mainly used for the collaboration with Alexey Kamenskiy
 |
 |  Author: Will Zhang
 |  Dependencies: None
 ----------------------------------------------------------------------------- */

  using namespace sim;


/* ******************************************************************************
 * Basic functions, e.g. calculating the residual, penalty, hysteresis etc.
 *
 * COMMENTS:
 * This is the same for all of the models, so they share these same codes
****************************************************************************** */
  double calculate_penalty(double pars[], double fiber[], double visco[])
  {

    double d_alpha = fiber[1] - ctv::ideal_alpha;


    return 1.0 + ctv::p_fiber * fiber[0]*fiber[0] + ctv::p_alpha * d_alpha*d_alpha + ctv::p_elastin * pars[4] * pars[4];
  }


  double calculate_residual_body(double stress[], double weights[],
    int index[], int select[], int dim, int nprot, int skip,
    double sims[]
  ){

    int strd_i, strd_k;
    double ds, res, eps;
    int kid, start, stop;

    res = 0;
    for (int k = 0; k < nprot; k++)
    {
      kid    = select[k];
      strd_k = dim*kid;
      start  = index[kid];
      stop   = index[kid + 1];
      for (int j = 0; j < dim; j+=3)
      {
        eps = 0;

        for (int i = start; i < stop; i+=skip)
        {
          strd_i = i*dim + j;
          ds = sims[strd_i] - stress[strd_i];
          eps = eps + ds * ds;

        }

        res = res + weights[strd_k + j] * eps;

      }
    }
    return res;
  }


  double calculate_weighted_hysteresis_body(double sims[], double deltaCG[], double hysteresis[], double weights[],
    int index[], int select[], int dim, int nprot, int skip)
  {

    int strd_i, strd_k;
    double ds, res, hyst;
    int kid, start, stop;

    res = 0;
    for (int k = 0; k < nprot; k++)
    {
      kid    = select[k];
      strd_k = dim*kid;
      start  = index[kid];
      stop   = index[kid + 1];
      for (int j = 0; j < dim; j+=3)
      {
        hyst = 0;

        for (int i = start; i < stop; i+=skip)
        {
          strd_i = i*dim + j;
          hyst   = hyst + sims[strd_i] * deltaCG[strd_i];
        }
        ds = hyst - hysteresis[strd_k + j];
        res = res + weights[strd_k + j] * ds*ds;
      }
    }

    return res;
  }


  double calculate_viscopart_body(double visco[], double alphas[])
  {

    double delta_alpha_C = visco[0] - alphas[0];
    double delta_alpha_CM = 0.5*(visco[0] + visco[1]) - alphas[1];


    return delta_alpha_C*delta_alpha_C + delta_alpha_CM*delta_alpha_CM;
  }

  double lowerbound(double var, double value)
  {

    double rat = var / value;


    return 1.0/(exp(rat) - exp(-rat));
  }

/*******************************************************************************
 * Calculating the residual for the hyperelastic functions
 *
 * COMMENTS:
 * The main form is predetermined, the forms differs by whether the model is
 * scaled.
*******************************************************************************/

  double calculate_hyperE_residual(double pars[], double fiber[], double visco[],
    double Tf,
    double args[], double stress[], double dt[], double weights[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {

    double * sims = new double[n*dim]();

    he_simulate(pars, fiber, visco, Tf, args, dt, &sims[0], n);

    double res = calculate_residual_body(stress, weights, index, select, dim, nprot, skip, sims);

    delete [] sims;

    return res * (1.0 + ctv::p_fiber * fiber[0]*fiber[0]);
  }


  double calculate_hyperE_residual_scaled(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {
    double * sims = new double[n*dim]();

    he_simulate_scaled(pars, fiber, visco, Tf, Cmax, args, dt, &sims[0], n);

    double res = calculate_residual_body(stress, weights, index, select, dim, nprot, skip, sims);

    delete [] sims;

    return res * (1.0 + ctv::p_fiber * fiber[0]*fiber[0]);
  }


/*******************************************************************************
 * Calculating the residual for the viscoelastic models
 *
 * COMMENTS:
 * The main form is predetermined, the forms differs by whether the model is
 * scaled.
*******************************************************************************/

  double calculate_viscoE_residual_C_M(double pars[], double fiber[], double visco[], double Tf,
    double args[], double stress[], double dt[], double weights[], int index[], int select[], int n, int dim,
    int nprot, int skip
  ){

    double * sims = new double[n*dim]();

    caputo_simulate_C_M(pars, fiber, visco, Tf, args, dt, &sims[0], n);

    double res = calculate_residual_body(stress, weights, index, select, dim, nprot, skip, sims);

    delete [] sims;

    return res;
  }


  double calculate_weighted_viscoE_residual_C_M(double pars[], double fiber[],
    double visco[], double Tf,
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {
    double * sims = new double[n*dim]();

    caputo_simulate_C_M(pars, fiber, visco, Tf, args, dt, &sims[0], n);

    double res = calculate_residual_body(stress, weights, index, select, dim, nprot, skip, sims);

    double hyst = calculate_weighted_hysteresis_body(sims, deltaCG, hysteresis, weights, index,
      select, dim, nprot, skip);

    double alp  = calculate_viscopart_body(visco, alphas);

    delete [] sims;

    double scale = 1.0 /std::max(weights[dim*select[nprot - 1]], weights[dim*select[nprot - 1] + 3]);

    return (res + ctv::w_hyst * scale * hyst) * (calculate_penalty(pars, fiber, visco) + ctv::w_visco * alp);
  }


  double calculate_weighted_viscoE_residual_C_M_scaled(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip)
  {
    double * sims = new double[n*dim]();

    caputo_simulate_C_M_scaled(pars, fiber, visco, Tf, Cmax, args, dt, &sims[0], n);

    double res = calculate_residual_body(stress, weights, index, select, dim, nprot, skip, sims);

    double hyst = calculate_weighted_hysteresis_body(sims, deltaCG, hysteresis, weights, index,
      select, dim, nprot, skip);

    double alp  = calculate_viscopart_body(visco, alphas);

    delete [] sims;

    double scale = 1.0 /std::max(weights[dim*select[nprot - 1]], weights[dim*select[nprot - 1] + 3]);

    // std::cout << res << "    " << scale * hyst << "    " << alp <<'\n';

    return (res + ctv::w_hyst * scale * hyst) * (calculate_penalty(pars, fiber, visco) + ctv::w_visco * alp);
  }


/*------------------------------------------------------------------------------
 |  THE END
 -----------------------------------------------------------------------------*/
}

