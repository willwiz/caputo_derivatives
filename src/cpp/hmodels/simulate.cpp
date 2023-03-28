#define _USE_MATH_DEFINES

#include <cmath>
#include <stdio.h>
#include <iostream>
#include "../CTvalues_optimization.hpp"
#include "../kinematics/tensor_algebra.hpp"
#include "../kinematics/kinematics.hpp"
#include "../constitutive/neohookean.hpp"
#include "../constitutive/HOG2D.hpp"
#include "../constitutive/caputo.hpp"
#include "simulate.hpp"

/*----------------------------------------------------------------------
 |  This file provides the definitions of the different model forms
 |  which combines the primative constitive model from the other cpp files
 |  in the constitutive_models namespace
 |
 |  These models are mainly used for the collaboration with Alexey Kamenskiy
 |
 |  Author: Will Zhang
 |  Dependencies: None
 -----------------------------------------------------------------------*/


using namespace constitutive_models;

namespace sim {


/*----------------------------------------------------------------------
 |  This provides the main models in the full constitutive model
 -----------------------------------------------------------------------*/

  // The Non-scaled version
  setup_models::setup_models(double pars[], double fiber[]):
      matrix(pars[0]),
      collagen(pars[1], pars[2], fiber[0], fiber[1], -fiber[1], ctv::M_kip, ctv::M_kop),
      elastin(pars[3], pars[4], fiber[0]),
      muscle(pars[5], pars[6], fiber[0] + M_PI_2)
  {}

  // Scaled for the maximum right Cauchy Green tensor
  setup_models::setup_models(double pars[], double fiber[], double Cmax[]):
      matrix(pars[0]),
      collagen(pars[1], pars[2], fiber[0], fiber[1], -fiber[1], ctv::M_kip, ctv::M_kop, Cmax),
      elastin(pars[3], pars[4], fiber[0], Cmax),
      muscle(pars[5], pars[6], fiber[0] + M_PI_2, Cmax)
  {}

  void setup_models::get_scaled_pars(double pars[])
  {
    pars[0] = matrix.mu;
    pars[1] = collagen.get_scaled_modulus();
    pars[2] = collagen.k2;
    pars[3] = elastin.get_scaled_modulus();
    pars[4] = elastin.k2;
    pars[5] = muscle.get_scaled_modulus();
    pars[6] = muscle.k2;
  }

  void get_model_parameters(double pars[], double fiber[], double visco[], double Tf,
    double pars_out[11])
  {
    setup_models psi(pars, fiber);
    psi.get_scaled_pars(&pars_out[0]);
    pars_out[7]  = fiber[0];
    pars_out[8]  = fiber[1];
    pars_out[9]  = visco[0];
    pars_out[10] = visco[1];
  }

  void get_model_parameters_scaled(double pars[], double fiber[], double visco[], double Tf,
    double Cmax[], double pars_out[11])
  {
    setup_models psi(pars, fiber, Cmax);
    psi.get_scaled_pars(&pars_out[0]);
    pars_out[7]  = fiber[0];
    pars_out[8]  = fiber[1];
    pars_out[9]  = visco[0];
    pars_out[10] = visco[1];
  }
/*----------------------------------------------------------------------
 |  The hyperelastic models
 -----------------------------------------------------------------------*/

  void he_simulate(double pars[], double fiber[], double caputo[], double Tf,
    double args[], double dt[], double stress[], int n
  ) {

    int strd_i;

    double p_g, p_c, p;
    double matrix[ctv::prob_dim], collagen[ctv::prob_dim], elastin[ctv::prob_dim], muscle[ctv::prob_dim];

    kinematics::deformation2D kin;
    setup_models psi(pars, fiber);

    for (int i = 1; i < n; i++)
    {
      strd_i = ctv::prob_dim * i;
      // Calculate deformation
      kin.precompute(&args[strd_i]);
      // Calculate hyper E
      p_g  = psi.matrix.stress(kin,   matrix);
      p_c  = psi.collagen.stress(kin, collagen);
      (void) psi.elastin.stress(kin,  elastin);
      (void) psi.muscle.stress(kin,   muscle);
      // Fractional Updates
      // Calculate pressure
      p = p_g + p_c;
      // Compute Final Stress
      for (int j = 0; j < ctv::prob_dim; j++)
      {
        stress[strd_i + j] = matrix[j] + collagen[j] + elastin[j] + muscle[j] - p * kin.C33Cinv[j];
      }
    }
  }


  void he_simulate_scaled(double pars[], double fiber[], double caputo[], double Tf, double Cmax[],
    double args[], double dt[], double stress[], int n)
  {

    int strd_i;

    double p_g, p_c, p;
    double matrix[ctv::prob_dim], collagen[ctv::prob_dim], elastin[ctv::prob_dim], muscle[ctv::prob_dim];

    kinematics::deformation2D kin;
    setup_models psi(pars, fiber, Cmax);

    for (int i = 1; i < n; i++)
    {
      strd_i = ctv::prob_dim * i;
      // Calculate deformation
      kin.precompute(&args[strd_i]);
      // Calculate hyper E
      p_g  = psi.matrix.stress(kin,   matrix);
      p_c  = psi.collagen.stress(kin, collagen);
      (void) psi.elastin.stress(kin,  elastin);
      (void) psi.muscle.stress(kin,   muscle);
      // Fractional Updates
      // Calculate pressure
      p = p_g + p_c;
      // Compute Final Stress
      for (int j = 0; j < ctv::prob_dim; j++)
      {
        stress[strd_i + j] = matrix[j] + collagen[j] + elastin[j] + muscle[j] - p * kin.C33Cinv[j];
      }
    }
  }


/*----------------------------------------------------------------------
 |  The viscoelastic models, with the collagen and S.M.C. being viscoelastic
 -----------------------------------------------------------------------*/

  void caputo_simulate_C_M(double pars[], double fiber[], double caputo[], double Tf,
    double args[], double dt[], double stress[], int n
  ) {

    int strd_i;

    double p_g, p_c, p;
    double matrix[ctv::prob_dim], collagen[ctv::prob_dim], elastin[ctv::prob_dim], muscle[ctv::prob_dim];

    kinematics::deformation2D kin;
    setup_models psi(pars, fiber);

    double f_c;
    double frac_c[ctv::prob_dim], frac_m[ctv::prob_dim];
    caputo::caputo_init_vec<4> carpC(caputo[0], Tf, 0.0);
    caputo::caputo_init_scl    carpC0(caputo[0], Tf, 0.0);
    caputo::caputo_init_vec<4> carpM(caputo[1], Tf, 0.0);

    for (int i = 1; i < n; i++)
    {
      strd_i = ctv::prob_dim * i;
      // Calculate deformation
      kin.precompute(&args[strd_i]);
      // Calculate hyper E
      p_g  = psi.matrix.stress(kin,   matrix);
      p_c  = psi.collagen.stress(kin, collagen);
      (void) psi.elastin.stress(kin,  elastin);
      (void) psi.muscle.stress(kin,   muscle);
      // Fractional Updates
      carpC.caputo_iter(collagen, dt[i], frac_c);
      carpM.caputo_iter(muscle,   dt[i], frac_m);
      f_c = carpC0.caputo_iter(p_c, dt[i]);
      // Calculate pressure
      p = p_g + f_c;
      // Compute Final Stress
      for (int j = 0; j < ctv::prob_dim; j++)
      {
        stress[strd_i + j] = matrix[j] + frac_c[j] + elastin[j] + frac_m[j] - p * kin.C33Cinv[j];
      }
    }
  }


  void caputo_simulate_C_M_scaled(double pars[], double fiber[], double caputo[], double Tf, double Cmax[],
    double args[], double dt[], double stress[], int n)
  {

    int strd_i;

    double p_g, p_c, p;
    double matrix[ctv::prob_dim], collagen[ctv::prob_dim], elastin[ctv::prob_dim], muscle[ctv::prob_dim];

    kinematics::deformation2D kin;
    setup_models psi(pars, fiber, Cmax);

    double f_c;
    double frac_c[ctv::prob_dim], frac_m[ctv::prob_dim];
    caputo::caputo_init_vec<4> carpC(caputo[0], Tf, 0.0);
    caputo::caputo_init_scl    carpC0(caputo[0], Tf, 0.0);
    caputo::caputo_init_vec<4> carpM(caputo[1], Tf, 0.0);

    for (int i = 1; i < n; i++)
    {
      strd_i = ctv::prob_dim * i;
      // Calculate deformation
      kin.precompute(&args[strd_i]);
      // Calculate hyper E
      p_g  = psi.matrix.stress(kin,   matrix);
      p_c  = psi.collagen.stress(kin, collagen);
      (void) psi.elastin.stress(kin,  elastin);
      (void) psi.muscle.stress(kin,   muscle);
      // Fractional Updates
      carpC.caputo_iter(collagen, dt[i], frac_c);
      carpM.caputo_iter(muscle,   dt[i], frac_m);
      f_c = carpC0.caputo_iter(p_c, dt[i]);
      // Calculate pressure
      p = p_g + f_c;
      // Compute Final Stress
      for (int j = 0; j < ctv::prob_dim; j++)
      {
        stress[strd_i + j] = matrix[j] + frac_c[j] + elastin[j] + frac_m[j] - p * kin.C33Cinv[j];
      }
    }
  }
}