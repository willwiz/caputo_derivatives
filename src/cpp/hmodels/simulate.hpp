#define _USE_MATH_DEFINES
#ifndef SIM_H
#define SIM_H

#include "../constitutive/neohookean.hpp"
#include "../constitutive/HOG2D.hpp"

namespace sim {

  class setup_models
  {
  public:

    constitutive_models::neohookean matrix;
    constitutive_models::StrucHOG2D collagen;
    constitutive_models::HOG2D      elastin;
    constitutive_models::HOG2D      muscle;

    setup_models() {}
    setup_models(double pars[], double fiber[]);
    setup_models(double pars[], double fiber[], double Cmax[]);
    ~setup_models() {}
    void get_scaled_pars(double pars[]);
  };

  void get_model_parameters(double pars[], double fiber[], double visco[], double Tf,
    double pars_out[11]);

  void get_model_parameters_scaled(double pars[], double fiber[], double visco[], double Tf,
    double Cmax[], double pars_out[11]);

  void he_simulate(double pars[], double fiber[], double caputo[],
    double Tf, double args[], double dt[], double stress[], int n);

  void he_simulate_scaled(double pars[], double fiber[], double caputo[],
    double Tf, double Cmax[], double args[], double dt[], double stress[], int n);

  void caputo_simulate_C_M(double pars[], double fiber[], double caputo[],
    double Tf, double args[], double dt[], double stress[], int n);

  void caputo_simulate_C_M_scaled(double pars[], double fiber[], double caputo[],
    double Tf, double Cmax[], double args[], double dt[], double stress[], int n);

}

#endif