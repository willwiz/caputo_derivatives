#include <cmath>
#include "../kinematics/tensor_algebra.hpp"
#include "../kinematics/kinematics.hpp"
#include "neohookean.hpp"

namespace constitutive_models {

/*----------------------------------------------------------------------
 |  Standard Neohookean, no much to say
 |
 |  Author: Will Zhang
 |  Dependencies: None
 -----------------------------------------------------------------------*/


  void neohookean::set_pars(double mu) {

    this -> mu = mu;
  }

  double neohookean::stress(const kinematics::deformation2D &kin, double stress[]){

    for (int i = 0; i < 4; i++)
    {
      stress[i] = mu*id2d[i];
    }
    return mu;
  }

  void neohookean::stress(double args[], double stress[]){

    kinematics::deformation2D kin(args);
    double p = this->stress(kin, stress);
    for (int i = 0; i < 4; i++)
    {
      stress[i] = stress[i] - p*kin.C33Cinv[i];
    }
  }

  // Fung form
  void fungIso::set_pars(double mu, double k) {
    this -> mu = mu;
    this -> k  = k;
  }

  double fungIso::stress(const kinematics::deformation2D &kin, double stress[]){

    exponent = mu*exp(k*kin.I_1m3);

    for (int i = 0; i < 4; i++)
    {
      stress[i] = exponent*id2d[i];
    }
    return exponent;
  }

  void fungIso::stress(double args[], double stress[]){

    kinematics::deformation2D kin(args);
    double p = this->stress(kin, stress);
    for (int i = 0; i < 4; i++)
    {
      stress[i] = stress[i] - p*kin.C33Cinv[i];
    }
  }
}
