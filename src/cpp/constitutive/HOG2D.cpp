#include <algorithm>
#include <cmath>
#include "../kinematics/tensor_algebra.hpp"
#include "../kinematics/kinematics.hpp"
#include "HOG2D.hpp"



/*----------------------------------------------------------------------
 |  This file provides the holzapfel class of models in 2D. Single and
 |  dual family models are available
 |
 |  Author: Will Zhang
 |  Dependencies: None
 -----------------------------------------------------------------------*/


namespace constitutive_models {

/*******************************************************************************
 * Structural tensor forms
 *
 * COMMENTS:
 *
*******************************************************************************/

  void StrucHOG2D::set_pars(double k1, double k2, double theta, double alpha, double beta,
    double kip, double kop)
  {

    this -> k1 = k1;
    this -> k2 = k2;
    this -> A = 2.0*kop*kip;
    this -> B = 2.0*kop*(1.0-2.0*kip);
    this -> C = 1.0 - 3.0*A - B;

    double ca4 = cos(theta + alpha);
    double sa4 = sin(theta + alpha);
    double ca6 = cos(theta + beta);
    double sa6 = sin(theta + beta);

    this -> m4[0] = ca4*ca4;
    this -> m4[1] = ca4*sa4;
    this -> m4[2] = m4[1];
    this -> m4[3] = sa4*sa4;

    this -> m6[0] = ca6*ca6;
    this -> m6[1] = ca6*sa6;
    this -> m6[2] = m6[1];
    this -> m6[3] = sa6*sa6;

    for (int i = 0; i < 4; i++)
    {
      this -> H4[i] = A*id2d[i] + B*m4[i];
      this -> H6[i] = A*id2d[i] + B*m6[i];
    }
  }

  // Scaled version
  void StrucHOG2D::set_pars(double k1, double k2, double theta, double alpha, double beta,
    double kip, double kop, double Cmax[])
  {
    this -> set_pars(k1, k2, theta, alpha, beta, kip, kop);

    double det = Cmax[0]*Cmax[3] - Cmax[1]*Cmax[1];
    double I_n = 1.0 / det;
    double I_1 = Cmax[0] + Cmax[3] + I_n;
    double I_4 = ddot(m4, Cmax);
    double I_6 = ddot(m6, Cmax);
    double E6  = A * I_1 + C * I_n - 1.0;
    double E4  = E6 + B*I_4;
    E6  = E6 + B*I_6;
    E1  = 0.5 * (E4 + E6);
    this->k1 = k1 / E1;
    this->E2 = E4 * E4;
  }


  double StrucHOG2D::get_scaled_modulus()
  {
    return k1 * exp(-k2*E2);
  }

  // Stress functions
  double StrucHOG2D::stress(const kinematics::deformation2D &kin, double stress[4]){

    double I_4 = ddot(m4, kin.C);
    double I_6 = ddot(m6, kin.C);
    double E6  = A * kin.I_1 + C * kin.I_n - 1.0;
    double E4  = E6 + B*I_4;
    E6  = E6 + B*I_6;
    double dWd4 = k1*E4*exp(k2*(E4*E4 - E2));
    double dWd6 = k1*E6*exp(k2*(E6*E6 - E2));

    for (int i = 0; i < 4; i++)
    {
      stress[i] = dWd4*H4[i] + dWd6*H6[i];
    }
    return C * (dWd4 + dWd6);
  }

  void StrucHOG2D::stress(double args[4], double stress[4]){

    kinematics::deformation2D kin(args);
    double p = this->stress(kin, stress);
    for (int i = 0; i < 4; i++)
    {
      stress[i] = stress[i] - p*kin.C33Cinv[i];
    }
  }


/*******************************************************************************
 * Standard holzapfel ogden with two fiber families
 *
 * COMMENTS:
 *
*******************************************************************************/

  void HOGDouble2D::set_pars(double k1, double k2, double theta, double alpha) {

    this -> k1 = k1;
    this -> k2 = k2;

    double ca4 = cos(theta + alpha);
    double sa4 = sin(theta + alpha);
    double ca6 = cos(theta - alpha);
    double sa6 = sin(theta - alpha);

    this -> m4[0] = ca4*ca4;
    this -> m4[1] = ca4*sa4;
    this -> m4[2] = m4[1];
    this -> m4[3] = sa4*sa4;

    this -> m6[0] = ca6*ca6;
    this -> m6[1] = ca6*sa6;
    this -> m6[2] = m6[1];
    this -> m6[3] = sa6*sa6;
  }

  // Stress functions
  double HOGDouble2D::stress(const kinematics::deformation2D &kin, double stress[4]){


    double I_4 = ddot(m4, kin.C) - 1;
    double I_6 = ddot(m6, kin.C) - 1;
    double dWd4 = k1*I_4*exp(k2*I_4 * I_4);
    double dWd6 = k1*I_6*exp(k2*I_6 * I_6);
    for (int i = 0; i < 4; i++)
    {
      stress[i] = dWd4*m4[i] + dWd6*m6[i];
    }
    return 0.0;
  }

  void HOGDouble2D::stress(double args[4], double stress[4]){

    kinematics::deformation2D kin(args);
    (void) this->stress(kin, stress);
    // for (int i = 0; i < 4; i++)
    // {
    //   stress[i] = stress[i] - p*kin.C33Cinv[i];
    // }
  }


/*******************************************************************************
 * Standard holzapfel ogden with one fiber family
 *
 * COMMENTS:
 *
*******************************************************************************/

  void HOG2D::set_pars(double k1, double k2, double theta) {
    this -> k1 = k1;
    this -> k2 = k2;
    double c = cos(theta);
    double s = sin(theta);
    this -> m[0] = c*c;
    this -> m[1] = c*s;
    this -> m[2] = this -> m[1];
    this -> m[3] = s*s;
  }

  // Scaled version
  void HOG2D::set_pars(double k1, double k2, double theta, double Cmax[]) {
    this-> set_pars(k1, k2, theta);
    E1 = ddot(m, Cmax) - 1;
    this-> k1 = k1 / E1;
    this-> E2 = E1*E1;
  }

  double HOG2D::get_scaled_modulus()
  {
    return k1 * exp(-k2*E2);
  }

  // Stress functions
  double HOG2D::stress(const kinematics::deformation2D &kin, double stress[4]){

    double I_4 = ddot(m, kin.C) - 1;
    double dWd4 = k1*I_4*exp(k2* (I_4 * I_4 - E2));

    for (int i = 0; i < 4; i++)
    {
      stress[i] = dWd4*m[i];
    }
    return 0.0;
  }

  void HOG2D::stress(double args[4], double stress[4]){

    kinematics::deformation2D kin(args);
    (void) this->stress(kin, stress);
  }

/*------------------------------------------------------------------------------
 |  THE END
 -----------------------------------------------------------------------------*/

}