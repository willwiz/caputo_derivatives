#ifndef CT_PARAMETERS_H
#define CT_PARAMETERS_H

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

namespace ctv
{
  extern const double p_fiber = 0.001;
  extern const double p_alpha = 0.1;
  extern const double p_elastin = 1e-3;
  extern const double b_visco = 0.01;
  extern const double b_modulus = 0.01;

  extern const double w_hyst = 0.001;
  extern const double w_visco = 100.0;

  extern const int    prob_dim = 4;

  extern const double M_kip = 0.155;
  extern const double M_kop = 0.424285714285714;

  extern const double ideal_alpha = M_PI_4;

} // namespace ct_pars

#endif