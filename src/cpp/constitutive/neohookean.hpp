#ifndef NEOHOOKEAN_H
#define NEOHOOKEAN_H

#include "../kinematics/kinematics.hpp"

namespace constitutive_models {

  class neohookean {
    public:
      double mu;

      neohookean () {};
      neohookean (double mu) {
        this->set_pars(mu);
      };
      ~neohookean () {};
      void set_pars(double mu);
      double stress(const kinematics::deformation2D &kin, double stress[]);
      void stress(double args[], double stress[]);
  };

  class fungIso {
    public:
      double mu;
      double k;

      double exponent;

      fungIso () {};
      fungIso(double mu, double k){
        this->set_pars(mu, k);
      };
      ~fungIso () {};
      void set_pars(double mu, double k);
      double stress(const kinematics::deformation2D &kin, double stress[]);
      void stress(double args[], double stress[]);
  };
}

#endif