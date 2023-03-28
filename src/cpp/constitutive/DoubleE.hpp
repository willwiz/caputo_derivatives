#ifndef DoubleE2D_H
#define DoubleE2D_H

#include "kinematics.hpp"

namespace constitutive_models {
  class DoubleE2D {
    public:
      double b1, b2;
      double k1, k2, k3;
      double mxm[4], nxn[4], mxn[4];

      DoubleE2D () {};
      DoubleE2D (double b1, double b2, double k1, double k2, double k3, double theta) {
        this->set_pars(b1, b2, k1, k2, k3, theta);
      };
      ~DoubleE2D () {};
      void set_pars(double b1, double b2, double k1, double k2, double k3, double theta);

      double stress(const kinematics::deformation2D &kin, double stress[]);
      void stress(double args[], double stress[]);
  };

}

#endif