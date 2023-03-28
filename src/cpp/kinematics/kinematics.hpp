#ifndef KINEMATICS_H
#define KINEMATICS_H

namespace kinematics {
  class deformation2D {
    public:
      double det;
      double I_n;
      double I_1;
      double I_1m3;
      double C[4];
      double Cinv[4];
      double C33Cinv[4];

      deformation2D();
      deformation2D(double vC[4]);
      ~deformation2D();
      void precompute(double vC[4]);
  };

}

#endif