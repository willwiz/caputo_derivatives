#ifndef HOG2D_H
#define HOG2D_H

#include "../kinematics/kinematics.hpp"

namespace constitutive_models {

  class StrucHOG2D
  {
  public:
    double k1, k2;
    double A, B, C;
    double m4[4], m6[4], H4[4], H6[4];
    double E1;
    double E2;

    StrucHOG2D (): E1{}, E2{} {};
    StrucHOG2D (double k1, double k2, double theta, double alpha, double beta, double kip,
      double kop) : E1{}, E2{}
    {
      this->set_pars(k1, k2, theta, alpha, beta, kip, kop);
    };
    StrucHOG2D (double k1, double k2, double theta, double alpha, double beta, double kip, double kop,
    double Cmax[])
    {
      this->set_pars(k1, k2, theta, alpha, beta, kip, kop, Cmax);
    };
    ~StrucHOG2D () {};
    void set_pars(double k1, double k2, double theta, double alpha, double beta, double kip,
      double kop);
    void set_pars(double k1, double k2, double theta, double alpha, double beta, double kip,
      double kop, double Cmax[]);
    double get_scaled_modulus();
    double stress(const kinematics::deformation2D &kin, double stress[4]);
    void stress(double args[4], double stress[4]);
  };

  class HOGDouble2D
  {
  public:
    double k1, k2;
    double m4[4], m6[4];

    HOGDouble2D () {};
    HOGDouble2D (double k1, double k2, double theta, double alpha) {
      this->set_pars(k1, k2, theta, alpha);
    };
    ~HOGDouble2D () {};
    void set_pars(double k1, double k2, double theta, double alpha);
    double stress(const kinematics::deformation2D &kin, double stress[4]);
    void stress(double args[4], double stress[4]);
  };



  class HOG2D
  {
  public:
    double k1, k2;
    double m[4];
    double E1;
    double E2;

    HOG2D () : E1{}, E2{} {};
    HOG2D (double k1, double k2, double theta) : E1{}, E2{}
    {
      this->set_pars(k1, k2, theta);
    };
    HOG2D (double k1, double k2, double theta, double Cmax[]) {
      this->set_pars(k1, k2, theta, Cmax);
    };
    ~HOG2D () {};
    void set_pars(double k1, double k2, double theta);
    void set_pars(double k1, double k2, double theta, double Cmax[]);
    double get_scaled_modulus();
    double stress(const kinematics::deformation2D &kin, double stress[4]);
    void stress(double args[4], double stress[4]);
  };

}

#endif