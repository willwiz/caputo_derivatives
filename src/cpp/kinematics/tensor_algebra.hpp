#ifndef TENSOR_ALGEBRA_H
#define TENSOR_ALGEBRA_H

namespace constitutive_models {

  extern const double id2d[4];

  double ddot(const double a[4], const double b[4]);

  void addto(const double a[], double b[], int dim);

}

#endif