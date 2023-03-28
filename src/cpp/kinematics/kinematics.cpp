#define _USE_MATH_DEFINES

#include <cmath>
#include <stdio.h>
#include <iostream>
#include "kinematics.hpp"

namespace kinematics {

/*----------------------------------------------------------------------
 |  This file provides the structure for storing the info for deformation
 |  gradient. It is helpful for reducing recomputation.
 |
 |  Author: Will Zhang
 |  Dependencies: None
 -----------------------------------------------------------------------*/

  deformation2D::deformation2D () : C{}, Cinv{}, C33Cinv{} {}

  deformation2D::deformation2D (double args[4]) : C{}, Cinv{}, C33Cinv{}
  {
    this -> precompute(args);
  }

  deformation2D::~deformation2D () {}

  void deformation2D::precompute(double args[4]) {

    this -> det = args[0]*args[3] - args[1]*args[1];
    this -> I_n = 1 / det;
    this -> Cinv[0] =  args[3] * I_n;
    this -> Cinv[1] = -args[1] * I_n;
    this -> Cinv[2] = -args[2] * I_n;
    this -> Cinv[3] =  args[0] * I_n;

    for (int i = 0; i < 4; i++)
    {
      this ->C[i] = args[i];
      this ->C33Cinv[i] = I_n*Cinv[i];
    }

    this -> I_1 = args[0] + args[3] + I_n;
    this -> I_1m3 = I_1 - 3.0;
  }


}