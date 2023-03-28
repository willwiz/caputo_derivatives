#define _USE_MATH_DEFINES

#include <cmath>
#include <stdio.h>
#include <iostream>
#include "../kinematics/tensor_algebra.hpp"
#include "../kinematics/kinematics.hpp"
#include "../constitutive/neohookean.hpp"
#include "../constitutive/HOG2D.hpp"
#include "../constitutive/caputo.hpp"
#include "least_squares.hpp"

using namespace constitutive_models;

namespace lsqreg {

/*----------------------------------------------------------------------
 |  Do part of the optimization using least squares regression, useful
 |  idea but doens't work due to allowing for negative parameters.
 |
 |  Author: Will Zhang
 |  Dependencies: None
 -----------------------------------------------------------------------*/

  void leastsquareM_he(
    double pars[], double fiber[], double caputo[], double Tf,
    double args[], double dt[], double weights[], int index[], int select[],
    int n, int dim, int nprot, int skip,
    double stress[]
  ) {

    const int nlaws = 4;

    int strd_i, strd_j, strd_m;

    double p;
    double * sims = new double[n*nlaws*dim]();

    kinematics::deformation2D kin;
    neohookean matrix(1.0);
    StrucHOG2D collagen(1.0, pars[0], fiber[0], fiber[1], -fiber[1], fiber[2], fiber[3]);
    HOG2D elastin(1.0, pars[1], fiber[0]);
    HOG2D muscle(1.0, pars[2], fiber[0] + M_PI_2);

    // Set initial state to 0
    for (int j = 0; j < nlaws; j++)
    {
      strd_j = dim * j;
      for (int k = 0; k < dim; k++)
      {
        stress[strd_j + k] = 0;
      }
    }


    for (int i = 0; i < n; i++)
    {
      strd_i = i * dim;

      kin.precompute(&args[strd_i]);

      strd_j = strd_i * nlaws;
      p = matrix.stress(kin, &sims[strd_j]);
      for (int k = 0; k < dim; k++)
      {
        sims[strd_j + k] = sims[strd_j + k] - p * kin.C33Cinv[k];
      }

      strd_j = strd_j + dim;
      p = collagen.stress(kin, &sims[strd_j]);
      for (int k = 0; k < dim; k++)
      {
        sims[strd_j + k] = sims[strd_j + k] - p * kin.C33Cinv[k];
      }

      strd_j = strd_j + dim;
      (void) elastin.stress(kin, &sims[strd_j]);
      strd_j = strd_j + dim;
      (void) muscle.stress(kin, &sims[12]);
    }

    for (int m = 0; m < nprot; m++)
    {
      strd_m = m*dim;
      for (int i = index[select[m]]; i < index[select[m] + 1]; i+=skip)
      {
        strd_i = i*dim*nlaws;
        for (int j = 0; j < nlaws; j++)
        {
          strd_j =j*dim;
          for (int k = 0; k < dim; k++)
          {
            stress[strd_i + k*nlaws + j] = sqrt(weights[strd_m + k]) * sims[strd_i + strd_j + k];
          }
        }
      }
    }

    delete [] sims;
  }

}