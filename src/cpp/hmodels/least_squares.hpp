#ifndef LSR_H
#define LSR_H

namespace lsqreg {

  void leastsquareM_he(
    double pars[], double fiber[], double caputo[], double Tf,
    double args[], double dt[], double weights[], int index[], int select[],
    int n, int dim, int nprot, int skip,
    double stress[]);

//   void leastsquareM_caputo_CM(
//     double pars[], double fiber[], double caputo[], double Tf,
//     double args[], double dt[], double weights[], int index[], int select[],
//     int n, int dim, int nprot, int skip,
//     double stress[]);

//   double leastsquare_residual(
//     double pars[], double fiber[], double caputo[], double Tf,
//     double args[], double dt[], double weights[], int index[], int select[],
//     int n, int dim, int nprot, int skip,
//     double stress[]);
}

#endif