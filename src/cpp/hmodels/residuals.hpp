#ifndef RESIDUALS_H
#define RESIDUALS_H

namespace residuals {

  double calculate_penalty(double pars[], double fiber[], double visco[]);

  double calculate_residual_body(double stress[], double weights[],
    int index[], int select[], int dim, int nprot, int skip,
    double sims[]);

  double calculate_weighted_hysteresis_body(double sims[], double deltaCG[], double hysteresis[], double weights[],
    int index[], int select[], int dim, int nprot, int skip);

  double calculate_hyperE_residual(double pars[], double fiber[],
    double visco[], double Tf,
    double args[], double stress[], double dt[], double weights[], int index[], int select[],
    int n, int dim, int nprot, int skip);

  double calculate_hyperE_residual_scaled(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[],
    int index[], int select[], int n, int dim, int nprot, int skip);

  double calculate_viscoE_residual_C_M(double pars[], double fiber[],
    double visco[], double Tf,
    double args[], double stress[], double dt[], double weights[], int index[], int select[],
    int n, int dim, int nprot, int skip);

  double calculate_weighted_viscoE_residual_C_M(double pars[], double fiber[],
    double visco[], double Tf,
    double args[], double stress[], double dt[], double weights[],
    double deltaCG[], double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip);

  double calculate_weighted_viscoE_residual_C_M_scaled(double pars[], double fiber[],
    double visco[], double Tf, double Cmax[],
    double args[], double stress[], double dt[], double weights[], double deltaCG[],
    double hysteresis[], double alphas[],
    int index[], int select[], int n, int dim, int nprot, int skip);

}

#endif