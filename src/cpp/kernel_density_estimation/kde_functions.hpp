#ifndef KDE_FUNC_H
#define KDE_FUNC_H

namespace kde_functions {

  void kde_gaussian_bounded_estimate(const double mean[], const int n_mean, double bandwidth,
    const double x[], const int n_x, double y[]);

  void kde_gaussian_estimate(const double mean[], const int n_mean, double bandwidth,
    const double x[], const int n_x, double y[]);
}

#endif