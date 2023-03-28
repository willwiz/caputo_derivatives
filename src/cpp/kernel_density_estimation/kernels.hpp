#ifndef KDE_KERNELS_H
#define KDE_KERNELS_H

namespace kernels {

  class Triweight_kernel
  {
  private:
    double m_mean;
    double m_1_sigma;
    double m_scaling;

  public:
    Triweight_kernel(/* args */);
    Triweight_kernel(double mean, double bandwidth, int n);
    ~Triweight_kernel();

    void set_mean(double mean, double bandwidth, int n);

    double pdf(double x);
  };

  class Gaussian_kernel
  {
  private:
    double m_mean;
    double m_1_sigma;
    double m_scaling;

  public:
    Gaussian_kernel(/* args */);
    Gaussian_kernel(double mean, double bandwidth, int n);
    ~Gaussian_kernel();

    void set_mean(double mean, double bandwidth, int n);

    double pdf(double x);
  };

  class Gamma_kernel
  {
  private:
    double m_a;
    double m_b;
    double m_mean;
    double m_1_sigma;
    double m_scaling;

  public:
    Gamma_kernel(/* args */);
    Gamma_kernel(double mean, double bandwidth, int n);
    ~Gamma_kernel();

    void set_mean(double mean, double bandwidth, int n);

    double pdf(double x);
  };

}

#endif