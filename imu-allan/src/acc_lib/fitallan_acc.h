#ifndef FitAllanAcc_H
#define FitAllanAcc_H

#include <ceres/ceres.h>
#include <cmath>
#include <eigen3/Eigen/Eigen>

namespace imu
{
// calculate coefficients to fit Allan curve line
class FitAllanAcc
{
  class AllanSigmaError
  {
   public:
    AllanSigmaError(const double &_sigma2, const double &_tau) :
      sigma2(_sigma2), tau(_tau)
    {
    }

    template <typename T>
    T calcLog10(T src) const
    {
      return (log(src)) / (log(10));
    }

    template <typename T>
    T calcSigma2(T _Q, T _N, T _B, T _K, T _R, T _tau) const
    {
      // clang-format off
      return  _Q * _Q / ( _tau * _tau )
            + _N * _N / _tau
            + _B * _B
            + _K * _K * _tau
            + _R * _R * _tau * _tau;
      // clang-format on
    }

    template <typename T>
    bool operator()(const T *const _param, T *residuals) const
    {
      T _Q   = T(_param[0]);
      T _N   = T(_param[1]);
      T _B   = T(_param[2]);
      T _K   = T(_param[3]);
      T _R   = T(_param[4]);
      T _tau = T(tau);

      T _sigma2    = calcSigma2(_Q, _N, _B, _K, _R, _tau);
      T _d_sigma2  = T(calcLog10(_sigma2)) - T(calcLog10(sigma2));
      residuals[0] = _d_sigma2;
      // std::cout << "_err " << T(sigma2) << " " << _sigma2 << std::endl;

      return true;
    }

    double sigma2; // unit is m^2/s^4=(m/s^2)^2
    double tau;    // unit is second
  };

 public:
  FitAllanAcc(std::vector<double> sigma2s, std::vector<double> taus, double _freq);
  std::vector<double> calcSimDeviation(const std::vector<double> taus) const;
  double getBiasInstability() const;
  double getWhiteNoise() const;

 private:
  std::vector<double> checkData(std::vector<double> sigma2s, std::vector<double> taus);

  std::vector<double> initValue(std::vector<double> sigma2s, std::vector<double> taus);
  double findMinNum(const std::vector<double> num) const;
  int findMinIndex(std::vector<double> num);
  double calcSigma2(double _Q, double _N, double _B, double _K, double _R, double _tau) const;

 public:
  /**
     * @brief getQ
     *          Quantization Noise
     * @unit: m/s
     * @return
     */
  double getQ() const;

  /**
     * @brief getN
     *          Velocity Random Walk
     * @unit: m/s/sqrt(s)
     * @return
     */
  double getN() const;

  /**
     * @brief getB
     *        Bias Instability
     * @unit: m/s^2
     * @return
     */
  double getB() const;

  /**
     * @brief getK
     *      Acceleration Random Walk
     * @unit: m/s^2/sqrt(s)
     * @return
     */
  double getK() const;

  /**
     * @brief getR
     *        Angle Rate Ramp
     * @unit: m/s^3
     * @return
     */
  double getR() const;

  double Q; // Quantization Noise, QN, slope=-1
  double N; // Velocity Random Walk, VRW, slope=-0.5
  double B; // Bias Instability, BI, slope=0
  double K; // Acceleration Random Walk, ARW, slope=0.5
  double R; // Angle Rate Ramp, AR

 private:
  std::vector<double> m_taus; // sample period in Allan
  double freq;
};
} // namespace imu

/*
Interpret the Allan Deviation Plot
  white noise(velocity random walk): -0.5 slope
  bias instability: flat
  flicker noise: positive slope at longer tau
  random walk(acceleration random walk): manifest at the shortest values and contributes to random function in acceleration readings

Derive Key parameters:
  velocity random walk: -0.5 slope
  bias instability: flat
  acceleration random walk: the short tau regions
*/

#endif // FitAllanAcc_H
