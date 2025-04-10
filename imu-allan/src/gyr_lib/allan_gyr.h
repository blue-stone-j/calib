#ifndef ALLAN_H
#define ALLAN_H

#include "../type.h"
#include <iostream>
#include <math.h>
#include <vector>

namespace imu
{
class AllanGyr
{
 public:
  AllanGyr(std::string name, int maxCluster = 10000);
  ~AllanGyr();

  // receive angular velocity and convert unit to degree/hour
  void pushRadPerSec(double data, double time);
  void pushDegreePerSec(double data, double time);
  void pushDegreePerHou(double data, double time);

  // entrance function to calculate all parameters of gyr
  void calc();

  // get result and won't make influence to result
  std::vector<double> getVariance() const;
  std::vector<double> getDeviation();
  std::vector<double> getTimes();
  std::vector<int> getFactors() const;
  double getAvgValue();
  double getFreq() const;

 private:
  std::vector<double> calcVariance(double period);

  std::vector<double> calcThetas(const double freq);
  void initStrides();
  std::vector<double> getLogSpace(float a, float b);
  // calculate average frequency, unit is Hz
  double getAvgFreq() { return 1.0 / getAvgDt(); }
  // calculate average interval between two consecutive imu msg, unit is second
  double getAvgDt();
  // calculate sample period of imu, unit is second
  double getAvgPeriod() { return getAvgDt(); }
  int getFactorsNum() { return numFactors; }

  std::string m_name;             // sensor name
  double m_freq;                  // average frequency of imu msg
  int numData;                    // num of buffered imu msg
  std::vector<GyrData> m_rawData; // buffered data, unit is degree/hour and second
  std::vector<double> m_thetas;   // accumulated angle from begining to current imu msg, unit is arcsecond
  int numCluster;                 // how many points of Allan curve line we want to obtain
  int numFactors;                 // num of all strides
  std::vector<int> mFactors;      // stride for different periods(virtual periods)

  std::vector<double> mVariance;
};
} // namespace imu

#endif // ALLAN_H
