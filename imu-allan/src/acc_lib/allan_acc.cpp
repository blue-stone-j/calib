#include "allan_acc.h"

imu::AllanAcc::AllanAcc(std::string name, int maxCluster) :
  m_name(name), numData(0), numCluster(maxCluster)
{
  std::cout << m_name << " "
            << " num of Cluster " << numCluster << std::endl;
}

imu::AllanAcc::~AllanAcc()
{
  m_rawData.clear();
  m_thetas.clear();
  mFactors.clear();
}

void imu::AllanAcc::pushRadPerSec(double data, double time)
{
  m_rawData.push_back(AccData(data * 57.3 * 3600, time));
  numData++;
}

void imu::AllanAcc::pushDegreePerSec(double data, double time)
{
  m_rawData.push_back(AccData(data * 3600, time));
  numData++;
}

void imu::AllanAcc::pushMPerSec2(double data, double time)
{
  m_rawData.push_back(AccData(data, time));
  numData++;
}

void imu::AllanAcc::calc()
{
  std::cout << m_name << " "
            << " numData " << numData << std::endl;
  if (numData < 10000)
    std::cout << m_name << " "
              << " Too few number" << std::endl;

  double start_t = m_rawData.begin()->t;
  double end_t   = m_rawData[numData - 1].t;
  std::cout << m_name << " "
            << " start_t " << start_t << std::endl;
  std::cout << m_name << " "
            << " end_t " << end_t << std::endl;
  std::cout << m_name << " "
            << "dt " << std::endl //
            << "-------------" << (end_t - start_t) << " s" << std::endl
            << "-------------" << (end_t - start_t) / 60 << " min" << std::endl
            << "-------------" << (end_t - start_t) / 3600 << " h" << std::endl;

  if ((end_t - start_t) / 60 < 10)
    std::cout << m_name << " "
              << " Too short time!!!!" << std::endl;

  m_freq = getAvgFreq(); // unit: hz
  std::cout << m_name << " "
            << " freq " << m_freq << std::endl;

  double period = getAvgPeriod(); // unit is second
  std::cout << m_name << " "
            << " period " << period << std::endl;

  m_thetas = calcThetas(m_freq); // unit is m/s

  initStrides();

  mVariance = calcVariance(period); // unit is m^2/s^4=(m/s^2)^2
}

std::vector<double>
imu::AllanAcc::getVariance() const
{
  return mVariance;
}

std::vector<double>
imu::AllanAcc::getDeviation()
{
  double period = getAvgPeriod();

  std::vector<double> sigma2 = calcVariance(period); // unit is m^2/s^4=(m/s^2)^2
  std::vector<double> sigma;

  for (auto &sig : sigma2)
  {
    sigma.push_back(sqrt(sig)); // unit is m/s^2
  }
  return sigma;
}

std::vector<double>
imu::AllanAcc::getTimes()
{
  double period = getAvgPeriod();
  std::vector<double> time(numFactors, 0.0);
  for (int i = 0; i < numFactors; i++)
  {
    int factor = mFactors[i];
    // double clusterPeriod2 = ( period * factor ) * ( period * factor );
    time[i] = period * factor;
  }
  return time;
}

std::vector<int>
imu::AllanAcc::getFactors() const
{
  return mFactors;
}

double
imu::AllanAcc::getFreq() const
{
  return m_freq;
}

std::vector<double>
imu::AllanAcc::calcVariance(double period)
{
  std::vector<double> sigma2(numFactors, 0.0);

  for (int i = 0; i < numFactors; i++)
  {
    int factor            = mFactors[i];
    double clusterPeriod2 = (period * factor) * (period * factor);       // unit is s^2
    double divided        = 2 * clusterPeriod2 * (numData - 2 * factor); // unit s^2
    int max               = numData - 2 * factor;

    for (int k = 0; k < max; k++)
    {
      double temp = (m_thetas[k + 2 * factor] - 2 * m_thetas[k + factor] + m_thetas[k]); // unit is m/s
      sigma2[i] += (temp * temp);                                                        // unit is m^2/s^2
    }

    //            std::cout << "sigma2 " << sigma2[i] << std::endl;
    sigma2[i] = sigma2[i] / divided;
  }

  //        for ( int index = 0; index < mVariance.size( ); ++index )
  //        {
  //            std::cout << sigma2[index] << std::endl;
  //        }
  return sigma2;
}

std::vector<double>
imu::AllanAcc::calcThetas(const double freq)
{
  std::vector<double> thetas;

  double sum = 0;
  for (auto &acc : m_rawData)
  {
    sum += acc.a;
    thetas.push_back(sum / freq);
  }
  return thetas;
}

void imu::AllanAcc::initStrides()
{
  // calculate the largest power of 2 (max_stride) that is less than or equal to (numdata/2)
  int mode                = numData / 2;
  unsigned int max_stride = 1; // unit: dimensionless
  int shift               = 0;
  while (mode)
  {
    mode       = mode >> 1;
    max_stride = 1 << shift;
    shift++;
  }
  /* You can use code below to replace this while loop
  ```
  shift = std::log2(mode);
  max_stride = 1 << shift;
  ```
  However, there are some disadvantages.
  1. Some system doesn't support std::log2
  2. counter error for edge cases such as mode is zero
  3. less flexcity
  4. depend on <cmath>
  */
  // std::cout << m_name << " " << " mode: " << mode << std::endl;
  // std::cout << m_name << " " << " shift: " << shift << std::endl;

  std::vector<double> avgFactors = getLogSpace(0, log10(max_stride));
  // std::cout << m_name << " " << " avgFactors " << avgFactors.size() << std::endl;

  for (int i = 0; i < numCluster; i++)
  {
    // std::cout << m_name << " " << " avgFactors " << i << " " << avgFactors[i] << std::endl;
    avgFactors[i] = ceil(avgFactors[i]);
  }

  std::vector<int> factors(numCluster, 0);

  numFactors = 1;
  factors[0] = (int)avgFactors[0];

  for (int i = 1; i < numCluster; i++)
  {
    // std::cout << m_name << " " << " avgFactors " << i << " " << avgFactors[i] << std::endl;
    if (avgFactors[i] != avgFactors[i - 1])
    {
      factors[numFactors] = (int)avgFactors[i];
      numFactors++;
    }
  }

  mFactors = factors;

  // std::cout << m_name << " " << " numFactors " << numFactors << std::endl;
}

std::vector<double>
imu::AllanAcc::getLogSpace(float a, float b)
{
  std::vector<double> logSpace;

  double start       = pow(10, a); // a always equals to 0 in this class --> start
  double end         = pow(10, b); // considering b=log10(max_stride) --> end=max_stride
  double progression = pow(end / start, (float)1 / (numCluster - 1));

  logSpace.push_back(start);
  for (int i = 1; i < numCluster; i++)
  {
    logSpace.push_back(logSpace[i - 1] * progression); // first is 10^a, last is (10^b)max_stride
  }
  // std::cout << m_name << " " << "size logSpace " << logSpace.size( ) << std::endl;
  return logSpace; // first element is 1, last element is max_stride
}

double
imu::AllanAcc::getAvgDt()
{
  double sum_dt  = 0.0;
  double start_t = m_rawData[0].t; // time of last data, unit is second
  bool first     = true;
  for (auto &acc : m_rawData)
  {
    if (!first)
    { sum_dt += (acc.t - start_t); }
    start_t = acc.t;
    first   = false;
  }
  return sum_dt / (numData - 1); // unit is second
}
