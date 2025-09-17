#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>

using namespace std;

#ifndef __Random__
#define __Random__

class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;
  int ParSize,ParRank;		//for parallel coding

protected:

public:
  // constructors
  Random(){ParSize=1,ParRank=0;}
  // destructor
  ~Random(){}
  // initializator
  void SetRandom(int* , int, int);
  void SetParallel(int size, int rank) {ParSize=size,ParRank=rank;}
  void SetRandom(string, string);
  //save l1,...,l4 for restart
  void SaveSeed();
  //gen from uniform distributions
  double Rannyu(void);
  double Rannyu(double min, double max);
  //gen from gaussian distribution (Box-Muller)
  double Gauss(double mean, double sigma);
  //gen from different distributions (inverse of cumulative)
  double Exp(double lambda);
  double Lorentz(double mu, double Gamma);
  //gen uniformly an angle in [-pi,pi]
  double Angle();
  //gen uniformly a solid angle, theta in [0,pi] and phi in [-pi,pi]
  void SolidAngle(double& theta, double& phi);
};

#endif // __Random__
