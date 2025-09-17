/*
This file defines some useful general functions that are used by the classes
    defined in Jacobi.h, Pinning.h, Superfluid.h and PDE.h
*/

#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include "Vector2.h"
#include <vector>
#include <fstream>

using std::cout;
using std::vector;


// returns p mod(2pi)
double Mod2pi(const double p);
// clockwise rotation of the Vector2 vec by the angle 'angle'
Vector2 Rotation(const Vector2 &vec, const double angle);

// compute the distance of the point x from the point x0 in PBC
//  (dist can not be higher in module than L/2)
double pbcDistance(const double x, const double x0, const double L);

// read input parameter and skip comment
template<typename T>
void ReadPar(std::ifstream& input, T& par, bool print=0) {
    std::string buffer;
    input >> par;
    std::getline(input, buffer);
    if (print) cout << endl << par << endl;
}


#endif  // FUNCTIONS_H_
