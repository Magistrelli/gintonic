/*
Class that defines the pinning potential as a scalar field
Implemented: cos, gaussian and gaussian rod pinning potentials with an
    arbitrary number of equally spaced, identical peaks
*/

#ifndef PINNING_H_
#define PINNING_H_

#include "./functions.h"

using std::string;
using std::cerr;
using std::endl;

class Pinning {
 private:
    string _type;           // type of the pinning function (cos, gauss, rod)
    bool _ifActivePart;     // if pinning is produced by active particles
    double _str, _sigma;    // pinning peaks: maximum values and widths
    Vector2 * _xp;          // position of the potential's peaks
    double _Lx, _Ly;        // fundamental cell dimensions
    unsigned int _nX, _nY;  // number of pinning peaks along the cartesian axis

 public:
    // constructors
    Pinning();
    Pinning(const string type, const bool ifActivePart,
            const vector<double> pinningParams, vector<unsigned int> npin);
    // destructor
    ~Pinning();
    // variable access
    double GetStr() const        {return _str;}
    double GetSigma() const      {return _sigma;}
    Vector2 * GetXp() const      {return _xp;}
    bool GetIfActivePart() const {return _ifActivePart;}
    void SetNewPos(const Vector2 * const xp);
    // constant and uniform background flow shift
    void BackgroundFlow(const Vector2 shift)  {*_xp += shift;}
    // pinning potential's evaluation
    double EvalPeak(const double x, const double y,
                    const unsigned int jpeak) const;
    double Eval(const double x, const double y) const;
};


#endif  // PINNING_H_
