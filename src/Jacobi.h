/*
Class that implements the Jacobi Theta1(z,q) function and the evaluation of
    the Wood-like (Jacobi-based) superfluid phase field thetaJ(x|x_v)
This phase field is employed by the class Superfluid (see './Superfluid.cpp')
    to define the initial condition of the phase of the superfluid's wave
    function (if at least one vortex is present)

Eqz (7.16) can be employed for this purpose only if the vortex lattice's first
    fundamental vector is l1=(pi,0)
In the general case, the coordinate system must be rescaled and rotated to
    restore this condition
The class Jacobi automatically implements this in their constructor and
    therefore it allows to calculate both the Jacobi theta function (7.17) and
    the phase (7.16) in the case of a general vortex lattice that has only one
    vortex inside the Wigner-Seitz cell

In the Superfluid class, each vortex of the fundamental cell will be considered
    as defining an independent vortex lattice
With this in mind, in Jacobi, we assume that the vortex lattice's fundamental
    vectors coincide with the fundamental cell's sides
*/


#ifndef JACOBI_H_
#define JACOBI_H_

#include "./functions.h"

class Jacobi {
 private:
    Vector2 _q;              // Jacobi Theta1 function's nome, q=exp(i pi tau)
    double _qMod, _qPhase;   // <_q>'s module and phase
    Vector2 _l1, _l2;        // unit cell sides (and vortex lattice fundamental
                             //    vector, see lines 24-32 of this file)
    Vector2 _rotl2;          // l2 in the l1=(pi,0) rescaling
    Vector2 _tau;            // Jacobi Theta's second complex variable,
                             //    Th_J1 = Th_J1(z, tau)
    unsigned int _upSum;     // last iteration of the series definition (7.17)
    double _angle, _lScale;  // <_l1>'s inclination, rescaling factor M_PI/|L1|
    Vector2 _posv;           // referral vortex position: first Jacobi Theta's
                             //    complex variable is z = (x-xv, y-yv)
    Vector2 _rotXv;          // vortex position in the l1=(pi,0) rescaling

    // compute q from tau
    void _qFromTau();
    // check that the referral vortex lies inside the fundamental cell
    void CheckXV() const;
    // find the cell where the position vector pos=(x,y) is
    vector<int> XYcell(const Vector2 &pos) const;

 public:
    // constructors
    Jacobi();
    Jacobi(const Vector2 &l1, const Vector2 &l2);
    Jacobi(const Vector2 &l1, const Vector2 &l2, const Vector2 &posv);
    // destructor
    ~Jacobi() {}
    // Jacobi Theta1(z,q) function
    Vector2 Eval(const Vector2 &rotPos) const;
    // Wood-like (Jacobi-based) superfluid phase thetaJ(x|x_v)
    double Theta(const Vector2 &pos) const;
    // return phase(x,y) in QPBC
    double QPBC(const Vector2 &pos) const;
};

#endif  // JACOBI_H_
