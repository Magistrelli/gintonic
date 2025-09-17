/*
Class for solving time Partial Differential Equations
The Runge-Kutta (4-th order), the Adam-Bashforth (3-rd order) and the Euler
    methods are implemented in the case of a generic Equation for the time
    derivative of the searched solution
The actual form of the Equation must be specified in each derived class

IMPORTANT! Quasi-periodic complex solution are considered: the solution is
    looked for on the extended grid (both the edges of the fundamental cell
    are considered independently, see also './Superfluid.*' for the
    introduction of the "extended" fundamental cell)

GPE implements the explicit form of the Gross-Pitaevskii Equation in the case
    of a single spinless superfluid
It inherits from PDE the time evolution algorithms and from Superfluid the
    physical quantities and the output-generating functions

In PDE.*, [2] refers to the Runge-Kutta reference
    https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#The_Runge%E2%80%93Kutta_method
and [3] to the Adam-Bashforth reference
    https://en.wikipedia.org/wiki/Linear_multistep_method#Adams%E2%80%93Bashforth_methods
*/


#ifndef PDE_H_
#define PDE_H_

#include "./Superfluid.h"

class PDE {
 protected:
    ComplexFunc _XY;        // PDE's complex solution
    double *_params;        // PDE's parameters
    unsigned int _npar;     // number of PDE's parameters
    int _llenx, _lleny;     // number of points of the spatial grid
    double _d_t;            // time step
    int _tSkip;             // out parameter: write some data every _tSkip steps
    unsigned int _tburn;    // burn-in before the actual time evolution
    unsigned int _ntburn;   // number of burn-in steps before the actual time evolution
    double _gamma;          // dissipation parameter
    bool _ifGammaDyn;       // wheter to keep dissipation active during the dynamics

    // 4-th order Runge-Kutta algorithm for dXY/dt=der
    void RKstep(ComplexFunc& XY, ComplexFunc& der, ComplexFunc* k);
    // 3-th order Adam-Bashforth algorithm for dXY/dt=der
    void ABstep(ComplexFunc& der, ComplexFunc* ff);
    // 3-th order Adam-Bashforth algorithm initialization
    void ABinit(ComplexFunc& der, ComplexFunc* ff, const unsigned int order,
                const unsigned int delta, const string identifier = "");
    // actual differential equation
    virtual void Equation(ComplexFunc& der, const ComplexFunc& XY) const = 0;
    // general time evolution operations (e.g. outputs, conserved quantities)
    virtual void EvStuff(const unsigned int t, const unsigned int delta,
                                               const string identifier = "") {}

 public:
    // constructors
    PDE();
    PDE(const int lenx, const int leny, const double dt);
    PDE(const int lenx, const int leny, const double dt,
        const unsigned int npar);
    // destructor
    virtual ~PDE() {delete[] _params;}
    // parameters initialization
    virtual void InitParams(double* params) {}
    // access to interval variables
    void SetTSkip(const int tSkip) {_tSkip = tSkip;}
    // methods that solves the PDE
    void Solver(const double tend, const unsigned int delta,
                const unsigned int method = 0, const string identifier = "");
};



class GPE: public PDE, public Superfluid {
 private:
    double _Omega0;         // mean superfluid angular velocity
    double *_KelX, *_KelY;  // coordinates of the points of the circuit used
                            //  to check Kelvin circulation theorem
    bool _ifKelvin;         // =1 for testing Kelvin theorem

    // common part for all the constructors
    void AllConstr();
    // Gross-Pitaevskii Equation
    virtual void Equation(ComplexFunc& der, const ComplexFunc& XY) const;
    // general time evolution operations (e.g. outputs, conserved quantities)
    virtual void EvStuff(const unsigned int t, const unsigned int delta,
                                               const string identifier = "");
    void PinningEvolve();
    // evolves a closed circuit to check the Kelvin theorem
    void KelvinCheck();

 public:
    // default constructor
    GPE();
    // constructors inherited from Superfluid
    //  (no _params needed for GPE, they are all saved in Superfluid):
    // constructor for a superfluid with a vortex lattice
    explicit GPE(const string inVort): PDE(), Superfluid(inVort) {AllConstr();}
    // constructor for a superfluid with a (couple of) soliton waves
    GPE(const double x01, const double x02, const double nfact): \
                        PDE(), Superfluid(x01, x02, nfact) {AllConstr();}
    // constructor for a superfluid with random density or phase
    GPE(const bool rndN, const bool rndPhase): \
                        PDE(), Superfluid(rndN, rndPhase) {AllConstr();}
    // destructor
    virtual ~GPE();

    // initialize the variables needed for Kelvin theorem's test
    void SetKelvin();
};


#endif  // PDE_H_
