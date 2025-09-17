/*
This library defines the following classes:
- SupCons     ->  Class that gather all the evaluated superfluid's conserved
                    quantities: total particle number, (mean) total energy and
                    total momentum along x and y
- Superfluid  ->  Class for the study of the configurations of a single spinless
                    superfluid
                  All the physical quantities and the timespace grid properties
                    are defined within this class

In main.C, if one wants to study the time evolution of the superfluid's
    configuration, one should use the class GPE, which inherits from Superfluid
    the physical variables and also implements PDE solvers
*/


#ifndef SUPERFLUID_H_
#define SUPERFLUID_H_

#include <fstream>
#include <iomanip>
#include <vector>
#include "./Jacobi.h"
#include "./Pinning.h"
#include "./Random.h"

#define hbar 6.58212e-22    // Planck reduced constant in units of Mev s
#define mass 1.87913e3      // 2.*neutron_mass*c^2 in units of Mev
#define c_light 2.99792e10  // speed of light in units of cm/s
#define PARAM_file std::string("parameters.dat")    // parameter file



class SupCons {
 public:
    // internal variables: total particle number, mean total particle energy,
    //  circulation along the cell boundaries and mean particle momentum
    double totN, partE, circ;
    Vector2 partP;
    // default constructor: initializes all the constants to 0.
    SupCons(): totN(0.), partE(0.), circ(0.), partP() {}
    // destructor
    ~SupCons() {}
    // assignment operator overload
    SupCons& operator= (const SupCons &obj) {
        totN = obj.totN;
        partE = obj.partE;
        partP = obj.partP;
        circ = obj.circ;
        return *this;
    }
};



class Superfluid {
 protected:
    ComplexFunc _Psi;       // superfluid wave function, real and imaginary parts: (_lenx+1)*(_leny+1) complex numbers
    double _nmean, _mu;     // mean (if no vortices) density and chemical potential
    double _t0, _csi;       // time scale (in s) and coherence length (in cm)

    Vector2 _vBg;           // background flow velocity, i.e. velocity shift of the background potential
    Vector2 _vBg0;          // order of magnitude of vBg
    double _vBgT;     // time scale of vBg prescribed evolution
    string _vBgEvol;        // type of prescribed evolution for vBg

    SupCons _cons;          // total density, energy, circulation along the cell edges and momentum
    int _Nv;                // number of vortices-antivortices (considering possible multiple charge)
    unsigned int _NvAll;    // number of initial vortices+antivortices (not considering possible multiple charge)
    Vector2 *_posv;         // vortices initial positions (do not evolve)
    int *_iXv, *_iYv;       // grid indexes for discretized vortex position (not bounded to the fundamental cell)
    int *_vcharge;          // list of the charges of the initial vortices (does not evolve)

    Pinning *_pinning;      // pinning landscape _lenx*_leny
    double * _pinEval;      // specific evaluation of the pinning potential
    double _mpart;          // mass of the active particles [Mev/c2], if present
    unsigned int _npart;    // number of active particles, if present
    string _forceType;      // type of the external force (elastic, constant, null)
    double _forceModule;    // elastic constant [Mev/csi], if active particles
    Vector2 * _xpEq;        // equilibrium positions of the active particles, if present 
    Vector2 * _xpart;       // positions of the active particles, if present
    Vector2 * _vpart;       // velocity of the active particles, if present
    Vector2 * _apart;       // acceleration of the active particles, if present
    Vector2 * _fpart;       // external force acting on the active particles, if present

    double _L1, _L2;        // (Cartesian) unit cell's dimensions
    int _lenx, _leny;       // number of points of the spatial grid
    double *_x, *_y;        // space grid (does not evolve)
    double _dx, _dy, _dt;   // space and time resolution

    double _c1, _c2;        // QPBC integration constants (they evolve)

    vector<int> _mmx;       // outputs-printed smaller cell's x-indices limits
    vector<int> _mmy;       // outputs-printed smaller cell's y-indices limits
    unsigned int _outGrid;  // 1/fraction of the spatial grid printed in ouputs

    // common part of all the constructors
    void AllConstr();
    // initialization of the pinning background
    void InitPinning(const string pinType, const bool ifActivePart,
                const vector<double> pinParams, vector<unsigned int> npin);
    // evaluate the pinning potential
    void PinningEval() const;
    // read the initial vortices positions and their charges
    void ReadVortices(const string inVort);
    // compute the QPBC's _c1 and _c2
    void ConstC();
    // update the values on the extra side with QPBC
    void UpdateSide();
    // update the background potential velocity
    void UpdateVBg(const double time);

    // read a real f(_x[i],_y[j]) of type "len" in true PBC
    double PBC(const double* f, const int i, const int j) const;
    // returns a complex psi(x[i],y[j]) in QPBC
    Vector2 RI_PBC(const ComplexFunc& psi, const int i, const int j,
                                              bool extraSide = 0) const;
    // brings i and j in back to the fundamental cell
    void CountPBC(int &countx, int &county, const int i, const int j,
                  int &ii, int &jj, bool extraSide = 0) const;

    // compute the <order>-th order partial derivative of the wave function
    //  <psi> in (x[<i>],y[<j>]) in the direction x (<axis>=0) or y (<axis>=1),
    //  3-points or 5-points stencil (set by <nStencil>)
    Vector2 Derivative(const ComplexFunc& psi, const int order,
                       const int axis, const int i, const int j,
                       const unsigned int nStencil = 3) const;
    // compute the Laplacian of the wave function <psi> in (x[<i>],y[<j>]),
    //  9-points stencil
    Vector2 Lap9(const ComplexFunc& psi, const int i, const int j) const;

    // compute the phase of the wave function in (x[<i>],y[<j>]), QPBC used
    double Phase(const int i, const int j) const \
                    {return Mod2pi(RI_PBC(_Psi, i, j).Phase());}
    // compute the superfluid density in (_x[i],_y[j]), QPBC used
    double Dens(const int i, const int j) const \
                    {return pow(RI_PBC(_Psi, i, j).Mod(), 2);}
    //returns real part in (x[<i>],y[<j>])
    double Real(const int i, const int j) const \
                    {return RI_PBC(_Psi, i, j).GetRe();}
    //returns imaginary part in (x[<i>],y[<j>])
    double Img(const int i, const int j) const \
                    {return RI_PBC(_Psi, i, j).GetIm();}
    // superfluid (no bg flow) current density in (x[<i>],y[<j>])
    Vector2 Current(const int i, const int j) const;
    // superfluid (no bg flow) velocity in (x[<i>],y[<j>]) from
    //  its current density
    Vector2 Velocity(const int i, const int j) const \
                    {return Current(i, j)/Dens(i, j);}
    // superfluid conserved quantities
    SupCons Conserved() const;
    // total standard deviation of the density field
    double NStdDev() const;

    // find the non-QPBC-ed k-th vortex position indices
    void VortexPos(const unsigned int k);

 public:
    // default constructor
    Superfluid();
    // constructor for a superfluid with a vortex lattice
    explicit Superfluid(const string inVort);
    // constructor for a superfluid with random density or phase
    Superfluid(const bool rndN, const bool rndPhase);
    // constructor for a superfluid with a (couple of) soliton waves
    Superfluid(const double x01, const double x02, const double nfact);
    // destructor
    ~Superfluid();
    // set the background flow velocity
    void SetVbg(const Vector2 &vBg) {_vBg = vBg;}
    void SetVbgT(const double tScale) {_vBgT = tScale;}
    void SetVbgEvol(const string evType) {_vBgEvol = evType;}
    // access to private variables
    double* Getn() const;
    double* GetPhase() const;
    double GetDt() const {return _dt;}
    // read the non-QPBC-ed k-th vortex position vector from _iXv, iYv
    Vector2 GetPosV(const unsigned int k) const;
    // read the k-th vortex starting position
    Vector2 GetIniPosV(const unsigned int k) const {return _posv[k];}

    // set the size of the output cell and its grid
    void SetOutCell(const double fracCell, const unsigned int outGrid);
    // print out the phase, density and total (vortex-generated + bg flow)
    //  velocity fields in natural units
    void OutputPhase(const string fname, double time) const;
    // output of the pinning landscape, physical units (MeV)
    void OutputPin(const string fname) const;
    // extends OutputPhase also outside the fundamental cell (<n>x<m> table of
    //  unit cells) employing QPBC
    void PrintAllGrid(const int n, const int m, const string fname) const;
    // extends the spatial grid of the fundamental cell to that one of a <n>x<m>
    //  table of unit cells centered on the fundamental one
    void LargeGrid(const int n, const int m, double* &x, double* &y) const;
    // print c1, c2 computed independently all over the fundamental cell edges
    void OutCs(const string fname) const;
    // output of the position, velocity, accelerations of the active particles
    // (if present) and force acting on them; natural units
    void OutputParticles(const string fname) const;
};


#endif  // SUPERFLUID_H_
