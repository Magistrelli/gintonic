#include "./Superfluid.h"


// CONSTRUCTORS AND DESTRUCTOR

// Default constructor: initializes all the numerical variables to 0 or 0.
//  and all the pointers to NULL
Superfluid::Superfluid(): _Psi(), _nmean(0.), _mu(0.), _t0(0.), _csi(0.), \
        _vBg(), _vBg0(), _vBgT(-1.), _vBgEvol(""), _cons(), _Nv(0), _NvAll(0), \
        _posv(NULL), _iXv(NULL), _iYv(NULL), _vcharge(NULL), _pinning(NULL), \
        _pinEval(NULL), _mpart(NAN), _npart(0), _forceType(""), \
        _forceModule(NAN), _xpEq(NULL), _xpart(NULL), _vpart(NULL), \
        _apart(NULL), _fpart(NULL), _L1(0.), _L2(0.), _lenx(0), _leny(0), \
        _x(NULL), _y(NULL), _dx(0.), _dy(0.), _dt(0.), _c1(0.), _c2(0.), \
        _mmx(2), _mmy(2), _outGrid(0) {}


void Superfluid::AllConstr() {
    /*
    Method that implements several operations that are common to all the
        superfluid's constructors
    It reads the global input file, allocates all the internal dynamic arrays
        that are not related to the vortices and defines the simulation grid
        and the pinning potential
    See './README.md' for the appropriate formatting of the global input file
    */

    ifstream input;
    string pinType, buffer;
    bool ifActivePart=0;
    vector<double> pinParams(6);    // pinning parameters
    vector<unsigned int> npin(2);   // number of the pinning potential's peaks
    double vxBg=0., vyBg=0.;        // background velocity
    double xmin, xmax, ymin, ymax;  // edges of the fundamental cell
    double CourNum;                 // CFL Courant number

    // read the global input file (information about pinning and spacetime grid)
    input.open(PARAM_file);
    ReadPar(input, buffer);         // comment
    ReadPar(input, _nmean);
    ReadPar(input, _mu);
    ReadPar(input, vxBg);
    ReadPar(input, vyBg);
    ReadPar(input, buffer);         // comment
    ReadPar(input, _L1);
    ReadPar(input, _L2);
    ReadPar(input, _lenx);
    ReadPar(input, _leny);
    ReadPar(input, buffer);         // comment
    ReadPar(input, _dt);
    ReadPar(input, buffer);         // tend
    ReadPar(input, buffer);         // tout
    ReadPar(input, buffer);         // comment
    ReadPar(input, pinType);
    ReadPar(input, ifActivePart);
    if (ifActivePart) {
        cout << endl << "ERROR! Active particles are NOT tested!" << endl;
        exit(-999);
    }
    for (unsigned int i = 0; i < 4; i++) {ReadPar(input, pinParams[i]);}
    for (unsigned int i = 0; i < 2; i++) {ReadPar(input, npin[i]);}
    input.close();
    
    _t0 = hbar/_mu;
    _csi = hbar*c_light/sqrt(mass*_mu);
    _vBg = Vector2(vxBg, vyBg);
    _vBg0 = _vBg;
    if ((pinParams[2] > _L1/2.) || (pinParams[3] > _L2/2.)) {
        cerr << \
        "center of the pinning potential outside of the fundamental cell!" \
        << endl;
        exit(-10);
    }
    pinParams[4] = _L1;
    pinParams[5] = _L2;

    // vectors dynamic allocation
    // _Psi is evaluated twice on the unit cell's edges (e.g. both in x=-L1/2
    //   and in x=+L1/2), while _pinning is only evaluated once (e.g. only in
    //   x=-L1/2): _pinEval(L1/2) = _pinEval(-L1/2) because of the (true) PBC,
    //   while the difference between _Psi(L1/2) and _Psi(-L1/2) will fix the
    //   integration constants _c1 and _c2 of the QPBC
    // Correspondingly, the +1 in the lengths of _x and _y is needed to define
    //   this "extended" cell
    _x = new double[_lenx+1];
    _y = new double[_leny+1];
    _Psi.SetLen((_lenx+1)*(_leny+1)); // on the extended unit cell

    // define the simulation grid, x=y=0 is the center of the unit cell
    xmax = 0.5*abs(_L1);
    ymax = 0.5*abs(_L2);
    xmin = -xmax;
    ymin = -ymax;
    _dx = (xmax-xmin)/_lenx;
    _dy = (ymax-ymin)/_leny;
    if (_dx-_dy < 1e-15) {
        _dy = _dx;  // removing numerical errors from input operations
        cout << endl << \
        "The 9-stencil Laplacian will be used for the spatial derivatives calculation" \
        << endl;
    } else {        // Laplacian needs squared grid elements
        cout << endl << "_dx-_dy = " << _dx-_dy << ": " << \
        "The Laplacial will be computed summing the 5-stencil 2-nd derivatives" \
        << endl;
    }
    // define the coordinate of each point of the spatial grid
    for (int i = 0; i < _lenx+1; i++) {_x[i] = xmin + i*_dx;}
    for (int j = 0; j < _leny+1; j++) {_y[j] = ymin + j*_dy;}

    InitPinning(pinType, ifActivePart, pinParams, npin);

    // by default, print in outputs the whole cell and grid
    SetOutCell(1., 1);

    // if dt < 0, set dt such that CourNum=0.1, otherwise check CFL condition
    CourNum = sqrt(_mu) * _dt * (1./_dx + 1./_dy);
    if (_dt < 0.) {
        if (_dt < -1.) {
            cerr << "\nERROR! CFL condition not satisfied, " \
                    << "Courant number = " << CourNum << " > 1!";
            exit(-12);
        }
        CourNum = -_dt;
        _dt = CourNum / sqrt(_mu) / (1./_dx + 1./_dy);
        CourNum = sqrt(_mu) * _dt * (1./_dx + 1./_dy);
        cout << "\nSet dt = " << _dt << " with CFL condition " \
                << "and Courant number Cnum = " << CourNum << endl;
    } else if (CourNum > 1.) {
        cerr << "\nERROR! CFL condition not satisfied, " \
                << "Courant number = " << CourNum << " > 1!";
        cerr << "  Decrease space resolution or increase time resolution!\n" \
                << endl;
        exit(-12);
    } else if (CourNum > 0.2) {
        cout << "\nWarning: CFL ok, but Courant number Cnum = " << CourNum \
                << " > 0.2 might cause numerical instability." << endl;
        cout << "  You can decrease Cnum by reducing space resolution " \
                << "or increasing time resolution\n" << endl;
    }

}


void Superfluid::InitPinning(const string pinType, const bool ifActivePart,
                const vector<double> pinParams, vector<unsigned int> npin) {
    // not generalized yet to the case of an arbitrary number of active particles
    Vector2 displ;  // initial displacement of particles from equilibrium
    Vector2 inivel; // initial velocity of the particles
    _pinEval = new double[_lenx*_leny]; // on the std unit cell
    double buffer;

    if (ifActivePart) {
        ifstream inPart;
        inPart.open("activePart.dat");
        inPart >> _npart;
        // check inputs
        if (_npart != npin[0]*npin[1]) {
            throw std::runtime_error("number of active particles (in \
                        activePart.dat) not compatible with the number of \
                        potential's peaks (in "+PARAM_file+")");
        }

        inPart >> _mpart;
        inPart >> _forceType;
        inPart >> _forceModule;
        // check input concistency
        if (_forceType == "elastic")
          cout << "An elastic force F = " << _forceModule << "* (x-x0) Mev is \
                acting on the active particles" << endl;
        else if (_forceType == "constant")
          cout << "A consntant force F = " << _forceModule << " MeV is acting \
                on the active particles" << endl;
        else if (_forceType == "null")
          cout << "No external force is acting on the active particles" << endl;
        else
            throw std::runtime_error("The required external force is not \
                                      implemented!");

        for (int comp = 0; comp < 2; comp++) {
            inPart >> buffer;
            displ.SetComp(comp, buffer);
            inPart >> buffer;
            inivel.SetComp(comp, buffer);
        }
        inPart.close();
    }

    _pinning = new Pinning(pinType, ifActivePart, pinParams, npin);

    if (ifActivePart) {
        _xpEq = new Vector2[_npart];
        _xpart = new Vector2[_npart];
        _vpart = new Vector2[_npart];
        _apart = new Vector2[_npart];
        _fpart = new Vector2[_npart];
        for (unsigned int ipart = 0; ipart < _npart; ipart++) {
            _xpEq[ipart] = _pinning->GetXp()[ipart];
            _xpart[ipart] = _xpEq[ipart] + displ;
            _vpart[ipart] = inivel;
        }
        _pinning->SetNewPos(_xpart);
    }

    PinningEval();
}


Superfluid::Superfluid(const string inVort): Superfluid() {
    /*
    Initializes all the physical quantities for a superlfuid system that shows
        a periodic array of vortices of arbitrary charges
    
    The initial density accounts for the deplition due to the pinning potential
        and the vortex cores
    The pinning-related deplition is accounted by corretcing the mean density
        with the multiplicative factor (1 - V/mu)
    
    Each vortex of the fundamental cell is considered as defining an
        independent vortex lattice and its contribution to the phase field is
        computed within a <Jacobi> instance separately
    The initial phase field is given by the sum of the contributions of the
        single vortices, computed generalizing the Wood-like phase to the case
        of arbitrary lattice's fundamental vectors (see './Jacobi.*' for the
        implementation of this operation)
    The inital phase, i.e. the imaginary part of the argument of the
        exponential of the wave function, does NOT include the contribution of
        the background flow (which is instead accounted for in the explicit
        form of the Gross-Pitaevskii equation, see './PDE.cpp.)
    
    :arg: <inVort>  name of the input file containing the initial position and
                    charges of the vortices
    */

    double* phase, *pet_f;      // phase and square modulus of _Psi(t=0)
    Jacobi *iniPhase = NULL;
    double depl, rho, appo;
    int gpos;
    double distx, disty;

    // read the global input file, dynamic allocate all the useful array,
    //  define the simulation grid and the pinning potential
    AllConstr();
    // dynamic allocate the phase field on the "extended" grid (see AllConstr())
    phase = new double[(_lenx+1)*(_leny+1)];
    // dynamic allocate the density field on the "usual" grid (see AllConstr())
    pet_f = new double[_lenx*_leny];
    // read the vortex input file
    ReadVortices(inVort);

    // define single-vortex lattice's fundamental vectors
    Vector2 l1(_L1, 0.);
    Vector2 l2(0., _L2);
    // compute the initial phase
    for (unsigned int k = 0; k < _NvAll; k++) {
        iniPhase = new Jacobi(l1, l2, _posv[k]);
        for (int i = 0; i < _lenx+1; i++) {
            for (int j = 0; j < _leny+1; j++) {
                gpos = (_leny+1)*i + j;  // (x[i],y[j]) on the extended grid
                if (k == 0) {phase[gpos] = 0.;}
                // compute and then sum the single vortex terms
                appo = iniPhase->Theta(Vector2(_x[i], _y[j]));
                phase[gpos] += _vcharge[k]*appo;
                if (k == _NvAll-1) {phase[gpos] = Mod2pi(phase[gpos]);}
            }
        }
        delete iniPhase;
    }

    // compute the initial density
    for (int i = 0; i < _lenx; i++) {
        for (int j = 0; j < _leny; j++) {
            gpos = _leny*i+j;   // (x[i],y[j]) on the usual grid
            pet_f[gpos] = 1.;
            for (unsigned int k = 0; k < _NvAll; k++) {
                // compute the distance in PBC centered on the vortex core
                distx = pbcDistance(_x[i], _posv[k][0], _L1);
                disty = pbcDistance(_y[j], _posv[k][1], _L2);
                rho = sqrt(pow(distx, 2) + pow(disty, 2));  // r of (7.19)
                pet_f[gpos] *= sqrt(_nmean) * \
                               rho/sqrt(2.*pow(_vcharge[k], 2) + pow(rho, 2));
            }
            // account for the pinning-related depletion
            depl = 1. - _pinEval[gpos]/_mu;
            //if(i ==0 and j ==0){
            //    cout << _pinEval[gpos];
            //}
            if (depl < 0.) {
                cerr << "\nThe pinning is too strong (total depletion!)" \
                     << endl;
                exit(-18);
            }
            pet_f[gpos] *= sqrt(depl);
        }
    }

    // copmute the intial wave function's real and imaginary parts
    for (int i = 0; i < _lenx+1; i++) {
        for (int j = 0; j < _leny+1; j++) {
            gpos = (_leny+1)*i + j;     // (x[i],y[j]) on the extended grid
            _Psi.SetValue(gpos, PBC(pet_f, i, j)*cos(phase[gpos]),
                                PBC(pet_f, i, j)*sin(phase[gpos]));
        }
    }
    delete[] phase;
    delete[] pet_f;

    // initialize _c1 and _c2 and the conserved quantities
    ConstC();
    _cons = Conserved();

    // print time and space physical scales
    cout << endl << "time scale: " << _t0 \
         << " s, \t coherence lenght: " << _csi << " cm" << endl;
}


Superfluid::Superfluid(const bool rndN, const bool rndPhase): Superfluid() {
    /*
    Initializes all the physical quantities for a superlfuid system that is
        generated with random or uniform density and random or uniform phase
    Since we are dealing with the factorized wave function psi^tilde (i.e. the
        wave function without the background flow term), the QPBC integration
        constants must always be 0 (no vortex is present)
    
    The initial density is computed by picking in each point of the (usual)
        unit cell's grid a random number in (0.9 n_mean, 1.1 n_mean)
    The pinning-related deplition is accounted by corretcing the mean density
        with the multiplicative factor (1 - V/mu)
    
    The initial phase field is computed by picking in each point of the (usual)
        unit cell's grid a random number in (0.9 M_PI, 1.1 M_PI)
    The inital phase does NOT include the contribution of the background flow
        (which is instead accounted for in the explicit form of the
        Gross-Pitaevskii equation, see './PDE.cpp.)
    
    The wave function on the extra sides of the extended unit cell is computed
        in QPBC with c_i = 0 (i.e. in true PBC) 
    
    :arg: <rndN>        = 0 to initialize the system with uniform density,
                        = 1 to initialize the system with random density
    :arg: <rndPhase>    = 0 to initialize the system with a uniform phase field,
                        = 1 to initialize the system with a random phase field
    */

    int gpos;
    double depl;
    Random rnd;
    rnd.SetRandom("src/Primes", "src/seed.in");  // initialize random_gen from files

    // read the global input file, dynamic allocate all the useful array,
    //  define the simulation grid and the pinning potential
    AllConstr();
    // dynamic allocate the phase and density on the "usual" grid
    double *dens = new double[_lenx*_leny];
    double *phase = new double[_lenx*_leny];
    _Nv = 0.;

    // compute the initial density on the "usual" unit cell
    if (rndN) {
        for (int i = 0; i < _lenx*_leny; i++) {
            dens[i] = rnd.Rannyu(0.9*_nmean, 1.1*_nmean);
        }
    } else {
        for (int i = 0; i < _lenx*_leny; i++) {
            dens[i] = _nmean;
        }
    }
    // correct for the pinning-related depletion guess
    for (int i = 0; i < _lenx; i++) {
        for (int j = 0; j < _leny; j++) {
            depl = 1. - _pinEval[_leny*i+j]/_mu;
            if (depl < 0.) {
                cerr << "\nThe pinning is too strong (total depletion!)" \
                << endl;
                exit(-18);
            }
            dens[_leny*i+j] *= depl;
        }
    }

    // compute the inital phase on the "usual" unit cell
    if (rndPhase) {
        for (int i = 0; i < _lenx*_leny; i++) {
            phase[i] = rnd.Rannyu(0.9*M_PI, 1.1*M_PI);
        }
    } else {
        for (int i = 0; i < _lenx*_leny; i++) {
            phase[i] = M_PI;
        }
    }

    // initialize the wave function's real and imaginary parts
    for (int i = 0; i < _lenx; i++) {
        for (int j = 0; j < _leny; j++) {
            gpos = _leny*i+j;
            _Psi.SetValue((_leny+1)*i+j, sqrt(dens[gpos]) * cos(phase[gpos]),
                                         sqrt(dens[gpos]) * sin(phase[gpos]));
        }
    }
    delete[] phase;
    delete[] dens;

    // extend the wave function on the extended unit cell using null cs
    //  (even with a background flow, as we are considering psi^tilde)
    _c1 = 0.;
    _c2 = 0.;
    UpdateSide();

    // compute the initial conserved quantities
    _cons = Conserved();
}


Superfluid::Superfluid(
        const double x01, const double x02, const double nfact): Superfluid() {
    /*
    Initializes all the physical quantities of a superlfuid system that is
        started with a couple of soliton waves propagating in the opposite
        horizontal directions
    Since we are dealing with the factorized wave function psi^tilde (i.e. the
        wave function without the background flow term), the QPBC integration
        constants must always be 0 (no vortex is present)
    
    The initial wave function is computed via Pethick's (7.159)
    The pinning-related deplition is accounted by corretcing the mean density
        with the multiplicative factor (1 - V/mu)
    The inital phase does NOT include the contribution of the background flow
        (which is instead accounted for in the explicit form of the
        Gross-Pitaevskii equation, see './PDE.cpp.)
    
    The wave function on the extra sides of the extended unit cell is computed
        in QPBC with c_i = 0 (i.e. in true PBC) 

    :arg: <x01>     initial position of the right-moving soliton wave
    :arg: <x02>     initial position of the left-moving soliton wave
    :arg: <nfact>   the maximum deplition produced by a single soliton gives a
                    minimum density of nmin = <nfact> * <_nmean>
    */

    double depl, nmin, deltan;
    double R1, R2, I1, I2;

    // read the global input file, dynamic allocate all the useful array,
    //  define the simulation grid and the pinning potential
    AllConstr();
    if ((x01 > _L1/2.) || (x01 < -_L1/2)) {
        cerr << "First wavefront outside the fundamental cell!" << endl;
        exit(-11);
    }
    if ((x02 > _L1/2.) || (x02 < -_L1/2)) {
        cerr << "Second wavefront outside the fundamental cell!" << endl;
        exit(-11);
    }
    nmin = _nmean*nfact;
    deltan = _nmean-nmin;
    _Nv = 0.;

    // compute the initial wave function on the "usual" unit cell
    for (int i = 0; i < _lenx; i++) {
        for (int j = 0; j < _leny; j++) {
            R1 = sqrt(deltan) * tanh((_x[i]-x01) * sqrt(0.5*deltan/_nmean));
            R2 = sqrt(deltan) * tanh((x02-_x[i]) * sqrt(0.5*deltan/_nmean));
            I1 = sqrt(nmin);
            I2 = sqrt(nmin);
            depl = 1. - PBC(_pinEval, i, j)/_mu;
            if (depl < 0.) {
                cerr << "\nThe pinning is too strong (total depletion!)" \
                << endl;
                exit(-18);
            }
            _Psi.SetValue((_leny+1)*i + j, (R1*R2 - I1*I2) * sqrt(depl),
                                           (R1*I2 + I1*R2) * sqrt(depl));
        }
    }
    // extend the wave function on the extended unit cell using null cs
    //  (even with a background flow, as we are considering psi^tilde)
    _c1 = 0.;
    _c2 = 0.;
    UpdateSide();

    // compute the initial conserved quantities
    _cons = Conserved();
}

// destructor
Superfluid::~Superfluid() {
    delete[] _posv;
    delete[] _iXv;
    delete[] _iYv;
    delete[] _vcharge;
    if (_pinning->GetIfActivePart()) {
        delete[] _xpEq;
        delete[] _xpart;
        delete[] _vpart;
        delete[] _apart;
        delete[] _fpart;
    }
    delete _pinning;
    delete[] _pinEval;
    delete[] _x;
    delete[] _y;
}





// PRIVATE FUNCTIONS, first chunck

void Superfluid::PinningEval() const {
    /*
    Evaluate a specific instance of the pinning potential <_pinning>
    Saves the result in the dynamic array pointed by the internal variable
        <_pinEval>
    Impose that the integral of the background potential is 0 over the whole
        unit cell
    */

    double intV = 0.;   // integral of the potential all over the unit cell
    for (int i = 0; i < _lenx; i++) {
        for (int j = 0; j < _leny; j++) {
            _pinEval[_leny*i + j] = _pinning->Eval(_x[i], _y[j]);
            intV += _pinEval[_leny*i + j];
        }
    }
    intV /= _lenx*_leny;
    for (int ij = 0; ij < _lenx*_leny; ij++) {_pinEval[ij] -= intV;}
}



void Superfluid::ReadVortices(const string fname) {
    /*
    Reads the vortex input file (containing the vortices' initial positions
        and charges) and allocates all the internal dynamic arrays that are
        related to the vortices
    See './README.md' for the appropriate formatting of the vortex input file
    
    :arg: <fname>   name of the input file containing the initial positions and
                    charges of the vortices
    */

    ifstream input;
    unsigned int nlines = 0;    // numbers of the input file's lines
    string oneline;
    double xv, yv, charge_val;  // components of the vortex positions

    // count the number of vortex centers
    input.open(fname);
    while (getline(input, oneline)) {nlines++;}
    input.close();
    if (nlines <= 2) {  // file's first line is the legend, second is empty
        cerr << "The vortices file is empty!" << endl;
        exit(-12);
    }

    _NvAll = nlines-2;
    // delete vortex vectors if they already exist
    if (_posv) {delete[] _posv;}
    if (_vcharge) {delete[] _vcharge;}
    // allocate the vortex vectors
    _posv = new Vector2[_NvAll];
    _iXv = new int[_NvAll];
    _iYv = new int[_NvAll];
    _vcharge = new int[_NvAll];

    // read the input file
    input.open(fname);
    getline(input, oneline);    // discard the first line (legend)
    for (unsigned int k = 0; k < _NvAll; k++) {
        // read input line with check 
        if (!(input >> xv >> yv >> charge_val)) {
            std::cerr << "Error: Could not read vortex data at index " << k << std::endl;
            exit(EXIT_FAILURE);
        }

        // check that the k-th vortex is inside the fundamental cell
        if ((xv < -0.5*abs(_L1)) || (xv >= 0.5*abs(_L1))) {
            cerr << to_string(k) << \
            "-th vortex out of the simulated cell (x-axis)!" << endl;
            exit(-13);
        }
        if ((yv < -0.5*abs(_L2)) || (yv >= 0.5*abs(_L2))) {
            cerr << to_string(k) << \
            "-th vortex out of the simulated cell (y-axis)!" << endl;
            exit(-13);
        }

        // save the k-th vortex position and charge
        _posv[k] = Vector2(xv, yv);
        _vcharge[k] = static_cast<int>(charge_val);
        // compute the indices of the spatial grid's point
        //  closest to the (real) vortex position
        _iXv[k] = static_cast<int>((xv + _L1/2)/_dx);
        _iYv[k] = static_cast<int>((yv + _L2/2)/_dy);
    }
    input.close();

    // compute the (charged) total number of vortices
    _Nv = 0;
    for (unsigned int k = 0; k < _NvAll; k++) {_Nv += _vcharge[k];}
}



void Superfluid::ConstC() {
    /*
    Computes the integration constants appearing in the expression of the
        Quasi-Periodic Boundary Conditions (QPBC)

    The constant c1 [c2] is obtained by comparing the values that the
        superfluid's phase assumes on the two vertical [horizontal] edges of
        the extended fundamental cell, identified by (i=0) and (i=_lenx) [(j=0)
        and (j=_leny)]
    The integration constants must not depend on the specific couple of points
        chosen on these edges: this method computes the values of the cs from
        each of these couples and then saves in the internal variables <_c1>
        and <_c2> their average
    This average is weighted with the superfluid's local density in order to
        avoid some shenanigans due to a possible proximity of the vortices to
        the cell edges
    
    For the purpose of this method, the dynamic evolution of the superfluid is
        computed explicitly within the GPE class (see './PDE.h') all over the
        "extended" spatial grid (this methods does NOT uses the QPBC)
    */

    double ThetaA, ThetaB;
    double dens;
    double cij, c0, c1 = 0., c2 = 0.;
    double w1 = 0., w2 = 0.;
    bool over;

    // compute c1
    for (int j = 0; j < _leny; j++) {
        ThetaA = _Psi[j].Phase();
        ThetaB = _Psi[(_leny+1)*_lenx + j].Phase();
        cij = Mod2pi(ThetaB - ThetaA - M_PI*_Nv*_y[j]/_L2);

        if (j == 0) {
            c0 = cij;
            if (c0 > M_PI) {
                over = 1;
            } else
                over = 0;
        } else if (abs(cij-c0) > M_PI) {
            if (over == 1) {
                cij += 2.*M_PI;
            } else
                cij -= 2.*M_PI;
        }

        dens = Dens(0, j);
        c1 += dens*cij;
        w1 += dens;
    }
    _c1 = Mod2pi(c1/w1);

    // compute c2
    for (int i = 0; i < _lenx; i++) {
        ThetaA = _Psi[(_leny+1)*i].Phase();
        ThetaB = _Psi[(_leny+1)*i + _leny].Phase();
        cij = Mod2pi(ThetaB - ThetaA + M_PI*_Nv*_x[i]/_L1);

        if (i == 0) {
            c0 = cij;
            if (c0 > M_PI) {
                over = 1;
            } else
                over = 0;
        } else if (abs(cij-c0) > M_PI) {
            if (over == 1) {
                cij += 2.*M_PI;
            } else
                cij -= 2.*M_PI;
        }

        dens = Dens(i, 0);
        c2 += dens*cij;
        w2 += dens;
    }
    _c2 = Mod2pi(c2/w2);
}



void Superfluid::UpdateSide() {
    /*
    compute the superfluid's wave function on the extra sides (extended grid)
        by employing QPBC
    */

    for (int j = 0; j < _leny+1; j++) {
        _Psi[(_leny+1)*_lenx + j] = RI_PBC(_Psi, _lenx, j, 1);
    }
    for (int i = 0; i < _lenx; i++) {   // i=_lenx already QPBC-ed
        _Psi[(_leny+1)*i + _leny] = RI_PBC(_Psi, i, _leny, 1);
    }
}



void Superfluid::UpdateVBg(const double time) {
    /*
    update the background potential velocity with a specified function for
      vBg = vBg(t)
    */

    if (_vBgT > 0.) {
        if (_vBgEvol == "sin") {
            if (time < 0.) {
                _vBg = Vector2(0., 0.);
            } else if (time < _vBgT) {
                _vBg = _vBg0 * (1. - cos(M_PI*time / (2.*_vBgT)));
            } else {
                _vBg = _vBg0 * (1. + M_PI/(2.*_vBgT)*(time - _vBgT));
            } 
        } else if (_vBgEvol == "cost"){
            if (time < 0.) {
                _vBg = Vector2(0., 0.);
            } else if (time < _vBgT) {
                _vBg = _vBg0 * (1. - cos(M_PI*time / _vBgT)) / 2 ;
            } else {
                _vBg = _vBg0;
            }
        } else if (_vBgEvol == "limit"){
            if (time < 0.) {
                _vBg = Vector2(0., 0.);
            } else if (time < _vBgT) {
                _vBg = _vBg0 * (1. - cos(M_PI*time / _vBgT)) / 2 ;
            } else if (time < _vBgT * 2){
                _vBg = _vBg0;
            } else if (time < _vBgT * 3) {
                _vBg = _vBg0 * (1 + cos(M_PI*time / _vBgT)) / 2 ;
            } else{
                _vBg = Vector2(0., 0.);
            }
        }
        
        else {
            cerr << "The required vBg evolution is not implemented\nThe implemented values are \"sin\" \"cost\" \"limit\"" << endl;
            exit(-19);
        }

    }
}



// BOUNDARY CONDITIONS

double Superfluid::PBC(const double* f, const int i, const int j) const {
    /*
    Reads the value of a perfectly periodic real function on a specified point
        of the (infinite) spatial grid
    This method assumes that (only) the values that this perfectly periodic
        function assumes on the usual grid are stored in a dynamic array
    
    :arg: <f>   pointer to the dynamic array storing the values of the function
                that has to be evaluated
    :arg: <i>   first index of the (infinite) grid's point where the function
                <f> has to be evaluated
    :arg: <j>   second index of the (infinite) grid's point where the function
                <f> has to be evaluated
    
    Returns <f>(x[<i>],y[<j>]) in true PBC, i.e. the value stored for the
        function <f> related to the correspondent point in the fundamental cell
    */

    int ii, jj, unused1 = 0, unused2 = 0;
    // brings i and j in back to the fundamental cell
    CountPBC(unused1, unused2, i, j, ii, jj, 1);
    // returns f(_x[i],_y[j]) = f(_x[ii],_y[jj])
    return f[_leny*ii+jj];
}


Vector2 Superfluid::RI_PBC(const ComplexFunc& psi, const int i, const int j,
                                                      bool extraSide) const {
    /*
    Reads the value that a quasi-periodic complex function assumes on a
        specified point of the (infinite) spatial grid by employing the
        Quasi-Periodic Boundary Conditions (QPBC)
    Usually these functions are described by a ComplexFunc that stores the
        values they assume on the "extended" fundamental cell; here, also the
        extra edges can be considered as external to the fundamental cell and
        then QPBC-ed
    
    :arg: <psi>         complex function that has to be QPBC-ed
    :arg: <i>           first index of the (infinite) grid's point where the
                        function <psi> has to be evaluated
    :arg: <j>           second index of the (infinite) grid's point where the
                        function <psi> has to be evaluated
    :arg: <extraSide>   = 0 to consider the extra side as part of the unit cell,
                        = 1 to calculate it with QPBC
    
    Returns <psi>(x[<i>], y[<j>]) in QPBC: its value is computed from the one
        assumed by <psi> on the correspondend point of the fundamental cell, by
        employing the QPBC
    */

    int ii, jj, countx = 0, county = 0;
    double A, vprod, real, im;
    // brings i and j in back to the fundamental cell
    CountPBC(countx, county, i, j, ii, jj, extraSide);
    // compute each term of (6.28) and (6.29)
    vprod = M_PI*_Nv*(countx*_y[jj]/_L2 \
                      - county*_x[ii]/_L1 \
                      - countx*county);     // first two terms of (6.30)
    A = vprod + countx*_c1 + county*_c2;
    real = psi[(_leny+1)*ii+jj][0]*cos(A) - psi[(_leny+1)*ii+jj][1]*sin(A);
    im = psi[(_leny+1)*ii+jj][1]*cos(A) + psi[(_leny+1)*ii+jj][0]*sin(A);
// test: impose PBC instead (if Nv=0, difference the is that now we impose ci=0)
//real = psi[(_leny+1)*ii+jj][0];
//im = psi[(_leny+1)*ii+jj][1];
    return Vector2(real, im);
}


void Superfluid::CountPBC(int &countx, int &county, const int i, const int j,
                          int &ii, int &jj, bool extraSide) const {
    /*
    Takes the grid indices of an arbitrary point of the space and finds: the
        indices of the cell where the point is, the grid indices of the
        corresponding point of the fundamental cell
    This function works with the "usual" fundamental cell (the extra edges of
        the extended cell are not considered as part of the fundamental cell)
    
    :arg: <countx>      will store the first cell index
    :arg: <county>      will store the second cell index
    :arg: <i>           first index of the given point of the (infinite) grid
    :arg: <j>           second index of the given point of the (infinite) grid
    :arg: <ii>          will store the first index of the QPBC-ed point
    :arg: <jj>          will store the second index of the QPBC-ed point
    :arg: <extraSide>   = 0 to consider the extra side as part of the unit cell,
                        = 1 to calculate it with QPBC
    */

    int limx = _lenx;
    int limy = _leny;
    if (extraSide) {
        limx--;
        limy--;
    }
    countx = 0;
    county = 0;
    ii = i;
    jj = j;
    while (ii > limx) {
        ii -= _lenx;
        countx++;
    }
    while (ii < 0) {
        ii += _lenx;
        countx--;
    }
    while (jj > limy) {
        jj -= _leny;
        county++;
    }
    while (jj < 0) {
        jj += _leny;
        county--;
    }
}





// SPATIAL DERIVATIVES

Vector2 Superfluid::Derivative(const ComplexFunc& psi, const int order,
                               const int axis, const int i, const int j,
                               const unsigned int nStencil) const {
    /*
    Compute the 1-st or 2-nd order partial derivatives of a quasi-periodic
        complex function in a specified point of the (infinite) spatial grid
    The QPBC implemented in RI_PBC are used to compute the values of the
        considered function outside from the "usual" fundamental cell
    A 3-points or 5-points stencil is used to compute the derivatives
    References:
        1-st order https://en.wikipedia.org/wiki/Compact_stencil#Two_Point_Stencil_Example
                   https://en.wikipedia.org/wiki/Five-point_stencil#1D_first_derivative
        2-nd order https://en.wikipedia.org/wiki/Compact_stencil#Three_Point_Stencil_Example
                   https://en.wikipedia.org/wiki/Five-point_stencil#1D_higher-order_derivatives)
    
    :arg: <psi>       complex function that has to be derivated
    :arg: <order>     order of the derivative (implemented 1-st and 2-nd)
    :arg: <axis>      direction of the partial derivative (0 for x, 1 for y)
    :arg: <i>         first index of the (infinite) grid's point of evaluation
    :arg: <j>         second index of the (infinite) grid's point of evaluation
    :arg: <nStencil>  number of points of the stencil used
    
    Returns the value that the <order>-order partial derivative, in the
        direction x (<axis>=0) or y (<axis>=1), of the complex wave function
        <psi> assumes in (x[<i>],y[<j>])
    */

    // 3-points stencil
    if (nStencil == 3) {
        if ((order == 1) && (axis == 0)) {
            return (RI_PBC(psi, i+1, j) - RI_PBC(psi, i-1, j))/(2.*_dx);
        } else if ((order == 1) && (axis == 1)) {
            return (RI_PBC(psi, i, j+1) - RI_PBC(psi, i, j-1))/(2.*_dy);
        } else if ((order == 2) && (axis == 0)) {
            return (RI_PBC(psi, i+1, j) - RI_PBC(psi, i, j)*2. \
                    + RI_PBC(psi, i-1, j)) /pow(_dx, 2);
        } else if ((order == 2) && (axis == 1)) {
            return (RI_PBC(psi, i, j+1) - RI_PBC(psi, i, j)*2. \
                    + RI_PBC(psi, i, j-1)) /pow(_dy, 2);
        } else {
            cerr << "The required derivative is not implemented" << endl;
            exit(-15);
        }
    // 5-points stencil
    } else if (nStencil == 5) {
        if ((order == 1) && (axis == 0)) {
            return (RI_PBC(psi, i+2, j) - RI_PBC(psi, i+1, j)*8. \
                    + RI_PBC(psi, i-1, j)*8. - RI_PBC(psi, i-2, j)) /(-12.*_dx);
        } else if ((order == 1) && (axis == 1)) {
            return (RI_PBC(psi, i, j+2) - RI_PBC(psi, i, j+1)*8. \
                    + RI_PBC(psi, i, j-1)*8. - RI_PBC(psi, i, j-2)) /(-12.*_dy);
        } else if ((order == 2) && (axis == 0)) {
            return (RI_PBC(psi, i+2, j) - RI_PBC(psi, i+1, j)*16. \
                    + RI_PBC(psi, i, j)*30. - RI_PBC(psi, i-1, j)*16. \
                    + RI_PBC(psi, i-2, j)) /(-12.*pow(_dx, 2));
        } else if ((order == 2) && (axis == 1)) {
            return (RI_PBC(psi, i, j+2) - RI_PBC(psi, i, j+1)*16. \
                    + RI_PBC(psi, i, j)*30. - RI_PBC(psi, i, j-1)*16. \
                    + RI_PBC(psi, i, j-2)) /(-12.*pow(_dy, 2));
        } else {
            cerr << "The required derivative is not implemented" << endl;
            exit(-15);
        }
    } else {
        cerr << "The required stencil is not implemented" << endl;
        exit(-15);
    }
}



Vector2 Superfluid::Lap9(const ComplexFunc& psi,
                         const int i, const int j) const {
    /*
    Compute the Laplacian of a quasi-periodic complex function in a specified
        point of the (infinite) spatial grid
    The QPBC implemented in RI_PBC are used to compute the values of the
        considered function outside from the "usual" fundamental cell
    A 9-points stencil is used to compute the Laplacian (reference:
        https://math.stackexchange.com/questions/2916234/how-to-obtain-the-9-point-laplacian-formula)
    
    :arg: <psi>     complex function that has to be derivated
    :arg: <i>       first index of the (infinite) grid's point of evaluation
    :arg: <j>       second index of the (infinite) grid's point of evaluation
    
    Returns the Laplacian of the wave function <psi> in (x[<i>],y[<j>])
    */

    if (_dx != _dy) {
        cerr << "Lap9 method works only if _dx==_dy!" << endl;
        exit(-16);
    }
    Vector2 u_00   = RI_PBC(psi, i,   j);
    Vector2 u_sx   = RI_PBC(psi, i-1, j);
    Vector2 u_dx   = RI_PBC(psi, i+1, j);
    Vector2 u_dw   = RI_PBC(psi, i,   j-1);
    Vector2 u_up   = RI_PBC(psi, i,   j+1);
    Vector2 u_sxdw = RI_PBC(psi, i-1, j-1);
    Vector2 u_sxup = RI_PBC(psi, i-1, j+1);
    Vector2 u_dxdw = RI_PBC(psi, i+1, j-1);
    Vector2 u_dxup = RI_PBC(psi, i+1, j+1);
    return ((u_sx + u_dx + u_dw + u_up)*4. \
            + u_sxdw + u_sxup + u_dxdw + u_dxup - u_00*20.) /(6.*pow(_dx, 2));
}





// SUPERFLUID QUANTITIES

Vector2 Superfluid::Current(const int i, const int j) const {
    /*
    Returns the superfluid internal (only vortex-generated, no bg flow) current
        density in a specific point
    
    :arg: <i>       first index of the (infinite) grid's point where the current
                    density has to be evaluated
    :arg: <j>       second index of the (infinite) grid's point where the current
                    density has to be evaluated
    */

    Vector2 dxPsi, dyPsi, psi;
    dxPsi = Derivative(_Psi, 1, 0, i, j);
    dyPsi = Derivative(_Psi, 1, 1, i, j);
    psi = RI_PBC(_Psi, i, j);
    return Vector2(psi.Cross(dxPsi), psi.Cross(dyPsi));
}


SupCons Superfluid::Conserved() const {
    /*
    Computes the superfluid quantities that are conserved during
        the time evolution of an isolated superfluid
    The total momentum is in natural units; the background velocities are
        NOT considered here
    The velocities circulation along the cell's boundaries is evaluated with
        its definition integral (sum of v dot dl along the extended cell's
        boundaries)
    
    Returns the ensamble of instant values of the quantities just listed as an
        instance of the SupCons class
    */

    SupCons cons;   // all the conserved quantities are initialized to 0.
    double dens;
    Vector2 Jcur, psi, vecRe, vecIm;

    // conserved surface densities
    for (int i = 0; i < _lenx; i++) {
      for (int j = 0; j < _leny; j++) {
        // superfluid density
        dens = Dens(i, j);
        cons.totN += dens;
        // particles linear momemtum
        psi = RI_PBC(_Psi, i, j);
        for (int ax = 0; ax < 2; ax++) {
          vecRe.SetComp(ax, Derivative(_Psi, 1, ax, i, j)[0]);
          vecIm.SetComp(ax, Derivative(_Psi, 1, ax, i, j)[1]);
        }
        cons.partE += 0.5*(pow(vecRe.Mod(), 2) + pow(vecIm.Mod(), 2)) \
                     + (_pinEval[_leny*i+j]/_mu-1.)*dens \
                     + 0.5*pow(dens, 2)/_nmean;
        // superfluid current density
        Jcur = Current(i, j);
        cons.partP += Jcur;
      }
    }
    cons.totN *= _dx*_dy;
    cons.partE *= _mu*_dx*_dy/cons.totN;
    cons.partP *= _dx*_dy/cons.totN;

    // conserved boundary's linear densities
    Vector2 vel;
    double apcirc = 0.;
    for (int i = 0; i < _lenx; i++) {
        vel = Velocity(i, 0);        // lower side contribution
        cons.circ += vel[0];
        vel = Velocity(i+1, _leny);  // upper side contribution
        cons.circ -= vel[0];
    }
    cons.circ *= _dx;
    for (int j = 0; j < _leny; j++) {
        vel = Velocity(_lenx, j);    // right side contribution
        apcirc += vel[1];
        vel = Velocity(0, j+1);      // left side contribution
        apcirc -= vel[1];
    }
    cons.circ += apcirc*_dy;
    cons.circ /= 2.*M_PI;           // units of 2pi

    return cons;
}


double Superfluid::NStdDev() const {
    /*
    Computes the total standard deviation of the density field (useful to
        track vortex formation)
    It is necessary to update the conserved quantities before using this
        function
    */

    double mean_n, delta_n, stdDev=0.;

    mean_n = _cons.totN/(_dx*_dy*_lenx*_leny);
    for (int i = 0; i < _lenx; i++) {
      for (int j = 0; j < _leny; j++) {
        delta_n = Dens(i, j) - mean_n;
        stdDev += pow(delta_n, 2);
      }
    }
    return sqrt(stdDev);
}


void Superfluid::VortexPos(const unsigned int k) {
    /*
    Updates the grid indices of the position of a specified vortex, which is
        found by looking, in a rectangle centered on its immediately previous
        position, for the point of the grid where the superfluid is locally
        maximally depleted (local minimum of the density)
    The grid indices, saved in the internal variables _iXv and _iYv, are not
        QPBC-ed: this allows to trace the "true" vortex motion and to correctly
        estimate the mean vortex velocity
    
    :arg: <k>       this method searches for the position of the <k>-th vortex
    */

    int newiXv, newiYv;
    double minDens, newDens;
    // define the semi-sides of the research rectangle
    int rangex = static_cast<int>(_lenx/10.);
    int rangey = static_cast<int>(_leny/10.);

    minDens = Dens(_iXv[k], _iYv[k]);
    newiXv = _iXv[k];
    newiYv = _iYv[k];
    for (int i = _iXv[k]-rangex; i < _iXv[k]+rangex; i++) {
        for (int j = _iYv[k]-rangey; j < _iYv[k]+rangey; j++) {
            newDens = Dens(i, j);
            if (newDens < minDens) {
                minDens = newDens;
                newiXv = i;
                newiYv = j;
            }
        }
    }
    _iXv[k] = newiXv;
    _iYv[k] = newiYv;
}





// ACCESS TO PRIVATE VARIABLES

double* Superfluid::Getn() const {
    /*
    Returns a pointer to a dynamic array of double of length _lenx*_leny that
        gives the superfluid density field all over the "usual" spatial grid
    */

    double* n = new double[_lenx*_leny];
    for (int i = 0; i < _lenx; i++) {
        for (int j = 0; j < _leny; j++) {
            n[_leny*i+j] = Dens(i, j);
        }
    }
    return n;
}


double* Superfluid::GetPhase() const {
    /*
    Returns a pointer to a dynamic array of double of length _lenx*_leny that
        gives the wave function phase field all over the "usual" spatial grid
    */

    double* phase = new double[_lenx*_leny];
    for (int i = 0; i < _lenx; i++) {
        for (int j = 0; j < _leny; j++) {
            phase[_leny*i+j] = Phase(i, j);
        }
    }
    return phase;
}


Vector2 Superfluid::GetPosV(const unsigned int k) const {
    /*
    Reads the non-PBC-ed position of a specified vortex from its (infinite)
        grid's indices 
    Avoiding to use the QPBC in the estimate of the vortex position allows to
        trace the actual vortex motion and to correctly compute the mean vortex
        velocity
    
    :arg: <k>   index of the vortex considered
    
    Returns the position of the <k>-th vortex
    */

    int countx, county, ii, jj;
    CountPBC(countx, county, _iXv[k], _iYv[k], ii, jj, 1);
    return Vector2(_x[ii] + countx*_L1, _y[jj] + county*_L2);
}





// OUTPUT FILES

void Superfluid::SetOutCell(const double fracCell, const unsigned int outGrid) {
    /*
    Computes the first and last indices of the smaller, central part of the
        fundamental cell that is actually printed in output
    Saves the input value for the internal variable <_outGrid> 
    
    :arg: <fracCell>    will be printed in outputs only the central region of
                        the central cell, large as (<fracCell>*<_lenx>,
                        <fracCell>*<_leny>)
    :arg: <outGrid>     will be printed in outputs only 1 point every <outGrid>
                        points of the spatial grid
    */

    _mmx[0] = static_cast<int>(0.5*_lenx*(1.-fracCell));
    _mmx[1] = static_cast<int>(0.5*_lenx*(1.+fracCell));
    _mmy[0] = static_cast<int>(0.5*_leny*(1.-fracCell));
    _mmy[1] = static_cast<int>(0.5*_leny*(1.+fracCell));
    _outGrid = outGrid;
}


void Superfluid::OutputPhase(const string fname, double time) const {
    /*
    Prints in a file the coordinates of the points of the "usual" space grid
        and the superfluid density, phase and total (vortex-generated + bg
        flow) velocity fields in natural units
    
    :arg: <fname>   the output will be printed in "out/<fname>"
    */

    ofstream output;
    output.open("out/"+fname);
    output << left << setw(15) << "#time" << setw(15) << "_c1" \
                        << setw(15) << "_c2" << endl;
    output << left << setw(15) << time << setw(15) << _c1 \
                        << setw(15) << _c2 << endl;
    
    output << left << setw(15) << "#x [csi]" << setw(15) << "y [csi]" \
                        << setw(18) << "n(x,y) [csi^-2]" \
                        << setw(15) << "phase(x,y)" \
                        << setw(18) << "vx(x,y) [csi/t0]" \
                        << setw(18) << "vy(x,y) [csi/t0]" \
                        << setw(15) << "Re(x,y)" \
                        << setw(15) << "Im(x,y)" \
                        << endl;
    output << setprecision(5);
    Vector2 vel;
    // print variables all over the smaller (output) space grid
    for (int i = _mmx[0]; i < _mmx[1]; i += _outGrid) {
        for (int j = _mmy[0]; j < _mmy[1]; j += _outGrid) {
            vel = Velocity(i, j);
            output << left << setw(15) << _x[i] << setw(15) << _y[j] \
                        << setw(18) << Dens(i, j) << setw(15) << Phase(i, j) \
                        << setw(18) << vel[0] << setw(18) << vel[1] \
                        << setw(15) << Real(i,j) << setw(15) << Img(i,j) \
                        << endl;
        }
    }
    output.close();
}


void Superfluid::OutputPin(const string fname) const {
    /*
    Prints in a file the coordinates of the points of the "usual" space grid
        and the pinning landscape in natural units
    
    :arg: <fname>   the output will be printed in "out/<fname>"
    */

    ofstream output;
    output.open("out/"+fname);
    output << left << setw(15) << "x [csi]" << setw(15) << "y [csi]" \
                   << setw(15) << "V(x,y) [mu]" << endl;
    output << setprecision(3);
    // print variables all over the "usual" space grid
    for (int i = _mmx[0]; i < _mmx[1]; i += _outGrid) {
        for (int j = _mmy[0]; j < _mmy[1]; j += _outGrid) {
            output << left << setw(15) << _x[i] << setw(15) << _y[j] \
                           << setw(15) << _pinEval[_leny*i+j] << endl;
        }
    }
    output.close();
}


void Superfluid::PrintAllGrid(const int n, const int m,
                              const string fname) const {
    /*
    Prints in a file the coordinates of the points of the spatial grid of a
        table of "usual" unit cells centered in the fundamental one and the
        superfluid density, phase, wave function's real and imaginary parts
        and total (vortex-generated + bg flow) velocity field in natural units
    The QPBC are employed to read the values of the superfluid's wave function
        outside of the fundamental cell
    
    :arg: <n>       number of rows of the table of unit cells, must be odd
    :arg: <m>       number of columns of the table of unit cells, must be odd
    :arg: <fname>   the output will be printed in "out/<fname>"
    */

    ofstream output;
    double *x = NULL;
    double *y = NULL;
    Vector2 psi, vel;
    // compute the indices of the lower left angle of the table of unit cells
    int ixmin = -_lenx*(m-1)/2;
    int jymin = -_leny*(n-1)/2;

    LargeGrid(n, m , x, y);  // generates the table of unit cells' spatial grid
    output.open("out/"+fname);
    output << left << setw(15) << "x [csi]" << setw(15) << "y [csi]" \
                        << setw(15) << "n(x,y) [csi^-2]" \
                        << setw(15) << "phase(x,y)" \
                        << setw(15) << "R [csi^-1]" \
                        << setw(15) << "I [csi^-1]" \
                        << setw(15) << "vx(x,y) [csi/t0]" \
                        << setw(15) << "vy(x,y)  [csi/t0]" \
                        << endl;
    output << setprecision(5);
    // print variables all over the table of unit cells' space grid
    for (int i = 0; i < _lenx*static_cast<int>(m); i += _outGrid) {
        for (int j = 0; j < _leny*static_cast<int>(n); j += _outGrid) {
            vel = Velocity(i+ixmin, j+jymin);
            psi = RI_PBC(_Psi, i+ixmin, j+jymin);
            output<< left << setw(15) << x[i] << setw(15) << y[j] \
                    << setw(15) << Dens(i+ixmin, j+jymin) \
                    << setw(15) << Phase(i+ixmin, j+jymin) \
                    << setw(15) << psi[0] << setw(15) << psi[1] \
                    << setw(15) << vel[0] << setw(15) << vel[1] \
                    << endl;
        }
    }
    output.close();
    delete[] x;
    delete[] y;
}


void Superfluid::LargeGrid(const int n, const int m,
                           double* &x, double* &y) const {
    /*
    Extends the "usual" spatial grid of the fundamental cell to that one of a
        table of unit cells centered on the fundamental cell
    
    :arg: <n>   number of rows of the table of unit cells, must be odd
    :arg: <m>   number of columns of the table of unit cells, must be odd
    :arg: <x>   will store the x-coordinates of the points of the spatial grid
                of the table of unit cells
    :arg: <y>   will store the y-coordinates of the points of the spatial grid
                of the table of unit cells
    */

    if ((n%2 != 1) || (m%2 != 1)) {  // for simplicity, n m assumed to be odd
        cerr << endl << "n and m must be odd!" << endl << endl;
        exit(-17);
    }

    x = new double[_lenx*m];
    y = new double[_leny*n];
    for (int nx = 0; nx < m; nx++) {
        for (int i = 0; i < _lenx; i++) {
            x[_lenx*nx+i] = _x[i] + (nx-(m-1)/2)*_L1;
        }
    }
    for (int ny = 0; ny < n; ny++) {
        for (int j = 0; j < _leny; j++) {
            y[_leny*ny+j] = _y[j] + (ny-(n-1)/2)*_L2;
        }
    }
}



void Superfluid::OutCs(const string fname) const {
    /*
    Computes the integration constants appearing in the expression of the
        Quasi-Periodic Boundary Conditions (QPBC) all over the fundamental
        cell's edges and prints the results in a file
    The constant c1 [c2] is obtained by comparing the values that the
        superfluid's phase assumes on the two vertical [horizontal] edges of
        the extended fundamental cell, identified by (i=0) and (i=_lenx) [(j=0)
        and (j=_leny)]
    This method is useful for debugging, as it allows to control that those
        quantities are actually constants over the cell's edges
    
    :arg: <fname>   the output will be printed in "out/<fname>"
    */

    ofstream output;
    double ThetaA, ThetaB;
    double c1, c2;

    output.open("out/"+fname);
    output << setprecision(8);
    for (int j = 0; j < _leny; j += _outGrid) {
        // compute c1
        ThetaA = Phase(0, j);
        ThetaB = _Psi[(_leny+1)*_lenx + j].Phase();
        c1 = Mod2pi(ThetaB - ThetaA - M_PI*_Nv*_y[j]/_L2);
        // print c1 mod2pi around _c1
        if (abs(c1-_c1) > M_PI) {c1 = 2.*M_PI - c1;}
        // print output
        output << left << setw(15) << c1;
    }
    output << endl;
    for (int i = 0; i < _lenx; i += _outGrid) {
        // compute c2
        ThetaA = Phase(i, 0);
        ThetaB = _Psi[(_leny+1)*i + _leny].Phase();
        c2 = Mod2pi(ThetaB - ThetaA + M_PI*_Nv*_x[i]/_L1);
        // print c2 mod2pi around _c2
        if (abs(c2-_c2) > M_PI) {c2 = 2.*M_PI - c2;}
        // print output
        output << left << setw(15) << c2;
    }
    output.close();
}



void Superfluid::OutputParticles(const string fname) const {
    /*
    Prints in a file the position, velocity and acceleration of the active
        particles and the force acting on them in natural units
    
    :arg: <fname>   the output will be printed in "out/<fname>"
    */

    ofstream output;
    output.open("out/"+fname);
    output<< left << setw(15) << "x_j [csi]" \
                  << setw(15) << "y_j [csi]" \
                  << setw(15) << "vx_j(x,y) [csi/t0]" \
                  << setw(15) << "vy_j(x,y) [csi/t0]" \
                  << setw(15) << "ax_j(x,y) [csi/t0^2]" \
                  << setw(15) << "ay_j(x,y) [csi/t0^2]" \
                  << setw(15) << "Fx_j(x,y) [MeV/csi]" \
                  << setw(15) << "Fy_j(x,y) [MeV/csi]" << endl;
    output << setprecision(5);
    for (unsigned int jpart = 0; jpart < _npart; jpart++) {
        output << left << setw(15) << _xpart[jpart][0] \
                       << setw(15) << _xpart[jpart][1] \
                       << setw(15) << _vpart[jpart][0] \
                       << setw(15) << _vpart[jpart][1] \
                       << setw(15) << _apart[jpart][0] \
                       << setw(15) << _apart[jpart][1] \
                       << setw(15) << _fpart[jpart][0] \
                       << setw(15) << _fpart[jpart][1] << endl;
    }
    output.close();
}

