#include "./PDE.h"


// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// PDE (Partial Differential Equation)
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


// CONSTRUCTORS AND DESTRUCTOR

// Default constructor: initializes all the numerical variables to 0 or 0.
//  (except for _tSkip) and all the pointers to NULL
PDE::PDE(): _XY(), _params(NULL), _npar(0), _llenx(0), _lleny(0), \
            _d_t(0), _tSkip(1), _tburn(0.), _ntburn(0), _gamma(0.), \
            _ifGammaDyn(false) {}


PDE::PDE(const int lenx, const int leny, const double dt): PDE() {
    /*
    Construct a PDE instance with specified spatial grid's dimensions and
        solving algorithm's time step
    The PDE's solution is a discrete complex field defined all over the
        "extended" fundamental cell (see also './Superfluid.*')
    
    :arg: <lenx>    number of grid's points along the horizontal side
    :arg: <leny>    number of grid's points along the vertical side
    :arg: <dt>      time resolution of the PDE-solving algorithm
    */

    _XY.SetLen((lenx+1)*(leny+1));
    _llenx = lenx;
    _lleny = leny;
    _d_t = dt;
}


PDE::PDE(const int lenx, const int leny, const double dt,
         const unsigned int npar): PDE(lenx, leny, dt) {
    /*
    Construct a PDE instance with specified spatial grid's dimensions, solving
        algorithm's time step and number of PDE's parameters
        
    :arg: <lenx>    number of grid's points along the horizontal side
    :arg: <leny>    number of grid's points along the vertical side
    :arg: <dt>      time resolution of the PDE-solving algorithm
    :arg: <npar>    number of PDE's parameters
    */

    _params = new double[npar];
    _npar = npar;
}





// PDE-SOLVING METHODS

void PDE::Solver(const double tend, const unsigned int delta,
                 const unsigned int method, const string identifier) {
    /*
    Solves the complex Partial Differential Equation dXY/dt=der by employing a
        specific numerical method (4-th order Runge-Kutta, 3-rd order
        Adam-Bashforth or Euler)
    Some operations on the time-dependent solution (X + iY) is (e.g. printing
        in outputs) could be implemented in the derivated' EvStuff
    The PDE-solving algorithm is runned for a specified number of time steps
    Some of the single step of the specific solving methods are implemented
        separately (see the methods RKstep, ABstep and ABinit)
    
    The PDE is solved independently for all the points of the spatial grid of
        the "extended" fundamental cell (see also './Superfluid.*'), i.e. the
        solution field is definied all over the fundamental cell and on all its
        edges
    
    Inside the implementation of the Runge-Kutta and Adam-Bashforth methods,
        the slopes coefficients k (reference [2]) and ff (reference [3]) are
        defined as an array of complex functions of length equal to the order
        of the algorithm
    
    :arg: <tend>        final time of the simulation
    :arg: <delta>       the output files will be updated every <delta> steps
    :arg: <method>      = 0 for the 4-th order Runge-Kutta explicit method
                        = 1 for the 3-rd order Adam-Bashforth explicit method
                        = 2 for the Euler explicit method
    :arg: <identifier>  some of the outputs produced by EvStuff will have added
                        in their name this string as '*<identifier>.out'
    */

    unsigned int len = (_llenx+1)*(_lleny+1);
    unsigned int order = -1;    // order of the method used
    ComplexFunc der(len);       // complex derivative of the solution XY

    // total number of time steps needed
    unsigned int ntend = int((tend + _tburn) / _d_t);

    // 4-th order Runge-Kutta method
    if (method == 0) {
        order = 4;
        ComplexFunc XY(len);    // support complex function
        ComplexFunc *k = new ComplexFunc[order];
        for (unsigned int ord = 0; ord < order; ord++) {k[ord].SetLen(len);}
        cout << \
        "The 4-th order Runge-Kutta method is used to solve the time GPE" \
        << endl << endl;

        for (unsigned int it = 0; it < ntend; it++) {
            RKstep(XY, der, k);
            EvStuff(it, delta, identifier);
        }

        delete[] k;

    // 3-rd order Adam-Bashforth method
    } else if (method == 1) {
        order = 3;
        ComplexFunc *ff = new ComplexFunc[order-1];
        for (unsigned int ord = 0; ord < order-1; ord++) {ff[ord].SetLen(len);}
        cout << \
        "The 3-th order Adams-Bashforth method is used to solve the time GPE" \
        << endl << endl;

        ABinit(der, ff, order, delta, identifier);
        for (unsigned int it = 2; it < ntend; it++) {
            ABstep(der, ff);
            EvStuff(it, delta, identifier);
        }

        delete[] ff;

    // basic Euler method
    } else if (method == 2) {
        order = 1;
        cout << \
        "The basic Euler method is used to solve the time GPE" \
        << endl << endl;

        for (unsigned int it = 0; it < ntend; it++) {
            Equation(der, _XY);
            _XY += der*_d_t;
            EvStuff(it, delta, identifier);
        }

    } else {
        cerr << "The selected PDE-solving method is not implemented!" << endl;
        exit(-50);
    }
}


void PDE::RKstep(ComplexFunc& XY, ComplexFunc& der, ComplexFunc* k) {
    /*
    Method that implements a single time step of the 4-th order Runge-Kutta
        algorithm for the solution of the the complex PDE dXY/dt=der
        (reference [2])
    The actual derivative of the solution is given and computed within the
        derived classes' method Equation
    
    :arg: <XY>    support complex function, passed as an argument to avoid
                  to dynamically allocate and deallocate it at each time step
    :arg: <der>   complex derivative of <XY>, passed as an argument to avoid
                  to dynamically allocate and deallocate it at each time step
    :arg: <k>     RK internal slope functions, passed as an argument to avoid
                  to dynamically allocate and deallocate it at each time step
    */

    unsigned int order = 4;  // order of the method used
    XY = _XY;
    for (unsigned int ord = 0; ord < order; ord++) {
        Equation(der, XY);
        k[ord] = der;
        if (ord == 2) {
            XY = _XY + der*_d_t;
        } else if (ord != 3) {
            XY = _XY + der*0.5*_d_t;
        }
    }
    XY = (k[0] + k[1]*2. + k[2]*2. + k[3])/6.;
    _XY += XY*_d_t;
}


void PDE::ABstep(ComplexFunc& der, ComplexFunc* ff) {
    /*
    Method that implements a single time step of the 3-rd order Adam-bashforth
        algorithm for the solution of the the complex PDE dXY/dt=der
        (reference [3])
    The actual derivative of the solution is given and computed within the
        derived classes' method Equation
    
    :arg: <der>   complex derivative of <XY>, passed as an argument to avoid
                  to dynamically allocate and deallocate it at each time step
    :arg: <ff>    AB internal slope functions, passed as an argument to avoid
                  to dynamically allocate and deallocate it at each time step
    */

    Equation(der, _XY);
    _XY += (der*23./12. - ff[1]*4./3. + ff[0]*5./12.)* _d_t;
    ff[0] = ff[1];
    ff[1] = der;
}

void PDE::ABinit(ComplexFunc& der, ComplexFunc* ff, const unsigned int order,
                 const unsigned int delta, const string identifier) {
    /*
    Method that initialize the Adam-bashforth of the specified order
        (reference [3])
    The actual derivative of the solution is given and computed within the
        derived classes' method Equation
    
    :arg: <der>         complex derivative of <XY>, passed as an argument to avoid
                        to dynamically allocate and deallocate it at each time step
    :arg: <ff>          AB internal slope functions, passed as an argument to avoid
                        to dynamically allocate and deallocate it at each time step
    :arg: <order>       order of the method used
    :arg: <delta>       the output files will be updated every <delta> steps
    :arg: <identifier>  some of the outputs produced by EvStuff will have added in
                        their name this string as '*<identifier>.out'
    */

    for (unsigned int nstp = 0; nstp < order-1; nstp++) {
        Equation(der, _XY);
        ff[nstp] = der;
        if (nstp == 0) {
            _XY += der*_d_t;
        } else {
            _XY += (der*1.5 - ff[0]*0.5)* _d_t;
        }
        EvStuff(nstp, delta, identifier);
    }
}





// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// GPE (Gross-Pitaevskii Equation)
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


// CONSTRUCTORS, DESTRUCTOR AND RELATED METHODS

GPE::GPE(): PDE(), Superfluid(), _Omega0(0.), \
                                 _KelX(NULL), _KelY(NULL), _ifKelvin(0) {
    /*
    Default constructor: initializes all the numerical variables to 0 or 0.,
        all the pointers to NULL and identifies the Superfluid's wave function
        with the solution of the GPE
    */

    _XY.SameObj(_Psi);
}


void GPE::AllConstr() {
    /*
    Method that implements the operations that are common to all the GPE's
        constructors
    It identifies the correspondent variables that are present both in PDE and
        in Superfluid classes
    This method assumes that the Kelvin theorem will not be checked; if it has
        to be checked, the method SetKelvin must be called in main after the
        GPE constructor
    */
    ifstream input;
    string buffer;

    _XY.SameObj(_Psi);
    _llenx = _lenx;
    _lleny = _leny;
    _d_t = _dt;
    input.open(PARAM_file);
    for (unsigned int ibuff=0; ibuff<36; ibuff++) {ReadPar(input, buffer);}
    ReadPar(input, _tburn);
    _ntburn = int(_tburn / _d_t);
    ReadPar(input, _gamma);
    ReadPar(input, _ifGammaDyn);
    input.close();
    _Omega0 = M_PI*_Nv/(_L1*_L2);
    _ifKelvin = 0;
    _KelX = NULL;
    _KelY = NULL;
}


// destructor
GPE::~GPE() {
    delete[] _KelX;
    delete[] _KelY;
    _XY.FreePointer();  // memory freed by Superfluid's destructor
}


void GPE::SetKelvin() {
    /*
    Initializes the variables needed for the Kelvin theorem's test by
        dynamically allocating the arrays which will store the coordinates of
        the points of a closed circuit inside the superfluid
    The number of points is arbitrarily taken equal to the number of points of
        the horizontal side of the spatial grid
    For simplicity, the initial closed curve is arbitrarily chosen as an
        horizontal line, and the coordinates of its point are printed, in n.u.,
        in 'out/Kelvin_0.out'
    Anyway, it must be a closed curve to be able to check the Kelvin theorem
    */

    ofstream circuit;
    _ifKelvin = 1;
    _KelX = new double[_lenx];  // arbitrary length
    _KelY = new double[_lenx];

    circuit.open("out/Kelvin_0.out");
    for (int i = 0; i < _lenx; i++) {
        _KelX[i] = i;
        _KelY[i] = _leny/2.5;
        circuit << left << setw(12) << _x[static_cast<int>(round(_KelX[i]))] \
                        << setw(12) << _y[static_cast<int>(round(_KelY[i]))] \
                        << endl;
    }
    circuit.close();
}





// GPE AND GPE-SOLVING METHODS

void GPE::Equation(ComplexFunc& der, const ComplexFunc& XY) const {
    /*
    Implements the Gross-Pitaevskii Equation by computing the time derivative
        of the superfluid wave function's real and imaginary parts
    Dissipation is included by solving i(1 + i gamma) dPsi = RHS instead
        of i dPsi = RHS
    Note that the GPE solution (equivalently the superfluid's wave function) is
        quasi-periodic and then it must be evolved on the "extended"
        fundamental cell in order to keep trace of the evolution of the
        integration constants c_i of the QPBC
    
    :arg: <der>  complex time-derivative of the superfluid wave function, passed
                 as an argument to avoid to dynamically allocate and deallocate
                 the correspondent Vector2 array at each time step
    :arg: <XY>   complex function which derivative is given by the GPE
    */

    double Op, dRe, dIm;
    Vector2 dxXY, dyXY, lap;
    int gpos;

    for (int i = 0; i < _lenx+1; i++) {
        for (int j = 0; j < _leny+1; j++) {
            gpos = (_leny+1)*i + j;
            Op = pow(XY[gpos].Mod(), 2)/_nmean \
                + PBC(_pinEval, i, j)/_mu - 1.;

            dxXY = Derivative(XY, 1, 0, i, j);
            dyXY = Derivative(XY, 1, 1, i, j);
            if (_dx == _dy) {
                lap = Lap9(XY, i, j);
            } else {
                lap = Derivative(XY, 2, 0, i, j) \
                     + Derivative(XY, 2, 1, i, j);
            }

            dRe = ( Op * (XY[gpos][1] - _gamma * XY[gpos][0]) \
                    - 0.5 * (lap[1] - _gamma * lap[0]) \
                    + _Omega0 * ( _x[i]*(dyXY[0] + _gamma*dyXY[1]) \
                                - _y[j]*(dxXY[0] + _gamma*dxXY[1]) ) ) \
                    / (1. + pow(_gamma, 2));
            dIm = ( - Op * (XY[gpos][0] + _gamma * XY[gpos][1]) \
                    + 0.5 * (lap[0] + _gamma * lap[1]) \
                    + _Omega0 * ( _x[i]*(dyXY[1] - _gamma*dyXY[0]) \
                                - _y[j]*(dxXY[1] - _gamma*dxXY[0]) ) ) \
                    / (1. + pow(_gamma, 2));

            der.SetValue(gpos, dRe, dIm);
        }
    }

}



void GPE::EvStuff(const unsigned int it, const unsigned int delta,
                                        const string identifier) {
    /*
    This method implements some generic operations that have to be repeated
        each time step
    
    It computes the complete time-evolution of the QPBC's integration constants
        c_i and, if required, of the shape of the circuit needed for the Kelvin
        theorem's check
    It constantly update the grid indices of the positions of the vortices,
        which are found with VortexPos
    
    Every <_tSkip> time steps, it computes and prints out in 'out/cons.out' the
        superfluid's conserved quantities
    Every chosen number of time steps, it prints on terminal the number of the
        steps already computed and prints out the istantaneous results for the
        c_i and, if required, for the Kelvin circuit respectively in
        'out/Cs_#.out' and 'out/Kelvin_#.out'
    It also prints out the instantaneous superfluid's configuration (phase,
        density and velocity) in 'out/phase_#.out' and update the istantaneous
        positions of the vortices, saved in 'out/vortPos.out', where is also
        saved the instantaneous velocity of the superfluid evaluated on the
        PBC-ed positions of the vortices
    
    :arg: <it>           number of the current time step
    :arg: <delta>       this method will update the intermediate outputs files
                        every <delta> steps and will print current the step on
                        screen every <delta>/5 steps
    :arg: <identifier>  some of the outputs will have added in their name this
                        string as '*<identifier>.out'
    */

    ofstream nout, circuit, outVort;
    string idOutName;
    int countx, county, ii, jj;
    int newt = it-_ntburn+1;
    int newtOverDelta;
    // time from the beginning of the simulation, n.u.
    double time = newt*_dt;

    // update background velocity
    UpdateVBg(time);

    // remove dissipation (only used to thermalize the initial condition)
    if ( (it > _ntburn) && (_gamma != 0.) && (!_ifGammaDyn) ) {_gamma = 0.;}

    // update the pinning landscape after a burn-in phase (let the sound waves
    // diffuse before starting with the actual simulation)
    if (it > _ntburn) {
        // account for active particles
        if (_pinning->GetIfActivePart())
            PinningEvolve();
        // account for background flow (fixed shift of the potential)
        _pinning->BackgroundFlow(_vBg * _dt);
        PinningEval();
    }

    // compute the new QPBC's integration constants
    ConstC();
    // intermediate outputs (only if the whole fundamental cell is printed)
    //if (((it+1)%delta == 0) && (_mmx[0] == 0) && (_mmy[0] == 0)) {
    //    OutCs("Cs_"+to_string((it+1)/delta)+".out");
    //}
    // update the values on the extra side with QPBC
    UpdateSide();

    // compute the superfluid's contants of motion
    if (it % _tSkip == 0) {_cons = Conserved();}
    // compute the new vortex indices
    for (unsigned int k = 0; k < _NvAll; k++) {VortexPos(k);}
    // update the Kelvin circuit
    if (_ifKelvin) {KelvinCheck();}

    // initialize the output files
    if (it == 0) {
        // conserved quantities
        nout.open("out/cons"+identifier+".out");
        nout.close();
        nout.open("out/cons"+identifier+".out", fstream::app);
        nout << left << setw(15) << "t [n.u.]" \
                     << setw(25) << "N" \
                     << setw(15) << "E [MeV]" \
                     << setw(15) << "Px [n.u.]" \
                     << setw(15) << "Py [n.u.]" \
                     << setw(15) << "circulation [2pi]" << endl;   // legend
        nout.close();
        
        // background flow velocity and total std dev of density
        if (_vBg0.Mod() > 0) {
            nout.open("out/vel_and_nStdDev"+identifier+".out");
            nout.close();
            nout.open("out/vel_and_nStdDev"+identifier+".out", fstream::app);
            nout << left << setw(15) << "# t [n.u]" \
                         << setw(15) << "vBgx [csi/t0]" \
                         << setw(15) << "vBgy [csi/t0]" \
                         << setw(15) << "n stdDev [n.u.]" << endl;   // legend
            nout.close();
        }

        // position of vortices
        if (_NvAll > 0) {
            outVort.open("out/vortPos"+identifier+".out");
            outVort.close();
            outVort.open("out/vortPos"+identifier+".out", fstream::app);
            outVort << left << setw(15) << "t [n.u.]" << flush;
            for (unsigned int k = 0; k < _NvAll; k++) {
                outVort << left << setw(3) << "xv_" << setw(6) << to_string(k) \
                                << setw(3) << "yv_" << setw(6) << to_string(k) \
                                << flush;
            }
            outVort << endl;
            outVort.close();
        }
    }

    if (it % _tSkip == 0) {
        // output of the superfluid's contants of motion
        nout.open("out/cons"+identifier+".out", fstream::app);
        nout << left << setw(15) << time \
                     << setw(25) << setprecision(15) << scientific << _cons.totN \
                     << setw(15) << setprecision(6) << _cons.partE \
                     << setw(15) << _cons.partP[0] \
                     << setw(15) << _cons.partP[1] \
                     << setw(15) << _cons.circ << endl;
        nout.close();
    }

    // intermediate outputs
    if (newt >= 0) {idOutName = "";}
    else {idOutName = "burn_";}
    if (newt%delta == 0) {
        newtOverDelta = newt/int(delta);
        // print the step number
        cout << endl << "    t = " << time \
                     << "\t_R[0] = " << _Psi[0][0] << endl;

        // phase, density and c_i
        OutputPhase("phase_"+idOutName+to_string(newtOverDelta)+".out", time);
        //pinning
        if ((_pinning->GetIfActivePart()) or (_vBg0.Mod()!=0.)) {
            OutputPin("pin_"+idOutName+to_string(newtOverDelta)+".out");
        }
        if (_pinning->GetIfActivePart()) {
            OutputParticles( \
                    "particles_"+idOutName+to_string(newtOverDelta)+".out");
        }

        // vortex positions
        if (_NvAll > 0) {
            outVort.open(\
                    "out/vortPos"+identifier+idOutName+".out", fstream::app);
            outVort << left << setw(15) << time << flush;
            for (unsigned int k = 0; k < _NvAll; k++) {
                CountPBC(countx, county, _iXv[k], _iYv[k], ii, jj, 1);
                outVort << left << setw(9) << _x[ii] << setw(9) << _y[jj] << flush;
            }
            outVort << endl;
            outVort.close();
        }

        // background flow velocity and n stdDev
        if (_vBg.Mod() > 0) {
            nout.open("out/vel_and_nStdDev" \
                        +identifier+idOutName+".out", fstream::app);
            nout << left << setw(15) << time \
                         << setw(15) << _vBg[0] \
                         << setw(15) << _vBg[1] \
                         << setw(15) << NStdDev() << endl;
            nout.close();
        }

        // Kelvin circuit
        if (_ifKelvin) {
            circuit.open("out/Kelvin_"+idOutName \
                        +to_string(newtOverDelta)+".out");
            for (int l = 0; l < _lenx; l++) {
                circuit << left \
                        << setw(15) << _x[static_cast<int>(round(_KelX[l]))] \
                        << setw(15) << _y[static_cast<int>(round(_KelY[l]))] \
                        << endl;
            }
            circuit.close();
        }
    } else if (newt%static_cast<int>(delta/5) == 0) {
        cout << "step: " << newt << ",  " << flush;
    }
}


void GPE::PinningEvolve() {
    /*
    Evolves the position of the active particles generating the pinning
        potential as dictated by their interaction with the superfluid
    */

    // Derivative requires a ComplexFunc; here, Vj = potential + i*0.
    ComplexFunc Vj(_lenx*_leny);
    Vector2 forcePsi;
    double forceBuffer;

    for (unsigned int jpart = 0; jpart < _npart; jpart++) {
        for (int i = 0; i < _lenx; i++) {
            for (int j = 0; j < _leny; j++) {
                Vj.SetValue(_leny*i + j,
                            _pinning->EvalPeak(_x[i], _y[j], jpart), 0.);
            }
        }
        for (int icomp = 0; icomp < 2; icomp++) {
            forceBuffer = 0.;
            for (int i = 0; i < _lenx; i++) {
                for (int j = 0; j < _leny; j++) {
                    forceBuffer += Dens(i,j)*Derivative(Vj, 1, icomp, i, j)[0];
                }
            }
            forcePsi.SetComp(icomp, forceBuffer*_dx*_dy);
        }
        if (_forceType == "elastic")
            _fpart[jpart] = (_xpart[jpart]-_xpEq[jpart])*(-_forceModule);
        // in the current version, only one component of the constant force
        //  can be assigned
        else if (_forceType == "constant")
            _fpart[jpart] = Vector2(_forceModule, 0.);
        _apart[jpart] = (_fpart[jpart] + forcePsi) /_mpart;
        _xpart[jpart] += _vpart[jpart]*_dt;
        _vpart[jpart] += _apart[jpart]*_dt;
    }

    _pinning->SetNewPos(_xpart);
    PinningEval();
}


void GPE::KelvinCheck() {
    /*
    Moves the points constituting the closed circuit chosen to test the validity
        of the Kelvin circulation theorem
    As imposed by the theorem, this points are moved with the local superfluid
        velocity (the circuit comoves with the superfluid)
    */

    Vector2 velCirc;
    double conv1 = _lenx/_L1;   // conversion factor from csi to n.u.
    double conv2 = _leny/_L2;
    for (int l = 0; l < _lenx; l++) {
        velCirc = Velocity(static_cast<int>(round(_KelX[l])),
                           static_cast<int>(round(_KelY[l])));
        _KelX[l] += _dt*velCirc[0]*conv1;
        _KelY[l] += _dt*velCirc[1]*conv2;
    }
}
