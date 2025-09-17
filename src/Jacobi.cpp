#include "./Jacobi.h"


// CONSTRUCTORS

// Default constructor: initializes all the numerical variables to 0 or 0.
//  except for _upSum = 30
Jacobi::Jacobi(): _q(), _qMod(0.), _qPhase(0.), _l1(), \
                  _l2(), _rotl2(), _tau(), _upSum(30), \
                  _angle(0.), _lScale(0.), _posv(), _rotXv() {}


Jacobi::Jacobi(const Vector2 &l1, const Vector2 &l2): Jacobi() {
    /*
    Constructs the Jacobi theta function related to a vortex lattice with
        specified fundamental vectors
    Because of the assumption described in the preamble of './Jacobi.h',
        this fundamental vectors are also the sides of the fundamental cell
    
    :arg: <l1>  first fundamental vector, the one that becomes (pi,0) after
                the rotation and rescaling of the vortex lattice
    :arg: <l2>  second fundamental vector
    */

    _l1 = l1;
    _l2 = l2;
    // copmute length and inclination of l1
    double modl1 = sqrt(pow(l1[0], 2.) + pow(l1[1], 2.));
    _angle = atan2(l1[1], l1[0]);
    // rescale and rotate l2 clockwise (coordinate system where l1 =(pi,0))
    _lScale = M_PI/modl1;
    _rotl2 = Rotation(l2*_lScale, _angle);
    // copmute tau and q after the rotation and rescaling
    _tau = _rotl2/M_PI;
    _qFromTau();
}


Jacobi::Jacobi(const Vector2 &l1, const Vector2 &l2,
               const Vector2 &posv): Jacobi(l1, l2) {
    /*
    Constructs the Jacobi theta function related to a vortex lattice with
        specified fundamental vectors and offset (position of a certain
        referral vortex)
    Because of the assumption described in the preamble of './Jacobi.h',
        the lattice's fundamental vectors are also the sides of the
        fundamental cell
    
    :arg: <l1>  first fundamental vector, the one that becomes (pi,0) after
                the rotation and rescaling of the vortex lattice
    :arg: <l2>  second fundamental vector
    :arg: <xv>  first component of the referral vortex position
    :arg: <yv>  second fcomponent of the referral vortex position
    */

    _posv = posv;
    // check that the referral vortex is inside the fundamental cell
    CheckXV();
    // rescale and rotate the referral vortex position
    _rotXv = Rotation(posv*_lScale, _angle);
}



// PRIVATE FUNCTIONS

void Jacobi::_qFromTau() {
    /*
    Computes the nome of the Jacobi Theta1 function as q = exp(i pi tau[0] +
        - pi tau[1]) and saves the result in the internal variable <_q>
    Computes also the nome's module and phase and saves the results in <_qMod>
        and <_qPhase>
    */

    double ept = exp(-M_PI*_tau[1]);
    double q1 = ept*cos(M_PI*_tau[0]);
    double q2 = ept*sin(M_PI*_tau[0]);
    _q = Vector2(q1, q2);
    _qMod = _q.Mod();
    _qPhase = _q.Phase();
}



void Jacobi::CheckXV() const {
    /*
    Method that verifies that the input position of the referral vortex lies
        inside the fundamental cell
    Exits the program with ValueError=-1 if this condition is not satisfied
    */

    // find the indices of the cell where the referral vortex position is
    vector<int> xycell = XYcell(_posv);
    // exit the program if those indices are not (0,0)
    if ((xycell[0] != 0) || (xycell[1] != 0)) {
        cerr << endl << "central vortex outside of the fundamental cell!!!" \
        <<endl;
        exit(-1);
    }
}



vector<int> Jacobi::XYcell(const Vector2 &pos) const {
    /*
    Finds the indices of the cell where a specified position lies
    The fundamental cell is (0,0)
    This method is written to work in the rotated and rescaled coordinate
        system where l1=(pi,0)
    
    :arg: <pos>  specified inestigated position
    
    Returns a vector<int> of length 2 storing the indices of the cell where
        the position pos=(x,y) lies: if this cell is xcell*pi + ycell*rotl2
        then the method returns (xcell, ycell)
    */

    vector<int> xycell(2);
    // rescale and rotate pos = (x,y)
    Vector2 rotx = Rotation(pos*_lScale, _angle);
    // distance from the bottom left angle of the fundamental cell
    double xpr = rotx[0] + 0.5*(M_PI+_rotl2[0]);  // M_PI=rotl1[0]
    double ypr = rotx[1] + 0.5*_rotl2[1];
    // compute the cell's indices
    xycell[0] = floor((xpr - ypr*(_rotl2[0]/_rotl2[1]))/M_PI);
    xycell[1] = floor(ypr/_rotl2[1]);
    return xycell;
}



// PUBLIC FUNCTIONS

Vector2 Jacobi::Eval(const Vector2 &rotPos) const {
    /*
    Evaluates the Jacobi Theta1(z,q) function function in a specified point
    This method assumes to be in the rotated and rescaled coordinate system
        where l1=(pi,0)
    
    :arg: <rotPos>  evaluation point's rotated and rescaled position
    
    Returns a Vector2 storing the real and imaginary parts of the complex value
        of the Jacobi Theta1 evaluated in <rotPos>
    */

    double Treal = 0., Tim = 0.;  // Jacobi Theta's real and imaginary parts
    double amp = 0., arg = 0.;    // the sin's coefficient in the Jacobi sum
                                  // will be written as amp*(cos(arg)+isin(arg))
    double csh = 0., sch = 0.;    // real and imaginary parts of sin((2n+1)z)
    Vector2 z = rotPos - _rotXv;  // rescaled and rotated z = z1 + i z2
    for (unsigned int n = 0; n < _upSum; n++) {
        amp = 2. * pow(-1., n) * pow(_qMod, pow(n+0.5, 2.));
        csh = sin((2.*n+1)*z[0]) * cosh((2.*n+1)*z[1]);
        sch = cos((2.*n+1)*z[0]) * sinh((2.*n+1)*z[1]);
        arg = pow(n+0.5, 2.)*_qPhase;
        Treal += amp*(cos(arg)*csh - sin(arg)*sch);
        Tim += amp*(sin(arg)*csh + cos(arg)*sch);
    }
    return Vector2(Treal, Tim);
}



double Jacobi::Theta(const Vector2 &pos) const {
    /*
    Rotates and rescale the position of a specified point following the
        coordinate transformation defined in the class constructors
    Then, this method evaluates the Wood-like (Jacobi-based) superfluid phase
        field thetaJ(x|x_v) in this transformed specified point
    Therefore, this method evaluates the actual value that the general
        (arbitrary lattice's fundamental vectors) Wood-like phase field assumes
        in the original (before rotation and rescaling) position specified by
        the arguments of this method
    
    :arg: <pos>  evaluation point's original position
    
    Returns a double storing the value of the general (l1 and l2 stored in the
        class' internal variables) Wood-like phase evaluated in <pos>
    */

    // rescale and rotate the position vector
    Vector2 rotPos = Rotation(pos*_lScale, _angle);
    // evaluate the Jacobi Theta1 function's argument
    double appo2 = Eval(rotPos).Phase();
    // evaluate the phase function (rotation around the origin)
    appo2 = rotPos[0]*(rotPos[1] - 2.*_rotXv[1])/(M_PI*_tau[1]) + appo2;
    return Mod2pi(appo2);
}



double Jacobi::QPBC(const Vector2 &pos) const {
    /*
    Evaluates the general (arbitrary lattice's fundamental vectors) Wood-like
        Jacobi-based phase field with the use of the Wood-like Quasi-Periodic
        Boundary Condition (QPBC)
    
    :arg: <pos>  evaluation point's
    
    Returns a double storing the QPBC-ed evaluation of the Wood-like phase
        field in <pos>
    */

    // compute the coefficients l_i of (6.21)
    vector<int> xycell = XYcell(pos);
    // compute the components of the translation vector L = l1 L1 + l2 L2
    Vector2 L(xycell[0]*_l1[0] + xycell[1]*_l2[0],
              xycell[0]*_l1[1] + xycell[1]*_l2[1]);
    // compute each term of (6.21), assume _Nv=1
    double vprod = (M_PI/(_l1[0]*_l2[1])) * L.Cross(pos);
    double c1 = Theta(_l1/2.) - Theta(_l1/(-2.));
    double c2 = Theta(_l2/2.) - Theta(_l2/(-2.));
    double sumi = xycell[0]*c1 + xycell[1]*c2;      // sum over i term
    return Mod2pi(Theta(pos-L) + vprod + sumi + xycell[0]*xycell[1]*M_PI);
}
