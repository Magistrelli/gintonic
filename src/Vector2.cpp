#include "./Vector2.h"


// &&&&&&&&&&&&&&&&
// Vector2
// &&&&&&&&&&&&&&&&


// CONSTRUCTORS

Vector2::Vector2() {
    /*
    Default constructor: initializes the vector's components to 0.
    */

    for (int r = 0; r < 2; r++) {_v[r] = 0.;}
}


Vector2::Vector2(const double comp1, const double comp2) {
    /*
    Constructs a Vector2 specifing its components
    
    :arg: <vx>  vector's first component
    :arg: <vy>  vector's second component
    */

    _v[0] = comp1;
    _v[1] = comp2;
}


Vector2::Vector2(const Vector2& vec) {
    /*
    Constructs a Vector2 instance as a copy of another one
    
    :arg: <vec>     instance copied
    */

    for (int r = 0; r < 2; r++) {_v[r] = vec._v[r];}
}


// OPERATOR OVERLOADING
Vector2& Vector2::operator= (const Vector2 &vec) {
    /*
    Overload of the assignment operator for two 2D vectors (or complex numbers):
        separately assigns both the components of the two Vector2
    
    :arg: <vec>  Vector2 instance copied
    
    Returns <*this> after have copied another specific Vector2 instance on it
    */

    // self-assignment guard
    if (this == &vec)
        return *this;
    // do the copy
    for (int r = 0; r < 2; r++) {_v[r] = vec[r];}
    return *this;
}


Vector2 Vector2::operator+ (const Vector2 &vec) const {
    /*
    Sum operator for two 2D vectors (or complex numbers): sums separately the
        two components of the two instances of Vector2
    
    :arg: <vec>  Vector2 instance summed to <*this>
    
    Returns the result of the operation <*this> + <vec>
    */

    Vector2 res;
    for (int r = 0; r < 2; r++) {res._v[r] = _v[r] + vec[r];}
    return res;
}

Vector2& Vector2::operator+= (const Vector2 &vec) {
    /*
    Overload of the operator += for 2D vectors (or complex numbers): sums
        separately the two components of the two instances of Vector2
    
    :arg: <vec>  Vector2 instance summed to <*this>
    
    Returns the result of the operation <*this> + <vec>
    */

    *this = *this + vec;
    return *this;
}


Vector2 Vector2::operator* (const double &scal) const {
    /*
    Product for a scalar coefficient: each component of Vector2 is multiplied
        by a specific scalar coefficient 
    
    :arg: <scal>    <*this> is multiplied for this scalar coefficient
    
    Returns the result of the operation <*this> * <scal>
    */

    return Vector2(_v[0]*scal, _v[1]*scal);
}

Vector2& Vector2::operator*= (const double &scal) {
    /*
    Overload of the operator *= for 2D vectors (or complex numbers): each
        component of Vector2 is multiplied by a specific scalar coefficient
    
    :arg: <scal>    <*this> is multiplied for this scalar coefficient
    
    Returns the result of the operation <*this> * <scal>
    */

    *this = *this * scal;
    return *this;
}


Vector2 Vector2::operator/ (const double &scal) const {
    /*
    Division for a scalar coefficient: each component of Vector2 is divided
        by a specific scalar coefficient
    
    :arg: <scal>    <*this> is divided for this scalar coefficient
    
    Returns the result of the operation <*this> / <scal>
    */

    return Vector2(_v[0]/scal, _v[1]/scal);
}


ostream& operator<< (ostream &os, const Vector2 &vec) {
    /*
    Overload of the operator << for 2D vectors (or complex numbers): print the
        two component of the considered instance of Vector2 by separating them
        with a comma
    
    :arg: <os>   output variable
    :arg: <vec>  Vector2 instance printed
    
    Returns the result of the output operations
    */

    return os << vec[0] << ", " << vec[1];
}




// OTHER FUNCTIONS


double Vector2::Cross(const Vector2 &vec) const {
    /*
    Vector product with a specific instance of Vector2
    
    :arg: <vec>  second argument of the cross product
    
    Returns the third component of the cross product's result
    */

    return _v[0]*vec[1] - _v[1]*vec[0];
}



double Vector2::Mod() const {
    /*
    Returns the modulus of Vector2
    */

    return sqrt(pow(_v[0], 2) + pow(_v[1], 2));
}


double Vector2::Phase() const {
    /*
    Returns the phase of Vector2 (angle with respect to the x-axes)
    */

    return atan2(_v[1], _v[0]);
}

double Vector2::GetRe() const {
    // Returns real part (component 1)
    
    return _v[0];
}

double Vector2::GetIm() const {
    // Returns imaginary part (component 2)
    
    return _v[1];
}





// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
// ComplexFunc
// &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


// CONSTRUCTORS

ComplexFunc::ComplexFunc(const unsigned int len): ComplexFunc() {
    /*
    Constructs a ComplexFunc specifing the number of evaluation points
    All the components of all the Vector2 costituing the ComplexFunc are
        set to 0
    
    :arg: <len>  number of points where the complex function is evaluated
    */

    _psi = new Vector2[len];
    _len = len;
}



// OPERATOR OVERLOADING
ComplexFunc& ComplexFunc::operator= (const ComplexFunc &func) {
    /*
    Assignment operator for two complex (or 2D) functions: separately assigns
        both the components of all the Vector2 constituting the two ComplexFunc
    
    :arg: <func>    ComplexFunc copied
    
    Returns <*this> after have copied another specific Vector2 instance on it
    */

    // self-assignment guard
    if (this == &func)
        return *this;
    // do the copy
    if (_len != func._len) {
        cerr << "Error! Assigment between two ComplexFunc of different len!" \
             << endl;
        exit(-90);
    }
    for (unsigned int i = 0; i < _len; i++) {_psi[i] = func[i];}
    return *this;
}


ComplexFunc ComplexFunc::operator+ (const ComplexFunc &func) const {
    /*
    Sum operator for two complex (or 2D) functions: sums separately the
        Vector2 instances constituting the ComplexFunc instances
    
    :arg: <func>    ComplexFunc instance summed to <*this>
    
    Returns the result of the operation <*this> + <func>
    */

    ComplexFunc res(_len);
    for (unsigned int i = 0; i < _len; i++) {res[i] = _psi[i] + func[i];}
    return res;
}

ComplexFunc& ComplexFunc::operator+= (const ComplexFunc &func) {
    /*
    Overload of the operator += for complex (or 2D) functions: sums separately
        the Vector2 instances constituting the ComplexFunc instances
    
    :arg: <func>    ComplexFunc instance summed to <*this>
    
    Returns the result of the operation <*this> + <func>
    */

    *this = *this + func;
    return *this;
}


ComplexFunc ComplexFunc::operator* (const double &scal) const {
    /*
    Product for a scalar coefficient: each component of each Vector2
        constituting the complex (or 2D) function is multiplied by a
        specific scalar coefficient 
    
    :arg: <scal>    <*this> is multiplied for this scalar coefficient
    
    Returns the result of the operation <*this> * <scal>
    */

    ComplexFunc res(_len);
    for (unsigned int i = 0; i < _len; i++) {res[i] = _psi[i]*scal;}
    return res;
}

ComplexFunc& ComplexFunc::operator*= (const double &scal) {
    /*
    Overload of the operator *= for complex (or 2D) functions: each Vector2
        constituting the ComplpexFunc is multiplied by a specific scalar
        coefficient
    
    :arg: <scal>    <*this> is multiplied for this scalar coefficient
    
    Returns the result of the operation <*this> * <scal>
    */

    *this = *this * scal;
    return *this;
}





// PUBLIC METHODS

void ComplexFunc::SetLen(const unsigned int len) {
    /*
    Deallocate the momory block pointed by <_psi> and allocate a new memory
        block of specified length
    Saves this specified length in the internal variable <_len>
    
    :arg: <len>  the pointer <_psi> will point to a memory block storing <len>
                 Vector2 instances
    */

    if (_psi) {delete[] _psi;}
    _psi = new Vector2[len];
    _len = len;
}

std::vector<double> ComplexFunc::GetComplex() const{
    std::vector<double> components(2 * _len);

    for (unsigned int i = 0; i < _len; i++) {
        components[2 * i] = _psi[i].GetRe();
        components[2 * i + 1] = _psi[i].GetIm();
    }
    return components; 
}

void ComplexFunc::SameObj(const ComplexFunc& func) {
    /*
    Completely identifies two ComplexFunc (indentify their pointers)
    In copied, saves the information that more than one pointers are pointing
        to the same memory block
    Before calling the destructor of this class (where <_copied> == 1), one MUST
        call FreePointer() to avoid double deallocations
    
    :arg: <func>    ComplexFunc instance that will be identified with <*this>
    */

    delete[] _psi;
    _len = func._len;
    _psi = func._psi;
    _copied = 1;
}


void ComplexFunc::FreePointer() {
    /*
    This method MUST be called before the destructor of a ComplexFunc where
        <_copied> == 1 to avoid double deletions (see also the method
        `SameObj`)
    Reinitializes the pointer (memory MUST be deleted explicitly somewhere,
        typically in the other ComplexFunc istance with the copy of the
        pointer _psi)
    */

    if (_copied) {
        _psi = new Vector2[1];
    } else {
        cerr \
        << "Your are lost the address of a NON-deallocated block of memory!!!" \
        << endl;
        exit(-99);
    }
}



void ComplexFunc::Print() const {
    /*
    Prints the ComplexFunc on terminal as a column of Vector2
    */

    for (unsigned int i = 0; i < _len; i++) {cout << _psi[i] << endl;}
}
