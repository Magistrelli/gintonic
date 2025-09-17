#include "./functions.h"


double Mod2pi(const double p) {
    /*
    This method does the modulo(2pi) operation on the input value
    (remember: the '%' of C++ is the remainder, not the modulo operator)
    
    :arg: <p>   argument of the mod(2pi) operation
    
    Returns a double containing the result of the operation: <p> mod(2pi)
    */

    return fmod( fmod(p, (2.*M_PI)) + (2.*M_PI), (2.*M_PI ) );
}


Vector2 Rotation(const Vector2 &vec, const double angle) {
    /*
    This method rotates clockwise an input 2D vector by a specified angle
    
    :arg: <vec>     input Vector2
    :arg: <angle>   value (in rad) of the rotation angle
    
    returns a Vector2 cointaing the components of the rotated vector
    */

    double comp1 = cos(angle)*vec[0] + sin(angle)*vec[1];
    double comp2 = -sin(angle)*vec[0] + cos(angle)*vec[1];
    return Vector2(comp1, comp2);
}



double pbcDistance(const double x, const double x0, const double L) {
    /*
    This method computes the 1-D distance (or a single component of the N-D
        distance) of a certain point from a reference one employing PBC
        (Periodic Boundary Conditions); i.e. the absolute value of this
        distance can not be larger than 1/2 the length of the cell's side)

    :arg: <x>   the position (or its i-th coordinate) of the point of which
                we want to calculate the distance from a given reference point
    :arg: <x0>  the position (or its i-th coordinate) of the reference point
    :arg: <L>   the length of the (i-th) side of the unit cell
    
    Returns a double giving (a certain component of) the distance of <x> from
        <x0> in PBC
    */

    double dist = x - x0;           // actual (component of the) distance
    while (dist > L/2.) {dist-=L;}  // reduce distance until it is less then L/2
    while (dist < -L/2.) {dist+=L;}
    return dist;
}
