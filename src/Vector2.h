/*
This library defines the following classes:
- Vector2     ->  represents a vector of dimension 2 (complex number or 2D
                  vector)
- ComplexFunc ->  represents a complex function, i.e. a discretized field
                  of complex number, i.e. instances of Vector2)
*/

#ifndef VECTOR2_H_
#define VECTOR2_H_

#include <cmath>
#include <iostream>
#include <vector>

using std::cout;
using std::cerr;
using std::endl;
using std::ostream;

class Vector2 {
 private:
    double _v[2];   // vector's components

 public:
    // constructors
    Vector2();
    Vector2(const double comp1, const double comp2);
    Vector2(const Vector2& vec);
    // destructor
    ~Vector2() {}

    // operators overloading
    double operator[](const unsigned int i) const {return _v[i];}
    Vector2& operator=(const Vector2 &vec);
    Vector2 operator+(const Vector2 &vec) const;
    Vector2& operator+=(const Vector2 &vec);
    Vector2 operator*(const double &scal) const;
    Vector2& operator*=(const double &scal);
    Vector2 operator/(const double &scal) const;
    Vector2 operator-(const Vector2 &v2) const {return *this + v2*(-1);}
    friend ostream& operator<<(ostream &os, const Vector2 &vec);

    // components acces
    void SetComp(const unsigned int r, const double comp) {_v[r] = comp;}
    // vector product
    double Cross(const Vector2 &vec) const;
    // modulus and phase
    double Mod() const;
    double Phase() const;
    
    double GetRe() const; //Real part
    double GetIm() const; //Imaginary part
    
};


class ComplexFunc {
 private:
    Vector2 *_psi;      // complex function (<_len> Vector2)
    unsigned int _len;  // number of the spatial grid points
    bool _copied;       // if another pointer that points where _psi does exists

 public:
    // constructors
    ComplexFunc(): _psi(NULL), _len(0), _copied(0) {}
    explicit ComplexFunc(const unsigned int len);
    // destructor
    ~ComplexFunc() {delete[] _psi;}

    // operators overloading
    Vector2& operator[](const int i) const {return _psi[i];}
    ComplexFunc& operator= (const ComplexFunc &func);
    ComplexFunc operator+(const ComplexFunc &func) const;
    ComplexFunc& operator+=(const ComplexFunc &func);
    ComplexFunc operator*(const double &scal) const;
    ComplexFunc& operator*=(const double &scal);
    ComplexFunc operator/(const double &scal) const {return *this * (1/scal);}
    ComplexFunc operator-(const ComplexFunc &func) const \
                                                   {return *this + func*(-1);}

    // variables acces
    void SetLen(const unsigned int len);
    void SetValue(const int gpos,
                  const double re,
                  const double im) {_psi[gpos] = Vector2(re, im);}
    std::vector<double> GetComplex() const;
    // identify two ComplexFunc (indentify their pointers)
    void SameObj(const ComplexFunc& func);
    // reinitialize the pointer (memory MUST be deleted explicitly somewhere)
    void FreePointer();
    // print ComplexFunc on terminal (debug method)
    void Print() const;
};


#endif  // VECTOR2_H_
