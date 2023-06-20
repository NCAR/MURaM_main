#ifndef __VECTOR3D_INCLUDED__
#define __VECTOR3D_INCLUDED__

/*using namespace std;
#include <iostream>*/
#include <cmath>
#include <assert.h>

template <class real> class Vector3D {

#ifndef WEAK_ENCAPSULATION
private:
#else
public:
#endif

  real x, y, z;

#ifndef WEAK_ENCAPSULATION
public:
#endif

  Vector3D() { x = 0.; y = 0.; z = 0.; }
  Vector3D(const real& xx, const real& yy, const real& zz): x(xx),y(yy),z(zz) {}

  static const Vector3D cylindrical(const real&, const real&, const real&);
  static const Vector3D spherical(const real&, const real&, const real &);

  void zero() { x = 0.; y = 0.; z = 0.; }
  real abs() const { return sqrt(x*x+y*y+z*z); }
  real sqr() const { return x*x+y*y+z*z; }

  real& operator[](int i) {
    assert( i>=0 && i<=2 );
    switch(i) {
    case 0 : return x;
    case 1 : return y;
    case 2 : return z;
    default: return x;
    };
  }
  const real& operator[](int i) const {
    assert(i >=0 && i<=2);
    switch(i) {
    case 0 : return x;
    case 1 : return y;
    case 2 : return z;
    default: return x;
    };
  }
  const Vector3D& operator+() const {
    return *this;
  }
  const Vector3D operator+(const Vector3D& A) const {
    return Vector3D(x+A.x,y+A.y,z+A.z);
  }
  const Vector3D operator-() const {
    return Vector3D(-x,-y,-z);
  }
  const Vector3D operator-(const Vector3D& A) const {
    return Vector3D(x-A.x,y-A.y,z-A.z);
  }
  const Vector3D operator*(const int a) const {
    return Vector3D(a*x,a*y,a*z);
  }
  const Vector3D operator*(const short a) const {
    return Vector3D(a*x,a*y,a*z);
  }
  const Vector3D operator*(const long a) const {
    return Vector3D(a*x,a*y,a*z);
  }
  const Vector3D operator*(const float a) const {
    return Vector3D(a*x,a*y,a*z);
  }
  const Vector3D operator*(const double& a) const {
    return Vector3D(a*x,a*y,a*z);
  }
  real operator*(const Vector3D& A) const {
    return x*A.x+y*A.y+z*A.z;
  }
  real dot(const Vector3D& A) const {
    return x*A.x+y*A.y+z*A.z;
  }
  // Vector product. Due to possible precedence conflicts
  // always enclose in brackets, i.e. write (A^B).
  const Vector3D operator^(const Vector3D& A) const {
    return Vector3D(y*A.z-z*A.y,z*A.x-x*A.z,x*A.y-y*A.x);
  }
  const Vector3D cross(const Vector3D& A) const {
    return Vector3D(y*A.z-z*A.y,z*A.x-x*A.z,x*A.y-y*A.x);
  }
  const Vector3D operator/(const int a) const {
    return Vector3D(x/a,y/a,z/a);
  }
  const Vector3D operator/(const short a) const {
    return Vector3D(x/a,y/a,z/a);
  }
  const Vector3D operator/(const long a) const {
    return Vector3D(x/a,y/a,z/a);
  }
  const Vector3D operator/(const float a) const {
    return Vector3D(x/a,y/a,z/a);
  }
  const Vector3D operator/(const double& a) const {
    return Vector3D(x/a,y/a,z/a);
  }
  const Vector3D& operator+=(const Vector3D& B) {
    x += B.x; y += B.y; z += B.z; return *this;
  }
  const Vector3D& operator-=(const Vector3D& B) {
    x -= B.x; y -= B.y; z -= B.z; return *this;
  }
  const Vector3D& operator*=(const int a) {
    x *= a; y *= a; z *= a; return *this;
  }
  const Vector3D& operator*=(const short a) {
    x *= a; y *= a; z *= a; return *this;
  }
  const Vector3D& operator*=(const long a) {
    x *= a; y *= a; z *= a; return *this;
  }
  const Vector3D& operator*=(const float a) {
    x *= a; y *= a; z *= a; return *this;
  }
  const Vector3D& operator*=(const double& a) {
    x *= a; y *= a; z *= a; return *this;
  }
  const Vector3D& operator/=(const int a) {
    x /= a; y /= a; z /= a; return *this;
  }
  const Vector3D& operator/=(const short a) {
    x /= a; y /= a; z /= a; return *this;
  }
  const Vector3D& operator/=(const long a) {
    x /= a; y /= a; z /= a; return *this;
  }
  const Vector3D& operator/=(const float a) {
    x /= a; y /= a; z /= a; return *this;
  }
  const Vector3D& operator/=(const double& a) {
    x /= a; y /= a; z /= a; return *this;
  }
  int operator==(const Vector3D& A) const {
    return x == A.x && y == A.y && z == A.z;
  }
  int operator!=(const Vector3D& A) const {
    return x != A.x || y != A.y || z != A.z;
  }
  int operator>=(const Vector3D& A) const {
    return x >= A.x && y >= A.y && z >= A.z;
  }
  int operator<=(const Vector3D& A) const {
    return x <= A.x && y <= A.y && z <= A.z;
  }
  int operator>(const Vector3D& A) const {
    return x > A.x && y > A.y && z > A.z;
  }
  int operator<(const Vector3D& A) const {
    return x < A.x && y < A.y && z < A.z;
  }
};

template <class real>
inline const Vector3D<real> Vector3D<real>::cylindrical(const real& r,
							const real& phi,
							const real& z) {
  return Vector3D(r*cos(phi),r*sin(phi),z);
}

template <class real>
inline const Vector3D<real> Vector3D<real>::spherical(const real& r,
						      const real& theta,
						      const real& phi) {
  real rsin = r*sin(theta);
  return Vector3D(rsin*cos(phi),rsin*sin(phi),r*cos(theta));
}

template <class real>
inline const Vector3D<real> operator*(const int a, const Vector3D<real>& V) {
#ifdef WEAK_ENCAPSULATION
  return Vector3D<real>(a*V.x,a*V.y,a*V.z);
#else
  return Vector3D<real>(a*V[0],a*V[1],a*V[2]);
#endif
}

template <class real>
inline const Vector3D<real> operator*(const short a, const Vector3D<real>& V) {
#ifdef WEAK_ENCAPSULATION
  return Vector3D<real>(a*V.x,a*V.y,a*V.z);
#else
  return Vector3D<real>(a*V[0],a*V[1],a*V[2]);
#endif
}

template <class real>
inline const Vector3D<real> operator*(const long a, const Vector3D<real>& V) {
#ifdef WEAK_ENCAPSULATION
  return Vector3D<real>(a*V.x,a*V.y,a*V.z);
#else
  return Vector3D<real>(a*V[0],a*V[1],a*V[2]);
#endif
}

template <class real>
inline const Vector3D<real> operator*(const float a, const Vector3D<real>& V) {
#ifdef WEAK_ENCAPSULATION
  return Vector3D<real>(a*V.x,a*V.y,a*V.z);
#else
  return Vector3D<real>(a*V[0],a*V[1],a*V[2]);
#endif
}

template <class real>
inline const Vector3D<real> operator*(const double& a, const Vector3D<real>& V) {
#ifdef WEAK_ENCAPSULATION
  return Vector3D<real>(a*V.x,a*V.y,a*V.z);
#else
  return Vector3D<real>(a*V[0],a*V[1],a*V[2]);
#endif
}

template <class real>
inline const Vector3D<real> min(const Vector3D<real>& A,
				const Vector3D<real>& B) {
#ifdef WEAK_ENCAPSULATION
  return Vector3D<real>((A.x < B.x) ? A.x : B.x,
			(A.y < B.y) ? A.y : B.y,
			(A.z < B.z) ? A.z : B.z);
#else
  return Vector3D<real>((A[0] < B[0]) ? A[0] : B[0],
			(A[1] < B[1]) ? A[1] : B[1],
			(A[2] < B[2]) ? A[2] : B[2]);
#endif
}

template <class real>
inline const Vector3D<real> max(const Vector3D<real>& A,
				const Vector3D<real>& B) {
#ifdef WEAK_ENCAPSULATION
  return Vector3D<real>((A.x > B.x) ? A.x : B.x,
			(A.y > B.y) ? A.y : B.y,
			(A.z > B.z) ? A.z : B.z);
#else
  return Vector3D<real>((A[0] > B[0]) ? A[0] : B[0],
			(A[1] > B[1]) ? A[1] : B[1],
			(A[2] > B[2]) ? A[2] : B[2]);
#endif
}

template <class real>
inline real dot(const Vector3D<real>& A, const Vector3D<real>& B) {
#ifdef WEAK_ENCAPSULATION
  return A.x*B.x+A.y*B.y+A.z*B.z;
#else
  return A[0]*B[0]+A[1]*B[1]+A[2]*B[2];
#endif
}

template <class real>
inline const Vector3D<real> cross(const Vector3D<real>& A,
				  const Vector3D<real>& B) {
#ifdef WEAK_ENCAPSULATION
  return Vector3D<real>(A.y*B.z-A.z*B.y,A.z*B.x-A.x*B.z,A.x*B.y-A.y*B.x);
#else
  return Vector3D<real>(A[1]*B[2]-A[2]*B[1],A[2]*B[0]-A[0]*B[2],A[0]*B[1]-A[1]*B[0]);
#endif
}

#endif // __VECTOR3D_INCLUDED__
