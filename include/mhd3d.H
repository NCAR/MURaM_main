#ifndef __MHD3D_INCLUDED__
#define __MHD3D_INCLUDED__

/*using namespace std;
#include <iostream>
#include <stdlib.h>
#include <cmath>*/
#include <vector3d.H>

namespace MHD {

  // Declare variable classes
  template <class real> class cState; // conserved

  // Define conserved variable class
  template <class real> class cState{

#ifndef WEAK_ENCAPSULATION
  private:
#else
  public:
#endif

    real           d; // Density
    Vector3D<real> M; // Momentum
    real           e; // Energy
    Vector3D<real> B; // Magnetic field
                      // population densities & time derivatives of population densities ...
#ifndef WEAK_ENCAPSULATION
  public:
#endif

    cState() { d = 0.; M.zero(); e = 0.; B.zero(); }
    cState(const real& dd, const Vector3D<real>& mm,
	   const real& ee, const Vector3D<real>& bb):
      d(dd),M(mm),e(ee),B(bb) {}
    cState(const real& dd, const real& mx, const real& my, const real& mz,
	   const real& ee, const real& bx, const real& by, const real& bz):
      d(dd),M(mx,my,mz),e(ee),B(bx,by,bz) {}

    // Utility functions
    void zero() { d = 0.; M.zero(); e = 0.; B.zero(); }
    int positive() const { return d>0. && (e-0.5*(M.sqr()/d+B.sqr()))>0.; }

    // Public interfaces for private members
    real& Density() { return d; }
    const real& Density() const { return d; }

    Vector3D<real>& Momentum() { return M; }
    const Vector3D<real>& Momentum() const { return M; }

    real& Energy() { return e; }
    const real& Energy() const { return e; }

    Vector3D<real>& MagneticField() { return B; }
    const Vector3D<real>& MagneticField() const { return B; }

    // Overloaded operators
    const cState& operator+() const {
      return *this;
    }
    const cState operator+(const cState& U) const {
      return cState(d+U.d,M+U.M,e+U.e,B+U.B);
    }
    const cState operator-(const cState& U) const {
      return cState(d-U.d,M-U.M,e-U.e,B-U.B);
    }
    const cState operator-() const {
      return cState(-d,-M,-e,-B);
    }
    const cState operator*(const real& a) const {
      return cState(d*a,M*a,e*a,B*a);
    }
    const cState operator/(const real& a) const {
      return cState(d/a,M/a,e/a,B/a);
    }
    const cState& operator+=(const cState& U) {
      d += U.d; M += U.M; e += U.e; B += U.B; return *this;
    }
    const cState& operator-=(const cState& U) {
      d -= U.d; M -= U.M; e -= U.e; B -= U.B; return *this;
    }
    const cState& operator*=(const real& a) {
      d *= a; M *= a; e *= a; B *= a; return *this;
    }
    const cState& operator/=(const real& a) {
      d /= a; M /= a; e /= a; B /= a; return *this;
    }
    int operator==(const cState& U) const {
      return d==U.d && M==U.M && e==U.e && B==U.B;
    }
    int operator!=(const cState& U) const {
      return d!=U.d || M!=U.M || e!=U.e || B!=U.B;
    }
    
    const Vector3D<real> Velocity() const { return M/d; }
  };
}

template <class real>
inline const MHD::cState<real> operator*(const int a,
					 const MHD::cState<real>& U) {
#ifdef WEAK_ENCAPSULATION
  return MHD::cState<real>(U.d*a,U.M*a,U.e*a,U.B*a);
#else
  return MHD::cState<real>(a*U.Density(),U.Momentum()*a,
			   a*U.Energy(),U.MagneticField()*a);
#endif
}

template <class real>
inline const MHD::cState<real> operator*(const float a,
					 const MHD::cState<real>& U) {
#ifdef WEAK_ENCAPSULATION
  return MHD::cState<real>(U.d*a,U.M*a,U.e*a,U.B*a);
#else
  return MHD::cState<real>(a*U.Density(),U.Momentum()*a,
			   a*U.Energy(),U.MagneticField()*a);
#endif
}

template <class real>
inline const MHD::cState<real> operator*(const double& a,
					 const MHD::cState<real>& U) {
#ifdef WEAK_ENCAPSULATION
  return MHD::cState<real>(U.d*a,U.M*a,U.e*a,U.B*a);
#else
  return MHD::cState<real>(a*U.Density(),U.Momentum()*a,
			   a*U.Energy(),U.MagneticField()*a);
#endif
}

#endif // __MHD3D_INCLUDED__
