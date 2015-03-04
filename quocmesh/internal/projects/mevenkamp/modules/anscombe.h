#ifndef ANSCOMBE_H
#define ANSCOMBE_H

#include <aol.h>
#include <op.h>

enum ANSCOMBE_INVERSE_TYPE {
  ANSCOMBE_INV_ASYMPTOTICALLY_UNBIASED,
  ANSCOMBE_INV_EXACT_UNBIASED
};


template<typename _RealType>
class AnscombeForward : public aol::Op<aol::Vector<_RealType> > {
  typedef _RealType RealType;
  
  const static int anscombeLastIdx;
  const static double anscombeEfz[];
  const static double anscombeEz[];
public:
  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    const RealType a = 3.0 / 8.0;
    for ( int k = 0; k < Arg.size ( ) ; ++k )
      Dest[k] = 2 * sqrt ( Arg[k] + a );
  }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::Vector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};

template<typename _RealType>
class AnscombeInverse : public aol::Op<aol::Vector<_RealType> > {
  typedef _RealType RealType;
public:
  const static int anscombeLastIdx;
  const static double anscombeEfz[];
  const static double anscombeEz[];
public:
  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    const RealType a = 3.0 / 8.0, b = 1.0 / 8.0;
    for ( int k = 0; k < Arg.size ( ) ; ++k ) {
      if ( Arg[k] < 2 * sqrt ( a ) )
        Dest[k] = 0;
      else if ( Arg[k] > anscombeEfz[anscombeLastIdx] )
        Dest[k] = aol::Sqr<RealType> ( 0.5 * Arg[k] ) - b;
      else
        Dest[k] = interpInverse ( Arg[k] );
    }
  }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::Vector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void applyInverseOfExactAlgebraicInverse ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) {
    const RealType b = 1.0 / 8.0;
    for ( int k = 0; k < Arg.size ( ) ; ++k ) {
      const RealType y = 2 * sqrt ( Arg[k] + b );
      if ( y  > anscombeEfz[anscombeLastIdx] )
        Dest[k] = y;
      else
        Dest[k] = interp ( Arg[k] );
    }
  }
  
private:
  RealType interpInverse ( const RealType Val ) const {
    // Return extrapolated value if D beyond range
    if ( Val > anscombeEfz[anscombeLastIdx] )
      return ( Val - anscombeEfz[anscombeLastIdx-1] ) * ( anscombeEz[anscombeLastIdx] - anscombeEz[anscombeLastIdx-1] )
           / ( anscombeEfz[anscombeLastIdx] - anscombeEfz[anscombeLastIdx-1] ) + anscombeEz[anscombeLastIdx-1];
    
    // Perform binary search for k s.t. D in [Efz[k], Efz[k+1])
    int a = 0, b = anscombeLastIdx, k = ( a + b ) / 2;
    while ( Val < anscombeEfz[k] || Val > anscombeEfz[k+1] ) {
      if ( Val < anscombeEfz[k] ) {
        b = k-1;
        k = ( b + a ) / 2;
      } else {
        a = k+1;
        k = ( b + a ) / 2;
      }
    }
    
    // Return interpolated value on [Efz[k], Efz[k+1]]
    return ( Val - anscombeEfz[k-1] ) * ( anscombeEz[k] - anscombeEz[k-1] ) / ( anscombeEfz[k] - anscombeEfz[k-1] ) + anscombeEz[k-1];
  }
  
  RealType interp ( const RealType Val ) const {
    // Return extrapolated value if D beyond range
    if ( Val > anscombeEz[anscombeLastIdx] )
      return ( Val - anscombeEz[anscombeLastIdx-1] ) * ( anscombeEfz[anscombeLastIdx] - anscombeEfz[anscombeLastIdx-1] )
      / ( anscombeEz[anscombeLastIdx] - anscombeEz[anscombeLastIdx-1] ) + anscombeEfz[anscombeLastIdx-1];
    
    // Perform binary search for k s.t. D in [Efz[k], Efz[k+1])
    int a = 0, b = anscombeLastIdx, k = ( a + b ) / 2;
    while ( Val < anscombeEz[k] || Val > anscombeEz[k+1] ) {
      if ( Val < anscombeEz[k] ) {
        b = k-1;
        k = ( b + a ) / 2;
      } else {
        a = k+1;
        k = ( b + a ) / 2;
      }
    }
    
    // Return interpolated value on [Efz[k], Efz[k+1]]
    return ( Val - anscombeEz[k-1] ) * ( anscombeEfz[k] - anscombeEfz[k-1] ) / ( anscombeEz[k] - anscombeEz[k-1] ) + anscombeEfz[k-1];
  }
};

#endif
