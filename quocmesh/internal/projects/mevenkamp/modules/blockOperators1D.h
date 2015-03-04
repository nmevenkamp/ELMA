#ifndef BLOCKOPERATORS1D_h
#define BLOCKOPERATORS1D_h

#include <aol.h>
#include <op.h>
#include <matrixInverse.h>
#include <preconditioner.h>
#include <wavelet.h>

enum TRANSFORM_TYPE {
  BIOR15,   // Bi-orthogonal (1.5) wavelets
  DCT,      // Discrete Cosine Transform
  HAAR      // Haar wavelet
};

class DWT {
public:
  static bool isDWT ( const std::string &DWT ) {
    for ( int i=0; i<numDWTs ; ++i ) {
      if ( DWT == dwts[i] ) return true;
    }
    return false;
  }
private:
  static const int numDWTs = 45;
  static const char* const dwts[];
};



template <typename _RealType>
class LinearUnitary1DTransformOp : aol::Op<aol::Vector<_RealType> > {
  typedef _RealType RealType;
  typedef aol::Vector<RealType> VectorType;
  typedef aol::FullMatrix<RealType> MatrixType;
protected:
  MatrixType _mTForward, _mTInverse;
public:
  LinearUnitary1DTransformOp ( const short SignalLength, const std::string &Transform ) {
    setTransformMatrices ( SignalLength, Transform );
  }
  
  virtual ~LinearUnitary1DTransformOp ( ) { }
  
  void apply ( const VectorType &Arg, VectorType &Dest ) const {
    _mTForward.apply ( Arg, Dest );
  }
  
  void applyInverse ( const VectorType &Arg, VectorType &Dest ) const {
    _mTInverse.apply ( Arg, Dest );
  }
  
  void applyAdd ( const aol::Vector<_RealType> &/*Arg*/, aol::Vector<_RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void setTransformMatrices ( const short SignalLength, const std::string &Transform ) {
    if ( SignalLength >= 2 ) {
      // Set matrix coefficients of forward transform
      _mTForward.reallocate ( SignalLength, SignalLength );
      if ( Transform == "dct" ) {
        for ( short row=0; row<SignalLength; ++row )
          for ( short col=0; col<SignalLength ; ++col )
            _mTForward.set ( row, col, cos ( aol::NumberTrait<RealType>::pi * ( 2 * col + 1 ) * row / ( 2 * SignalLength ) ) );
      } else if ( DWT::isDWT ( Transform ) ) {
        for ( short col=0; col<SignalLength; ++col ) {
          std::vector<double> signal ( SignalLength ), dwtOutput, flags;
          signal[col] = 1.0;
          dwt ( signal, log2 ( SignalLength ), Transform, dwtOutput, flags );
          for ( short row=0; row<SignalLength ; ++row )
            _mTForward.set ( row, col, dwtOutput[row] );
        }
      } else throw aol::Exception ( "Did not recognize transformation!", __FILE__, __LINE__ );
      
      // Normalize basis elements (normalize 2-Norm of the rows of the forward transform)
      RealType rowNormInv;
      for ( short row=0; row<SignalLength ; ++row ) {
        rowNormInv = 0;
        for ( short col=0; col<SignalLength ; ++col ) rowNormInv += aol::Sqr<RealType> ( _mTForward.get ( row, col ) );
        rowNormInv = 1 / sqrt ( rowNormInv );
        for ( short col=0; col<SignalLength ; ++col ) _mTForward.set ( row, col, _mTForward.get ( row, col ) * rowNormInv );
      }
      
      // Compute inverse transform matrix
      _mTInverse.reallocate ( SignalLength, SignalLength );
      aol::QRInverse<RealType> qrInverse ( _mTForward );
      _mTInverse = qrInverse.getFM ( );
    } else {
      _mTForward.reallocate ( 1, 1 );
      _mTForward.set ( 0, 0, 1.0 );
      _mTInverse.reallocate ( 1, 1 );
      _mTInverse.set ( 0, 0, 1.0 );
    }
  }
};


//template <typename _RealType>
//class WaveletTransformOp : aol::Op<aol::Vector<_RealType> > {
//  typedef _RealType RealType;
//  typedef aol::Vector<RealType> VectorType;
//protected:
//  std::string _transform;
//public:
//  WaveletTransformOp ( const std::string &Transform )
//    : _transform ( Transform ) {
//    if ( !DWT::isDWT ( Transform ) ) throw aol::Exception ( "Specified transform could not be recognized!", __FILE__, __LINE__ );
//  }
//  
//  virtual ~WaveletTransformOp ( ) { }
//  
//  void apply ( const VectorType &Arg, VectorType &Dest ) const {
//    std::vector<double> signal ( Arg.size ( ) ), dwtOutput, flags;
//    for ( short k=0; k<Arg.size ( ) ; ++k ) signal[k] = Arg[k];
//    dwt ( signal, log2 ( Arg.size ( ) ) - 1, _transform, dwtOutput, flags );
//    std::cerr << dwtOutput << std::endl;
//    std::cerr << flags << std::endl;
//    for ( short k=0; k<Dest.size ( ) ; ++k ) Dest[k] = dwtOutput[k];
//  }
//  
//  void applyInverse ( const VectorType &Arg, VectorType &Dest ) const {
//    const short pow2Size = nextPow2Size ( Arg.size ( ) );
//    std::vector<double> dwtOutput ( pow2Size ), flags ( 2 ), signal;
//    for ( short k=0; k<Arg.size ( ) ; ++k ) dwtOutput[k] = Arg[k];
//    for ( short k=Arg.size ( ) ; k<dwtOutput.size ( ) ; ++ k ) dwtOutput[k] = 0.0;
//    flags[0] = pow2Size - Arg.size ( ); flags[1] = log2 ( pow2Size ) - 2;
//    std::cerr << dwtOutput << std::endl;
//    idwt ( dwtOutput, flags, _transform, signal );
//    for ( short k=0; k<Dest.size ( ) ; ++k ) Dest[k] = signal[k];
//  }
//  
//  void applyAdd ( const aol::Vector<_RealType> &/*Arg*/, aol::Vector<_RealType> &/*Dest*/ ) const {
//    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
//  }
//  
//  void setTransform ( const std::string &Transform ) {
//    if ( DWT::isDWT ( Transform ) ) _transform = Transform;
//    else throw aol::Exception ( "Specified transform could not be recognized!", __FILE__, __LINE__ );
//  }
//
//private:
//  short nextPow2Size ( short Size ) const {
//    short res = 1;
//    while ( res < Size ) res *= 2;
//    return res;
//  }
//};


#endif
