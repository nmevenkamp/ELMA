#ifndef BLOCKOPERATORS2D_h
#define BLOCKOPERATORS2D_h

#include <array.h>
#include <convolution.h>
#include "blockOperators1D.h"
#ifdef USE_BOOST
#ifndef Q_MOC_RUN // See: https://bugreports.qt-project.org/browse/QTBUG-22829
#include <boost/math/special_functions/bessel.hpp>
#endif
#endif



// beta = pi * alpha
template <typename _RealType>
void setKaiserWindow ( qc::ScalarArray<_RealType, qc::QC_2D> &Window, const short Size, const _RealType Beta ) {
  aol::Vector<_RealType> kaiser ( Size );
  for ( int i=0; i<Size ; ++i )
#ifdef USE_BOOST
    kaiser[i] = static_cast<_RealType> ( boost::math::cyl_bessel_i<int, _RealType> ( 0, Beta * sqrt ( 1 - aol::Sqr<_RealType> ( 2.0 * i / ( Size - 1 ) - 1 ) ) )
                                       / boost::math::cyl_bessel_i<int, _RealType> ( 0, Beta ) );
#else
    throw aol::Exception ( "Boost required! Compile with -DUSE_BOOST=1", __FILE__, __LINE__ );
#endif
  
  Window.reallocate ( Size, Size );
  for ( int x=0; x<Size ; ++x )
    for ( int y=0; y<Size ; ++y )
      Window.set ( x, y, kaiser[x] * kaiser[y] );
}


template <typename _RealType>
void setGaussianWindow ( qc::ScalarArray<_RealType, qc::QC_2D> &Window, const short Size ) {
  if ( Size % 2 == 0 )
    throw aol::Exception ( "Window size must be odd!", __FILE__, __LINE__ );
  
  Window.reallocate ( Size, Size );
  const int offset = ( Size - 1 ) / 2;
  const _RealType stdDev = 0.5 * Size, c1 = 1.0 / ( stdDev * sqrt ( 2 * aol::NumberTrait<_RealType>::pi ) ), c2 = 2 * aol::Sqr<_RealType> ( stdDev );
  for ( int i=-offset; i<=offset ; ++i )
    for ( int j=-offset; j<=offset ; ++j )
      Window.set ( i + offset, j + offset, c1 * exp ( -( aol::Pow ( i, 2 ) + aol::Pow ( j, 2 ) ) / c2 ) );
}


template <typename _RealType>
class HardThresholding2DOp : aol::Op<qc::ScalarArray<_RealType, qc::QC_2D> > {
  typedef _RealType RealType;
  typedef qc::ScalarArray<RealType, qc::QC_2D> BlockType;
protected:
  BlockType _thresholds;
public:
  HardThresholding2DOp ( ) : _thresholds ( BlockType ( 1, 1 ) ) { }
  
  HardThresholding2DOp ( const BlockType &Thresholds ) : _thresholds ( Thresholds ) { }
  
  void apply ( const BlockType &Arg, BlockType &Dest ) const {
    for ( short x=0; x<Arg.getNumX ( ) ; ++x )
      for ( short y=0; y<Arg.getNumY ( ) ; ++y )
        Dest.set ( x, y, aol::Abs<RealType> ( Arg.get ( x, y ) ) > _thresholds.get ( x, y ) ? Arg.get ( x, y ) : 0.0 );
  }
  
  virtual void applyAdd ( const BlockType &/*Arg*/, BlockType &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void setThresholds ( const BlockType &Thresholds ) {
    _thresholds.reallocate ( Thresholds.getNumX ( ), Thresholds.getNumY ( ) );
    _thresholds = Thresholds;
  }
};


template <typename _RealType>
class LinearUnitary2DTransformOp : aol::Op<qc::ScalarArray<_RealType, qc::QC_2D> > {
  typedef _RealType RealType;
  typedef qc::ScalarArray<RealType, qc::QC_2D> BlockType;
protected:
  LinearUnitary1DTransformOp<RealType> _transform1D;
public:
  LinearUnitary2DTransformOp ( const short BlockSize, const std::string &Transform )
    : _transform1D ( BlockSize, Transform ) { }
  
  void apply ( const BlockType &Arg, BlockType &Dest ) const {
    // Apply row-column algorithm
    qc::Array<RealType> tmpArg ( Arg.getNumXYZ ( ), 1, 1 ), tmpDest ( Dest.getNumXYZ ( ), 1, 1 );
    
    // Transform rows
    for ( short y=0; y<Arg.getNumXYZ ( ) ; ++y ) {
      Arg.getSlice ( qc::QC_X, y, tmpArg );      
      _transform1D.apply ( tmpArg, tmpDest );
      Dest.putSlice ( qc::QC_X, y, tmpDest );
    }
    
    // Transform columns
    for ( short x=0; x<Arg.getNumXYZ ( ) ; ++x ) {
      Dest.getSlice ( qc::QC_Y, x, tmpArg );
      _transform1D.apply ( tmpArg, tmpDest );
      Dest.putSlice ( qc::QC_Y, x, tmpDest );
    }
  }
  
  void applyInverse ( const BlockType &Arg, BlockType &Dest ) const {
    // Apply row-column algorithm
    qc::Array<RealType> tmpArg ( Arg.getNumXYZ ( ), 1, 1 ), tmpDest ( Dest.getNumXYZ ( ), 1, 1 );
    
    // Transform rows
    for ( short y=0; y<Arg.getNumXYZ ( ) ; ++y ) {
      Arg.getSlice ( qc::QC_X, y, tmpArg );
      _transform1D.applyInverse ( tmpArg, tmpDest );
      Dest.putSlice ( qc::QC_X, y, tmpDest );
    }
    
    // Transform columns
    for ( short x=0; x<Arg.getNumXYZ ( ) ; ++x ) {
      Dest.getSlice ( qc::QC_Y, x, tmpArg );
      _transform1D.applyInverse ( tmpArg, tmpDest );
      Dest.putSlice ( qc::QC_Y, x, tmpDest );
    }
  }
  
  void applyAdd ( const BlockType &/*Arg*/, BlockType &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void setTransform ( const short BlockSize, const std::string &Transform ) {
    _transform1D.setTransformMatrices ( BlockSize, Transform );
  }
};


template <typename _RealType>
class BlockDenoisingOp : public aol::Op<qc::ScalarArray<_RealType, qc::QC_2D> > {
  typedef qc::ScalarArray<_RealType, qc::QC_2D> BlockType;
public:
  BlockDenoisingOp ( ) { }
  
  virtual ~BlockDenoisingOp ( ) { }
  
  virtual void apply ( const BlockType &/*Arg*/, BlockType &/*Dest*/ ) const = 0;
  
  void applyAdd ( const BlockType &/*Arg*/, BlockType &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};


template <typename _RealType>
class HTDenoisingInLocal2DTransformDomainOp : public BlockDenoisingOp<_RealType> {
  typedef _RealType RealType;
  typedef qc::ScalarArray<RealType, qc::QC_2D> BlockType;
protected:
  LinearUnitary2DTransformOp<RealType> _2DTransformOp;
  HardThresholding2DOp<RealType> _HTOp;
public:
  HTDenoisingInLocal2DTransformDomainOp ( )
    : _2DTransformOp ( 1, "dct" ), _HTOp ( BlockType ( 1, 1 ) ) { }
  
  HTDenoisingInLocal2DTransformDomainOp ( const short BlockSize, const std::string &Transform, const BlockType &Thresholds )
    : _2DTransformOp ( BlockSize, Transform ), _HTOp ( Thresholds ) { }
  
  void apply ( const BlockType &Arg, BlockType &Dest ) const {
    _2DTransformOp.apply ( Arg, Dest );
    _HTOp.apply ( Dest, Dest );
  }
  
  void setTransform ( const short BlockSize, const std::string &Transform ) {
    _2DTransformOp.setTransform ( BlockSize, Transform );
  }
  
  void setThresholds ( const BlockType &Thresholds ) {
    _HTOp.setThresholds ( Thresholds );
  }
};


template <typename _RealType>
class HTDenoisingInLocal2DFFTDomainOp : public BlockDenoisingOp<_RealType> {
  typedef _RealType RealType;
  typedef qc::ScalarArray<RealType, qc::QC_2D> BlockType;
protected:
  RealType _threshold;
public:
  HTDenoisingInLocal2DFFTDomainOp ( )
    : _threshold ( 0.0 ) { }
  
  HTDenoisingInLocal2DFFTDomainOp ( const RealType Threshold )
    : _threshold ( Threshold ) { }
  
  void apply ( const BlockType &Arg, BlockType &Dest ) const {
    qc::MultiArray<RealType, 2, 2> _function ( Arg, aol::STRUCT_COPY ), _transform ( Arg, aol::STRUCT_COPY );
    _function[0] = Arg;
    _function[1].setZero ( );
    qc::FourierTransform<RealType>	(	_function, _transform );
    for ( qc::RectangularIterator<qc::QC_2D> it ( Arg ); it.notAtEnd ( ) ; ++it ) {
      if ( _transform.get ( *it ).norm ( ) < _threshold )
        _transform.set ( *it, aol::Vec2<RealType> ( 0, 0 ) );
      Dest.set ( *it, _transform.get ( *it ).norm ( ) );
    }
  }
  
  void setThreshold ( const RealType Threshold ) {
    _threshold = Threshold;
  }
};


template <typename _RealType>
class ZeroMeanOp : public BlockDenoisingOp<_RealType> {
  typedef _RealType RealType;
  typedef aol::FullMatrix<RealType> MatrixType;
  typedef qc::ScalarArray<_RealType, qc::QC_2D> BlockType;
public:
  ZeroMeanOp ( ) { }
  
  void apply ( const BlockType &Arg, BlockType &Dest ) const {
    RealType mean = 0;
    for ( short x=0; x<Arg.getNumX ( ) ; ++x )
      for ( short y=0; y<Arg.getNumY ( ) ; ++y )
        mean += Arg.get ( x, y );
    mean /= Arg.getNumX( ) * Arg.getNumY ( );
    for ( short x=0; x<Dest.getNumX ( ) ; ++x )
      for ( short y=0; y<Dest.getNumY ( ) ; ++y )
        Dest.set ( x, y, Arg.get ( x, y ) - mean );
  }
};

#endif
