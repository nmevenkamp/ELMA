#ifndef BUMPFIT_H_
#define BUMPFIT_H_


#include "nonLinearRegression.h"


template <typename RealType>
class AsymmetricGaussianBumpFunction {
public:
  static const int NumberOfParameters = 7;

private:
  const aol::Vec2<RealType> _center;
  const RealType _height;
  const RealType _sigmaX, _sigmaY;
  const RealType _rotation;
  const RealType _offset;
  const RealType _c1, _c2;

public:
  AsymmetricGaussianBumpFunction (  const aol::Vec2<RealType> &Center,
                                    const RealType Height,
                                    const RealType SigmaX, const RealType SigmaY,
                                    const RealType Rotation,
                                    const RealType Offset )
   : _center ( Center ),
     _height ( Height ),
     _sigmaX ( SigmaX ), _sigmaY ( SigmaY ),
     _rotation ( Rotation ),
     _offset ( Offset ),
     _c1 ( - 1 / ( 2 * ( 1 - aol::Sqr<RealType> ( _rotation ) ) ) ),
     _c2 ( 2 * _rotation / ( _sigmaX * _sigmaY ) ) {}

  AsymmetricGaussianBumpFunction ( const aol::Vector<RealType> &Parameters )
   : _center ( Parameters.getData() ),
     _height ( Parameters[2] ),
     _sigmaX ( Parameters[3] ), _sigmaY ( Parameters[4] ),
     _rotation ( Parameters[5] ),
     _offset ( Parameters[6] ),
     _c1 ( - 1 / ( 2 * ( 1 - aol::Sqr<RealType> ( _rotation ) ) ) ),
     _c2 ( 2 * _rotation / ( _sigmaX * _sigmaY ) ) { }

  RealType evaluate ( const aol::Vec2<RealType> &X ) const {
    aol::Vec2<RealType> pos ( X );
    pos -= _center;

    return _height * exp ( _c1 * ( aol::Sqr<RealType> ( pos[0] / _sigmaX ) + aol::Sqr<RealType> ( pos[1] / _sigmaY ) - _c2 * pos[0] * pos[1] ) ) + _offset;
  }

  void evaluateParameterGradient ( const aol::Vec2<RealType> &X, aol::Vector<RealType> &Grad ) const {
    aol::Vec2<RealType> pos ( X );
    pos -= _center;
    
    const RealType f = evaluate ( X ) - _offset;
    Grad[0] = f * _c1 * ( - 2 * pos[0] / aol::Sqr<RealType> ( _sigmaX ) + _c2 * pos[1] );
    Grad[1] = f * _c1 * ( - 2 * pos[1] / aol::Sqr<RealType> ( _sigmaY ) + _c2 * pos[0] );
    Grad[2] = f / _height;
    Grad[3] = f * _c1 * ( - 2 * aol::Sqr<RealType> ( pos[0] / _sigmaX ) / _sigmaX + _c2 / _sigmaX * pos[0] * pos[1] );
    Grad[4] = f * _c1 * ( - 2 * aol::Sqr<RealType> ( pos[1] / _sigmaY ) / _sigmaY + _c2 / _sigmaY * pos[0] * pos[1] );
    Grad[5] = f * ( pos[0] * pos[1] / ( _sigmaX * _sigmaY * ( 1 - aol::Sqr<RealType> ( _rotation ) ) )
              - _rotation / aol::Sqr<RealType> ( ( 1 - aol::Sqr<RealType> ( _rotation ) ) ) * ( aol::Sqr<RealType> ( pos[0] / _sigmaX )
              + aol::Sqr<RealType> ( pos[1] / _sigmaY ) - _c2 * pos[0] * pos[1] ) );
    Grad[6] = 1;
  }

  const aol::Vec2<RealType>& getCenter ( ) const {
    return _center;
  }
};


template <typename RealType>
class SingleAsymmetricBumpFitTargetFunctional : public aol::Op<aol::Vector<RealType> > {
  const qc::ScalarArray<RealType, qc::QC_2D> &_data;
  const qc::FastILexMapper<qc::QC_2D> _mapper;
public:
  SingleAsymmetricBumpFitTargetFunctional ( const qc::ScalarArray<RealType, qc::QC_2D> &Data )
    : _data ( Data ), _mapper ( Data ) {}

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    AsymmetricGaussianBumpFunction<RealType> bump ( Arg );
    for ( qc::RectangularIterator<qc::QC_2D> it ( _data ); it.notAtEnd ( ) ; ++it ) {
      aol::Vec<qc::QC_2D, RealType> pos;
      pos.readFromBuffer ( (*it).getData ( ) );
      const RealType g = _data.get ( *it );
      Dest[_mapper.getGlobalIndex ( *it )] += ( !aol::isNaN ( g ) ? bump.evaluate ( pos ) - g : 0 );
    }
  }
};


template <typename RealType>
class SingleAsymmetricBumpFitTargetJacobian : public aol::Op<aol::Vector<RealType>, aol::FullMatrix<RealType> > {
  const qc::ScalarArray<RealType, qc::QC_2D> &_data;
  const qc::FastILexMapper<qc::QC_2D> _mapper;
public:
  SingleAsymmetricBumpFitTargetJacobian ( const qc::ScalarArray<RealType, qc::QC_2D> &Data )
    : _data ( Data ), _mapper ( Data ) { }

  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::FullMatrix<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__);
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::FullMatrix<RealType> &Dest ) const {
    AsymmetricGaussianBumpFunction<RealType> bump ( Arg );
    aol::Vector<RealType> gradient ( bump.NumberOfParameters );
    int i = 0;
    for ( qc::RectangularIterator<qc::QC_2D> it ( _data ); it.notAtEnd ( ) ; ++it ) {
      aol::Vec<qc::QC_2D, RealType> pos;
      pos.readFromBuffer ( (*it).getData ( ) );
      if ( aol::isFinite ( _data.get ( *it ) ) )
        bump.evaluateParameterGradient ( pos, gradient );
      else {
        gradient.setZero ( );
        gradient[bump.NumberOfParameters - 1] = 1;
      }
      for ( int j=0; j<bump.NumberOfParameters ; ++j )
        Dest.set ( i, j, gradient[j] );
      ++i;
    }
  }
};


template <typename RealType>
class AsymmetricGaussianDoubleBumpFunction {
public:
  static const int NumberOfParameters = 13;
  
private:
  const aol::Vec2<RealType> _center1, _center2;
  const RealType _height1, _height2;
  const RealType _sigmaX1, _sigmaY1, _sigmaX2, _sigmaY2;
  const RealType _rotation1, _rotation2;
  const RealType _offset;
  const RealType _c11, _c21, _c12, _c22;
  
public:
  AsymmetricGaussianDoubleBumpFunction ( const aol::Vec2<RealType> &Center1, const aol::Vec2<RealType> &Center2,
                                         const RealType Height1, const RealType Height2,
                                         const RealType SigmaX1, const RealType SigmaY1, const RealType SigmaX2, const RealType SigmaY2,
                                         const RealType Rotation1, const RealType Rotation2,
                                         const RealType Offset )
  : _center1 ( Center1 ), _center2 ( Center2 ),
    _height1 ( Height1 ), _height2 ( Height2 ),
    _sigmaX1 ( SigmaX1 ), _sigmaY1 ( SigmaY1 ), _sigmaX2 ( SigmaX2 ), _sigmaY2 ( SigmaY2 ),
    _rotation1 ( Rotation1 ), _rotation2 ( Rotation2 ),
    _offset ( Offset ),
    _c11 ( - 1 / ( 2 * ( 1 - aol::Sqr<RealType> ( _rotation2 ) ) ) ),
    _c21 ( 2 * _rotation1 / ( _sigmaX1 * _sigmaY1 ) ),
    _c12 ( - 1 / ( 2 * ( 1 - aol::Sqr<RealType> ( _rotation2 ) ) ) ),
    _c22 ( 2 * _rotation2 / ( _sigmaX2 * _sigmaY2 ) ) { }
  
  AsymmetricGaussianDoubleBumpFunction ( const aol::Vector<RealType> &Parameters )
    : _center1 ( Parameters[0], Parameters[1] ), _center2 ( Parameters[2], Parameters[3] ),
      _height1 ( Parameters[4] ), _height2 ( Parameters[5] ),
      _sigmaX1 ( Parameters[6] ), _sigmaY1 ( Parameters[7] ), _sigmaX2 ( Parameters[8] ), _sigmaY2 ( Parameters[9] ),
      _rotation1 ( Parameters[10] ), _rotation2 ( Parameters[11] ),
      _offset ( Parameters[12] ),
      _c11 ( - 1 / ( 2 * ( 1 - pow ( _rotation1, 2 ) ) ) ),
      _c21 ( 2 * _rotation1 / ( _sigmaX1 * _sigmaY1 ) ),
      _c12 ( - 1 / ( 2 * ( 1 - pow ( _rotation2, 2 ) ) ) ),
      _c22 ( 2 * _rotation2 / ( _sigmaX2 * _sigmaY2 ) ) { }
  
  RealType evaluate ( const aol::Vec2<RealType> &X ) const {
    aol::Vec2<RealType> pos1 ( X ), pos2 ( X );
    pos1 -= _center1;
    pos2 -= _center2;
    
    return _height1 * exp ( _c11 * ( aol::Sqr<RealType> ( pos1[0] / _sigmaX1 )
             + aol::Sqr<RealType> ( pos1[1] / _sigmaY1 ) - _c21 * pos1[0] * pos1[1] ) )
           + _height2 * exp ( _c12 * ( aol::Sqr<RealType> ( pos2[0] / _sigmaX2 )
             + aol::Sqr<RealType> ( pos2[1] / _sigmaY2 ) - _c22 * pos2[0] * pos2[1] ) )
           + _offset;
  }
  
  void evaluateParameterGradient ( const aol::Vec2<RealType> &X, aol::Vector<RealType> &Grad ) const {
    aol::Vec2<RealType> pos1 ( X ), pos2 ( X );
    pos1 -= _center1; pos2 -= _center2;
    
    const RealType c11 = aol::Sqr<RealType> ( pos1[0] / _sigmaX1 ), c21 = aol::Sqr<RealType> ( pos1[1] / _sigmaY1 ),
                   c12 = aol::Sqr<RealType> ( pos2[0] / _sigmaX2 ), c22 = aol::Sqr<RealType> ( pos2[1] / _sigmaY2 ),
                   f1 = _height1 * exp ( _c11 * ( pow ( pos1[0] / _sigmaX1, 2 ) + pow ( pos1[1] / _sigmaY1, 2 ) - _c21 * pos1[0] * pos1[1] ) ),
                   f2 = _height2 * exp ( _c12 * ( pow ( pos2[0] / _sigmaX2, 2 ) + pow ( pos2[1] / _sigmaY2, 2 ) - _c22 * pos2[0] * pos2[1] ) );
    Grad[0] = f1 * _c11 * ( -2 * pos1[0] / aol::Sqr<RealType> ( _sigmaX1 ) + _c21 * pos1[1] );
    Grad[1] = f1 * _c11 * ( -2 * pos1[1] / aol::Sqr<RealType> ( _sigmaY1 ) + _c21 * pos1[0] );
    Grad[2] = f2 * _c12 * ( -2 * pos2[0] / aol::Sqr<RealType> ( _sigmaX2 ) + _c22 * pos2[1] );
    Grad[3] = f2 * _c12 * ( -2 * pos2[1] / aol::Sqr<RealType> ( _sigmaY2 ) + _c22 * pos2[0] );
    Grad[4] = f1 / _height1;
    Grad[5] = f2 / _height2;
    Grad[6] = f1 * _c11 * ( -2 * c11 / _sigmaX1 + _c21 / _sigmaX1 * pos1[0] * pos1[1] );
    Grad[7] = f1 * _c11 * ( -2 * c21 / _sigmaY1 + _c21 / _sigmaY1 * pos1[0] * pos1[1] );
    Grad[8] = f2 * _c12 * ( -2 * c12 / _sigmaX2 + _c22 / _sigmaX2 * pos2[0] * pos2[1] );
    Grad[9] = f2 * _c12 * ( -2 * c22 / _sigmaY2 + _c22 / _sigmaY2 * pos2[0] * pos2[1] );
    Grad[10] = f1 * ( pos1[0] * pos1[1] / ( _sigmaX1 * _sigmaY1 * ( 1 - aol::Sqr<RealType> ( _rotation1 ) ) )
               - _rotation1 / aol::Sqr<RealType> ( ( 1 - aol::Sqr<RealType> ( _rotation1 ) ) ) * ( c11 + c21 - _c21 * pos1[0] * pos1[1] ) );
    Grad[11] = f2 * ( pos2[0] * pos2[1] / ( _sigmaX2 * _sigmaY2 * ( 1 - aol::Sqr<RealType> ( _rotation2 ) ) )
               - _rotation2 / aol::Sqr<RealType> ( ( 1 - aol::Sqr<RealType> ( _rotation2 ) ) ) * ( c12 + c22 - _c22 * pos2[0] * pos2[1] ) );
    Grad[12] = 1;
  }
};


template <typename RealType>
class SingleAsymmetricDoubleBumpFitTargetFunctional : public aol::Op<aol::Vector<RealType> > {
  const qc::ScalarArray<RealType, qc::QC_2D> &_data;
  const qc::FastILexMapper<qc::QC_2D> _mapper;
public:
  SingleAsymmetricDoubleBumpFitTargetFunctional ( const qc::ScalarArray<RealType, qc::QC_2D> &Data )
    : _data ( Data ), _mapper ( Data ) {}
  
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    AsymmetricGaussianDoubleBumpFunction<RealType> bump ( Arg );
    for ( qc::RectangularIterator<qc::QC_2D> it ( _data ); it.notAtEnd ( ) ; ++it ) {
      aol::Vec<qc::QC_2D, RealType> pos;
      pos.readFromBuffer ( (*it).getData ( ) );
      const RealType g = _data.get ( *it );
      Dest[_mapper.getGlobalIndex ( *it )] += ( !aol::isNaN ( g ) ? bump.evaluate ( pos ) - g : 0 );
    }
  }
};


template <typename RealType>
class SingleAsymmetricDoubleBumpFitTargetJacobian : public aol::Op<aol::Vector<RealType>, aol::FullMatrix<RealType> > {
  const qc::ScalarArray<RealType, qc::QC_2D> &_data;
  const qc::FastILexMapper<qc::QC_2D> _mapper;
public:
  SingleAsymmetricDoubleBumpFitTargetJacobian ( const qc::ScalarArray<RealType, qc::QC_2D> &Data )
  : _data ( Data ), _mapper ( Data ) { }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::FullMatrix<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__);
  }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::FullMatrix<RealType> &Dest ) const {
    AsymmetricGaussianDoubleBumpFunction<RealType> bump ( Arg );
    aol::Vector<RealType> gradient ( bump.NumberOfParameters );
    int i = 0;
    for ( qc::RectangularIterator<qc::QC_2D> it ( _data ); it.notAtEnd ( ) ; ++it ) {
      aol::Vec<qc::QC_2D, RealType> pos;
      pos.readFromBuffer ( (*it).getData ( ) );
      if ( aol::isFinite ( _data.get ( *it ) ) )
        bump.evaluateParameterGradient ( pos, gradient );
      else {
        gradient.setZero ( );
        gradient[bump.NumberOfParameters - 1] = 1;
      }
      for ( int j=0; j<bump.NumberOfParameters ; ++j )
        Dest.set ( i, j, gradient[j] );
      ++i;
    }
  }
};


//void testAsymmetricGaussianBumpFunction ( int argc, char **argv ) {
//  aol::ParameterParser parser ( argc, argv, "bumpFunctionTest.par" );
//  aol::Vector<RealType> parameters;
//  parser.getRealVec( "parameters", parameters );
//
//  aol::Vector<RealType> positionTmp;
//  parser.getRealVec ( "position", positionTmp );
//  aol::Vec2<RealType> x ( positionTmp[0], positionTmp[1] );
//
//  AsymmetricGaussianBumpFunction<RealType> bump ( parameters );
//  aol::Vector<RealType> grad ( 6 );
//  bump.evaluateParameterGradient ( x, grad );
//
//  std::cerr << bump.evaluate ( x ) << std::endl;
//  std::cerr << grad << std::endl;
//}

#endif /* BUMPFIT_H_ */
