#ifndef BLOCKREGULARIZATION_H_
#define BLOCKREGULARIZATION_H_

#include "scalarArrayExtensions.h"
#include <indexMapper.h>
#include "linearRegression.h"
#include "nonLinearRegression.h"


static const double DEFAULT_EPSILON = 0.001;


template <typename _RealType, typename _PictureType>
class BlockRegularityTargetFunctional : public aol::Op<aol::Vector<_RealType> > {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  const PictureType &_data;
  const qc::FastILexMapper<qc::QC_2D> _mapper;
  const aol::Vec2<short> _pos;
  const short _width, _height;
  const int _numPixels;
public:
  BlockRegularityTargetFunctional ( const PictureType &Data, const aol::Vec2<short> &Pos, const aol::Vec2<short> &Size )
    : _data ( Data ), _mapper ( Size[0], Size[1] ), _pos ( Pos ), _width ( Size[0] ), _height ( Size[1] ), _numPixels ( Size[0] * Size[1] ) { }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( Arg.size ( ) != _height )
      throw aol::Exception ( "Arguments size does not match patch height!", __FILE__, __LINE__ );

    if ( Dest.size ( ) != _numPixels )
      throw aol::Exception ( "Destination vector does not match patch size!", __FILE__, __LINE__ );

    for ( short xBlock=0; xBlock<_width ; ++xBlock ) {
      for ( short yBlock=0; yBlock<_height ; ++yBlock ) {
        const short x = _pos[0] + xBlock, y = _pos[1] + yBlock;
        Dest[_mapper.getGlobalIndex ( xBlock, yBlock )] += _data.interpolate ( x + ( ( yBlock + 1 < _height ) ? Arg[yBlock + 1] : 0 ), y + 1 ) - _data.interpolate ( x + Arg[yBlock], y );
      }
    }
  }

  RealType dxData ( const RealType X, const short Y ) const {
    const short i = floor ( X ), j = Y;
    return _data.get ( i + 1, j ) - _data.get ( i, j );
  }

  int getGlobalIndex ( const short X, const short Y ) const {
    return _mapper.getGlobalIndex ( X, Y );
  }

  short getWidth ( ) const {
    return _width;
  }

  short getHeight ( ) const {
    return _height;
  }

  short getNumArgs ( ) const {
    return _height;
  }

  int getNumPixels ( ) const {
    return _numPixels;
  }
};


template <typename _RealType, typename _PictureType>
class BlockRegularityTargetJacobian : public aol::Op<aol::Vector<_RealType>, aol::FullMatrix<_RealType> > {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  const PictureType &_data;
  const qc::FastILexMapper<qc::QC_2D> _mapper;
  const aol::Vec2<short> _pos;
  const short _width, _height;
  const int _numPixels;
public:
  BlockRegularityTargetJacobian ( const PictureType &Data, const aol::Vec2<short> &Pos, const aol::Vec2<short> &Size )
    : _data ( Data ), _mapper ( Size[0], Size[1] ), _pos ( Pos ), _width ( Size[0] ), _height ( Size[1] ), _numPixels ( Size[0] * Size[1] ) { }

  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::FullMatrix<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::FullMatrix<RealType> &Dest ) const {
    if ( Arg.size ( ) != _height )
      throw aol::Exception ( "Arguments size does not match patch height!", __FILE__, __LINE__ );

    if ( Dest.getNumRows ( ) != _numPixels || Dest.getNumCols ( ) != getNumArgs ( ) )
      throw aol::Exception ( "Destination dimensions do not fit patch size and parameters!", __FILE__, __LINE__ );

    for ( short k=0; k<_height ; ++k ) {
      for ( short xBlock=0; xBlock<_width ; ++xBlock ) {
        for ( short yBlock=0; yBlock<_height ; ++yBlock ) {
          const short x = _pos[0] + xBlock, y = _pos[1] + yBlock;
          const bool dkj1 = ( k == yBlock + 1 ), dkj = ( k == yBlock );
          Dest.set ( _mapper.getGlobalIndex ( xBlock, yBlock ), k, dkj1 * dxData ( x + ( ( yBlock + 1 < _height ) ? Arg[yBlock + 1] : 0 ), y + 1 ) - dkj * dxData ( x + Arg[yBlock], y ) );
        }
      }
    }
  }

  RealType dxData ( const RealType X, const short Y ) const {
    const short i = floor ( X ), j = Y;
    return _data.get ( i + 1, j ) - _data.get ( i, j );
  }

  int getGlobalIndex ( const short X, const short Y ) const {
    return _mapper.getGlobalIndex ( X, Y );
  }

  short getWidth ( ) const {
    return _width;
  }

  short getHeight ( ) const {
    return _height;
  }

  short getNumArgs ( ) const {
    return _height;
  }

  int getNumPixels ( ) const {
    return this->_numPixels;
  }
};


template <typename _RealType, typename _PictureType>
class BlockRegularityTargetFunctionalFixedCenter : public BlockRegularityTargetFunctional<_RealType, _PictureType> {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  const short _fixedLine;
  const RealType _epsilon;
public:
  BlockRegularityTargetFunctionalFixedCenter ( const PictureType &Data, const aol::Vec2<short> &Pos, const aol::Vec2<short> &Size, const short FixedLine,
                                               const RealType Epsilon = DEFAULT_EPSILON )
    : BlockRegularityTargetFunctional<RealType, PictureType> ( Data, Pos, Size ), _fixedLine ( FixedLine ), _epsilon ( Epsilon ) {
    if ( FixedLine < 0 || FixedLine >= this->_height )
      throw aol::Exception ( "Specified fixed line exceeds patch dimensions!", __FILE__, __LINE__ );
  }

  virtual ~BlockRegularityTargetFunctionalFixedCenter ( ) { }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( Arg.size ( ) != getNumArgs ( ) )
      throw aol::Exception ( "Arguments size does not match patch height!", __FILE__, __LINE__ );

    if ( Dest.size ( ) != this->_numPixels )
      throw aol::Exception ( "Destination vector does not match patch size!", __FILE__, __LINE__ );

    aol::Vector<RealType> shifts ( Arg );
    transformVariablesToCanonicalBasis ( shifts );

    for ( short xBlock=0; xBlock<this->_width ; ++xBlock ) {
      for ( short yBlock=0; yBlock<this->_height ; ++yBlock ) {
        Dest[this->_mapper.getGlobalIndex ( xBlock , yBlock )] = RegularizedL2Norm ( xBlock, yBlock, shifts );
      }
    }
  }

  RealType L2NormSqr ( const short X, const short Y, const aol::Vector<RealType> &Shifts ) const {
    RealType normSqr = 0;
    const short i = X, j = Y;
    const RealType sjp1 = ( ( j+1 < Shifts.size ( ) ) ? Shifts[j+1] : 0 ), sj = Shifts[j], fracShiftjp1 = ceil ( sjp1 ) - sjp1, fracShiftj = ceil ( sj ) - sj;
    aol::Vector<RealType> xl ( 4 );
    xl[0] = i;
    xl[1] = i + aol::Min ( fracShiftj, fracShiftjp1 );
    xl[2] = i + aol::Max ( fracShiftj, fracShiftjp1 );
    xl[3] = i+1;
    short h = 0;
    do {
      if ( xl[h] == xl[h+1] )
        xl.erase ( h+1 );
      else
        ++h;
    } while ( h < xl.size ( )-1 );
    for ( short l=0; l<xl.size ( )-1 ; ++l ) {
      RealType x = ( xl[l+1] + xl[l] ) / 2;
      RealType mlj = 0, nlj = 0, x0lj = 0, mljp1 = 0, nljp1 = 0, x0ljp1 = 0;
      setCoefficients ( x + sj, j, mlj, nlj, x0lj );
      setCoefficients ( x + sjp1, j+1, mljp1, nljp1, x0ljp1 );
      for ( short lp=0; lp<=1 ; ++lp ) {
        if ( mljp1 == mlj )
          normSqr += ( ( lp == 1 ) ? 1 : -1 ) * xl[l+lp] * pow ( nljp1 - nlj + mljp1 * ( sjp1 - x0lj - ( sj - x0lj ) ), 2 );
        else
          normSqr += ( ( lp == 1 ) ? 1 : -1 ) * pow ( mljp1 * ( xl[l+lp] + sjp1 - x0ljp1 ) + nljp1 - ( mlj * ( xl[l+lp] + sj - x0lj ) + nlj ), 3 ) / ( 3 * ( mljp1 - mlj ) );
      }
    }
    return normSqr;
  }

  RealType L2Norm ( const short X, const short Y, const aol::Vector<RealType> &Shifts ) const {
    return sqrt ( L2NormSqr ( X, Y, Shifts ) );
  }

  RealType RegularizedL2Norm ( const short X, const short Y, const aol::Vector<RealType> &Shifts ) const {
    return sqrt ( L2NormSqr ( X, Y, Shifts ) + pow ( _epsilon, 2 ) );
  }

  short getNumArgs ( ) const {
    return this->_height-1;
  }

  void transformVariablesToCanonicalBasis ( aol::Vector<RealType> &Shifts ) const {
    addCenterZeroShift ( Shifts );
    const aol::Vector<RealType> temp ( Shifts );
    for ( int k=0; k<_fixedLine ; ++k ) {
      Shifts[k] = 0;
      for ( int l=k; l<_fixedLine ; ++l )
        Shifts[k] += temp[l];
    }
    for ( int k=_fixedLine+1; k<this->_height ; ++k ) {
      Shifts[k] = 0;
      for ( int l=_fixedLine+1; l<=k ; ++l )
        Shifts[k] += temp[l];
    }
  }

private:
  virtual void setCoefficients ( const RealType X, const short Y, RealType &Slope, RealType &Intercept, RealType &X0 ) const {
    const short i = floor ( X ), j = Y;
    Slope = this->_data.get ( this->_pos[0] + i + 1, this->_pos[1] + j ) - this->_data.get ( this->_pos[0] + i, this-> _pos[1] + j );
    Intercept = this->_data.get ( this->_pos[0] + i, this->_pos[1] + j );
    X0 = i;
  }

  void addCenterZeroShift ( aol::Vector<RealType> &Shifts ) const {
    const aol::Vector<RealType> temp ( Shifts );
    Shifts.resize ( Shifts.size ( ) + 1 );
    for ( int k=0; k<_fixedLine ; ++k )
      Shifts[k] = temp[k];
    Shifts[_fixedLine] = 0;
    for ( int k=_fixedLine+1; k<this->_height ; ++k )
      Shifts[k] = temp[k-1];
  }
};


template <typename _RealType, typename _PictureType>
class BlockRegularityTargetJacobianFixedCenter : public BlockRegularityTargetJacobian<_RealType, _PictureType> {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  const short _fixedLine;
  const RealType _epsilon;
public:
  BlockRegularityTargetJacobianFixedCenter ( const PictureType &Data, const aol::Vec2<short> &Pos, const aol::Vec2<short> &Size, const short FixedLine,
                                             const RealType Epsilon=DEFAULT_EPSILON )
    : BlockRegularityTargetJacobian<RealType, PictureType> ( Data, Pos, Size ), _fixedLine ( FixedLine ), _epsilon ( Epsilon ) {
    if ( FixedLine < 0 || FixedLine >= this->_height )
      throw aol::Exception ( "Specified fixed line exceeds patch dimensions!", __FILE__, __LINE__ );
  }

  virtual ~BlockRegularityTargetJacobianFixedCenter ( ) { }

  void apply ( const aol::Vector<RealType> &Arg, aol::FullMatrix<RealType> &Dest ) const {
    if ( Arg.size ( ) != getNumArgs ( ) )
      throw aol::Exception ( "Arguments size does not match patch height!", __FILE__, __LINE__ );

    if ( Dest.getNumRows ( ) != this->_numPixels || Dest.getNumCols ( ) != this->getNumArgs ( ) )
      throw aol::Exception ( "Destination dimensions do not fit patch size and parameters!", __FILE__, __LINE__ );

    aol::Vector<RealType> shifts ( Arg );
    transformVariablesToCanonicalBasis ( shifts );

    for ( short xBlock=0; xBlock<this->_width ; ++xBlock ) {
      for ( short yBlock=0; yBlock<this->_height ; ++yBlock ) {
        for ( short k=0; k<_fixedLine ; ++k )
          Dest.set ( this->_mapper.getGlobalIndex ( xBlock, yBlock ), k, DL2NormSqr ( k, xBlock, yBlock, shifts ) / ( 2 * RegularizedL2Norm ( xBlock, yBlock, shifts ) ) );
        for ( short k=_fixedLine+1; k<this->_height ; ++k )
          Dest.set ( this->_mapper.getGlobalIndex ( xBlock, yBlock ), k-1, DL2NormSqr ( k, xBlock, yBlock, shifts ) / ( 2 * RegularizedL2Norm ( xBlock, yBlock, shifts ) ) );
      }
    }
  }

  short getNumArgs ( ) const {
    return this->_height-1;
  }

  void transformVariablesToCanonicalBasis ( aol::Vector<RealType> &Shifts ) const {
    addCenterZeroShift ( Shifts );
    const aol::Vector<RealType> temp ( Shifts );
    for ( int k=0; k<_fixedLine ; ++k ) {
      Shifts[k] = 0;
      for ( int l=k; l<_fixedLine ; ++l )
        Shifts[k] += temp[l];
    }
    for ( int k=_fixedLine+1; k<this->_height ; ++k ) {
      Shifts[k] = 0;
      for ( int l=_fixedLine+1; l<=k ; ++l )
        Shifts[k] += temp[l];
    }
  }

private:
  virtual void setCoefficients ( const RealType X, const short Y, RealType &Slope, RealType &Intercept, RealType &X0 ) const {
    const short i = floor ( X ), j = Y;
    Slope = this->_data.get ( this->_pos[0] + i + 1, this->_pos[1] + j ) - this->_data.get ( this->_pos[0] + i, this-> _pos[1] + j );
    Intercept = this->_data.get ( this->_pos[0] + i, this->_pos[1] + j );
    X0 = i;
  }

  void addCenterZeroShift ( aol::Vector<RealType> &Shifts ) const {
    const aol::Vector<RealType> temp ( Shifts );
    Shifts.resize ( Shifts.size ( ) + 1 );
    for ( int k=0; k<_fixedLine ; ++k )
      Shifts[k] = temp[k];
    Shifts[_fixedLine] = 0;
    for ( int k=_fixedLine+1; k<this->_height ; ++k )
      Shifts[k] = temp[k-1];
  }

  RealType L2NormSqr ( const short X, const short Y, const aol::Vector<RealType> &Shifts ) const {
    RealType normSqr = 0;
    const short i = X, j = Y;
    const RealType sjp1 = ( ( j+1 < Shifts.size ( ) ) ? Shifts[j+1] : 0 ), sj = Shifts[j], fracShiftjp1 = ceil ( sjp1 ) - sjp1, fracShiftj = ceil ( sj ) - sj;
    aol::Vector<RealType> xl ( 4 );
    xl[0] = i;
    xl[1] = i + aol::Min ( fracShiftj, fracShiftjp1 );
    xl[2] = i + aol::Max ( fracShiftj, fracShiftjp1 );
    xl[3] = i+1;
    short h = 0;
    do {
      if ( xl[h] == xl[h+1] )
        xl.erase ( h+1 );
      else
        ++h;
    } while ( h < xl.size ( )-1 );
    for ( short l=0; l<xl.size ( )-1 ; ++l ) {
      RealType x = ( xl[l+1] + xl[l] ) / 2;
      RealType mlj = 0, nlj = 0, x0lj = 0, mljp1 = 0, nljp1 = 0, x0ljp1 = 0;
      setCoefficients ( x + sj, j, mlj, nlj, x0lj );
      setCoefficients ( x + sjp1, j+1, mljp1, nljp1, x0ljp1 );
      for ( short lp=0; lp<=1 ; ++lp ) {
        if ( mljp1 == mlj )
          normSqr += ( ( lp == 1 ) ? 1 : -1 ) * xl[l+lp] * pow ( nljp1 - nlj + mljp1 * ( sjp1 - x0lj - ( sj - x0lj ) ), 2 );
        else
          normSqr += ( ( lp == 1 ) ? 1 : -1 ) * pow ( mljp1 * ( xl[l+lp] + sjp1 - x0ljp1 ) + nljp1 - ( mlj * ( xl[l+lp] + sj - x0lj ) + nlj ), 3 ) / ( 3 * ( mljp1 - mlj ) );
      }
    }
    return normSqr;
  }

  RealType L2Norm ( const short X, const short Y, const aol::Vector<RealType> &Shifts ) const {
    return sqrt ( L2NormSqr ( X, Y, Shifts ) );
  }

  RealType RegularizedL2Norm ( const short X, const short Y, const aol::Vector<RealType> &Shifts ) const {
    return sqrt ( L2NormSqr ( X, Y, Shifts ) + pow( _epsilon, 2 ) );
  }

  RealType DL2NormSqr ( const short Direction, const short X, const short Y, const aol::Vector<RealType> &Shifts ) const {
    RealType dNormSqr = 0;
    const short i = X, j = Y, k = Direction;
    const bool dkjp1 = ( ( k<_fixedLine ) ? ( k >= j+1 ) : ( k <= j+1 ) ), dkj = ( ( k<_fixedLine ) ? ( k >= j ) : ( k <= j ) );
    RealType sjp1 = ( ( j+1 < this->_height ) ? Shifts[j+1] : 0 ), sj = Shifts[j], fracShiftjp1 = ceil ( sjp1 ) - sjp1, fracShiftj = ceil ( sj ) - sj;
    aol::Vector<RealType> xl ( 4 );
    xl[0] = i;
    xl[1] = i + aol::Min ( fracShiftj, fracShiftjp1 );
    xl[2] = i + aol::Max ( fracShiftj, fracShiftjp1 );
    xl[3] = i+1;
    short h = 0;
    do {
      if ( xl[h] == xl[h+1] )
        xl.erase ( h+1 );
      else
        ++h;
    } while ( h < xl.size ( )-1 );
    for ( short l=0; l<xl.size ( )-1 ; ++l ) {
      RealType x = ( xl[l+1] + xl[l] ) / 2;
      RealType mlj = 0, nlj = 0, x0lj = 0, mljp1 = 0, nljp1 = 0, x0ljp1 = 0;
      setCoefficients ( x + sj, j, mlj, nlj, x0lj );
      setCoefficients ( x + sjp1, j+1, mljp1, nljp1, x0ljp1 );
      RealType dFl = ( mljp1 * dkjp1 - mlj * dkj );
      RealType Fl = 0;
      for ( short lp=0; lp<=1 ; ++lp ) {
        Fl += ( ( lp == 1 ) ? 1 : -1 ) * xl[l+lp] * ( mljp1 * ( 0.5 * xl[l+lp] + sjp1 - x0ljp1 ) + nljp1 - ( mlj * ( 0.5 * xl[l+lp] + sj - x0lj ) + nlj ) );
      }
      dNormSqr += 2 * dFl * Fl;
    }
    return dNormSqr;
  }
};


template <typename _RealType, typename _PictureType>
class ScalarBlockRegularityTargetFunctionalFixedCenter : public aol::Op<aol::Vector<_RealType>, aol::Scalar<_RealType> > {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  BlockRegularityTargetFunctionalFixedCenter<RealType, PictureType> _func;
  const short _xBlock, _yBlock;
public:
  ScalarBlockRegularityTargetFunctionalFixedCenter ( const PictureType &Data, const aol::Vec2<short> &Pos, const aol::Vec2<short> &Size, const short FixedLine, const short XBlock, const short YBlock )
    : _func ( Data, Pos, Size, FixedLine ), _xBlock ( XBlock ), _yBlock ( YBlock ) { }

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    aol::Vector<RealType> tmpDest ( _func.getNumPixels ( ) );
    _func.apply ( Arg, tmpDest );

    Dest = tmpDest[_func.getGlobalIndex ( _xBlock, _yBlock )];
  }
};


template <typename _RealType, typename _PictureType>
class ScalarBlockRegularityTargetJacobianFixedCenter : public aol::Op<aol::Vector<_RealType> > {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  BlockRegularityTargetJacobianFixedCenter<RealType, PictureType> _jacobian;
  const short _xBlock, _yBlock;
public:
  ScalarBlockRegularityTargetJacobianFixedCenter ( const PictureType &Data, const aol::Vec2<short> &Pos, const aol::Vec2<short> &Size, const short FixedLine, const short XBlock, const short YBlock )
    : _jacobian ( Data, Pos, Size, FixedLine ), _xBlock ( XBlock ), _yBlock ( YBlock ) { }

  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::Vector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( Dest.size ( ) != _jacobian.getNumArgs ( ) )
      throw aol::Exception ( "Destination vector does not match arguments size!", __FILE__, __LINE__ );

    aol::FullMatrix<RealType> tmpDest ( _jacobian.getNumPixels ( ), _jacobian.getNumArgs ( ) );
    _jacobian.apply ( Arg, tmpDest );

    for ( short k=0; k<_jacobian.getNumArgs ( ) ; ++k )
      Dest[k] = tmpDest.get ( _jacobian.getGlobalIndex ( _xBlock, _yBlock ), k );
  }
};


template <typename _RealType, typename _MatrixType, typename _PictureType, typename _LinearRegressionType>
class BlockRegularizerFixedCenter {
  typedef _RealType RealType;
  typedef _MatrixType MatrixType;
  typedef _PictureType PictureType;
  typedef _LinearRegressionType LinearRegressionType;
  typedef qc::ScalarArray<RealType, qc::QC_2D> BlockType;
protected:
  PictureType _data;
  const int _maxIterations;
  const RealType _mu0, _rho0, _rho1, _epsDeltaX, _epsF, _epsGradRealTargetFunc;
  const bool _verbose;
public:
  BlockRegularizerFixedCenter (  const PictureType &Data,
                                 const int MaxIterations = 50,
                                 const RealType Mu0 = 1,
                                 const RealType Rho0 = 0.2,
                                 const RealType Rho1 = 0.8,
                                 const RealType EpsDeltaX = 1e-6,
                                 const RealType EpsF = 1e-6,
                                 const RealType EpsGradRealTargetFunc = 1e-6,
                                 const bool Verbose = true )
    : _data ( Data ), _maxIterations ( MaxIterations ), _mu0 ( Mu0 ), _rho0 ( Rho0 ), _rho1 ( Rho1 ), _epsDeltaX ( EpsDeltaX ), _epsF ( EpsF ), _epsGradRealTargetFunc ( EpsGradRealTargetFunc ), _verbose ( Verbose ) { }

  virtual ~BlockRegularizerFixedCenter ( ) { }

  virtual void apply ( BlockType &Block, const aol::Vec2<short> &Pos, const short CenterLineIdx ) const {
    BlockRegularityTargetFunctionalFixedCenter<RealType, PictureType> func ( _data, Pos, aol::Vec2<short> ( Block.getNumX ( ), Block.getNumY ( ) ), CenterLineIdx );
    BlockRegularityTargetJacobianFixedCenter<RealType, PictureType> jacobian ( _data, Pos, aol::Vec2<short> ( Block.getNumX ( ), Block.getNumY ( ) ), CenterLineIdx );

    aol::Vector<RealType> dest ( Block.getNumY ( ) - 1 );
    optimizeShifts ( dest, func, jacobian, Block.size ( ), Pos, aol::Vec2<short> ( Block.getNumX ( ), Block.getNumY ( ) ) );
    applyShifts ( Block, Pos, dest );
  }

  virtual void apply ( ArbitraryAnchorBlockCollection<RealType, PictureType> &Blocks ) const {
    BlockType block ( Blocks.getBlockSize ( ), Blocks.getBlockSize ( ) );
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for ( short xAnchor=Blocks.getX0 ( ); xAnchor<Blocks.getXEnd ( ) ; ++xAnchor ) {
      for ( short yAnchor=Blocks.getY0 ( ); yAnchor<Blocks.getYEnd ( ) ; ++yAnchor ) {
        apply ( block, aol::Vec2<short> ( xAnchor - Blocks.getBlockAnchor ( )[0], yAnchor - Blocks.getBlockAnchor ( )[1] ), Blocks.getBlockAnchor ( )[1] );
        Blocks.copyFrom ( xAnchor, yAnchor, block );
      }
    }
  }

  virtual void applyShifts ( BlockType &Block, const aol::Vec2<short> &Pos, const aol::Vector<RealType> &Shifts ) const {
    for ( short xBlock=0; xBlock<Block.getNumX ( ) ; ++xBlock )
      for ( short yBlock=0; yBlock<Block.getNumY ( ) ; ++yBlock )
        Block.set ( xBlock, yBlock, _data.get ( Pos[0] + xBlock + round ( Shifts[yBlock] ), Pos[1] + yBlock ) );
  }
  
  void regularizeRows ( const short WindowSize ) {
    const short WindowOffset = ( WindowSize - 1 ) / 2;
    aol::Vector<RealType> est ( _data.getNumX ( ) );
    for ( int y=0; y<_data.getNumY ( ) ; ++y ) {
      for ( int x=WindowOffset; x<_data.getNumX ( ) - WindowOffset ; ++x ) {
        RealType avg = 0;
        for ( int dx=-WindowOffset; dx<=WindowOffset ; ++dx )
          avg += _data.get ( x + dx, y );
        est[x] = avg / WindowSize;
      }
      for ( int x=WindowOffset; x<_data.getNumX ( ) - WindowOffset ; ++x )
        _data.set ( y, x, est[x] );
    }
  }

protected:
  void optimizeShifts ( aol::Vector<RealType> &Dest,
                        BlockRegularityTargetFunctionalFixedCenter<RealType, PictureType> &Func,
                        BlockRegularityTargetJacobianFixedCenter<RealType, PictureType> &Jacobian,
                        const int DimRangeF, const aol::Vec2<short> &Pos, const aol::Vec2<short> &Size ) const {
    LevenbergMarquardtAlgorithm<RealType, MatrixType, LinearRegressionType> levenbergMarquardtAlg ( DimRangeF, Func, Jacobian,
                                                                                                    _maxIterations, _mu0, _rho0, _rho1, _epsDeltaX, _epsF, _epsGradRealTargetFunc, _verbose );
    aol::Vector<RealType> arg ( Size[1] - 1 );
    levenbergMarquardtAlg.apply ( arg, Dest );
    Func.transformVariablesToCanonicalBasis ( Dest );

    if ( _verbose && Dest.norm ( ) > Dest.size ( ) * Size[1] ) {
      std::cerr << "Warning: Norm of shifts is greater than the patch size times the number of shifts." << std::endl;
      std::cerr << "Pos=" << Pos << "; Size=" << Size << std::endl;
    }
  }
};


template <typename _RealType, typename _PictureType>
class BlockRegularityTargetFunctionalFixedCenterPenalized : public BlockRegularityTargetFunctionalFixedCenter<_RealType, _PictureType> {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  const RealType _maxVal;
public:
  BlockRegularityTargetFunctionalFixedCenterPenalized ( const PictureType &Data, const aol::Vec2<short> &Pos, const aol::Vec2<short> &Size, const short FixedLine, const RealType MaxVal )
    : BlockRegularityTargetFunctionalFixedCenter<RealType, PictureType> ( Data, Pos, Size, FixedLine ), _maxVal ( MaxVal ) { }
  
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( Dest.size ( ) != this->_numPixels + this->getHeight ( ) )
      throw aol::Exception ( "Destination vector does not match patch size plus patch height (used for penalty)!", __FILE__, __LINE__ );
    
    aol::Vector<RealType> tmpDest ( this->getNumPixels ( ) );
    BlockRegularityTargetFunctionalFixedCenter<RealType, PictureType>::applyAdd ( Arg, tmpDest );
    for ( short i=0; i<this->getNumPixels ( ) ; ++i )
      Dest[i] = tmpDest[i];
    
    // Add penalty function entries
    aol::Vector<RealType> shifts ( Arg );
    this->transformVariablesToCanonicalBasis ( shifts );
    RealType lambda = _maxVal / ( this->_width * sqrt ( this->_width ) );
    for ( short i=0; i<this->getHeight ( ) ; ++i )
      Dest[i + this->getNumPixels ( )] = ( abs ( shifts[i] ) > this->_width ) ? lambda * pow ( abs ( shifts[i] ) - this->_width, 2 ) : 0;
  }
};


template <typename _RealType, typename _PictureType>
class BlockRegularityTargetJacobianFixedCenterPenalized : public BlockRegularityTargetJacobianFixedCenter<_RealType, _PictureType> {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  const RealType _maxVal;
public:
  BlockRegularityTargetJacobianFixedCenterPenalized ( const PictureType &Data, const aol::Vec2<short> &Pos, const aol::Vec2<short> &Size, const short FixedLine, const RealType MaxVal )
    : BlockRegularityTargetJacobianFixedCenter<RealType, PictureType> ( Data, Pos, Size, FixedLine ), _maxVal ( MaxVal ) { }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::FullMatrix<RealType> &Dest ) const {
    if ( Dest.getNumRows ( ) != this->_numPixels + this->getHeight ( ) || Dest.getNumCols ( ) != this->getNumArgs ( ) )
      throw aol::Exception ( "Destination dimensions do not fit patch size plus patch height (used for penalty) and parameters!", __FILE__, __LINE__ );
    
    Dest.setZero ( );
    aol::FullMatrix<RealType> tmpDest ( this->_numPixels, this->getNumArgs ( ) );
    BlockRegularityTargetJacobianFixedCenter<RealType, PictureType>::apply ( Arg, tmpDest );
    for ( short i=0; i<this->_numPixels ; ++i )
      for ( short k=0; k<this->getNumArgs ( ) ; ++k )
        Dest.set ( i, k, tmpDest.get ( i, k ) );
    
    RealType lambda = _maxVal / ( this->_width * sqrt ( this->_width ) );
    aol::Vector<RealType> shifts ( Arg );
    this->transformVariablesToCanonicalBasis ( shifts );
    for ( short i=0; i<this->_height ; ++i ) {
      for ( short k=i; k<this->_fixedLine ; ++k )
        Dest.set ( i + this->_numPixels, k, ( abs ( shifts[i] ) > this->_width ) ? ( ( shifts[i] > 0 ) ? 1 : -1 ) * lambda * 2 * ( abs ( shifts[i] ) - this->_width ) : 0 );
      for ( short k=this->_fixedLine+1; k<=i ; ++k )
        Dest.set ( i + this->_numPixels, k-1, ( abs ( shifts[i] ) > this->_width ) ? ( ( shifts[i] > 0 ) ? 1 : -1 ) * lambda * 2 * ( abs ( shifts[i] ) - this->_width ) : 0 );
    }
  }
};


template <typename _RealType, typename _PictureType>
class ScalarBlockRegularityTargetFunctionalFixedCenterPenalized : public aol::Op<aol::Vector<_RealType>, aol::Scalar<_RealType> > {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  BlockRegularityTargetFunctionalFixedCenterPenalized<RealType, PictureType> _func;
  const short _penalizedComponent;
public:
  ScalarBlockRegularityTargetFunctionalFixedCenterPenalized ( const PictureType &Data, const aol::Vec2<short> &Pos, const aol::Vec2<short> &Size, const short FixedLine, const RealType MaxVal,
                                                              const short PenalizedComponent )
    : _func ( Data, Pos, Size, FixedLine, MaxVal ), _penalizedComponent ( PenalizedComponent ) { }
  
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    aol::Vector<RealType> tmpDest ( _func.getNumPixels ( ) + _func.getHeight ( ) );
    _func.apply ( Arg, tmpDest );
    
    Dest = tmpDest[_func.getNumPixels ( ) + _penalizedComponent];
  }
};


template <typename _RealType, typename _PictureType>
class ScalarBlockRegularityTargetJacobianFixedCenterPenalized : public aol::Op<aol::Vector<_RealType> > {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  BlockRegularityTargetJacobianFixedCenterPenalized<RealType, PictureType> _jacobian;
  const short _penalizedComponent;
public:
  ScalarBlockRegularityTargetJacobianFixedCenterPenalized ( const PictureType &Data, const aol::Vec2<short> &Pos, const aol::Vec2<short> &Size, const short FixedLine, const RealType MaxVal,
                                                            const short PenalizedComponent )
    : _jacobian ( Data, Pos, Size, FixedLine, MaxVal ), _penalizedComponent ( PenalizedComponent ) { }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::Vector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( Dest.size ( ) != _jacobian.getNumArgs ( ) )
      throw aol::Exception ( "Destination vector does not match arguments size!", __FILE__, __LINE__ );
    
    aol::FullMatrix<RealType> tmpDest ( _jacobian.getNumPixels ( ) + _jacobian.getHeight ( ), _jacobian.getNumArgs ( ) );
    _jacobian.apply ( Arg, tmpDest );
    
    for ( short k=0; k<_jacobian.getNumArgs ( ) ; ++k )
      Dest[k] = tmpDest.get ( _jacobian.getNumPixels ( ) + _penalizedComponent, k );
  }
};


template <typename _RealType, typename _MatrixType, typename _PictureType, typename _LinearRegressionType>
class BlockRegularizerFixedCenterPenalized : public BlockRegularizerFixedCenter<_RealType, _MatrixType, _PictureType, _LinearRegressionType> {
  typedef _RealType RealType;
  typedef _MatrixType MatrixType;
  typedef _PictureType PictureType;
  typedef _LinearRegressionType LinearRegressionType;
  typedef qc::ScalarArray<RealType, qc::QC_2D> BlockType;
protected:
  const RealType _maxVal;
public:
  BlockRegularizerFixedCenterPenalized ( const PictureType &Data,
                                         const RealType MaxVal,
                                         const int MaxIterations = 50,
                                         const RealType Mu0 = 1,
                                         const RealType Rho0 = 0.2,
                                         const RealType Rho1 = 0.8,
                                         const RealType EpsDeltaX = 1e-6,
                                         const RealType EpsF = 1e-6,
                                         const RealType EpsGradRealTargetFunc = 1e-6,
                                         const bool Verbose = true )
    : BlockRegularizerFixedCenter<RealType, MatrixType, PictureType, LinearRegressionType> ( Data, MaxIterations, Mu0, Rho0, Rho1, EpsDeltaX, EpsF, EpsGradRealTargetFunc, Verbose ),
      _maxVal ( MaxVal ) { }
  
  void apply ( BlockType &Block, const aol::Vec2<short> &Pos, const short CenterLineIndex ) const {
    aol::Vector<RealType> shifts ( Block.getNumY ( ) - 1 );
    setOptimalShifts ( shifts, Pos, aol::Vec2<short> ( Block.getNumX ( ), Block.getNumY ( ) ), CenterLineIndex );
    this->applyShifts ( Block, Pos, shifts );
  }
  
  void apply ( ArbitraryAnchorBlockCollection<RealType, PictureType> &Blocks ) const {
    BlockType block ( Blocks.getBlockSize ( ), Blocks.getBlockSize ( ) );
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for ( short xAnchor=Blocks.getX0 ( ); xAnchor<Blocks.getX0 ( ) + Blocks.getNumXEff ( ) ; ++xAnchor ) {
      for ( short yAnchor=Blocks.getY0 ( ); yAnchor<Blocks.getY0 ( ) + Blocks.getNumYEff ( ) ; ++yAnchor ) {
        apply ( block, aol::Vec2<short> ( xAnchor - Blocks.getX0 ( ), yAnchor - Blocks.getY0 ( ) ), Blocks.getBlockAnchor ( )[1] );
        Blocks.copyFrom ( xAnchor, yAnchor, block );
      }
    }
  }
  
  void setOptimalShifts ( aol::Vector<RealType> &Dest, const aol::Vec2<short> &Pos, const aol::Vec2<short> &Size, const short FixedLine ) const {
    if ( Dest.size ( ) != Size[1] - 1 )
      throw aol::Exception ( "Destination vector does not match patch height!", __FILE__, __LINE__ );
    
    BlockRegularityTargetFunctionalFixedCenterPenalized<RealType, PictureType> func ( this->_data, Pos, Size, FixedLine, _maxVal );
    BlockRegularityTargetJacobianFixedCenterPenalized<RealType, PictureType> jacobian ( this->_data, Pos, Size, FixedLine, _maxVal );
    
    this->optimizeShifts ( Dest, func, jacobian, Size[0] * Size[1] + Size[1], Pos, Size );
  }
  
  const PictureType& getDataConstRef ( ) const {
    return this->_data;
  }
};


template <typename _RealType, typename _MatrixType, typename _PictureType, typename _LinearRegressionType>
class MultiBlockRegularizer {
  typedef _RealType RealType;
  typedef _MatrixType MatrixType;
  typedef _PictureType PictureType;
  typedef qc::ScalarArray<RealType, qc::QC_2D> BlockType;
  typedef _LinearRegressionType LinearRegressionType;
protected:
  BlockRegularizerFixedCenterPenalized<RealType, MatrixType, PictureType, LinearRegressionType> _singleBlockRegularizer;
  const PictureType &_data;
  bool _verbose;
  qc::FastILexMapper<qc::QC_2D> _mapper;
public:
  MultiBlockRegularizer ( const PictureType &Data,
                          const RealType MaxVal,
                          const int MaxIterations = 20,
                          const RealType Mu0 = 1,
                          const RealType Rho0 = 0.2,
                          const RealType Rho1 = 0.8,
                          const RealType EpsDeltaX = 1e-6,
                          const RealType EpsF = 1e-6,
                          const RealType EpsGradRealTargetFunc = 1e-6,
                          const bool Verbose = true )
    : _singleBlockRegularizer ( Data, MaxVal, MaxIterations, Mu0, Rho0, Rho1, EpsDeltaX, EpsF, EpsGradRealTargetFunc, false ),
      _data ( Data ), _verbose ( Verbose ) { }
  
  void apply ( ArbitraryAnchorBlockCollection<RealType, PictureType> &Dest, aol::MultiVector<RealType> &Shifts,
               const bool RegularizeShifts = true, const bool PreRegularizeRows = true,
               const bool SkipBorderBlocks = false ) {
    if ( PreRegularizeRows )
      _singleBlockRegularizer.regularizeRows ( Dest.getBlockSize ( ) );
    
    // Initialize shift container
    _mapper.resize ( Dest.getNumX ( ), Dest.getNumY ( ) );
    Shifts.reallocate ( Dest.size ( ), Dest.getBlockSize ( ) );
    
    // Calculate individual shifts for block regularization
    int numIt = 0;
    const aol::Vec2<short> size ( Dest.getBlockSize ( ), Dest.getBlockSize ( ) );
    const short centerLineIdx = Dest.getBlockAnchor ( )[1];
    const short y0 = ( SkipBorderBlocks ) ? Dest.getY0 ( ) + 1 : Dest.getY0 ( );
    const short yEnd = ( SkipBorderBlocks ) ? Dest.getYEnd ( ) - 1 : Dest.getYEnd ( );
    aol::Vec2<short> position;
    #ifdef _OPENMP
    #pragma omp parallel for firstprivate ( numIt )
    #endif
    for ( short xAnchor=Dest.getX0 ( ); xAnchor<Dest.getX0 ( ) + Dest.getNumXEff ( ) ; ++xAnchor ) {
      for ( short yAnchor=y0; yAnchor<yEnd; ++yAnchor ) {
        aol::Vector<RealType> shifts ( Dest.getBlockSize ( ) - 1  );
        position.set ( xAnchor - Dest.getBlockAnchor ( )[0], yAnchor - Dest.getBlockAnchor ( )[1] );
        _singleBlockRegularizer.setOptimalShifts ( shifts, position, size, centerLineIdx );
        Shifts[_mapper.getGlobalIndex ( xAnchor, yAnchor )] = shifts;
        numIt++;
        if ( _verbose )
          std::cerr << "Regularized " << numIt << "/" << ( Dest.getNumXEff ( ) * Dest.getNumYEff ( ) ) << " blocks." << std::endl;
      }
    }
    
    applyShifts ( Dest, Shifts );
    if ( RegularizeShifts ) {
      regularizeShifts ( Dest, Shifts );
      applyShifts ( Dest, Shifts );
    }
  }
  
  void setQuietMode ( const bool Verbose = false ) {
    _verbose = Verbose;
  }
  
private:
  void regularizeShifts ( ArbitraryAnchorBlockCollection<RealType, PictureType> &Dest, aol::MultiVector<RealType> &Shifts ) const {
    aol::RandomAccessContainer<aol::Vector<RealType> > regularizedShifts ( Dest.getNumX ( ) * Dest.getNumY ( ) );
    BlockType compBlock;
    for ( short xAnchor=Dest.getX0 ( ); xAnchor<Dest.getXEnd ( ) ; ++xAnchor ) {
      for ( short yAnchor=Dest.getX0 ( ); yAnchor<Dest.getYEnd ( ) ; ++yAnchor ) {
        regularizedShifts[_mapper.getGlobalIndex ( xAnchor, yAnchor )].resize ( Dest.getBlockSize ( ) );
        const short blockOffset = Dest.getBlockSize ( ) / 2;
        RealType normFactor = 0;
        for ( short dx=-blockOffset; dx<=blockOffset; ++dx ) {
          if ( xAnchor+dx >= Dest.getX0 ( ) && xAnchor+dx < Dest.getX0 ( ) + Dest.getNumXEff ( ) ) {
            Dest.copyTo ( xAnchor+dx, yAnchor, compBlock );
            const RealType weight = std::exp ( -abs ( dx ) / compBlock.getNumY ( ) ) * std::exp ( compBlock.getMeanValue ( ) );
            regularizedShifts[_mapper.getGlobalIndex ( xAnchor, yAnchor )].addMultiple ( Shifts[_mapper.getGlobalIndex ( xAnchor+dx, yAnchor )], weight );
            normFactor += weight;
          }
        }
        regularizedShifts[_mapper.getGlobalIndex ( xAnchor, yAnchor )] /= normFactor;
      }
    }
    for ( int k=0; k<Shifts.numComponents ( ) ; ++k ) {
      if ( regularizedShifts[k].size ( ) > 0 )
        Shifts[k] = regularizedShifts[k];
    }
  }
  
  void applyShifts ( ArbitraryAnchorBlockCollection<RealType, PictureType> &Dest, aol::MultiVector<RealType> &Shifts ) {
    BlockType block ( Dest.getBlockSize ( ), Dest.getBlockSize ( ) );
    for ( short xAnchor=Dest.getX0 ( ); xAnchor<Dest.getXEnd ( ) ; ++xAnchor ) {
      for ( short yAnchor=Dest.getY0 ( ); yAnchor<Dest.getYEnd ( ) ; ++yAnchor ) {
        for ( short xBlock=0; xBlock<block.getNumX ( ) ; ++xBlock )
          for ( short yBlock=0; yBlock<block.getNumY ( ) ; ++yBlock )
            block.set ( xBlock, yBlock, Dest.get ( xAnchor - Dest.getX0 ( ) + xBlock + round ( Shifts[_mapper.getGlobalIndex ( xAnchor, yAnchor )][yBlock] ), yAnchor - Dest.getY0 ( ) + yBlock ) );
        Dest.copyFrom ( xAnchor, yAnchor, block );
      }
    }
  }
};


#endif /* BLOCKREGULARIZATION_H_ */
