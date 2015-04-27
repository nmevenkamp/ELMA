#ifndef NONLINEARREGRESSION_H_
#define NONLINEARREGRESSION_H_

#include <matrixInverse.h>
#include <linearRegression.h>
#include <projectors.h>
#include <statistics.h>

template <typename _RealType, typename _MatrixType, typename _LinearRegressionType>
class LevenbergMarquardtAlgorithm {
  typedef _RealType RealType;
  typedef _MatrixType MatrixType;
  typedef _LinearRegressionType LinearRegressionType;

protected:
  //! stores the evaluation of DF and is deleted in the destructor.
  mutable aol::DeleteFlagPointer<MatrixType> _jacobian;

  const int _dimRangeF;
  const aol::Op<aol::Vector<RealType> > &_F;
  const aol::Op<aol::Vector<RealType>, MatrixType> &_DF;
  const int _maxIterations;
  const RealType _mu0, _rho0, _rho1;
  const RealType _epsDeltaX, _epsF, _epsGradRealTargetFunc;
  const bool _verbose;

public:
  LevenbergMarquardtAlgorithm ( const int DimRangeF,
                                const aol::Op<aol::Vector<RealType> > &F,
                                const aol::Op<aol::Vector<RealType>, MatrixType> &DF,
                                const int MaxIterations = 50,
                                const RealType Mu0 = 1,
                                const RealType Rho0 = 0.2,
                                const RealType Rho1 = 0.8,
                                const RealType EpsDeltaX = 1e-6,
                                const RealType EpsF = 1e-6,
                                const RealType EpsGradRealTargetFunc = 1e-6,
                                const bool Verbose = true )
    : _jacobian ( NULL ),
      _dimRangeF ( DimRangeF ),
      _F ( F ),
      _DF ( DF ),
      _maxIterations ( MaxIterations ),
      _mu0 ( Mu0 ),
      _rho0 ( Rho0 ), _rho1 ( Rho1 ),
      _epsDeltaX ( EpsDeltaX ), _epsF ( EpsF ), _epsGradRealTargetFunc ( EpsGradRealTargetFunc ),
      _verbose ( Verbose ) { }

  virtual ~LevenbergMarquardtAlgorithm ( ) { }

  void apply( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( _jacobian.get ( ) == NULL )
      _jacobian.reset ( new MatrixType ( _dimRangeF, Arg.size ( ) ), true );

    if ( _verbose )
      std::cerr << "LevenbergMarquardt: Optimizing..." << std::endl;

    RealType mu = _mu0;
    int iteration = 0;
    aol::Vector<RealType> gradientRealValuedTargetFunction ( _jacobian->getNumCols ( ) );

    aol::Vector<RealType> x0 ( Arg ), x1 ( Arg );
    aol::Vector<RealType> f0 ( _dimRangeF ), f1 ( _dimRangeF );
    aol::Vector<RealType> direction ( Dest, aol::STRUCT_COPY );

    _F.apply ( x0, f0 );
    f1 = f0;

    if ( _verbose )
      std::cerr << "LevenbergMarquardt: Iteration " << iteration << ": x=" << x1 << "; res=" << f1.norm ( ) << std::endl;

    do {
      x0 = x1;
      f0 = f1;

      _DF.apply ( x0, *_jacobian );

      // Setup right-hand side
      aol::Vector<RealType> rhs ( _jacobian->getNumRows ( ) + _jacobian->getNumCols ( ) );
      for ( int i=0; i<f0.size ( ) ; ++i )
        rhs[i] = -f0[i];

      do {
        // Setup system matrix
        aol::FullMatrix<RealType> systemMatrix ( _jacobian->getNumRows ( ) + _jacobian->getNumCols ( ), _jacobian->getNumCols ( ) );
        for ( int i=0; i<_jacobian->getNumRows ( ) ; ++i ) {
          for ( int j=0; j<_jacobian->getNumCols ( ) ; ++j ) {
            systemMatrix.set ( i, j, _jacobian->get ( i, j ) );
            systemMatrix.set ( _jacobian->getNumRows ( ) + j, j, mu );
          }
        }

        LinearRegressionType linearRegression ( systemMatrix );
        linearRegression.apply ( rhs, direction );

        step ( x0, direction, x1 );
        _F.apply ( x1, f1 );

        // Test the step (and possibly correct mu)
        aol::Vector<RealType> linEstimate ( f0 );
        _jacobian->applyAdd ( direction, linEstimate );
        RealType linDiffNormSqr = f0.normSqr ( ) - linEstimate.normSqr ( );
        if ( linDiffNormSqr < 0 ) {
          if ( _verbose )
            std::cerr << "LevenbergMarquardt: Linear residual is negative. Breaking inner iteration." << std::endl;
          break;
        }
        RealType rhoMu = aol::NumberTrait<RealType>::Inf;
        rhoMu = ( f0.normSqr ( ) - f1.normSqr ( ) ) / linDiffNormSqr;

        if ( !aol::isFinite ( rhoMu ) ) {
          if ( _verbose )
            std::cerr << "LevenbergMarquardt: Quotient between linear and non-linear Residual is NaN or Inf. Breaking inner iteration." << std::endl;
          x1 = x0;
          _F.apply ( x1, f1 );
          break;
        }

        if ( rhoMu <= _rho0 ) {
          mu *= 10;
          if ( _verbose )
            std::cerr << "LevenbergMarquardt: Adapting damping parameter: Increasing to mu=" << mu << " and recalculating step." << std::endl;
        } else if ( _rho0 < rhoMu && rhoMu < _rho1 ) {
          if ( _verbose )
            std::cerr << "LevenbergMarquardt: Accepting step without changing the damping parameter." << std::endl;
          break;
        } else {
          mu /= 10;
          if ( _verbose )
            std::cerr << "LevenbergMarquardt: Adapting damping parameter: Decreasing to mu=" << mu << " and accepting step." << std::endl;
          break;
        }
      } while ( true );

      ++iteration;
      if ( _verbose )
        std::cerr << "LevenbergMarquardt: Iteration " << iteration << ": x=" << x1 << "; res=" << f1.norm ( ) << std::endl;

      // Calculate gradient of scalar target function phi(x) = 1/2 * ( F(x) ).normSqr ( )
      // Use: grad phi(x) = ( F'(x) )^T F(x)
      for ( int i=0; i<_jacobian->getNumCols ( ) ; ++i ) {
        gradientRealValuedTargetFunction[i] = 0;
        for ( int k=0; k<_jacobian->getNumRows ( ) ; ++k )
          gradientRealValuedTargetFunction[i] += _jacobian->get ( k, i ) * f1[k];
      }
    } while ( iteration < _maxIterations && f1.norm ( ) >= _epsF
              && gradientRealValuedTargetFunction.norm ( ) >= _epsGradRealTargetFunc
              && direction.norm ( ) / x0.norm ( ) >= _epsDeltaX );

    if ( _verbose ) {
      std::cerr << "LevenebrgMarquardt: Terminating iteration. Reason(s):";
      if ( iteration == _maxIterations )
        std::cerr << " maximum number of iterations reached";
      if ( f1.norm ( ) < _epsF )
        std::cerr << " residual below threshold";
      if ( gradientRealValuedTargetFunction.norm ( ) < _epsGradRealTargetFunc )
        std::cerr << " gradient of target function below threshold";
      if ( direction.norm ( ) / x0.norm ( ) < _epsDeltaX )
        std::cerr << " step size below threshold";
      std::cerr << "." << std::endl;
    }

    if ( _verbose ) {
      std::cerr << "LevenbergMarquardt: Optimization finished." << std::endl;
      std::cerr << "LevenbergMarquardt: Residual: " << f1.norm ( ) << std::endl;
      std::cerr << "LevenbergMarquardt: Optimal parameters: " << x1 << std::endl;
    }

    Dest = x1;
  }

  void applyAdd( const aol::Vector<RealType> &/*Arg*/, aol::Vector<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__);
  }

protected:
  virtual void step ( const aol::Vector<RealType> &X0, const aol::Vector<RealType> &Direction, aol::Vector<RealType> &X1 ) const {
    X1 = X0;
    X1 += Direction;
  }
};


// In each iteration solves unconstrained linearized regression Problem ||DF(x_k) s_k + F(x_k)||^2 + mu^2 ||x_k||^2 -> min.
// The projection of x_k+1 = x_k + s_k onto the convex space X = { x in R^n | A x <= B } is then used as the update
template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ProjectorType>
class ConstrainedLevenbergMarquardtAlgorithm : public LevenbergMarquardtAlgorithm<_RealType, _MatrixType, _LinearRegressionType> {
  typedef _RealType RealType;
  typedef aol::Vector<RealType> VectorType;
  typedef _MatrixType MatrixType;
  typedef _LinearRegressionType LinearRegressionType;
  typedef _ProjectorType ProjectorType;
protected:
  const ProjectorType _projector;
public:
  ConstrainedLevenbergMarquardtAlgorithm (  const int DimRangeF,
                                            const aol::Op<aol::Vector<RealType> > &F,
                                            const aol::Op<aol::Vector<RealType>, MatrixType> &DF,
                                            const ProjectorType &Projector,
                                            const int MaxIterations = 50,
                                            const RealType Mu0 = 1,
                                            const RealType Rho0 = 0.2,
                                            const RealType Rho1 = 0.8,
                                            const RealType EpsDeltaX = 1e-6,
                                            const RealType EpsF = 1e-6,
                                            const RealType EpsGradRealTargetFunc = 1e-6,
                                            const bool Verbose = true )
    : LevenbergMarquardtAlgorithm<RealType, MatrixType, LinearRegressionType> ( DimRangeF, F, DF, MaxIterations, Mu0, Rho0, Rho1, EpsDeltaX, EpsF, EpsGradRealTargetFunc, Verbose ),
      _projector ( Projector ) { }
protected:
  void step ( const aol::Vector<RealType> &X0, const aol::Vector<RealType> &Direction, aol::Vector<RealType> &X1 ) const {
    LevenbergMarquardtAlgorithm<RealType, MatrixType, LinearRegressionType>::step ( X0, Direction, X1 );
    _projector.apply ( X1, X1 );
  }
};


template <typename _RealType>
class SumOfSinesLeastSquaresEnergy : public aol::Op<aol::Vector<_RealType>, aol::Scalar<_RealType> > {
  typedef _RealType RealType;
protected:
  const std::vector<std::pair<RealType, RealType> > &_data;
  const int _numTerms;
public:
  SumOfSinesLeastSquaresEnergy ( const std::vector<std::pair<RealType, RealType> > &Data, const int NumTerms )
    : _data ( Data ), _numTerms ( NumTerms ) { }
  
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Scalar<RealType> &Dest ) const {
    if ( Arg.size ( ) != 3 * _numTerms )
      throw aol::Exception ( "Arguments size does not match 3 * number of sin terms!", __FILE__, __LINE__ );
    
    for ( int i=0; i<_data.size ( ) ; ++i )
      Dest += aol::Sqr<RealType> ( sumOfSines ( _data[i].first, Arg ) - _data[i].second );
    Dest *= 0.5;
  }
  
  RealType sumOfSines ( const RealType X, const aol::Vector<RealType> &Arg ) const {
    RealType res = 0;
    for ( short k=0; k<_numTerms ; ++k )
      res += Arg[3*k] * std::sin ( Arg[3*k+1] * X + Arg[3*k+2] );
    
    return res;
  }
  
  static RealType sumOfSines ( const RealType X, const aol::Vector<RealType> &Arg, const short NumTerms ) {
    RealType res = 0;
    for ( short k=0; k<NumTerms ; ++k )
      res += Arg[3*k] * std::sin ( Arg[3*k+1] * X + Arg[3*k+2] );
    
    return res;
  }
};

template <typename _RealType>
class SumOfSinesLeastSquaresEnergyDerivative : public aol::Op<aol::Vector<_RealType> > {
  typedef _RealType RealType;
protected:
  const std::vector<std::pair<RealType, RealType> > &_data;
  const int _numTerms;
public:
  SumOfSinesLeastSquaresEnergyDerivative ( const std::vector<std::pair<RealType, RealType> > &Data, const int NumTerms )
  : _data ( Data ), _numTerms ( NumTerms ) { }
  
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( Arg.size ( ) != 3 * _numTerms )
      throw aol::Exception ( "Arguments size does not match 3 * number of sin terms!", __FILE__, __LINE__ );
    
    if ( Dest.size ( ) != Arg.size ( ) )
      throw aol::Exception ( "Destination vector (gradient) size does not match argument size!", __FILE__, __LINE__ );
    
    for ( int i=0; i<_data.size ( ) ; ++i ) {
      const RealType x = _data[i].first;
      for ( short k=0; k<_numTerms ; ++k ) {
        Dest[3*k] += ( sumOfSines ( _data[i].first, Arg ) - _data[i].second ) * std::sin ( Arg[3*k+1] * x + Arg[3*k+2] );
        Dest[3*k+1] += ( sumOfSines ( _data[i].first, Arg ) - _data[i].second ) * Arg[3*k] * x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] );
        Dest[3*k+2] += ( sumOfSines ( _data[i].first, Arg ) - _data[i].second ) * Arg[3*k] * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] );
      }
    }
  }
  
  RealType sumOfSines ( const RealType X, const aol::Vector<RealType> &Arg ) const {
    RealType res = 0;
    for ( short k=0; k<_numTerms ; ++k )
      res += Arg[3*k] * std::sin ( Arg[3*k+1] * X + Arg[3*k+2] );
    
    return res;
  }
};

template <typename _RealType>
class SumOfSinesLeastSquaresEnergySecondDerivative : public aol::Op<aol::Vector<_RealType>, aol::FullMatrix<_RealType> > {
  typedef _RealType RealType;
protected:
  const std::vector<std::pair<RealType, RealType> > &_data;
  const int _numTerms;
public:
  SumOfSinesLeastSquaresEnergySecondDerivative ( const std::vector<std::pair<RealType, RealType> > &Data, const int NumTerms )
  : _data ( Data ), _numTerms ( NumTerms ) { }
  
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::FullMatrix<RealType> &Dest ) const {
    if ( Arg.size ( ) != 3 * _numTerms )
      throw aol::Exception ( "Arguments size does not match 3 * number of sin terms!", __FILE__, __LINE__ );
    
    if ( Dest.getNumRows ( ) != Arg.size ( ) || Dest.getNumCols ( ) != Arg.size ( ) )
      throw aol::Exception ( "Destination matrix (Hessian) row and/or column size does not match arguments size!", __FILE__, __LINE__ );
    
    for ( int i=0; i<_data.size ( ) ; ++i ) {
      const RealType x = _data[i].first;
      for ( short k=0; k<_numTerms ; ++k ) {
        for ( short l=0; l<_numTerms ; ++l ) {
          if ( k == l ) {
            Dest.ref ( 3*k, 3*k ) += aol::Sqr<RealType> ( std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) );
            Dest.ref ( 3*k, 3*k+1 ) += std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) * Arg[3*k] * x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] )
                                       + x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * ( sumOfSines ( _data[i].first, Arg ) - _data[i].second );
            Dest.ref ( 3*k, 3*k+2 ) += std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) * Arg[3*k] * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] )
                                       + std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * ( sumOfSines ( _data[i].first, Arg ) - _data[i].second );
            Dest.ref ( 3*k+1, 3*k ) += Arg[3*k] * x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * std::sin ( Arg[3*k+1] * x + Arg[3*k+2] )
                                       + x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * ( sumOfSines ( _data[i].first, Arg ) - _data[i].second );
            Dest.ref ( 3*k+1, 3*k+1 ) += aol::Sqr<RealType> ( Arg[3*k] * x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) )
                                         - Arg[3*k] * aol::Sqr<RealType> ( x ) * std::sin ( Arg[3*k+1] * x + Arg[3*k+2] )
                                         * ( sumOfSines ( _data[i].first, Arg ) - _data[i].second );
            Dest.ref ( 3*k+1, 3*k+2 ) += x * aol::Sqr<RealType> ( Arg[3*k] *  std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) )
                                         - Arg[3*k] * x * std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) * ( sumOfSines ( _data[i].first, Arg ) - _data[i].second );
            Dest.ref ( 3*k+2, 3*k ) += Arg[3*k] * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * std::sin ( Arg[3*k+1] * x + Arg[3*k+2] )
                                       + std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * ( sumOfSines ( _data[i].first, Arg ) - _data[i].second );
            Dest.ref ( 3*k+2, 3*k+1 ) += x * aol::Sqr<RealType> ( Arg[3*k] * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) )
                                         - Arg[3*k] * x * std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) * ( sumOfSines ( _data[i].first, Arg ) - _data[i].second );
            Dest.ref ( 3*k+2, 3*k+2 ) += aol::Sqr<RealType> ( Arg[3*k] * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) )
                                         - Arg[3*k] * std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) * ( sumOfSines ( _data[i].first, Arg ) - _data[i].second );
          } else {
            Dest.ref ( 3*k, 3*l ) += std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) * std::sin ( Arg[3*l+1] * x + Arg[3*l+2] );
            Dest.ref ( 3*k, 3*l+1 ) += std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) * Arg[3*l] * x * std::cos ( Arg[3*l+1] * x + Arg[3*l+2] );
            Dest.ref ( 3*k, 3*l+2 ) += std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) * Arg[3*l] * std::cos ( Arg[3*l+1] * x + Arg[3*l+2] );
            Dest.ref ( 3*k+1, 3*l ) += Arg[3*k] * x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * std::sin ( Arg[3*l+1] * x + Arg[3*l+2] );
            Dest.ref ( 3*k+1, 3*l+1 ) += Arg[3*k] * x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * Arg[3*l] * x * std::cos ( Arg[3*l+1] * x + Arg[3*l+2] );
            Dest.ref ( 3*k+1, 3*l+2 ) += Arg[3*k] * x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * Arg[3*l] * std::cos ( Arg[3*l+1] * x + Arg[3*l+2] );
            Dest.ref ( 3*k+2, 3*l ) += Arg[3*k] * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * std::sin ( Arg[3*l+1] * x + Arg[3*l+2] );
            Dest.ref ( 3*k+2, 3*l+1 ) += Arg[3*k] * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * Arg[3*l] * x * std::cos ( Arg[3*l+1] * x + Arg[3*l+2] );
            Dest.ref ( 3*k+2, 3*l+2 ) += Arg[3*k] * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) * Arg[3*l] * std::cos ( Arg[3*l+1] * x + Arg[3*l+2] );
          }
        }
      }
    }
  }
  
  RealType sumOfSines ( const RealType X, const aol::Vector<RealType> &Arg ) const {
    RealType res = 0;
    for ( short k=0; k<_numTerms ; ++k )
      res += Arg[3*k] * std::sin ( Arg[3*k+1] * X + Arg[3*k+2] );
    
    return res;
  }
};


template <typename _RealType>
class SumOfSinesTargetFunctional : public aol::Op<aol::Vector<_RealType> > {
  typedef _RealType RealType;
protected:
  const std::vector<std::pair<RealType, RealType> > &_data;
  const int _numTerms;
public:
  SumOfSinesTargetFunctional ( const std::vector<std::pair<RealType, RealType> > &Data, const int NumTerms )
    : _data ( Data ), _numTerms ( NumTerms ) { }
  
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( Arg.size ( ) != 3 * _numTerms )
      throw aol::Exception ( "Arguments size does not match 3 * number of sin terms!", __FILE__, __LINE__ );
    
    if ( Dest.size ( ) != _data.size ( ) )
      throw aol::Exception ( "Destination vector does not match number of data points!", __FILE__, __LINE__ );
    
    for ( int i=0; i<_data.size ( ) ; ++i )
      Dest[i] += sumOfSines ( _data[i].first, Arg ) - _data[i].second;
  }
  
  RealType sumOfSines ( const RealType X, const aol::Vector<RealType> &Arg ) const {
    RealType res = 0;
    for ( short k=0; k<_numTerms ; ++k )
      res += Arg[3*k] * std::sin ( Arg[3*k+1] * X + Arg[3*k+2] );
    
    return res;
  }
};

template <typename _RealType, typename _MatrixType>
class SumOfSinesTargetJacobian : public aol::Op<aol::Vector<_RealType>, _MatrixType> {
  typedef _RealType RealType;
protected:
  const std::vector<std::pair<RealType, RealType> > &_data;
  const int _numTerms;
public:
  SumOfSinesTargetJacobian ( const std::vector<std::pair<RealType, RealType> > &Data, const int NumTerms )
    : _data ( Data ), _numTerms ( NumTerms ) { }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::FullMatrix<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::FullMatrix<RealType> &Dest ) const {
    if ( Arg.size ( ) != 3 * _numTerms )
      throw aol::Exception ( "Arguments size does not match 3 * number of sin terms!", __FILE__, __LINE__ );
    
    if ( Dest.getNumRows ( ) != _data.size ( ) || Dest.getNumCols ( ) != 3 * _numTerms )
      throw aol::Exception ( "Destination dimensions do not fit number of data points and parameters!", __FILE__, __LINE__ );
  
    for ( int i=0; i<_data.size ( ) ; ++i ) {
       const RealType x = _data[i].first;
       for ( short k=0; k<_numTerms ; ++k ) {
        Dest.set ( i, 3*k, std::sin ( Arg[3*k+1] * x + Arg[3*k+2] ) );
        Dest.set ( i, 3*k+1, Arg[3*k] * x * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) );
        Dest.set ( i, 3*k+2, Arg[3*k] * std::cos ( Arg[3*k+1] * x + Arg[3*k+2] ) );
      }
    }
  }
};


template <typename _RealType>
class Gaussian1DTargetFunctional : public aol::Op<aol::Vector<_RealType> > {
  typedef _RealType RealType;
protected:
  const std::vector<std::pair<RealType, RealType> > &_data;
public:
  Gaussian1DTargetFunctional ( const std::vector<std::pair<RealType, RealType> > &Data )
    : _data ( Data ) { }
  
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( Arg.size ( ) != 3 )
      throw aol::Exception ( "Arguments size is inequal to two (mean, variance)!", __FILE__, __LINE__ );
    
    if ( Dest.size ( ) != _data.size ( ) )
      throw aol::Exception ( "Destination vector does not match number of data points!", __FILE__, __LINE__ );
    
    for ( int i=0; i<_data.size ( ) ; ++i )
      Dest[i] += Arg[2] * NormalDistribution<RealType>::PDF ( _data[i].first, Arg[0], Arg[1] ) - _data[i].second;
  }
};

template <typename _RealType, typename _MatrixType>
class Gaussian1DTargetJacobian : public aol::Op<aol::Vector<_RealType>, _MatrixType> {
  typedef _RealType RealType;
protected:
  const std::vector<std::pair<RealType, RealType> > &_data;
public:
  Gaussian1DTargetJacobian ( const std::vector<std::pair<RealType, RealType> > &Data )
    : _data ( Data ) { }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::FullMatrix<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::FullMatrix<RealType> &Dest ) const {
    if ( Arg.size ( ) != 3 )
      throw aol::Exception ( "Arguments size is inequal to two (mean, variance)!", __FILE__, __LINE__ );
    
    if ( Dest.getNumRows ( ) != _data.size ( ) || Dest.getNumCols ( ) != 3 )
      throw aol::Exception ( "Destination dimensions do not fit number of data points and parameters!", __FILE__, __LINE__ );
    
    for ( int i=0; i<_data.size ( ) ; ++i ) {
      const RealType x = _data[i].first;
      Dest.set ( i, 0, Arg[2] * ( ( x - Arg[0] ) / Arg[1] * NormalDistribution<RealType>::PDF ( _data[i].first, Arg[0], Arg[1] ) ) );
      Dest.set ( i, 1, Arg[2] * ( 0.5 * ( aol::Sqr<RealType> ( x - Arg[0] ) * NormalDistribution<RealType>::PDF ( _data[i].first, Arg[0], Arg[1] ) / aol::Sqr<RealType> ( Arg[1] )
                                          - NormalDistribution<RealType>::PDF ( _data[i].first, Arg[0], Arg[1] ) / Arg[1] ) ) );
      Dest.set ( i, 2, NormalDistribution<RealType>::PDF ( _data[i].first, Arg[0], Arg[1] ) );
      
    }
  }
};

#endif /* NONLINEARREGRESSION_H_ */
