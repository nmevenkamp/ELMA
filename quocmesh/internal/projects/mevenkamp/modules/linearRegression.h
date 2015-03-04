#ifndef LINEARREGRESSION_H_
#define LINEARREGRESSION_H_


#include <matrixInverse.h>
#include <preconditioner.h>


template <typename _RealType>
class LinearRegression {
  typedef _RealType RealType;
  typedef aol::Vector<_RealType> VectorType;
  typedef aol::FullMatrix<_RealType> MatrixType;
protected:
  MatrixType _A;
public:
  LinearRegression ( const MatrixType &SystemMatrix )
    : _A ( SystemMatrix ) { }

  virtual ~LinearRegression ( ) { }

  virtual void apply ( const VectorType &/*RHS*/, VectorType &/*X*/ ) = 0;
};


template <typename _RealType>
class LinearRegressionQR : public LinearRegression<_RealType> {
  typedef _RealType RealType;
  typedef aol::Vector<_RealType> VectorType;
  typedef aol::FullMatrix<_RealType> MatrixType;
protected:
  const aol::QRInverse<RealType> _qrInv;
public:
  LinearRegressionQR ( const MatrixType &SystemMatrix )
    : LinearRegression<RealType> ( SystemMatrix ), _qrInv ( SystemMatrix ) { }

  void apply ( const VectorType &RHS, VectorType &X ) {
    if ( RHS.size ( ) != this->_A.getNumRows ( ) )
      throw aol::Exception ( "Number of equations does not match between system matrix and right-hand side!", __FILE__, __LINE__ );

    if ( X.size ( ) != this->_A.getNumCols ( ) )
      throw aol::Exception ( "Vector size does not match number of columns of the system matrix!", __FILE__, __LINE__ );

    _qrInv ( RHS, X );

    if ( X.checkForNANsAndINFs() )
      throw aol::Exception ( "Linear regression failed! NaN / Inf entries detected in the solution!", __FILE__, __LINE__ );
  }
};


template <typename _RealType>
class LinearRegressionNormalEquations : public LinearRegression<_RealType> {
  typedef _RealType RealType;
  typedef aol::Vector<_RealType> VectorType;
  typedef aol::FullMatrix<_RealType> MatrixType;
protected:
  MatrixType _ATA;
  VectorType _ATb;
  aol::DiagonalMatrix<RealType> _P;
public:
  LinearRegressionNormalEquations ( const MatrixType &SystemMatrix )
    : LinearRegression<RealType> ( SystemMatrix ), _ATA ( SystemMatrix.getNumCols ( ), SystemMatrix.getNumCols ( ) ), _ATb ( SystemMatrix.getNumCols ( ) ),
      _P ( SystemMatrix.getNumCols ( ) ) {
    // Calculate _A^T * _A
    for ( int i=0; i<this->_A.getNumCols ( ) ; ++i ) {
      for ( int j=0; j<this->_A.getNumCols ( ) ; ++j ) {
        _ATA.set ( i, j, 0 );
        for ( int k=0; k<this->_A.getNumRows ( ) ; ++k )
          _ATA.add ( i, j, this->_A.get ( k, i ) * this->_A.get ( k, j ) );
      }
    }

    // Apply preconditioner to system matrix
    aol::DiagonalPreconditioner<aol::Vector<RealType> > preconditioner ( _ATA );
    _P.resize ( _ATA.getNumRows ( ) );
    preconditioner.getPreconditionerMatrix ( _P );
    for ( int i=0; i<_ATA.getNumRows ( ) ; ++i ) {
      for ( int j=0; j<_ATA.getNumCols ( ) ; ++j )
        _ATA.set ( i, j, _ATA.get ( i, j ) * _P.get ( i, i ) );
    }
  }

  void apply ( const VectorType &RHS, VectorType &X ) {
    if ( RHS.size ( ) != this->_A.getNumRows ( ) )
      throw aol::Exception ( "Number of equations does not match between system matrix and right-hand side!", __FILE__, __LINE__ );

    if ( X.size ( ) != this->_A.getNumCols ( ) )
      throw aol::Exception ( "Vector size does not match number of columns of the system matrix!", __FILE__, __LINE__ );

    // Calculate _A^T * _b and apply preconditioner
    _ATb.setZero ( );
    this->_A.applyAddTranspose ( RHS, _ATb );
    this->_P.apply ( _ATb, _ATb );

    // Solve P ATA X = P ATb using QR decomposition
    aol::QRInverse<RealType> qrInverse ( _ATA );
    qrInverse.apply ( _ATb, X );

    if ( X.checkForNANsAndINFs() )
      throw aol::Exception ( "Linear regression failed! NaN / Inf entries detected in the solution!", __FILE__, __LINE__ );
  }
};


#endif /* LINEARREGRESSION_H_ */
