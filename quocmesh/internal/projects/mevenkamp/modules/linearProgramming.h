#ifndef LINEARPROGRAMMING_H_
#define LINEARPROGRAMMING_H_


template <typename _RealType, typename _VectorType, typename _MatrixType>
class Projector {
  typedef _RealType RealType;
  typedef _VectorType VectorType;
  typedef _MatrixType MatrixType;
public:
  Projector ( ) { }

  virtual ~Projector ( ) { }

  virtual void apply ( const VectorType &/*Arg*/, VectorType &/*Dest*/ ) const = 0;

  virtual bool isFeasible ( const VectorType &/*X*/ ) const = 0;
};


// Given a constraints matrix A and vector b, and a point p this class calculates the
// projection x of p onto X = { y in R^n, A y <= b }
// TODO: find / create algorithm to solve the projection problem
template <typename _RealType, typename _VectorType, typename _MatrixType>
class PolyhedronProjector : public Projector<_RealType, _VectorType, _MatrixType> {
  typedef _RealType RealType;
  typedef _VectorType VectorType;
  typedef _MatrixType MatrixType;
protected:
  MatrixType _A;
  VectorType _b;
public:
  PolyhedronProjector ( const MatrixType &A, const VectorType &B )
    : Projector<RealType, VectorType, MatrixType> ( ), _A ( A ), _b ( B ) { }

  void apply ( const VectorType &Arg, VectorType &Dest ) const {
    if ( isFeasible ( Arg ) )
      Dest = Arg;
    // TODO: implement efficient projection onto polyhedron
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__);
  }

  bool isFeasible ( const VectorType &X ) const {
    VectorType test ( _b );
    _A.apply ( X, test );
    for ( int i=0; i<_b.size ( ) ; ++i ) {
      if ( test[i] > _b[i] )
        return false;
    }
    return true;
  }
};


template <typename _RealType, typename _VectorType, typename _MatrixType>
class BoxProjector : public Projector<_RealType, _VectorType, _MatrixType> {
  typedef _RealType RealType;
  typedef _VectorType VectorType;
  typedef _MatrixType MatrixType;
protected:
  const VectorType _lowerBounds, _upperBounds;
  const aol::BitVector _constrainedDirections;
public:
  BoxProjector ( const VectorType &LowerBounds, const VectorType &UpperBounds )
    : Projector<RealType, VectorType, MatrixType> ( ), _lowerBounds ( LowerBounds ), _upperBounds ( UpperBounds ),
      _constrainedDirections ( LowerBounds.Size ( ), true ) {
    if ( LowerBounds.size ( ) != UpperBounds.size ( ) )
      throw aol::Exception ( "Lower and upper bounds dimensions do not match!", __FILE__, __LINE__ );
  }

  BoxProjector ( const VectorType &LowerBounds, const VectorType &UpperBounds, const aol::BitVector &ConstrainedDirections )
      : Projector<RealType, VectorType, MatrixType> ( ), _lowerBounds ( LowerBounds ), _upperBounds ( UpperBounds ),
        _constrainedDirections ( ConstrainedDirections ) {
    if ( LowerBounds.size ( ) != UpperBounds.size ( ) )
      throw aol::Exception ( "Lower and upper bounds dimensions do not match!", __FILE__, __LINE__ );

    if ( ConstrainedDirections.size ( ) != LowerBounds.size ( ) )
      throw aol::Exception ( "Dimension of vector specifying the constrained directions does not match lower/upper variable bounds!", __FILE__, __LINE__ );
  }

  void apply ( const VectorType &Arg, VectorType &Dest ) const {
    if ( Arg.size ( ) != Dest.size ( ) || Arg.size ( ) != _lowerBounds.size ( ) )
      throw aol::Exception ( "Given variable dimensions do not match the dimensions of previously specified bounds!", __FILE__, __LINE__ );

    Dest = Arg;
    for ( int i=0; i<Arg.size ( ) ; ++i ) {
      if ( _constrainedDirections[i] ) {
        if ( Arg[i] < _lowerBounds[i] )
          Dest[i] = _lowerBounds[i];
        else if ( Arg[i] > _upperBounds[i] )
          Dest[i] = _upperBounds[i];
      }
    }
  }

  bool isFeasible ( const VectorType &Arg ) const {
    if ( Arg.size ( ) != _lowerBounds.size( ) )
      throw aol::Exception ( "Given variable dimension does not match the dimensions of previously specified bounds!", __FILE__, __LINE__ );

    for ( int i=0; i<Arg.size ( ) ; ++i ) {
      if ( Arg[i] < _lowerBounds[i] || Arg[i] > _upperBounds[i] )
        return false;
    }
    return true;
  }
};


#endif /* LINEARPROGRAMMING_H_ */
