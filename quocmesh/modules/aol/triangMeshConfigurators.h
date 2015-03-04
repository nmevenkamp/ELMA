#ifndef __TRIANGMESHCONFIGURATORS_H
#define __TRIANGMESHCONFIGURATORS_H

#include <aol.h>
#include <configurators.h>
#include <FEOpInterface.h>
#include <solver.h>
#include <triangMesh.h>
#include <triangGaussQuadrature.h>


namespace aol {



//! BasefunctionSet for linear finite element functions on triangular meshes
//! NOTE: in order to use the gradient of the basis functions, 
//!       TriangleType has to provide a function that returns the inverted metric (first fund. form)
/*!
 * \author Droske, vDeylen, Heeren
 */

template <typename RealType, typename QuadRuleType, typename TriangleType>
class TriangMeshBaseFunctionSet : public aol::BaseFunctionSetInterface<RealType, aol::Vec3<RealType>,
  aol::Vec2<RealType>, 3, QuadRuleType, TriangMeshBaseFunctionSet<RealType, QuadRuleType, TriangleType> > {


  static RealType _b1   ( const aol::Vec2<RealType> &c ) { return 1. - c[0] - c[1]; }
  static RealType _b2   ( const aol::Vec2<RealType> &c ) { return c[0]; }
  static RealType _b3   ( const aol::Vec2<RealType> &c ) { return c[1]; }

  static RealType _d1_b1   ( const aol::Vec2<RealType> & ) { return - 1.; }
  static RealType _d1_b2   ( const aol::Vec2<RealType> & ) { return 1.; }
  static RealType _d1_b3   ( const aol::Vec2<RealType> & ) { return 0.; }

  static RealType _d2_b1   ( const aol::Vec2<RealType> & ) { return - 1.; }
  static RealType _d2_b2   ( const aol::Vec2<RealType> & ) { return 0.; }
  static RealType _d2_b3   ( const aol::Vec2<RealType> & ) { return 1.; }

  typedef RealType ( *BASIS_FUNC_TYPE ) ( const aol::Vec2<RealType> &RefCoord );
  BASIS_FUNC_TYPE _dbasis[2][3];
  BASIS_FUNC_TYPE _basis[3];

  const TriangleType *_triangle;

public:
  TriangMeshBaseFunctionSet(  ) : _triangle ( NULL ) {
    _basis[0] = _b1;
    _basis[1] = _b2;
    _basis[2] = _b3;

    _dbasis[0][0] = _d1_b1;
    _dbasis[0][1] = _d1_b2;
    _dbasis[0][2] = _d1_b3;

    _dbasis[1][0] = _d2_b1;
    _dbasis[1][1] = _d2_b2;
    _dbasis[1][2] = _d2_b3;
  }

  enum { numBaseFuncs = 3 };

  void setTriangle ( const TriangleType &T ) {
    _triangle = &T;
  }

  void evaluateGradient ( int BaseFuncNum, const aol::Vec2<RealType> &RefCoord, aol::Vec3<RealType> &Gradient ) const {
    // initialize vectors
    aol::Vec2<RealType> tmp, tmp2;
    // gradient at quad point in barycentric coords
    tmp[0] = _dbasis[0][BaseFuncNum] ( RefCoord );
    tmp[1] = _dbasis[1][BaseFuncNum] ( RefCoord );
    // change coordinate system
    _triangle->ginv(  ).mult ( tmp, tmp2 );
    // compute tangent vectors
    const aol::Vec3<RealType> dir0 = _triangle->edge(0,1);
    const aol::Vec3<RealType> dir1 = _triangle->edge(0,2);
    for ( int i = 0; i < 3; i++ ) {
      Gradient[i] = tmp2[0] * ( dir0[i] ) + tmp2[1] * ( dir1[i] );
    }
  }

  inline aol::Vec3<RealType> evaluateGradient ( int BaseFuncNum, int QuadPoint ) const {
    aol::Vec3<RealType> g;
    evaluateGradient ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ), g );
    return g;
  }

  RealType evaluate ( int BaseFuncNum, const aol::Vec2<RealType> &RefCoord ) const {
    return _basis[BaseFuncNum] ( RefCoord );
  }

  inline const RealType evaluate ( int BaseFuncNum, int QuadPoint ) const {
    return evaluate ( BaseFuncNum, this->_quadRule.getRefCoord ( QuadPoint ) );
  }
protected:

};


//! Mesh configurator for general triangular meshes, e.g. aol::TriangMesh<> or om::Trimesh<>.
/*!
 * \author Droske, vDeylen, Heeren
 */

template <typename _RealType, typename MeshType, typename _QuadType>
class TriangMeshConfigurator {
protected:
  const MeshType &_mesh;
public:
  typedef _RealType RealType;
  typedef _QuadType QuadType;

  static const qc::Dimension Dim = qc::QC_3D;
  enum { DomDim = 2 }; // ???

  static const qc::Dimension DimOfWorld = qc::QC_3D; // ??? currently not used, but should be correct anyway ...

  typedef RealType Real;
  typedef MeshType                                 InitType;               //!< that's the type that is needed by the constructor of the configurator
  typedef aol::Vec<2, RealType>                        DomVecType;
  typedef aol::Vec3<RealType>                          VecType;
  typedef aol::Mat<3, 3, RealType>                     MatType;
  typedef aol::Vector<RealType>                        VectorType;
  typedef aol::Vector<RealType>                        ArrayType;
  typedef aol::SparseMatrix<RealType>                  MatrixType;
  typedef aol::BitVector                               MaskType;
  typedef typename MeshType::ElementType           ElementType;
  typedef TriangMeshBaseFunctionSet<RealType, QuadType,
            ElementType>                               BaseFuncSetType;

  typedef typename MeshType::ElementIteratorType   ElementIteratorType;
  typedef typename MeshType::NodeIteratorType      NodeIteratorType;

  class DOFIterator : public NodeIteratorType {
    typedef TriangMeshConfigurator<_RealType, MeshType, _QuadType> ConfType;
  public:
    DOFIterator (const ConfType& conf ) : NodeIteratorType( conf.getInitializer() ){}
  };

  typedef DOFIterator DOFIteratorType;


  TriangMeshConfigurator ( const InitType &Mesh ) :
      _mesh ( Mesh ) {}

  //! returns the begin iterator of the grid
  inline const MeshType & begin( ) const {
    return _mesh;
  }

  //! returns the end iterator of the grid
  inline const MeshType & end( ) const {
    return _mesh;
  }

  const InitType& getInitializer( ) const { return this->_mesh; }



  mutable BaseFuncSetType _baseFuncSet;

  static const int maxNumLocalDofs = 3;

  inline int getNumLocalDofs ( const ElementType & ) const {
    return 3;
  }

  int getNumGlobalDofs( ) const {
    return this->_mesh.getNumVertices();
  }

  int maxNumQuadPoints( ) const {
    return QuadType::numQuadPoints;
  }


  const BaseFuncSetType& getBaseFunctionSet ( const ElementType &T ) const {
    _baseFuncSet.setTriangle ( T );
    return _baseFuncSet;
  }

  //! returns global index of the dof with number localIndex
  inline int localToGlobal ( const ElementType &T, int localIndex ) const {
    return T.globNodeIdx( localIndex );
  }

  //! \warning copied and pasted from QuocConfiguratorTraitMultiLin.
  //! I have absolutely no idea if this works or not.
  inline void localToGlobal ( const ElementType &T, const int localIndex0, const int localIndex1, aol::Vec2<int> &glob ) const {
    glob[0] = localToGlobal ( T, localIndex0 );
    glob[1] = localToGlobal ( T, localIndex1 );
  }

  RealType vol ( const ElementType &T ) const {
    return T.area();
  }

  //! create a new, clean matrix
  MatrixType* createNewMatrix( ) const {
    int num = getNumGlobalDofs();
    MatrixType *mat = new MatrixType ( num, num );
    return mat;
  }

  void fillBoundaryMask( MaskType& mask ) const {
    _mesh.fillBoundaryMask( mask );
  }

};

} // namespace aol

#endif // __TRIANGMESHCONFIGURATORS_H
