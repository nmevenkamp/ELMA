#ifndef __PROLONGATION_H
#define __PROLONGATION_H

#include <gridBase.h>
#include <sparseMatrices.h>
#include <scalarArray.h>

namespace qc {

template <typename RealType>
class ProlongOp : public aol::BiOp< aol::Vector<RealType> > {
protected:
  const GridDefinition &_coarseGrid;  //!< Coarse grid
  const GridDefinition &_fineGrid;    //!< Fine grid
  aol::OperatorType _opType;
  mutable aol::SparseMatrix< RealType > *mat; // may want to use specially designed matrix for this purpose to reduce computational and memory overhead

public:
  /** Construct the operator to prolongate between the given grids
   */
  ProlongOp ( const GridDefinition &Coarse,
              const GridDefinition &Fine,
              aol::OperatorType OpType = aol::ONTHEFLY ) :
      _coarseGrid ( Coarse ), _fineGrid ( Fine ), _opType ( OpType ), mat ( NULL ) {}

  ~ProlongOp() {
    if ( mat )
      delete ( mat );
  }

  void apply ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( _opType == aol::ONTHEFLY ) {
      std_mg_prolongate ( _fineGrid, Dest, _coarseGrid, Arg );
    }

    if ( _opType == aol::ASSEMBLED ) {
      if ( !mat ) {
        assembleMatrix();
      }
      mat->apply ( Arg, Dest );
    }
  }

  using aol::BiOp<aol::Vector< RealType> >::apply;

  virtual void assembleAddMatrix ( aol::SparseMatrix<RealType> &Mat ) const {
    // const is complete nonsense here. unfortunately, apply and apply add are const, thus require constness here.
    // could allow more general matrix, but want consistency with restriction operator

    assemble_std_prolong_matrix ( Mat );

  }

protected:
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    aol::Vector<RealType> tmp ( Dest.size() );
    apply ( Arg, tmp );
    Dest += tmp;
  }

  using aol::BiOp<aol::Vector< RealType> >::applyAdd;

  void assembleMatrix() const {
    // const is complete nonsense here. unfortunately, apply and apply add are const, thus require constness here.
    if ( mat != NULL )
      delete ( mat );
    mat = new aol::SparseMatrix<RealType> ( _fineGrid.getNumberOfNodes(), _coarseGrid.getNumberOfNodes() );
    assembleAddMatrix ( *mat );
  }

  void assemble_std_prolong_matrix ( aol::SparseMatrix<RealType> &mat ) const ;

  static void std_mg_prolongate ( const GridDefinition &FineGrid,
                                  aol::Vector<RealType> &FineVector,
                                  const GridDefinition &CoarseGrid,
                                  const aol::Vector<RealType> &CoarseVector );

  ProlongOp ( const ProlongOp<RealType>& ); // do not implement
  ProlongOp<RealType>& operator= ( const ProlongOp<RealType>& );
 // end class ProlongOp
};


namespace simplex {

//! Class for prolongation of data on a simplicial grid consistent with qc::simplex::TopologyLookup.
template <typename RealType, qc::Dimension Dim>
class ProlongOp :
  public aol::Op< aol::Vector<RealType> > {

protected:
  const GridDefinition &_coarseGrid;
  const GridDefinition &_fineGrid;

public:
  ProlongOp ( const GridDefinition &CoarseGrid,
              const GridDefinition &FineGrid ) :
    _coarseGrid ( CoarseGrid ),
    _fineGrid ( FineGrid ) {}

  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    qc::ScalarArray<RealType,Dim> arg( Arg, _coarseGrid, aol::FLAT_COPY ), dest( Dest, _fineGrid, aol::FLAT_COPY );

    for ( int z1 = 0, z2 = 0; z1 < _coarseGrid.getNumZ(); z1++, z2 += 2 )
      for ( int y1 = 0, y2 = 0; y1 < _coarseGrid.getNumY(); y1++, y2 += 2 )
        for ( int x1 = 0, x2 = 0; x1 < _coarseGrid.getNumX(); x1++, x2 += 2 ) {
          qc::CoordType pos1( x1, y1, z1 ), auxPos1, pos2( x2, y2, z2 ), auxPos2, aux;
          // copy values at old nodes
          dest.add( pos2, arg.get( pos1 ) );
          // interpolate values at new nodes on edges
          for ( int i = 0; i < Dim; ++i ) {
            auxPos1 = pos1;
            auxPos1[i] += 1;
            auxPos2 = pos2;
            auxPos2[i] += 1;
            if ( auxPos2[i] < _fineGrid.getSize()[i] )
              dest.add( auxPos2, .5 * ( arg.get( auxPos1 ) + arg.get( pos1 ) ) );
          }
          // interpolate values at new nodes on xy-faces
          auxPos1 = pos1;
          auxPos1[0] += 1;
          aux = pos1;
          aux[1] += 1;
          auxPos2 = pos2;
          auxPos2[0] += 1;
          auxPos2[1] += 1;
          if ( auxPos2[0] < _fineGrid.getSize()[0] && auxPos2[1] < _fineGrid.getSize()[1] )
            dest.add( auxPos2, .5 * ( arg.get( auxPos1 ) + arg.get( aux ) ) );
          // in 3d, interpolate values at new nodes on xz- and yz-faces and in the cube-middle
          if ( Dim == qc::QC_3D ) {
            auxPos1 = pos1;
            auxPos1[0] += 1;
            auxPos1[2] += 1;
            auxPos2 = pos2;
            auxPos2[0] += 1;
            auxPos2[2] += 1;
            if ( auxPos2[0] < _fineGrid.getSize()[0] && auxPos2[2] < _fineGrid.getSize()[2] )
              dest.add( auxPos2, .5 * ( arg.get( auxPos1 ) + arg.get( pos1 ) ) );
            auxPos1 = pos1;
            auxPos1[1] += 1;
            auxPos1[2] += 1;
            auxPos2 = pos2;
            auxPos2[1] += 1;
            auxPos2[2] += 1;
            if ( auxPos2[1] < _fineGrid.getSize()[1] && auxPos2[2] < _fineGrid.getSize()[2] )
              dest.add( auxPos2, .5 * ( arg.get( auxPos1 ) + arg.get( pos1 ) ) );
            auxPos1 = pos1;
            auxPos1[0] += 1;
            aux = pos1;
            aux[1] += 1;
            aux[2] += 1;
            auxPos2 = pos2;
            for ( int i = 0; i < Dim; ++i )
              auxPos2[i] += 1;
            if ( auxPos2[0] < _fineGrid.getSize()[0] && auxPos2[1] < _fineGrid.getSize()[1] && auxPos2[2] < _fineGrid.getSize()[2] )
              dest.add( auxPos2, .5 * ( arg.get( auxPos1 ) + arg.get( aux ) ) );
          }
        }
  }
}; // end class ProlongOpSimplex

} // end of namespace simplex

} // end of namespace qc

#endif
