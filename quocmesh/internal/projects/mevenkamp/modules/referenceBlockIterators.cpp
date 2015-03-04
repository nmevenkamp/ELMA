#include "referenceBlockIterators.h"
#include "neighborhoodFilter.h"
#include "nonLocalMeansFilter.h"
#include "bm3dFilter.h"
#include "emNonLocalMeansFilter.h"
#include "emBM3DFilter.h"


template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait>
ReferenceBlockIterator<_RealType, _PictureType, _NeighborhoodFilterTrait>::ReferenceBlockIterator ( NeighborhoodFilter<_RealType, _PictureType, _NeighborhoodFilterTrait> &Filter ) : _filter ( Filter ) { }


/*
 *  GlobalReferenceBlockIterator
 */
template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait>
GlobalReferenceBlockIterator<_RealType, _PictureType, _NeighborhoodFilterTrait>::GlobalReferenceBlockIterator ( NeighborhoodFilter<_RealType, _PictureType, _NeighborhoodFilterTrait> &Filter )
  : ReferenceBlockIterator<_RealType, _PictureType, _NeighborhoodFilterTrait> ( Filter ),
    _xL ( Filter._X0 ), _xR ( Filter._XEnd - 1 ), _yL ( Filter._Y0 ), _yR ( Filter._YEnd - 1 ),
    _refStep ( 1 ), _cur ( _xL, _yL ) {
  if ( this->_filter._progressBar != NULL )
    this->_filter._progressBar->start ( this->_filter._input.size ( ) );
}


template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait>
GlobalReferenceBlockIterator<_RealType, _PictureType, _NeighborhoodFilterTrait>::GlobalReferenceBlockIterator ( CollaborativeNeighborhoodFilter<_RealType, _PictureType, _NeighborhoodFilterTrait> &Filter )
: ReferenceBlockIterator<_RealType, _PictureType, _NeighborhoodFilterTrait> ( Filter ),
  _xL ( Filter._X0 ), _xR ( Filter._XEnd - 1 ), _yL ( Filter._Y0 ), _yR ( Filter._YEnd - 1 ),
  _refStep ( Filter._refStep ), _cur ( _xL, _yL ) {
  if ( this->_filter._progressBar != NULL )
    this->_filter._progressBar->start ( this->_filter._input.size ( ) / ( _refStep * _refStep ) );
}


template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait>
bool GlobalReferenceBlockIterator<_RealType, _PictureType, _NeighborhoodFilterTrait>::notAtEnd ( ) const {
  if ( _cur[0] <= _xR && _cur[1] <= _yR ) return true;
  else {
    if ( this->_filter._progressBar != NULL )
      this->_filter._progressBar->finish ( );
    return false;
  }
}


template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait>
aol::Vec2<short>& GlobalReferenceBlockIterator<_RealType, _PictureType, _NeighborhoodFilterTrait>::operator* ( ) {
  return this->_cur;
}


template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait>
GlobalReferenceBlockIterator<_RealType, _PictureType, _NeighborhoodFilterTrait>& GlobalReferenceBlockIterator<_RealType, _PictureType, _NeighborhoodFilterTrait>::operator++ ( ) {
  _cur[0] += _refStep;
  if ( _cur[0] > _xR ) {
    if ( _cur[0] - _refStep < _xR ) _cur[0] = _xR;
    else {
      _cur[0] = _xL;
      _cur[1] += _refStep;
      if ( _cur[1] > _yR && _cur[1] - _refStep < _yR ) _cur[1] = _yR;
    }
  }
  ReferenceBlockIterator<_RealType, _PictureType, _NeighborhoodFilterTrait>::operator++ ( );
  
  return *this;
}


/*
 *  PixelReferenceBlockIterator
 */
template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait>
PixelReferenceBlockIterator<_RealType, _PictureType, _NeighborhoodFilterTrait>::PixelReferenceBlockIterator ( NeighborhoodFilter<_RealType, _PictureType, _NeighborhoodFilterTrait> &Filter )
  : ReferenceBlockIterator<_RealType, _PictureType, _NeighborhoodFilterTrait> ( Filter ), _curIdx ( 0 ), _curPos ( Filter._pixels[0][0], Filter._pixels[0][1] ), _pixels ( Filter._pixels ) {
  this->_filter._progressBar->start ( _pixels.numComponents ( ) );
}


template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait>
bool PixelReferenceBlockIterator<_RealType, _PictureType, _NeighborhoodFilterTrait>::notAtEnd ( ) const {
  if ( _curIdx < _pixels.numComponents ( ) ) return true;
  else {
    this->_filter._progressBar->finish ( );
    return false;
  }
}


template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait>
aol::Vec2<short>& PixelReferenceBlockIterator<_RealType, _PictureType, _NeighborhoodFilterTrait>::operator* ( ) {
  return _curPos;
}


template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait>
PixelReferenceBlockIterator<_RealType, _PictureType, _NeighborhoodFilterTrait>& PixelReferenceBlockIterator<_RealType, _PictureType, _NeighborhoodFilterTrait>::operator++ ( ) {
  ++_curIdx;
  if ( notAtEnd ( ) )
    _curPos.set ( _pixels[_curIdx][0], _pixels[_curIdx][1] );
  ReferenceBlockIterator<_RealType, _PictureType, _NeighborhoodFilterTrait>::operator++ ( );
  
  return *this;
}



template class ReferenceBlockIterator<double, qc::ScalarArray<double, qc::QC_2D>, NonLocalMeansFilterTrait<double, qc::ScalarArray<double, qc::QC_2D> > >;
template class ReferenceBlockIterator<double, qc::ScalarArray<double, qc::QC_2D>, BM3DFilterTrait<double, qc::ScalarArray<double, qc::QC_2D> > >;
template class ReferenceBlockIterator<double, qc::ScalarArray<double, qc::QC_2D>, EMNonLocalMeansFilterTrait<double, qc::ScalarArray<double, qc::QC_2D> > >;
template class ReferenceBlockIterator<double, qc::ScalarArray<double, qc::QC_2D>, EMBM3DFilterTrait<double, qc::ScalarArray<double, qc::QC_2D> > >;

template class GlobalReferenceBlockIterator<double, qc::ScalarArray<double, qc::QC_2D>, NonLocalMeansFilterTrait<double, qc::ScalarArray<double, qc::QC_2D> > >;
template class GlobalReferenceBlockIterator<double, qc::ScalarArray<double, qc::QC_2D>, BM3DFilterTrait<double, qc::ScalarArray<double, qc::QC_2D> > >;
template class GlobalReferenceBlockIterator<double, qc::ScalarArray<double, qc::QC_2D>, EMNonLocalMeansFilterTrait<double, qc::ScalarArray<double, qc::QC_2D> > >;
template class GlobalReferenceBlockIterator<double, qc::ScalarArray<double, qc::QC_2D>, EMBM3DFilterTrait<double, qc::ScalarArray<double, qc::QC_2D> > >;

template class PixelReferenceBlockIterator<double, qc::ScalarArray<double, qc::QC_2D>, NonLocalMeansFilterTrait<double, qc::ScalarArray<double, qc::QC_2D> > >;
template class PixelReferenceBlockIterator<double, qc::ScalarArray<double, qc::QC_2D>, BM3DFilterTrait<double, qc::ScalarArray<double, qc::QC_2D> > >;
template class PixelReferenceBlockIterator<double, qc::ScalarArray<double, qc::QC_2D>, EMNonLocalMeansFilterTrait<double, qc::ScalarArray<double, qc::QC_2D> > >;
template class PixelReferenceBlockIterator<double, qc::ScalarArray<double, qc::QC_2D>, EMBM3DFilterTrait<double, qc::ScalarArray<double, qc::QC_2D> > >;