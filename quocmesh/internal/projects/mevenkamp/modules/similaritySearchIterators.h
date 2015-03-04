#ifndef SIMILARITYSEARCHITERATORS_H_
#define SIMILARITYSEARCHITERATORS_H_

#include <aol.h>
#include <vec.h>
#include <geom.h>

#include "patternAnalysis.h"


template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait>
class NeighborhoodFilter;

template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait, typename _BaseClass>
class EMNeighborhoodFilter;


template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait>
class SimilaritySearchIterator {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  NeighborhoodFilter<RealType, PictureType, _NeighborhoodFilterTrait> &_filter;
  short _xL, _xR, _yL, _yR;
  aol::Vec2<short> _cur;
public:
  SimilaritySearchIterator ( NeighborhoodFilter<RealType, PictureType, _NeighborhoodFilterTrait> &Filter )
    : _filter ( Filter ) { }
  
  virtual ~SimilaritySearchIterator ( ) { };
  
  virtual void reset ( const aol::Vec2<short> &XRef ) = 0;
  
  virtual void update ( const RealType Distance ) = 0;
  
  virtual bool notAtEnd ( ) const = 0;
  
  virtual SimilaritySearchIterator& operator++ ( ) = 0;
  
  aol::Vec2<short>& operator* ( ) {
    return _cur;
  }
  
  virtual bool isCurFinal ( ) const {
    return true;
  }
  
  virtual short getNumNonLocalWindows ( ) const {
    return 0;
  }
};


template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait>
class LocalSimilaritySearchIterator : public SimilaritySearchIterator<_RealType, _PictureType, _NeighborhoodFilterTrait> {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
public:
  LocalSimilaritySearchIterator ( NeighborhoodFilter<RealType, PictureType, _NeighborhoodFilterTrait> &Filter )
    : SimilaritySearchIterator<_RealType, _PictureType, _NeighborhoodFilterTrait> ( Filter ) { }
  
  virtual void reset ( const aol::Vec2<short> &XRef ) {
    short searchWindowOffset = ( this->_filter._searchWindowSize - 1 ) / 2;
    this->_xL = aol::Max<short> ( XRef[0] - searchWindowOffset, this->_filter._X0 );
    this->_xR = aol::Min<short> ( XRef[0] + searchWindowOffset, this->_filter._XEnd - 1 );
    this->_yL = aol::Max<short> ( XRef[1] - searchWindowOffset, this->_filter._Y0 );
    this->_yR = aol::Min<short> ( XRef[1] + searchWindowOffset, this->_filter._YEnd - 1 );
    this->_cur.set ( this->_xL, this->_yL );
  }
  
  void update ( const RealType /*Distance*/ ) { }
  
  bool notAtEnd ( ) const {
    return ( this->_cur[0] <= this->_xR && this->_cur[1] <= this->_yR );
  }
  
  LocalSimilaritySearchIterator& operator++ ( ) {
    ++this->_cur[0];
    if ( this->_cur[0] > this->_xR ) {
      this->_cur[0] = this->_xL;
      ++this->_cur[1];
    }
    
    return *this;
  }
};


template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait>
class GlobalSimilaritySearchIterator : public LocalSimilaritySearchIterator<_RealType, _PictureType, _NeighborhoodFilterTrait> {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
public:
  GlobalSimilaritySearchIterator ( NeighborhoodFilter<RealType, PictureType, _NeighborhoodFilterTrait> &Filter )
    : LocalSimilaritySearchIterator<_RealType, _PictureType, _NeighborhoodFilterTrait> ( Filter ) { }
  
  void reset ( const aol::Vec2<short> &/*XRef*/ ) {
    this->_xL = this->_filter._X0;
    this->_xR = this->_filter._XEnd - 1;
    this->_yL = this->_filter._Y0;
    this->_yR = this->_filter._YEnd -1;
    this->_cur.set ( this->_xL, this->_yL );
  }
};


template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait, typename _BaseClass>
class PeriodicSimilaritySearchIterator : public SimilaritySearchIterator<_RealType, _PictureType, _NeighborhoodFilterTrait> {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  EMNeighborhoodFilter<RealType, PictureType, _NeighborhoodFilterTrait, _BaseClass> &_filter;
  
  const short _localNs = 5;
  const short _localNsOffset;
  
  aol::Vector<RealType> _angles, _periods;
  aol::Vec2<short> _localMinCoords, _curAxisDirection, _curFirstAxisCenter, _curSecondAxisCenter;
  RealType _localMinDist;
  short _localXL, _localXR, _localYL, _localYR;
  
  const aol::AlignedQuad<RealType> _c;
  
  short _numNonLocalWindows;
public:
  PeriodicSimilaritySearchIterator ( EMNeighborhoodFilter<RealType, PictureType, _NeighborhoodFilterTrait, _BaseClass> &Filter )
    : SimilaritySearchIterator<RealType, PictureType, _NeighborhoodFilterTrait> ( Filter ), _filter ( Filter ),
      _localNsOffset ( ( _localNs - 1 ) / 2 ),
      _angles ( 2 ), _periods ( 2 ),
      _c ( aol::Vec2<RealType> ( 0, 0 ), aol::Vec2<RealType> ( Filter._input.getNumX ( )-1, Filter._input.getNumY ( )-1 ) ) {
    PatternAnalyzer<RealType, PictureType> patternAnalyzer ( _filter._input, _filter._outputDir, true, Filter._progressBar );
    _angles = patternAnalyzer.getPeriodicityAnglesRadians ( );
    _periods = patternAnalyzer.getPeriodicitySpacingsPixels ( );
    if ( !_filter._quietMode )
      std::cerr << "Periodicity analysis: angles=" << _angles << "; spacings=" << _periods << std::endl;
  }
  
  void reset ( const aol::Vec2<short> &XRef ) {
    this->_xL = _filter._X0;
    this->_xR = _filter._XEnd - 1;
    this->_yL = _filter._Y0;
    this->_yR = _filter._YEnd - 1;
    periodicityReset ( XRef );
    localReset ( XRef );
    _numNonLocalWindows = 0;
  }
  
  void update ( const RealType Distance ) {
    if ( Distance < _localMinDist ) {
      _localMinDist = Distance;
      _localMinCoords.set ( this->_cur );
    }
  }
  
  bool notAtEnd ( ) const {
    return ( this->_cur[0] >= this->_xL && this->_cur[0] <= this->_xR && this->_cur[1] >= this->_yL && this->_cur[1] <= this->_yR );
  }
  
  PeriodicSimilaritySearchIterator& operator++ ( ) {
    if ( this->_cur[0] < _localXR )
      ++this->_cur[0];
    else if ( this->_cur[1] < _localYR) {
      this->_cur[0] = _localXL;
      ++this->_cur[1];
    } else {
      this->_cur.set ( _localMinCoords );
      advanceAlongAxis ( 1 );
      if ( !notAtEnd ( ) ) {
        if ( _curAxisDirection[1] == 1 ) {
          this->_cur.set ( _curSecondAxisCenter );
          _curAxisDirection[1] = -1;
          advanceAlongAxis ( 1 );
        }
        if ( _curAxisDirection[1] == -1 && !notAtEnd ( ) ) {
          this->_cur.set ( _curSecondAxisCenter );
          advanceAlongAxis ( 0 );
          _curSecondAxisCenter.set ( this->_cur );
          _curAxisDirection[1] = 1;
          if ( _curAxisDirection[0] == 1 && !notAtEnd ( ) ) {
            this->_cur.set ( _curFirstAxisCenter );
            _curAxisDirection[0] = -1;
            advanceAlongAxis ( 0 );
            _curSecondAxisCenter.set ( this->_cur );
            _curAxisDirection[1] = 1;
          }
        }
      }
      if ( notAtEnd ( ) ) {
        localReset ( this->_cur );
        ++_numNonLocalWindows;
      }
    }
    
    return *this;
  }
  
  short getNumNonLocalWindows ( ) const {
    return _numNonLocalWindows;
  }
private:
  void periodicityReset ( const aol::Vec2<short> &Center ) {
    _curFirstAxisCenter.set ( Center );
    _curSecondAxisCenter.set ( Center );
    _curAxisDirection.set ( 1, 1 );
    
    // Swap angles and periods such that the intersection of the line through Center with angle _angles[0] with the image is the greater one
    aol::RandomAccessContainer<aol::Vec2<RealType> > pts ( 4, Center );
    for ( int axis=0; axis<2 ; ++axis ) {
      pts[2*axis] += 2 * ( _c._max[0] - _c._min[0] ) * aol::Vec2<RealType> ( cos ( _angles[axis] ), sin ( _angles[axis] ) );
      pts[2*axis+1] -= 2 * ( _c._max[0] - _c._min[0] ) * aol::Vec2<RealType> ( cos ( _angles[axis] ), sin ( _angles[axis] ) );
    }
    aol::LineSegment<RealType, qc::QC_2D> axis1 ( pts[0], pts[1] );
    aol::LineSegment<RealType, qc::QC_2D> axis2 ( pts[2], pts[3] );
    if ( aol::Intersector<aol::AlignedQuad<RealType>, aol::LineSegment<RealType, qc::QC_2D> >::cutLength ( _c, axis1 )
        <  aol::Intersector<aol::AlignedQuad<RealType>, aol::LineSegment<RealType, qc::QC_2D> >::cutLength ( _c, axis2 ) ) {
      RealType tmp = _angles[0];
      _angles[0] = _angles[1];
      _angles[1] = tmp;
      tmp = _periods[0];
      _periods[0] = _periods[1];
      _periods[1] = tmp;
    }
  }
  
  void localReset ( const aol::Vec2<short> &Center ) {
    _localMinCoords.set ( Center );
    _localMinDist = aol::NumberTrait<RealType>::getInf ( );
    _localXL = aol::Max<short> ( Center[0] - _localNsOffset, _filter._X0 );
    _localXR = aol::Min<short> ( Center[0] + _localNsOffset, _filter._XEnd - 1 );
    _localYL = aol::Max<short> ( Center[1] - _localNsOffset, _filter._Y0 );
    _localYR = aol::Min<short> ( Center[1] + _localNsOffset, _filter._YEnd - 1 );
    this->_cur.set ( _localXL, _localYL );
  }
  
  void advanceAlongAxis ( const short Axis ) {
    this->_cur[0] += _curAxisDirection[Axis] * _periods[Axis] * cos ( _angles[Axis] );
    this->_cur[1] += _curAxisDirection[Axis] * _periods[Axis] * sin ( _angles[Axis] );
  }
};


#endif /* SIMILARITYSEARCH_H_ */
