#ifndef REFERENCEBLOCKITERATORS_h
#define REFERENCEBLOCKITERATORS_h

#include <aol.h>
#include <vec.h>
#include <multiVector.h>

template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait>
class NeighborhoodFilter;

template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait>
class CollaborativeNeighborhoodFilter;


template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait>
class ReferenceBlockIterator {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  NeighborhoodFilter<RealType, PictureType, _NeighborhoodFilterTrait> &_filter;
public:
  ReferenceBlockIterator ( NeighborhoodFilter<RealType, PictureType, _NeighborhoodFilterTrait> &Filter );
  
  virtual bool notAtEnd ( ) const = 0;
  
  virtual ReferenceBlockIterator& operator++ ( ) {
    if ( _filter._progressBar != NULL ) (*_filter._progressBar)++;
    return *this;
  };

  virtual aol::Vec2<short>& operator* ( ) = 0;
};


template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait>
class GlobalReferenceBlockIterator : public ReferenceBlockIterator<_RealType, _PictureType, _NeighborhoodFilterTrait> {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  short _xL, _xR, _yL, _yR, _refStep;
  aol::Vec2<short> _cur;
public:
  GlobalReferenceBlockIterator ( NeighborhoodFilter<RealType, PictureType, _NeighborhoodFilterTrait> &Filter );
  
  GlobalReferenceBlockIterator ( CollaborativeNeighborhoodFilter<RealType, PictureType, _NeighborhoodFilterTrait> &Filter );
  
  bool notAtEnd ( ) const;
  
  GlobalReferenceBlockIterator& operator++ ( );
  
  aol::Vec2<short>& operator* ( );
};


template <typename _RealType, typename _PictureType, typename _NeighborhoodFilterTrait>
class PixelReferenceBlockIterator : public ReferenceBlockIterator<_RealType, _PictureType, _NeighborhoodFilterTrait> {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  int _curIdx;
  aol::Vec2<short> _curPos;
  aol::MultiVector<short> _pixels;
public:
  PixelReferenceBlockIterator ( NeighborhoodFilter<RealType, PictureType, _NeighborhoodFilterTrait> &Filter );
  
  bool notAtEnd ( ) const;
  
  PixelReferenceBlockIterator& operator++ ( );
  
  aol::Vec2<short>& operator* ( );
};


#endif /* REFERENCEBLOCKITERATORS_H_ */
