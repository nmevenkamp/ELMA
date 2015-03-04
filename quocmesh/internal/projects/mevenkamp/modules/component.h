#ifndef COMPONENT_H_
#define COMPONENT_H_


#include <aol.h>
#include <smallVec.h>
#include <vectorExtensions.h>
#include <geom.h>

template<typename _DataType>
class RandomAccessStencil : protected aol::RandomAccessContainer<_DataType> {
public:
  typedef _DataType DataType;

  RandomAccessStencil ( ) {
  }

  RandomAccessStencil ( const DataType &N, const DataType &E, const DataType &S, const DataType &W ) {
    set ( N, E, S, W );
  }

  void set ( const DataType &N, const DataType &E, const DataType &S, const DataType &W ) {
    this->clear ( );
    this->pushBack ( N );
    this->pushBack ( E );
    this->pushBack ( S );
    this->pushBack ( W );
  }

  const DataType& N ( ) const {
    return (*this)[0];
  }

  const DataType& E ( ) const {
    return (*this)[1];
  }

  const DataType& S ( ) const {
    return (*this)[2];
  }

  const DataType& W ( ) const {
    return (*this)[3];
  }

  bool empty ( ) const {
    return ( this->_data.size ( ) < 4 );
  }
};


template<typename _RealType>
class Component : public std::set<aol::Vec2<short> > {
public:
  typedef _RealType RealType;

protected:
  RandomAccessStencil<short> _boundaryBox;

public:
  Component ( ) { }

  const aol::Vec2<short> getGeometricCenter ( ) const {
    // short is too small to store the distance.
    aol::Vec2<RealType> center ( 0.0, 0.0 );
    for ( std::set<aol::Vec2<short> >::const_iterator it=this->begin ( ); it != this->end ( ) ; ++it ) {
      center[0] += (*it)[0];
      center[1] += (*it)[1];
    }
    center /= this->size ( );
    center[0] = round ( center[0] );
    center[1] = round ( center[1] );
    const aol::Vec2<short> res ( center );
    return res;
  }
  
  const aol::Vec<4, RealType> getMinMaxDiameterAndAngle ( ) {
    const aol::Vec2<short> geomCenter = getGeometricCenter ( );
    const Component bndry = boundary ( );
    const RealType maxD = 2 * aol::Max ( width ( ), height ( ) );
    RealType maxDiameter = 0.0, minDiameter = maxD, diameter = 0.0, maxAngle = 0.0, minAngle = 0.0;
    for ( int deg=0; deg<360 ; deg+=10 ) {
      const RealType rad = deg * aol::NumberTrait<RealType>::pi / 180.0;
      diameter = 0;
      for ( short sign=-1; sign<=1 ; sign+=2 ) {
        Component boundaryIntersections;
        for ( std::set<aol::Vec2<short> >::const_iterator it=bndry.begin ( ); it != bndry.end ( ) ; ++it ) {
          aol::LineSegment<RealType, qc::QC_2D> line ( aol::Vec2<RealType> ( -maxD * cos ( rad ), -maxD * sin ( rad ) ),
                                                       aol::Vec2<RealType> ( maxD * cos ( rad ), maxD * sin ( rad ) ) );
          aol::Vec2<RealType> pt ( (*it)[0] - geomCenter[0], (*it)[1] - geomCenter[1] ), projPt;
          line.projectTo ( pt, projPt );
          pt -= projPt;
          if ( pt.norm ( ) < 1 )
            boundaryIntersections.insert ( *it );
        }
        RealType intersectionDist = 0.0, maxIntersectionDist = 0.0;
        for ( std::set<aol::Vec2<short> >::const_iterator it=boundaryIntersections.begin ( ); it != boundaryIntersections.end ( ) ; ++it ) {
          intersectionDist = aol::Vec2<RealType> ( (*it)[0] - geomCenter[0], (*it)[1] - geomCenter[1] ).norm ( );
          if ( intersectionDist > maxIntersectionDist )
            maxIntersectionDist = intersectionDist;
        }
        diameter += maxIntersectionDist;
      }
      if ( diameter > maxDiameter ) {
        maxDiameter = diameter;
        maxAngle = rad;
      }
      if ( diameter < minDiameter ) {
        minDiameter = diameter;
        minAngle = rad;
      }
    }
    
    aol::Vec<4, RealType> res;
    res[0] = maxDiameter;
    res[1] = maxAngle;
    res[2] = minDiameter;
    res[3] = minAngle;
    return res;
  }
  
  const Component boundary ( ) const {
    Component bndry;
    aol::Vec2<short> pos;
    for ( std::set<aol::Vec2<short> >::const_iterator it=this->begin ( ); it != this->end ( ) ; ++it ) {
      bool isInBoundary = false;
      for ( int dx=-1; dx<=1 && !isInBoundary ; dx+=2 ) {
        pos.set ( *it );
        pos[0] += dx;
        if ( this->find ( pos ) == this->end ( ) ) {
          bndry.insert ( *it );
          isInBoundary = true;
        }
      }
      for ( int dy=-1; dy<=1 && !isInBoundary ; dy+=2 ) {
        pos.set ( *it );
        pos[1] += dy;
        if ( this->find ( pos ) == this->end ( ) ) {
          bndry.insert ( *it );
          isInBoundary = true;
        }
      }
    }
    
    return bndry;
  }

  RealType naiveNearestNeighborSearchDistance ( const aol::Vec2<short> &Point ) const {
    RealType distance = aol::NumberTrait<RealType>::Inf;
    for ( std::set<aol::Vec2<short> >::const_iterator it=this->begin ( ); it != this->end ( ) ; ++it ) {
      // short is too small to store the distance.
      const aol::Vec2<RealType> distanceVec ( (*it)[0] - Point[0], (*it)[1] - Point[1] );
      distance = aol::Min ( distance, distanceVec.normSqr ( ) );
    }
    return sqrt ( distance );
  }

  const RandomAccessStencil<short>& getBoundaryIndices ( ) {
    if ( _boundaryBox.empty ( ) ) {
      short N, E, S, W;
      // find left-most X coordinate
      W = 1000; // TODO: find more appropriate definition of upper bound
      for ( std::set<aol::Vec2<short> >::const_iterator it=this->begin ( ); it != this->end ( ) ; ++it )
        if ( (*it)[0] < W )
          W = (*it)[0];
      // find right-most X coordinate
      E = -1000;
      for ( std::set<aol::Vec2<short> >::const_iterator it=this->begin ( ); it != this->end ( ) ; ++it )
        if ( (*it)[0] > E )
          E = (*it)[0];
      // find bottom-most Y coordinate
      S = 1000;
      for ( std::set<aol::Vec2<short> >::const_iterator it=this->begin ( ); it != this->end ( ) ; ++it )
        if ( (*it)[1] < S )
          S = (*it)[1];
      // find top-most Y coordinate
      N = -1000;
      for ( std::set<aol::Vec2<short> >::const_iterator it=this->begin ( ); it != this->end ( ) ; ++it )
        if ( (*it)[1] > N )
          N = (*it)[1];
      _boundaryBox.set ( N, E, S, W );
    }
    return _boundaryBox;
  }

  short width ( ) {
    RandomAccessStencil<short> boundaryIndices = getBoundaryIndices ( );
    return boundaryIndices.E ( ) - boundaryIndices.W ( );
  }

  short height ( ) {
    RandomAccessStencil<short> boundaryIndices = getBoundaryIndices ( );
    return boundaryIndices.N ( ) - boundaryIndices.S ( );
  }
};


#endif /* COMPONENT_H_ */
