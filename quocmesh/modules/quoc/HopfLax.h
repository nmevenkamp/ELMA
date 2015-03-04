#ifndef __HOPFLAX_H
#define __HOPFLAX_H

#include<gridBase.h>
#include<scalarArray.h>

namespace qc {

template <typename RealType>
class HopfLax {
protected:
  const qc::GridDefinition &_grid;
public:
  HopfLax( const qc::GridDefinition &grid ) : _grid( grid ) {
  }

  bool boundary( const CoordType &c ) {
    return c.x() == 0;
    //return ( c.x() == (_grid.getWidth()>>1)  && c.x() == c.y() );
    //return !( c.x() > 0 && c.x() < _grid.getWidth() - 1 && c.y() > 0 && c.y() < _grid.getWidth() - 1 );
  }

  void jacobiUpdate( const qc::ScalarArray<RealType, qc::QC_2D> &func ) {
    qc::ScalarArray<RealType, qc::QC_2D> oldfunc( _grid );
    oldfunc = func;
    qc::GridDefinition::OldFullNodeIterator it;
    for ( it = _grid.begin(); it != _grid.end(); ++it ) {
      if ( !boundary( *it ) ) {
  func.set( *it, getLocalHopfLaxUpdate( *it, oldfunc ) );
      }
    }
  }

  RealType getLocalHopfLaxUpdate( const CoordType &c, const qc::ScalarArray<RealType, qc::QC_2D> &oldfunc ) const {
    bool pos = ( oldfunc.get( c ) >= 0. );
    //RealType update = pos ? 1e20 : -1e20;
    RealType update = 1e20;

    vector<aol::Vec<2,CoordType> > segments;
    findNeighborSegments( c, segments );
    for ( vector<aol::Vec<2,CoordType> >::iterator sit=segments.begin(); sit!=segments.end(); sit++ ) {
      RealType Delta = (aol::Abs(oldfunc.get( (*sit)[1] )) - aol::Abs(oldfunc.get( (*sit)[0] ))) / _grid.H();
      RealType dist1, dist2;
      RealType cosalpha, cosbeta;
      if ( c.x() == (*sit)[0][0] || c.y() == (*sit)[0][1] ) {
  cosalpha = 0.;
  cosbeta =  1. / sqrt( 2. );
  dist1 = _grid.H();
  dist2 = _grid.H() * sqrt( 2. );
      } else {
  cosalpha = 1. / sqrt( 2. );
  cosbeta = 0.;
  dist2 = _grid.H();
  dist1 = _grid.H() * sqrt( 2. );
      }

      // here interpolation takes place
      if ( cosalpha < Delta ) {
  RealType v = oldfunc.get( (*sit)[0] );
  update = aol::Min( update, aol::Abs(v) + dist1 );
      } else if ( Delta <= -cosbeta ) {
  RealType v = oldfunc.get( (*sit)[1] );
  update = aol::Min( update, aol::Abs(v) + dist2 );
      } else {
  RealType v = oldfunc.get( (*sit)[0] );
  update = aol::Min( update, aol::Abs(v) + ( cosalpha * Delta + sqrt( ( 1. - aol::Sqr(cosalpha))*( 1. - aol::Sqr( Delta )))) * dist1 );
      }
    }
    return pos ? update : -update;
  }

  RealType getLocalHopfLaxUpdateAndExtend( const CoordType &c, const qc::ScalarArray<RealType, qc::QC_2D> &oldfunc, const qc::ScalarArray<RealType, qc::QC_2D> &toExtend, RealType &extendValue, RealType &distExtend, aol::Vec3<RealType> &extendFrom ) const {
    bool pos = ( oldfunc.get( c ) >= 0. );
    //RealType update = pos ? 1e20 : -1e20;
    RealType update = 1e20;
    RealType e = 0.;

    vector<aol::Vec<2,CoordType> > segments;
    findNeighborSegments( c, segments );
    int n=0;
    for ( vector<aol::Vec<2,CoordType> >::iterator sit=segments.begin(); sit!=segments.end(); sit++, n++ ) {
      RealType Delta = (aol::Abs(oldfunc.get( (*sit)[1] )) - aol::Abs(oldfunc.get( (*sit)[0] ))) / _grid.H();
      RealType s0 = toExtend.get( (*sit)[0] );
      RealType s1 = toExtend.get( (*sit)[1] );
      RealType dist1, dist2;
      RealType cosalpha, cosbeta;
      if ( c.x() == (*sit)[0][0] || c.y() == (*sit)[0][1] ) {
  cosalpha = 0.;
  cosbeta =  1. / sqrt( 2. );
  dist1 = _grid.H();
  dist2 = _grid.H() * sqrt( 2. );
      } else {
  cosalpha = 1. / sqrt( 2. );
  cosbeta = 0.;
  dist2 = _grid.H();
  dist1 = _grid.H() * sqrt( 2. );
      }

      // here interpolation takes place
      if ( cosalpha < Delta ) {
  RealType v = oldfunc.get( (*sit)[0] );
  RealType newV = aol::Abs(v) + dist1;
  if ( newV < update ) {
    update = newV;
    e = s0;
    extendFrom[0] = (*sit)[0][0];
    extendFrom[1] = (*sit)[0][1];
    distExtend = dist1;
  }
      } else if ( Delta <= -cosbeta ) {
  RealType v = oldfunc.get( (*sit)[1] );
  RealType newV = aol::Abs(v) + dist2;
  if ( newV < update ) {
    update = newV;
    e = s1;
    extendFrom[0] = (*sit)[1][0];
    extendFrom[1] = (*sit)[1][1];
    distExtend = dist2;
  }
      } else {
  RealType v = oldfunc.get( (*sit)[0] );
  RealType weight = ( cosalpha * Delta + sqrt( ( 1. - aol::Sqr(cosalpha))*( 1. - aol::Sqr( Delta ))));
  RealType newV = aol::Abs(v) + weight * dist1;
  if ( newV < update ) {
    update = newV;

    const RealType lambda = ( cosalpha - Delta * sqrt(   ( 1. - aol::Sqr(cosalpha)) / ( 1. - aol::Sqr( Delta )) ) ) * dist1 / _grid.H();
    extendFrom[0] = (*sit)[0][0] * ( 1. - lambda ) + (*sit)[1][0] * lambda;
    extendFrom[1] = (*sit)[0][1] * ( 1. - lambda ) + (*sit)[1][1] * lambda;
    e = s0 +  lambda * ( s1 - s0 );
    distExtend = weight * dist1;
  }
      }
    }


    extendValue = e;
    return pos ? update : -update;
  }

  void findNeighborSegments( const CoordType &c,
           vector<aol::Vec<2,CoordType> > &segments ) const {
    short minX = aol::Max( c.x() - 1, 0 );
    short maxX = aol::Min( c.x() + 1, _grid.getWidth()-1 );
    short minY = aol::Max( c.y() - 1, 0 );
    short maxY = aol::Min( c.y() + 1, _grid.getWidth()-1 );

    segments.resize( 8 );
    int n=0;
    for ( short Y = minY; Y < maxY; Y++ ) {
      segments[n][0][0]=minX;
      segments[n][0][1]=Y;
      segments[n][1][0]=minX;
      segments[n][1][1]=Y+1;

      n++;
      segments[n][0][0]=maxX;
      segments[n][0][1]=Y;
      segments[n][1][0]=maxX;
      segments[n][1][1]=Y+1;
      n++;
    }
    for ( short X = minX; X < maxX; X++ ) {
      segments[n][0][0]=X;
      segments[n][0][1]=minY;
      segments[n][1][0]=X+1;
      segments[n][1][1]=minY;
      n++;

      segments[n][0][0]=X;
      segments[n][0][1]=maxY;
      segments[n][1][0]=X+1;
      segments[n][1][1]=maxY;
      n++;
    }
    segments.resize( n );
  }
};

} // end namespace qc

#endif
