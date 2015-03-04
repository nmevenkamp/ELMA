#ifndef COMPONENTSCOLLECTION_H_
#define COMPONENTSCOLLECTION_H_

// standard
#include <cmath>

// quocmesh
#include "component.h"
#include <imageTools.h>
#include <vectorExtensions.h>

template<typename _RealType>
class ComponentsCollection {
public:
  typedef _RealType RealType;

protected:
  aol::RandomAccessContainer<Component<RealType> > _components;
  qc::ScalarArray<int, qc::QC_2D> _labelArr;
  aol::RandomAccessContainer<Component<RealType> > _neighborhoods;
  qc::ScalarArray<int, qc::QC_2D> _neighborhoodArr;


public:
  ComponentsCollection ( ) {
  }

  ComponentsCollection ( const qc::ScalarArray<int, qc::QC_2D> &Labels )
    : _labelArr ( Labels ), _neighborhoodArr ( Labels ) {
    init( );
  }

  ComponentsCollection ( const qc::BitArray<qc::QC_2D> &Mask )
    : _labelArr ( qc::GridSize<qc::QC_2D>::createFrom ( Mask ) ), _neighborhoodArr ( qc::GridSize<qc::QC_2D>::createFrom ( Mask ) ) {
    qc::ConnectedComponentsLabeler::doLabel ( Mask, _labelArr );
    init ( );
  }
  
  void initializeFrom ( const qc::BitArray<qc::QC_2D> &Mask ) {
    _labelArr.reallocate ( qc::GridSize<qc::QC_2D>::createFrom ( Mask ) );
    _neighborhoodArr.reallocate ( qc::GridSize<qc::QC_2D>::createFrom ( Mask ) );
    qc::ConnectedComponentsLabeler::doLabel ( Mask, _labelArr );
    init ( );
  }

private:
  void init ( ) {
    for ( int i=0; i<=_labelArr.getMaxValue ( ) ; ++i )
      _components.pushBack ( Component<RealType> ( ) );
    for ( qc::RectangularIterator<qc::QC_2D, aol::Vec2<short> > it ( _labelArr ); it.notAtEnd ( ) ; ++it )
      if ( _labelArr.get ( *it ) )
        _components[_labelArr.get ( *it )].insert ( *it );
    for ( qc::RectangularIterator<qc::QC_2D, aol::Vec2<short> > it ( _neighborhoodArr ); it.notAtEnd ( ) ; ++it )
      _neighborhoodArr.set ( *it, 0 );
  }

public:
  void clearBoundaryComponents ( ) {
    qc::GridStructure grid ( qc::GridSize2d::createFrom ( _labelArr ) );
    for ( qc::GridStructure::AllBoundaryNodeIterator itArr ( grid ); itArr.notAtEnd ( ) ; ++itArr ) {
      const int componentNr = _labelArr.get ( *itArr );
      if ( componentNr ) {
        for ( typename Component<RealType>::const_iterator itSet = _components[componentNr].begin ( ); itSet != _components[componentNr].end ( ) ; ++itSet )
          _labelArr.set ( *itSet, 0 );
        _components[componentNr].clear ( );
      }
    }
  }

  void cropComponentsByGeometricCenter ( const aol::Vec2<short> &Origin, const aol::Vec2<short> &Size ) {
    const std::set<aol::Vec2<short> > centers = getGeometricCenters ( );
    for ( std::set<aol::Vec2<short> >::const_iterator it = centers.begin ( ); it != centers.end ( ) ; ++it ) {
      if ( (*it)[0] < Origin[0] || (*it)[0] >= Origin[0] + Size[0] || (*it)[1] < Origin[1] || (*it)[1] >= Origin[1] + Size[1] ) {
        const int componentNr = getComponentNumberByPosition ( *it );
        for ( typename Component<RealType>::const_iterator itSet = _components[componentNr].begin ( ); itSet != _components[componentNr].end ( ) ; ++itSet )
          _labelArr.set ( *itSet, 0 );
        _components[componentNr].clear ( );
      }
    }
  }

  void createPictureBWComponents ( qc::ScalarArray<RealType, qc::QC_2D> &Picture ) const {
    for ( qc::RectangularIterator<qc::QC_2D, aol::Vec2<short> > it ( _labelArr ); it.notAtEnd ( ) ; ++it )
      Picture.set ( *it, ( _labelArr.get ( *it ) > 0 ) );
    Picture.setOverflowHandlingToCurrentValueRange ( );
  }

  void createPictureBWComponentsRedGeometricCenters ( qc::MultiArray<RealType, qc::QC_2D, 3> &Picture ) {
    qc::GridStructure grid ( qc::GridSize2d::createFrom ( _labelArr ) );
    qc::ScalarArray<RealType, qc::QC_2D> bwComponents ( grid );
    createPictureBWComponents ( bwComponents );
    for ( qc::RectangularIterator<qc::QC_2D, aol::Vec2<short> > it ( grid ); it.notAtEnd ( ) ; ++it )
      for ( int k=0; k<3 ; ++k )
        Picture[k].set ( *it, bwComponents.get ( *it ) );
    const std::set<aol::Vec2<short> > geometricCenters = getGeometricCenters ( );
    for ( std::set<aol::Vec2<short> >::const_iterator it = geometricCenters.begin ( ); it != geometricCenters.end ( ) ; ++it ) {
      Picture[0].set ( *it, 1 );
      Picture[1].set ( *it, 0 );
      Picture[2].set ( *it, 0 );
    }
  }

  void createPictureBWComponentsRGBNeighborhoods ( qc::MultiArray<RealType, qc::QC_2D, 3> &Picture ) {
    int color = 0;
    for ( ComponentsCollection<RealType>::NonEmptyComponentsIterator it ( *this ); it.notAtEnd ( ) ; ++it ) {
      const Component<RealType> neighborhood = getNeighborhoodByComponentNumber ( it.getCurrentComponentNumber ( ) );
      for ( std::set<aol::Vec2<short> >::const_iterator nhIt=neighborhood.begin ( ); nhIt != neighborhood.end ( ) ; ++nhIt ) {
        for ( int k=0; k<3 ; ++k )
          Picture[k].set ( *nhIt, 0 );
        Picture[color].set ( *nhIt, 1 );
      }
      color++;
      color %= 3;
    }
    qc::GridStructure grid ( qc::GridSize2d::createFrom ( _labelArr ) );
    qc::ScalarArray<RealType, qc::QC_2D> bwComponents ( grid );
    createPictureBWComponents ( bwComponents );
    for ( qc::RectangularIterator<qc::QC_2D, aol::Vec2<short> > it ( grid ); it.notAtEnd ( ) ; ++it )
      if ( bwComponents.get( *it ) )
        for ( int k=0; k<3 ; ++k )
          Picture[k].set ( *it, bwComponents.get ( *it ) );
  }

  int getComponentNumberByPosition ( const aol::Vec2<short int> &Coord ) const {
    return _labelArr.get ( Coord );
  }

  int getComponentNumberByPosition ( const int X, const int Y ) const {
    return _labelArr.get ( X, Y );
  }

  Component<RealType>& getComponentByNumber ( const int I ) {
    return _components[I];
  }

  Component<RealType>& getComponentByPosition ( const aol::Vec2<short int> &Coord ) {
    return _components[getComponentNumberByPosition( Coord )];
  }

  Component<RealType>& getComponentByPosition ( const int X, const int Y ) {
    return _components[getComponentNumberByPosition( X , Y )];
  }

  std::set<aol::Vec2<short> > getGeometricCenters ( ) {
    std::set<aol::Vec2<short> > res;
    for ( NonEmptyComponentsIterator it ( *this ); it.notAtEnd ( ) ; ++it )
      res.insert ( (*it).getGeometricCenter( ) );
    return res;
  }

  int getNumNonEmptyComponents ( ) {
    int n = 0;
    for ( NonEmptyComponentsIterator it ( *this ); it.notAtEnd ( ) ; ++it )
      ++n;
    return n;
  }

  int getComponentNumberByNeighborhoodPosition ( const aol::Vec2<short int> &Coord ) {
    if ( _neighborhoodArr.size ( ) == 0 )
      naiveCalculateNeighborhoods ( );
    return _neighborhoodArr.get ( Coord );
  }

  int getComponentNumberByNeighborhoodPosition ( const int X, const int Y ) {
    if ( _neighborhoods.size ( ) == 0 )
      naiveCalculateNeighborhoods ( );
    return _neighborhoodArr.get ( X, Y );
  }

  const Component<RealType>& getNeighborhoodByComponentNumber ( const int I ) {
    if ( _neighborhoods.size ( ) == 0 )
      naiveCalculateNeighborhoods ( );
    return _neighborhoods[I];
  }

  void setSquareNeighborhoods ( const int Size ) {
    aol::Vec2<short> pos;
    const short neighborhoodOffset = ( Size - 1 ) / 2;
    for ( int i=0; i<=_labelArr.getMaxValue ( ) ; ++i )
      _neighborhoods.pushBack ( Component<RealType> ( ) );
    const std::set<aol::Vec2<short> > centers = getGeometricCenters ( );
    for ( std::set<aol::Vec2<short> >::const_iterator it = centers.begin ( ); it != centers.end ( ) ; ++it ) {
      const int componentNr = getComponentNumberByPosition ( *it );
      for ( int i=-neighborhoodOffset; i<=neighborhoodOffset ; ++i )
        for ( int j=-neighborhoodOffset; j<=neighborhoodOffset ; ++j )
          if ( (*it)[0] + i >= 0 && (*it)[0] + i < _labelArr.getNumXYZ ( ) && (*it)[1] + j >= 0 && (*it)[1] + j < _labelArr.getNumXYZ ( ) ) {
          pos.set ( i, j );
          pos += *it;
          _neighborhoods[componentNr].insert ( pos );
          _neighborhoodArr.set ( pos, componentNr );
        }
    }
  }

  void naiveCalculateNeighborhoods ( ) {
    RealType distance, minDist;
    for ( qc::RectangularIterator<qc::QC_2D, aol::Vec2<short> > arrIt ( _neighborhoodArr ); arrIt.notAtEnd ( ) ; ++arrIt ) {
      minDist = aol::NumberTrait<RealType>::Inf;
      for ( NonEmptyComponentsIterator compIt ( *this ); compIt.notAtEnd ( ) ; ++compIt ) {
        distance = (*compIt).naiveNearestNeighborSearchDistance ( *arrIt );
        if ( distance < minDist ) {
          minDist = distance;
          _neighborhoodArr.set ( *arrIt, compIt.getCurrentComponentNumber ( ) );
        }
      }
    }
    for ( int i=0; i<=_labelArr.getMaxValue ( ) ; ++i )
      _neighborhoods.pushBack ( Component<RealType> ( ) );
    for ( qc::RectangularIterator<qc::QC_2D, aol::Vec2<short> > it ( _neighborhoodArr ); it.notAtEnd ( ) ; ++it )
      if ( _neighborhoodArr.get ( *it ) )
        _neighborhoods[_neighborhoodArr.get ( *it )].insert ( *it );
  }

  class NonEmptyComponentsIterator {
  private:
    ComponentsCollection &_compCollection;
    int _curCompNr, _lastNonEmptyCompNr;
  public:
    typedef NonEmptyComponentsIterator Self;

    NonEmptyComponentsIterator ( ComponentsCollection &CompCollection ) : _compCollection ( CompCollection ) {
      _lastNonEmptyCompNr = 0;
      for ( int i=1; i<_compCollection._components.size ( ) ; ++i ) {
        if ( !_compCollection.getComponentByNumber ( i ).empty ( ) )
          _lastNonEmptyCompNr = i;
      }
      
      _curCompNr = 0;
      ++(*this);
    }

    bool notAtEnd ( ) const {
      return ( _curCompNr <= _lastNonEmptyCompNr );
    }

    Self& operator++ ( ) {
      ++_curCompNr;
      while ( notAtEnd ( ) && _compCollection.getComponentByNumber ( _curCompNr ).empty ( ) )
        ++_curCompNr;
      return *this;
    }

    Component<RealType>& operator* ( ) const {
      return _compCollection.getComponentByNumber ( _curCompNr );
    }

    int getCurrentComponentNumber ( ) const {
      return _curCompNr;
    }
  };
};


#endif /* COMPONENTSCOLLECTION_H_ */
