#include "atomFinder.h"


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
                _ColoredPictureType>::getComponentsCollection ( ComponentsCollection<RealType> &ComponentsCollection,
                                                                qc::BitArray<qc::QC_2D> &Segmented,
                                                                const PictureType &Data,
                                                                const RealType Gamma, const int MaxIt, const RealType Epsilon ) const {
  if ( !_quietMode )
    std::cerr << "Approximating atom positions as geometric centers of connected components.." << std::endl;
  
  PictureType u ( Data );
  u /= u.getMaxValue ( );
  u.setOverflowHandling( aol::CLIP_THEN_SCALE, 0, 1 );
  
  InitType grid ( qc::GridSize<ConfiguratorType::Dim> ( u.getSize ( ) ) );
  
  ArrayType uArr ( u, aol::FLAT_COPY );
  qc::PiecewiseConstantTwoPhaseMSSegmentor<ConfiguratorType, 1, qc::FirstOrderPrimalTwoPhaseMSSegmentor<ConfiguratorType> > segmentor ( grid, Gamma, uArr );
  segmentor.setQuietMode ( this->_quietMode );
  segmentor.setCatchCtrlC ( true );
  PictureType segmentation ( grid );
  segmentor.setMaxIterations ( MaxIt );
  segmentor.setStopEpsilon ( Epsilon );
  segmentor.segmentAndAdjustGrayValues ( segmentation );
  
  if ( _diskOutput ) {
    std::stringstream ss;
    ss << _outputDir << "/segmented.png";
    segmentation.setOverflowHandlingToCurrentValueRange ( );
    segmentation.savePNG ( ss.str ( ).c_str ( ) );
  }

  segmentation.threshold ( 0.5, 0, 1 );
  
  if ( _diskOutput ) {
    std::stringstream ss;
    ss << _outputDir << "/segmentedAndThresholded.png";
    segmentation.setOverflowHandlingToCurrentValueRange ( );
    segmentation.savePNG ( ss.str ( ).c_str ( ) );
  }
  
  qc::BitArray<qc::QC_2D> mask ( grid );
  for ( int i=0; i<segmentation.getNumXYZ ( ); ++i )
    for ( int j=0; j<segmentation.getNumXYZ ( ); ++j )
      mask.set ( i, j, ( segmentation.get ( i, j ) == 1 ) );

  // We assume that atoms intensities are higher than the void intensities.
  // So we have to make sure that the mask corresponds to regions of higher intensities.
  if ( segmentor.getMeanValuesReference()[0][0] > segmentor.getMeanValuesReference()[1][0] )
    mask.invert();

  Segmented.reallocate ( mask.getNumX ( ), mask.getNumY ( ) );
  Segmented = mask;
  
  ComponentsCollection.initializeFrom ( mask );
  
  if ( !_quietMode )
    std::cerr << "# of atoms after segmentation: " << ComponentsCollection.getNumNonEmptyComponents ( ) << std::endl;
  
  if ( _diskOutput ) {
    PictureType connectedRegions ( grid );
    ComponentsCollection.createPictureBWComponents ( connectedRegions );
    std::stringstream ss;
    ss << _outputDir << "/connectedRegions" << getDefaultArraySuffix ( qc::QC_2D );
    connectedRegions.save ( ss.str ( ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    
    ColoredPictureType geometricCentersImg ( grid );
    ComponentsCollection.createPictureBWComponentsRedGeometricCenters ( geometricCentersImg );
    ss.str ( std::string ( ) ); // clear stringstream
    ss << _outputDir << "/geometricCenters.png";
    geometricCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
    geometricCentersImg.savePNG ( ss.str ( ).c_str ( ) );
  }
}


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
_ColoredPictureType>::getApproximateAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<int> &Dimensions,
                                                    aol::MultiVector<RealType> &DumbbellCenters, aol::MultiVector<int> &DumbbellDimensions,
                                                    aol::Vector<RealType> &DumbbellOrientations, aol::Vector<RealType> &DumbbellSeparations,
                                                    qc::BitArray<qc::QC_2D> &Segmented,
                                                    const PictureType &Data,
                                                    const RealType Gamma, const int MaxIt, const RealType Epsilon ) const {
  ComponentsCollection<RealType> componentsCollection;
  getComponentsCollection ( componentsCollection, Segmented, Data, Gamma, MaxIt, Epsilon );
  
  // Split components into dumbbells and single atoms and collect orientation and separation for dumbbells
  int numNonEmptyComps = componentsCollection.getNumNonEmptyComponents ( );
  Centers.reallocate ( numNonEmptyComps, 2 );
  Dimensions.reallocate ( numNonEmptyComps, 2 );
  DumbbellCenters.reallocate ( numNonEmptyComps, 2 );
  DumbbellDimensions.reallocate ( numNonEmptyComps, 2 );
  DumbbellOrientations.reallocate ( numNonEmptyComps );
  DumbbellSeparations.reallocate ( numNonEmptyComps );
  int kSingle = 0, kDumbbell = 0;
  for ( NonEmptyComponentsIterator it ( componentsCollection ); it.notAtEnd ( ) ; ++it ) {
    const RandomAccessStencil<short> boundaryIndices = (*it).getBoundaryIndices ( );
    aol::Vec2<short> geometricCenter = (*it).getGeometricCenter ( );
    const short width = 1 + 2 * round ( 1.25 * aol::Max<short> ( aol::Abs<short> ( boundaryIndices.W ( ) - geometricCenter[0] ), aol::Abs<short> ( boundaryIndices.E ( ) - geometricCenter[0] ) ) );
    const short height = 1 + 2 * round ( 1.25 * aol::Max<short> ( aol::Abs<short> ( boundaryIndices.N ( ) - geometricCenter[1] ), aol::Abs<short> (  boundaryIndices.S ( ) - geometricCenter[1] ) ) );
    aol::Vec<4,RealType> minMaxDiamAngle = (*it).getMinMaxDiameterAndAngle ( );
    if ( minMaxDiamAngle[0] > 1.5 * ( minMaxDiamAngle[2] ) ) {
      DumbbellCenters[kDumbbell][0] = geometricCenter[0];
      DumbbellCenters[kDumbbell][1] = geometricCenter[1];
      DumbbellDimensions[kDumbbell][0] = width;
      DumbbellDimensions[kDumbbell][1] = height;
      DumbbellOrientations[kDumbbell] = minMaxDiamAngle[1];
      DumbbellSeparations[kDumbbell] = 0.5 * minMaxDiamAngle[0];
      ++kDumbbell;
    } else {
      Centers[kSingle][0] = geometricCenter[0];
      Centers[kSingle][1] = geometricCenter[1];
      Dimensions[kSingle][0] = width;
      Dimensions[kSingle][1] = height;
      ++kSingle;
    }
  }
  DumbbellCenters.resize ( kDumbbell, 2 );
  DumbbellDimensions.resize ( kDumbbell, 2 );
  DumbbellOrientations.resize ( kDumbbell );
  DumbbellSeparations.resize ( kDumbbell );
  Centers.resize ( kSingle, 2 );
  Dimensions.resize ( kSingle, 2 );
  
  if ( _diskOutput ) {
    // Add orientation lines to the dumbbell atoms
    std::stringstream ss;
    ss << _outputDir << "/geometricCenters.png";
    ColoredPictureType geometricCentersImg ( ss.str ( ).c_str ( ) );
    aol::Vec2<short> pos;
    RealType maxVal = geometricCentersImg.getMaxValue ( );
    for ( int k=0; k<kDumbbell ; ++k ) {
      for ( short i=-DumbbellSeparations[k]; i<=DumbbellSeparations[k] ; ++i ) {
        pos.set ( round ( DumbbellCenters[k][0] + cos ( DumbbellOrientations[k] ) * i ), round ( DumbbellCenters[k][1] + sin ( DumbbellOrientations[k] ) * i ) );
        if ( pos[0] >=0 && pos[0] < geometricCentersImg[0].getNumX ( ) && pos[1] >=0 && pos[1] < geometricCentersImg[0].getNumY ( ) ) {
          geometricCentersImg[0].set ( pos, 0 );
          geometricCentersImg[1].set ( pos, maxVal );
          geometricCentersImg[2].set ( pos, 0 );
        }
      }
      pos.set ( DumbbellCenters[k][0], DumbbellCenters[k][1] );
      geometricCentersImg[0].set ( pos, maxVal );
      geometricCentersImg[1].set ( pos, 0 );
    }
    geometricCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
    geometricCentersImg.savePNG ( ss.str ( ).c_str ( ) );
  }
}


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
                _ColoredPictureType>::getApproximateSingleAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<int> &Dimensions,
                                                                          qc::BitArray<qc::QC_2D> &Segmented,
                                                                          const PictureType &Data,
                                                                          const RealType Gamma, const int MaxIt, const RealType Epsilon ) const {
  ComponentsCollection<RealType> componentsCollection;
  getComponentsCollection ( componentsCollection, Segmented, Data, Gamma, MaxIt, Epsilon );
  
  // Split components into dumbbells and single atoms and collect orientation and separation for dumbbells
  int numNonEmptyComps = componentsCollection.getNumNonEmptyComponents ( );
  Centers.reallocate ( numNonEmptyComps, 2 );
  Dimensions.reallocate ( numNonEmptyComps, 2 );
  int k = 0;
  for ( NonEmptyComponentsIterator it ( componentsCollection ); it.notAtEnd ( ) ; ++it ) {
    const RandomAccessStencil<short> boundaryIndices = (*it).getBoundaryIndices ( );
    aol::Vec2<short> geometricCenter = (*it).getGeometricCenter ( );
    const short width = 1 + 2 * round ( 1.25 * aol::Max<short> ( aol::Abs<short> ( boundaryIndices.W ( ) - geometricCenter[0] ), aol::Abs<short> ( boundaryIndices.E ( ) - geometricCenter[0] ) ) );
    const short height = 1 + 2 * round ( 1.25 * aol::Max<short> ( aol::Abs<short> ( boundaryIndices.N ( ) - geometricCenter[1] ), aol::Abs<short> (  boundaryIndices.S ( ) - geometricCenter[1] ) ) );
    Centers[k][0] = geometricCenter[0];
    Centers[k][1] = geometricCenter[1];
    Dimensions[k][0] = width;
    Dimensions[k][1] = height;
    ++k;
  }
  Centers.resize ( k, 2 );
  Dimensions.resize ( k, 2 );
}


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
_ColoredPictureType>::getApproximateDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<int> &Dimensions,
                                                            aol::Vector<RealType> &DumbbellOrientations, aol::Vector<RealType> &DumbbellSeparations,
                                                            qc::BitArray<qc::QC_2D> &Segmented,
                                                            const PictureType &Data,
                                                            const RealType Gamma, const int MaxIt, const RealType Epsilon ) const {
  ComponentsCollection<RealType> componentsCollection;
  getComponentsCollection ( componentsCollection, Segmented, Data, Gamma, MaxIt, Epsilon );
  
  // Split components into dumbbells and single atoms and collect orientation and separation for dumbbells
  int numNonEmptyComps = componentsCollection.getNumNonEmptyComponents ( );
  Centers.reallocate ( numNonEmptyComps, 2 );
  Dimensions.reallocate ( numNonEmptyComps, 2 );
  DumbbellOrientations.reallocate ( numNonEmptyComps );
  DumbbellSeparations.reallocate ( numNonEmptyComps );
  int k = 0;
  for ( NonEmptyComponentsIterator it ( componentsCollection ); it.notAtEnd ( ) ; ++it ) {
    const RandomAccessStencil<short> boundaryIndices = (*it).getBoundaryIndices ( );
    aol::Vec2<short> geometricCenter = (*it).getGeometricCenter ( );
    const short width = 1 + 2 * round ( 1.25 * aol::Max<short> ( aol::Abs<short> ( boundaryIndices.W ( ) - geometricCenter[0] ), aol::Abs<short> ( boundaryIndices.E ( ) - geometricCenter[0] ) ) );
    const short height = 1 + 2 * round ( 1.25 * aol::Max<short> ( aol::Abs<short> ( boundaryIndices.N ( ) - geometricCenter[1] ), aol::Abs<short> (  boundaryIndices.S ( ) - geometricCenter[1] ) ) );
    aol::Vec<4,RealType> minMaxDiamAngle = (*it).getMinMaxDiameterAndAngle ( );
    Centers[k][0] = geometricCenter[0];
    Centers[k][1] = geometricCenter[1];
    Dimensions[k][0] = width;
    Dimensions[k][1] = height;
    DumbbellOrientations[k] = minMaxDiamAngle[1];
    DumbbellSeparations[k] = 0.5 * minMaxDiamAngle[0];
    ++k;
  }
  Centers.resize ( k, 2 );
  Dimensions.resize ( k, 2 );
  DumbbellOrientations.resize ( k );
  DumbbellSeparations.resize ( k );
  
  if ( _diskOutput ) {
    // Add orientation lines to the dumbbell atoms
    std::stringstream ss;
    ss << _outputDir << "/geometricCenters.png";
    ColoredPictureType geometricCentersImg ( ss.str ( ).c_str ( ) );
    aol::Vec2<short> pos;
    RealType maxVal = geometricCentersImg.getMaxValue ( );
    for ( int j=0; j<k ; ++j ) {
      for ( short i=-DumbbellSeparations[j]; i<=DumbbellSeparations[j] ; ++i ) {
        pos.set ( round ( Centers[j][0] + cos ( DumbbellOrientations[j] ) * i ), round ( Centers[j][1] + sin ( DumbbellOrientations[j] ) * i ) );
        if ( pos[0] >=0 && pos[0] < geometricCentersImg[0].getNumX ( ) && pos[1] >=0 && pos[1] < geometricCentersImg[0].getNumY ( ) ) {
          geometricCentersImg[0].set ( pos, 0 );
          geometricCentersImg[1].set ( pos, maxVal );
          geometricCentersImg[2].set ( pos, 0 );
        }
      }
      pos.set ( Centers[j][0], Centers[j][1] );
      geometricCentersImg[0].set ( pos, maxVal );
      geometricCentersImg[1].set ( pos, 0 );
    }
    geometricCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
    geometricCentersImg.savePNG ( ss.str ( ).c_str ( ) );
  }
}


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
                _ColoredPictureType>::getApproximateAtomPositions ( aol::MultiVector<RealType> &Centers, qc::BitArray<qc::QC_2D> &Segmented,
                                                                    const PictureType &Data,
                                                                    const RealType Gamma, const int MaxIt, const RealType Epsilon ) const {
  ComponentsCollection<RealType> componentsCollection;
  getComponentsCollection ( componentsCollection, Segmented, Data, Gamma, MaxIt, Epsilon );
  
  if ( !_quietMode )
    std::cerr << "# of atoms after segmentation: " << componentsCollection.getNumNonEmptyComponents ( ) << std::endl;
  
  if ( _diskOutput ) {
    InitType grid ( qc::GridSize<ConfiguratorType::Dim> ( Data.getSize ( ) ) );
    
    PictureType connectedRegions ( grid );
    componentsCollection.createPictureBWComponents ( connectedRegions );
    std::stringstream ss;
    ss << _outputDir << "/connectedRegions" << getDefaultArraySuffix ( qc::QC_2D );
    connectedRegions.save ( ss.str ( ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    
    ColoredPictureType geometricCentersImg ( grid );
    componentsCollection.createPictureBWComponentsRedGeometricCenters ( geometricCentersImg );
    ss.str ( std::string ( ) ); // clear stringstream
    ss << _outputDir << "/geometricCenters.png";
    geometricCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
    geometricCentersImg.savePNG ( ss.str ( ).c_str ( ) );
  }
  
  short k = 0;
  Centers.reallocate ( componentsCollection.getNumNonEmptyComponents ( ), 2 );
  for ( NonEmptyComponentsIterator it ( componentsCollection ); it.notAtEnd ( ) ; ++it ) {
    aol::Vec2<short> geometricCenter = (*it).getGeometricCenter ( );
    Centers[k][0] = geometricCenter[0];
    Centers[k][1] = geometricCenter[1];
    ++k;
  }
}


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
_ColoredPictureType>::getRefinedAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                                                const aol::MultiVector<RealType> &ApproximateCenters, const aol::MultiVector<int> &ApproximateDimensions,
                                                const PictureType &Data ) const {
  PictureType u ( Data );
  u /= u.getMaxValue ( );
  u.setOverflowHandling( aol::CLIP_THEN_SCALE, 0, 1 );
  
  if ( !_quietMode )
    std::cerr << "Refining single atom positions via 2D Gaussian fit.." << std::endl;
  
  Centers.reallocate ( ApproximateCenters.numComponents ( ), 2 );
  const short numVar = AsymmetricGaussianBumpFunction<RealType>::NumberOfParameters;
  GaussianParams.reallocate ( ApproximateCenters.numComponents ( ), numVar );
  PictureType bumpFunctionImg ( u.getNumX ( ), u.getNumY ( ) );
  RealType approxAtomRadius = 0;
  if ( _progressBar != NULL ) _progressBar->start ( ApproximateCenters.numComponents ( ) );
  for ( short k=0; k<ApproximateCenters.numComponents ( ) && !wantsInterrupt ( ) ; ++k ) {
    aol::Vec2<short> maxAtomOffset ( ( ApproximateDimensions[k][0] - 1 ) / 2, ( ApproximateDimensions[k][1] - 1 ) / 2 ), atomOffset, atomSize, approxCenterInt;
    approxCenterInt.set ( round ( ApproximateCenters[k][0] ), round ( ApproximateCenters[k][1] ) );
    atomOffset.set ( aol::Min<short> ( aol::Min<short> ( maxAtomOffset[0], approxCenterInt[0] ), u.getNumX ( ) - 1 - approxCenterInt[0] ),
                    aol::Min<short> ( aol::Min<short> ( maxAtomOffset[1], approxCenterInt[1] ), u.getNumY ( ) - 1 - approxCenterInt[1] ) );
    atomSize.set ( 2 * atomOffset[0] + 1, 2 * atomOffset[1] + 1 );
    aol::Vec2<short> corner ( approxCenterInt[0] - atomOffset[0], approxCenterInt[1] - atomOffset[1] );
    PictureType atomData ( atomSize[0], atomSize[1] );
    for ( short dx=0; dx<atomSize[0] ; ++dx )
      for ( short dy=0; dy<atomSize[1] ; ++dy )
        atomData.set ( dx, dy, u.get ( corner[0] + dx, corner[1] + dy ) );
    SingleAsymmetricBumpFitTargetFunctional<RealType> gaussNewtonF ( atomData );
    SingleAsymmetricBumpFitTargetJacobian<RealType> gaussNewtonDF ( atomData );
    aol::Vector<RealType> lowerBounds ( numVar );
    aol::Vector<RealType> upperBounds ( numVar );
    lowerBounds[0] = 0;
    upperBounds[0] = atomSize[0];
    lowerBounds[1] = 0;
    upperBounds[1] = atomSize[1];
    aol::BitVector constrainedDirections ( numVar );
    constrainedDirections.set ( 0, true );
    constrainedDirections.set ( 1, true );
    ProjectorType boxProjector ( lowerBounds, upperBounds, constrainedDirections );
    ConstrainedLevenbergMarquardtAlgorithm<RealType, MatrixType, LinearRegressionType, ProjectorType> levenbergMarquardtAlg ( atomData.size ( ), gaussNewtonF, gaussNewtonDF, boxProjector, 50, 1, 0.2, 0.8, 1e-6, 1e-6, 1e-6, false );
    aol::Vector<RealType> Arg ( numVar ), Dest ( numVar );
    Arg[0] = 0.5 * atomSize[0];   // Center x
    Arg[1] = 0.5 * atomSize[1];   // Center y
    Arg[2] = 0.5;                 // Intensity
    Arg[3] = 0.5 * atomSize[0];   // Width
    Arg[4] = 0.5 * atomSize[1];   // Height
    Arg[5] = 0;                   // Rotation
    Arg[6] = 0;                   // Offset
    levenbergMarquardtAlg.apply ( Arg, Dest );
    Dest[0] += corner[0];
    Dest[1] += corner[1];
    Centers[k][0] = Dest[0];
    Centers[k][1] = Dest[1];
    
    for ( int j=0; j<numVar ; ++j )
      GaussianParams[k][j] = Dest[j];
    
    approxAtomRadius += Dest[3] + Dest[4];
    
    if ( _diskOutput ) {
      Dest[6] = 0;
      const AsymmetricGaussianBumpFunction<RealType> bumpFunc ( Dest );
      for ( int x=corner[0]; x<corner[0]+atomSize[0] ; ++x )
        for ( int y=corner[1]; y<corner[1]+atomSize[1] ; ++y )
          bumpFunctionImg.add ( x, y, bumpFunc.evaluate ( aol::Vec2<RealType> ( x, y ) ) );
    }
    if ( _progressBar != NULL ) (*_progressBar)++;
  }
  if ( _progressBar != NULL ) _progressBar->finish ( );
  
  approxAtomRadius /= 2 * ApproximateCenters.numComponents ( );
  
  // Clear boundary atoms
  short k = 0;
  while ( k < Centers.numComponents ( ) ) {
    if ( Centers[k][0] - approxAtomRadius < 0 || Centers[k][0] + approxAtomRadius >= Data.getNumX ( )
        || Centers[k][1] - approxAtomRadius < 0 || Centers[k][1] + approxAtomRadius >= Data.getNumY ( ) ) {
      Centers.eraseComponent ( k );
      GaussianParams.eraseComponent ( k );
    } else ++k;
  }
  
  if ( !_quietMode )
    std::cerr << "# of (inner) atoms after bump fitting: " << Centers.numComponents ( ) << std::endl;
  
  if ( _diskOutput ) {
    std::stringstream ss;
    ss << _outputDir << "/bumpFunctions_single" << getDefaultArraySuffix ( qc::QC_2D );
    bumpFunctionImg.save ( ss.str ( ).c_str ( ), qc::PGM_DOUBLE_BINARY );
  }
  
  if ( _diskOutput ) {
    InitType grid ( qc::GridSize<ConfiguratorType::Dim> ( u.getSize ( ) ) );
    ColoredPictureType atomicCentersImg ( grid );
    for ( short i=0; i<3 ; ++ i ) atomicCentersImg[i] = u;
    aol::Vec2<short> pos;
    for ( int i=0; i<Centers.numComponents ( ) ; ++i ) {
      pos.set ( round ( Centers[i][0] ), round ( Centers[i][1] ) );
      atomicCentersImg[0].set ( pos, 1 );
      atomicCentersImg[1].set ( pos, 0 );
      atomicCentersImg[2].set ( pos, 0 );
    }
    std::stringstream ss;
    ss << _outputDir << "/atomicCenters_single.png";
    atomicCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
    atomicCentersImg.savePNG ( ss.str ( ).c_str ( ) );
  }
}


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
                _ColoredPictureType>::getRefinedAtomPositions ( aol::MultiVector<RealType> &Centers, const aol::MultiVector<RealType> &ApproximateCenters, aol::MultiVector<RealType> &GaussianParams,
                                                                const PictureType &Data, const aol::Vec2<short> &AtomSize ) const {
  PictureType u ( Data );
  u /= u.getMaxValue ( );
  u.setOverflowHandling( aol::CLIP_THEN_SCALE, 0, 1 );
                  
  if ( !_quietMode )
    std::cerr << "Refining atom positions via 2D Gaussian fit.." << std::endl;
                  
  Centers.reallocate ( ApproximateCenters.numComponents ( ), 2 );
  const short numVar = AsymmetricGaussianBumpFunction<RealType>::NumberOfParameters;
  GaussianParams.reallocate ( ApproximateCenters.numComponents ( ), numVar );
  PictureType bumpFunctionImg ( u.getNumX ( ), u.getNumY ( ) );
  aol::Vec2<short> maxAtomOffset ( ( AtomSize[0] - 1 ) / 2, ( AtomSize[1] - 1 ) / 2 ), atomOffset, atomSize, approxCenterInt;
  RealType approxAtomRadius = 0;
  for ( short k=0; k<ApproximateCenters.numComponents ( ) ; ++k ) {
    approxCenterInt.set ( round ( ApproximateCenters[k][0] ), round ( ApproximateCenters[k][1] ) );
    atomOffset.set ( aol::Min<short> ( aol::Min<short> ( maxAtomOffset[0], approxCenterInt[0] ), u.getNumX ( ) - 1 - approxCenterInt[0] ),
                    aol::Min<short> ( aol::Min<short> ( maxAtomOffset[1], approxCenterInt[1] ), u.getNumY ( ) - 1 - approxCenterInt[1] ) );
    atomSize.set ( 2 * atomOffset[0] + 1, 2 * atomOffset[1] + 1 );
    aol::Vec2<short> corner ( approxCenterInt[0] - atomOffset[0], approxCenterInt[1] - atomOffset[1] );
    PictureType atomData ( atomSize[0], atomSize[1] );
    for ( short dx=0; dx<atomSize[0] ; ++dx )
      for ( short dy=0; dy<atomSize[1] ; ++dy )
        atomData.set ( dx, dy, u.get ( corner[0] + dx, corner[1] + dy ) );
    SingleAsymmetricBumpFitTargetFunctional<RealType> gaussNewtonF ( atomData );
    SingleAsymmetricBumpFitTargetJacobian<RealType> gaussNewtonDF ( atomData );
    aol::Vector<RealType> lowerBounds ( numVar );
    aol::Vector<RealType> upperBounds ( numVar );
    lowerBounds[0] = 0;
    upperBounds[0] = atomSize[0];
    lowerBounds[1] = 0;
    upperBounds[1] = atomSize[1];
    aol::BitVector constrainedDirections ( numVar );
    constrainedDirections.set ( 0, true );
    constrainedDirections.set ( 1, true );
    ProjectorType boxProjector ( lowerBounds, upperBounds, constrainedDirections );
    ConstrainedLevenbergMarquardtAlgorithm<RealType, MatrixType, LinearRegressionType, ProjectorType> levenbergMarquardtAlg ( atomData.size ( ), gaussNewtonF, gaussNewtonDF, boxProjector, 50, 1, 0.2, 0.8, 1e-6, 1e-6, 1e-6, false );
    aol::Vector<RealType> Arg ( numVar ), Dest ( numVar );
    Arg[0] = 0.5 * atomSize[0];   // Center x
    Arg[1] = 0.5 * atomSize[1];   // Center y
    Arg[2] = 0.5;                 // Intensity
    Arg[3] = 0.5 * atomSize[0];   // Width
    Arg[4] = 0.5 * atomSize[1];   // Height
    Arg[5] = 0;                   // Rotation
    Arg[6] = 0;                   // Offset
    levenbergMarquardtAlg.apply ( Arg, Dest );
    Dest[0] += corner[0];
    Dest[1] += corner[1];
    Centers[k][0] = Dest[0];
    Centers[k][1] = Dest[1];
    
    for ( int j=0; j<numVar ; ++j )
      GaussianParams[k][j] = Dest[j];
    
    approxAtomRadius += Dest[3] + Dest[4];
    
    if ( _diskOutput ) {
      Dest[6] = 0;
      const AsymmetricGaussianBumpFunction<RealType> bumpFunc ( Dest );
      for ( int x=corner[0]; x<corner[0]+atomSize[0] ; ++x )
        for ( int y=corner[1]; y<corner[1]+atomSize[1] ; ++y )
          bumpFunctionImg.add ( x, y, bumpFunc.evaluate ( aol::Vec2<RealType> ( x, y ) ) );
    }
  }
                  
  approxAtomRadius /= 2 * ApproximateCenters.numComponents ( );
  
  // Clear boundary atoms
  short k = 0;
  while ( k < Centers.numComponents ( ) ) {
    if ( Centers[k][0] - approxAtomRadius < 0 || Centers[k][0] + approxAtomRadius >= Data.getNumX ( )
        || Centers[k][1] - approxAtomRadius < 0 || Centers[k][1] + approxAtomRadius >= Data.getNumY ( ) ) {
      Centers.eraseComponent ( k );
      GaussianParams.eraseComponent ( k );
    } else ++k;
  }
                  
  if ( !_quietMode )
    std::cerr << "# of atoms after bump fitting: " << Centers.numComponents ( ) << std::endl;
  
  if ( _diskOutput ) {
    std::stringstream ss;
    ss << _outputDir << "/bumpFunctions" << getDefaultArraySuffix ( qc::QC_2D );
    bumpFunctionImg.save ( ss.str ( ).c_str ( ), qc::PGM_DOUBLE_BINARY );
  }
  
  if ( _diskOutput ) {
    InitType grid ( qc::GridSize<ConfiguratorType::Dim> ( u.getSize ( ) ) );
    ColoredPictureType atomicCentersImg ( grid );
    for ( short i=0; i<3 ; ++ i ) atomicCentersImg[i] = u;
    aol::Vec2<short> pos;
    for ( int i=0; i<Centers.numComponents ( ) ; ++i ) {
      pos.set ( round ( Centers[i][0] ), round ( Centers[i][1] ) );
      atomicCentersImg[0].set ( pos, 1 );
      atomicCentersImg[1].set ( pos, 0 );
      atomicCentersImg[2].set ( pos, 0 );
    }
    std::stringstream ss;
    ss << _outputDir << "/atomicCenters.png";
    atomicCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
    atomicCentersImg.savePNG ( ss.str ( ).c_str ( ) );
  }
}


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
_ColoredPictureType>::getRefinedDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                                                        const aol::MultiVector<RealType> &ApproximateCenters, const aol::MultiVector<int> &ApproximateDimensions,
                                                        const aol::Vector<RealType> &Orientations, const aol::Vector<RealType> &Separations,
                                                        const PictureType &Data ) const {
  PictureType u ( Data );
  u /= u.getMaxValue ( );
  u.setOverflowHandling( aol::CLIP_THEN_SCALE, 0, 1 );
  
  if ( !_quietMode )
    std::cerr << "Refining positions within dumbbells via 2D Gaussian fit.." << std::endl;
  
  Centers.reallocate ( 2 * ApproximateCenters.numComponents ( ), 2 );
  const short numVar = AsymmetricGaussianDoubleBumpFunction<RealType>::NumberOfParameters;
  GaussianParams.reallocate ( ApproximateCenters.numComponents ( ), numVar );
  PictureType bumpFunctionImg ( u.getNumX ( ), u.getNumY ( ) );
  RealType approxAtomRadius = 0;
  if ( _progressBar != NULL ) _progressBar->start ( ApproximateCenters.numComponents ( ) );
  for ( short k=0; k<ApproximateCenters.numComponents ( ) && !wantsInterrupt ( ) ; ++k ) {
    const RealType atomAngleRads = Orientations[k], atomSeparation = Separations[k];
    aol::Vec2<short> maxAtomOffset ( ( ApproximateDimensions[k][0] - 1 ) / 2, ( ApproximateDimensions[k][1] - 1 ) / 2 ), atomOffset, atomSize, approxCenterInt;
    approxCenterInt.set ( round ( ApproximateCenters[k][0] ), round ( ApproximateCenters[k][1] ) );
    atomOffset.set ( aol::Min<short> ( aol::Min<short> ( maxAtomOffset[0], approxCenterInt[0] ), u.getNumX ( ) - 1 - approxCenterInt[0] ),
                    aol::Min<short> ( aol::Min<short> ( maxAtomOffset[1], approxCenterInt[1] ), u.getNumY ( ) - 1 - approxCenterInt[1] ) );
    atomSize.set ( 2 * atomOffset[0] + 1, 2 * atomOffset[1] + 1 );
    aol::Vec2<short> corner ( approxCenterInt[0] - atomOffset[0], approxCenterInt[1] - atomOffset[1] );
    PictureType atomData ( atomSize[0], atomSize[1] );
    for ( short dx=0; dx<atomSize[0] ; ++dx )
      for ( short dy=0; dy<atomSize[1] ; ++dy )
        atomData.set ( dx, dy, u.get ( corner[0] + dx, corner[1] + dy ) );
    SingleAsymmetricDoubleBumpFitTargetFunctional<RealType> gaussNewtonF ( atomData );
    SingleAsymmetricDoubleBumpFitTargetJacobian<RealType> gaussNewtonDF ( atomData );
    aol::Vector<RealType> lowerBounds ( numVar );
    aol::Vector<RealType> upperBounds ( numVar );
    lowerBounds[0] = 0;
    upperBounds[0] = atomSize[0];
    lowerBounds[1] = 0;
    upperBounds[1] = atomSize[1];
    lowerBounds[2] = 0;
    upperBounds[2] = atomSize[0];
    lowerBounds[3] = 0;
    upperBounds[3] = atomSize[1];
    aol::BitVector constrainedDirections ( numVar );
    constrainedDirections.set ( 0, true );
    constrainedDirections.set ( 1, true );
    constrainedDirections.set ( 2, true );
    constrainedDirections.set ( 3, true );
    ProjectorType boxProjector ( lowerBounds, upperBounds, constrainedDirections );
    ConstrainedLevenbergMarquardtAlgorithm<RealType, MatrixType, LinearRegressionType, ProjectorType> levenbergMarquardtAlg ( atomData.size ( ), gaussNewtonF, gaussNewtonDF, boxProjector, 50, 1, 0.2, 0.8, 1e-6, 1e-6, 1e-6, false );
    aol::Vector<RealType> Arg ( numVar ), Dest ( numVar );
    Arg[0] = 0.5 * atomSize[0] - std::cos ( atomAngleRads ) * 0.5 * atomSeparation;    // Center 1 x
    Arg[1] = 0.5 * atomSize[1] - std::sin ( atomAngleRads ) * 0.5 * atomSeparation;    // Center 1 y
    Arg[2] = 0.5 * atomSize[0] + std::cos ( atomAngleRads ) * 0.5 * atomSeparation;    // Center 2 x
    Arg[3] = 0.5 * atomSize[1] + std::sin ( atomAngleRads ) * 0.5 * atomSeparation;    // Center 2 y
    Arg[4] = 0.5;                   // Intensity 1
    Arg[5] = 0.5;                   // Intensity 2
    Arg[6] = 0.25 * atomSize.getMinValue ( );     // Width 1
    Arg[7] = 0.25 * atomSize.getMinValue ( );     // Height 1
    Arg[8] = 0.25 * atomSize.getMinValue ( );     // Width 2
    Arg[9] = 0.25 * atomSize.getMinValue ( );     // Height 2
    Arg[10] = 0;                    // Rotation 1
    Arg[11] = 0;                    // Rotation 2
    Arg[12] = 0;                    // Offset
    levenbergMarquardtAlg.apply ( Arg, Dest );
    Dest[0] += corner[0];
    Dest[1] += corner[1];
    Dest[2] += corner[0];
    Dest[3] += corner[1];
    Centers[2*k][0] = Dest[0];
    Centers[2*k][1] = Dest[1];
    Centers[2*k+1][0] = Dest[2];
    Centers[2*k+1][1] = Dest[3];
    
    for ( int j=0; j<numVar ; ++j )
      GaussianParams[k][j] = Dest[j];
    
    approxAtomRadius += Dest[6] + Dest[7] + Dest[8] + Dest[9];
    
    if ( _diskOutput ) {
      Dest[12] = 0;
      const AsymmetricGaussianDoubleBumpFunction<RealType> bumpFunc ( Dest );
      for ( int x=corner[0]; x<corner[0]+atomSize[0] ; ++x )
        for ( int y=corner[1]; y<corner[1]+atomSize[1] ; ++y )
          bumpFunctionImg.add ( x, y, bumpFunc.evaluate ( aol::Vec2<RealType> ( x, y ) ) );
    }
    
    if ( _progressBar != NULL ) (*_progressBar)++;
  }
  if ( _progressBar != NULL ) _progressBar->finish ( );
  
  approxAtomRadius /= 4 * ApproximateCenters.numComponents ( );
  // Clear boundary atoms
  short k = 0;
  while ( k < Centers.numComponents ( ) - 1 ) {
    if ( Centers[k][0] - approxAtomRadius < 0 || Centers[k][0] + approxAtomRadius >= Data.getNumX ( )
        || Centers[k][1] - approxAtomRadius < 0 || Centers[k][1] + approxAtomRadius >= Data.getNumY ( )
        || Centers[k+1][0] - approxAtomRadius < 0 || Centers[k+1][0] + approxAtomRadius >= Data.getNumX ( )
        || Centers[k+1][1] - approxAtomRadius < 0 || Centers[k+1][1] + approxAtomRadius >= Data.getNumY ( ) ) {
      Centers.eraseComponent ( k );
      Centers.eraseComponent ( k ); // component k+1 is now component k
      GaussianParams.eraseComponent ( k / 2 );
    } else k += 2;
  }
  
  if ( !_quietMode )
    std::cerr << "# of (inner) dumbbells after bump fitting: " << Centers.numComponents ( ) / 2 << std::endl;
  
  if ( _diskOutput ) {
    std::stringstream ss;
    ss << _outputDir << "/bumpFunctions_dumbbell" << getDefaultArraySuffix ( qc::QC_2D );
    bumpFunctionImg.save ( ss.str ( ).c_str ( ), qc::PGM_DOUBLE_BINARY );
  }
  
  if ( _diskOutput ) {
    InitType grid ( qc::GridSize<ConfiguratorType::Dim> ( u.getSize ( ) ) );
    ColoredPictureType atomicCentersImg ( grid );
    for ( short i=0; i<3 ; ++ i ) atomicCentersImg[i] = u;
    aol::Vec2<short> pos;
    for ( int i=0; i<Centers.numComponents ( ) ; ++i ) {
      pos.set ( round ( Centers[i][0] ), round ( Centers[i][1] ) );
      atomicCentersImg[0].set ( pos, 1 );
      atomicCentersImg[1].set ( pos, 0 );
      atomicCentersImg[2].set ( pos, 0 );
    }
    std::stringstream ss;
    ss << _outputDir << "/atomicCenters_dumbbell.png";
    atomicCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
    atomicCentersImg.savePNG ( ss.str ( ).c_str ( ) );
  }
}


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
_ColoredPictureType>::getRefinedDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers, const aol::MultiVector<RealType> &ApproximateCenters, aol::MultiVector<RealType> &GaussianParams,
                                                        const PictureType &Data, const aol::Vec2<short> &AtomSize,
                                                        const RealType AtomSeparation, const RealType AtomAngle ) const {
  PictureType u ( Data );
  u /= u.getMaxValue ( );
  u.setOverflowHandling( aol::CLIP_THEN_SCALE, 0, 1 );
  
  if ( !_quietMode )
    std::cerr << "Refining atom positions via 2D Gaussian fit.." << std::endl;
  
  Centers.reallocate ( 2 * ApproximateCenters.numComponents ( ), 2 );
  const short numVar = AsymmetricGaussianDoubleBumpFunction<RealType>::NumberOfParameters;
  GaussianParams.reallocate ( ApproximateCenters.numComponents ( ), numVar );
  PictureType bumpFunctionImg ( u.getNumX ( ), u.getNumY ( ) );
  const RealType atomAngleRads = AtomAngle * aol::NumberTrait<RealType>::pi / 180;
  aol::Vec2<short> maxAtomOffset ( ( AtomSize[0] - 1 ) / 2, ( AtomSize[1] - 1 ) / 2 ), atomOffset, atomSize, approxCenterInt;
  RealType approxAtomRadius = 0;
  for ( short k=0; k<ApproximateCenters.numComponents ( ) ; ++k ) {
    approxCenterInt.set ( round ( ApproximateCenters[k][0] ), round ( ApproximateCenters[k][1] ) );
    atomOffset.set ( aol::Min<short> ( aol::Min<short> ( maxAtomOffset[0], approxCenterInt[0] ), u.getNumX ( ) - 1 - approxCenterInt[0] ),
                     aol::Min<short> ( aol::Min<short> ( maxAtomOffset[1], approxCenterInt[1] ), u.getNumY ( ) - 1 - approxCenterInt[1] ) );
    atomSize.set ( 2 * atomOffset[0] + 1, 2 * atomOffset[1] + 1 );
    aol::Vec2<short> corner ( approxCenterInt[0] - atomOffset[0], approxCenterInt[1] - atomOffset[1] );
    PictureType atomData ( atomSize[0], atomSize[1] );
    for ( short dx=0; dx<atomSize[0] ; ++dx )
      for ( short dy=0; dy<atomSize[1] ; ++dy )
        atomData.set ( dx, dy, u.get ( corner[0] + dx, corner[1] + dy ) );
    SingleAsymmetricDoubleBumpFitTargetFunctional<RealType> gaussNewtonF ( atomData );
    SingleAsymmetricDoubleBumpFitTargetJacobian<RealType> gaussNewtonDF ( atomData );
    aol::Vector<RealType> lowerBounds ( numVar );
    aol::Vector<RealType> upperBounds ( numVar );
    lowerBounds[0] = 0;
    upperBounds[0] = atomSize[0];
    lowerBounds[1] = 0;
    upperBounds[1] = atomSize[1];
    lowerBounds[2] = 0;
    upperBounds[2] = atomSize[0];
    lowerBounds[3] = 0;
    upperBounds[3] = atomSize[1];
    aol::BitVector constrainedDirections ( numVar );
    constrainedDirections.set ( 0, true );
    constrainedDirections.set ( 1, true );
    constrainedDirections.set ( 2, true );
    constrainedDirections.set ( 3, true );
    ProjectorType boxProjector ( lowerBounds, upperBounds, constrainedDirections );
    ConstrainedLevenbergMarquardtAlgorithm<RealType, MatrixType, LinearRegressionType, ProjectorType> levenbergMarquardtAlg ( atomData.size ( ), gaussNewtonF, gaussNewtonDF, boxProjector, 50, 1, 0.2, 0.8, 1e-6, 1e-6, 1e-6, false );
    aol::Vector<RealType> Arg ( numVar ), Dest ( numVar );
    Arg[0] = 0.5 * atomSize[0] - std::cos ( atomAngleRads ) * 0.5 * AtomSeparation;    // Center 1 x
    Arg[1] = 0.5 * atomSize[1] - std::sin ( atomAngleRads ) * 0.5 * AtomSeparation;    // Center 1 y
    Arg[2] = 0.5 * atomSize[0] + std::cos ( atomAngleRads ) * 0.5 * AtomSeparation;    // Center 2 x
    Arg[3] = 0.5 * atomSize[1] + std::sin ( atomAngleRads ) * 0.5 * AtomSeparation;    // Center 2 y
    Arg[4] = 0.5;                   // Intensity 1
    Arg[5] = 0.5;                   // Intensity 2
    Arg[6] = 0.25 * atomSize.getMinValue ( );     // Width 1
    Arg[7] = 0.25 * atomSize.getMinValue ( );     // Height 1
    Arg[8] = 0.25 * atomSize.getMinValue ( );     // Width 2
    Arg[9] = 0.25 * atomSize.getMinValue ( );     // Height 2
    Arg[10] = 0;                    // Rotation 1
    Arg[11] = 0;                    // Rotation 2
    Arg[12] = 0;                    // Offset
    levenbergMarquardtAlg.apply ( Arg, Dest );
    Dest[0] += corner[0];
    Dest[1] += corner[1];
    Dest[2] += corner[0];
    Dest[3] += corner[1];
    Centers[2*k][0] = Dest[0];
    Centers[2*k][1] = Dest[1];
    Centers[2*k+1][0] = Dest[2];
    Centers[2*k+1][1] = Dest[3];
    
    for ( int j=0; j<numVar ; ++j )
      GaussianParams[k][j] = Dest[j];
    
    approxAtomRadius += Dest[6] + Dest[7] + Dest[8] + Dest[9];
    
    if ( _diskOutput ) {
      Dest[12] = 0;
      const AsymmetricGaussianDoubleBumpFunction<RealType> bumpFunc ( Dest );
      for ( int x=corner[0]; x<corner[0]+atomSize[0] ; ++x )
        for ( int y=corner[1]; y<corner[1]+atomSize[1] ; ++y )
          bumpFunctionImg.add ( x, y, bumpFunc.evaluate ( aol::Vec2<RealType> ( x, y ) ) );
    }
  }
  
  approxAtomRadius /= 4 * ApproximateCenters.numComponents ( );
  // Clear boundary atoms
  short k = 0;
  while ( k < Centers.numComponents ( ) - 1 ) {
    if ( Centers[k][0] - approxAtomRadius < 0 || Centers[k][0] + approxAtomRadius >= Data.getNumX ( )
      || Centers[k][1] - approxAtomRadius < 0 || Centers[k][1] + approxAtomRadius >= Data.getNumY ( )
      || Centers[k+1][0] - approxAtomRadius < 0 || Centers[k+1][0] + approxAtomRadius >= Data.getNumX ( )
      || Centers[k+1][1] - approxAtomRadius < 0 || Centers[k+1][1] + approxAtomRadius >= Data.getNumY ( ) ) {
      Centers.eraseComponent ( k );
      Centers.eraseComponent ( k ); // component k+1 is now component k
      GaussianParams.eraseComponent ( k / 2 );
    } else k += 2;
  }
  
  if ( !_quietMode )
    std::cerr << "# of (dumbbell-) atoms after bump fitting: " << Centers.numComponents ( ) / 2 << std::endl;
  
  if ( _diskOutput ) {
    std::stringstream ss;
    ss << _outputDir << "/bumpFunctions" << getDefaultArraySuffix ( qc::QC_2D );
    bumpFunctionImg.save ( ss.str ( ).c_str ( ), qc::PGM_DOUBLE_BINARY );
  }
  
  if ( _diskOutput ) {
    InitType grid ( qc::GridSize<ConfiguratorType::Dim> ( u.getSize ( ) ) );
    ColoredPictureType atomicCentersImg ( grid );
    for ( short i=0; i<3 ; ++ i ) atomicCentersImg[i] = u;
    aol::Vec2<short> pos;
    for ( int i=0; i<Centers.numComponents ( ) ; ++i ) {
      pos.set ( round ( Centers[i][0] ), round ( Centers[i][1] ) );
      atomicCentersImg[0].set ( pos, 1 );
      atomicCentersImg[1].set ( pos, 0 );
      atomicCentersImg[2].set ( pos, 0 );
    }
    std::stringstream ss;
    ss << _outputDir << "/atomicCenters.png";
    atomicCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
    atomicCentersImg.savePNG ( ss.str ( ).c_str ( ) );
  }
}


template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
void AtomFinder<_RealType, _MatrixType, _LinearRegressionType, _ScalarPictureType,
                _ColoredPictureType>::readCentersFromCSV ( aol::MultiVector<RealType> &Centers, const std::string &Path ) const {
  std::ifstream file ( Path.c_str() );
  std::string value;
  std::list<std::string> values;
  while ( file.good() ) {
    getline ( file, value, ',' );
    if ( value.find('\n') != std::string::npos ) {
      size_t pos = 0;
      while ( ( pos = value.find ( "\n", pos + 1 ) ) != std::string::npos ) {
        std::string p = value.substr(0, pos);
        values.push_back(p);
        value = value.substr(pos + 1);
      }
      if (!value.empty()) {
        values.push_back(value);
      }
    } else {
      values.push_back(value);
    }
  }
  for ( int i=0; i<2; ++i ) // skip header
    values.pop_front ( );
  RealType xPos = 0, yPos = 0, xMax = 0, yMax = 0;
  int k = 0;
  for ( std::list<std::string>::const_iterator it = values.begin ( ); it != values.end ( ) ; it++ ) {
    const RealType val = strtod ( (*it).c_str(), NULL );
    if ( k % 2 == 0 ) {
      xPos = val;
      if ( xPos > xMax )
        xMax = xPos;
    }
    if ( k % 2 == 1 ) {
      yPos = val;
      if ( yPos > yMax )
        yMax = yPos;
      
      if ( !_quietMode )
        std::cerr << xPos << "; " << yPos << std::endl;
      
      Centers.resize ( Centers.numComponents ( ) + 1, 2 );
      Centers[Centers.numComponents ( )-1][0] = xPos;
      Centers[Centers.numComponents ( )-1][1] = yPos;
    }
    ++k;
  }
  file.close ( );
  
  if ( _diskOutput ) {
    short size = aol::Max ( round ( xMax * 1.1 ), round ( yMax * 1.1 ) );
    ColoredPictureType atomicCentersImg ( size, size );
    aol::Vec2<short> pos;
    for ( int i=0; i<Centers.numComponents ( ) ; ++i ) {
      pos.set ( round ( Centers[i][0] ), round ( Centers[i][1] ) );
      atomicCentersImg[0].set ( pos, 1 );
      atomicCentersImg[1].set ( pos, 0 );
      atomicCentersImg[2].set ( pos, 0 );
    }
    std::stringstream ss;
    ss << _outputDir << "/atomicCenters.png";
    atomicCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
    atomicCentersImg.savePNG ( ss.str ( ).c_str ( ) );
  }
}

template class AtomFinder<double, aol::FullMatrix<double>, LinearRegressionQR<double>,
                          qc::ScalarArray<double, qc::QC_2D>, qc::MultiArray<double, qc::QC_2D, 3> >;