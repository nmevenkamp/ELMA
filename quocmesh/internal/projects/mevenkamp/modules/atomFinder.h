#ifndef ATOMFINDER_H_
#define ATOMFINDER_H_


// standard
#include <iostream>
#include <fstream>
#include <cmath>

// quocmesh
#include <aol.h>
#include "component.h"
#include "componentsCollection.h"
#include <configurators.h>
#include <finiteDifferences.h>
#include <firstOrderTVAlgos.h>
#include <imageTools.h>
#include <multiArray.h>
#include "scalarArrayExtensions.h"
#include "nonLocalMeansFilter.h"
#include <parameterParser.h>
#include <scalarArray.h>
#include "bumpFit.h"
#include "linearRegression.h"
#include "nonLinearRegression.h"

const static double GAMMA = 1e-4, EPSILON = 1e-3;
const static int MAXIT = 10000;

template <typename _RealType, typename _MatrixType, typename _LinearRegressionType, typename _ScalarPictureType, typename _ColoredPictureType>
class AtomFinder {
  typedef _RealType RealType;
  typedef _MatrixType MatrixType;
  typedef _LinearRegressionType LinearRegressionType;
  typedef _ScalarPictureType PictureType;
  typedef _ColoredPictureType ColoredPictureType;
  typedef qc::MultiArray<RealType, qc::QC_2D, 1> ArrayType;
  typedef qc::RectangularGridConfigurator<RealType, qc::QC_2D, aol::GaussQuadrature<RealType, qc::QC_2D, 3> > ConfiguratorType;
  typedef typename ConfiguratorType::InitType InitType;
  typedef typename ComponentsCollection<RealType>::NonEmptyComponentsIterator NonEmptyComponentsIterator;
  typedef BoxProjector<RealType, aol::Vector<RealType>, MatrixType> ProjectorType;
protected:
  std::string _outputDir;
  bool _quietMode, _diskOutput;
  aol::ProgressBar<> *_progressBar;
  bool _catchCtrlC;
  mutable sigfunc _previousCtrlCHandler;
public:
  AtomFinder ( const std::string &OutputDir = "", const bool Quiet = true )
    : _outputDir ( OutputDir ), _quietMode ( Quiet ), _diskOutput ( OutputDir != "" ), _progressBar ( NULL ), _catchCtrlC ( false ) { }
  
  void getAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams, qc::BitArray<qc::QC_2D> &Segmented,
                          const PictureType &Data,
                          const RealType Gamma = GAMMA, const int MaxIt = MAXIT, const RealType Epsilon = EPSILON ) {
    setCtrlCHandler ( );
    aol::MultiVector<RealType> approximateCenters, approximateDumbbellCenters;
    aol::MultiVector<int> approximateDimensions, approximateDumbbellDimensions;
    aol::Vector<RealType> dumbbellOrientations, dumbbellSeparations;
    if ( _progressBar != NULL ) _progressBar->setText ( "AtomFinder: segmenting and finding approximate centers (step 1/3)" );
    getApproximateAtomPositions ( approximateCenters, approximateDimensions, approximateDumbbellCenters, approximateDumbbellDimensions, dumbbellOrientations, dumbbellSeparations,
                                  Segmented, Data, Gamma, MaxIt, Epsilon );
    aol::MultiVector<RealType> *dumbbellCenters = new aol::MultiVector<RealType> ( 0, 0 ), *dumbbellGaussianParams = new aol::MultiVector<RealType> ( 0, 0 );
    if ( approximateCenters.numComponents ( ) > 0 ) {
      if ( _progressBar != NULL ) _progressBar->setText ( "AtomFinder: Gaussian fitting single atoms (step 2/3)" );
      getRefinedAtomPositions ( Centers, GaussianParams, approximateCenters, approximateDimensions, Data );
    }
    if ( approximateDumbbellCenters.numComponents ( ) > 0 ) {
      if ( _progressBar != NULL ) _progressBar->setText ( "AtomFinder: Gaussian fitting dumbbells (step 3/3)" );
      getRefinedDumbbellAtomPositions ( *dumbbellCenters, *dumbbellGaussianParams, approximateDumbbellCenters, approximateDumbbellDimensions, dumbbellOrientations, dumbbellSeparations, Data );
    }
    Centers.appendReference ( *dumbbellCenters, true );
    GaussianParams.appendReference ( *dumbbellGaussianParams, true );
    
    if ( _diskOutput ) {
      PictureType bumpFunctionImg ( Data.getNumX ( ), Data.getNumY ( ) );
      getBumpFunctionImage ( GaussianParams, bumpFunctionImg );
      std::stringstream ss;
      ss << _outputDir << "/bumpFunctions" << getDefaultArraySuffix ( qc::QC_2D );
      bumpFunctionImg.save ( ss.str ( ).c_str ( ), qc::PGM_DOUBLE_BINARY );
      
      ColoredPictureType atomicCentersImg ( Data.getNumX ( ), Data.getNumY ( ) );
      for ( short c=0; c<3 ; ++c ) atomicCentersImg[c] = Data;
      atomicCentersImg.scaleValuesTo01 ( );
      aol::Vec2<short> pos;
      for ( int i=0; i<Centers.numComponents ( ) ; ++i ) {
        pos.set ( round ( Centers[i][0] ), round ( Centers[i][1] ) );
        atomicCentersImg[0].set ( pos, 1 );
        atomicCentersImg[1].set ( pos, 0 );
        atomicCentersImg[2].set ( pos, 0 );
      }
      ss.str ( std::string ( ) ); // clear stringstream
      ss << _outputDir << "/atomicCenters.png";
      atomicCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
      atomicCentersImg.savePNG ( ss.str ( ).c_str ( ) );
    }
    unsetCtrlCHandler ( );
  }
  
  void getSingleAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams, qc::BitArray<qc::QC_2D> &Segmented,
                                 const PictureType &Data,
                                 const RealType Gamma = GAMMA, const int MaxIt = MAXIT, const RealType Epsilon = EPSILON ) {
    setCtrlCHandler ( );
    aol::MultiVector<RealType> approximateCenters;
    aol::MultiVector<int> approximateDimensions;
    if ( _progressBar != NULL ) _progressBar->setText ( "AtomFinder: segmenting and finding approximate centers (step 1/2)" );
    getApproximateSingleAtomPositions ( approximateCenters, approximateDimensions, Segmented, Data, Gamma, MaxIt, Epsilon );
    if ( _progressBar != NULL ) _progressBar->setText ( "AtomFinder: Gaussian fitting single atoms (step 2/2)" );
    getRefinedAtomPositions ( Centers, GaussianParams, approximateCenters, approximateDimensions, Data );
    
    if ( _diskOutput ) {
      PictureType bumpFunctionImg ( Data.getNumX ( ), Data.getNumY ( ) );
      getBumpFunctionImage ( GaussianParams, bumpFunctionImg );
      std::stringstream ss;
      ss << _outputDir << "/bumpFunctions" << getDefaultArraySuffix ( qc::QC_2D );
      bumpFunctionImg.save ( ss.str ( ).c_str ( ), qc::PGM_DOUBLE_BINARY );
      
      ColoredPictureType atomicCentersImg ( Data.getNumX ( ), Data.getNumY ( ) );
      for ( short c=0; c<3 ; ++c ) atomicCentersImg[c] = Data;
      atomicCentersImg.scaleValuesTo01 ( );
      aol::Vec2<short> pos;
      for ( int i=0; i<Centers.numComponents ( ) ; ++i ) {
        pos.set ( round ( Centers[i][0] ), round ( Centers[i][1] ) );
        atomicCentersImg[0].set ( pos, 1 );
        atomicCentersImg[1].set ( pos, 0 );
        atomicCentersImg[2].set ( pos, 0 );
      }
      ss.str ( std::string ( ) ); // clear stringstream
      ss << _outputDir << "/atomicCenters.png";
      atomicCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
      atomicCentersImg.savePNG ( ss.str ( ).c_str ( ) );
    }
    unsetCtrlCHandler ( );
  }
  
  void getDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams, qc::BitArray<qc::QC_2D> &Segmented,
                                  const PictureType &Data,
                                  const RealType Gamma = GAMMA, const int MaxIt = MAXIT, const RealType Epsilon = EPSILON ) {
    setCtrlCHandler ( );
    aol::MultiVector<RealType> approximateCenters;
    aol::MultiVector<int> approximateDimensions;
    aol::Vector<RealType> dumbbellOrientations, dumbbellSeparations;
    if ( _progressBar != NULL ) _progressBar->setText ( "AtomFinder: segmenting and finding approximate centers (step 1/2)" );
    getApproximateDumbbellAtomPositions ( approximateCenters, approximateDimensions, dumbbellOrientations, dumbbellSeparations, Segmented, Data, Gamma, MaxIt, Epsilon );
    if ( _progressBar != NULL ) _progressBar->setText ( "AtomFinder: Gaussian fitting dumbbells (step 2/2)" );
    getRefinedDumbbellAtomPositions ( Centers, GaussianParams, approximateCenters, approximateDimensions, dumbbellOrientations, dumbbellSeparations, Data );
    
    if ( _diskOutput ) {
      PictureType bumpFunctionImg ( Data.getNumX ( ), Data.getNumY ( ) );
      getBumpFunctionImage ( GaussianParams, bumpFunctionImg );
      std::stringstream ss;
      ss << _outputDir << "/bumpFunctions" << getDefaultArraySuffix ( qc::QC_2D );
      bumpFunctionImg.save ( ss.str ( ).c_str ( ), qc::PGM_DOUBLE_BINARY );
      
      ColoredPictureType atomicCentersImg ( Data.getNumX ( ), Data.getNumY ( ) );
      for ( short c=0; c<3 ; ++c ) atomicCentersImg[c] = Data;
      atomicCentersImg.scaleValuesTo01 ( );
      aol::Vec2<short> pos;
      for ( int i=0; i<Centers.numComponents ( ) ; ++i ) {
        pos.set ( round ( Centers[i][0] ), round ( Centers[i][1] ) );
        atomicCentersImg[0].set ( pos, 1 );
        atomicCentersImg[1].set ( pos, 0 );
        atomicCentersImg[2].set ( pos, 0 );
      }
      ss.str ( std::string ( ) ); // clear stringstream
      ss << _outputDir << "/atomicCenters.png";
      atomicCentersImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0, 1 );
      atomicCentersImg.savePNG ( ss.str ( ).c_str ( ) );
    }
    unsetCtrlCHandler ( );
  }
  
  void getAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams, qc::BitArray<qc::QC_2D> &Segmented,
                          const PictureType &Data, const aol::ParameterParser &Parser ) {
    aol::Vec2<short> atomSize;
    RealType atomSeparation, atomAngle, gamma, epsilon;
    int maxIt;
    bool singleFit, dumbbellFit;
    readParameters ( Parser, gamma, maxIt, epsilon, atomSize, atomSeparation, atomAngle, singleFit, dumbbellFit );
    if ( singleFit ) getAtomPositions ( Centers, GaussianParams, Segmented, Data, atomSize, gamma, maxIt, epsilon );
    else if ( dumbbellFit ) getDumbbellAtomPositions ( Centers, GaussianParams, Segmented, Data, atomSize, atomSeparation, atomAngle, gamma, maxIt, epsilon );
    else getAtomPositions ( Centers, GaussianParams, Segmented, Data, gamma, maxIt, epsilon );
  }
  
  void getAtomPositions ( aol::MultiVector<RealType> &Centers, const PictureType &Data, const aol::ParameterParser &Parser ) {
    aol::MultiVector<RealType> gaussianParams;
    qc::BitArray<qc::QC_2D> segmented;
    getAtomPositions ( Centers, gaussianParams, segmented, Data, Parser );
  }
  
  void getAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams, qc::BitArray<qc::QC_2D> &Segmented,
                          const PictureType &Data, const aol::Vec2<short> &AtomSize,
                          const RealType Gamma = GAMMA, const int MaxIt = MAXIT, const RealType Epsilon = EPSILON ) const {
    aol::MultiVector<RealType> approximateCenters;
    getApproximateAtomPositions ( approximateCenters, Segmented, Data, Gamma, MaxIt, Epsilon );
    getRefinedAtomPositions ( Centers, approximateCenters, GaussianParams, Data, AtomSize );
  }
  
  void getAtomPositions ( aol::MultiVector<RealType> &Centers,
                         const PictureType &Data, const aol::Vec2<short> &AtomSize,
                         const RealType Gamma = GAMMA, const int MaxIt = MAXIT, const RealType Epsilon = EPSILON ) const {
    aol::MultiVector<RealType> gaussianParams;
    qc::BitArray<qc::QC_2D> segmented;
    getAtomPositions ( Centers, gaussianParams, segmented, Data, AtomSize, Gamma, MaxIt, Epsilon );
  }
  
  void getDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams, qc::BitArray<qc::QC_2D> &Segmented,
                                 const PictureType &Data, const aol::Vec2<short> &AtomSize,
                                 const RealType AtomSeparation, const RealType AtomAngle,
                                 const RealType Gamma = GAMMA, const int MaxIt = MAXIT, const RealType Epsilon = EPSILON ) const {
    aol::MultiVector<RealType> approximateCenters;
    getApproximateAtomPositions ( approximateCenters, Segmented, Data, Gamma, MaxIt, Epsilon );
    getRefinedDumbbellAtomPositions ( Centers, approximateCenters, GaussianParams, Data, AtomSize, AtomSeparation, AtomAngle );
  }
  
  void getDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers,
                                 const PictureType &Data, const aol::Vec2<short> &AtomSize,
                                 const RealType AtomSeparation, const RealType AtomAngle,
                                 const RealType Gamma = GAMMA, const int MaxIt = MAXIT, const RealType Epsilon = EPSILON ) const {
    aol::MultiVector<RealType> gaussianParams;
    qc::BitArray<qc::QC_2D> segmented;
    getDumbbellAtomPositions ( Centers, gaussianParams, segmented, Data, AtomSize, AtomSeparation, AtomAngle, Gamma, MaxIt, Epsilon );
  }
  
  void getApproximateAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<int> &Dimensions,
                                     aol::MultiVector<RealType> &DumbbellCenters, aol::MultiVector<int> &DumbbellDimensions,
                                     aol::Vector<RealType> &DumbbellOrientations, aol::Vector<RealType> &DumbbellSeparations,
                                     qc::BitArray<qc::QC_2D> &Segmented,
                                     const PictureType &Data,
                                     const RealType Gamma = GAMMA, const int MaxIt = MAXIT, const RealType Epsilon = EPSILON ) const;
  
  void getApproximateSingleAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<int> &Dimensions,
                                           qc::BitArray<qc::QC_2D> &Segmented,
                                           const PictureType &Data,
                                           const RealType Gamma = GAMMA, const int MaxIt = MAXIT, const RealType Epsilon = EPSILON ) const;
  
  void getApproximateDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<int> &Dimensions,
                                             aol::Vector<RealType> &DumbbellOrientations, aol::Vector<RealType> &DumbbellSeparations,
                                             qc::BitArray<qc::QC_2D> &Segmented,
                                             const PictureType &Data,
                                             const RealType Gamma = GAMMA, const int MaxIt = MAXIT, const RealType Epsilon = EPSILON ) const;
  
  void getApproximateAtomPositions ( aol::MultiVector<RealType> &Centers, qc::BitArray<qc::QC_2D> &Segmented,
                                    const PictureType &Data,
                                    const RealType Gamma = GAMMA, const int MaxIt = MAXIT, const RealType Epsilon = EPSILON ) const;
  
  void getRefinedAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                                 const aol::MultiVector<RealType> &ApproximateCenters, const aol::MultiVector<int> &ApproximateDimensions,
                                 const PictureType &Data ) const;
  
  void getRefinedAtomPositions ( aol::MultiVector<RealType> &Centers, const aol::MultiVector<RealType> &ApproximateCenters, aol::MultiVector<RealType> &GaussianParams,
                                 const PictureType &Data, const aol::Vec2<short> &AtomSize ) const;
  
  void getRefinedDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers, aol::MultiVector<RealType> &GaussianParams,
                                         const aol::MultiVector<RealType> &ApproximateCenters, const aol::MultiVector<int> &ApproximateDimensions,
                                         const aol::Vector<RealType> &Orientations, const aol::Vector<RealType> &Separations,
                                         const PictureType &Data ) const;
  
  void getRefinedDumbbellAtomPositions ( aol::MultiVector<RealType> &Centers, const aol::MultiVector<RealType> &ApproximateCenters, aol::MultiVector<RealType> &GaussianParams,
                                         const PictureType &Data, const aol::Vec2<short> &AtomSize,
                                         const RealType AtomSeparation, const RealType AtomAngle ) const;

  void readCentersFromCSV ( aol::MultiVector<RealType> &Centers, const std::string &Path ) const;
  
  void setOutputDir ( const std::string &OutputDir ) {
    _outputDir = OutputDir;
    _diskOutput = OutputDir != "";
  }
  
  void setQuietMode ( const bool Quiet = true ) {
    _quietMode = Quiet;
  }
  
  static void getBumpFunctionImage ( const aol::MultiVector<RealType> &GaussianParameters, PictureType &BumpFunctionImg ) {
    BumpFunctionImg.setZero ( );
    for ( int k=0; k<GaussianParameters.numComponents ( ) ; ++k ) {
      if ( GaussianParameters[k].size ( ) == AsymmetricGaussianBumpFunction<RealType>::NumberOfParameters ) {
        const AsymmetricGaussianBumpFunction<RealType> bumpFunc ( GaussianParameters[k] );
        for ( int x=0; x<BumpFunctionImg.getNumX ( ) ; ++x )
          for ( int y=0; y<BumpFunctionImg.getNumY ( ) ; ++y )
            BumpFunctionImg.add ( x, y, bumpFunc.evaluate ( aol::Vec2<RealType> ( x, y ) ) );
      } else if ( GaussianParameters[k].size ( ) == AsymmetricGaussianDoubleBumpFunction<RealType>::NumberOfParameters ) {
        const AsymmetricGaussianDoubleBumpFunction<RealType> bumpFunc ( GaussianParameters[k] );
        for ( int x=0; x<BumpFunctionImg.getNumX ( ) ; ++x )
          for ( int y=0; y<BumpFunctionImg.getNumY ( ) ; ++y )
            BumpFunctionImg.add ( x, y, bumpFunc.evaluate ( aol::Vec2<RealType> ( x, y ) ) );
      }
    }
  }
  
  void setProgressBar ( aol::ProgressBar<> *ProgressBar ) {
    _progressBar = ProgressBar;
  }
  
  void setCatchCtrlC ( bool catchCtrlC ) {
    _catchCtrlC = catchCtrlC;
  }
protected:
  void getComponentsCollection ( ComponentsCollection<RealType> &ComponentsCollection,
                                 qc::BitArray<qc::QC_2D> &Segmented,
                                 const PictureType &Data,
                                 const RealType Gamma = GAMMA, const int MaxIt = MAXIT, const RealType Epsilon = EPSILON ) const;
  
  void readParameters ( const aol::ParameterParser &Parser, RealType &Gamma, int &MaxIt, RealType &Epsilon,
                        aol::Vec2<short> &AtomSize, RealType &AtomSeparation, RealType &AtomAngle,
                        bool &SingleFit, bool &DumbbellFit ) const {
    Gamma = Parser.getDoubleOrDefault ( "gamma", GAMMA );
    MaxIt = Parser.getIntOrDefault ( "maxIt", MAXIT );
    Epsilon = Parser.getDoubleOrDefault ( "epsilon", EPSILON );
    AtomSize.set ( Parser.getIntOrDefault ( "atomWidth", 0 ), Parser.getIntOrDefault ( "atomHeight", 0 ) );
    AtomSeparation = Parser.getDoubleOrDefault ( "atomSeparation", 0 );
    AtomAngle = Parser.getDoubleOrDefault ( "atomAngle", 0 );
    SingleFit = ( AtomSize[0] > 0 && AtomSize[1] > 0 && AtomSize[0] % 2 != 0 && AtomSize[1] % 2 != 0 );
    DumbbellFit = ( Parser.checkAndGetBool ( "dumbbellFit" ) && SingleFit && AtomSeparation > 0 && Parser.hasVariable ( "atomAngle" ) );
    if ( DumbbellFit ) SingleFit = false;
  }
  
  void setCtrlCHandler () const {
    if (_catchCtrlC)
      _previousCtrlCHandler = signal ( InterruptSignal, aol::ctrlCHandler );
  }
  
  void unsetCtrlCHandler () const {
    if (_catchCtrlC)
      signal ( InterruptSignal, _previousCtrlCHandler );
  }
  
  bool wantsInterrupt() const {
    if (!_catchCtrlC || !aol::getCtrlCState())
      return false;
    else
      return true;
  }
};


#endif /* ATOMFINDER_H_ */
