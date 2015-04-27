#ifndef NEIGHBORHOODFILTER_H_
#define NEIGHBORHOODFILTER_H_

// standard
#define _USE_MATH_DEFINES
#define _OPEN_SYS
#include <climits>
#include <cmath>
#include <dirent.h>
#include <iostream>
#include <vector>
#include <progressBar.h>

// quocmesh
#include <aol.h>
#include <ctrlCCatcher.h>
#include <multiArray.h>
#include <op.h>
#include <scalarArray.h>
#include "scalarArrayExtensions.h"
#include <qmException.h>
#include <bzipiostream.h>
#include "statistics.h"
#include "anscombe.h"
#include "similaritySearchIterators.h"
#include "referenceBlockIterators.h"


template <typename _RealType, typename _PictureType, typename _OptionsType>
class ReferenceBlockIterator;
template <typename _RealType, typename _PictureType, typename _OptionsType>
class GlobalReferenceBlockIterator;
template <typename _RealType, typename _PictureType, typename _OptionsType>
class SparseReferenceBlockIterator;

template <typename _RealType, typename _PictureType, typename _OptionsType>
class SimilaritySearchIterator;
template <typename _RealType, typename _PictureType, typename _OptionsType>
class LocalBSimilaritySearchIterator;
template <typename _RealType, typename _PictureType, typename _OptionsType>
class GlobalSimilaritySearchIterator;


class NOISE_TYPE {
public:
  static const int GAUSSIAN = 0;  // Remove additive Gaussian white noise (AGWN)
  static const int POISSON  = 1;  // Remove Poisson noise
  
  static const int NUM = 2;
  
  static std::string getIdentifier ( const int NoiseType ) {
    if ( NoiseType == GAUSSIAN ) return "G";
    else if ( NoiseType == POISSON ) return "P";
    else throw aol::Exception ( "Did not recognize noise type!", __FILE__, __LINE__ );
  }
  
  static std::string toString ( const int NoiseType ) {
    if ( NoiseType == GAUSSIAN ) return "Gaussian noise";
    else if ( NoiseType == POISSON ) return "Poisson noise";
    else throw aol::Exception ( "Did not recognize noise type!", __FILE__, __LINE__ );
  }
};

class POISSON_NOISE_ADAPTATION {
public:
  static const int ANSCOMBE           = 0;  // Three step procedure: 1. Anscombe transform, 2. Remove Gaussian noise, 3. Inverse Anscombe transform
  static const int MAXIMUMLIKELIHOOD  = 1;  // Directly apply denoising procedures on Poisson noise statistics using maximum-likelihood ratios as intensity distance measure
  
  static const int NUM = 2;
  
  static std::string getIdentifier ( const int PoissonNoiseAdaptation ) {
    if ( PoissonNoiseAdaptation == ANSCOMBE ) return "A";
    else if ( PoissonNoiseAdaptation == MAXIMUMLIKELIHOOD ) return "ML";
    else throw aol::Exception ( "Did not recognize Poisson noise adaptation!", __FILE__, __LINE__ );
  }
  
  static std::string toString ( const int PoissonNoiseAdaptation ) {
    if ( PoissonNoiseAdaptation == ANSCOMBE ) return "Anscombe transform";
    else if ( PoissonNoiseAdaptation == MAXIMUMLIKELIHOOD ) return "Maximum-likelihood ratios";
    else throw aol::Exception ( "Did not recognize Poisson noise adaptation!", __FILE__, __LINE__ );
  }
};

class SIMILARITYSEARCH_METHOD {
public:
  static const int LOCAL  = 0;  // A local square is used as a search window
  static const int GLOBAL = 1;  // Similarity search is performed on the whole image
  
  static const int NUM = 2;
  
  static std::string getIdentifier ( const int SimilaritySearchMethod ) {
    if ( SimilaritySearchMethod == LOCAL ) return "l";
    else if ( SimilaritySearchMethod == GLOBAL ) return "gl";
    else throw aol::Exception ( "NH similarity search method: Did not recognize similarity search method!", __FILE__, __LINE__ );
  }
  
  static std::string toString ( const int SimilaritySearchMethod ) {
    if ( SimilaritySearchMethod == LOCAL ) return "Local";
    else if ( SimilaritySearchMethod == GLOBAL ) return "Global";
    else throw aol::Exception ( "NH similarity search method: Did not recognize similarity search method!", __FILE__, __LINE__ );
  }
};

struct REFERENCEITERATION_METHOD {
  static const int GLOBAL = 0;  // The whole image is denoised (possibly except some boundary pixels)
  static const int PIXELS = 1;  // Only a specified set of pixels is denoised
  
  static const int NUM = 2;
};


template<typename _RealType, typename _PictureType>
struct NHFilterOptions {
  _PictureType input, preprocessedInput, *groundTruth;
  std::string outputDir;
  bool quietMode;
  aol::ProgressBar<> *progressBar;
  _RealType stdDev;
  short blockSize;
  int noiseType;
  int poissonNoiseAdaptation;
  int similaritySearchMethod;
  
  NHFilterOptions ( _PictureType &Input, const std::string OutputDir = "", const bool Quiet = true )
    : input ( Input ), preprocessedInput ( Input ), groundTruth ( NULL ),
      outputDir ( OutputDir ), quietMode ( Quiet ), progressBar ( NULL ),
      stdDev ( 0 ), blockSize ( 0 ), noiseType ( NOISE_TYPE::GAUSSIAN ), poissonNoiseAdaptation ( POISSON_NOISE_ADAPTATION::ANSCOMBE ),
      similaritySearchMethod ( SIMILARITYSEARCH_METHOD::LOCAL ) { }
  
  NHFilterOptions ( const NHFilterOptions<_RealType, _PictureType> &Options )
    : input ( Options.input ), preprocessedInput ( Options.preprocessedInput ), groundTruth ( Options.groundTruth ),
      outputDir ( Options.outputDir ), quietMode ( Options.quietMode ), progressBar ( Options.progressBar ),
      stdDev ( Options.stdDev ), blockSize ( Options.blockSize ), noiseType ( Options.noiseType ), poissonNoiseAdaptation ( Options.poissonNoiseAdaptation ),
      similaritySearchMethod ( Options.similaritySearchMethod ) { }
};


template <typename _RealType, typename _PictureType>
class NHFilterTrait {
public:
  typedef NHFilterOptions<_RealType, _PictureType> OptionsType;
  typedef ArbitraryAnchorBlockCollection<_RealType, _PictureType> BlockCollectionType;
  typedef SIMILARITYSEARCH_METHOD SimilaritySearchMethodType;
};


template<typename _RealType, typename _PictureType, typename NeighborhoodFilterTrait>
class NeighborhoodFilter {
  typedef _RealType RealType;
  typedef aol::FullMatrix<RealType> MatrixType;
  typedef _PictureType PictureType;
  typedef qc::MultiArray<RealType, qc::QC_2D, 3> ColoredPictureType;
public:
  typedef typename NeighborhoodFilterTrait::OptionsType OptionsType;
  typedef typename NeighborhoodFilterTrait::BlockCollectionType BlockCollectionType;
  typedef typename NeighborhoodFilterTrait::SimilaritySearchMethodType SimilaritySearchMethodType;
  typedef void ( NeighborhoodFilter<RealType, PictureType, NeighborhoodFilterTrait>::*BlockDistanceFunctionType )( RealType&, const aol::Vec2<short>&, const aol::Vec2<short>& ) const;
protected:
  // Input
  PictureType _input, _preprocessedInput;                                 // Noisy image
  PictureType _estimate;
  RealType _stdDev;                                                       // Standard deviation of the noise (if Gaussian) of the input image (for [0,255] scale values)
  
  // Output
  bool _quietMode;                                                        // In quiet mode no console output is generated
  std::string _outputDir;                                                 // If a path to an output directory is specified, detailed results of all steps are stored to disk
  aol::ProgressBar<> *_progressBar;
  bool _catchCtrlC;
  mutable sigfunc _previousCtrlCHandler;
  
  // Similarity search
  short _blockSize, _blockSizeSqr;                                        // Size of the square neighborhoods
  int _X0, _Y0, _NX, _NY, _XEnd, _YEnd, _effSize;                       // Parameters defining the dimensions of the subregion of the image containing well-defined block anchors
  short _searchWindowSize;                                                // Size of the window used for similarity search (in case of a local similarity search)
  aol::MultiVector<short> _pixels;                                        // Array of reference pixels for which similarity search is performed (in case of a sparse reconstruction)
  SimilaritySearchIterator<RealType, PictureType, NeighborhoodFilterTrait> *_searchIt;  // Pointer to iterator over pixels used within similarity search
  
  // Blocks
  BlockCollectionType *_blocks;         // Pointer to currently used block collection
  BlockCollectionType _blocksInitial;   // Block collection containing values from input image
  BlockDistanceFunctionType _setBlockDistance;                            // Pointer to currently used block distance function
  
  // Internal variables used for Anscombe transformation
  const RealType _scaleRange, _scaleShift, _minTransformed;
  RealType _maxTransformed;
public:
  NeighborhoodFilter ( ) : _catchCtrlC ( false ), _searchIt ( NULL ),
    _scaleRange ( 0.7 ),
    _scaleShift ( 0.5 * ( 1.0 - _scaleRange ) ),
    _minTransformed ( 2 * sqrt( 0.0 + 3.0 / 8.0 ) ) { }
  
  void setCatchCtrlC ( bool catchCtrlC ) {
    _catchCtrlC = catchCtrlC;
  }
  
  virtual ~NeighborhoodFilter ( ) { }
  
  virtual void apply ( const OptionsType &Options, PictureType &Dest ) = 0;
  
  virtual void setSimilaritySearchGrid ( const OptionsType &Options, ColoredPictureType &Dest, const aol::Vec2<short> &XRef ) = 0;
  
  static std::string getMethodIdentifier ( const OptionsType &Options ) {
    return getMethodPrefix ( Options ) + "-" + name ( );
  }
  
  static std::string getMethodPrefix ( const OptionsType &Options ) {
    std::stringstream ss;
    ss << NOISE_TYPE::getIdentifier ( Options.noiseType );
    if ( Options.noiseType == NOISE_TYPE::POISSON ) ss << "-" + POISSON_NOISE_ADAPTATION::getIdentifier ( Options.poissonNoiseAdaptation );
    ss << "-" << SimilaritySearchMethodType::getIdentifier ( Options.similaritySearchMethod );
    ss << "-" << Options.blockSize;
    return ss.str ( );
  }
  
  static const std::string name ( ) { return "NHF"; }
protected:
  void initialize ( const OptionsType &Options ) {
    _outputDir = Options.outputDir;
    _quietMode = Options.quietMode;
    _progressBar = Options.progressBar;
    _NX = Options.input.getNumX ( );
    _NY = Options.input.getNumY ( );
    _input.reallocate ( _NX, _NY );
    _input = Options.input;
    _preprocessedInput.reallocate ( _NX, _NY );
    _preprocessedInput = Options.preprocessedInput;
    _estimate.reallocate ( _NX, _NY );
    _estimate = _preprocessedInput;
    
    setSimilaritySearchIterator ( Options );
    setBlockDistanceFunction ( Options );
  }
  
  void preprocess ( OptionsType &Options ) {
    Options.preprocessedInput.reallocate ( Options.input.getNumX ( ), Options.input.getNumY ( ) );
    
    if ( Options.noiseType == NOISE_TYPE::GAUSSIAN ) {
      // Scale noisy image to [0,1]
      for ( int k=0; k<Options.preprocessedInput.size ( ) ; ++k )
        Options.preprocessedInput[k] = Options.input[k] / 255.0;
      
      // Set noise standard deviation
      if ( Options.stdDev > 0 ) _stdDev = Options.stdDev / 255.0;
      else throw aol::Exception ( "Noise standard deviation estimation is not supported yet!", __FILE__, __LINE__ );
    } else if ( Options.noiseType == NOISE_TYPE::POISSON && Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::ANSCOMBE ) {
      applyAnscombeForwardAndScale ( Options );
    } else
      Options.preprocessedInput = Options.input;
  }
  
  void postprocess ( const OptionsType &Options ) {
    if ( Options.noiseType == NOISE_TYPE::GAUSSIAN ) {
      // Scale estimate back to [0,255]
      for ( int k=0; k<_estimate.size ( ) ; ++k )
        _estimate[k] *= 255.0;
    } else if ( Options.noiseType == NOISE_TYPE::POISSON && Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::ANSCOMBE ) {
      scaleBackAndApplyAnscombeInverse ( Options );
    }
  }
  
  void applyAnscombeForwardAndScale ( OptionsType &Options ) {
    AnscombeForward<RealType> anscombeFW;
    anscombeFW.apply ( Options.input, Options.preprocessedInput );
    
    // Scale image to [0.15,0.85]
    _maxTransformed = Options.preprocessedInput.getMaxValue ( );
    for ( int k=0; k<Options.preprocessedInput.size ( ) ; ++k ) {
      Options.preprocessedInput[k] = ( Options.preprocessedInput[k] - _minTransformed ) / ( _maxTransformed - _minTransformed );
      Options.preprocessedInput[k] = Options.preprocessedInput[k] * _scaleRange + _scaleShift;
    }
    
    // Set noise standard deviation
    _stdDev = 1.0 / ( _maxTransformed - _minTransformed ) * _scaleRange;
  }
  
  void scaleBackAndApplyAnscombeInverse ( const OptionsType &/*Options*/ ) {
    for ( int k=0; k<_estimate.size ( ) ; ++k ) {
      _estimate[k] = ( _estimate[k] - _scaleShift ) / _scaleRange;
      _estimate[k] = _estimate[k] * ( _maxTransformed - _minTransformed ) + _minTransformed;
    }
    
    // Apply inverse Anscombe transformation
    AnscombeInverse<RealType> anscombeInverse;
    anscombeInverse.apply ( _estimate, _estimate );
  }
  
  virtual void setSimilaritySearchIterator ( const OptionsType &Options ) {
    if ( Options.similaritySearchMethod == SIMILARITYSEARCH_METHOD::GLOBAL )
      _searchIt = new GlobalSimilaritySearchIterator<RealType, PictureType, NeighborhoodFilterTrait> ( *this );
    else if ( Options.similaritySearchMethod == SIMILARITYSEARCH_METHOD::LOCAL )
      _searchIt = new LocalSimilaritySearchIterator<RealType, PictureType, NeighborhoodFilterTrait> ( *this );
    else
      throw aol::Exception ( "Neighborhood filter: Did not recognize similarity search method!", __FILE__, __LINE__ );
  }
  
  void setBlockDistanceFunction ( const OptionsType &Options ) {
    if ( Options.noiseType == NOISE_TYPE::GAUSSIAN || ( Options.noiseType == NOISE_TYPE::POISSON && Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::ANSCOMBE ) )
      _setBlockDistance = &NeighborhoodFilter<RealType, PictureType, NeighborhoodFilterTrait>::setBlockDistanceL2NormSqr;
    else if ( Options.noiseType == NOISE_TYPE::POISSON )
      _setBlockDistance = &NeighborhoodFilter<RealType, PictureType, NeighborhoodFilterTrait>::setBlockDistancePoissonLikelihoodRatio;
    else
      throw aol::Exception ( "Could not determine appropriate block distance measure!", __FILE__, __LINE__ );
  }
  
  virtual inline void setBlockDistanceL2NormSqr ( RealType &Dist, const aol::Vec2<short> &XRef, const aol::Vec2<short> &X ) const {
    Dist = 0.0;
    for ( short k=0; k<_blockSizeSqr; ++k )
      Dist += aol::Sqr<RealType> ( (*this->_blocks).get ( XRef, k ) - (*this->_blocks).get ( X, k ) );
  }
  
  virtual inline void setBlockDistancePoissonLikelihoodRatio ( RealType &Dist, const aol::Vec2<short> &XRef, const aol::Vec2<short> &X ) const {
    Dist = 0.0;
    for ( short k=0; k<_blockSizeSqr; ++k ) {
      const int k1 = (*_blocks).get ( XRef, k ), k2 = (*_blocks).get ( X, k );
      Dist += ( ( k1 > 0 ) ? k1 * log ( k1 ) : 0 ) + ( ( k2 > 0 ) ? k2 * log ( k2 ) : 0 ) - ( ( k1 + k2 > 0 ) ? ( k1 + k2 ) * log ( 0.5 * ( k1 + k2 ) ) : 0 );
    }
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
  
  
  friend class ReferenceBlockIterator<RealType, PictureType, NeighborhoodFilterTrait>;
  friend class GlobalReferenceBlockIterator<RealType, PictureType, NeighborhoodFilterTrait>;
  friend class PixelReferenceBlockIterator<RealType, PictureType, NeighborhoodFilterTrait>;
  
  friend class SimilaritySearchIterator<RealType, PictureType, NeighborhoodFilterTrait>;
  friend class LocalSimilaritySearchIterator<RealType, PictureType, NeighborhoodFilterTrait>;
  friend class GlobalSimilaritySearchIterator<RealType, PictureType, NeighborhoodFilterTrait>;
};


template<typename _RealType, typename _PictureType, typename NeighborhoodFilterTrait>
class CollaborativeNeighborhoodFilter : public NeighborhoodFilter<_RealType, _PictureType, NeighborhoodFilterTrait> {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  short _refStep;
  
  
  friend class ReferenceBlockIterator<RealType, PictureType, NeighborhoodFilterTrait>;
  friend class GlobalReferenceBlockIterator<RealType, PictureType, NeighborhoodFilterTrait>;
  friend class PixelReferenceBlockIterator<RealType, PictureType, NeighborhoodFilterTrait>;
  
  friend class SimilaritySearchIterator<RealType, PictureType, NeighborhoodFilterTrait>;
  friend class LocalSimilaritySearchIterator<RealType, PictureType, NeighborhoodFilterTrait>;
  friend class GlobalSimilaritySearchIterator<RealType, PictureType, NeighborhoodFilterTrait>;
};

#endif /* NEIGHBORHOODFILTER_H_ */