#ifndef EMNEIGHBORHOODFILTER_H_
#define EMNEIGHBORHOODFILTER_H_

#include "neighborhoodFilter.h"
#include "blockRegularization.h"
#include "globalScanDistortionCorrection.h"


template <typename _RealType, typename _PictureType>
class PeriodicBlockMatchingIterator;


class EM_SIMILARITYSEARCH_METHOD : public SIMILARITYSEARCH_METHOD {
public:
  static const int PERIODIC = 2; // Periodic search grid for a perfect crystal lattice is used
  
  static const int NUM = 3;
  
  static std::string getIdentifier ( const int SimilaritySearchMethod ) {
    if ( SimilaritySearchMethod < PERIODIC ) return SIMILARITYSEARCH_METHOD::getIdentifier ( SimilaritySearchMethod );
    else if ( SimilaritySearchMethod == PERIODIC ) return "pi";
    else throw aol::Exception ( "EM similarity search method: Did not recognize similarity search method!", __FILE__, __LINE__ );
  }
  
  static std::string toString ( const int SimilaritySearchMethod ) {
    if ( SimilaritySearchMethod < PERIODIC ) return SIMILARITYSEARCH_METHOD::toString ( SimilaritySearchMethod );
    else if ( SimilaritySearchMethod == PERIODIC ) return "Periodic";
    else throw aol::Exception ( "EM similarity search method: Did not recognize similarity search method!", __FILE__, __LINE__ );
  }
};


class EM_SCANDISTORTIONCORRECTION_METHOD {
public:
  static const int NONE = 0;  // No correction of scan distortion is performed
  
  static const int BR_NORMAL                = 1; // Patches are regularized by minimizing the vertical pixel gradient over a vector of shifts per patch
  static const int BR_REGSHIFTS             = 2; // Same as normal, but additionally the shifts are post-processed by horizontal regularization
  static const int BR_REGSHIFTS_PREREGROWS  = 3; // Same as REGSHIFTS, but additionally the rows are regularized before the block regularization
  
  static const int JITTERBUG = 4;  // Localized registration of scan line segments via normalized cross-correlation (Jitterbug)
  
  
  static const int NUM = 5;
  
  static const int SHIFTMODE_NONE = 0;
  static const int SHIFTMODE_BLOCKROWSHIFTS = 1;
  static const int SHIFTMODE_GLOBALHORIZONTALSHIFTS = 2;
  
  static std::string getIdentifier ( const int ScanDistortionCorrectionMethod ) {
    if ( ScanDistortionCorrectionMethod == NONE ) return "NR";
    else if ( ScanDistortionCorrectionMethod == BR_NORMAL ) return "R";
    else if ( ScanDistortionCorrectionMethod == BR_REGSHIFTS ) return "RR";
    else if ( ScanDistortionCorrectionMethod == BR_REGSHIFTS_PREREGROWS ) return "RRR";
    else if ( ScanDistortionCorrectionMethod == JITTERBUG ) return "JB";
    else throw aol::Exception ( "Did not recognize scan distortion correction method!", __FILE__, __LINE__ );
  }
  
  static std::string toString ( const int ScanDistortionCorrectionMethod ) {
    if ( ScanDistortionCorrectionMethod == NONE ) return "None";
    else if ( ScanDistortionCorrectionMethod == BR_NORMAL ) return "Block regularization";
    else if ( ScanDistortionCorrectionMethod == BR_REGSHIFTS ) return "Block regularization + shift reg.";
    else if ( ScanDistortionCorrectionMethod == BR_REGSHIFTS_PREREGROWS ) return "Block regularization + shift reg. + row reg.";
    else if ( ScanDistortionCorrectionMethod == JITTERBUG ) return "Jitterbug";
    else throw aol::Exception ( "Did not recognize scan distortion correction method!", __FILE__, __LINE__ );
  }
  
  static int getShiftMode ( const int ScanDistortionCorrectionMethod ) {
    if ( ScanDistortionCorrectionMethod >= BR_NORMAL && ScanDistortionCorrectionMethod <= BR_REGSHIFTS_PREREGROWS ) return SHIFTMODE_BLOCKROWSHIFTS;
    else if ( ScanDistortionCorrectionMethod == JITTERBUG ) return SHIFTMODE_GLOBALHORIZONTALSHIFTS;
    else return 0;
  }
};


template<typename _RealType, typename _PictureType, typename BaseClass>
struct EMNHFilterOptions : BaseClass {
  const int noiseType = NOISE_TYPE::POISSON;
  int scanDistortionCorrectionMethod;
  
  aol::MultiVector<_RealType> *blockRowShifts;
  _PictureType *globalHorizontalShifts;
  bool useGroundTruthForScanDistortionCorrection;
  
  EMNHFilterOptions ( _PictureType &Input, const std::string OutputDir = "", const bool Quiet = true )
    : BaseClass ( Input, OutputDir, Quiet ), scanDistortionCorrectionMethod ( EM_SCANDISTORTIONCORRECTION_METHOD::NONE ),
      blockRowShifts ( NULL ), globalHorizontalShifts ( NULL ), useGroundTruthForScanDistortionCorrection ( false ) { }
  
  EMNHFilterOptions ( const EMNHFilterOptions<_RealType, _PictureType, BaseClass> &Options )
    : BaseClass ( Options ), scanDistortionCorrectionMethod ( Options.scanDistortionCorrectionMethod ),
      blockRowShifts ( Options.blockRowShifts ), globalHorizontalShifts ( Options.globalHorizontalShifts ), useGroundTruthForScanDistortionCorrection ( Options.useGroundTruthForScanDistortionCorrection ) { }
};


template<typename _RealType, typename _PictureType, typename NeighborhoodFilterTrait, typename BaseClass>
class EMNeighborhoodFilter : public BaseClass {
  typedef _RealType RealType;
  typedef aol::FullMatrix<RealType> MatrixType;
  typedef _PictureType PictureType;
  typedef LinearRegressionQR<RealType> LinearRegressionType;
  typedef typename NeighborhoodFilterTrait::OptionsType OptionsType;
  typedef typename NeighborhoodFilterTrait::BlockCollectionType BlockCollectionType;
  typedef typename NeighborhoodFilterTrait::SimilaritySearchMethodType SimilaritySearchMethodType;
public:
  static std::string getMethodIdentifier ( const OptionsType &Options ) {
    return getMethodPrefix ( Options ) + "-" + name ( );
  }
  
  static std::string getMethodPrefix ( const OptionsType &Options ) {
    std::string methodPrefix = BaseClass::getMethodPrefix ( Options );
    methodPrefix += "-" + EM_SCANDISTORTIONCORRECTION_METHOD::getIdentifier ( Options.scanDistortionCorrectionMethod );
    return methodPrefix;
  }
  
  static std::string name ( ) { return "EMNHF"; }
protected:
  aol::MultiVector<RealType> _blockRowShifts; // Shift vectors applied to each block (used )
  PictureType _globalHorizontalShifts;        // A single horizontal shift is applied to each pixel
  
  using BaseClass::setSimilaritySearchIterator;
  virtual void setSimilaritySearchIterator ( const OptionsType &Options ) {
    if ( Options.similaritySearchMethod < EM_SIMILARITYSEARCH_METHOD::PERIODIC )
      BaseClass::setSimilaritySearchIterator ( Options );
    else {
      if ( Options.similaritySearchMethod == EM_SIMILARITYSEARCH_METHOD::PERIODIC )
        this->_searchIt = new PeriodicSimilaritySearchIterator<RealType, PictureType, NeighborhoodFilterTrait, BaseClass> ( *this );
      else
        throw aol::Exception ( "EMNHFilter: Did not recognize similarity search method!", __FILE__, __LINE__ );
    }
  }
  
  void setScanDistortionCorrectionParameters ( const OptionsType &Options ) {
    if ( Options.scanDistortionCorrectionMethod != EM_SCANDISTORTIONCORRECTION_METHOD::NONE ) {
      if ( Options.blockRowShifts != NULL && Options.blockRowShifts->numComponents ( ) == Options.input.size ( ) ) {
        _blockRowShifts.reallocate ( Options.blockRowShifts->numComponents ( ), this->_blockSize );
        for ( int k=0; k<_blockRowShifts.numComponents ( ) ; ++k )
          _blockRowShifts[k] = (*Options.blockRowShifts)[k];
      } else if ( Options.globalHorizontalShifts != NULL && Options.globalHorizontalShifts->getNumX ( ) == this->_NX && Options.globalHorizontalShifts->getNumY ( ) == this->_NY ) {
        _globalHorizontalShifts.reallocate ( this->_NX, this->_NY );
        _globalHorizontalShifts = *Options.globalHorizontalShifts;
      } else {
        PictureType blockRegularizerInput ( Options.preprocessedInput );
        if ( Options.useGroundTruthForScanDistortionCorrection ) blockRegularizerInput = *Options.groundTruth;
        
        if ( EM_SCANDISTORTIONCORRECTION_METHOD::getShiftMode ( Options.scanDistortionCorrectionMethod ) == EM_SCANDISTORTIONCORRECTION_METHOD::SHIFTMODE_BLOCKROWSHIFTS ) {
          this->_blocksInitial.fillFrom ( Options.preprocessedInput, this->_blockSize );
          MultiBlockRegularizer<RealType, MatrixType, PictureType, LinearRegressionType> multiBlockRegularizer ( blockRegularizerInput, 1, 10, 1, 0.2, 0.8, 1e-6, 1e-6, 1e-6, true );
          multiBlockRegularizer.apply ( this->_blocksInitial, _blockRowShifts,
                                        Options.scanDistortionCorrectionMethod == EM_SCANDISTORTIONCORRECTION_METHOD::BR_REGSHIFTS
                                     || Options.scanDistortionCorrectionMethod == EM_SCANDISTORTIONCORRECTION_METHOD::BR_REGSHIFTS_PREREGROWS,
                                        Options.scanDistortionCorrectionMethod == EM_SCANDISTORTIONCORRECTION_METHOD::BR_REGSHIFTS_PREREGROWS,
                                        true );
        } else if ( EM_SCANDISTORTIONCORRECTION_METHOD::getShiftMode ( Options.scanDistortionCorrectionMethod ) == EM_SCANDISTORTIONCORRECTION_METHOD::SHIFTMODE_GLOBALHORIZONTALSHIFTS ) {
          if ( Options.scanDistortionCorrectionMethod == EM_SCANDISTORTIONCORRECTION_METHOD::JITTERBUG )
            jitterbug<RealType, PictureType> ( blockRegularizerInput, _globalHorizontalShifts );
        }
      }
    }
  }
  
  void fillBlocks ( BlockCollectionType &Blocks, const PictureType &Input, const short BlockSize, const OptionsType &Options ) {
    if ( EM_SCANDISTORTIONCORRECTION_METHOD::getShiftMode ( Options.scanDistortionCorrectionMethod ) == EM_SCANDISTORTIONCORRECTION_METHOD::SHIFTMODE_NONE )
      Blocks.fillFrom ( Input, BlockSize );
    else if ( EM_SCANDISTORTIONCORRECTION_METHOD::getShiftMode ( Options.scanDistortionCorrectionMethod ) == EM_SCANDISTORTIONCORRECTION_METHOD::SHIFTMODE_BLOCKROWSHIFTS )
      Blocks.fillFrom ( Input, BlockSize, _blockRowShifts );
    else if ( EM_SCANDISTORTIONCORRECTION_METHOD::getShiftMode ( Options.scanDistortionCorrectionMethod ) == EM_SCANDISTORTIONCORRECTION_METHOD::SHIFTMODE_GLOBALHORIZONTALSHIFTS )
      Blocks.fillFrom ( Input, BlockSize, _globalHorizontalShifts );
    else
      throw aol::Exception ( "Did not recognize specified scan distortion correction method!", __FILE__, __LINE__ );
  }
  
  
  friend class ReferenceBlockIterator<RealType, PictureType, NeighborhoodFilterTrait>;
  friend class GlobalReferenceBlockIterator<RealType, PictureType, NeighborhoodFilterTrait>;
  friend class PixelReferenceBlockIterator<RealType, PictureType, NeighborhoodFilterTrait>;
  
  friend class SimilaritySearchIterator<RealType, PictureType, NeighborhoodFilterTrait>;
  friend class LocalSimilaritySearchIterator<RealType, PictureType, NeighborhoodFilterTrait>;
  friend class GlobalSimilaritySearchIterator<RealType, PictureType, NeighborhoodFilterTrait>;
  friend class PeriodicSimilaritySearchIterator<RealType, PictureType, NeighborhoodFilterTrait, BaseClass>;
};

#endif /* EMNEIGHBORHOODFILTER_H_ */