#ifndef BM3DFILTER_H_
#define BM3DFILTER_H_

#include "neighborhoodFilter.h"
#include "scalarArrayExtensions.h"
#include "blockOperators3D.h"


class BM3D_SIMILARITYSEARCH_METHOD : public SIMILARITYSEARCH_METHOD {
public:
  static std::string getIdentifier ( const int SimilaritySearchMethod ) {
    return SIMILARITYSEARCH_METHOD::getIdentifier ( SimilaritySearchMethod );
  }
  
  static std::string toString ( const int SimilaritySearchMethod ) {
    return SIMILARITYSEARCH_METHOD::toString ( SimilaritySearchMethod );
  }
};

class BM3D_PROFILE {
public:
  static const int NP = 0;      // Normal Profile (default, balanced quality)
  static const int LC = 1;      // Low Complexity Profile (fast, lower quality)
  static const int HIGH = 2;    // High Profile (high quality, not documented in [1])
  
  static const int NUM = 3;
  
  static std::string getIdentifier ( const int Profile ) {
    if ( Profile == NP ) return "np";
    else if ( Profile == LC ) return "lc";
    else if ( Profile == HIGH ) return "high";
    else throw aol::Exception ( "Did not recognize profile!", __FILE__, __LINE__ );
  }
  
  static std::string toString ( const int Profile ) {
    if ( Profile == NP ) return "Normal";
    else if ( Profile == LC ) return "Low complexity";
    else if ( Profile == HIGH ) return "High quality";
    else throw aol::Exception ( "Did not recognize profile!", __FILE__, __LINE__ );
  }
};

class BM3D_ESTIMATION_METHOD {
public:
  static const int HT = 0;      // Parameters are set for hard-thresholding (HT) in 2D and 3D transform domain
  static const int WIENER = 1;  // Parameters are set for Wiener filtering in 3D transform domain
  
  static const int NUM = 2;
  
  static std::string getIdentifier ( const int EstimationMethod ) {
    if ( EstimationMethod == HT ) return "ht";
    if ( EstimationMethod == WIENER ) return "wie";
    else throw aol::Exception ( "Did not recognize estimation method!", __FILE__, __LINE__ );
  }
  
  static std::string toString ( const int EstimationMethod ) {
    if ( EstimationMethod == HT ) return "Hard-thresholding";
    if ( EstimationMethod == WIENER ) return "Wiener filtering";
    else throw aol::Exception ( "Did not recognize estimation method!", __FILE__, __LINE__ );
  }
};

template<typename _RealType, typename _PictureType>
struct BM3DOptions : public NHFilterOptions<_RealType, _PictureType> {
  int profile;
  int method;
  
  BM3DOptions ( _PictureType &Input, const std::string OutputDir = "", const bool Quiet = true )
    : NHFilterOptions<_RealType, _PictureType> ( Input, OutputDir, Quiet ), profile ( BM3D_PROFILE::NP ), method ( BM3D_ESTIMATION_METHOD::HT ) {
    this->blockSize = 8;
  }
  
  BM3DOptions ( const BM3DOptions<_RealType, _PictureType> &Options )
    : NHFilterOptions<_RealType, _PictureType> ( Options ), profile ( Options.profile ), method ( Options.method ) { }
};


template <typename _RealType, typename _PictureType>
class BM3DFilterTrait {
public:
  typedef BM3DOptions<_RealType, _PictureType> OptionsType;
  typedef BlockCollection<_RealType, _PictureType> BlockCollectionType;
  typedef BM3D_SIMILARITYSEARCH_METHOD SimilaritySearchMethodType;
};


//  The following class BM3DFilter has been implemented by
//
//    Niklas Mevenkamp, email: mevenkamp@aices.rwth-aachen.de
//
//  based on the description of the BM3D algorithm as published in:
//
//      K. Dabov, A. Foi, V. Katkovnik, and K. Egiazarian, "Image Denoising
//      by Sparse 3D Transform-Domain Collaborative Filtering,"
//      IEEE Transactions on Image Processing, vol. 16, no. 8, August, 2007.
//      preprint at http://www.cs.tut.fi/~foi/GCF-BM3D
//
//  AUTHORS:
//      Kostadin Dabov, email: dabov _at_ cs.tut.fi
template<typename _RealType, typename _PictureType, typename FilterTrait = BM3DFilterTrait<_RealType, _PictureType> >
class BM3DFilter : public CollaborativeNeighborhoodFilter<_RealType, _PictureType, FilterTrait> {
  typedef _RealType RealType;
  typedef aol::Vector<RealType> VectorType;
  typedef aol::FullMatrix<RealType> MatrixType;
  typedef _PictureType PictureType;
  typedef qc::MultiArray<RealType, qc::QC_2D, 3> ColoredPictureType;
  typedef qc::ScalarArray<RealType, qc::QC_2D> BlockType;
  typedef qc::ScalarArray<RealType, qc::QC_3D> BlockStackType;
  typedef typename FilterTrait::OptionsType OptionsType;
  typedef typename FilterTrait::BlockCollectionType BlockCollectionType;
  typedef typename FilterTrait::SimilaritySearchMethodType SimilaritySearchMethodType;
protected:
  short _N2;                                           // Maximum number of similar blocks for 3D stacks (must be a power of 2, if wavelet denoising is used along 3rd dimension)
  RealType _tauMatch;                                  // Threshold for the block distance (d-distance)
  
  // Reconstruction
  RealType _beta;                                      // Parameter for Bessel functions used to create Kaiser window
  BlockType _WWin2D;                                   // 2D Kaiser window used in the aggregation of the HT and Wiener filtering weights
  
  
  // Internal variables
  qc::FastILexMapper<qc::QC_2D> _mapper;               // Class that transforms between flat (linear) and global (2D) indices
  PictureType _eBuff, _wBuff;                          // Buffers to accumulate estimates and their weights
  PictureType _numMatchedBlocks, _numAggregates, _numAggregatesHT, _numAggregatesWiener, _distances, _cumulatedDistancesHT, _cumulatedDistancesWiener;           // FOR TESTING ONLY
  BlockCollectionType _blocksEstimate, _blocksInitialTransformed, _blocksEstimateTransformed;
  
  // 2D Operators
  LinearUnitary2DTransformOp<RealType> _unitaryTransform2D;
  
  // 3D Operators
  BlockStackDenoisingOp<RealType> *_op3D;              // Operator applied to 3D block stacks (usually: forward 3D transform -> denoising -> inverse 3D transform)
  HTDenoisingInLocal2D1DTransformDomainOp<RealType> _htDenoisingInLocal2D1DTransformDomainOp;
  WienerDenoisingInLocal2D1DTransformDomainOp<RealType> _wienerDenoisingInLocal2D1DTransformDomainOp;
public:
  BM3DFilter ( );
  
  void apply ( const OptionsType &Options, PictureType &Dest );
  
  void setSimilaritySearchGrid ( const OptionsType &Options, ColoredPictureType &Dest, const aol::Vec2<short> &XRef );
  
  static std::string getMethodIdentifier ( const OptionsType &Options ) {
    return getMethodPrefix ( Options ) + "-" + name ( );
  }
  
  static std::string getMethodPrefix ( const OptionsType &Options );
  
  static std::string name ( ) { return "BM3D"; }
  
  PictureType& getNumAggregatesHT ( ) {
    return _numAggregatesHT;
  }
  
  PictureType& getNumAggregatesWiener ( ) {
    return _numAggregatesWiener;
  }
  
  PictureType& getCumulatedDistancesHT ( ) {
    return _cumulatedDistancesHT;
  }
  
  PictureType& getCumulatedDistancesWiener ( ) {
    return _cumulatedDistancesWiener;
  }

  
protected:
  void preprocess ( OptionsType &Options ) {
    CollaborativeNeighborhoodFilter<RealType, PictureType, FilterTrait>::preprocess ( Options );
    
    // For BM3D with Poisson Maximum-Likelihood ratios, the Anscombe transformed image is still required to fill the 3D block stacks
    if ( Options.noiseType == NOISE_TYPE::POISSON && Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::MAXIMUMLIKELIHOOD ) {
      this->applyAnscombeForwardAndScale ( Options );
    }
  }
  
  virtual void initializeOnce ( const OptionsType &Options );
  
  void initializeIteration ( const OptionsType &Options );
  
  virtual void fillBlocksFromPreviousEstimate ( const OptionsType &Options );
  
  void applyLinearUnitary2DTransformsToBlocks ( const OptionsType &Options );

  void denoise ( const OptionsType &Options );
  
  virtual void blockMatching ( const OptionsType &Options, const aol::Vec2<short> &XRef, aol::MultiVector<short> &Sx, VectorType &BlockDistances, int &Nz );
  
  void blockStacking ( const OptionsType &Options, const aol::MultiVector<short> &Sx, const short Nz,
                       Blocks3DInitialAndEstimate<RealType> &Blocks3DInitialAndEstimate );
  
  virtual void aggregate ( const OptionsType &/*Options*/, const Blocks3DAndWeight<RealType> &Blocks3D, const aol::MultiVector<short> &Sx, const short Nz );
  
  virtual void setEstimateFromBuffers ( const OptionsType &/*Options*/ );
  
  void setParameters ( const OptionsType &Options );
  
  void setOperators ( const OptionsType &Options );
  
  int numAggregatesSum ( const short X, const short Y ) {
    int res = 0;
    for ( short xBlock=0; xBlock<this->_blockSize ; ++xBlock )
      for ( short yBlock=0; yBlock<this->_blockSize ; ++yBlock )
        res += _numAggregates.get ( X + xBlock, Y + yBlock );
    return res;
  }
    

  friend class ReferenceBlockIterator<RealType, PictureType, OptionsType>;
  friend class GlobalReferenceBlockIterator<RealType, PictureType, OptionsType>;
  
  friend class SimilaritySearchIterator<RealType, PictureType, OptionsType>;
  friend class LocalSimilaritySearchIterator<RealType, PictureType, OptionsType>;
  friend class GlobalSimilaritySearchIterator<RealType, PictureType, OptionsType>;
};

#endif /* BM3DFILTER_H_ */
