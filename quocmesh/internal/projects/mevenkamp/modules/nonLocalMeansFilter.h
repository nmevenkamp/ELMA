#ifndef NONLOCALMEANSFILTER_H_
#define NONLOCALMEANSFILTER_H_


#include "neighborhoodFilter.h"
#include "scalarArrayExtensions.h"


class NLM_SIMILARITYSEARCH_METHOD : public SIMILARITYSEARCH_METHOD {
public:
  static std::string getIdentifier ( const int SimilaritySearchMethod ) {
    return SIMILARITYSEARCH_METHOD::getIdentifier ( SimilaritySearchMethod );
  }
  
  static std::string toString ( const int SimilaritySearchMethod ) {
    return SIMILARITYSEARCH_METHOD::toString ( SimilaritySearchMethod );
  }
};


template <typename _RealType, typename _PictureType>
struct NLMOptions : NHFilterOptions<_RealType, _PictureType> {
  short searchWindowSize;
  _RealType filterParameter;
  aol::MultiVector<short> pixels;
  
  NLMOptions ( _PictureType &Input, const std::string OutputDir = "", const bool Quiet = true )
    : NHFilterOptions<_RealType, _PictureType> ( Input, OutputDir, Quiet ), searchWindowSize ( 0 ), filterParameter ( 0 ) { }
};


template <typename _RealType, typename _PictureType>
class NonLocalMeansFilterTrait {
public:
  typedef NLMOptions<_RealType, _PictureType> OptionsType;
  typedef CenteredBlockCollection<_RealType, _PictureType> BlockCollectionType;
  typedef NLM_SIMILARITYSEARCH_METHOD SimilaritySearchMethodType;
};


//  The following class NonLocalMeansFilter has been implemented by
//
//    Niklas Mevenkamp, email: mevenkamp@aices.rwth-aachen.de
//
//  based on the description of the Non-local Means filter as published in:
//
//      A. Buades, B. Coll, J.M. Morel "A review of image denoising methods, with a new one"
//      Multiscale Modeling and Simulation, Vol. 4 (2), pp: 490-530, 2006. DOI: 10.1137/040616024
//      available at http://epubs.siam.org/doi/abs/10.1137/040616024
//
//  Online applet and source code of patchwise implemenation: http://www.ipol.im/pub/art/2011/bcm_nlm/
//
//  AUTHORS:
//      Antoni Buades toni.buades@uib.es, CNRS-Paris Descartes
//      Bartomeu Coll tomeu.coll@uib.es, Universitat Illes Balears
//      Jean-Michel Morel morel@cmla.ens-cachan.fr, CMLA, ENS-Cachan
template<typename _RealType, typename _PictureType, typename FilterTrait = NonLocalMeansFilterTrait<_RealType, _PictureType> >
class NonLocalMeansFilter : public NeighborhoodFilter<_RealType, _PictureType, FilterTrait> {
  typedef _RealType RealType;
  typedef aol::FullMatrix<RealType> MatrixType;
  typedef _PictureType PictureType;
  typedef qc::MultiArray<RealType, qc::QC_2D, 3> ColoredPictureType;
  typedef typename FilterTrait::OptionsType OptionsType;
  typedef typename FilterTrait::BlockCollectionType BlockCollectionType;
protected:
  short _blockOffset;
  
  RealType _filterParameter, _filterParameterSqr, _stdDevSqr;
public:
  void apply ( const OptionsType &Options, PictureType &Dest );
  
  void setSimilaritySearchGrid ( const OptionsType &Options, ColoredPictureType &Dest, const aol::Vec2<short> &XRef ) { }
  
  void setWeights ( const OptionsType &Options, const aol::Vec2<short> &XRef, PictureType &Weights, const bool ReInitialize = true );
  
  void cropBorders ( PictureType &Picture ) {
    Picture.crop ( aol::Vec2<int> ( _blockOffset ), aol::Vec2<int> ( Picture.getNumXYZ ( ) - 2 * _blockOffset, Picture.getNumY ( ) - 2 * _blockOffset ) );
  }
  
  static std::string getMethodIdentifier ( const OptionsType &Options ) {
    return getMethodPrefix ( Options ) + "-" + name ( );
  }
  
  static std::string getMethodPrefix ( const OptionsType &Options );
  
  static std::string name ( ) { return "NLM"; }
protected:
  void denoise ( const OptionsType &Options );
  
  void initialize ( const OptionsType &Options );
  
  void initializeBlocks ( const OptionsType &/*Options*/ );
  
  void setParameters ( const OptionsType &Options );
  
  void setOperators ( const OptionsType &Options );
  
  inline RealType getWeight ( const RealType Dist ) {
    return exp ( -aol::Max<RealType> ( Dist / this->_blockSizeSqr - 2 * _stdDevSqr, 0.0 ) / _filterParameterSqr  );
  }
  
  RealType getPUniformityOptimalFilterParameter ( const OptionsType &Options );

  
  friend class SimilaritySearchIterator<RealType, PictureType, OptionsType>;
  friend class LocalSimilaritySearchIterator<RealType, PictureType, OptionsType>;
  friend class GlobalSimilaritySearchIterator<RealType, PictureType, OptionsType>;
};


#endif /* NONLOCALMEANSFILTER_H_ */