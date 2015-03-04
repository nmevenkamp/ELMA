#ifndef EMNONLOCALMEANSFILTER_H_
#define EMNONLOCALMEANSFILTER_H_


#include "emNeighborhoodFilter.h"
#include "nonLocalMeansFilter.h"


class EMNLM_SIMILARITYSEARCH_METHOD : public EM_SIMILARITYSEARCH_METHOD {
public:
  static std::string getIdentifier ( const int SimilaritySearchMethod ) {
    return EM_SIMILARITYSEARCH_METHOD::getIdentifier ( SimilaritySearchMethod );
  }
  
  static std::string toString ( const int SimilaritySearchMethod ) {
    return EM_SIMILARITYSEARCH_METHOD::toString ( SimilaritySearchMethod );
  }
};

class EMNLM_SCANDISTORTIONCORRECTION_METHOD : public EM_SCANDISTORTIONCORRECTION_METHOD { };


template <typename _RealType, typename _PictureType>
struct EMNLMOptions : public EMNHFilterOptions<_RealType, _PictureType, NLMOptions<_RealType, _PictureType> > {
  EMNLMOptions ( _PictureType &Input, const std::string OutputDir = "", const bool Quiet = true )
    : EMNHFilterOptions<_RealType, _PictureType, NLMOptions<_RealType, _PictureType> > ( Input, OutputDir, Quiet ) { }
};


template <typename _RealType, typename _PictureType>
class EMNonLocalMeansFilterTrait : public NHFilterTrait<_RealType, _PictureType> {
public:
  typedef EMNLMOptions<_RealType, _PictureType> OptionsType;
  typedef CenteredShiftedBlockCollection<_RealType, _PictureType> BlockCollectionType;
  typedef EMNLM_SIMILARITYSEARCH_METHOD SimilaritySearchMethodType;
};


template<typename _RealType, typename _PictureType>
class EMNonLocalMeansFilter : public EMNeighborhoodFilter<_RealType, _PictureType, EMNonLocalMeansFilterTrait<_RealType, _PictureType>,
                                                                                   NonLocalMeansFilter<_RealType, _PictureType, EMNonLocalMeansFilterTrait<_RealType, _PictureType> > > {
  typedef _RealType RealType;
  typedef aol::FullMatrix<RealType> MatrixType;
  typedef _PictureType PictureType;
  typedef EMNonLocalMeansFilterTrait<_RealType, _PictureType> FilterTrait;
  typedef typename FilterTrait::OptionsType OptionsType;
  typedef typename FilterTrait::BlockCollectionType BlockCollectionType;
public:
  
  static std::string getMethodIdentifier ( const OptionsType &Options ) {
    return getMethodPrefix ( Options ) + "-" + name ( );
  }
  
  static std::string getMethodPrefix ( const OptionsType &Options ) {
    return EMNeighborhoodFilter<RealType, PictureType, FilterTrait, NonLocalMeansFilter<RealType, PictureType, FilterTrait> >::getMethodPrefix ( Options );
  }
  
  static std::string name ( ) { return "EMNLM"; }
protected:
  void initializeBlocks ( const OptionsType &Options );
  
  friend class SimilaritySearchIterator<RealType, PictureType, FilterTrait>;
  friend class LocalSimilaritySearchIterator<RealType, PictureType, FilterTrait>;
  friend class GlobalSimilaritySearchIterator<RealType, PictureType, FilterTrait>;
  friend class PeriodicSimilaritySearchIterator<RealType, PictureType, FilterTrait, EMNonLocalMeansFilterTrait<_RealType, _PictureType> >;
};


#endif /* EMNONLOCALMEANSFILTER_H_ */