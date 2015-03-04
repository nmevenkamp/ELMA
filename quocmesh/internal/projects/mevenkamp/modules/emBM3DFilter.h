#ifndef EMBM3DFILTER_H_
#define EMBM3DFILTER_H_

#include "emNeighborhoodFilter.h"
#include "bm3dFilter.h"


class EMBM3D_SIMILARITYSEARCH_METHOD : public EM_SIMILARITYSEARCH_METHOD {
public:
  static const int UNIFORM_PERIODIC = 3;
  
  static const int NUM = 4;
  
  static std::string getIdentifier ( const int SimilaritySearchMethod ) {
    if ( SimilaritySearchMethod < UNIFORM_PERIODIC ) return EM_SIMILARITYSEARCH_METHOD::getIdentifier ( SimilaritySearchMethod );
    else if ( SimilaritySearchMethod == UNIFORM_PERIODIC ) return "piU";
    else throw aol::Exception ( "EMBM3D: Did not recognize similarity search method", __FILE__, __LINE__ );
  }
  
  static std::string toString ( const int SimilaritySearchMethod ) {
    if ( SimilaritySearchMethod < UNIFORM_PERIODIC ) return EM_SIMILARITYSEARCH_METHOD::toString ( SimilaritySearchMethod );
    else if ( SimilaritySearchMethod == UNIFORM_PERIODIC ) return "Uniform periodic";
    else throw aol::Exception ( "EMBM3D: Did not recognize similarity search method", __FILE__, __LINE__ );
  }
};

class EMBM3D_PROFILE : public BM3D_PROFILE {
public:
  static std::string getIdentifier ( const int Profile ) {
    return BM3D_PROFILE::getIdentifier ( Profile );
  }
  
  static std::string toString ( const int Profile ) {
    return BM3D_PROFILE::toString ( Profile );
  }
};

class EMBM3D_SCANDISTORTIONCORRECTION_METHOD : public EM_SCANDISTORTIONCORRECTION_METHOD {
public:
  static std::string getIdentifier ( const int ScanDistortionCorrectionMethod ) {
    return EM_SCANDISTORTIONCORRECTION_METHOD::getIdentifier ( ScanDistortionCorrectionMethod );
  }
};

class EMBM3D_ESTIMATION_METHOD : public BM3D_ESTIMATION_METHOD {
public:
  static std::string getIdentifier ( const int EstimationMethod ) {
    return BM3D_ESTIMATION_METHOD::getIdentifier ( EstimationMethod );
  }
};


template <typename _RealType, typename _PictureType>
struct EMBM3DOptions : public EMNHFilterOptions<_RealType, _PictureType, BM3DOptions<_RealType, _PictureType> > {
  EMBM3DOptions ( _PictureType &Input, const std::string OutputDir = "", const bool Quiet = true )
    : EMNHFilterOptions<_RealType, _PictureType, BM3DOptions<_RealType, _PictureType> > ( Input, OutputDir, Quiet ) { }
};


template <typename _RealType, typename _PictureType>
class EMBM3DFilterTrait : public NHFilterTrait<_RealType, _PictureType> {
public:
  typedef EMBM3DOptions<_RealType, _PictureType> OptionsType;
  typedef ShiftedBlockCollection<_RealType, _PictureType> BlockCollectionType;
  typedef EMBM3D_SIMILARITYSEARCH_METHOD SimilaritySearchMethodType;
};


template <typename _RealType, typename _PictureType>
class EMBM3DFilter : public EMNeighborhoodFilter<_RealType, _PictureType, EMBM3DFilterTrait<_RealType, _PictureType>,
                                                                          BM3DFilter<_RealType, _PictureType, EMBM3DFilterTrait<_RealType, _PictureType> > > {
  typedef _RealType RealType;
  typedef aol::Vector<RealType> VectorType;
  typedef aol::FullMatrix<RealType> MatrixType;
  typedef _PictureType PictureType;
  typedef qc::MultiArray<RealType, qc::QC_2D, 3> ColoredPictureType;
  typedef EMBM3DFilterTrait<RealType, PictureType> FilterTrait;
  typedef typename FilterTrait::OptionsType OptionsType;
  typedef typename FilterTrait::BlockCollectionType BlockCollectionType;
public:
  void apply ( const OptionsType &Options, PictureType &Dest ) {
    BM3DFilter<RealType, PictureType, FilterTrait>::apply ( Options, Dest );
  }
  
  void setSimilaritySearchGrid ( const OptionsType &Options, ColoredPictureType &Dest,
                                 const aol::Vec2<short> &XRef = aol::Vec2<short> ( -1, -1 ) );
  
  static std::string getMethodIdentifier ( const OptionsType &Options ) {
    return getMethodPrefix ( Options ) + "-" + name ( );
  }
  
  static std::string getMethodPrefix ( const OptionsType &Options ) {
    return EMNeighborhoodFilter<RealType, PictureType, FilterTrait, BM3DFilter<RealType, PictureType, FilterTrait> >::getMethodPrefix ( Options );
  }
  
  static std::string name ( ) { return "EMBM3D"; }
protected:
  void initializeOnce ( const OptionsType &Options ) {
    BM3DFilter<RealType, PictureType, FilterTrait>::initializeOnce ( Options );
    this->setScanDistortionCorrectionParameters ( Options );
  }
  
  void fillBlocksFromPreviousEstimate ( const OptionsType &Options );
  
  void aggregate ( const OptionsType &Options, const Blocks3DAndWeight<RealType> &Blocks3D, const aol::MultiVector<short> &Sx, const short Nz );
  
  void setEstimateFromBuffers ( const OptionsType &Options );
  
  void inpaint ( const aol::Vec2<short> &X );
  
  void blockMatching ( const OptionsType &Options, const aol::Vec2<short> &XRef, aol::MultiVector<short> &Sx, VectorType &BlockDistances, int &Nz );
                                                                            
  void setSimilaritySearchIterator ( const OptionsType &Options ) {
    if ( Options.similaritySearchMethod < EMBM3D_SIMILARITYSEARCH_METHOD::UNIFORM_PERIODIC )
      EMNeighborhoodFilter<_RealType, _PictureType, EMBM3DFilterTrait<_RealType, _PictureType>,
                                                    BM3DFilter<_RealType, _PictureType, EMBM3DFilterTrait<_RealType, _PictureType> > >::setSimilaritySearchIterator ( Options );
    else {
      if ( Options.similaritySearchMethod == EMBM3D_SIMILARITYSEARCH_METHOD::UNIFORM_PERIODIC )
        this->_searchIt = new PeriodicSimilaritySearchIterator<RealType, PictureType, EMBM3DFilterTrait<_RealType, _PictureType>,
                                                                                      BM3DFilter<_RealType, _PictureType, EMBM3DFilterTrait<_RealType, _PictureType> > > ( *this );
      else
        throw aol::Exception ( "EMBM3DFilter: Did not recognize similarity search method!", __FILE__, __LINE__ );
    }
  }

  
  friend class SimilaritySearchIterator<RealType, PictureType, FilterTrait>;
  friend class LocalSimilaritySearchIterator<RealType, PictureType, FilterTrait>;
  friend class GlobalSimilaritySearchIterator<RealType, PictureType, FilterTrait>;
  friend class PeriodicSimilaritySearchIterator<RealType, PictureType, FilterTrait, BM3DFilter<_RealType, _PictureType, EMBM3DFilterTrait<_RealType, _PictureType> > >;
};


#endif /* BM3DFILTER_H_ */