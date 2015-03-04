#include "nonLocalMeansFilter.h"
#include "emNonLocalMeansFilter.h"


template <typename _RealType, typename _PictureType, typename FilterTrait>
void NonLocalMeansFilter<_RealType, _PictureType, FilterTrait>::apply ( const OptionsType &Options, PictureType &Dest ) {
  OptionsType options ( Options );
  this->preprocess ( options );
  
  initialize ( options );
  denoise ( options );
  
  this->postprocess ( options );
  
  // Store estimate in Dest argument
  Dest.reallocate ( this->_NX, this->_NY );
  Dest = this->_estimate;
  
  // (Optional) console output
  if ( !this->_quietMode && options.groundTruth != NULL && (*options.groundTruth).getNumX ( ) == this->_estimate.getNumX ( ) && (*options.groundTruth).getNumY ( ) == this->_estimate.getNumY ( ) ) {
    if ( options.noiseType == NOISE_TYPE::GAUSSIAN ) std::cerr << "PSNR = " << aol::PSNR<RealType> ( *options.groundTruth, this->_estimate, 255 ) << " dB" << std::endl;
    else if ( options.noiseType == NOISE_TYPE::POISSON ) std::cerr << "PSNR = " << aol::PSNR<RealType> ( *options.groundTruth, this->_estimate ) << " dB" << std::endl;
  }
}


template <typename _RealType, typename _PictureType, typename FilterTrait>
void NonLocalMeansFilter<_RealType, _PictureType, FilterTrait>::initialize ( const OptionsType &Options ) {
  NeighborhoodFilter<RealType, PictureType, FilterTrait>::initialize ( Options );
  setParameters ( Options );
  initializeBlocks ( Options );
}


template <typename _RealType, typename _PictureType, typename FilterTrait>
void NonLocalMeansFilter<_RealType, _PictureType, FilterTrait>::initializeBlocks ( const OptionsType &/*Options*/ ) {
  this->_blocksInitial.fillFrom ( this->_preprocessedInput, this->_blockSize );
  this->_blocks = &this->_blocksInitial;
}


template <typename _RealType, typename _PictureType, typename FilterTrait>
void NonLocalMeansFilter<_RealType, _PictureType, FilterTrait>::setParameters ( const OptionsType &Options ) {
  if ( Options.blockSize > 0 ) {
    if ( Options.blockSize % 2 == 1 ) this->_blockSize = Options.blockSize;
    else throw aol::Exception ( "Specified block size must be odd!", __FILE__, __LINE__ );
  } else if ( Options.noiseType == NOISE_TYPE::GAUSSIAN && this->_stdDev > 0 ) {
    if ( Options.stdDev <= 15 ) this->_blockSize = 3;
    else if ( Options.stdDev <= 30 ) this->_blockSize = 5;
    else if ( Options.stdDev <= 45 ) this->_blockSize = 7;
    else if ( Options.stdDev <= 75 ) this->_blockSize = 9;
    else this->_blockSize = 11;
  } else this->_blockSize = 7;
  _blockOffset = ( this->_blockSize - 1 ) / 2;
  this->_blockSizeSqr = this->_blockSize * this->_blockSize;
  this->_X0 = _blockOffset;
  this->_Y0 = _blockOffset;
  this->_XEnd = this->_X0 + this->_NX - this->_blockSize + 1;
  this->_YEnd = this->_Y0 + this->_NY - this->_blockSize + 1;
  this->_effSize = ( this->_NX - this->_blockSize + 1 ) * ( this->_NY - this->_blockSize + 1 );
  
  if ( Options.similaritySearchMethod == NLM_SIMILARITYSEARCH_METHOD::LOCAL ) {
    if ( Options.searchWindowSize > 0 ) {
      if ( Options.searchWindowSize % 2 == 1 ) this->_searchWindowSize = Options.searchWindowSize;
      else throw aol::Exception ( "Specified search window size must be odd!", __FILE__, __LINE__ );
    } else if ( this->_stdDev > 0 && Options.noiseType == NOISE_TYPE::GAUSSIAN ) {
      if ( Options.stdDev <= 30 ) this->_searchWindowSize = 21;
      else this->_searchWindowSize = 35;
    } else this->_searchWindowSize = 35;
  }
  
  if ( Options.filterParameter > 0 ) _filterParameter = Options.filterParameter;
  else if ( this->_stdDev > 0 ) {
    _filterParameter = this->_stdDev;
    if ( Options.noiseType == NOISE_TYPE::GAUSSIAN ) {
      if ( Options.stdDev <= 30 ) _filterParameter *= 0.4;
      else if ( Options.stdDev <= 75 ) _filterParameter *= 0.35;
      else if ( Options.stdDev <= 100 ) _filterParameter *= 0.3;
      else throw aol::Exception ( "Removing Gaussian noise with standard deviation exceeding 100 not feasible!", __FILE__, __LINE__ );
    } else if ( Options.noiseType == NOISE_TYPE::POISSON && Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::ANSCOMBE ) {
      _filterParameter *= 0.35;
    } else throw aol::Exception ( "Automatic filter parameter settings are only provided for Gaussian noise (or Anscombe transformed Poisson noise)!", __FILE__, __LINE__ );
  } else throw aol::Exception ( "Neither filter parameter nor noise standard deviation was specified!", __FILE__, __LINE__ );
  _filterParameterSqr = aol::Sqr<RealType> ( _filterParameter );
  _stdDevSqr = aol::Sqr<RealType> ( this->_stdDev );
}


template <typename _RealType, typename _PictureType, typename FilterTrait>
void NonLocalMeansFilter<_RealType, _PictureType, FilterTrait>::denoise ( const OptionsType &Options ) {
  RealType dist, weight, normFactor;
  for ( GlobalReferenceBlockIterator<RealType, PictureType, FilterTrait> it ( *this ); it.notAtEnd ( ) ; ++it ) {
    this->_estimate.set ( *it, 0 );
    normFactor = 0;
    
    for ( this->_searchIt->reset ( *it ); this->_searchIt->notAtEnd ( ) ; ++(*this->_searchIt) ) {      
      (this->*this->_setBlockDistance) ( dist, *it, *(*this->_searchIt) );
      weight = getWeight ( dist );
    
      normFactor += weight;
      this->_estimate.set ( *it, this->_estimate.get ( *it ) + weight * this->_preprocessedInput.get ( *(*this->_searchIt) ) );
    }
    
    this->_estimate.set ( *it, this->_estimate.get ( *it ) / normFactor );
  }
}


template <typename _RealType, typename _PictureType, typename FilterTrait>
void NonLocalMeansFilter<_RealType, _PictureType, FilterTrait>::setWeights ( const OptionsType &Options, const aol::Vec2<short> &XRef, PictureType &Weights, const bool ReInitialize ) {
  if ( ReInitialize ) {
    OptionsType options ( Options );
    this->preprocess ( options );
    initialize ( options );
  }
  
  if ( !this->_quietMode )
    std::cerr << "NonLocalMeansFilter: Calculating weights for pixel (" << XRef[0] << "," << XRef[1] << ").." << std::endl;
  
  RealType dist;
  for ( this->_searchIt->reset ( XRef ); this->_searchIt->notAtEnd ( ) ; ++(*this->_searchIt) ) {
    (this->*this->_setBlockDistance) ( dist, XRef, *(*this->_searchIt) );
    Weights.set ( *(*this->_searchIt), getWeight ( dist ) );
    this->_searchIt->update ( dist );
  }
}


template <typename _RealType, typename _PictureType, typename FilterTrait>
_RealType NonLocalMeansFilter<_RealType, _PictureType, FilterTrait>::getPUniformityOptimalFilterParameter ( const OptionsType &Options ) {
  if ( Options.noiseType == NOISE_TYPE::POISSON )
    throw aol::Exception ( "P-values are not properly defined for Poisson noise statistics!", __FILE__, __LINE__ );
  
  if ( !this->_quietMode )
    std::cerr << "NonLocalMeansFilter: searching for p-uniformity optimal filter parameter.." << std::endl;
  
  const int numBins = 100;
  const RealType p = 0.5, alpha = 0.05;
  RealType lBound = 0, rBound = aol::Sqr<RealType> ( this->_input.getMaxValue ( ) ), h0 = 0.5 * ( lBound + rBound ), h1 = h0, chiSqr0 = numBins * 2, chiSqr1 = chiSqr0;
  int numPixels = 128, numItInner = 0, numItOuter = 1, pValuesQuantileThreshold;
  aol::RandomGenerator randomGenerator;
  randomGenerator.randomize ( );
  
  if ( !this->_quietMode )
    std::cerr << "Iteration " << numItOuter << "." << numItInner << ": h0=" << h0 << " in [" << lBound << ", " << rBound << "]" << std::endl;
  
  do {
    chiSqr0 = chiSqr1;
    
    aol::MultiVector<short> pixels ( numPixels, 2 );
    short numGeneratedPixels = 0;
    do {
      aol::Vec2<short> pixel ( randomGenerator.rInt ( this->_blockOffset, this->_NX - this->_blockOffset ),
                               randomGenerator.rInt ( this->_blockOffset, this->_NY - this->_blockOffset ) );
      bool isNewPixel = true;
      for ( int k=0; k<numGeneratedPixels ; ++k ) {
        if ( pixels[k][0] == pixel[0] && pixels[k][1] == pixel[1] ) {
          isNewPixel = false;
          break;
        }
      }
      if ( isNewPixel ) {
        pixels.set ( numGeneratedPixels, pixel );
        ++numGeneratedPixels;
      }
    } while ( numGeneratedPixels < numPixels );
    pValuesQuantileThreshold = DistributionTester<RealType>::getPValuesQuantileThreshold ( numPixels, p, alpha );
    do {
      h0 = h1;
      
      PictureType means ( this->_NX, this->_NY );
      OptionsType options ( Options );
      options.pixels = pixels;
      options.filterParameter = h1;
      apply ( options, means );
      
      GaussianNoiseImageAnalyzer<RealType, PictureType> noiseImageAnalyzer ( this->_input, means, numBins, this->_outputDir );
      const std::pair<int, int> pValuesInSymmetricQuantiles ( noiseImageAnalyzer.getNumPValuesInSymmetricQuantiles ( p, pixels ) );
      chiSqr1 = noiseImageAnalyzer.getChiSquareOfUniformDistributionOfPValues ( pixels );
      
      if ( pValuesInSymmetricQuantiles.first > pValuesQuantileThreshold ) {
        h1 = 0.5 * ( lBound + h0 );
        rBound = h0;
      } else if ( pValuesInSymmetricQuantiles.second > pValuesQuantileThreshold ) {
        h1 = 0.5 * ( h0 + rBound );
        lBound = h0;
      }
      numItInner++;
      
      if ( !this->_quietMode )
        std::cerr << "Iteration " << numItOuter << "." << numItInner << ": h1=" << h1 << " in [" << lBound << ", " << rBound << "]; (pv<, pv>)=(" << pValuesInSymmetricQuantiles.first << ", " << pValuesInSymmetricQuantiles.second << "); pv_thres=" << pValuesQuantileThreshold << "; chi^2=" << chiSqr1 << std::endl;
    } while ( aol::Abs<RealType> ( h1 - h0 ) / h1 > 1e-3 );
    numItInner = 0;
    numItOuter++;
    numPixels *= 2;
  } while ( numPixels < this->_input.size ( ) / 4 );
  
  if ( !this->_quietMode )
    std::cerr << "Finished binary search. h_opt=" << h1 << std::endl;
  
  return h1;
}


template <typename _RealType, typename _PictureType, typename FilterTrait>
std::string NonLocalMeansFilter<_RealType, _PictureType, FilterTrait>::getMethodPrefix ( const OptionsType &Options ) {
  std::stringstream ss;
  ss << NeighborhoodFilter<RealType, PictureType, FilterTrait>::getMethodPrefix ( Options );
  ss << "-" << Options.filterParameter;
  return ss.str ( );
}


template class NonLocalMeansFilter<double, qc::ScalarArray<double, qc::QC_2D> >;
template class NonLocalMeansFilter<double, qc::ScalarArray<double, qc::QC_2D>, EMNonLocalMeansFilterTrait<double, qc::ScalarArray<double, qc::QC_2D> > >;
