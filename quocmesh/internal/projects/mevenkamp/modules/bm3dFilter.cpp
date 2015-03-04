#include "bm3dFilter.h"
#include "referenceBlockIterators.h"
#include "similaritySearchIterators.h"
#include "emBM3DFilter.h"


template <typename _RealType, typename _PictureType, typename FilterTrait>
BM3DFilter<_RealType, _PictureType, FilterTrait>::BM3DFilter ( )
  : _unitaryTransform2D ( 1, "dct" ),
    _htDenoisingInLocal2D1DTransformDomainOp ( ), _wienerDenoisingInLocal2D1DTransformDomainOp ( ) { }

template <typename _RealType, typename _PictureType, typename FilterTrait>
void BM3DFilter<_RealType, _PictureType, FilterTrait>::apply ( const OptionsType &Options, PictureType &Dest ) {
  this->setCtrlCHandler();
  
  OptionsType options ( Options );
  this->preprocess ( options );
  
  // Denoise
  initializeOnce ( options );
  
  options.method = BM3D_ESTIMATION_METHOD::HT;
  if ( !this->wantsInterrupt ( ) ) denoise ( options );
  
    _numAggregatesHT = _numAggregates;
    _cumulatedDistancesHT = _distances;
  
  options.method = BM3D_ESTIMATION_METHOD::WIENER;
  if ( !this->wantsInterrupt ( ) ) denoise ( options );

    _numAggregatesWiener = _numAggregates;
    _cumulatedDistancesWiener = _distances;
  
  
  this->postprocess ( options );
  
  // Store estimate in Dest argument
  Dest.reallocate ( this->_NX, this->_NY );
  Dest = this->_estimate;
  
  // (Optional) console output
  if ( !this->_quietMode && options.groundTruth != NULL && (*options.groundTruth).getNumX ( ) == this->_estimate.getNumX ( ) && (*options.groundTruth).getNumY ( ) == this->_estimate.getNumY ( ) ) {
    if ( options.noiseType == NOISE_TYPE::GAUSSIAN ) std::cerr << "PSNR = " << aol::PSNR<RealType> ( *options.groundTruth, this->_estimate, 255 ) << " dB" << std::endl;
    else if ( options.noiseType == NOISE_TYPE::POISSON ) std::cerr << "PSNR = " << aol::PSNR<RealType> ( *options.groundTruth, this->_estimate ) << " dB" << std::endl;
  }
  
  this->unsetCtrlCHandler();
}


template <typename _RealType, typename _PictureType, typename FilterTrait>
void BM3DFilter<_RealType, _PictureType, FilterTrait>::initializeOnce ( const OptionsType &Options ) {
  this->initialize ( Options );
  
  this->_X0 = 0;
  this->_Y0 = 0;
  _mapper.resize ( this->_NX, this->_NY );
  _eBuff.reallocate ( this->_NX, this->_NY );
  _wBuff.reallocate ( this->_NX, this->_NY );
  
  // For testing purposes
  _numMatchedBlocks.reallocate ( this->_NX, this->_NY );
  _numAggregates.reallocate ( this->_NX, this->_NY );
  _numAggregatesHT.reallocate ( this->_NX, this->_NY );
  _numAggregatesWiener.reallocate ( this->_NX, this->_NY );
  _distances.reallocate ( this->_NX, this->_NY );
  _cumulatedDistancesHT.reallocate ( this->_NX, this->_NY );
  _cumulatedDistancesWiener.reallocate ( this->_NX, this->_NY );
}

template <typename _RealType, typename _PictureType, typename FilterTrait>
void BM3DFilter<_RealType, _PictureType, FilterTrait>::initializeIteration ( const OptionsType &Options ) {
  _wBuff.setZero ( );
  _eBuff.setZero ( );
  _numMatchedBlocks.setZero ( );
  _numAggregates.setZero ( );
  _distances.setZero ( );
  setParameters ( Options );
  setOperators ( Options );
}


template <typename _RealType, typename _PictureType, typename FilterTrait>
void BM3DFilter<_RealType, _PictureType, FilterTrait>::fillBlocksFromPreviousEstimate ( const OptionsType &Options ) {
  if ( Options.noiseType == NOISE_TYPE::POISSON && Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::MAXIMUMLIKELIHOOD && Options.method == BM3D_ESTIMATION_METHOD::HT ) {
    this->_blocksInitial.fillFrom ( this->_input, this->_blockSize );
    this->_blocks = &this->_blocksInitial;
  } else {
    if ( Options.method == BM3D_ESTIMATION_METHOD::HT ) {
      this->_blocksInitial.fillFrom ( this->_preprocessedInput, this->_blockSize );
      this->_blocks = &this->_blocksInitial;
    } else {
      this->_blocksEstimate.fillFrom ( this->_estimate, this->_blockSize );
      this->_blocks = &this->_blocksEstimate;
    }
  }
}


template <typename _RealType, typename _PictureType, typename FilterTrait>
void BM3DFilter<_RealType, _PictureType, FilterTrait>::applyLinearUnitary2DTransformsToBlocks ( const OptionsType &Options ) {
  BlockType tmpBlockArg ( this->_blockSize, this->_blockSize ), tmpBlockDest ( this->_blockSize, this->_blockSize );
  _blocksInitialTransformed.resize ( this->_NX, this->_NY, this->_blockSize );
  if ( Options.method == BM3D_ESTIMATION_METHOD::HT ) this->_progressBar->setText ( "BM3D: Transforming initial 2D blocks (step 1/5)" );
  else if ( Options.method == BM3D_ESTIMATION_METHOD::WIENER ) this->_progressBar->setText ( "BM3D: Transforming initial 2D blocks (step 3/5)" );
  this->_progressBar->start ( this->_XEnd );
  for ( short x=0; x<this->_XEnd && !this->wantsInterrupt ( ); ++x ) {
    for ( short y=0; y<this->_YEnd && !this->wantsInterrupt ( ); ++y ) {
      if ( Options.noiseType == NOISE_TYPE::POISSON && Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::MAXIMUMLIKELIHOOD
        && Options.method == BM3D_ESTIMATION_METHOD::HT ) {
        for ( short xBlock=0; xBlock<this->_blockSize ; ++xBlock )
          for ( short yBlock=0; yBlock<this->_blockSize ; ++yBlock )
            tmpBlockArg.set ( xBlock, yBlock, this->_preprocessedInput.get ( x + xBlock, y + yBlock ) );
      } else this->_blocksInitial.copyTo ( x, y, tmpBlockArg );
      _unitaryTransform2D.apply ( tmpBlockArg, tmpBlockDest );
      _blocksInitialTransformed.copyFrom ( x, y, tmpBlockDest );
    }
    (*this->_progressBar)++;
  }
  this->_progressBar->finish ( );
  if ( Options.method == BM3D_ESTIMATION_METHOD::WIENER ) {
    this->_progressBar->setText ( "BM3D: Transforming estimate 2D blocks (step 4/5)" );
    this->_progressBar->start ( this->_XEnd );
    _blocksEstimateTransformed.resize ( this->_NX, this->_NY, this->_blockSize );
    for ( short x=0; x<this->_XEnd && !this->wantsInterrupt ( ) ; ++x ) {
      for ( short y=0; y<this->_YEnd && !this->wantsInterrupt ( ) ; ++y ) {
        _blocksEstimate.copyTo ( x, y, tmpBlockArg );
        _unitaryTransform2D.apply ( tmpBlockArg, tmpBlockDest );
        _blocksEstimateTransformed.copyFrom ( x, y, tmpBlockDest );
      }
      (*this->_progressBar)++;
    }
  }
  this->_progressBar->finish ( );
}


template <typename _RealType, typename _PictureType, typename FilterTrait>
void BM3DFilter<_RealType, _PictureType, FilterTrait>::denoise ( const OptionsType &Options ) {
  initializeIteration ( Options );
  fillBlocksFromPreviousEstimate ( Options );
  applyLinearUnitary2DTransformsToBlocks ( Options );
  
  // Process all reference blocks in a sliding manner
  int Nz;
  aol::MultiVector<short> Sx ( this->_effSize, 2 );
  VectorType blockDistances ( this->_effSize );
  Blocks3DInitialAndEstimate<RealType> blocks3DInitialAndEstimate ( this->_blockSize, this->_blockSize, this->_N2 );
  Blocks3DAndWeight<RealType> blocks3DAndWeight ( this->_blockSize, this->_blockSize, this->_N2 );
  
  if ( Options.method == BM3D_ESTIMATION_METHOD::HT ) this->_progressBar->setText ( "BM3D: Computing initial estimate (step 2/5)" );
  else if ( Options.method == BM3D_ESTIMATION_METHOD::WIENER ) this->_progressBar->setText ( "BM3D: Computing final estimate (step 5/5)" );
  for ( GlobalReferenceBlockIterator<RealType, PictureType, FilterTrait> it ( *this ); it.notAtEnd ( ) && !this->wantsInterrupt ( ) ; ++it ) {
    blockMatching ( Options, *it, Sx, blockDistances, Nz );
    blockStacking ( Options, Sx, Nz, blocks3DInitialAndEstimate );
    _op3D->applyUsingOnly1DTransform ( blocks3DInitialAndEstimate, blocks3DAndWeight );    
    aggregate ( Options, blocks3DAndWeight, Sx, Nz );
  }
  
  setEstimateFromBuffers ( Options );
}


template <typename _RealType, typename _PictureType, typename FilterTrait>
void BM3DFilter<_RealType, _PictureType, FilterTrait>::blockMatching ( const OptionsType &Options, const aol::Vec2<short> &XRef,
                                                          aol::MultiVector<short> &Sx, VectorType &BlockDistances, int &Nz ) {
  // Calculate distances to all blocks in the search region
  int z = 0;
  RealType dist;
  for ( this->_searchIt->reset ( XRef ); this->_searchIt->notAtEnd ( ) ; ++(*this->_searchIt) ) {
    (this->*this->_setBlockDistance) ( dist, XRef, *(*this->_searchIt) );
    if ( this->_searchIt->isCurFinal ( ) && dist <= _tauMatch ) {
      BlockDistances[z] = dist;
      Sx[z][0] = (*(*this->_searchIt))[0]; Sx[z][1] = (*(*this->_searchIt))[1];
      ++z;
    }
    this->_searchIt->update ( dist );
  }
  Nz = z;
  _numMatchedBlocks.set ( XRef, Nz );
  
  // If Nz > _N2 or Nz is not a power of 2 drop coordinates with highest distances until |Sx| = Nz = _N2 or |Sx| = Nz is a power of 2
  short finalNz = finalNz = aol::Min<int> ( aol::Pow ( 2, floor ( log2 ( Nz ) ) ), _N2 );
  if ( Nz > finalNz ){
    aol::Vec2<short> tmpPos;
    RealType tmpDist;
    int zMin;
    for ( int z=0; z<finalNz ; ++z ) {
      zMin = z;
      for ( int z2=z+1; z2<Nz ; ++z2 )
        if ( BlockDistances[z2] < BlockDistances[zMin] ) zMin = z2;
      tmpPos.set ( Sx[z][0], Sx[z][1] );
      tmpDist = BlockDistances[z];
      Sx[z][0] = Sx[zMin][0]; Sx[z][1] = Sx[zMin][1];
      BlockDistances[z] = BlockDistances[zMin];
      Sx[zMin][0] = tmpPos[0]; Sx[zMin][1] = tmpPos[1];
      BlockDistances[zMin] = tmpDist;
    }
    Nz = finalNz;
  }
  
  for ( int z=0; z<finalNz ; ++z )
    _distances[_mapper.getGlobalIndex ( XRef[0], XRef[1] )] += BlockDistances[z];
}


template <typename _RealType, typename _PictureType, typename FilterTrait>
void BM3DFilter<_RealType, _PictureType, FilterTrait>::blockStacking ( const OptionsType &Options, const aol::MultiVector<short> &Sx, const short Nz,
                                                          Blocks3DInitialAndEstimate<RealType> &Blocks3DInitialAndEstimate ) {
  Blocks3DInitialAndEstimate.getInitialReference ( ).reallocate ( this->_blockSize, this->_blockSize, Nz );
  Blocks3DInitialAndEstimate.getEstimateReference ( ).reallocate ( this->_blockSize, this->_blockSize, Nz );
  for ( int z=0; z<Nz ; ++z )
    for ( short xBlock=0; xBlock<this->_blockSize ; ++xBlock )
      for ( short yBlock=0; yBlock<this->_blockSize ; ++yBlock )
        Blocks3DInitialAndEstimate.setInitial ( xBlock, yBlock, z, _blocksInitialTransformed.get ( Sx[z][0], Sx[z][1], xBlock, yBlock ) );
  
  if ( Options.method == BM3D_ESTIMATION_METHOD::WIENER ) {
    for ( int z=0; z<Nz ; ++z )
      for ( short xBlock=0; xBlock<this->_blockSize ; ++xBlock )
        for ( short yBlock=0; yBlock<this->_blockSize ; ++yBlock )
          Blocks3DInitialAndEstimate.setEstimate ( xBlock, yBlock, z, _blocksEstimateTransformed.get ( Sx[z][0], Sx[z][1], xBlock, yBlock ) );
  }
}


template <typename _RealType, typename _PictureType, typename FilterTrait>
void BM3DFilter<_RealType, _PictureType, FilterTrait>::aggregate ( const OptionsType &/*Options*/, const Blocks3DAndWeight<RealType> &Blocks3D, const aol::MultiVector<short> &Sx, const short Nz ) {
  RealType weight;
  int flatIdx;
  for ( short xBlock=0; xBlock<this->_blockSize ; ++xBlock ) {
    for ( short yBlock=0; yBlock<this->_blockSize ; ++yBlock ) {
      weight = Blocks3D.weight ( ) * _WWin2D.get ( xBlock, yBlock );
      for ( short z=0; z<Nz ; ++z ) {        
        flatIdx = _mapper.getGlobalIndex ( Sx[z][0] + xBlock, Sx[z][1] + yBlock );
        _eBuff[flatIdx] += weight * Blocks3D.get ( xBlock, yBlock, z );
        _wBuff[flatIdx] += weight;
        ++_numAggregates[flatIdx];
      }
    }
  }
}


template <typename _RealType, typename _PictureType, typename FilterTrait>
void BM3DFilter<_RealType, _PictureType, FilterTrait>::setEstimateFromBuffers ( const OptionsType &/*Options*/ ) {
  for ( int k=0; k<this->_estimate.size ( ) ; ++k ) {
    if ( _wBuff[k] == 0 ) this->_estimate[k] = 0;
    else this->_estimate[k] = _eBuff[k] / _wBuff[k];
  }
}


template <typename _RealType, typename _PictureType, typename FilterTrait>
void BM3DFilter<_RealType, _PictureType, FilterTrait>::setParameters ( const OptionsType &Options ) {
  this->_blockSize = Options.blockSize;
  this->_blockSizeSqr = this->_blockSize * this->_blockSize;
  
  // Initially set all parameters for Normal Profile
  this->_searchWindowSize = 39;
  this->_refStep = 3;
  _beta = 2.0;
  if ( Options.method == BM3D_ESTIMATION_METHOD::HT ) {
    _N2 = 16;
    _tauMatch = 3000;
  } else if ( Options.method == BM3D_ESTIMATION_METHOD::WIENER ) {
    _N2 = 32;
    _tauMatch = 400;
  }
  
  // Adjust parameters according to actually chosen profile
  if ( Options.profile == BM3D_PROFILE::LC ) {
    if ( Options.method == BM3D_ESTIMATION_METHOD::HT )
      this->_refStep = 6;
    else if ( Options.method == BM3D_ESTIMATION_METHOD::WIENER ) {
      this->_refStep = 5;
      _N2 = 16;
    }
    this->_searchWindowSize = 25;
  }
  
  if ( Options.profile == BM3D_PROFILE::HIGH ) {
    if ( Options.method == BM3D_ESTIMATION_METHOD::HT )
      _beta = 2.5;
    else if ( Options.method == BM3D_ESTIMATION_METHOD::WIENER )
      _beta = 1.5;
    this->_refStep = 2;
  }
  
  if ( Options.method == BM3D_ESTIMATION_METHOD::HT && Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::MAXIMUMLIKELIHOOD ) {
    _tauMatch = -log(0.55); // TODO
    _tauMatch *= this->_blockSizeSqr;
  } else
    _tauMatch *= this->_blockSizeSqr / aol::Sqr<RealType> ( 255.0 );
  
  this->_XEnd = this->_NX - this->_blockSize + 1;
  this->_YEnd = this->_NY - this->_blockSize + 1;
  this->_effSize = ( this->_NX - this->_blockSize + 1 ) * ( this->_NY - this->_blockSize + 1 );
  
  setKaiserWindow ( _WWin2D, this->_blockSize, _beta );
}


template <typename _RealType, typename _PictureType, typename FilterTrait>
void BM3DFilter<_RealType, _PictureType, FilterTrait>::setOperators ( const OptionsType &Options ) {
  if ( Options.method == BM3D_ESTIMATION_METHOD::HT ) {
    if ( Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::MAXIMUMLIKELIHOOD )
      this->_setBlockDistance = &BM3DFilter<RealType, PictureType, FilterTrait>::setBlockDistancePoissonLikelihoodRatio;
    else
      this->_setBlockDistance = &BM3DFilter<RealType, PictureType, FilterTrait>::setBlockDistanceL2NormSqr;
    
    // Set 3D operator
    _unitaryTransform2D.setTransform ( this->_blockSize, "bior1.5" );
    _htDenoisingInLocal2D1DTransformDomainOp.setTransforms ( this->_blockSize, "bior1.5", _N2, "haar" );

    RealType lambdaThr3D = 2.7;
    if ( Options.profile == BM3D_PROFILE::HIGH ) lambdaThr3D = 2.5;
    const RealType threshold3D = lambdaThr3D * this->_stdDev;
    _htDenoisingInLocal2D1DTransformDomainOp.setThreshold ( threshold3D );
    
    _op3D = &_htDenoisingInLocal2D1DTransformDomainOp;
  } else if ( Options.method == BM3D_ESTIMATION_METHOD::WIENER ) {
    this->_setBlockDistance = &BM3DFilter<RealType, PictureType, FilterTrait>::setBlockDistanceL2NormSqr;

    // Set 3D operator
    _unitaryTransform2D.setTransform ( this->_blockSize, "dct" );
    _wienerDenoisingInLocal2D1DTransformDomainOp.setTransforms ( this->_blockSize, "dct", _N2, "haar" );
    _wienerDenoisingInLocal2D1DTransformDomainOp.setNoiseStandardDeviation ( this->_stdDev );
    _op3D = &_wienerDenoisingInLocal2D1DTransformDomainOp;
  }
}


template <typename _RealType, typename _PictureType, typename FilterTrait>
void BM3DFilter<_RealType, _PictureType, FilterTrait>::setSimilaritySearchGrid ( const OptionsType &Options, ColoredPictureType &Dest, const aol::Vec2<short> &XRef ) {
  OptionsType options ( Options );
  
  // Initialize PNG with intensities from input
  Dest.setZero ( );
  Dest[0] = Options.input; Dest[1] = Options.input; Dest[2] = Options.input;
  RealType maxVal = Options.input.getMaxValue ( );
  
  this->preprocess ( options );

  initializeOnce ( options );
  options.method = BM3D_ESTIMATION_METHOD::HT;
  initializeIteration ( options );
  fillBlocksFromPreviousEstimate ( Options );

  PictureType tmp ( Options.input );
  tmp.set ( XRef, 0 );
  
  // Perform block-matching and in output PNG mark all corners of iterated blocks in the search region red
  int Nz;
  aol::MultiVector<short> Sx ( this->_effSize, 2 );
  int z = 0;
  VectorType blockDistances ( this->_effSize );
  RealType dist;
  for ( this->_searchIt->reset ( XRef ); this->_searchIt->notAtEnd ( ) ; ++(*this->_searchIt) ) {
    (this->*this->_setBlockDistance) ( dist, XRef, *(*this->_searchIt) );
    _distances.set ( *(*this->_searchIt), dist );
    
    Dest[1].set ( *(*this->_searchIt), 0 );
    if ( this->_searchIt->isCurFinal ( ) ) {
      Dest[0].set ( *(*this->_searchIt), maxVal );
      Dest[2].set ( *(*this->_searchIt), 0 );
    } else {
      Dest[0].set ( *(*this->_searchIt), 0 );
      Dest[2].set ( *(*this->_searchIt), maxVal );
    }
    
    if ( this->_searchIt->isCurFinal ( ) && dist <= _tauMatch ) {
      blockDistances[z] = dist;
      Sx[z][0] = (*(*this->_searchIt))[0]; Sx[z][1] = (*(*this->_searchIt))[1];
      ++z;
    }
  
    this->_searchIt->update ( dist );
  }
  Nz = z;
  
  if ( !this->_quietMode && this->_outputDir != "" ) {
    aol::Vector<RealType> blockDists;
    for ( int z=0; z<Nz ; ++z )
      blockDists.pushBack ( blockDistances[z] );
    std::stringstream ss;
    ss << this->_outputDir << "/blockDistances_HT_" << this->getMethodIdentifier ( Options );
    std::vector<std::pair<RealType, int> > histo;
    blockDists.createHistogramOfValues ( histo, 256 );
    aol::plotHistogram<RealType> ( histo, ss.str ( ).c_str ( ), true );
  }
  
  // If Nz > _N2 or Nz is not a power of 2 drop coordinates with highest distances until |Sx| = Nz = _N2 or |Sx| = Nz is a power of 2
  const short finalZ = aol::Min<int> ( aol::Pow ( 2, floor ( log2 ( Nz ) ) ), _N2 );
  if ( Nz > finalZ ) {
    aol::Vec2<short> tmpPos;
    RealType tmpDist;
    int zMin;
    for ( int z=0; z<finalZ ; ++z ) {
      zMin = z;
      for ( int z2=z+1; z2<Nz ; ++z2 )
        if ( blockDistances[z2] < blockDistances[zMin] ) zMin = z2;
      tmpPos.set ( Sx[z][0], Sx[z][1] );
      tmpDist = blockDistances[z];
      Sx[z][0] = Sx[zMin][0]; Sx[z][1] = Sx[zMin][1];
      blockDistances[z] = blockDistances[zMin];
      Sx[zMin][0] = tmpPos[0]; Sx[zMin][1] = tmpPos[1];
      blockDistances[zMin] = tmpDist;
    }
    Nz = finalZ;
  }
  
  if ( !this->_quietMode && this->_outputDir != "" ) {
    RealType finalBlocksCumulativeDist = 0;
    for ( short z=0; z<Nz ; ++z ) {
      Dest[0].set ( Sx[z][0], Sx[z][1], maxVal );
      Dest[1].set ( Sx[z][0], Sx[z][1], maxVal );
      Dest[2].set ( Sx[z][0], Sx[z][1], 0 );
      finalBlocksCumulativeDist += blockDistances[z];
      
      PictureType block ( this->_blockSize, this->_blockSize );
      for ( short x=0; x<this->_blockSize ; ++x )
        for ( short y=0; y<this->_blockSize ; ++y )
          block.set ( x, y, this->_input.get ( Sx[z][0] + x, Sx[z][1] + y ) );
      std::stringstream ss;
      ss << this->_outputDir << "/blocks_" << this->getMethodIdentifier ( Options );
      aol::makeDirectory ( ss.str ( ).c_str ( ) );
      ss << "/blocks_" << z << ".q2bz";
      block.save ( ss.str ( ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    }
    std::cerr << "Final cumulative block distances: " << finalBlocksCumulativeDist << std::endl;
    
    // In output PNG mark corner of reference block green
    Dest[0].set ( XRef, 0 );
    Dest[1].set ( XRef, maxVal );
    Dest[2].set ( XRef, 0 );
  
    std::stringstream ss;
    ss << this->_outputDir << "/distances_" << this->getMethodIdentifier ( Options ) << ".q2bz";
    _distances.save ( ss.str ( ).c_str ( ), qc::PGM_DOUBLE_BINARY );
  }
}


template <typename _RealType, typename _PictureType, typename FilterTrait>
std::string BM3DFilter<_RealType, _PictureType, FilterTrait>::getMethodPrefix ( const OptionsType &Options ) {
  std::string methodPrefix = BM3D_PROFILE::getIdentifier ( Options.profile );
  methodPrefix += "-" + NeighborhoodFilter<RealType, PictureType, FilterTrait>::getMethodPrefix ( Options );
  return methodPrefix;
}


template class BM3DFilter<double, qc::ScalarArray<double, qc::QC_2D> >;
template class BM3DFilter<double, qc::ScalarArray<double, qc::QC_2D>, EMBM3DFilterTrait<double, qc::ScalarArray<double, qc::QC_2D> > >;