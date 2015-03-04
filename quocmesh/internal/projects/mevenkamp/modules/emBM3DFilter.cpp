#include "emBM3DFilter.h"


template <typename _RealType, typename _PictureType>
void EMBM3DFilter<_RealType, _PictureType>::fillBlocksFromPreviousEstimate ( const OptionsType &Options ) {
  if ( Options.noiseType == NOISE_TYPE::POISSON && Options.poissonNoiseAdaptation == POISSON_NOISE_ADAPTATION::MAXIMUMLIKELIHOOD && Options.method == BM3D_ESTIMATION_METHOD::HT ) {
    this->fillBlocks ( this->_blocksInitial, this->_input, this->_blockSize, Options );
    this->_blocks = &this->_blocksInitial;
  } else {
    if ( Options.method == BM3D_ESTIMATION_METHOD::HT ) {
      this->fillBlocks ( this->_blocksInitial, this->_preprocessedInput, this->_blockSize, Options );
      this->_blocks = &this->_blocksInitial;
    } else {
      this->fillBlocks ( this->_blocksEstimate, this->_estimate, this->_blockSize, Options );
      this->_blocks = &this->_blocksEstimate;
    }
  }
}


template <typename _RealType, typename _PictureType>
void EMBM3DFilter<_RealType, _PictureType>::aggregate ( const OptionsType &Options, const Blocks3DAndWeight<RealType> &Blocks3D, const aol::MultiVector<short> &Sx, const short Nz ) {
  short x;
  RealType weight;
  int flatIdx;
  for ( short xBlock=0; xBlock<this->_blockSize ; ++xBlock ) {
    for ( short yBlock=0; yBlock<this->_blockSize ; ++yBlock ) {
      weight = Blocks3D.weight ( ) * this->_WWin2D.get ( xBlock, yBlock );
      for ( short z=0; z<Nz ; ++z ) {
        if ( EM_SCANDISTORTIONCORRECTION_METHOD::getShiftMode ( Options.scanDistortionCorrectionMethod ) == EM_SCANDISTORTIONCORRECTION_METHOD::SHIFTMODE_NONE )
          x = Sx[z][0] + xBlock;
        else if ( EM_SCANDISTORTIONCORRECTION_METHOD::getShiftMode ( Options.scanDistortionCorrectionMethod ) == EM_SCANDISTORTIONCORRECTION_METHOD::SHIFTMODE_BLOCKROWSHIFTS ) {
          x = Sx[z][0] + xBlock + round ( this->_blockRowShifts[this->_mapper.getGlobalIndex ( Sx[z][0], Sx[z][1] )][yBlock] );
          x = x % this->_estimate.getNumX ( );
          if ( x < 0 ) x += this->_estimate.getNumX ( );
        } else if ( EM_SCANDISTORTIONCORRECTION_METHOD::getShiftMode ( Options.scanDistortionCorrectionMethod ) == EM_SCANDISTORTIONCORRECTION_METHOD::SHIFTMODE_GLOBALHORIZONTALSHIFTS )
          x = Sx[z][0] + xBlock + round ( this->_globalHorizontalShifts.get ( Sx[z][0] + xBlock, Sx[z][1] + yBlock ) );
        else
          throw aol::Exception ( "Did not recognize scan distortion correction method!", __FILE__, __LINE__ );
        
        flatIdx = this->_mapper.getGlobalIndex ( x, Sx[z][1] + yBlock );
        this->_eBuff[flatIdx] += weight * Blocks3D.get ( xBlock, yBlock, z );
        this->_wBuff[flatIdx] += weight;
        ++this->_numAggregates[flatIdx];
      }
    }
  }
}


template <typename _RealType, typename _PictureType>
void EMBM3DFilter<_RealType, _PictureType>::setEstimateFromBuffers ( const OptionsType &Options ) {
  for ( int k=0; k<this->_estimate.size ( ) ; ++k ) {
    if ( this->_wBuff[k] == 0 ) {
      if ( Options.scanDistortionCorrectionMethod == EM_SCANDISTORTIONCORRECTION_METHOD::JITTERBUG )
        inpaint ( aol::Vec2<short> ( this->_mapper.splitGlobalIndex ( k )[0], this->_mapper.splitGlobalIndex ( k )[1] ) );
      else {
        this->_estimate[k] = 0;
        std::cerr << "WARNING: no weights accumulated for pixel " << this->_mapper.splitGlobalIndex ( k ) << "!" << std::endl;
      }
    } else
      this->_estimate[k] = this->_eBuff[k] / this->_wBuff[k];
  }
}


template <typename _RealType, typename _PictureType>
void EMBM3DFilter<_RealType, _PictureType>::inpaint ( const aol::Vec2<short> &X ) {
  if ( this->_wBuff.get ( X[0]-1, X[1] ) > 0 && this->_wBuff.get ( X[0]+1, X[1] ) > 0 && X[0]-1 >= 0 && X[0]+1 < this->_estimate.getNumX ( ) )
    this->_estimate.set ( X, 0.5 * ( this->_eBuff.get ( X[0]-1, X[1] ) / this->_wBuff.get ( X[0]-1, X[1] ) + this->_eBuff.get ( X[0]+1, X[1] ) / this->_wBuff.get ( X[0]+1, X[1] ) ) );
  else {
    short dx=-1, dy=-1;
    while ( X[0]+dx < 0 || X[0]+dx >= this->_estimate.getNumX ( ) || X[1]+dy < 0 || X[1]+dy >= this->_estimate.getNumY ( ) || this->_wBuff.get ( X[0]+dx, X[1]+dy ) == 0 ) {
      ++dx;
      if ( dx == 1 ) {
        dx = -1;
        ++dy;
      }
    }
    this->_estimate.set ( X, this->_eBuff.get ( X[0]+dx, X[1]+dy ) / this->_wBuff.get ( X[0]+dx, X[1]+dy ) );
  }
}


template <typename _RealType, typename _PictureType>
void EMBM3DFilter<_RealType, _PictureType>::blockMatching ( const OptionsType &Options, const aol::Vec2<short> &XRef,
                                                            aol::MultiVector<short> &Sx, VectorType &BlockDistances, int &Nz ) {
  // Calculate distances to all blocks in the search region
  int z = 0;
  RealType dist;
  for ( this->_searchIt->reset ( XRef ); this->_searchIt->notAtEnd ( ) ; ++(*this->_searchIt) ) {
    (this->*this->_setBlockDistance) ( dist, XRef, *(*this->_searchIt) );
    if ( this->_searchIt->isCurFinal ( ) && dist <= this->_tauMatch ) {
      BlockDistances[z] = dist;
      Sx[z][0] = (*(*this->_searchIt))[0]; Sx[z][1] = (*(*this->_searchIt))[1];
      ++z;
    }
    this->_searchIt->update ( dist );
  }
  Nz = z;
  this->_numMatchedBlocks.set ( XRef, Nz );
  
  // If Nz > _N2 or Nz is not a power of 2 drop coordinates with highest distances until |Sx| = Nz = _N2 or |Sx| = Nz is a power of 2
  // If uniform block-matching is performed, drop all coordinates with highest distances until |Sx| has reached a certain threshold  
  short finalNz;
  if ( Options.similaritySearchMethod == EMBM3D_SIMILARITYSEARCH_METHOD::UNIFORM_PERIODIC )
    finalNz = aol::Min<int> ( aol::Pow ( 2, floor ( log2 ( Nz ) ) ), aol::Max<int> ( this->_N2, this->_searchIt->getNumNonLocalWindows ( ) ) );
  else finalNz = aol::Min<int> ( aol::Pow ( 2, floor ( log2 ( Nz ) ) ), this->_N2 );
  
  if ( Nz > finalNz || Options.similaritySearchMethod == EMBM3D_SIMILARITYSEARCH_METHOD::UNIFORM_PERIODIC ) {
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
  
  // From all remaining blocks choose Nz blocks that satisfy a certain criterion (e.g. uniformity of aggregates)
  if ( Options.similaritySearchMethod == EMBM3D_SIMILARITYSEARCH_METHOD::UNIFORM_PERIODIC ) {
    finalNz = aol::Min<int> ( aol::Pow ( 2, floor ( log2 ( Nz ) ) ), this->_N2 );
    if ( Nz > finalNz ) {
      aol::Vec2<short> tmpPos;
      RealType tmpDist;
      int zMin;
      int min, sum;
      for ( int z=1; z<finalNz ; ++z ) {
        zMin = z;
        min = this->numAggregatesSum ( Sx[zMin][0], Sx[zMin][1] );
        for ( int z2=z+1; z2<Nz && min > 0 ; ++z2 ) {
          sum = this->numAggregatesSum ( Sx[z2][0], Sx[z2][1] );
          if ( sum < min ) {
            zMin = z2;
            min = sum;
          }
        }
        for ( int z2=zMin; z2>z ; --z2 ) {
          tmpPos.set ( Sx[z2][0], Sx[z2][1] );
          tmpDist = BlockDistances[z2];
          Sx[z2][0] = Sx[z2-1][0]; Sx[z2][1] = Sx[z2-1][1];
          BlockDistances[z2] = BlockDistances[z2-1];
          Sx[z2-1][0] = tmpPos[0]; Sx[z2-1][1] = tmpPos[1];
          BlockDistances[z2-1] = tmpDist;
        }
      }
      Nz = finalNz;
    }
  }
  
  for ( int z=0; z<finalNz ; ++z )
    this->_distances[this->_mapper.getGlobalIndex ( XRef[0], XRef[1] )] += BlockDistances[z];
}


template <typename _RealType, typename _PictureType>
void EMBM3DFilter<_RealType, _PictureType>::setSimilaritySearchGrid ( const OptionsType &Options, ColoredPictureType &Dest, const aol::Vec2<short> &XRef ) {
  aol::Vec2<short> xRef ( XRef );
  while ( xRef[0] < 0 || xRef[0] >= this->_XEnd || xRef[1] < 0 || xRef[1] >= this->_YEnd ) {
    const qc::FastILexMapper<qc::QC_2D> mapper ( Options.input.getNumX ( ), Options.input.getNumY ( ) );
    xRef.set ( mapper.splitGlobalIndex ( Options.input.getMaxIndexAndValue ( ).first )[0], mapper.splitGlobalIndex ( Options.input.getMaxIndexAndValue ( ).first )[1] );
  }
  
  BM3DFilter<RealType, PictureType, FilterTrait>::setSimilaritySearchGrid ( Options, Dest, xRef );
}

template class EMBM3DFilter<double, qc::ScalarArray<double, qc::QC_2D> >;

