#include "emNonLocalMeansFilter.h"


template <typename _RealType, typename _PictureType>
void EMNonLocalMeansFilter<_RealType, _PictureType>::initializeBlocks ( const OptionsType &Options ) {
  this->setScanDistortionCorrectionParameters ( Options );
  this->fillBlocks ( this->_blocksInitial, this->_preprocessedInput, this->_blockSize, Options );
  this->_blocks = &this->_blockInitial;
}