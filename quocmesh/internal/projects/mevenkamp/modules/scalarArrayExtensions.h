#ifndef SCALARARRAYEXTENSIONS_H_
#define SCALARARRAYEXTENSIONS_H_


#include <scalarArray.h>


//template <typename _DataType>
//class ScalarArrayPeriodic : public qc::ScalarArray<_DataType, qc::QC_2D> {
//  typedef _DataType DataType;
//  typedef aol::Vec3<short> CoordType;
//public:
//  ScalarArrayPeriodic ( )
//    : qc::ScalarArray<DataType, qc::QC_2D> ( ) { }
//
//  ScalarArrayPeriodic ( const int Size )
//    : qc::ScalarArray<DataType, qc::QC_2D> ( Size ) { }
//
//  ScalarArrayPeriodic ( const int SizeX, const int SizeY )
//    : qc::ScalarArray<DataType, qc::QC_2D> ( SizeX, SizeY ) { }
//
//  ScalarArrayPeriodic ( const qc::ScalarArray<DataType, qc::QC_2D> &Org )
//    : qc::ScalarArray<DataType, qc::QC_2D> ( Org ) { }
//
//  ScalarArrayPeriodic ( const qc::GridStructure &Grid )
//    : qc::ScalarArray<DataType, qc::QC_2D> ( Grid ) { }
//
//  ScalarArrayPeriodic ( const qc::GridSize<qc::QC_2D> &GridSize )
//    : qc::ScalarArray<DataType, qc::QC_2D> ( GridSize ) { }
//
//  ScalarArrayPeriodic ( const char* Src )
//    : qc::ScalarArray<DataType, qc::QC_2D> ( Src ) { }
//
//  DataType get ( const CoordType &Coord ) const {
//    return getPeriodic ( Coord[0], Coord[1] );
//  }
//
//  DataType get ( const aol::Vec2<short> &Coord ) const {
//    return getPeriodic ( Coord[0], Coord[1] );
//  }
//
//  DataType get ( int X, int Y ) const {
//    return getPeriodic ( X, Y );
//  }
//
//  void add ( const CoordType &Coord, const DataType Val ) {
//    addPeriodic ( Coord[0], Coord[1], Val );
//  }
//
//  void add ( const aol::Vec2<short> &Coord, const DataType Val ) {
//    addPeriodic ( Coord[0], Coord[1], Val );
//  }
//
//  void add ( int X, int Y, const DataType Val ) {
//    addPeriodic ( X, Y, Val );
//  }
//
//private:
//  DataType getPeriodic ( int X, int Y ) const {
//    mapIndicesToDomain ( X, Y );
//    return qc::ScalarArray<DataType, qc::QC_2D>::get ( X, Y );
//  }
//
//  void addPeriodic ( int X, int Y, const DataType Val ) {
//    mapIndicesToDomain ( X, Y );
//    qc::ScalarArray<DataType, qc::QC_2D>::add ( X, Y, Val );
//  }
//
//  void mapIndicesToDomain ( int &X, int &Y ) const {
//    if ( this->getNumX ( ) > 0 )
//      X = X % this->getNumX ( );
//    if ( X < 0 )
//      X += this->getNumX ( );
//    if ( this->getNumY ( ) > 0 )
//      Y = Y % this->getNumY ( );
//    if ( Y < 0 )
//      Y += this->getNumY ( );
//  }
//};


template <typename _RealType, typename _PictureType>
class ArbitraryAnchorBlockCollection {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
  typedef qc::ScalarArray<RealType, qc::QC_2D> BlockType;
protected:
  PictureType _data;
  short _blockSize;
  int _blockSizeSqr;
  aol::Vec2<short> _blockAnchor;
  aol::Vector<RealType> _blocks;
  std::vector<int> _YCornerLookup, _BlockLookup, _YLookup;
public:
  void fillFrom ( const PictureType &Data, const short BlockSize, const aol::Vec2<short> &BlockAnchor ) {
    resize ( Data.getNumX ( ), Data.getNumY ( ), BlockSize, BlockAnchor );
    _data = Data;
    fillBlocks ( );
  }
  
  void fillFrom ( const PictureType &Data, const short BlockSize, const short BlockOffset ) {
    resize ( Data.getNumX ( ), Data.getNumY ( ), BlockSize, BlockOffset );
    _data = Data;
    fillBlocks ( );
  }
  
  void resize ( const short NX, const short NY, const short BlockSize, const aol::Vec2<short> &BlockAnchor ) {
    _data.reallocate ( NX, NY );
    _blockSize = BlockSize;
    _blockSizeSqr = BlockSize * BlockSize;
    _blocks.reallocate ( size ( ) * _blockSizeSqr );
    _blockAnchor = BlockAnchor;
    
    _YCornerLookup.resize ( NY );
    _BlockLookup.resize ( NX * NY );
    _YLookup.resize ( BlockSize );
    initializeLookups ( );
  }
  
  void resize ( const short NX, const short NY, const short BlockSize, const short BlockOffset ) {
    resize ( NX, NY, BlockSize, aol::Vec2<short> ( BlockOffset, BlockOffset ) );
  }
  
  void copyFrom ( const short X, const short Y, const BlockType &Block ) {
    for ( short xBlock=0; xBlock<_blockSize ; ++xBlock )
      for ( short yBlock=0; yBlock<_blockSize ; ++yBlock )
        _blocks[getLinearIndex ( X, Y, xBlock, yBlock )] = Block.get ( xBlock, yBlock );
  }
  
  void copyTo ( const short X, const short Y, BlockType &Block ) {
    for ( short xBlock=0; xBlock<_blockSize ; ++xBlock )
      for ( short yBlock=0; yBlock<_blockSize ; ++yBlock )
        Block.set ( xBlock, yBlock, _blocks[getLinearIndex ( X, Y, xBlock, yBlock )] );
  }
  
  RealType get ( const aol::Vec2<short> &Corner, const short XBlock, const short YBlock ) const {
    return _blocks[getLinearIndex ( Corner[0], Corner[1], XBlock, YBlock )];
  }
  
  RealType get ( const short XCorner, const short YCorner, const short XBlock, const short YBlock ) const {
    return _blocks[getLinearIndex ( XCorner, YCorner, XBlock, YBlock )];
  }
  
  RealType get ( const aol::Vec2<short> &Corner, const short K ) const {
    return _blocks[getLinearIndex ( Corner[0], Corner[1], K )];
  }
  
  RealType get ( const short X, const short Y ) const {
    return _data.get ( X, Y );
  }
  
  RealType get ( const aol::Vec2<short> &Pos ) const {
    return get ( Pos[0], Pos[1] );
  }
  
  void set ( const aol::Vec2<short> &Corner, const short XBlock, const short YBlock, const RealType Val ) {
    _blocks[getLinearIndex ( Corner[0], Corner[1], XBlock, YBlock )] = Val;
  }
  
  void set ( const short XCorner, const short YCorner, const short XBlock, const short YBlock, const RealType Val ) {
    _blocks[getLinearIndex ( XCorner, YCorner, XBlock, YBlock )] = Val;
  }
  
  void set ( const aol::Vec2<short> &Corner, const short K, const RealType Val ) {
    _blocks[getLinearIndex ( Corner[0], Corner[1], K )] = Val;
  }
  
  short getNumX ( ) const {
    return _data.getNumX ( );
  }
  
  short getNumY ( ) const {
    return _data.getNumY ( );
  }
  
  short getNumXEff ( ) const {
    return _data.getNumX ( ) - _blockSize + 1;
  }
  
  short getNumYEff ( ) const {
    return _data.getNumY ( ) - _blockSize + 1;
  }
  
  short getX0 ( ) const {
    return _blockAnchor[0];
  }
  
  short getY0 ( ) const {
    return _blockAnchor[1];
  }
  
  short getXEnd ( ) const {
    return getX0 ( ) + getNumXEff ( );
  }
  
  short getYEnd ( ) const {
    return getY0 ( ) + getNumYEff ( );
  }
  
  int size ( ) const {
    return _data.size ( );
  }
  
  int effSize ( ) const {
    return ( getNumXEff ( ) ) * ( getNumYEff ( ) );
  }
  
  short getBlockSize ( ) const {
    return _blockSize;
  }
  
  const aol::Vec2<short>& getBlockAnchor ( ) const {
    return _blockAnchor;
  }

  short getBlockAnchorX ( ) const {
	  return _blockAnchor[0];
  }

  short getBlockAnchorY ( ) const {
    return _blockAnchor[1];
  }
  
protected:
  inline int getLinearIndex ( const int XCorner, const int YCorner, const int XBlock, const int YBlock ) const {
    return _BlockLookup[_YCornerLookup[YCorner] + XCorner] + _YLookup[YBlock] + XBlock;
  }
  
  inline int getLinearIndex ( const int XCorner, const int YCorner, const int BlockFlatIdx ) const {
    return _BlockLookup[_YCornerLookup[YCorner] + XCorner] + BlockFlatIdx;
  }
  
  void initializeLookups ( ) {
    _YCornerLookup.resize ( getNumY ( ) );
    for ( short yCorner=0; yCorner<getNumY ( ) ; ++yCorner )
      _YCornerLookup[yCorner] = yCorner * getNumX ( );
    
    _BlockLookup.resize ( size ( ) );
    for ( short yCorner=0; yCorner<getNumY ( ) ; ++yCorner )
      for ( short xCorner=0; xCorner<getNumX ( ) ; ++xCorner )
        _BlockLookup[_YCornerLookup[yCorner] + xCorner] = ( _YCornerLookup[yCorner] + xCorner ) * _blockSizeSqr;
    
    _YLookup.resize ( _blockSize );
    for ( short yBlock=0; yBlock<_blockSize ; ++yBlock )
      _YLookup[yBlock] = yBlock * _blockSize;
  }
    
  void fillBlocks ( ) {
    for ( short xCorner=_blockAnchor[0]; xCorner<=getNumX ( )-_blockSize+_blockAnchor[0] ; ++xCorner )
      for ( short yCorner=_blockAnchor[1]; yCorner<=getNumY ( )-_blockSize+_blockAnchor[1] ; ++yCorner )
        for ( short xBlock=0; xBlock<_blockSize ; ++xBlock )
          for ( short yBlock=0; yBlock<_blockSize ; ++yBlock )
            _blocks[getLinearIndex ( xCorner, yCorner, xBlock, yBlock )] = _data.get ( xCorner + xBlock - _blockAnchor[0], yCorner + yBlock - _blockAnchor[1] );
  }
};
  
  
template <typename _RealType, typename _PictureType>
class ArbitraryAnchorShiftedBlockCollection : public ArbitraryAnchorBlockCollection<_RealType, _PictureType> {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
  typedef qc::ScalarArray<RealType, qc::QC_2D> BlockType;
protected:
  aol::MultiVector<RealType> _blockShifts;
  PictureType _globalShifts;
public:
  void fillFrom ( const PictureType &Data, const short BlockSize, const aol::Vec2<short> &BlockAnchor,
                  const aol::MultiVector<RealType> &Shifts ) {
    this->resize ( Data.getNumX ( ), Data.getNumY ( ), BlockSize, BlockAnchor );
    this->_data = Data;
    fillBlocks ( Shifts );
  }
  
  void fillFrom ( const PictureType &Data, const short BlockSize, const short BlockOffset,
                  const aol::MultiVector<RealType> &Shifts ) {
    this->resize ( Data.getNumX ( ), Data.getNumY ( ), BlockSize, aol::Vec2<short> ( BlockOffset, BlockOffset ) );
    this->_data = Data;
    fillBlocks ( Shifts );
  }
  
  void fillFrom ( const PictureType &Data, const short BlockSize, const aol::Vec2<short> &BlockAnchor,
                  const PictureType &Shifts ) {
    this->resize ( Data.getNumX ( ), Data.getNumY ( ), BlockSize, BlockAnchor );
    this->_data = Data;
    fillBlocks ( Shifts );
  }
  
  void fillFrom ( const PictureType &Data, const short BlockSize, const short BlockOffset,
                  const PictureType &Shifts ) {
    this->resize ( Data.getNumX ( ), Data.getNumY ( ), BlockSize, aol::Vec2<short> ( BlockOffset, BlockOffset ) );
    this->_data = Data;
    fillBlocks ( Shifts );
  }
protected:
  void fillBlocks ( const aol::MultiVector<RealType> &Shifts ) {
    int flatIdx;
    qc::FastILexMapper<qc::QC_2D> mapper ( this->getNumX ( ), this->getNumY ( ) );
    for ( short xCorner=this->_blockAnchor[0]; xCorner<=this->getNumX ( )-this->_blockSize+this->_blockAnchor[0] ; ++xCorner ) {
      for ( short yCorner=this->_blockAnchor[1]; yCorner<=this->getNumY ( )-this->_blockSize+this->_blockAnchor[1] ; ++yCorner ) {
        flatIdx = mapper.getGlobalIndex ( xCorner, yCorner );
        for ( short xBlock=0; xBlock<this->_blockSize ; ++xBlock )
          for ( short yBlock=0; yBlock<this->_blockSize ; ++yBlock )
            this->_blocks[this->getLinearIndex ( xCorner, yCorner, xBlock, yBlock )]
              = this->_data.get ( xCorner + xBlock - this->_blockAnchor[0] + round ( Shifts[flatIdx][yBlock] ), yCorner + yBlock - this->_blockAnchor[1] );
      }
    }
  }
  
  void fillBlocks ( const PictureType &Shifts ) {
    for ( short xCorner=this->_blockAnchor[0]; xCorner<=this->getNumX ( )-this->_blockSize+this->_blockAnchor[0] ; ++xCorner )
      for ( short yCorner=this->_blockAnchor[1]; yCorner<=this->getNumY ( )-this->_blockSize+this->_blockAnchor[1] ; ++yCorner )
        for ( short xBlock=0; xBlock<this->_blockSize ; ++xBlock )
          for ( short yBlock=0; yBlock<this->_blockSize ; ++yBlock )
            this->_blocks[this->getLinearIndex ( xBlock, yBlock, xCorner, yCorner )]
              = this->_data.get ( xCorner + xBlock - this->_blockAnchor[0] + round ( Shifts.get ( xCorner + xBlock - this->_blockAnchor[0], yCorner + yBlock - this->_blockAnchor[1] ) ),
                                  yCorner + yBlock - this->_blockAnchor[1] );
  }
};
  

template <typename _RealType, typename _PictureType>
class BlockCollection : public ArbitraryAnchorBlockCollection<_RealType, _PictureType> {
public:
  void fillFrom ( const _PictureType &Data, const short BlockSize ) {
    ArbitraryAnchorBlockCollection<_RealType, _PictureType>::fillFrom ( Data, BlockSize, 0 );
  }
  
  void resize ( const short NX, const short NY, const short BlockSize ) {
    ArbitraryAnchorBlockCollection<_RealType, _PictureType>::resize ( NX, NY, BlockSize, 0 );
  }
};
  
  
template <typename _RealType, typename _PictureType>
class ShiftedBlockCollection : public ArbitraryAnchorShiftedBlockCollection<_RealType, _PictureType> {
public:
  void fillFrom ( const _PictureType &Data, const short BlockSize ) {
    ArbitraryAnchorBlockCollection<_RealType, _PictureType>::fillFrom ( Data, BlockSize, 0 );
  }
  
  void fillFrom ( const _PictureType &Data, const short BlockSize,
                  const aol::MultiVector<_RealType> &Shifts ) {
    ArbitraryAnchorShiftedBlockCollection<_RealType, _PictureType>::fillFrom ( Data, BlockSize, 0, Shifts );
  }
  
  void fillFrom ( const _PictureType &Data, const short BlockSize,
                  const _PictureType &Shifts ) {
    ArbitraryAnchorShiftedBlockCollection<_RealType, _PictureType>::fillFrom ( Data, BlockSize, 0, Shifts );
  }
  
  void resize ( const short NX, const short NY, const short BlockSize ) {
    ArbitraryAnchorShiftedBlockCollection<_RealType, _PictureType>::resize ( NX, NY, BlockSize, 0 );
  }
};
  
  
template <typename _RealType, typename _PictureType>
class CenteredBlockCollection : public ArbitraryAnchorBlockCollection<_RealType, _PictureType> {
public:
  void fillFrom ( const _PictureType &Data, const short BlockSize ) {
    if ( BlockSize % 2 == 0 )
      throw aol::Exception ( "Block size must be odd!", __FILE__, __LINE__ );

    ArbitraryAnchorBlockCollection<_RealType, _PictureType>::fillFrom ( Data, BlockSize, ( BlockSize - 1 ) / 2 );
  }
  
  void resize ( const short NX, const short NY, const short BlockSize ) {
    ArbitraryAnchorBlockCollection<_RealType, _PictureType>::resize ( NX, NY, BlockSize, ( BlockSize - 1 ) / 2 );
  }
};
  
  
template <typename _RealType, typename _PictureType>
class CenteredShiftedBlockCollection : public ArbitraryAnchorShiftedBlockCollection<_RealType, _PictureType> {
public:
  void fillFrom ( const _PictureType &Data, const short BlockSize ) {
    if ( BlockSize % 2 == 0 )
      throw aol::Exception ( "Block size must be odd!", __FILE__, __LINE__ );
    
    ArbitraryAnchorBlockCollection<_RealType, _PictureType>::fillFrom ( Data, BlockSize, ( BlockSize - 1 ) / 2 );
  }
  
  void fillFrom ( const _PictureType &Data, const short BlockSize,
                  const aol::MultiVector<_RealType> &Shifts ) {
    if ( BlockSize % 2 == 0 )
      throw aol::Exception ( "Block size must be odd!", __FILE__, __LINE__ );
    
    ArbitraryAnchorShiftedBlockCollection<_RealType, _PictureType>::fillFrom ( Data, BlockSize, ( BlockSize - 1 ) / 2, Shifts );
  }
  
  void fillFrom ( const _PictureType &Data, const short BlockSize,
                  const _PictureType &Shifts ) {
    if ( BlockSize % 2 == 0 )
      throw aol::Exception ( "Block size must be odd!", __FILE__, __LINE__ );
    
    ArbitraryAnchorShiftedBlockCollection<_RealType, _PictureType>::fillFrom ( Data, BlockSize, ( BlockSize - 1 ) / 2, Shifts );
  }
  
  void resize ( const short NX, const short NY, const short BlockSize ) {
    ArbitraryAnchorShiftedBlockCollection<_RealType, _PictureType>::resize ( NX, NY, BlockSize, ( BlockSize - 1 ) / 2 );
  }
};


#endif /* SCALARARRAYEXTENSIONS_H_ */
