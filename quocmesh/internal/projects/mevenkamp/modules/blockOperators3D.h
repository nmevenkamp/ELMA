#ifndef BLOCKOPERATORS3D_h
#define BLOCKOPERATORS3D_h

#include "blockOperators2D.h"


template <typename _RealType>
class Blocks3DInitialAndEstimate : public aol::MultiVector<_RealType> {
  typedef _RealType RealType;
  typedef qc::ScalarArray<RealType, qc::QC_3D> BlockStackType;
protected:
  BlockStackType _initial, _estimate;
public:
  Blocks3DInitialAndEstimate ( const int NX, const int NY, const int NZ )
    : _initial ( NX, NY, NZ ), _estimate ( NX, NY, NZ ) {
    appendReferences ( );
  }
  
  Blocks3DInitialAndEstimate ( const BlockStackType &Initial, const BlockStackType &Estimate )
    : _initial ( Initial ), _estimate ( Estimate ) {
    appendReferences ( );
  }
  
  RealType getInitial ( const int X, const int Y, const int Z ) const {
    return _initial.get ( X, Y, Z );
  }
  
  RealType getEstimate ( const int X, const int Y, const int Z ) const {
    return _estimate.get ( X, Y, Z );
  }
  
  BlockStackType& getInitialReference ( ) {
    return _initial;
  }
  
  const BlockStackType& getInitialConstReference ( ) const {
    return _initial;
  }
  
  BlockStackType& getEstimateReference ( ) {
    return _estimate;
  }
  
  const BlockStackType& getEstimateConstReference ( ) const {
    return _estimate;
  }
  
  void setInitial ( const int X, const int Y, const int Z, const RealType Val ) {
    _initial.set ( X, Y, Z, Val );
  }
  
  void setEstimate ( const int X, const int Y, const int Z, const RealType Val ) {
    _estimate.set ( X, Y, Z, Val );
  }
  
  short getNumX ( ) const {
    return _initial.getNumX ( );
  }
  
  short getNumY ( ) const {
    return _initial.getNumY ( );
  }
  
  short getNumZ ( ) const {
    return _initial.getNumZ ( );
  }
private:
  void appendReferences ( ) {
    this->appendReference ( _initial );
    this->appendReference ( _estimate );
  }
};


template <typename _RealType>
class Blocks3DAndWeight : public aol::MultiVector<_RealType> {
  typedef _RealType RealType;
  typedef qc::ScalarArray<RealType, qc::QC_3D> BlockStackType;
protected:
  BlockStackType _blockStack;
  aol::Vector<RealType> _weight;
public:
  Blocks3DAndWeight ( const int NX, const int NY, const int NZ )
  : _blockStack ( NX, NY, NZ ), _weight ( 1 ) {
    appendReferences ( );
  }
  
  Blocks3DAndWeight ( const BlockStackType &BlockStack, const RealType Weight )
  : _blockStack ( BlockStack ), _weight ( 1 ) {
    _weight[0] = Weight;
    appendReferences ( );
  }
  
  BlockStackType& getReference ( ) {
    return _blockStack;
  }
  
  const BlockStackType& getConstReference ( ) const {
    return _blockStack;
  }
  
  RealType get ( const int X, const int Y, const int Z ) const {
    return _blockStack.get ( X, Y , Z );
  }
  
  RealType weight ( ) const {
    return _weight[0];
  }
  
  void set ( const int X, const int Y, const int Z, RealType Val ) {
    _blockStack.set ( X, Y, Z, Val );
  }
  
  void setWeight ( RealType Weight ) {
    _weight[0] = Weight;
  }
  
  short getNumX ( ) const {
    return _blockStack.getNumX ( );
  }
  
  short getNumY ( ) const {
    return _blockStack.getNumY ( );
  }
  
  short getNumZ ( ) const {
    return _blockStack.getNumZ ( );
  }
private:
  void appendReferences ( ) {
    this->appendReference ( _blockStack );
    this->appendReference ( _weight );
  }
};


template <typename _RealType>
class HardThresholding3DOp : aol::Op<qc::ScalarArray<_RealType, qc::QC_3D>, Blocks3DAndWeight<_RealType> > {
  typedef _RealType RealType;
  typedef qc::ScalarArray<RealType, qc::QC_3D> BlockStackType;
protected:
  RealType _threshold;
public:
  HardThresholding3DOp ( ) : _threshold ( 0.0 ) { }
  
  HardThresholding3DOp ( const RealType Threshold ) : _threshold ( Threshold ) { }
  
  void apply ( const BlockStackType &Arg, Blocks3DAndWeight<_RealType> &Dest ) const {
    int Nhar = 0;
    for ( int k=0; k<Arg.size ( ) ; ++k ) {
      if ( aol::Abs<RealType> ( Arg[k] ) > _threshold ) {
        Dest[0][k] = Arg[k];
        ++Nhar;
      } else {
        Dest[0][k] = 0.0;
      }
    }
    Dest[1][0] = ( Nhar >= 1 ) ? 1.0 / Nhar : 1.0;
  }
    
  void applyAdd ( const BlockStackType &/*Arg*/, Blocks3DAndWeight<_RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void setThreshold ( const RealType Threshold ) {
    _threshold = Threshold;
  }
};


template <typename _RealType>
class WienerFiltering3DOp : aol::Op<Blocks3DInitialAndEstimate<_RealType>, Blocks3DAndWeight<_RealType> > {
  typedef _RealType RealType;
protected:
  RealType _stdDev, _stdDevSqr;
public:
  WienerFiltering3DOp ( ) : _stdDev ( 0.0 ), _stdDevSqr ( 0.0 ) { }
  
  WienerFiltering3DOp ( const RealType NoiseStandardDeviation ) : _stdDev ( NoiseStandardDeviation ), _stdDevSqr ( aol::Sqr<RealType> ( NoiseStandardDeviation ) ) { }
  
  void apply ( const Blocks3DInitialAndEstimate<_RealType> &Arg, Blocks3DAndWeight<_RealType> &Dest ) const {
    Dest[1][0] = 0;
    RealType weight, estimateBlockEntrySqr;
    for ( int k=0; k<Arg[0].size ( ) ; ++k ) {
      estimateBlockEntrySqr = aol::Sqr<RealType> ( Arg[1][k] );
      weight = estimateBlockEntrySqr / ( estimateBlockEntrySqr + _stdDevSqr );
      Dest[0][k] = weight * Arg[0][k];
      Dest[1][0] += aol::Sqr<RealType> ( weight );
    }
    Dest[1][0] = 1.0 / Dest[1][0];
  }
    
  void applyAdd ( const Blocks3DInitialAndEstimate<_RealType> &/*Arg*/, Blocks3DAndWeight<_RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void setNoiseStandardDeviation ( const RealType NoiseStandardDeviation ) {
    _stdDev = NoiseStandardDeviation;
    _stdDevSqr = aol::Sqr<RealType> ( NoiseStandardDeviation );
  }
};


template <typename _RealType>
class LinearUnitary2D1DTransformOp : aol::Op<qc::ScalarArray<_RealType, qc::QC_3D> > {
  typedef _RealType RealType;
  typedef aol::FullMatrix<RealType> MatrixType;
  typedef qc::ScalarArray<RealType, qc::QC_2D> BlockType;
  typedef qc::ScalarArray<RealType, qc::QC_3D> BlockStackType;
protected:
  LinearUnitary2DTransformOp<RealType> _transformXY;
  aol::RandomAccessContainer<LinearUnitary1DTransformOp<RealType> > _transformsZ;
public:
  LinearUnitary2D1DTransformOp ( const short BlockSizeXY, const std::string &TransformXY, const short MaxBlockSizeZ, const std::string &TransformZ )
    : _transformXY ( BlockSizeXY, TransformXY ), _transformsZ ( ) {
    setTransformsZ ( MaxBlockSizeZ, TransformZ );
  }
  
  void apply ( const BlockStackType &Arg, BlockStackType &Dest ) const {
    // Apply 2D transform to all blocks (x-y slices along Z-direction)
    BlockType tmpBlockArg ( Arg.getNumX ( ), Arg.getNumY ( ) ), tmpBlockDest ( Arg.getNumX ( ), Arg.getNumY ( ) );
    for ( short z=0; z<Arg.getNumZ ( ) ; ++z ) {
      Arg.getSlice ( qc::QC_Z, z, tmpBlockArg );
      _transformXY.apply ( tmpBlockArg, tmpBlockDest );
      Dest.putSlice ( qc::QC_Z, z, tmpBlockDest );
    }
    
    // Apply 1D transform along Z-axis in each point (x,y)
    aol::Vector<RealType> tmpVecArg ( Arg.getNumZ ( ) ), tmpVecDest ( Arg.getNumZ ( ) );
    for ( short x=0; x<Arg.getNumX ( ) ; ++x )
      for ( short y=0; y<Arg.getNumY ( ) ; ++y ) {
        for ( short z=0; z<Arg.getNumZ ( ) ; ++z ) tmpVecArg[z] = Dest.get ( x, y, z );
        _transformsZ[floor ( log2 ( Arg.getNumZ ( ) ) )].apply ( tmpVecArg, tmpVecDest );
        for ( short z=0; z<Arg.getNumZ ( ) ; ++z ) Dest.set ( x, y, z, tmpVecDest[z] );
      }
  }
    
  void applyInverse ( const BlockStackType &Arg, BlockStackType &Dest ) const {
    // Apply inverse 1D transform along Z-axis in each point (x,y)
    aol::Vector<RealType> tmpVecArg ( Arg.getNumZ ( ) ), tmpVecDest ( Arg.getNumZ ( ) );
    for ( short x=0; x<Arg.getNumX ( ) ; ++x )
      for ( short y=0; y<Arg.getNumY ( ) ; ++y ) {
        for ( short z=0; z<Arg.getNumZ ( ) ; ++z ) tmpVecArg[z] = Arg.get ( x, y, z );
        _transformsZ[floor ( log2 ( Arg.getNumZ ( ) ) )].applyInverse ( tmpVecArg, tmpVecDest );
        for ( short z=0; z<Arg.getNumZ ( ) ; ++z ) Dest.set ( x, y, z, tmpVecDest[z] );
      }
    
    // Apply inverse 2D transform to all blocks (x-y slices along Z-direction)
    BlockType tmpBlockArg ( Arg.getNumX ( ), Arg.getNumY ( ) ), tmpBlockDest ( Arg.getNumX ( ), Arg.getNumY ( ) );
    for ( short z=0; z<Arg.getNumZ ( ) ; ++z ) {
      Dest.getSlice ( qc::QC_Z, z, tmpBlockArg );
      _transformXY.applyInverse ( tmpBlockArg, tmpBlockDest );
      Dest.putSlice ( qc::QC_Z, z, tmpBlockDest );
    }
  }
  
  void applyOnly1DTransform ( const BlockStackType &Arg, BlockStackType &Dest ) const {
    aol::Vector<RealType> tmpVecArg ( Arg.getNumZ ( ) ), tmpVecDest ( Arg.getNumZ ( ) );
    for ( short x=0; x<Arg.getNumX ( ) ; ++x )
      for ( short y=0; y<Arg.getNumY ( ) ; ++y ) {
        for ( short z=0; z<Arg.getNumZ ( ) ; ++z ) tmpVecArg[z] = Arg.get ( x, y, z );
        _transformsZ[floor ( log2 ( Arg.getNumZ ( ) ) )].apply ( tmpVecArg, tmpVecDest );
        for ( short z=0; z<Arg.getNumZ ( ) ; ++z ) Dest.set ( x, y, z, tmpVecDest[z] );
      }
  }
  
  void applyAdd ( const BlockStackType &/*Arg*/, BlockStackType &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void setTransformXY ( const short BlockSizeXY, const std::string &TransformXY ) {
    _transformXY.setTransform ( BlockSizeXY, TransformXY );
  }
  
  void setTransformsZ ( const short MaxBlockSizeZ, const std::string &TransformZ ) {
    _transformsZ.clear ( );
    for ( int i=0; i<=log2 ( MaxBlockSizeZ ) ; ++ i )
      _transformsZ.pushBack ( LinearUnitary1DTransformOp<RealType> ( aol::Pow ( 2, i ), TransformZ ) );
  }
  
  void setTransforms ( const short BlockSizeXY, const std::string &TransformXY, const short MaxBlockSizeZ, const std::string &TransformZ ) {
    setTransformXY ( BlockSizeXY, TransformXY );
    setTransformsZ ( MaxBlockSizeZ, TransformZ );
  }
};


template <typename _RealType>
class BlockStackDenoisingOp : aol::Op<Blocks3DInitialAndEstimate<_RealType>, Blocks3DAndWeight<_RealType> > {
public:
  BlockStackDenoisingOp ( ) { }
  
  virtual ~BlockStackDenoisingOp ( ) { }
  
  virtual void apply ( const Blocks3DInitialAndEstimate<_RealType> &/*Arg*/, Blocks3DAndWeight<_RealType> &/*Dest*/ ) const = 0;
  
  virtual void applyUsingOnly1DTransform ( const Blocks3DInitialAndEstimate<_RealType> &/*Arg*/, Blocks3DAndWeight<_RealType> &/*Dest*/ ) const = 0;
  
  void applyAdd ( const Blocks3DInitialAndEstimate<_RealType> &/*Arg*/, Blocks3DAndWeight<_RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
};


template <typename _RealType>
class HTDenoisingInLocal2D1DTransformDomainOp : public BlockStackDenoisingOp<_RealType> {
  typedef _RealType RealType;
  typedef qc::ScalarArray<RealType, qc::QC_3D> BlockStackType;
protected:
  LinearUnitary2D1DTransformOp<RealType> _2D1DTransformOp;
  HardThresholding3DOp<RealType> _HT3DOp;
public:
  HTDenoisingInLocal2D1DTransformDomainOp ( )
    : _2D1DTransformOp ( 1, "dct", 1, "haar" ), _HT3DOp ( 0.0 ) { }
  
  HTDenoisingInLocal2D1DTransformDomainOp ( const short BlockSizeXY, const std::string &TransformXY, const short MaxBlockSizeZ, const std::string &TransformZ,
                                            const RealType Threshold )
    : _2D1DTransformOp ( BlockSizeXY, TransformXY, MaxBlockSizeZ, TransformZ ), _HT3DOp ( Threshold ) { }
  
  void apply ( const Blocks3DInitialAndEstimate<_RealType> &Arg, Blocks3DAndWeight<_RealType> &Dest ) const {
    // Apply 2D-1D transform on original 3D block stack
    BlockStackType transformedBlockStack ( Arg.getNumX ( ), Arg.getNumY ( ), Arg.getNumZ ( ) );
    _2D1DTransformOp.apply ( Arg.getInitialConstReference ( ), transformedBlockStack );
    
    // Apply hard-thresholding filter on 2D-1D transform of original 3D block stack
    Blocks3DAndWeight<RealType> filteredBlockStackAndWeight ( Arg.getNumX ( ), Arg.getNumY ( ), Arg.getNumZ ( ) );
    _HT3DOp.apply ( transformedBlockStack, filteredBlockStackAndWeight );
    Dest.setWeight ( filteredBlockStackAndWeight.weight ( ) );
    
    // Apply inverse 2D-1D transform on hard-thresholded and transformed original 3D block stack
    _2D1DTransformOp.applyInverse ( filteredBlockStackAndWeight.getConstReference ( ), Dest.getReference ( ) );
  }
  
  void applyUsingOnly1DTransform ( const Blocks3DInitialAndEstimate<_RealType> &Arg, Blocks3DAndWeight<_RealType> &Dest ) const {
    // Apply 2D-1D transform on original 3D block stack
    BlockStackType transformedBlockStack ( Arg.getNumX ( ), Arg.getNumY ( ), Arg.getNumZ ( ) );
    _2D1DTransformOp.applyOnly1DTransform ( Arg.getInitialConstReference ( ), transformedBlockStack );
    
    // Apply hard-thresholding filter on 2D-1D transform of original 3D block stack
    Blocks3DAndWeight<RealType> filteredBlockStackAndWeight ( Arg.getNumX ( ), Arg.getNumY ( ), Arg.getNumZ ( ) );
    _HT3DOp.apply ( transformedBlockStack, filteredBlockStackAndWeight );
    Dest.setWeight ( filteredBlockStackAndWeight.weight ( ) );
    
    // Apply inverse 2D-1D transform on hard-thresholded and transformed original 3D block stack
    _2D1DTransformOp.applyInverse ( filteredBlockStackAndWeight.getConstReference ( ), Dest.getReference ( ) );
  }
  
  void setTransforms ( const short BlockSizeXY, const std::string &TransformXY, const short MaxBlockSizeZ, const std::string &TransformZ ) {
    _2D1DTransformOp.setTransforms ( BlockSizeXY, TransformXY, MaxBlockSizeZ, TransformZ );
  }
  
  void setThreshold ( const RealType Threshold ) {
    _HT3DOp.setThreshold ( Threshold );
  }
};


template <typename _RealType>
class WienerDenoisingInLocal2D1DTransformDomainOp : public BlockStackDenoisingOp<_RealType> {
  typedef _RealType RealType;
  typedef qc::ScalarArray<RealType, qc::QC_3D> BlockStackType;
protected:
  LinearUnitary2D1DTransformOp<RealType> _2D1DTransformOp;
  WienerFiltering3DOp<RealType> _Wiener3DOp;
public:
  WienerDenoisingInLocal2D1DTransformDomainOp ( )
    : _2D1DTransformOp ( 1, "dct", 1, "haar" ), _Wiener3DOp ( 0.0 ) { }
  
  WienerDenoisingInLocal2D1DTransformDomainOp ( const short BlockSizeXY, const std::string &TransformXY, const short MaxBlockSizeZ, const std::string &TransformZ,
                                                const RealType NoiseStandardDeviation )
    : _2D1DTransformOp ( BlockSizeXY, TransformXY, MaxBlockSizeZ, TransformZ ), _Wiener3DOp ( NoiseStandardDeviation ) { }
  
  void apply ( const Blocks3DInitialAndEstimate<_RealType> &Arg, Blocks3DAndWeight<_RealType> &Dest ) const {
    // Apply 2D-1D transforms on original and estimate 3D block stacks
    Blocks3DInitialAndEstimate<RealType> transformedBlockStacks ( Arg.getNumX ( ), Arg.getNumY ( ), Arg.getNumZ ( ) );
    _2D1DTransformOp.apply ( Arg.getInitialConstReference ( ), transformedBlockStacks.getInitialReference ( ) );
    _2D1DTransformOp.apply ( Arg.getEstimateConstReference ( ), transformedBlockStacks.getEstimateReference ( ) );
    
    // Apply Wiener filter on 2D-1D transform of original 3D block stack
    Blocks3DAndWeight<RealType> filteredBlockStackAndWeight ( Arg.getNumX ( ), Arg.getNumY ( ), Arg.getNumZ ( ) );
    _Wiener3DOp.apply ( transformedBlockStacks, filteredBlockStackAndWeight );
    Dest.setWeight ( filteredBlockStackAndWeight.weight ( ) );
    
    // Apply inverse 2D-1D transform on Wiener filtered and transformed original 3D block stack
    _2D1DTransformOp.applyInverse ( filteredBlockStackAndWeight.getConstReference ( ), Dest.getReference ( ) );
  }
  
  void applyUsingOnly1DTransform ( const Blocks3DInitialAndEstimate<_RealType> &Arg, Blocks3DAndWeight<_RealType> &Dest ) const {
    // Apply 2D-1D transforms on original and estimate 3D block stacks
    Blocks3DInitialAndEstimate<RealType> transformedBlockStacks ( Arg.getNumX ( ), Arg.getNumY ( ), Arg.getNumZ ( ) );
    _2D1DTransformOp.applyOnly1DTransform ( Arg.getInitialConstReference ( ), transformedBlockStacks.getInitialReference ( ) );
    _2D1DTransformOp.applyOnly1DTransform ( Arg.getEstimateConstReference ( ), transformedBlockStacks.getEstimateReference ( ) );
    
    // Apply Wiener filter on 2D-1D transform of original 3D block stack
    Blocks3DAndWeight<RealType> filteredBlockStackAndWeight ( Arg.getNumX ( ), Arg.getNumY ( ), Arg.getNumZ ( ) );
    _Wiener3DOp.apply ( transformedBlockStacks, filteredBlockStackAndWeight );
    Dest.setWeight ( filteredBlockStackAndWeight.weight ( ) );
    
    // Apply inverse 2D-1D transform on Wiener filtered and transformed original 3D block stack
    _2D1DTransformOp.applyInverse ( filteredBlockStackAndWeight.getConstReference ( ), Dest.getReference ( ) );
  }
  
  void setTransforms ( const short BlockSizeXY, const std::string &TransformXY, const short MaxBlockSizeZ, const std::string &TransformZ ) {
    _2D1DTransformOp.setTransforms ( BlockSizeXY, TransformXY, MaxBlockSizeZ, TransformZ );
  }
  
  void setNoiseStandardDeviation ( const RealType &NoiseStandardDeviation ) {
    _Wiener3DOp.setNoiseStandardDeviation ( NoiseStandardDeviation );
  }
};


//template <typename _RealType>
//class HTDenoisingInLocal3DFFTDomainOp : public BlockStackDenoisingOp<_RealType> {
//  typedef _RealType RealType;
//  typedef qc::ScalarArray<RealType, qc::QC_3D> BlockStackType;
//protected:
//  RealType _threshold;
//public:
//  HTDenoisingInLocal3DFFTDomainOp ( )
//    : _threshold ( 0.0 ) { }
//  
//  HTDenoisingInLocal3DFFTDomainOp ( const RealType Threshold )
//    : _threshold ( Threshold ) { }
//  
//  void apply ( const Blocks3DInitialAndEstimate<_RealType> &/*Arg*/, Blocks3DAndWeight<_RealType> &/*Dest*/ ) const {
//    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
//    
//    const qc::GridSize<qc::QC_3D> size ( Arg.getNumX ( ), Arg.getNumY ( ), Arg.getNumZ ( ) );
//    qc::MultiArray<RealType, 3, 2> function ( size ), transform ( size );
//    function[0] = Arg[0];
//    function[1].setZero ( );
//    qc::FourierTransform<RealType> ( function, transform, qc::FTForward );
//    
//    int Nhar = 0;
//    for ( qc::RectangularIterator<qc::QC_3D> it ( aol::Vec3<short> ( 0, 0, 0 ), aol::Vec3<short> ( Arg.getNumX ( ), Arg.getNumY ( ), Arg.getNumZ ( ) ) ); it.notAtEnd ( ) ; ++it ) {
//      if ( transform.get ( *it ).norm ( ) > _threshold ) {
//        ++Nhar;
//      } else {
//        transform.set ( *it, aol::Vec2<RealType> ( 0, 0 ) );
//      }
//    }
//    Dest[1][0] = ( Nhar >= 1 ) ? 1.0 / Nhar : 1.0;
//
//    qc::FourierTransform<RealType> ( transform, function, qc::FTBackward );
//    Dest[0] = function[0];
//  }
//  
//  void setThreshold ( const RealType Threshold ) {
//    _threshold = Threshold;
//  }
//};


//template <typename _RealType>
//class WienerDenoisingInLocal3DFFTDomainOp : public BlockStackDenoisingOp<_RealType> {
//  typedef _RealType RealType;
//  typedef qc::ScalarArray<RealType, qc::QC_3D> BlockStackType;
//protected:
//  RealType _noiseStdDev;
//public:
//  WienerDenoisingInLocal3DFFTDomainOp ( )
//    : _noiseStdDev ( 0.0 ) { }
//  
//  WienerDenoisingInLocal3DFFTDomainOp ( const RealType NoiseStandardDeviation )
//  : _noiseStdDev ( NoiseStandardDeviation ) { }
//  
//  void apply ( const Blocks3DInitialAndEstimate<_RealType> &/*Arg*/, Blocks3DAndWeight<_RealType> &/*Dest*/ ) const {
//    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
//    
//    // TODO: 3D FFT on Arg.getInitialConstReference ( ) and Arg.getEstimateConstReference ( )
//    
//    // TODO: apply Wiener Filter to 3D FFT of Arg.getInitialConstReference ( ) based on the values of the 3D FFT of Arg.getEstimateConstReference ( )
//    
//    // TODO: inverse 3D FFT on Wiener filtered 3D FFT of Arg.getInitialConstReference ( )
//  }
//  
//  void setNoiseStandardDeviation ( const RealType &NoiseStandardDeviation ) {
//    _noiseStdDev = NoiseStandardDeviation;
//  }
//};


#endif
