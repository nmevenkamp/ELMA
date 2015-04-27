#include <convolution.h>

namespace qc {
  
#ifdef USE_LIB_FFTW

//! 1D Fourier transform (complex-to-complex)
//! Unnecessarily complicated because ScalarArray<complex> is not possible
//! (otherwise copying would not be necessary)
template <typename RealType>
void FourierTransform ( const aol::MultiVector<RealType>& function, aol::MultiVector<RealType>& transform, enum FourierTransformDirection direction ) {
  typedef ConvolutionTrait<qc::QC_1D,RealType> ConvType;
  typedef typename ConvType::FFTWComplex FFTWComplex;
  
  
  // Prepare transformation
  if ( function.numComponents ( ) != 2 || transform.numComponents ( ) != 2 )
    throw aol::Exception ( "MultiVector components not equal to two in FourierTransform", __FILE__, __LINE__ );
  int size = function [0].size ();
  if ( function [1].size () != size || transform [0].size () != size ||  transform [1].size () != size )
    throw aol::Exception ( "Array sizes not equal in FourierTransform", __FILE__, __LINE__ );
  FFTWComplex* f = static_cast<FFTWComplex*> ( ConvType::fftwMalloc ( sizeof ( FFTWComplex ) * size ) );
  FFTWComplex* t = static_cast<FFTWComplex*> ( ConvType::fftwMalloc ( sizeof ( FFTWComplex ) * size ) );
  typename ConvType::FFTWPlan plan = ConvType::fftwPlan_dft ( size, f, t, direction, FFTW_ESTIMATE );
  
  // Copy data, transform and copy back
  for ( int i = 0; i < size; ++i )
    for ( int k = 0; k < 2; ++k ) f [i] [k] = function [k] [i];
  ConvType::fftwExecute ( plan );
  for ( int i = 0; i < size; ++i )
    for ( int k = 0; k < 2; ++k ) transform [k] [i] = t [i] [k];
  
  // Cleanup
  ConvType::fftwDestroy_plan ( plan );
  ConvType::fftwFree ( f );
  ConvType::fftwFree ( t );
}
  
//! 2D Fourier transform (complex-to-complex)./
//! Unnecessarily complicated because ScalarArray<complex> is not possible
//! (otherwise copying would not be necessary)
template <typename RealType>
void FourierTransform ( const qc::MultiArray<RealType, 2, 2>& function, qc::MultiArray<RealType, 2, 2>& transform, enum FourierTransformDirection direction ) {
  typedef ConvolutionTrait<qc::QC_2D,RealType> ConvType;
  typedef typename ConvType::FFTWComplex FFTWComplex;


  // Prepare transformation
  int numX = function [0].getNumX (), numY = function [0].getNumY ();
  if ( function [1].getNumX () != numX || transform [0].getNumX () != numX ||  transform [1].getNumX () != numX ||
       function [1].getNumY () != numY || transform [0].getNumY () != numY ||  transform [1].getNumY () != numY )
    throw aol::Exception ( "Array sizes not equal in FourierTransform", __FILE__, __LINE__ );
  FFTWComplex* f = static_cast<FFTWComplex*> ( ConvType::fftwMalloc ( sizeof ( FFTWComplex ) * numX * numY ) );
  FFTWComplex* t = static_cast<FFTWComplex*> ( ConvType::fftwMalloc ( sizeof ( FFTWComplex ) * numX * numY ) );
  typename ConvType::FFTWPlan plan = ConvType::fftwPlan_dft ( aol::Vec2<int> ( numX, numY ), f, t, direction, FFTW_ESTIMATE );

  // Copy data, transform and copy back
  for ( int j = 0; j < numY; ++j ) {
    for ( int i = 0; i < numX; ++i ) {
      const int ind = qc::ILexCombine2 ( i, j, numX );
      const aol::Vec2<short int> pos ( i, j );
      aol::Vec2<RealType> val = function.get ( pos );
      for ( int k = 0; k < 2; ++k ) f [ind] [k] = val [k];
    }
  }
  ConvType::fftwExecute ( plan );
  for ( int j = 0; j < numY; ++j ) {
    for ( int i = 0; i < numX; ++i ) {
      const int ind = qc::ILexCombine2 ( i, j, numX );
      const aol::Vec2<short int> pos ( i, j );
      aol::Vec2<RealType> val ( t [ind] [0], t [ind] [1] );
      transform.set ( pos, val );
    }
  }

  // Cleanup
  ConvType::fftwDestroy_plan ( plan );
  ConvType::fftwFree ( f );
  ConvType::fftwFree ( t );
}
  
#elif defined ( USE_KISSFFT )

//! 1D Fourier transform (complex-to-complex)
//! Unnecessarily complicated because ScalarArray<complex> is not possible
//! (otherwise copying would not be necessary)
template<>
void FourierTransform<double> ( const aol::MultiVector<double>& function, aol::MultiVector<double>& transform, enum FourierTransformDirection direction ) {
  typedef ConvolutionTrait<qc::QC_1D,double> ConvType;
  typedef ConvType::KISSFFTComplex KISSFFTComplex;
  
  
  // Prepare transformation
  if ( function.numComponents ( ) != 2 || transform.numComponents ( ) != 2 )
    throw aol::Exception ( "MultiVector components not equal to two in FourierTransform", __FILE__, __LINE__ );
  int size = function [0].size ();
  if ( function [1].size () != size || transform [0].size () != size ||  transform [1].size () != size )
    throw aol::Exception ( "Array sizes not equal in FourierTransform", __FILE__, __LINE__ );
  KISSFFTComplex* f = static_cast<KISSFFTComplex*> ( ConvType::kissfftMalloc ( sizeof ( KISSFFTComplex ) * size ) );
  KISSFFTComplex* t = static_cast<KISSFFTComplex*> ( ConvType::kissfftMalloc ( sizeof ( KISSFFTComplex ) * size ) );
  ConvType::KISSFFTPlan plan = ConvType::kissfftPlan_dft ( size, direction );
  
  // Copy data, transform and copy back
  for ( int i = 0; i < size; ++i ) {
    f [i].r = function [0] [i];
    f [i].i = function [1] [i];
  }
  ConvType::kissfftExecute ( plan, f, t );
  for ( int i = 0; i < size; ++i ) {
    transform [0] [i] = t [i].r;
    transform [1] [i] = t [i].i;
  }
  
  // Cleanup
  ConvType::kissfftDestroy_plan ( plan );
  ConvType::kissfftFree ( f );
  ConvType::kissfftFree ( t );
}
  
//! 2D Fourier transform (complex-to-complex)
//! Unnecessarily complicated because ScalarArray<complex> is not possible
//! (otherwise copying would not be necessary)
template<>
void FourierTransform<double> ( const qc::MultiArray<double, 2, 2>& function, qc::MultiArray<double, 2, 2>& transform, enum FourierTransformDirection direction ) {
  typedef ConvolutionTrait<qc::QC_2D,double> ConvType;
  typedef ConvType::KISSFFTComplex KISSFFTComplex;
  
  
  // Prepare transformation
  int numX = function [0].getNumX (), numY = function [0].getNumY ();
  if ( function [1].getNumX () != numX || transform [0].getNumX () != numX ||  transform [1].getNumX () != numX ||
      function [1].getNumY () != numY || transform [0].getNumY () != numY ||  transform [1].getNumY () != numY )
    throw aol::Exception ( "Array sizes not equal in FourierTransform", __FILE__, __LINE__ );
  KISSFFTComplex* f = static_cast<KISSFFTComplex*> ( ConvType::kissfftMalloc ( sizeof ( KISSFFTComplex ) * numX * numY ) );
  KISSFFTComplex* t = static_cast<KISSFFTComplex*> ( ConvType::kissfftMalloc ( sizeof ( KISSFFTComplex ) * numX * numY ) );
  ConvType::KISSFFTPlan plan = ConvType::kissfftPlan_dft ( aol::Vec2<int> ( numX, numY ), direction );
  
  // Copy data, transform and copy back
  for ( int j = 0; j < numY; ++j ) {
    for ( int i = 0; i < numX; ++i ) {
      const int ind = qc::ILexCombine2 ( i, j, numX );
      const aol::Vec2<short int> pos ( i, j );
      aol::Vec2<double> val = function.get ( pos );
      f [ind].r = val [0];
      f [ind].i = val [1];
    }
  }
  ConvType::kissfftExecute ( plan, f, t );
  for ( int j = 0; j < numY; ++j ) {
    for ( int i = 0; i < numX; ++i ) {
      const int ind = qc::ILexCombine2 ( i, j, numX );
      const aol::Vec2<short int> pos ( i, j );
      aol::Vec2<double> val ( t [ind].r, t [ind].i );
      transform.set ( pos, val );
    }
  }
  
  // Cleanup
  ConvType::kissfftDestroy_plan ( plan );
  ConvType::kissfftFree ( f );
  ConvType::kissfftFree ( t );
}

template<>
void FourierTransform<float> ( const aol::MultiVector<float>& /*function*/, aol::MultiVector<float>& /*transform*/, enum FourierTransformDirection /*direction*/ )
{
  throw aol::Exception ( "libffts does not support float type! Compile with -DUSE_LIB_FFTW", __FILE__, __LINE__ );
}
  
template<>
void FourierTransform<float> ( const qc::MultiArray<float, 2, 2>& /*function*/, qc::MultiArray<float, 2, 2>& /*transform*/, enum FourierTransformDirection /*direction*/ )
{
  throw aol::Exception ( "libffts does not support float type! Compile with -DUSE_LIB_FFTW", __FILE__, __LINE__ );
}
  
#else

//! 1D Fourier transform (complex-to-complex)
template <typename RealType>
void FourierTransform ( const aol::MultiVector<RealType>& /*function*/, aol::MultiVector<RealType>& /*transform*/, enum FourierTransformDirection /*direction*/ )
// this (unnamed parameter but default value) seems to work with gcc-4.2.0 and VC++ 2005. if it does not work with other compilers, write separate function with two parameters.
{
  throw aol::Exception ( "FourierTransform needs libfftw! Compile with -DUSE_LIB_FFTW", __FILE__, __LINE__ );
}
  
//! 2D Fourier transform (complex-to-complex)
template <typename RealType>
void FourierTransform ( const qc::MultiArray<RealType, 2, 2>& /*function*/, qc::MultiArray<RealType, 2, 2>& /*transform*/, enum FourierTransformDirection /*direction*/ )
// this (unnamed parameter but default value) seems to work with gcc-4.2.0 and VC++ 2005. if it does not work with other compilers, write separate function with two parameters.
{
  throw aol::Exception ( "FourierTransform needs libfftw! Compile with -DUSE_LIB_FFTW", __FILE__, __LINE__ );
}

#endif

template void FourierTransform<float> ( const aol::MultiVector<float>&, aol::MultiVector<float>&, enum FourierTransformDirection );
template void FourierTransform<double> ( const aol::MultiVector<double>&, aol::MultiVector<double>&, enum FourierTransformDirection );
template void FourierTransform<float> ( const qc::MultiArray<float, 2, 2>&, qc::MultiArray<float, 2, 2>&, enum FourierTransformDirection );
template void FourierTransform<double> ( const qc::MultiArray<double, 2, 2>&, qc::MultiArray<double, 2, 2>&, enum FourierTransformDirection );
// Since fftw3l seems to be missing on our Linux machines, we can't use the long double version.
//template void FourierTransform<long double> ( const qc::MultiArray<long double, 2, 2>&, qc::MultiArray<long double, 2, 2>&, enum FourierTransformDirection );

  
void addMotionBlurToArray ( const aol::Vec2<double> &Velocity, const qc::ScalarArray<double, qc::QC_2D> &Arg, qc::ScalarArray<double, qc::QC_2D> &Dest ) {
  qc::Convolution<qc::QC_2D> conv ( aol::Vec2<int>( Arg.getNumX(), Arg.getNumY() ) );
  qc::ScalarArray<double, qc::QC_2D> kernel ( Arg, aol::STRUCT_COPY );
  qc::generateMotionBlurKernel<double> ( Velocity, kernel );
  conv.convolve ( Arg, kernel, Dest );
}


} // end namespace
