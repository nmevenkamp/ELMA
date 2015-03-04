#include <scalarArray.h>
#include <multiArray.h>
#include <pngInterface.h>

namespace qc {


template <typename DataType, int imagedim>
#ifdef USE_LIB_PNG
void loadPNGToMArray ( qc::MultiArray<DataType, 2, imagedim> &MArg, const char *FileName, const bool StripAlpha ) {
  PNGLoadInterface pngInterface( FileName, false, StripAlpha );
  MArg.reallocate ( pngInterface.getWidth(), pngInterface.getHeight() );
  pngInterface.writeDataToMultiArray<DataType>( MArg );
#else
void loadPNGToMArray ( qc::MultiArray<DataType, 2, imagedim> &, const char *, const bool ) {
  throw aol::Exception ( "Reading PNG files requires libpng. Compile with -DUSE_LIB_PNG.", __FILE__, __LINE__ );
#endif //USE_LIB_PNG
}

template <typename DataType, int imagedim>
#ifdef USE_LIB_PNG
void savePNGFromMArray ( const qc::MultiArray<DataType, 2, imagedim> &MArg, const char *FileName ) {
  PNGSaveInterface pngInterface( FileName );
  pngInterface.loadDataFromMultiArray<DataType>( MArg );
#else
void savePNGFromMArray ( const qc::MultiArray<DataType, 2, imagedim> &, const char * ) {
  throw aol::Exception ( "Writing PNG files requires libpng. Compile with -DUSE_LIB_PNG.", __FILE__, __LINE__ );
#endif //USE_LIB_PNG
}

/**
 * This struct + static function construction is a workaround to C++ limitations.
 *
 * \author Berkels
 */
template <typename DataType, int rangedim, int imagedim>
struct doPNGLoadSave {};

template <typename DataType>
struct doPNGLoadSave<DataType, 2, 2> {
  static void loadPNG ( qc::MultiArray<DataType, 2, 2> &MArg, const char *FileName ) {
    loadPNGToMArray ( MArg, FileName, false );
  }

  static void savePNG ( const qc::MultiArray<DataType, 2, 2> &/*MArg*/, const char * /*FileName*/ ) {
    throw aol::UnimplementedCodeException( "doPNGLoadSave<DataType, 2, 2> PNGSaveInterface::savePNG is not implemented.", __FILE__, __LINE__);
  }
};

template <typename DataType>
struct doPNGLoadSave<DataType, 2, 3> {
  static void loadPNG ( qc::MultiArray<DataType, 2, 3> &MArg, const char *FileName ) {
    loadPNGToMArray ( MArg, FileName, true );
  }

  static void savePNG ( const qc::MultiArray<DataType, 2, 3> &MArg, const char *FileName ) {
    savePNGFromMArray ( MArg, FileName );
  }
};

template <typename DataType>
struct doPNGLoadSave<DataType, 2, 4> {
  static void loadPNG ( qc::MultiArray<DataType, 2, 4> &MArg, const char *FileName ) {
    loadPNGToMArray ( MArg, FileName, false );
  }

  static void savePNG ( const qc::MultiArray<DataType, 2, 4> &MArg, const char *FileName ) {
    savePNGFromMArray ( MArg, FileName );
  }
};

template <typename DataType>
struct doPNGLoadSave<DataType, 2, 5> {
  static void loadPNG ( qc::MultiArray<DataType, 2, 5> &MArg, const char *FileName ) {
    loadPNGToMArray ( MArg, FileName, false );
  }

  static void savePNG ( const qc::MultiArray<DataType, 2, 5> &MArg, const char *FileName ) {
    savePNGFromMArray ( MArg, FileName );
  }
};

template <typename DataType>
struct doPNGLoadSave<DataType, 2, 1> {
  static void loadPNG ( qc::MultiArray<DataType, 2, 1> &MArg, const char *FileName ) {
    MArg[0].loadPNG ( FileName );
  }

  static void savePNG ( const qc::MultiArray<DataType, 2, 1> &MArg, const char *FileName ) {
    MArg[0].savePNG ( FileName );
  }
};

template <class DataType, int rangedim, int imagedim>
void qc::MultiArray<DataType, rangedim, imagedim>::loadPNG ( const char *FileName ) {
  doPNGLoadSave<DataType, rangedim, imagedim>::loadPNG ( *this, FileName );
}

template <class DataType, int rangedim, int imagedim>
void qc::MultiArray<DataType, rangedim, imagedim>::savePNG ( const char *FileName ) const {
  doPNGLoadSave<DataType, rangedim, imagedim>::savePNG ( *this, FileName );
}


} // namespace qc

template class qc::MultiArray<double, 2, 1>;
template class qc::MultiArray<double, 2, 2>;
template class qc::MultiArray<float, 2, 3>;
template class qc::MultiArray<double, 2, 3>;
template class qc::MultiArray<long double, 2, 3>;
template class qc::MultiArray<unsigned char, 2, 3>;
template class qc::MultiArray<float, 2, 4>;
template class qc::MultiArray<double, 2, 4>;
template class qc::MultiArray<long double, 2, 4>;

template class qc::MultiArray<double, 2, 5>;
