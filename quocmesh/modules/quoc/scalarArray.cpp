#include <scalarArray.h>
#include <aol.h>

#include <pngInterface.h>
#include <bzipiostream.h>

#include <quoc.h>

#include <auxiliary.h>
#include <qmException.h>
#include <indexMapper.h>
#include <imageTools.h>
#include <dm3Import.h>

#ifdef USE_LIB_TIFF
#include <tiffio.h>
#endif

#ifdef USE_EXTERNAL_CIMG
#include <cimgIncludes.h>
#endif

char FILE_TYPE[][14] = { "", "",
                         "ASCII",
                         "", "",
                         "BINARY",
                         "",
                         "ASCII FLOAT",
                         "RAW FLOAT",
                         "RAW DOUBLE"
                       };

template < typename DataType > const aol::Format & qc::ScalarArray < DataType, qc::QC_2D > ::format = aol::mixedFormat;
template < typename DataType > bool qc::ScalarArray < DataType, qc::QC_2D > ::prettyFormat = true;

template <typename _DataType>
qc::ScalarArray<_DataType, qc::QC_1D>::ScalarArray ( const string &filename ) : Array<_DataType> ( ) {
  aol::Bzipifstream file ( filename.c_str () );
  ArrayHeader header;
  ReadArrayHeader ( file, header );
  reallocate ( header.numX );
  load ( filename.c_str () );
  // apparently no init() call in 1D
}


template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_1D>::load ( istream &in ) {
  char          tmp[256];
  int           type, inWidth, inMaxColor;

  if ( ! ( this->quietMode ) ) {
    cerr << "qc::ScalarArray<DataType, qc::QC_1D>::load( istream &in )\n";
  }

  in.get ( tmp, 255, '\n' );
  if ( tmp[0] != 'O' )
    throw aol::TypeException ( "qc::ScalarArray<DataType, qc::QC_1D>::load: wrong file format",
                               __FILE__, __LINE__ );

  switch ( tmp[1] ) {
    case '2':
      type = PGM_UNSIGNED_CHAR_ASCII;
      break;
    case '5':
      type = PGM_UNSIGNED_CHAR_BINARY;
      break;
    case '7':
      type = PGM_FLOAT_ASCII;
      break;
    case '8':
      type = PGM_FLOAT_BINARY;
      break;
    case '9':
      type = PGM_DOUBLE_BINARY;
      break;

    default:
      throw aol::TypeException ( "qc::ScalarArray<DataType, qc::QC_1D>::load: wrong magic number",
                                 __FILE__, __LINE__ );
  }

  aol::READ_COMMENTS ( in );
  in >> inWidth;
  aol::READ_COMMENTS ( in );
  in >> inMaxColor;
  in.ignore();

#ifdef VERBOSE
  cerr << "Read from header: " << inWidth << " " << inMaxColor << endl;
#endif

  this->loadRaw ( in, type, inWidth );
}

template < typename _DataType >
void qc::ScalarArray<_DataType, qc::QC_1D>::loadRaw ( istream &in, const int Type, const int InWidth ) {
  if ( InWidth != this->numX ) {
    char errorMsg[1024];
    sprintf ( errorMsg,
              "qc::ScalarArray<DataType, qc::QC_1D>::loadRaw: dimension does not match\n"
              "loading size (%d) but expecting (%d)", InWidth, this->numX );
    throw aol::TypeException ( errorMsg, __FILE__, __LINE__ );
  }

  this->aol::Vector<_DataType>::loadRaw ( in, Type );

  if ( in.good() ) {
    if ( !this->quietMode ) cerr << "Success. 1D Stream (" << InWidth << ") was of type " << Type << ". " << endl;
  } else throw aol::IOException ( "qc::ScalarArray<DataType, qc::QC_1D>::loadRaw: Error reading from file", __FILE__, __LINE__ );
}


template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_1D>::load ( const char *fileName ) {
  if ( !this->quietMode ) {
    cerr << "loading from file " << fileName << endl;
  };

  if ( aol::fileNameEndsWith ( fileName, ".png" ) ) {
    throw aol::Exception ( "qc::ScalarArray<DataType, qc::QC_1D>::loadPNG not yet implemented", __FILE__, __LINE__ );
    //loadPNG ( fileName );
  } else {

    // Decompression is done by aol::Bzipifstream based on filename extension
    aol::Bzipifstream in ( fileName );
    load ( in );
  }
  if ( !this->quietMode ) cerr << "done." << endl;
}

template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_1D>::save ( ostream &out, qc::SaveType type, const char *comment ) const {

  if ( type == PNG_2D )
    throw aol::TypeException ( "qc::ScalarArray<DataType, qc::QC_1D>::save( ostream&, qc::SaveType, const char* ): impossible with type PNG_2D", __FILE__, __LINE__ );

  DataType  maximum = this->getMaxValue(), minimum = this->getMinValue();

  char      info[4096];

  if ( !comment ) {
    sprintf ( info, "# This is a QuOcMesh file of type %d (=%s) written %s", type, FILE_TYPE[type],
              aol::generateCurrentTimeAndDateString().c_str() );
    comment = info;
  }

  // make sure comment starts with '#'
  if ( comment[0] != '#' ) {
    info[0] = '#';
    strncpy ( info + 1, comment, 4095 );
    info[strlen ( comment ) ] = 0;
    comment = info;
  }

  out << "O" << type << endl <<  comment << endl << this->numX << endl;
  switch ( type ) {
    case PGM_UNSIGNED_CHAR_ASCII:
      out << "255" << endl;
      break;
    case PGM_UNSIGNED_CHAR_BINARY:
      out << "255\n";
      break;
    case PGM_FLOAT_ASCII:
      out << static_cast<int> ( maximum ) << "\n";
      break;
    case PGM_FLOAT_BINARY:
      out <<  static_cast<int> ( maximum ) << "\n";
      break;
    case PGM_DOUBLE_BINARY:
      out << static_cast<int> ( maximum ) << "\n";
      break;
    default:
      throw aol::TypeException ( "qc::ScalarArray<DataType, qc::QC_1D>::save: illegal type", __FILE__, __LINE__ );
  }

  this->saveRaw ( out, type, minimum, maximum );
}

template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_1D>::save ( const char *fileName,
                                                   qc::SaveType type,
                                                   const char *comment  ) const {

//   if ( type == PNG_2D )
//     if ( comment != NULL )
//       throw aol::TypeException ( "qc::ScalarArray<DataType, qc::QC_2D>::save( const char*, qc::SaveType, const char* ) with type PNG allows no comment", __FILE__, __LINE__ );
//     else
//       savePNG( fileName );
//   else {
  if ( type != 2 && type != 5 && ! ( type >= 7 && type <= 9 ) )
    throw aol::TypeException ( "qc::ScalarArray<DataType, qc::QC_1D>::save: Unknown magic number",
                               __FILE__, __LINE__ );

  aol::Bzipofstream *out = new aol::Bzipofstream ( fileName );

  if ( !out->good() ) {
    delete out;
    throw aol::FileException ( "qc::ScalarArray<DataType, qc::QC_1D>::save: Cannot open file for writing",
                               __FILE__, __LINE__ );
  }

  save ( *out, type, comment );
  delete out;
  if ( !this->quietMode ) cerr << "done.\n";
//  }
}


template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_1D>::saveToFile ( const char *filename ) const {
  aol::Bzipofstream out ( filename );
  out << aol::VectorFileMagicChar::ScalarArray_1D << aol::FileFormatMagicNumber<DataType>::FFType << endl;
  out << "# This is a QuOcMesh file of type " << aol::VectorFileMagicChar::ScalarArray_1D << aol::FileFormatMagicNumber<DataType>::FFType
      << " storing a qc::ScalarArray<" << aol::FileFormatMagicNumber<DataType>::FFContainedName << ", QC_1D>" << endl;
  out << this->getNumX() << endl;
  out << this->getMaxValue() << endl; // for compatibility with PGM image format in some cases, else can be ignored
  const char* buffer = reinterpret_cast<char*> ( this->getData() );
  out.write ( buffer, this->size() * sizeof ( DataType ) );
}


template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_1D>::loadFromFile ( const char *filename ) {
  aol::Bzipifstream reader ( filename );
  char M;
  int ident;
  reader >> M;
  reader >> ident;
  if ( ( M != aol::VectorFileMagicChar::ScalarArray_1D ) || ( ident != aol::FileFormatMagicNumber<DataType>::FFType ) ) {
    cerr << M << ident << ", should be " << aol::VectorFileMagicChar::ScalarArray_1D << aol::FileFormatMagicNumber<DataType>::FFType
         << " (" << aol::FileFormatMagicNumber<DataType>::FFContainedName << ")" << endl;
    throw aol::Exception ( "Illegal magic number for qc::ScalarArray1D", __FILE__, __LINE__ );
  }
  aol::READ_COMMENTS ( reader );
  int sizeX;
  reader >> sizeX;
  char buffer[1024];
  reader.getline ( buffer, 1 );
  reader.getline ( buffer, 1024 ); // read and ignore maximum value

  this->reallocate ( sizeX );
  reader.read ( reinterpret_cast<char*>( this->getData() ), this->size() * sizeof( DataType ) );
}


// do not use bool and char
template class qc::ScalarArray<unsigned char,  qc::QC_1D>;
template class qc::ScalarArray<signed char,    qc::QC_1D>;
template class qc::ScalarArray<short,          qc::QC_1D>;
template class qc::ScalarArray<unsigned short, qc::QC_1D>;
template class qc::ScalarArray<int,            qc::QC_1D>;
template class qc::ScalarArray<unsigned int,   qc::QC_1D>;
template class qc::ScalarArray<float,          qc::QC_1D>;
template class qc::ScalarArray<double,         qc::QC_1D>;
template class qc::ScalarArray<long double,    qc::QC_1D>;

//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------


template <typename _DataType>
qc::ScalarArray<_DataType, qc::QC_2D>:: ScalarArray ( const string &filename ) : Array<DataType> ( ) {
  load ( filename.c_str () );
  init ();
}

template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_2D>::save ( ostream &out, qc::SaveType type, const char *comment ) const {

  if ( type == PNG_2D )
    throw aol::TypeException ( "qc::ScalarArray<DataType, qc::QC_2D>::save( ostream&, qc::SaveType, const char* ): impossible with type PNG_2D", __FILE__, __LINE__ );

  DataType  maximum = this->getMaxValue(), minimum = this->getMinValue();

  char      info[4096];

  if ( !comment ) {
    sprintf ( info, "# This is a QuOcMesh file of type %d (=%s) written %s", type, FILE_TYPE[type],
              aol::generateCurrentTimeAndDateString().c_str() );
    comment = info;
  }

  // make sure comment starts with '#'
  if ( comment[0] != '#' ) {
    info[0] = '#';
    strncpy ( info + 1, comment, 4095 );
    info[strlen ( comment ) ] = 0;
    comment = info;
  }

  out << "P" << type << endl <<  comment << endl << this->numX << " " << this->numY << endl;
  switch ( type ) {
    case PGM_UNSIGNED_CHAR_ASCII:
    case PGM_UNSIGNED_CHAR_BINARY:
      out << "255\n";
      break;
    case PGM_FLOAT_ASCII:
    case PGM_FLOAT_BINARY:
    case PGM_DOUBLE_BINARY:
      out << static_cast<int> ( maximum ) << "\n";
      break;
    default:
      throw aol::TypeException ( "qc::ScalarArray<DataType, qc::QC_2D>::save: illegal type", __FILE__, __LINE__ );
  }

  this->saveRaw ( out, type, minimum, maximum );
}

template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_2D>::save ( const char *fileName,
                                                   qc::SaveType type,
                                                   const char *comment ) const {

  if ( type == PNG_2D )
    if ( comment != NULL )
      throw aol::TypeException ( "qc::ScalarArray<DataType, qc::QC_2D>::save( const char*, qc::SaveType, const char* ) with type PNG allows no comment", __FILE__, __LINE__ );
    else
      savePNG ( fileName );
  else {
    if ( type != 2 && type != 5 && ! ( type >= 7 && type <= 9 ) )
      throw aol::TypeException ( "qc::ScalarArray<DataType, qc::QC_2D>::save: Unknown magic number",
                                 __FILE__, __LINE__ );

    aol::Bzipofstream *out = new aol::Bzipofstream ( fileName );

    if ( !out->good() ) {
      delete out;
      throw aol::FileException ( "qc::ScalarArray<DataType, qc::QC_2D>::save: Cannot open file for writing",
                                 __FILE__, __LINE__ );
    }

    save ( *out, type, comment );
    delete out;
    if ( !this->quietMode ) cerr << "done.\n";
  }
}


template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_2D>::saveToFile ( const char *filename ) const {
  aol::Bzipofstream out ( filename );
  out << aol::VectorFileMagicChar::ScalarArray_2D << aol::FileFormatMagicNumber<DataType>::FFType << endl;
  out << "# This is a QuOcMesh file of type " << aol::VectorFileMagicChar::ScalarArray_2D << aol::FileFormatMagicNumber<DataType>::FFType
      << " storing a qc::ScalarArray<" << aol::FileFormatMagicNumber<DataType>::FFContainedName << ", QC_2D>" << endl;
  out << this->getNumX() << " " << this->getNumY() << endl;
  out << this->getMaxValue() << endl; // for compatibility with PGM image format in some cases, else can be ignored
  const char* buffer = reinterpret_cast<char*> ( this->getData() );
  out.write ( buffer, this->size() * sizeof ( DataType ) );
}


template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_2D>::loadFromFile ( const char *filename ) {
  aol::Bzipifstream reader ( filename );
  char M;
  int ident;
  reader >> M;
  reader >> ident;
  if ( ( M != aol::VectorFileMagicChar::ScalarArray_2D ) || ( ident != aol::FileFormatMagicNumber<DataType>::FFType ) ) {
    cerr << M << ident << ", should be " << aol::VectorFileMagicChar::ScalarArray_2D << aol::FileFormatMagicNumber<DataType>::FFType
         << " (" << aol::FileFormatMagicNumber<DataType>::FFContainedName << ")" << endl;
    throw aol::Exception ( "Illegal magic number for qc::ScalarArray2D", __FILE__, __LINE__ );
  }
  aol::READ_COMMENTS ( reader );
  int sizeX, sizeY;
  reader >> sizeX;
  reader >> sizeY;
  char buffer[1024];
  reader.getline ( buffer, 1 );
  reader.getline ( buffer, 1024 ); // read and ignore maximum value

  this->reallocate ( sizeX, sizeY );
  reader.read ( reinterpret_cast<char*>( this->getData() ), this->size() * sizeof( DataType ) );
}


namespace {

template <typename DataType, bool add>
void pasteHelper ( qc::ScalarArray<DataType, qc::QC_2D> &dest, const qc::ScalarArray<DataType, qc::QC_2D> &image, const int xPosition, const int yPosition ) {
  if ( dest.getNumX() < image.getNumX() + xPosition || dest.getNumY() < image.getNumY() + yPosition ) {
    throw aol::Exception ( "ScalarArray<QC_2D>::pasteFrom(): Image should fit into the current array!\n", __FILE__, __LINE__ );
  }

  const int sRanX[2] = { 0, image.getNumX() };
  const int sRanY[2] = { 0, image.getNumY() };
  const int tRanX[2] = { xPosition, image.getNumX() + xPosition };
  const int tRanY[2] = { yPosition, image.getNumY() + yPosition };

  if ( add ) {
    typename qc::ScalarArray<DataType, qc::QC_2D>::AddCopier aCopier ( dest );
    dest.fillValuesFrom ( image, sRanX, sRanY, tRanX, tRanY, aCopier );
  } else {
    typename qc::ScalarArray<DataType, qc::QC_2D>::StandardCopier sCopier;
    dest.fillValuesFrom ( image, sRanX, sRanY, tRanX, tRanY, sCopier );
  }
}

} // end of nameless namespace

template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_2D>::pasteFrom ( const qc::ScalarArray<DataType, qc::QC_2D> &Image, const int XPosition, const int YPosition ) {
  pasteHelper<_DataType, false> ( *this, Image, XPosition, YPosition );
}

template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_2D>::pasteAddFrom ( const qc::ScalarArray<DataType, qc::QC_2D> &Image, const int XPosition, const int YPosition ) {
  pasteHelper<_DataType, true> ( *this, Image, XPosition, YPosition );
}

template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_2D>::load ( istream &in ) {
  char          tmp[256];
  int           type, inWidth, inHeight, inMaxColor;

  if ( ! ( this->quietMode ) ) {
    cerr << "qc::ScalarArray<DataType, qc::QC_2D>::load( istream &in )\n";
  }

  in.get ( tmp, 255, '\n' );
  if ( tmp[0] != 'P' )
    throw aol::TypeException ( "qc::ScalarArray<DataType, qc::QC_2D>::load: wrong file format",
                               __FILE__, __LINE__ );

  switch ( tmp[1] ) {
    case '2':
      type = PGM_UNSIGNED_CHAR_ASCII;
      break;
    case '5':
      type = PGM_UNSIGNED_CHAR_BINARY;
      break;
    case '7':
      type = PGM_FLOAT_ASCII;
      break;
    case '8':
      type = PGM_FLOAT_BINARY;
      break;
    case '9':
      type = PGM_DOUBLE_BINARY;
      break;
    case 'a':
      type = PGM_UNSIGNED_SHORT_BINARY;
      break;
    case 'b':
      type = PGM_SHORT_BINARY;
      break;
    default:
      throw aol::TypeException ( "qc::ScalarArray<DataType, qc::QC_2D>::load: wrong magic number",
                                 __FILE__, __LINE__ );
  }

  aol::READ_COMMENTS ( in );
  in >> inWidth;
  if ( in.fail() )
    throw aol::FileFormatException ( "qc::ScalarArray<DataType, qc::QC_2D>::load: error reading width", __FILE__, __LINE__ );

  aol::READ_COMMENTS ( in );
  in >> inHeight;
  if ( in.fail() )
    throw aol::FileFormatException ( "qc::ScalarArray<DataType, qc::QC_2D>::load: error reading height", __FILE__, __LINE__ );

  aol::READ_COMMENTS ( in );
  in >> inMaxColor;
  if ( in.fail() )
    throw aol::FileFormatException ( "qc::ScalarArray<DataType, qc::QC_2D>::load: error reading max color", __FILE__, __LINE__ );

  in.ignore();

#ifdef VERBOSE
  cerr << "Read from header: " << inWidth << " " << inHeight << " " << inMaxColor << endl;
#endif

  this->loadRaw ( in, type, inWidth, inHeight );
}

template < typename _DataType >
void qc::ScalarArray<_DataType, qc::QC_2D>::loadRaw ( istream &in, const int Type, const int InWidth, const int InHeight ) {
  if ( InWidth != this->numX || InHeight != this->numY )
    reallocate ( InWidth, InHeight );

  this->aol::Vector<_DataType>::loadRaw ( in, Type );

  if ( in.good() ) {
    if ( !this->quietMode ) cerr << "Success. 2D Stream (" << InWidth << "x" << InHeight
                                   << ") was of type " << Type << ". " << endl;
  } else throw aol::IOException ( "qc::ScalarArray<DataType, qc::QC_2D>::loadRaw: Error reading from file", __FILE__, __LINE__ );
}


template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_2D>::load ( const char *fileName ) {
  if ( !this->quietMode ) {
    cerr << "loading from file " << fileName << endl;
  };

  if ( aol::fileNameEndsWith ( fileName, ".png" ) ) {
    loadPNG ( fileName );
  } else if ( aol::fileNameEndsWith ( fileName, ".dm3" ) || aol::fileNameEndsWith ( fileName, ".dm4" ) ) {
    loadDM3 ( fileName );
  } else if ( aol::fileNameEndsWith ( fileName, ".tif" ) || aol::fileNameEndsWith ( fileName, ".tiff" ) ) {
    loadTIFF ( fileName );
  } else {

    // Decompression is done by aol::Bzipifstream based on filename extension
    aol::Bzipifstream in ( fileName );
    load ( in );
  }
  if ( !this->quietMode ) cerr << "done." << endl;
}

// code for loadPNG and savePNG based on libpng.txt from libpng-1.2.12
template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_2D>::loadPNG ( const char *fileName ) {
#ifdef USE_LIB_PNG

  PNGLoadInterface pngInterface ( fileName, this->quietMode );

  if ( pngInterface.getWidth() != this->numX || pngInterface.getHeight() != this->numY )
    reallocate ( pngInterface.getWidth(), pngInterface.getHeight() );

  pngInterface.writeDataToScalarArray2d<DataType> ( *this );


#else
  aol::doNothingWithArgumentToPreventUnusedParameterWarning ( fileName );
  throw aol::Exception ( "Reading PNG files requires libpng. Compile with -DUSE_LIB_PNG.", __FILE__, __LINE__ );
#endif //USE_LIB_PNG
}

template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_2D>::loadTIFF ( const char *fileName ) {
#if defined ( USE_LIB_TIFF )
  TIFF* tif = TIFFOpen ( fileName, "r" );
  if ( tif ) {
    uint32 width, height;
    unsigned short bits, format;
    TIFFGetField ( tif, TIFFTAG_IMAGEWIDTH, &width );
    TIFFGetField ( tif, TIFFTAG_IMAGELENGTH, &height );
    TIFFGetField ( tif, TIFFTAG_BITSPERSAMPLE, &bits );
    TIFFGetField ( tif, TIFFTAG_SAMPLEFORMAT, &format );

    this->reallocate ( width, height );

    tsize_t scanline = TIFFScanlineSize(tif);
    tdata_t buf = _TIFFmalloc(scanline);
    for ( uint32 row = 0; row < height; ++row ) {
      TIFFReadScanline(tif, buf, row);
      for ( uint32 col = 0; col < width; ++col ) {
        if ( format == SAMPLEFORMAT_INT ) {
          if ( bits == 64 )
            this->set ( col, row, static_cast<TIFF_INT64_T*>(buf)[col] );
          else if ( bits == 32 )
            this->set ( col, row, static_cast<TIFF_INT32_T*>(buf)[col] );
          else if ( bits == 16 )
            this->set ( col, row, static_cast<TIFF_INT16_T*>(buf)[col] );
          else if ( bits == 8 )
            this->set ( col, row, static_cast<TIFF_INT8_T*>(buf)[col] );
          else
            throw aol::UnimplementedCodeException ( aol::strprintf ( "%d bits uint per sample unimplemented.", bits ).c_str(), __FILE__, __LINE__ );
        } else if ( format == SAMPLEFORMAT_IEEEFP ) {
          if ( bits == 64 )
            this->set ( col, row, static_cast<double*>(buf)[col] );
          else if ( bits == 32 )
            this->set ( col, row, static_cast<float*>(buf)[col] );
          else
            throw aol::UnimplementedCodeException ( aol::strprintf ( "%d bits uint per sample unimplemented.", bits ).c_str(), __FILE__, __LINE__ );
        } else {
          // If the format was not defined by the TIFF writer, one usually assumes UINT.
          if ( bits == 64 )
            this->set ( col, row, static_cast<TIFF_UINT64_T*>(buf)[col] );
          else if ( bits == 32 )
            this->set ( col, row, static_cast<TIFF_UINT32_T*>(buf)[col] );
          else if ( bits == 16 )
            this->set ( col, row, static_cast<TIFF_UINT16_T*>(buf)[col] );
          else if ( bits == 8 )
            this->set ( col, row, static_cast<TIFF_UINT8_T*>(buf)[col] );
          else
            throw aol::UnimplementedCodeException ( aol::strprintf ( "%d bits uint per sample unimplemented.", bits ).c_str(), __FILE__, __LINE__ );
        }
      }
    }
    _TIFFfree ( buf );
    TIFFClose ( tif );
  }
#elif defined ( USE_EXTERNAL_CIMG )
  cimg_library::CImg<_DataType> cimg;
  cimg.load_tiff ( fileName );

  this->reallocate ( cimg.width(), cimg.height() );
  this->readFromBuffer ( cimg.data() );
#else
  aol::doNothingWithArgumentToPreventUnusedParameterWarning ( fileName );
  throw aol::Exception ( "Reading TIFF files requires libtiff or external cimg.", __FILE__, __LINE__ );
#endif
}

template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_2D>::saveTIFF ( const char *fileName ) const {
#if defined ( USE_LIB_TIFF )
  TIFF* tif = TIFFOpen ( fileName, "w" );
  if ( tif ) {
    TIFFSetField ( tif, TIFFTAG_IMAGEWIDTH, this->getNumX ( ) );
    TIFFSetField ( tif, TIFFTAG_IMAGELENGTH, this->getNumY ( ) );
    TIFFSetField ( tif, TIFFTAG_SAMPLESPERPIXEL, 1 );
    TIFFSetField ( tif, TIFFTAG_SAMPLEFORMAT, SAMPLEFORMAT_IEEEFP );
    TIFFSetField ( tif, TIFFTAG_BITSPERSAMPLE, 32);
    TIFFSetField ( tif, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT );
    TIFFSetField ( tif, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG );
    TIFFSetField ( tif, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_MINISBLACK );
    
    tsize_t scanline = TIFFScanlineSize(tif);
    tdata_t buf = _TIFFmalloc(scanline);
    TIFFSetField ( tif, TIFFTAG_ROWSPERSTRIP, TIFFDefaultStripSize ( tif, scanline ) );
    for ( uint32 row = 0; row < static_cast<uint32>(this->getNumY ( )) ; ++row ) {
      for ( uint32 col = 0; col < static_cast<uint32>(this->getNumY ( )); ++col ) {
        static_cast<float*>(buf)[col] = static_cast<float> ( this->get ( col, row ) );
      }
      TIFFWriteScanline ( tif, buf, row, 0 );
    }
    
    _TIFFfree ( buf );
    TIFFClose ( tif );
  }
#elif defined ( USE_EXTERNAL_CIMG )
  cimg_library::CImg<_DataType> cimg ( this->getNumX ( ), this->getNumY ( ) );
  for ( int row = 0; row < this->getNumY ( ) ; ++row )
    for ( int col = 0; col < this->getNumX ( ) ; ++col )
      cimg(col,row) = this->get ( col, row );
  cimg.save_tiff ( fileName );
#else
  aol::doNothingWithArgumentToPreventUnusedParameterWarning ( fileName );
  throw aol::Exception ( "Writing TIFF files requires libtiff or external cimg.", __FILE__, __LINE__ );
#endif
}

template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_2D>::savePNG ( const char *fileName ) const {
#ifdef USE_LIB_PNG
  PNGSaveInterface pngInterface ( fileName );
  pngInterface.loadDataFromScalarArray2d<DataType> ( *this );
#else
  aol::doNothingWithArgumentToPreventUnusedParameterWarning ( fileName );
  throw aol::Exception ( "Writing PNG files requires libpng. Compile with -DUSE_LIB_PNG.", __FILE__, __LINE__ );
#endif //USE_LIB_PNG
}

template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_2D>::loadDM3 ( const char *FileName ) {
  qc::DM3Reader dmreader( FileName );
  dmreader.exportDataToScalarArray ( *this );
}

template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_2D>::saveVTK ( const char *fileName ) const {

  ofstream file ( fileName, fstream::out );
  if ( file.is_open() ) {
    file << "# vtk DataFile Version 2.0" << endl
         << "QuocMesh ScalarArray2D VTK export, written " <<  aol::generateCurrentTimeAndDateString().c_str() << endl
         << "ASCII" << endl
         << "DATASET STRUCTURED_POINTS" << endl
         << "DIMENSIONS " << this->getNumX() + 1 << " " << this->getNumY() + 1 << " 1" << endl
         << "ORIGIN 0 0 0" << endl
         << "SPACING 1 1 1" << endl
         << "CELL_DATA " << this->size() << endl
         << "SCALARS scalarValue float" << endl
         << "LOOKUP_TABLE default" << endl;
    for ( int y = 0; y < this->getNumY(); ++y )
      for ( int x = 0; x < this->getNumX(); ++x )
        file << get( x, y ) << endl;
  }
  else {
    throw aol::Exception ( "qc::ScalarArray<DataType,qc::QC_2D>::saveVTK: Cannot open file for writing", __FILE__, __LINE__ );
  }
}

// needed for element saturation
template <typename _DataType>
_DataType qc::ScalarArray<_DataType, qc::QC_2D>::getElementSaturationMax ( int X, int Y ) {
  DataType val, max = this->get ( X, Y );

  /* the following block is used for saturation to ensure
     1-level-transitions to neighbours connected only by a node */
  int XMin, XMax, YMin, YMax;
  if ( X > 0 ) XMin = X - 1;
  else XMin = 0;
  if ( Y > 0 ) YMin = Y - 1;
  else YMin = 0;
  if ( X < this->numX - 2 ) XMax = X + 2;
  else XMax = X + 1;
  if ( Y < this->numY - 2 ) YMax = Y + 2;
  else YMax = Y + 1;

  for ( int x = XMin; x <= XMax; ++x ) {
    for ( int y = YMin; y <= YMax; ++y ) {
      val = this->get ( x, y );
      if ( val > max ) max = val;
    }
  }
  return max;

  val = this->get ( X + 1 , Y );
  if ( val > max ) max = val;

  val = this->get ( X     , Y + 1 );
  if ( val > max ) max = val;

  val = this->get ( X + 1 , Y + 1 );
  if ( val > max ) max = val;


  return max;
}

template <typename _DataType>
_DataType qc::ScalarArray<_DataType, qc::QC_2D>::getElementSaturationMin ( int X, int Y ) {

  DataType val, min = this->get ( X, Y );

  val = this->get ( X + 1 , Y );
  if ( val < min ) min = val;

  val = this->get ( X     , Y + 1 );
  if ( val < min ) min = val;

  val = this->get ( X + 1 , Y + 1 );
  if ( val < min ) min = val;

  if ( X > 0 ) {
    val = this->get ( X - 1 , Y );
    if ( val < min ) min = val;

    val = this->get ( X - 1 , Y + 1 );
    if ( val < min ) min = val;
  }

  if ( X < this->numX - 2 ) {
    val = this->get ( X + 2 , Y );
    if ( val < min ) min = val;

    val = this->get ( X + 2 , Y + 1 );
    if ( val < min ) min = val;
  }

  if ( Y > 0 ) {
    val = this->get ( X     , Y - 1 );
    if ( val < min ) min = val;

    val = this->get ( X + 1 , Y - 1 );
    if ( val < min ) min = val;
  }

  if ( Y < this->numY - 2 ) {
    val = this->get ( X     , Y + 2 );
    if ( val < min ) min = val;

    val = this->get ( X + 1 , Y + 2 );
    if ( val < min ) min = val;
  }

  return min;
}

template <typename _DataType>
typename qc::ScalarArray<_DataType, qc::QC_2D>::DataType
qc::ScalarArray<_DataType, qc::QC_2D>::getMedianFilterValue ( int X, int Y ) const {
  int XMin, XMax, YMin, YMax;
  int Off = ( this->MED_FILTER_WIDTH - 1 ) >> 1;

  XMin = aol::Max ( 0, X - Off );
  YMin = aol::Max ( 0, Y - Off );
  XMax = aol::Min ( this->numX - 1, X + Off );
  YMax = aol::Min ( this->numY - 1, Y + Off );

  aol::Vector<DataType> medianSortVec ( ( XMax - XMin + 1 ) * ( YMax - YMin + 1 ) );
  int index = 0;
  for ( int x = XMin; x <= XMax; ++x ) {
    for ( int y = YMin; y <= YMax; ++y ) {
      medianSortVec[index++] = this->get ( x, y );
    }
  }
  medianSortVec.sortValues();
  return medianSortVec[ ( medianSortVec.size() >> 1 ) ];
}

template <typename _DataType>
typename qc::ScalarArray<_DataType, qc::QC_2D>::DataType
qc::ScalarArray<_DataType, qc::QC_2D>::getConvolveValue ( int X, int Y,
                                                          const qc::Kernel2d<RealType> &Kernel ) const {
  int ksize = Kernel.getSize();
  int offset = ksize >> 1;

  RealType val = 0.0;

  if ( X >= offset && Y >= offset &&
       X < this->numX - offset && Y < this->numY - offset ) {

    int ind = qc::Array<DataType>::index ( X - offset, Y - offset );
    int Yoff = ( this->numX - ksize );
    int kindmax = ksize * ksize;

    for ( int kind = 0; kind < kindmax; ) {
      val += static_cast<RealType> ( Kernel.get ( kind++ ) * this->get ( ind++ ) );
      if ( kind % ksize == 0 ) ind += Yoff;
    }
  } else {
    for ( int OffY = -offset; OffY <= offset; ++OffY ) {
      for ( int OffX = -offset; OffX <= offset; ++OffX ) {
        val += static_cast<RealType> ( Kernel.getValue ( OffX, OffY )
                                       * getPeriodic ( X + OffX, Y + OffY ) );
      }
    }
  }
  return static_cast< DataType > ( val );
}

template <typename _DataType>
typename qc::ScalarArray<_DataType, qc::QC_2D>::DataType
qc::ScalarArray<_DataType, qc::QC_2D>::getWeightedConvolveValue ( int X, int Y,
                                                                  const qc::Kernel2d<RealType> &Kernel,
                                                                  const qc::ScalarArray<RealType, qc::QC_2D> &Weight ) const {
  int ksize = Kernel.getSize();
  int offset = ksize >> 1;

  RealType val = 0.0, weight = 0.0, inc;

  if ( X >= offset && Y >= offset &&
       X < this->numX - offset && Y < this->numY - offset ) {

    int ind = qc::Array<DataType>::index ( X - offset, Y - offset );
    int Yoff = ( this->numX - ksize );
    int kindmax = ksize * ksize;

    for ( int kind = 0; kind < kindmax; ) {
      inc = Kernel.get ( kind++ ) * Weight.get ( ind );
      val += static_cast<RealType> ( inc * this->get ( ind++ ) );
      weight += inc;
      if ( kind % ksize == 0 ) ind += Yoff;
    }
  } else {
    for ( int OffY = -offset; OffY <= offset; ++OffY ) {
      for ( int OffX = -offset; OffX <= offset; ++OffX ) {
        inc = Kernel.getValue ( OffX, OffY ) * Weight.getPeriodic ( X + OffX, Y + OffY );
        val += static_cast<RealType> ( inc * getPeriodic ( X + OffX, Y + OffY ) );
        weight += inc;
      }
    }
  }
  return static_cast< DataType > ( ( aol::Abs ( weight ) > 1e-10 ) ? ( val / weight ) : 0 );
}

template <typename _DataType>
typename qc::ScalarArray<_DataType, qc::QC_2D>::DataType
qc::ScalarArray<_DataType, qc::QC_2D>::getWeightedConvolveValue ( int X, int Y, int /*Z*/,
                                                                  const qc::Kernel2d<RealType> &Kernel,
                                                                  const qc::ScalarArray<RealType, qc::QC_2D> &Weight ) const {
  return getWeightedConvolveValue ( X, Y, Kernel, Weight );
}

template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_2D>::applyMedianFilter() {
  qc::ScalarArray<DataType, qc::QC_2D> tmp ( *this );

  for ( int X = 0; X < this->numX; ++X ) {
    for ( int Y = 0; Y < this->numY; ++Y ) {
      set ( X, Y, tmp.getMedianFilterValue ( X, Y ) );
    }
  }
}

template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_2D>::saltAndPepperNoise ( double Frac ) {
  for ( int X = 0; X < this->numX; ++X ) {
    for ( int Y = 0; Y < this->numY; ++Y ) {
      double r = static_cast< double > ( rand() ) / ( static_cast< double > ( RAND_MAX ) );
      if ( r < Frac ) {
        if ( r < Frac / 2.0 ) {
          set ( X, Y, static_cast< DataType > ( 1.0 ) );
        } else {
          set ( X, Y, static_cast< DataType > ( 0.0 ) );
        }
      }
    }
  }
}

template <typename _DataType>
typename qc::ScalarArray<_DataType, qc::QC_2D>::DataType
qc::ScalarArray<_DataType, qc::QC_2D>::getCannyEdgeValue ( int X, int Y ) const {
  RealType dX, dY, dXX, dXY, dYY;

  RealType gradNorm;

  dXX = static_cast<RealType> ( getConvolveValue ( X, Y, *diffKernelXX ) );
  dYY = static_cast<RealType> ( getConvolveValue ( X, Y, *diffKernelYY ) );
  dXY = static_cast<RealType> ( getConvolveValue ( X, Y, *diffKernelXY ) );

  dX = static_cast<RealType> ( getConvolveValue ( X, Y, *diffKernelX ) );
  dY = static_cast<RealType> ( getConvolveValue ( X, Y, *diffKernelY ) );

  gradNorm = sqrt ( dX * dX + dY * dY );

  dX /= gradNorm;
  dY /= gradNorm;

  RealType c = dXX * dX * dX + 2 * dXY * dY * dX + dYY * dY * dY;

  return static_cast< DataType > ( c );
}

template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_2D>::setDiffSigma ( RealType DiffSigma ) {
  if ( !diffKernelX || diffKernelX->getSigma() != DiffSigma ) {
    if ( !diffKernelX ) {
      diffKernelX = new qc::GaussDiffKernel2d<RealType> ( this->DIFF_STENCIL, DiffSigma, DIFF_X );
    } else {
      diffKernelX->setSigma ( DiffSigma );
    }
  }

  if ( !diffKernelY || diffKernelY->getSigma() != DiffSigma ) {
    if ( !diffKernelY ) {
      diffKernelY = new qc::GaussDiffKernel2d<RealType> ( this->DIFF_STENCIL, DiffSigma, DIFF_Y );
    } else {
      diffKernelY->setSigma ( DiffSigma );
    }
  }

  if ( !diffKernelXX || diffKernelXX->getSigma() != DiffSigma ) {
    if ( !diffKernelXX ) {
      diffKernelXX = new qc::GaussDiffKernel2d<RealType> ( this->DIFF_STENCIL, DiffSigma, DIFF_XX );
    } else {
      diffKernelXX->setSigma ( DiffSigma );
    }
  }

  if ( !diffKernelYY || diffKernelYY->getSigma() != DiffSigma ) {
    if ( !diffKernelYY ) {
      diffKernelYY = new qc::GaussDiffKernel2d<RealType> ( this->DIFF_STENCIL, DiffSigma, DIFF_YY );
    } else {
      diffKernelYY->setSigma ( DiffSigma );
    }
  }

  if ( !diffKernelXY || diffKernelXY->getSigma() != DiffSigma ) {
    if ( !diffKernelXY ) {
      diffKernelXY = new qc::GaussDiffKernel2d<RealType> ( this->DIFF_STENCIL, DiffSigma, DIFF_XY );
    } else {
      diffKernelXY->setSigma ( DiffSigma );
    }
  }
}

template <typename _DataType>
ostream & qc::ScalarArray<_DataType, qc::QC_2D>::print ( ostream &os ) const {
  if ( prettyFormat ) os << "ScalarArray2D (" << this->numX << " x " << this->numY << "):" << endl;

  os << format ( ' ' ) << ( prettyFormat ? "|" : " " );
  for ( int j = 0 ; j < this->numY ; ++j ) {
    os << format ( j ) << " ";
  }
  os << endl;

  if ( prettyFormat ) {
    os << format ( '-' ) << "+";
    for ( int j = 0 ; j < this->numY ; ++j ) {
      os << format ( '-' ) << "-";
    }
    os << endl;
  }

  for ( int i = 0 ; i < this->numX ; ++i ) {
    os << format ( i ) << ( prettyFormat ? "|" : " " );
    for ( int j = 0 ; j < this->numY ; ++j )
      os << format ( this->get ( i, j ) ) << " ";
    os << endl;
  }
  return os;
}

template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_2D>::bicubicResampleFrom ( const ScalarArray<_DataType, qc::QC_2D> &Image ) {
  qc::CachedBicubicInterpolator<RealType> interpolator;
  // p stores the 16 values of Image necessary for the bicubic interpolation at the current position
  aol::Matrix44<RealType> p;
  int xPosCached = -1;
  int yPosCached = -1;
  const int imageNumXMinus1 = Image.getNumX() - 1;
  const int imageNumYMinus1 = Image.getNumY() - 1;
  const RealType xFac = ( imageNumXMinus1 ) / static_cast<RealType> ( this->numX - 1 );
  const RealType yFac = ( imageNumYMinus1 ) / static_cast<RealType> ( this->numY - 1 );
  for ( int i = 0; i < this->numX; ++i ) {
    for ( int j = 0; j < this->numY; ++j ) {
      const RealType x = i * xFac;
      const RealType y = j * yFac;
      const int xPos = static_cast<int> ( x );
      const int yPos = static_cast<int> ( y );

      if ( ( xPosCached != xPos ) || ( yPosCached != yPos ) ) {
        for ( int k = 0; k < 4; ++k )
          for ( int l = 0; l < 4; ++l )
            p[k][l] = Image.get ( aol::Clamp ( xPos + k - 1, 0, imageNumXMinus1 ), aol::Clamp ( yPos + l - 1, 0, imageNumYMinus1 ) );
        interpolator.updateCoefficients ( p );
        xPosCached = xPos;
        yPosCached = yPos;
      }
      set ( i, j, static_cast<_DataType> ( interpolator.getValue ( x - xPos, y - yPos ) ) );
    }
  }
}

// template class qc::ScalarArray<bool, qc::QC_2D>; // qc::ScalarArray<bool, qc::QC_2D> should not be used ( += doesn't make sense with bool ), use qc::BitArray<qc::QC_2D> instead.
// template class qc::ScalarArray<char, qc::QC_2D>; // DO NOT USE. char is signed or unsigned, depending on your platform. use signed and unsigned char instead.
template class qc::ScalarArray<signed char,    qc::QC_2D>;
template class qc::ScalarArray<unsigned char,  qc::QC_2D>;
template class qc::ScalarArray<short,          qc::QC_2D>;
template class qc::ScalarArray<unsigned short, qc::QC_2D>;
template class qc::ScalarArray<int,            qc::QC_2D>;
template class qc::ScalarArray<unsigned int,   qc::QC_2D>;
template class qc::ScalarArray<float,          qc::QC_2D>;
template class qc::ScalarArray<double,         qc::QC_2D>;
template class qc::ScalarArray<long double,    qc::QC_2D>;



//--------------------------------------------------------------------------------
//--------------------------------------------------------------------------------

template <typename _DataType>
qc::ScalarArray<_DataType, qc::QC_3D>::ScalarArray ( const string &filename ) : Array<DataType> ( ) {
  load ( filename.c_str () );
  initPolyBases ();
  init ();
}


template <typename _DataType> qc::ScalarArray<_DataType, qc::QC_3D>::~ScalarArray() {
  if ( diffKernelX )
    delete diffKernelX;

  if ( diffKernelY )
    delete diffKernelY;

  if ( diffKernelZ )
    delete diffKernelZ;
}

template < typename _DataType >
void qc::ScalarArray < _DataType, qc::QC_3D >::save ( ostream &out, qc::SaveType type, const char *comment ) const {
  if ( this->size() == 0 )
    throw aol::FileFormatException ( "qc::ScalarArray<DataType, qc::QC_3D>::save: Can't save an array with size() == 0", __FILE__, __LINE__ );

  DataType  maximum = this->getMaxValue();
  DataType  minimum = this->getMinValue();
  char      info[4096];

  if ( type == PNG_2D )
    throw aol::TypeException ( "qc::ScalarArray<DataType, qc::QC_3D>::save: impossible with type PNG_2D", __FILE__, __LINE__ );

  if ( !comment ) {
    sprintf ( info, "# This is a QuOcMesh file of type %d (=%s) written %s", type, FILE_TYPE[type], aol::generateCurrentTimeAndDateString().c_str() );
    comment = info;
  }

  if ( comment[0] != '#' ) {
    info[0] = '#';
    strncpy ( info + 1, comment, 4095 );
    info[strlen ( comment ) ] = 0;
  }

  out << "Q";
  if ( type == 10 )
    out << "a";
  else
    out << type;
  out << "\n"
      <<  comment << "\n" << this->numX << " " << this->numY << " " << this->numZ;

  switch ( type ) {
    case PGM_UNSIGNED_CHAR_ASCII:
      out << std::endl << static_cast<int> ( maximum ) << std::endl;
      break;
    case PGM_UNSIGNED_CHAR_BINARY:
      out << "\n255\n";
      break;
    case PGM_FLOAT_ASCII:
    case PGM_FLOAT_BINARY:
    case PGM_DOUBLE_BINARY:
    case PGM_UNSIGNED_SHORT_BINARY:
      out << std::endl << static_cast<int> ( maximum ) << std::endl;
      break;
    default:
      throw aol::TypeException ( "qc::ScalarArray<QC_3D>::save: invalid SaveType specified", __FILE__, __LINE__ );
  }

  this->saveRaw ( out, type, minimum, maximum );
}


template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_3D>::saveToFile ( const char *filename ) const {
  aol::Bzipofstream out ( filename );
  out << aol::VectorFileMagicChar::ScalarArray_3D << aol::FileFormatMagicNumber<DataType>::FFType << endl;
  out << "# This is a QuOcMesh file of type " << aol::VectorFileMagicChar::ScalarArray_3D << aol::FileFormatMagicNumber<DataType>::FFType
      << " storing a qc::ScalarArray<" << aol::FileFormatMagicNumber<DataType>::FFContainedName << ", QC_3D>" << endl;
  out << this->getNumX() << " " << this->getNumY() << " " << this->getNumZ() << endl;
  out << this->getMaxValue() << endl; // for compatibility with PGM image format in some cases, else can be ignored
  const char* buffer = reinterpret_cast<char*> ( this->getData() );
  out.write ( buffer, this->size() * sizeof ( DataType ) );
}


template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_3D>::loadFromFile ( const char *filename ) {
  aol::Bzipifstream reader ( filename );
  char M;
  int ident;
  reader >> M;
  reader >> ident;
  if ( ( M != aol::VectorFileMagicChar::ScalarArray_3D ) || ( ident != aol::FileFormatMagicNumber<DataType>::FFType ) ) {
    cerr << M << ident << ", should be " << aol::VectorFileMagicChar::ScalarArray_3D << aol::FileFormatMagicNumber<DataType>::FFType
         << " (" << aol::FileFormatMagicNumber<DataType>::FFContainedName << ")" << endl;
    throw aol::Exception ( "Illegal magic number for qc::ScalarArray3D", __FILE__, __LINE__ );
  }
  aol::READ_COMMENTS ( reader );
  int sizeX, sizeY, sizeZ;
  reader >> sizeX;
  reader >> sizeY;
  reader >> sizeZ;
  char buffer[1024];
  reader.getline ( buffer, 1 );
  reader.getline ( buffer, 1024 ); // read and ignore maximum value

  this->reallocate ( sizeX, sizeY, sizeZ );
  reader.read ( reinterpret_cast<char*>( this->getData() ), this->size() * sizeof( DataType ) );
}

template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_3D>::loadSlices ( const char *fileNameMask, qc::Comp dir, int begin, int end ) {
  int inWidth = 0, inHeight = 0, depth = 0;
  switch ( dir ) {
    case QC_X:
      inWidth = this->numY;
      inHeight = this->numZ;
      depth = this->numX;
      break;
    case QC_Y:
      inWidth = this->numX;
      inHeight = this->numZ;
      depth = this->numY;
      break;
    case QC_Z:
      inWidth = this->numX;
      inHeight = this->numY;
      depth = this->numZ;
      break;
    default:
        throw aol::UnimplementedCodeException ( "Unknown qc::Comp", __FILE__, __LINE__ );
  }
  qc::ScalarArray<DataType, qc::QC_2D> slice ( inWidth, inHeight );
  slice.setQuietMode ( this->quietMode );
  if ( end - begin + 1 > depth ) {
    throw aol::Exception ( "too many slices to fit in the ScalarArray", __FILE__, __LINE__ );
  }
  // center the slices in the array:
  const int offset = ( depth - ( end - begin ) ) / 2;
  char filename[1024];
  for ( int s = begin; s <= end; ++s ) {
    sprintf ( filename, fileNameMask, s );
    slice.load ( filename );
    this->putSlice ( dir, s - begin + offset, slice );
  }
}


template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_3D>::
saveSlices ( const char *fileNameMask,
             qc::Comp dir,
             qc::SaveType type,
             const char *comment,
             aol::OverflowHandlingType oh,
             typename qc::ScalarArray<_DataType, qc::QC_3D>::DataType mi,
             typename qc::ScalarArray<_DataType, qc::QC_3D>::DataType ma ) const {


  int tmpNumX = 0, tmpNumY = 0, numOfSlices = 0;
  switch ( dir ) {
    case QC_X:
      tmpNumX = this->getNumY();
      tmpNumY = this->getNumZ();
      numOfSlices = this->getNumX();
      break;
    case QC_Y:
      tmpNumX = this->getNumX();
      tmpNumY = this->getNumZ();
      numOfSlices = this->getNumY();
      break;
    case qc::QC_Z:
      tmpNumX = this->getNumX();
      tmpNumY = this->getNumY();
      numOfSlices = this->getNumZ();
      break;
    default:
      throw aol::UnimplementedCodeException ( "Unknown qc::Comp", __FILE__, __LINE__ );
  }

  qc::ScalarArray<DataType, qc::QC_2D> tmp ( tmpNumX, tmpNumY );
  tmp.setOverflowHandling ( oh, mi, ma );

  tmp.quietMode = this->quietMode;

  for ( int z = 0; z < numOfSlices; ++z ) {
    char fileName[1024];
    sprintf ( fileName, fileNameMask, z );
    this->getSlice ( dir, z, tmp );
    tmp.save ( fileName, type, comment );
  }
  if ( !this->quietMode )
    cerr << "Successfully wrote to 2D PGM " << type << " stream slices\n";

}


template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_3D>::saveMetaImageFile ( const char *BaseFileName, SaveType Type ) const {
  if ( Type != PGM_UNSIGNED_CHAR_BINARY &&
       Type != PGM_UNSIGNED_SHORT_BINARY  )
    throw aol::Exception ( "qc::ScalarArray<DataType, qc::QC_3D>::saveMetaImageFile: Unsupported SaveType",
                           __FILE__, __LINE__ );

  string headerFileName = aol::strprintf ( "%s.mhd", BaseFileName );
  string dataFileName = aol::strprintf ( "%s.raw", BaseFileName );
  aol::Bzipofstream out ( headerFileName.c_str() );

  if ( !out.good() )
    throw aol::Exception ( "qc::ScalarArray<DataType, qc::QC_3D>::saveMetaImageFile: Cannot open file for writing", __FILE__, __LINE__ );

  out << "NDims = 3\n";
  out << "DimSize = " << this->numX << " " << this->numY << " " << this->numZ << endl;
  out << "ElementType = ";
  if ( Type == PGM_UNSIGNED_CHAR_BINARY )
    out << "MET_UCHAR\n";
  else if ( Type == PGM_UNSIGNED_SHORT_BINARY )
    out << "MET_USHORT\n";
  else
    throw aol::Exception ( "qc::ScalarArray<DataType, qc::QC_3D>::saveMetaImageFile: Anomalous error encountered!", __FILE__, __LINE__ );

  out << "ElementSpacing = 1.0 1.0 1.0\n";
  out << "ElementByteOrderMSB = False\n";
  out << "ElementDataFile = " << dataFileName << endl;

  this->saveRaw ( dataFileName.c_str(), Type, this->getMinValue(), this->getMaxValue() );
}


template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_3D>::save ( const char *fileName,
                                                   qc::SaveType type, const char *comment ) const {

  if ( type == PNG_2D )
    throw aol::TypeException ( "qc::ScalarArray<DataType, qc::QC_3D>::save: impossible with type PNG_2D", __FILE__, __LINE__ );

  if ( type != PGM_UNSIGNED_CHAR_ASCII &&
       type != PGM_UNSIGNED_CHAR_BINARY &&
       type != PGM_FLOAT_ASCII &&
       type != PGM_FLOAT_BINARY &&
       type != PGM_DOUBLE_BINARY &&
       type != PGM_UNSIGNED_SHORT_BINARY  )
    throw aol::Exception ( "qc::ScalarArray<DataType, qc::QC_3D>::save: Unknown SaveType",
                           __FILE__, __LINE__ );

  aol::Bzipofstream *out = new aol::Bzipofstream ( fileName );

  if ( !out->good() ) {
    delete out;
    throw aol::Exception ( "qc::ScalarArray<DataType, qc::QC_3D>::save: Cannot open file for writing", __FILE__, __LINE__ );
  }

  save ( *out, type, comment );

  delete out;

  if ( !this->quietMode ) cerr << "done.\n";
}

//-----------------------------------------------------------------------------------------------------
template < typename _DataType >
void qc::ScalarArray<_DataType, qc::QC_3D>::load ( istream &in ) {
  char tmp[256];
  int  type, inWidth, inHeight, inDepth, inMaxColor;

  if ( !in ) {
    throw aol::Exception ( "qc::ScalarArray<DataType, qc::QC_3D>::load: instream not open. ", __FILE__, __LINE__ );
  };

  in.get ( tmp, 255, '\n' );
  if ( tmp[0] != 'Q' ) {
    cerr << "first two bytes: " << tmp[0] << tmp[1] << endl;
    throw aol::Exception ( "qc::ScalarArray<DataType, qc::QC_3D>::load: wrong file format", __FILE__, __LINE__ );
  }

  switch ( tmp[1] ) {
    case '2':
      type = PGM_UNSIGNED_CHAR_ASCII;
      break;
    case '5':
      type = PGM_UNSIGNED_CHAR_BINARY;
      break;
    case '7':
      type = PGM_FLOAT_ASCII;
      break;
    case '8':
      type = PGM_FLOAT_BINARY;
      break;
    case '9':
      type = PGM_DOUBLE_BINARY;
      break;
    case 'a':
      type = PGM_UNSIGNED_SHORT_BINARY;
      break;
    case 'b':
      type = PGM_SHORT_BINARY;
      break;
    default:
      throw aol::TypeException ( "qc::ScalarArray<DataType, qc::QC_3D>::load: wrong magic number", __FILE__, __LINE__ );
  }

  aol::READ_COMMENTS ( in );
  in >> inWidth;
  aol::READ_COMMENTS ( in );
  in >> inHeight;
  aol::READ_COMMENTS ( in );
  in >> inDepth;
  aol::READ_COMMENTS ( in );
  in >> inMaxColor;
  in.ignore();

  if ( !this->quietMode ) cerr << "\twidth = " << inWidth << endl
                                 << "\theight = " << inHeight << endl
                                 << "\tdepth = " << inDepth << endl
                                 << "\tmaxColor = " << inMaxColor << endl;

  if ( !this->quietMode ) cerr << this->getNumX() << " " << this->getNumY() << " " << this->getNumZ() << endl;

  this->loadRaw ( in, type, inWidth, inHeight, inDepth );
}

template < typename _DataType >
void qc::ScalarArray<_DataType, qc::QC_3D>::loadRaw ( istream &in, const int Type, const int InWidth, const int InHeight, const int InDepth ) {
  if ( InWidth != this->getNumX() || InHeight != this->getNumY() || InDepth != this->getNumZ() )
    reallocate ( InWidth, InHeight, InDepth );

  this->aol::Vector<_DataType>::loadRaw ( in, Type );

  if ( in.good() ) {
    if ( !this->quietMode ) cerr << "Success. 3D Stream (" << InWidth << "x" << InHeight << "x" << InDepth
                                   << ") was of type " << Type << ". " << endl;
  } else throw aol::Exception ( "qc::ScalarArray<DataType, qc::QC_3D>::loadRaw: Error reading from file", __FILE__, __LINE__ );
}


// This method is exactly implemented as in qc::ScalarArray<QC_2D>.
template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_3D>::load ( const char *fileName ) {
#ifdef DEBUG
  cerr << "loading " << fileName << endl;
#endif
  if ( aol::fileNameEndsWith ( fileName, ".mrc" ) ) {
    loadMRC ( fileName );
  } else {
    aol::Bzipifstream in ( fileName );
    load ( in );
  }
}


template < typename _DataType >
typename qc::ScalarArray<_DataType, qc::QC_3D>::RealType
qc::ScalarArray < _DataType, qc::QC_3D >::integrateOnSquare ( RealType OffX, RealType OffY, RealType OffZ, int BasisFuncNum, const Hexahedron &Hex ) const {
  const RealType alpha[2] = {
    0.21132487, 0.78867513
  };
  RealType temp = 0.0;

  for ( int i = 0; i < 2; ++i ) {
    for ( int j = 0; j < 2; ++j ) {
      for ( int k = 0; k < 2; ++k ) {
        temp += Hex.localEvaluate ( alpha[i], alpha[j], alpha[k] ) *
                polyBases[BasisFuncNum] ( OffX + alpha[i], OffY + alpha[j], OffZ + alpha[k] );
      }
    }
  }
  return static_cast<RealType> ( temp / 8.0 );
}


template < typename _DataType >
typename qc::ScalarArray<_DataType, qc::QC_3D>::RealType qc::ScalarArray < _DataType, qc::QC_3D >
::getConvolveValue ( int X, int Y, int Z, const qc::Kernel3d < RealType > & Kernel ) const {
  int ksize = Kernel.getSize();
  int offset = ksize >> 1;

  RealType val = 0.0;

  if ( X >= offset && Y >= offset && Z >= offset && X < this->getNumX() - offset && Y < this->getNumY() - offset &&
       Z < this->getNumZ() - offset ) {

    int ind = this->index ( X - offset, Y - offset, Z - offset );
    int Yoff = this->getNumX() - ksize;
    int Zoff = this->getNumX() * ( this->getNumY() - ksize );
    int kindmax = ksize * ksize * ksize;
    int ksizesqr = ksize * ksize;

    for ( int kind = 0; kind < kindmax; ) {
      val += static_cast<RealType> ( Kernel.get ( kind++ ) * ( this->get ( ind++ ) ) );
      if ( kind % ksize == 0 ) ind += Yoff;
      if ( kind % ksizesqr == 0 ) ind += Zoff;
    }
  } else {
    for ( int OffZ = -offset; OffZ <= offset; ++OffZ ) {
      for ( int OffY = -offset; OffY <= offset; ++OffY ) {
        for ( int OffX = -offset; OffX <= offset; ++OffX ) {
          val += static_cast<RealType> ( Kernel.getValue ( OffX, OffY, OffZ ) * this->getPeriodic ( X + OffX, Y + OffY, Z + OffZ ) );
        }
      }
    }
  }

#if 0
  cerr << "template <typename DataType> DataType qc::ScalarArray<DataType, qc::QC_3D>" <<
       "::getConvolveValue( int X, int Y, int Z, const qc::Kernel3d<RealType> &Kernel ) const" << "Uncompilable code!!!\n";
#endif
  return val;
}


template < typename _DataType >
typename qc::ScalarArray<_DataType, qc::QC_3D>::RealType qc::ScalarArray < _DataType, qc::QC_3D >
::getWeightedConvolveValue ( int X, int Y, int Z, const qc::Kernel3d < RealType > & Kernel, const qc::ScalarArray < RealType, qc::QC_3D > & Weight ) const {
  int ksize = Kernel.getSize();
  int offset = ksize >> 1;

  RealType val = 0.0, weight = 0.0, inc;

  if ( X >= offset && Y >= offset && Z >= offset && X < this->getNumX() - offset && Y < this->getNumY() - offset &&
       Z < this->getNumZ() - offset ) {

    int ind = this->index ( X - offset, Y - offset, Z - offset );
    int Yoff = this->getNumX() - ksize;
    int Zoff = this->getNumX() * ( this->getNumY() - ksize );
    int kindmax = ksize * ksize * ksize;
    int ksizesqr = ksize * ksize;

    for ( int kind = 0; kind < kindmax; ) {
      inc = Kernel.get ( kind++ ) * Weight.get ( ind );
      val += static_cast<RealType> ( inc * this->get ( ind++ ) );
      weight += inc;
      if ( kind % ksize == 0 ) ind += Yoff;
      if ( kind % ksizesqr == 0 ) ind += Zoff;
    }
  } else {
    for ( int OffZ = -offset; OffZ <= offset; ++OffZ ) {
      for ( int OffY = -offset; OffY <= offset; ++OffY ) {
        for ( int OffX = -offset; OffX <= offset; ++OffX ) {
          inc = Kernel.getValue ( OffX, OffY, OffZ ) * Weight.getPeriodic ( X + OffX, Y + OffY, Z + OffZ );
          val += static_cast<RealType> ( inc * this->getPeriodic ( X + OffX, Y + OffY, Z + OffZ ) );
          weight += inc;
        }
      }
    }
  }
  return static_cast< DataType > ( aol::Abs ( weight ) > 1e-10 ) ? ( val / weight ) : 0;
}

//! \brief Resizes the array. If the array is made smaller, values are omitted.
//! \author toelkes
template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_3D>::resize ( const int newNumX, const int newNumY, const int newNumZ ) {
  qc::ScalarArray<DataType, qc::QC_3D> copyOfThis ( *this, aol::DEEP_COPY );
  this->reallocate ( newNumX, newNumY, newNumZ );

  const int sRanX[2] = { 0, aol::Min ( copyOfThis.getNumX(), newNumX ) };
  const int sRanY[2] = { 0, aol::Min ( copyOfThis.getNumY(), newNumY ) };
  const int sRanZ[2] = { 0, aol::Min ( copyOfThis.getNumZ(), newNumZ ) };

  // target ranges are the same, as at most the whole old array has to be copied into the new (resized) one
  const int tRanX[2] = { sRanX[0], sRanX[1] };
  const int tRanY[2] = { sRanY[0], sRanY[1] };
  const int tRanZ[2] = { sRanZ[0], sRanZ[1] };

  StandardCopier copier;
  this->fillValuesFrom ( copyOfThis, sRanX, sRanY, sRanZ, tRanX, tRanY, tRanZ, copier );
}


template<typename _DataType>
_DataType qc::ScalarArray<_DataType, qc::QC_3D>::getElementSaturationMax ( int X, int Y, int Z ) const {
  DataType val,
           max = this->get ( X, Y, Z );

  int Xmin, Xmax, Ymin, Ymax, Zmin, Zmax;

  if ( X > 0 ) Xmin = X - 1;
  else Xmin = 0;
  if ( Y > 0 ) Ymin = Y - 1;
  else Ymin = 0;
  if ( Z > 0 ) Zmin = Z - 1;
  else Zmin = 0;
  if ( X + 2 < this->getNumX() ) Xmax = X + 2;
  else Xmax = this->getNumX() - 1;
  if ( Y + 2 < this->getNumY() ) Ymax = Y + 2;
  else Ymax = this->getNumY() - 1;
  if ( Z + 2 < this->getNumZ() ) Zmax = Z + 2;
  else Zmax = this->getNumZ() - 1;

  for ( int z = Zmin; z <= Zmax; ++z ) {
    for ( int y = Ymin; y <= Ymax; ++y ) {
      for ( int x = Xmin; x <= Xmax; ++x ) {
        if ( ( val = this->get ( x, y, z ) ) > max ) max = val;
      }
    }
  }

  // is this comment block still needed?

  /*val = get( X + 1 , Y     , Z );
    if ( val > max ) max = val;

    val = get(  X     , Y + 1 , Z );
    if ( val > max ) max = val;

    val = get(  X + 1 , Y + 1 , Z );
    if ( val > max ) max = val;


    val = get(  X + 1 , Y     , Z + 1 );
    if ( val > max ) max = val;

    val = get(  X     , Y     , Z + 1 );
    if ( val > max ) max = val;

    val = get(  X + 1 , Y + 1 , Z + 1 );
    if ( val > max ) max = val;

    val = get(  X     , Y + 1 , Z + 1 );
    if ( val > max ) max = val;


    if ( X > 0 ) {
    val = get(  X - 1 , Y     , Z );
    if ( val > max ) max = val;

    val = get(  X - 1 , Y + 1 , Z );
    if ( val > max ) max = val;

    val = get(  X - 1 , Y     , Z + 1 );
    if ( val > max ) max = val;

    val = get(  X - 1 , Y + 1 , Z + 1 );
    if ( val > max ) max = val;

    }

    if ( X < getNumX() - 2 ) {
    val = get(  X + 2 , Y     , Z );
    if ( val > max ) max = val;

    val = get(  X + 2 , Y + 1 , Z );
    if ( val > max ) max = val;

    val = get(  X + 2 , Y     , Z + 1 );
    if ( val > max ) max = val;

    val = get(  X + 2 , Y + 1 , Z + 1 );
    if ( val > max ) max = val;

    }

    if ( Y > 0 ) {
    val = get(  X     , Y - 1 , Z );
    if ( val > max ) max = val;

    val = get(  X + 1 , Y - 1 , Z );
    if ( val > max ) max = val;

    val = get(  X     , Y - 1 , Z + 1 );
    if ( val > max ) max = val;

    val = get(  X + 1 , Y - 1 , Z + 1 );
    if ( val > max ) max = val;

    }

    if ( Y < getNumY() - 2 ) {
    val = get(  X     , Y + 2 , Z );
    if ( val > max ) max = val;

    val = get(  X + 1 , Y + 2 , Z );
    if ( val > max ) max = val;

    val = get(  X     , Y + 2 , Z + 1 );
    if ( val > max ) max = val;

    val = get(  X + 1 , Y + 2 , Z + 1 );
    if ( val > max ) max = val;

    }

    if ( Z > 0 ) {
    val = get(  X     , Y     , Z - 1 );
    if ( val > max ) max = val;

    val = get(  X + 1 , Y     , Z - 1 );
    if ( val > max ) max = val;

    val = get(  X     , Y + 1 , Z - 1 );
    if ( val > max ) max = val;

    val = get(  X + 1 , Y + 1 , Z - 1 );
    if ( val > max ) max = val;

    }

    if ( Z < getNumZ() - 2 ) {
    val = get( X     , Y     , Z + 2 );
    if ( val > max ) max = val;

    val = get( X + 1 , Y     , Z + 2 );
    if ( val > max ) max = val;

    val = get( X     , Y + 1 , Z + 2 );
    if ( val > max ) max = val;

    val = get( X + 1 , Y + 1 , Z + 2 );
    if ( val > max ) max = val;

    } */
  return max;
}

// needed for element saturation
template<typename _DataType>
_DataType qc::ScalarArray<_DataType, qc::QC_3D>::getElementSaturationMin ( int X, int Y, int Z ) const {
  DataType val, min = this->get ( X, Y, Z );

  // can this be implemented in the analogous way as above?

  val = this->get ( X + 1 , Y     , Z );
  if ( val < min ) min = val;

  val = this->get ( X     , Y + 1 , Z );
  if ( val < min ) min = val;

  val = this->get ( X + 1 , Y + 1 , Z );
  if ( val < min ) min = val;


  val = this->get ( X + 1 , Y     , Z + 1 );
  if ( val < min ) min = val;

  val = this->get ( X     , Y     , Z + 1 );
  if ( val < min ) min = val;

  val = this->get ( X + 1 , Y + 1 , Z + 1 );
  if ( val < min ) min = val;

  val = this->get ( X     , Y + 1 , Z + 1 );
  if ( val < min ) min = val;


  if ( X > 0 ) {
    val = this->get ( X - 1 , Y     , Z );
    if ( val < min ) min = val;

    val = this->get ( X - 1 , Y + 1 , Z );
    if ( val < min ) min = val;

    val = this->get ( X - 1 , Y     , Z + 1 );
    if ( val < min ) min = val;

    val = this->get ( X - 1 , Y + 1 , Z + 1 );
    if ( val < min ) min = val;

  }

  if ( X < this->getNumX() - 2 ) {
    val = this->get ( X + 2 , Y     , Z );
    if ( val < min ) min = val;

    val = this->get ( X + 2 , Y + 1 , Z );
    if ( val < min ) min = val;

    val = this->get ( X + 2 , Y     , Z + 1 );
    if ( val < min ) min = val;

    val = this->get ( X + 2 , Y + 1 , Z + 1 );
    if ( val < min ) min = val;

  }

  if ( Y > 0 ) {
    val = this->get ( X     , Y - 1 , Z );
    if ( val < min ) min = val;

    val = this->get ( X + 1 , Y - 1 , Z );
    if ( val < min ) min = val;

    val = this->get ( X     , Y - 1 , Z + 1 );
    if ( val < min ) min = val;

    val = this->get ( X + 1 , Y - 1 , Z + 1 );
    if ( val < min ) min = val;

  }

  if ( Y < this->getNumY() - 2 ) {
    val = this->get ( X     , Y + 2 , Z );
    if ( val < min ) min = val;

    val = this->get ( X + 1 , Y + 2 , Z );
    if ( val < min ) min = val;

    val = this->get ( X     , Y + 2 , Z + 1 );
    if ( val < min ) min = val;

    val = this->get ( X + 1 , Y + 2 , Z + 1 );
    if ( val < min ) min = val;

  }

  if ( Z > 0 ) {
    val = this->get ( X     , Y     , Z - 1 );
    if ( val < min ) min = val;

    val = this->get ( X + 1 , Y     , Z - 1 );
    if ( val < min ) min = val;

    val = this->get ( X     , Y + 1 , Z - 1 );
    if ( val < min ) min = val;

    val = this->get ( X + 1 , Y + 1 , Z - 1 );
    if ( val < min ) min = val;

  }

  if ( Z < this->getNumZ() - 2 ) {
    val = this->get ( X     , Y     , Z + 2 );
    if ( val < min ) min = val;

    val = this->get ( X + 1 , Y     , Z + 2 );
    if ( val < min ) min = val;

    val = this->get ( X     , Y + 1 , Z + 2 );
    if ( val < min ) min = val;

    val = this->get ( X + 1 , Y + 1 , Z + 2 );
    if ( val < min ) min = val;

  }
  return min;
}


template <typename _DataType>
void qc::ScalarArray<_DataType, qc::QC_3D>::writePovrayVolumeHeader ( std::ostream &Out ) const {
  const signed short int sX = this->getNumX(), sY = this->getNumY(), sZ = this->getNumZ();
  // cerr << "size is " << aol::intFormat ( sX ) << " "  << aol::intFormat ( sY ) << " "  << aol::intFormat ( sZ ) << endl;
  const signed char
    sXh = static_cast<signed char> ( sX >> 8 ),
    sYh = static_cast<signed char> ( sY >> 8 ),
    sZh = static_cast<signed char> ( sZ >> 8 ),
    sXl = static_cast<signed char> ( sX & 255 ),
    sYl = static_cast<signed char> ( sY & 255 ),
    sZl = static_cast<signed char> ( sZ & 255 );

  Out.write ( reinterpret_cast< const char* > ( &sXh ), 1 );   Out.write ( reinterpret_cast< const char* > ( &sXl ), 1 );
  Out.write ( reinterpret_cast< const char* > ( &sYh ), 1 );   Out.write ( reinterpret_cast< const char* > ( &sYl ), 1 );
  Out.write ( reinterpret_cast< const char* > ( &sZh ), 1 );   Out.write ( reinterpret_cast< const char* > ( &sZl ), 1 );
}

// template class qc::ScalarArray <bool, qc::QC_3D>; // qc::ScalarArray<bool, qc::QC_3D> should not be used ( += doesn't make sense with bool ), use qc::BitArray<qc::QC_3D> instead.
// template class qc::ScalarArray <char, qc::QC_3D>; // DO NOT USE. depending on your platform, char is unsigned or signed. use signed char or unsigned char.
template class qc::ScalarArray<signed char,    qc::QC_3D>;
template class qc::ScalarArray<unsigned char,  qc::QC_3D>;
template class qc::ScalarArray<short,          qc::QC_3D>;
template class qc::ScalarArray<unsigned short, qc::QC_3D>;
template class qc::ScalarArray<int,            qc::QC_3D>;
template class qc::ScalarArray<unsigned int,   qc::QC_3D>;
template class qc::ScalarArray<float,          qc::QC_3D>;
template class qc::ScalarArray<double,         qc::QC_3D>;
template class qc::ScalarArray<long double,    qc::QC_3D>;
