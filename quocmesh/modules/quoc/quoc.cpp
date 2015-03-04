#include <quoc.h>

namespace qc {

int getSizeOfSaveType ( const qc::SaveType Type ) {
  switch ( Type ) {
  case PGM_UNSIGNED_CHAR_ASCII:
  case PGM_UNSIGNED_CHAR_BINARY:
    return 1;
    break;
  case PGM_UNSIGNED_SHORT_BINARY:
  case PGM_SHORT_BINARY:
    return 2;
    break;
  case PGM_FLOAT_ASCII:
  case PGM_FLOAT_BINARY:
  case PGM_UNSIGNED_INT_BINARY:
  case PGM_SIGNED_INT_BINARY:
    return 4;
    break;
  case PGM_DOUBLE_BINARY:
    return 8;
    break;
  default:
    throw aol::Exception ( "qc::getSizeOfSaveType: Unsupported SaveType", __FILE__, __LINE__ );
    break;
  }
}

const char *getDefaulSuffixOfSaveType ( const qc::SaveType Type ) {
  switch ( Type ) {
    case PGM_UNSIGNED_CHAR_ASCII:
    case PGM_UNSIGNED_CHAR_BINARY:
      return ".pgm";
      break;
    case PGM_UNSIGNED_SHORT_BINARY:
    case PGM_SHORT_BINARY:
    case PGM_FLOAT_ASCII:
    case PGM_FLOAT_BINARY:
    case PGM_UNSIGNED_INT_BINARY:
    case PGM_SIGNED_INT_BINARY:
    case PGM_DOUBLE_BINARY:
      return ".dat.bz2";
      break;
    case PNG_2D:
      return ".png";
      break;
    default:
      throw aol::Exception ( "qc::getDefaulSuffixOfSaveType: Unsupported SaveType", __FILE__, __LINE__ );
      break;
  }
}

template<> const qc::SaveType qc::SaveTypeTrait< unsigned char  >::AsciiSaveType = qc::PGM_UNSIGNED_CHAR_ASCII;
template<> const qc::SaveType qc::SaveTypeTrait< float          >::AsciiSaveType = qc::PGM_FLOAT_ASCII;

template<> const qc::SaveType qc::SaveTypeTrait< unsigned char  >::BinarySaveType = qc::PGM_UNSIGNED_CHAR_BINARY;
template<> const qc::SaveType qc::SaveTypeTrait< short          >::BinarySaveType = qc::PGM_SHORT_BINARY;
template<> const qc::SaveType qc::SaveTypeTrait< unsigned short >::BinarySaveType = qc::PGM_UNSIGNED_SHORT_BINARY;
template<> const qc::SaveType qc::SaveTypeTrait< float          >::BinarySaveType = qc::PGM_FLOAT_BINARY;
template<> const qc::SaveType qc::SaveTypeTrait< double         >::BinarySaveType = qc::PGM_DOUBLE_BINARY;

}
