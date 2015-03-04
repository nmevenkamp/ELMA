#ifndef __PPMWRITER_H
#define __PPMWRITER_H

#include <scalarArray.h>

inline void set_colortrans_256 ( aol::Vector<float>* red, aol::Vector<float>* green, aol::Vector<float>* blue ) {
  aol::Vector<float> redv ( 5 ), greenv ( 5 ), bluev ( 5 );

  // hsv blue to red
  redv[ 0] = 0.0;
  greenv[ 0] = 0.0;
  bluev [ 0] = 1.0; // blue
  redv[ 1] = 0.0;
  greenv[ 1] = 1.0;
  bluev [ 1] = 1.0; // cyan
  redv[ 2] = 0.0;
  greenv[ 2] = 1.0;
  bluev [ 2] = 0.0; // green
  redv[ 3] = 1.0;
  greenv[ 3] = 1.0;
  bluev [ 3] = 0.0; // yellow
  redv[ 4] = 1.0;
  greenv[ 4] = 0.0;
  bluev [ 4] = 0.0; // red

  red->resize ( 256 );
  green->resize ( 256 );
  blue->resize ( 256 );

  for ( int i = 0; i < 256; i++ ) {
    red->set (   i, redv.interpolateInRange   ( i / 255.f, 0.0f, 1.0f ) );
    green->set ( i, greenv.interpolateInRange ( i / 255.f, 0.0f, 1.0f ) );
    blue->set (  i, bluev.interpolateInRange  ( i / 255.f, 0.0f, 1.0f ) );
  }
}

inline void set_morecolortrans_256 ( aol::Vector<float>* red, aol::Vector<float>* green, aol::Vector<float>* blue ) {
  aol::Vector<float> redv ( 6 ), greenv ( 6 ), bluev ( 6 );

  // hsv blue to red
  redv[ 5] = 1.0;
  greenv[ 5] = 0.0;
  bluev [ 5] = 0.0; // red
  redv[ 4] = 0.0;
  greenv[ 4] = 1.0;
  bluev [ 4] = 1.0; // cyan
  redv[ 3] = 0.0;
  greenv[ 3] = 0.0;
  bluev [ 3] = 1.0; // blue
  redv[ 2] = 1.0;
  greenv[ 2] = 1.0;
  bluev [ 2] = 0.0; // yellow
  redv[ 1] = 0.0;
  greenv[ 1] = 1.0;
  bluev [ 1] = 0.0; // green
  redv[ 0] = 1.0;
  greenv[ 0] = 0.0;
  bluev [ 0] = 1.0; // mangenta

  red->resize ( 256 );
  green->resize ( 256 );
  blue->resize ( 256 );

  for ( int i = 0; i < 256; i++ ) {
    red->set (   i, redv.interpolateInRange   ( i / 255.f, 0.0f, 1.0f ) );
    green->set ( i, greenv.interpolateInRange ( i / 255.f, 0.0f, 1.0f ) );
    blue->set (  i, bluev.interpolateInRange  ( i / 255.f, 0.0f, 1.0f ) );
  }
}

template <typename RealType>
inline void img_scaled_save_ctrscl ( const qc::ScalarArray<RealType, qc::QC_2D> &img, const char* filename_color, const RealType ctr, RealType scl, const aol::Vector<float>* red, const aol::Vector<float>* green, const aol::Vector<float>* blue ) {
  qc::ScalarArray<RealType, qc::QC_2D> tmpimg ( img, aol::DEEP_COPY );
  tmpimg.addToAll ( -ctr );

  if ( scl == 0.0 ) { // this means no explicit scaling. use full color width.
    const RealType imin = tmpimg.getMinValue(), imax = tmpimg.getMaxValue();
    scl = 127.5 / ( max ( imax, -imin ) );
#ifdef VERBOSE
    cerr << "saving image: min = " << imin << ", max = " << imax << ", scalefactor = " << scl << endl;
#endif
  } else {
    scl *= 127.5;
  };
  tmpimg *= scl;
  tmpimg.addToAll ( 127.5 );

  const int magicNumber = 6;

  ofstream outdat ( filename_color, ios::binary );
  outdat << "P" << magicNumber << "\n" << tmpimg.getNumX() << " " << tmpimg.getNumY() << "\n255\n";

  const int size = 3 * tmpimg.getNumX() * tmpimg.getNumY();
  const int numX = tmpimg.getNumX();
  unsigned char *dummyRGB = new unsigned char[size];

  int clippingNeeded = 0;
  for ( int j = 0; j < tmpimg.getNumX(); j++ ) {
    for ( int i = 0; i < tmpimg.getNumY(); i++ ) {
      int val = static_cast<int> ( tmpimg.get ( i, j ) );
      if ( val <   0 ) {
        clippingNeeded++;
        val =   0;
      };
      if ( val > 255 ) {
        clippingNeeded++;
        val = 255;
      };

      if ( magicNumber == 6 ) {
        const int index = 3 * ( i + j * numX );
        dummyRGB[index] = static_cast<unsigned char> ( 255.0 *   red->get ( val ) );
        dummyRGB[index+1] = static_cast<unsigned char> ( 255.0 * green->get ( val ) );
        dummyRGB[index+2] = static_cast<unsigned char> ( 255.0 *  blue->get ( val ) );
      }
      if ( magicNumber == 3 ) {
        const int rval = static_cast<int> ( 255.0 *   red->get ( val ) );
        const int gval = static_cast<int> ( 255.0 * green->get ( val ) );
        const int bval = static_cast<int> ( 255.0 *  blue->get ( val ) );

        outdat << rval << " " << gval << " " << bval << endl;
      }
    };
  };
  if ( magicNumber == 6 ) {
    outdat.write ( reinterpret_cast< char* > ( dummyRGB ), size );
  }
  delete[] dummyRGB;
  if ( clippingNeeded > 0 )
    cerr << clippingNeeded << " values had to been clipped into 0-255.\n";
  outdat.close();
}

template <typename RealType>
inline void img_scaled_save ( const qc::ScalarArray<RealType, qc::QC_2D> &img, const char* filename_color, const aol::Vector<float>* red, const aol::Vector<float>* green, const aol::Vector<float>* blue ) {
  img_scaled_save_ctrscl<RealType> ( img, filename_color, 0.5 * ( img.getMinValue() + img.getMaxValue() ), 0.0, red, green, blue );
}

template <typename RealType>
inline void img_scaled_save_scl ( const qc::ScalarArray<RealType, qc::QC_2D> &img, const char* filename_color, const RealType scl, const aol::Vector<float>* red, const aol::Vector<float>* green, const aol::Vector<float>* blue ) {
  img_scaled_save_ctrscl<RealType> ( img, filename_color, 0.5 * ( img.getMinValue() + img.getMaxValue() ), scl, red, green, blue );
}

template <typename RealType>
inline void img_scaled_save_ctr ( const qc::ScalarArray<RealType, qc::QC_2D> &img, const char* filename_color, const RealType ctr, const aol::Vector<float>* red, const aol::Vector<float>* green, const aol::Vector<float>* blue ) {
  img_scaled_save_ctrscl<RealType> ( img, filename_color, ctr, 0.0, red, green, blue );
}

#endif
