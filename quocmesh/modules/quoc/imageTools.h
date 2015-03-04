#ifndef __IMAGETOOLS_H
#define __IMAGETOOLS_H

#include <array.h>
#include <matrix.h>
#include <kernel3d.h>
#include <auxiliary.h>
#include <bzipiostream.h>
#include <multiArray.h>
#include <rectangularGrid.h>
#include <gnuplotter.h>
#include <parameterParser.h>

namespace qc {

enum ConcatType { LEFT_RIGHT, TOP_BOTTOM };

/** Class that allows the concatenation of ScalarArray<QC_2D>-objects for the purpose of output as a PGM image.
 *  It the ScalarArrays are of different sizes the images will be aligned top and right in the output image.
 *  Any background will be filled with the specified color.
 *  Please note that the recursive use of ConcatImage is possible. I.e. you can concatenate ConcatImages with
 *  this class as well.
 *  Please note further that this class is meaningful only for output as PGM image. In fact to save memory the
 *  implementation does some nasty things in terms of memory management.
 */

template<typename DataType, typename RealType = float>
class ConcatImage : public ScalarArray<DataType, qc::QC_2D> {
protected:
  const ScalarArray<DataType, qc::QC_2D> &_image1, &_image2;
  const ConcatType          _type;
  const int                 _spacing;
  const unsigned char       _bgColor;
public:
  ConcatImage ( ScalarArray<DataType, qc::QC_2D> &image1,
                ScalarArray<DataType, qc::QC_2D> &image2,
                const ConcatType type = LEFT_RIGHT,
                const int spacing = 1, const unsigned char bgColor = 255 ) :
      ScalarArray<DataType, qc::QC_2D> ( 1, 1 ),
      _image1 ( image1 ), _image2 ( image2 ), _type ( type ), _spacing ( spacing ), _bgColor ( bgColor ) {
    switch ( type ) {
    case LEFT_RIGHT :
      this->numX = image1.getNumX() + spacing + image2.getNumX();
      this->numY = aol::Max ( image1.getNumY(), image2.getNumY() );
      break;
    case TOP_BOTTOM:
      this->numX = aol::Max ( image1.getNumX(), image2.getNumX() );
      this->numY = image1.getNumY() + spacing + image2.getNumY();
      break;
    }
  }

  virtual bool createOverflowHandeledData ( unsigned char *tmp, DataType * /*min*/ = NULL, DataType * /*max*/ = NULL ) const {
    const int w1 = _image1.getNumX();
    const int w2 = _image2.getNumX();
    const int h1 = _image1.getNumY();
    const int h2 = _image2.getNumY();

    const int w = this->getNumX();
    const int h = this->getNumY();

    typedef unsigned char UCHAR;

    UCHAR *img1 = new UCHAR[w1*h1];
    UCHAR *img2 = new UCHAR[w2*h2];

    _image1.createOverflowHandeledData ( img1, NULL, NULL );
    _image2.createOverflowHandeledData ( img2, NULL, NULL );

    memset ( tmp, _bgColor, w*h );

    int x, y;

    switch ( _type ) {
    case LEFT_RIGHT: {
      UCHAR *ptr = tmp, *ptr1 = img1, *ptr2 = img2;
      const int c = aol::Min ( h1, h2 );
      const int d = aol::Max ( h1, h2 );
      for ( y = 0; y < c; ++y ) {
        for ( x = 0; x < w1; ++x ) *ptr++ = *ptr1++;
        ptr += _spacing;
        for ( x = 0; x < w2; ++x ) *ptr++ = *ptr2++;
      }
      if ( c == d ) break;
      if ( c == h1 ) {
        for ( y = c; y < d; ++y ) {
          ptr += ( w1 + _spacing );
          for ( x = 0; x < w2; ++x ) *ptr++ = *ptr2++;
        }
      } else {
        for ( y = c; y < d; ++y ) {
          for ( x = 0; x < w1; ++x ) *ptr++ = *ptr1++;
          ptr += ( w2 + _spacing );
        }
      }
    }
    break;
    case TOP_BOTTOM: {
      UCHAR *ptr = tmp, *ptr1 = img1, *ptr2 = img2;
      const int c = w - w1;
      const int d = w - w2;
      for ( y = 0; y < h1; ++y ) {
        for ( x = 0; x < w1; ++x ) *ptr++ = *ptr1++;
        ptr += c;
      }
      ptr += ( _spacing * w );
      for ( y = 0; y < h2; ++y ) {
        for ( x = 0; x < w2; ++x ) *ptr++ = *ptr2++;
        ptr += d;
      }
    }
    break;
    };

    delete[] img1;
    delete[] img2;

    return false;
  }

};

/** Class for the joint output of 1-3 ScalarArray<QC_2D> objects in one single PPM file. Thereby the different
 *  ScalarArray<QC_2D>-objects correspond to the colors red, green and blue in the output file.
 *  If a color is set to be NULL, it is considered zero in the output. At least one color must be non NULL.
 */

template<typename REAL>
class PPMImage {
public:
  typedef qc::ScalarArray<REAL, qc::QC_2D> VType;
protected:
  int _numX, _numY, size;
  VType *_red, *_green, *_blue;

public:
  PPMImage ( VType *red, VType *green = NULL, VType *blue = NULL ) : _red ( red ), _green ( green ), _blue ( blue ) {
    VType *ptr = NULL;
    if ( red ) ptr = red;
    else if ( green ) ptr = green;
    else if ( blue ) ptr = blue;
    if ( !ptr ) throw aol::Exception ( "Cannot create PPMImage with all colors set to NULL", __FILE__, __LINE__ );

    _numX = ptr->getNumX();
    _numY = ptr->getNumY();
    if ( ( red   && ( _numX != red->getNumX() || _numY != red->getNumY() ) ) ||
         ( green && ( _numX != green->getNumX() || _numY != green->getNumY() ) ) ||
         ( blue  && ( _numX != blue->getNumX() || _numY != blue->getNumY() ) ) ) throw aol::Exception ( "Dimensions of color components do not match", __FILE__, __LINE__ );

    size = _numX * _numY;
  }

  void save ( const char *fileName ) const {
    ofstream out ( fileName );
    save ( out );
  }

  void save ( ostream &out ) const {
    typedef unsigned char UCHAR;

    // Prepare the image data
    aol::Vector<UCHAR> red ( size ), green ( size ), blue( size );

    if ( _red ) _red->createOverflowHandledData ( red );
    if ( _green ) _green->createOverflowHandledData ( green );
    if ( _blue ) _blue->createOverflowHandledData ( blue );

    UCHAR *rgb = new UCHAR[size*3];
    UCHAR *ptr = rgb;
    for ( int i = 0; i < size; i++ ) {
      *ptr++ = red[i];
      *ptr++ = green[i];
      *ptr++ = blue[i];
    }

    // Write image header
    out << "P6" << endl << _numX << " " << _numY << endl << "255" << endl;
    out.write ( reinterpret_cast< char* > ( rgb ), size*3 );

    delete[] rgb;
  }

};

template <typename REAL>
std::ostream &operator<< ( std::ostream &out, PPMImage<REAL> &img ) {
  img.save ( out );
  return out;
}

template <typename RealType, typename GridType>
void write2dImage ( const GridType &Grid, const aol::Vector< RealType > &Image, const char* baseFileName, const bool /*WriteSlices*/ = 0, const Comp /*Direction*/ = qc::QC_X ) {
  if ( Grid.getDimOfWorld() == 2 ) {
    char fileName[1024];
    sprintf ( fileName, "%s.png", baseFileName );
    ScalarArray<RealType, qc::QC_2D> ImageArray ( Image, Grid.getNumX(), Grid.getNumY() );
    ImageArray.setQuietMode ( true );
    ImageArray.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0., 1. );
    ImageArray.savePNG ( fileName );
  }
  else
    throw aol::Exception ( "You can't execute write2dImage on a 3D grid!\n" );
}

template <typename RealType, typename GridType>
void write3dImage ( const GridType &Grid, const aol::Vector< RealType > &Image, const char* baseFileName, const bool WriteSlices = 0, const Comp Direction = qc::QC_X ) {
  if ( Grid.getDimOfWorld() == 3 ) {
    char fileName[1024];
    qc::ScalarArray<RealType, qc::QC_3D> *ImageArray = NULL;
    RealType minValue = Image.getMinValue();
    if ( ( minValue < 0 ) && WriteSlices ) {
      ImageArray = new qc::ScalarArray<RealType, qc::QC_3D> ( Image, Grid.getNumX(), Grid.getNumY(), Grid.getNumZ(), aol::DEEP_COPY );
      ImageArray->addToAll ( - minValue );
    } else {
      ImageArray = new qc::ScalarArray<RealType, qc::QC_3D> ( Image, Grid.getNumX(), Grid.getNumY(), Grid.getNumZ() );
    }
    ImageArray->setQuietMode ( true );
    if ( WriteSlices == 0 ) {
      sprintf ( fileName, "%s.dat.bz2", baseFileName );
      ImageArray->save ( fileName, PGM_DOUBLE_BINARY );
    } else {
      sprintf ( fileName, "%s_%%03d.pgm", baseFileName );
      ImageArray->saveSlices ( fileName, Direction, qc::PGM_UNSIGNED_CHAR_BINARY, NULL, aol::CLIP_THEN_SCALE, 0, 1 );
    }
    delete ImageArray;
  }
  else
    throw aol::Exception ( "You can't execute write3dImage on a 2D grid!\n" );
}

template <typename RealType, typename GridType>
void writeImage ( const GridType &Grid, const aol::Vector< RealType > &Image, const char* baseFileName, const bool WriteSlices = 0, const Comp Direction = qc::QC_X ) {
  if ( Grid.getDimOfWorld() == 2 ) {
    write2dImage<RealType, GridType> ( Grid, Image, baseFileName, WriteSlices, Direction );
  }
  if ( Grid.getDimOfWorld() == 3 ) {
    write3dImage<RealType, GridType> ( Grid, Image, baseFileName, WriteSlices, Direction );
  }
}

template <typename RealType, typename GridType>
void writeImage ( const GridType &Grid, const aol::MultiVector< RealType > &MImage, const char* baseFileName, const bool WriteSlices = 0, const Comp Direction = qc::QC_X ) {
  char filename[1024];
  for ( int i = 0; i < MImage.numComponents(); i++ ) {
    sprintf ( filename, "%s_%d", baseFileName, i );
    qc::writeImage<RealType, GridType> ( Grid, MImage[i], filename, WriteSlices, Direction );
  }
}

template <typename RealType, typename GridType>
void writeImage ( const GridType &Grid, const qc::MultiArray< RealType, 2, 1 > &MImage, const char* baseFileName, const bool WriteSlices = 0, const Comp Direction = qc::QC_X ) {
  qc::writeImage<RealType, GridType> ( Grid, MImage[0], baseFileName, WriteSlices, Direction );
}

template <typename RealType, typename GridType>
void writeImage ( const GridType &/*Grid*/, const qc::MultiArray< RealType, 2, 3 > &MImage, const char* baseFileName, const bool /*WriteSlices*/ = 0, const Comp /*Direction*/ = qc::QC_X ) {
  qc::MultiArray< RealType, 2, 3 > image ( MImage, aol::FLAT_COPY );
  image.setOverflowHandling ( aol::CLIP_THEN_SCALE, 0., 1. );
  image.savePNG ( aol::strprintf ( "%s.png", baseFileName ).c_str() );
}

/**
 * \author Berkels
 */
template <typename RealType>
void plotLineOf2DArray ( const ScalarArray<RealType, qc::QC_2D> &Array, const int LineNumber, const char *BaseName, const aol::PlotOutFileType OutType ){
  const int numX = Array.getNumX();
  aol::Vector<RealType> slice ( numX );
  for ( int i = 0; i < numX; ++i )
    slice[i] = Array.get ( i, LineNumber );

  aol::Plotter<RealType> plotter;
  plotter.set_outfile_base_name( BaseName );
  aol::PlotDataFileHandler<RealType> plotHandler;
  plotHandler.generateFunctionPlot( slice );
  plotter.addPlotCommandsFromHandler( plotHandler );
  plotter.genPlot( OutType );
}

/**
 * Maps the scalar values of a ScalarArray from [0,1] to color values in a MultiArray
 * using the standard HSV hue scale.
 *
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
void convertScalarToColor ( const qc::ScalarArray<RealType, Dim> &Scalar, qc::MultiArray<unsigned char, Dim, 3> &BufArray ) {

  RealType rgb[3] = {0, 0, 0}, hsv[3];

  for ( int i = 0; i < Scalar.size(); ++i ) {
    RealType b = Scalar[i];

    while ( b > 1. ) b -= 1.;
    while ( b < 0. ) b += 1.;

    // A weighting to alter the brightness could be added here.
    const RealType w = 1;

    hsv[0] = b;
    hsv[1] = w;
    hsv[2] = w;

    aol::RGBColorMap<RealType>::hsv2rgb ( hsv, rgb );

    BufArray[0][i] = static_cast< unsigned char> ( rgb[0] * 255 );
    BufArray[1][i] = static_cast< unsigned char> ( rgb[1] * 255 );
    BufArray[2][i] = static_cast< unsigned char> ( rgb[2] * 255 );
  }
}

/**
 * Algorithm for bicubic interpolation.
 *
 * \author Paul Breeuwsma
 *
 * Converted from Java by Berkels.
 * Original Java source available at http://www.paulinternet.nl/?page=bicubic
 */
template <typename RealType>
class CachedBicubicInterpolator {
  aol::Mat<4,4, RealType> _a;

public:
  void updateCoefficients ( const aol::Mat<4,4, RealType> &P ) {
    _a[0][0] = P[1][1];
    _a[0][1] = -.5*P[1][0] + .5*P[1][2];
    _a[0][2] = P[1][0] - 2.5*P[1][1] + 2*P[1][2] - .5*P[1][3];
    _a[0][3] = -.5*P[1][0] + 1.5*P[1][1] - 1.5*P[1][2] + .5*P[1][3];
    _a[1][0] = -.5*P[0][1] + .5*P[2][1];
    _a[1][1] = .25*P[0][0] - .25*P[0][2] - .25*P[2][0] + .25*P[2][2];
    _a[1][2] = -.5*P[0][0] + 1.25*P[0][1] - P[0][2] + .25*P[0][3] + .5*P[2][0] - 1.25*P[2][1] + P[2][2] - .25*P[2][3];
    _a[1][3] = .25*P[0][0] - .75*P[0][1] + .75*P[0][2] - .25*P[0][3] - .25*P[2][0] + .75*P[2][1] - .75*P[2][2] + .25*P[2][3];
    _a[2][0] = P[0][1] - 2.5*P[1][1] + 2*P[2][1] - .5*P[3][1];
    _a[2][1] = -.5*P[0][0] + .5*P[0][2] + 1.25*P[1][0] - 1.25*P[1][2] - P[2][0] + P[2][2] + .25*P[3][0] - .25*P[3][2];
    _a[2][2] = P[0][0] - 2.5*P[0][1] + 2*P[0][2] - .5*P[0][3] - 2.5*P[1][0] + 6.25*P[1][1] - 5*P[1][2] + 1.25*P[1][3] + 2*P[2][0] - 5*P[2][1] + 4*P[2][2] - P[2][3] - .5*P[3][0] + 1.25*P[3][1] - P[3][2] + .25*P[3][3];
    _a[2][3] = -.5*P[0][0] + 1.5*P[0][1] - 1.5*P[0][2] + .5*P[0][3] + 1.25*P[1][0] - 3.75*P[1][1] + 3.75*P[1][2] - 1.25*P[1][3] - P[2][0] + 3*P[2][1] - 3*P[2][2] + P[2][3] + .25*P[3][0] - .75*P[3][1] + .75*P[3][2] - .25*P[3][3];
    _a[3][0] = -.5*P[0][1] + 1.5*P[1][1] - 1.5*P[2][1] + .5*P[3][1];
    _a[3][1] = .25*P[0][0] - .25*P[0][2] - .75*P[1][0] + .75*P[1][2] + .75*P[2][0] - .75*P[2][2] - .25*P[3][0] + .25*P[3][2];
    _a[3][2] = -.5*P[0][0] + 1.25*P[0][1] - P[0][2] + .25*P[0][3] + 1.5*P[1][0] - 3.75*P[1][1] + 3*P[1][2] - .75*P[1][3] - 1.5*P[2][0] + 3.75*P[2][1] - 3*P[2][2] + .75*P[2][3] + .5*P[3][0] - 1.25*P[3][1] + P[3][2] - .25*P[3][3];
    _a[3][3] = .25*P[0][0] - .75*P[0][1] + .75*P[0][2] - .25*P[0][3] - .75*P[1][0] + 2.25*P[1][1] - 2.25*P[1][2] + .75*P[1][3] + .75*P[2][0] - 2.25*P[2][1] + 2.25*P[2][2] - .75*P[2][3] - .25*P[3][0] + .75*P[3][1] - .75*P[3][2] + .25*P[3][3];
  }

  RealType getValue ( const RealType X, const RealType Y ) const {
    RealType x2 = aol::Sqr( X );
    RealType x3 = x2 * X;
    RealType y2 = aol::Sqr( Y );
    RealType y3 = y2 * Y;

    return ( _a[0][0] + _a[0][1] * Y + _a[0][2] * y2 + _a[0][3] * y3 ) +
           ( _a[1][0] + _a[1][1] * Y + _a[1][2] * y2 + _a[1][3] * y3 ) * X +
           ( _a[2][0] + _a[2][1] * Y + _a[2][2] * y2 + _a[2][3] * y3 ) * x2 +
           ( _a[3][0] + _a[3][1] * Y + _a[3][2] * y2 + _a[3][3] * y3 ) * x3;
  }
};

/**
 * ConnectedComponentsLabeler is an algorithm that applies the Connected Component Labeling
 * alogrithm to an input qc::BitArray<qc::QC_2D>.
 *
 * \author Neil Brown, DAI
 * \author Judy Robertson, SELLIC OnLine
 *
 * Converted from Java by Berkels, cf.
 * http://homepages.inf.ed.ac.uk/rbf/HIPR2/label.htm
 * http://homepages.inf.ed.ac.uk/rbf/HIPR2/flatjavasrc/ImageLabel.java
 */
class ConnectedComponentsLabeler {
private:
  /**
   * Associate(equivalence) A with B.
   * A should be less than B to give some ordering (sorting)
   * if B is already associated with some other value, then propagate
   * down the list.
    */
  static void associate( aol::Vector<unsigned int> &Labels, unsigned int A, unsigned int B ) {
    if( A > B ) {
      associate( Labels, B, A );
      return;
    }
    if( ( A == B ) || ( Labels[ B ] == A ) ) return;
    if( Labels[ B ] == B ) {
      Labels[ B ] = A;
    } else {
      associate( Labels, Labels[ B ], A );
      if (Labels[ B ] > A) {             //***rbf new
        Labels[ B ] = A;
      }
    }
  }
  /**
   * Reduces the number of Labels.
   */
  static unsigned int reduce( aol::Vector<unsigned int> &Labels, unsigned int A ){
    if( Labels[ A ] == A ){
      return A;
    } else {
      return reduce( Labels, Labels[ A ] );
    }
  }
public:

  /**
   * doLabel applies the Labeling alogrithm
   *
   * \param Mask The input pixel array
   * \param LabelArray A pixel array containing the labelled image
   *
   * Returns the number of labels used.
   *
   * NB For images  0,0 is the top left corner.
   */
  static int doLabel ( const qc::BitArray<qc::QC_2D> &Mask, qc::ScalarArray<int, qc::QC_2D> &LabelArray ) {

    unsigned int nextlabel = 1;
    aol::Vec<4, bool> nbs;
    aol::Vec<4, unsigned int> nbls;
    unsigned int result = 0;
    const int numX = Mask.getNumX();
    const int numY = Mask.getNumY();
    // the most labels there can be is 1/2 of the points in checkerboard
    aol::Vector<unsigned int> labels ( numX * numY / 2 );
    //initialise labels
    for (int i=0; i<labels.size(); i++) labels[ i ] = i;

    //now Label the image
    for ( int y = 0; y < numY; ++y ) {
      for ( int x = 0; x < numX; ++x ) {
        if ( Mask.get( x, y ) == false ) {
          result = 0;  //nothing here
        } else {

          //The 4 visited neighbours
          nbs[ 0 ] = ( x > 0 ) && Mask.get ( x-1, y );
          nbs[ 1 ] = ( y > 0 ) && Mask.get ( x, y-1 );
          nbs[ 2 ] = ( x > 0 ) && ( y > 0 ) && Mask.get ( x-1, y-1 );
          nbs[ 3 ] = ( x < numX - 1 ) && ( y > 0 ) && Mask.get ( x+1, y-1 );

          //Their corresponding labels
          nbls[ 0 ] = ( x > 0 ) ? LabelArray.get ( x-1, y ) : 0;
          nbls[ 1 ] = ( y > 0 ) ? LabelArray.get ( x, y-1 ) : 0;
          nbls[ 2 ] = ( ( x > 0 ) && ( y > 0 ) ) ? LabelArray.get ( x-1, y-1 ) : 0;
          nbls[ 3 ] = ( ( x < numX - 1 ) && ( y > 0 ) ) ? LabelArray.get ( x+1, y-1 ) : 0;

          //label the point
          if( (nbs[0] == nbs[1]) && (nbs[1] == nbs[2]) && (nbs[2] == nbs[3])
            && (nbs[0] == 0 )) {
              // all neighbours are 0 so gives this point a new label
              result = nextlabel;
              nextlabel++;
          } else { //one or more neighbours have already got labels
            int count = 0;
            int found = -1;
            for( int j=0; j<4; j++){
              if( nbs[ j ] != 0 ){
                count +=1;
                found = j;
              }
            }
            if( count == 1 ) {
              // only one neighbour has a label, so assign the same label to this.
              result = nbls[ found ];
            } else {
              // more than 1 neighbour has a label
              result = nbls[ found ];
              // Equivalence the connected points
              for( int j=0; j<4; j++){
                if( ( nbls[ j ] != 0 ) && (nbls[ j ] != result ) ){
                  associate( labels, nbls[ j ], result );
                }
              }
            }
          }
        }
        LabelArray.set ( x, y, result );
      }
    }
    //reduce labels ie 76=23=22=3 -> 76=3
    //done in reverse order to preserve sorting
    for( unsigned int i= labels.size()-1; i > 0; i-- ){
      labels[ i ] = reduce( labels, i );
    }

    /*now labels will look something like 1=1 2=2 3=2 4=2 5=5.. 76=5 77=5
      this needs to be condensed down again, so that there is no wasted
      space eg in the above, the labels 3 and 4 are not used instead it jumps
      to 5.
      */
    aol::Vector<unsigned int> condensed ( nextlabel ); // cant be more than nextlabel labels

    unsigned int count = 0;
    for (unsigned int i=0; i< nextlabel; i++){
      if( i == labels[ i ] ) condensed[ i ] = count++;
    }
    // Record the number of labels
    const int numberOfLabels = count - 1;

    // now run back through our preliminary results, replacing the raw label
    // with the reduced and condensed one, and do the scaling and offsets too
    for (int i=0; i< LabelArray.size(); i++)
      LabelArray[i] = condensed[ labels[ LabelArray[ i ] ] ];

    return numberOfLabels;
  }

  /**
   * Extension of doLabel to 3D. Can also handle periodic domains.
   * \author Paul Springer (AICES)
   * \param[out] labelSize Stores the number of pixels belonging to each connected component. Label-size of label i is stored at position i (i.e., the last label-size is stored at labelSize[numberOfLabels]).
   * \param connectivity 0: 6-way connectivity (i.e. left, right, top, bottom, front, back) or 1: full connectivity(all 26 neighbors (i.e. with diagonal elements))
   * \return numberOfLabels Number of labels found (excluding label 0)
   */
  template<int PBC_X, int PBC_Y, int PBC_Z, int connectivity>
  static int doLabel3D ( const qc::BitArray<qc::QC_3D> &Mask, qc::ScalarArray<int, qc::QC_3D> &LabelArray, aol::Vector<unsigned int> &labelSize ) {

    unsigned int nextlabel = 1;
    unsigned int result = 0;
    const int numX = Mask.getNumX();
    const int numY = Mask.getNumY();
    const int numZ = Mask.getNumZ();
    // the most labels there can be is 1/2 of the points in checkerboard
    aol::Vector<unsigned int> labels ( numX * numY * numZ / 2 );
    //initialise labels
    for (int i=0; i<labels.size(); i++) labels[ i ] = i;

    //now Label the image
    for ( int z = 0; z < numZ; ++z ) {
      for ( int y = 0; y < numY; ++y ) {
        for ( int x = 0; x < numX; ++x ) {

          if ( Mask.get( x, y, z ) == false ) {
            result = 0;  //nothing here
          } else if ( LabelArray.get( x, y, z ) == 0 ) {

            int found = 0;
            int myLabel = 0;
            //loop over neighbors
            //1) look for a neighbor which already has a label
            //2) associate myLable with all my neihgbors
            for(int xx = -1; xx <= 1; ++xx){
              int tmpX;
              if(PBC_X){
                tmpX = (x + xx+numX)%numX;
              }else {
                if( x+xx >= numX ||  x+xx < 0) continue;
                tmpX = x + xx;
              }

              for(int yy = -1; yy <= 1; ++yy){

                //if yy*yy + xx*xx > 1 we are looking at diagonal elements or corners, so we skip these.
                if(connectivity == 0 && xx*xx + yy*yy > 1)
                  continue;

                int tmpY;
                if(PBC_Y){
                  tmpY = (y + yy+numY)%numY;
                }else{
                  if( y+yy >= numY ||  y+yy < 0) continue;
                  tmpY = y+yy;
                }

                for(int zz = -1; zz <= 1; ++zz){
                  int tmpZ;
                  if(PBC_Z){
                    tmpZ = (z + zz+numZ)%numZ;
                  }else{
                    if( z+zz >= numZ ||  z+zz < 0) continue;
                    tmpZ = z+zz;
                  }

                  //if zz*zz + yy*yy + xx*xx > 1 we are looking at diagonal elements or corners, so we skip these.
                  if(connectivity == 0 && zz*zz + yy*yy + xx*xx > 1)
                    continue;

                  if(Mask.get(tmpX, tmpY, tmpZ)){
                    int nbrLabel = LabelArray.get ( tmpX, tmpY, tmpZ );
                    if(nbrLabel){
                      if(found == 0){
                        found = 1;
                        myLabel = nbrLabel;
                        LabelArray.set ( x, y, z, myLabel);
                      }else{
                        associate( labels, nbrLabel, myLabel);
                      }
                    }
                  }
                }
              }
            }

            //none of my neighbors has a label yet
            if ( found == 0 ){
              myLabel = nextlabel;
              LabelArray.set ( x, y, z, myLabel);
              nextlabel++;
            }

            //loop over neighbors
            //assign same label to all neighbors which don't have a label yet
            for(int xx = -1; xx <= 1; ++xx){
              int tmpX;
              if(PBC_X){
                tmpX = (x + xx+numX)%numX;
              }else {
                if( x+xx >= numX ||  x+xx < 0) continue;
                tmpX = x + xx;
              }

              for(int yy = -1; yy <= 1; ++yy){

                //if yy*yy + xx*xx > 1 we are looking at diagonal elements or corners, so we skip these.
                if(connectivity == 0 && xx*xx + yy*yy > 1)
                  continue;

                int tmpY;
                if(PBC_Y){
                  tmpY = (y + yy+numY)%numY;
                }else{
                  if( y+yy >= numY ||  y+yy < 0) continue;
                  tmpY = y+yy;
                }

                for(int zz = -1; zz <= 1; ++zz){
                  int tmpZ;
                  if(PBC_Z){
                    tmpZ = (z + zz+numZ)%numZ;
                  }else{
                    if( z+zz >= numZ ||  z+zz < 0) continue;
                    tmpZ = z+zz;
                  }

                  //if zz*zz + yy*yy + xx*xx > 1 we are looking at diagonal elements or corners, so we skip these.
                  if(connectivity == 0 && zz*zz + yy*yy + xx*xx > 1)
                    continue;

                  if(Mask.get(tmpX, tmpY, tmpZ) && LabelArray.get ( tmpX, tmpY, tmpZ ) == 0){
                    LabelArray.set ( tmpX, tmpY, tmpZ, myLabel );
                  }
                }
              }
            }

          }
        }
      }
    }

    //reduce labels ie 76=23=22=3 -> 76=3
    //done in reverse order to preserve sorting
    for( unsigned int i= labels.size()-1; i > 0; i-- ){
      labels[ i ] = reduce( labels, i );
    }

    /*now labels will look something like 1=1 2=2 3=2 4=2 5=5.. 76=5 77=5
      this needs to be condensed down again, so that there is no wasted
      space eg in the above, the labels 3 and 4 are not used instead it jumps
      to 5.
      */
    aol::Vector<unsigned int> condensed ( nextlabel ); // cant be more than nextlabel labels

    unsigned int count = 0;
    for (unsigned int i=0; i< nextlabel; i++){
      if( i == labels[ i ] ) condensed[ i ] = count++;
    }
    // Record the number of labels
    const int numberOfLabels = count - 1;

    // now run back through our preliminary results, replacing the raw label
    // with the reduced and condensed one, and do the scaling and offsets too
    for (int i=0; i< LabelArray.size(); i++){
      LabelArray[i] = condensed[ labels[ LabelArray[ i ] ] ];
      labelSize[LabelArray[i]]++;
    }

    return numberOfLabels;
  }
};


/**
 * Very simple volume rendering using orthogonal projection and a basic emission / absorption model.
 *
 * \author Berkels
 */
template <typename RealType>
void renderVolume ( const qc::ScalarArray<RealType, qc::QC_3D> &Volume, const qc::Comp Direction, qc::ScalarArray<RealType, qc::QC_2D> &Projection ){
  int tmpNumX = 0, tmpNumY = 0, numOfSlices = 0;

  switch ( Direction ) {
    case qc::QC_X:
      tmpNumX = Volume.getNumY();
      tmpNumY = Volume.getNumZ();
      numOfSlices = Volume.getNumX();
      break;
    case qc::QC_Y:
      tmpNumX = Volume.getNumX();
      tmpNumY = Volume.getNumZ();
      numOfSlices = Volume.getNumY();
      break;
    case qc::QC_Z:
      tmpNumX = Volume.getNumX();
      tmpNumY = Volume.getNumY();
      numOfSlices = Volume.getNumZ();
      break;
    default:
      throw aol::UnimplementedCodeException ( "Unknown qc::Comp", __FILE__, __LINE__ );
  }

  Projection.reallocate ( tmpNumX, tmpNumY );

  // The following may look like an awful example of code duplication, but
  // the order of the loops is important for the execution speed. Furthermore,
  // it can't hurt not to have the switch in the inner most loop.
  switch ( Direction ) {
    case qc::QC_X:
      for ( int y = 0; y < tmpNumY; ++y ) {
        for ( int x = 0; x < tmpNumX; ++x ) {
          for ( int z = 0; z <numOfSlices; ++z ) {
            const RealType value = Volume.get ( z, x, y );
            const RealType opacity = value;
            Projection.set ( x, y, opacity * value + ( 1 - opacity ) * Projection.get ( x, y ) );
          }
        }
      }
      break;
    case qc::QC_Y:
      for ( int y = 0; y < tmpNumY; ++y ) {
        for ( int z = 0; z <numOfSlices; ++z ) {
          for ( int x = 0; x < tmpNumX; ++x ) {
            const RealType value = Volume.get ( x, z, y );
            const RealType opacity = value;
            Projection.set ( x, y, opacity * value + ( 1 - opacity ) * Projection.get ( x, y ) );
          }
        }
      }
      break;
    case qc::QC_Z:
      for ( int z = 0; z <numOfSlices; ++z ) {
        for ( int y = 0; y < tmpNumY; ++y ) {
          for ( int x = 0; x < tmpNumX; ++x ) {
            const RealType value = Volume.get ( x, y, z );
            const RealType opacity = value;
            Projection.set ( x, y, opacity * value + ( 1 - opacity ) * Projection.get ( x, y ) );
          }
        }
      }
      break;
    default:
      throw aol::UnimplementedCodeException ( "Unknown qc::Comp", __FILE__, __LINE__ );
  }
}

/**
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
void resampleArrayToQuadraticArrayKeepingAspectRatio ( const qc::ScalarArray<RealType, Dim> &InputArray, qc::ScalarArray<RealType, Dim> &OutputArray, const RealType FillValue = 0 ) {
  // OutputArray needs to be quadatic, otherwise the code below won't keep the aspect ratio.
  GridSize<Dim> sizeChecker ( OutputArray );
  sizeChecker.quadraticOrDie ();

  GridSize<Dim> paddedSize ( static_cast<short> ( aol::Max ( InputArray.getNumX(), InputArray.getNumY(), InputArray.getNumZ() ) ) );
  qc::ScalarArray<RealType, Dim> paddedArray ( paddedSize );
  paddedArray.padFrom ( InputArray, FillValue );
  OutputArray.resampleFrom ( paddedArray );
}

/**
 * \author Berkels
 */
template <typename RealType>
void joinTwo2DArraysVertically( const qc::ScalarArray<RealType, qc::QC_2D> &Arg1, const qc::ScalarArray<RealType, qc::QC_2D> &Arg2, qc::ScalarArray<RealType, qc::QC_2D> &Dest, const int Spacing = 0 ) {
  if( Arg1.getNumX() != Arg2.getNumX() )
    throw aol::Exception( "Arg1.getNumX() != Arg2.getNumX() !", __FILE__, __LINE__);

  Dest.reallocate( Arg1.getNumX(), Arg1.getNumY() + Spacing + Arg2.getNumY() );
  Dest.pasteFrom ( Arg1, 0, 0 );
  Dest.pasteFrom ( Arg2, 0, Arg1.getNumY() + Spacing );

  for( int j = 0; j < Spacing; j++ ) {
    for( int i = 0; i < Arg1.getNumX(); i++ ) {
      Dest.set( i, j + Arg1.getNumY(), 0.5 );
    }
  }
}

/**
 * \author Berkels
 */
template <typename RealType>
void shrinkAndPad2DArray ( qc::ScalarArray<RealType, qc::QC_2D> &Array, const RealType ScaleFactor, const RealType FillValue = 0 ) {
  qc::ScalarArray<RealType, qc::QC_2D> arrayScaled ( static_cast<int> ( Array.getNumX() * ScaleFactor ), static_cast<int> ( Array.getNumY() * ScaleFactor ) );
  arrayScaled.resampleFrom ( Array );
  Array.padFrom ( arrayScaled, FillValue );
}

/**
 * \author Berkels
 */
template <typename RealType>
void crop2DImage ( qc::ScalarArray<RealType, qc::QC_2D> &Image, const aol::ParameterParser &Parser ) {
  if ( Parser.checkAndGetBool ( "cropInput" ) ) {
    qc::ScalarArray<RealType, qc::QC_2D> imageCropped ( Parser.getInt ( "cropSizeX" ), Parser.getInt ( "cropSizeY" ) );
    Image.copyBlockTo ( Parser.getInt ( "cropStartX" ), Parser.getInt ( "cropStartY" ) , imageCropped );
    Image.reallocate ( imageCropped );
    Image = imageCropped;
  }
}

} // end of namespace qc.

#endif
