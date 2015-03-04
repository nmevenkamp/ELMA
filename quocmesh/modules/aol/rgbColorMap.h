#ifndef __RGBCOLORMAP_H
#define __RGBCOLORMAP_H

#include <aol.h>
#include <multiVector.h>
#include <smallVec.h>

namespace aol {

enum colorTrans { BLACK_TO_WHITE,
                  WHITE_TO_BLACK,
                  BLUE_TO_RED,
                  HSV_BLUE_TO_RED,
                  HSV_BLUE_TO_MAGENTA,
                  HSV_RED_TO_BLUE,
                  TOPOGRAPHIC,
                  WHITE_TO_GREEN_TO_BLUE_TO_RED,
                  WHITE_TO_BLUE,
                  UNDEFINED };

/** Routines associated with color mapping
 *
 * <ul><li>Mapping from scalar values in range [Min, Max] to three-channel color with predefined color maps</li>
 *     <li>rgb to hsv</li>
 *     <li>hsv to rgb</li>
 *     <li>rgb color string</li></ul>
 *
 *  \author Schwen
 */
template< typename RealType >
class RGBColorMap {
protected:
  RealType _min, _max;
  aol::MultiVector<RealType> _colorValues;

public:
  RGBColorMap ( ) : _min ( aol::NumberTrait<RealType>::NaN ), _max ( aol::NumberTrait<RealType>::NaN ), _colorValues() {
  }

  // copy constructor and assignment operator should do the correct thing

  RGBColorMap ( const RealType Min, const RealType Max, const colorTrans ColorTrans = BLACK_TO_WHITE );

  //! map scalar value to color triple, returning color value (slow due to temporary object)
  inline aol::Vec3<RealType> scalarToColor ( const RealType val ) const {
    aol::Vec3<RealType> col;
    scalarToColor ( val, col );
    return ( col );
  }

  //! map scalar value to color triple
  inline void scalarToColor ( const RealType val, aol::Vec3<RealType> &color ) const {
    for ( short i = 0; i < 3; ++i )
      color[i] = _colorValues[i].interpolateInRange ( val, _min, _max );
  }

  //! Simple routine which converts from a HSV triple to an RGB triple. Both have to be in the range 0..1,0..1,0..1
  inline static void hsv2rgb ( const RealType* hsv, RealType* rgb )  {
    const int hueCase = static_cast< int > ( hsv[0] * 6 );
    const RealType frac = 6 * hsv[0] - hueCase;
    const RealType lx  = hsv[2] * ( 1.0 - hsv[1] );
    const RealType ly  = hsv[2] * ( 1.0 - hsv[1] * frac );
    const RealType lz  = hsv[2] * ( 1.0 - hsv[1] * ( 1.0 - frac ) );

    switch ( hueCase ) {
    case 0:
    case 6:
      rgb[0] = hsv[2];
      rgb[1] = lz;
      rgb[2] = lx;
      break;  /* 0<hue<1/6   */
    case 1:
      rgb[0] = ly;
      rgb[1] = hsv[2];
      rgb[2] = lx;
      break;  /* 1/6<hue<2/6 */
    case 2:
      rgb[0] = lx;
      rgb[1] = hsv[2];
      rgb[2] = lz;
      break;  /* 2/6<hue<3/6 */
    case 3:
      rgb[0] = lx;
      rgb[1] = ly;
      rgb[2] = hsv[2];
      break;  /* 3/6<hue/4/6 */
    case 4:
      rgb[0] = lz;
      rgb[1] = lx;
      rgb[2] = hsv[2];
      break;  /* 4/6<hue<5/6 */
    case 5:
      rgb[0] = hsv[2];
      rgb[1] = lx;
      rgb[2] = ly;
      break;  /* 5/6<hue<1   */
    default:
      throw aol::UnimplementedCodeException ( "hsv2rgb: unhandled switch case!", __FILE__, __LINE__ );
    }
  }

  //! Simple routine which converts from an RGB triple to an HSV triple. Both have to be in the range 0..1,0..1,0..1
  inline static void rgb2hsv ( const RealType* rgb, RealType* hsv )  {
    const RealType M = aol::Max<RealType> ( rgb[0], rgb[1], rgb[2] );
    const RealType m = aol::Min<RealType> ( rgb[0], rgb[1], rgb[2] );
    const RealType d = M - m;

    hsv[2] = M;    //value == max(r,g,b)
    hsv[1] = ( M > 0.00001 ) ? d / M : 0; //saturation

    if ( hsv[1] == 0 )
      hsv[0] = 0;  //achromatic case, hue is 0 by convention
    else     //chromatic case
    {
      if ( rgb[0] == M )
        hsv[0] = ( rgb[1] - rgb[2] ) / d;
      else if ( rgb[1] == M )
        hsv[0] = 2 + ( rgb[2] - rgb[0] ) / d;
      else
        hsv[0] = 4 + ( rgb[0] - rgb[1] ) / d;

      hsv[0] /= 6;
      if ( hsv[0] < 0 )
        hsv[0] += 1;
    }
  }

  //! \brief Returns the color values.
  inline const aol::MultiVector<RealType>& getColorValues () const {
    return _colorValues;
  }

  inline const RealType getMin () const {
    return _min;
  }

  inline const RealType getMax () const {
    return _max;
  }

  //! Generate rgb hex notation used in html, gnuplot, ...
  static string gnuplotColorString( const RealType* rgb );
  static string gnuplotColorString( const aol::MultiVector < RealType >& rgb, int i );

};

}

#endif
