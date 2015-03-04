#ifndef GLOBALSCANDISTORTIONCORRECTION_H
#define GLOBALSCANDISTORTIONCORRECTION_H

#include <aol.h>


template <typename _RealType, typename _PictureType>
_RealType normXCorr ( const _PictureType &Data, const short x, const short y, const short dy, const short L, const short s ) {
  aol::Vector<_RealType> band1 ( 2 * L + 1 ), band2 ( 2 * L + 1 );
  for ( int dx=-L; dx<=L ; ++dx ) {
    band1[dx+L] = Data.get ( x+dx+s, y );
    band2[dx+L] = Data.get ( x+dx, y+dy );
  }
  _RealType mean0 = band1.getMeanValue ( ), meanDY = band2.getMeanValue ( ), sum0Sqr = 0, sumDYSqr = 0, res = 0, val0, valDY;
  for ( int dx=-L; dx<=L ; ++dx ) {
    val0 = Data.get ( x+dx+s, y ) - mean0;
    valDY = Data.get ( x+dx, y+dy ) - meanDY;
    res += val0 * valDY;
    sum0Sqr += aol::Sqr<_RealType> ( val0 );
    sumDYSqr += aol::Sqr<_RealType> ( valDY );
  }
  _RealType sumProd = sum0Sqr * sumDYSqr;
  if ( sumProd > 0 )
    res /= sqrt ( sumProd );
  
  return res;
}


template <typename _RealType, typename _PictureType>
void jitterbug ( const _PictureType &Data, _PictureType &Shifts ) {
  const short NX = Data.getNumX ( ), NY = Data.getNumY ( );
  
  const short L = 12, S = 11;
  
  _PictureType est ( Data );
  Shifts.reallocate ( NX, NY );
  
  for ( short y=1; y<NY ; ++y ) {
    // Determine optimal shifts in row y
    for ( short x=(L+S); x<(NX-S-L) ; ++x ) {
      // Determine maximum of X-Corr. between a horizontal band of size (2L+1) around pixel (x,y) and the pixels (x,y+1) and (x,y-1) within a range (-S,S)
      _RealType optNormXCorr = 0, optS = 0, normXCorrTop, normXCorrBot, normXCorrAvg;
      for ( short s=-S; s<=S ; ++s ) {
        normXCorrTop = normXCorr<_RealType, _PictureType> ( est, x, y, -1, L, s );
        normXCorrBot = normXCorr<_RealType, _PictureType> ( est, x, y, -1, L, s );
        normXCorrAvg = 0.5 * ( normXCorrTop + normXCorrBot );
        if ( normXCorrAvg > optNormXCorr ) {
          optNormXCorr = normXCorrAvg;
          optS = s;
        }
      }
      if ( optNormXCorr >= 0.2 )
        Shifts.set ( x, y, ( optNormXCorr - 0.2 ) / 0.8 * optS );
    }
    
    // Smooth shifts
    aol::Vector<_RealType> shifts ( NX );
    for ( short x=0; x<NX ; ++x )
      shifts[x] = Shifts.get ( x, y );
    
    // TODO: replace moving averages with gaussian filter
    const short WindowSize = 7, WindowOffset = ( WindowSize - 1 ) / 2;
    aol::Vector<_RealType> smoothedShifts ( shifts );
    for ( short x=WindowOffset; x<shifts.size ( ) - WindowOffset ; ++x ) {
      _RealType avg = 0;
      for ( short dx=-WindowOffset; dx<=WindowOffset ; ++dx )
        avg += shifts[x + dx];
      smoothedShifts[x] = avg / WindowSize;
    }
    
    // Resample row y
    for ( short x=0; x<NX ; ++x )
      est.set ( x, y, Data.get ( x + round ( smoothedShifts[x] ), y ) );
  }
}


#endif