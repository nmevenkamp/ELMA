#ifndef PATTERNANALYSIS_H_
#define PATTERNANALYSIS_H_

#include <aol.h>
#include <multiArray.h>
#include <convolution.h>
#include <gnuplotter.h>
#include <FEOpInterface.h>
#include <scalarArray.h>
#include <solver.h>
#include <fastUniformGridMatrix.h>
#include <linearSmoothOp.h>
#include <configurators.h>
#include <generator.h>
#include "nonLinearRegression.h"
#include "linearRegression.h"
#include <cellCenteredGrid.h>
#include <linearSmoothOp.h>


template <typename _RealType, typename _PictureType>
class PatternAnalyzer {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
  typedef qc::MultiArray<RealType, qc::QC_2D, 3> ColoredPictureType;
  typedef qc::RectangularGridConfigurator<RealType, qc::QC_1D, aol::GaussQuadrature<RealType,qc::QC_1D,3> > ConfType;
protected:
  const PictureType &_data;
  const aol::Vec2<RealType> _center;
  const std::string &_outputDir;
  const bool _verbose;
  aol::ProgressBar<> *_progressBar;
  mutable sigfunc _previousCtrlCHandler;
  const qc::FastILexMapper<qc::QC_2D> _mapper;
  PictureType _fourierPowerCoefficients, _fourierPowerPeaks;
  aol::RandomAccessContainer<aol::Vec2<RealType> > _peaks;
  aol::Vector<RealType> _periodicityAnglesRadians, _periodicityAnglesDegrees, _periodicitySpacingsPixels;
public:
  PatternAnalyzer ( const PictureType &Data, const std::string &OutputDir = "", const bool Verbose = false,
                    aol::ProgressBar<> *ProgressBar = NULL )
    : _data ( Data ), _center ( Data.getNumX ( ) / 2, Data.getNumY ( ) / 2 ), _outputDir ( OutputDir ), _verbose ( Verbose ),
      _progressBar ( ProgressBar ),
      _mapper ( Data.getNumX ( ), Data.getNumY ( ) ),
      _fourierPowerCoefficients ( Data.getNumX ( ), Data.getNumY ( ) ), _fourierPowerPeaks ( Data.getNumX ( ), Data.getNumY ( ) ),
      _peaks ( ), _periodicityAnglesRadians ( ), _periodicityAnglesDegrees ( ), _periodicitySpacingsPixels ( )
  {
    // Calculate fourier power coefficients
    const RealType mean = _data.getMeanValue ( );
    qc::MultiArray<RealType, 2, 2> dataEliminatedMean ( _data.getNumX ( ), _data.getNumY ( ) );
    for ( short i=0; i<_data.getNumX ( ) ; ++i )
      for ( short j=0; j<_data.getNumY ( ) ; ++j )
        dataEliminatedMean[0].set ( i, j, _data.get ( i, j ) - mean );
    qc::ScalarArray<RealType, qc::QC_2D> modulus ( _data.getNumX ( ), _data.getNumY ( ) );
    qc::computeLogFFTModulus<RealType> ( dataEliminatedMean[0], modulus, 0, false );
    _fourierPowerCoefficients.rotate90From ( modulus );
    _fourierPowerCoefficients.scaleValuesTo01 ( );

    // Find fourier power peaks
    _fourierPowerPeaks = _fourierPowerCoefficients;
    std::pair<int, RealType> maxIndVal;
    aol::Vec2<RealType> peakPos;
    while ( _peaks.size ( ) < 2 ) {
      maxIndVal = _fourierPowerPeaks.getMaxIndexAndValue ( );
      for ( short dx=-1; dx<=1 ; ++dx )
        for ( short dy=-1; dy<=1 ; ++dy )
          _fourierPowerPeaks.set ( _mapper.splitGlobalIndex ( maxIndVal.first )[0] + dx, _mapper.splitGlobalIndex ( maxIndVal.first )[1] + dy, 0 );
      peakPos.set ( _mapper.splitGlobalIndex ( maxIndVal.first )[0], _mapper.splitGlobalIndex ( maxIndVal.first )[1] + 1 ); // TODO: why is the center incorrect?
      peakPos -= _center;
      bool angleTooSmall = false;
      for ( short k=0; k<_peaks.size ( ) ; ++k ) {
        const RealType p1dotp2 = peakPos.dotProduct ( _peaks[k] ) / ( peakPos.norm ( ) * _peaks[k].norm ( ) );
        if ( p1dotp2 < -0.9 || p1dotp2 > 0.9 )
          angleTooSmall = true;
      }
      if ( !angleTooSmall ) {
        _peaks.pushBack ( peakPos );
        _fourierPowerPeaks[maxIndVal.first] = -1;
      }
    }
    _fourierPowerPeaks.clamp ( -1, 0 );
    _fourierPowerPeaks *= -1;

    // Calculate angles of main axes of periodicity from the peak positions
    _periodicityAnglesRadians.resize ( 2 );
    _periodicityAnglesDegrees.resize ( 2 );
    for ( short i=0; i<2 ; ++i ) {
      _periodicityAnglesRadians[i] = atan ( _peaks[i][1] / _peaks[i][0] );
      _periodicityAnglesDegrees[i] = _periodicityAnglesRadians[i] * 180 / aol::NumberTrait<RealType>::pi;
    }
    
    /*
     * BEGIN: Calculate periodicity spacings
     */
    // Preprocess data (to make peak finding and sine fitting more robust)
    const short filterSize = ( _data.getMaxValue ( ) <= 11 ) ? 7 : 5, filterOffset = ( filterSize - 1 ) / 2;
    PictureType preprocessedData ( _data.getNumX ( ) - filterSize + 1, _data.getNumY ( ) - filterSize + 1 );
    PictureType block ( filterSize, filterSize );
    for ( int x=filterOffset; x<_data.getNumX ( )-filterOffset ; ++x ) {
      for ( int y=filterOffset; y<_data.getNumY ( )-filterOffset ; ++ y) {
        _data.copyBlockTo ( x-filterOffset, y-filterOffset, block );
        preprocessedData.set ( x-filterOffset, y-filterOffset, block.getMeanValue ( ) );
      }
    }
    
    if ( _verbose && _outputDir.size ( ) > 0 ) {
      std::stringstream ss;
      ss << _outputDir << "/preprocessedData.png";
      preprocessedData.setOverflowHandling ( aol::CLIP_THEN_SCALE, preprocessedData.getMinValue ( ), preprocessedData.getMaxValue ( ) );
      preprocessedData.savePNG ( ss.str ( ).c_str ( ) );
    }
    
    
    _periodicitySpacingsPixels.resize ( 2 );
    
    // Find brightest peak and see which of the two axes offers the larger intersection with the image
    const qc::FastILexMapper<qc::QC_2D> mapper ( preprocessedData.getNumX ( ), preprocessedData.getNumY ( ) );
    aol::Vec2<short> peak ( mapper.splitGlobalIndex ( preprocessedData.getMaxIndexAndValue ( ).first )[0],
                            mapper.splitGlobalIndex ( preprocessedData.getMaxIndexAndValue ( ).first )[1] );
    aol::RandomAccessContainer<std::vector<std::pair<RealType, RealType> > > intensitiesContainer ( 2 );
    for ( int k=0; k<2 ; ++k ) setIntensitiesAlongAxis ( intensitiesContainer[k], preprocessedData, _periodicityAnglesRadians[k], peak );
    const short kFirst = ( intensitiesContainer[0].size ( ) > intensitiesContainer[1].size ( ) ) ? 0 : 1;

    // Get periodicity spacing along primary axis kFirst
    RealType energy;
    _periodicitySpacingsPixels[kFirst] = getPeriodicitySpacing ( intensitiesContainer[kFirst], energy, 1 );
    
    // Move origin of secondary axis along primary axis from brightest peak in periodicity steps,
    // and extract intensities from where the intersection of secondary axis with the image is largest
    std::vector<std::pair<RealType, RealType> > intensities;
    setIntensitiesAlongAxis ( intensities, preprocessedData, _periodicityAnglesRadians[1-kFirst], peak );
    short direction = 1;
    aol::Vec2<RealType> pos ( peak[0], peak[1] );
    aol::Vec2<RealType> stepVector ( direction * _periodicitySpacingsPixels[kFirst] * cos ( _periodicityAnglesRadians[kFirst] ),
                                     direction * _periodicitySpacingsPixels[kFirst] * sin ( _periodicityAnglesRadians[kFirst] ) );
    short maxIntensitiesSize = 0;
    aol::Vec2<short> maxIntersectionOrigin;
    while ( true ) {
      pos += stepVector;
      if ( pos[0] < 1 || pos[0] >= preprocessedData.getNumX ( ) -1 || pos[1] < 1 || pos[1] >= preprocessedData.getNumY ( ) - 1 ) {
        if ( direction == -1 ) break;
        else {
          direction = -1;
          stepVector *= direction;
          pos.set ( peak[0], peak[1] );
          pos += stepVector;
        }
      }
      RealType localMax = 0;
      aol::Vec2<short> localPos, localMaxPos;
      for ( int dx=-1; dx<=1 ; ++dx ) {
        for ( int dy=-1; dy<=1 ; ++dy ) {
          localPos.set ( round ( pos[0] ) + dx, round ( pos[1] ) + dy );
          if ( localPos[0] >= 0 && localPos[0] < preprocessedData.getNumX ( ) && localPos[1] >= 0 && localPos[1] < preprocessedData.getNumY ( )
            && preprocessedData.get ( localPos ) > localMax ) {
            localMax = preprocessedData.get ( round ( pos[0] ) + dx, round ( pos[1] ) + dy );
            localMaxPos.set ( localPos );
          }
        }
      }
      pos.set ( localMaxPos[0], localMaxPos[1] );
      setIntensitiesAlongAxis ( intensities, preprocessedData, _periodicityAnglesRadians[1-kFirst], localMaxPos );
      if ( intensities.size ( ) > maxIntensitiesSize ) {
        maxIntensitiesSize = intensities.size ( );
        maxIntersectionOrigin.set ( localMaxPos );
      }
    }
    setIntensitiesAlongAxis ( intensities, preprocessedData, _periodicityAnglesRadians[1-kFirst], maxIntersectionOrigin );
    
    // Get periodicity spacing along secondary axis 1-kFirst
    _periodicitySpacingsPixels[1-kFirst] = getPeriodicitySpacing ( intensities, energy, 2 );
    
    /*
     * END: Calculate periodicity spacings
     */
    
    
    // Create some output for debugging/analysis (if requested)
    if ( _verbose && _outputDir.size ( ) > 0 ) {
      std::stringstream ss;
      ss << _outputDir << "/periodicityFourierPeaks.png";
      saveFourierCoefficients ( ss.str ( ).c_str ( ), true, true );
      
      ss.str ( std::string ( ) ); // clear stringstream
      ss << _outputDir << "/periodicityAxes.png";
      saveDataPeriodicityAxesImg ( ss.str ( ).c_str ( ) );
      
      ss.str ( std::string ( ) ); // clear stringstream
      ss << _outputDir << "/periodicityPattern.png";
      saveDataPeriodicPatternImg ( ss.str ( ).c_str ( ) );
    }
    
    if ( _verbose )
      std::cerr << _periodicitySpacingsPixels << std::endl;
    
    if ( _verbose && _outputDir.size ( ) > 0 ) {
      std::stringstream ss;
      ss << _outputDir << "/periodicityAnalysis.txt";
      ofstream txtFile;
      txtFile.open ( ss.str ( ).c_str ( ) );
      txtFile << "Estimated grid parameters" << std::endl;
      txtFile << std::endl;
      txtFile << "Delta x_1 = " << _periodicitySpacingsPixels[0] << " pixels" << std::endl;
      txtFile << "Delta x_2 = " << _periodicitySpacingsPixels[1] << " pixels" << std::endl;
      txtFile << std::endl;
      txtFile << "alpha_1 = " << _periodicityAnglesDegrees[0] << " degrees" << std::endl;
      txtFile << "alpha_2 = " << _periodicityAnglesDegrees[1] << " degrees" << std::endl;
      txtFile << "alpha_1 = " << _periodicityAnglesRadians[0] << " radians" << std::endl;
      txtFile << "alpha_2 = " << _periodicityAnglesRadians[1] << " radians" << std::endl;
      txtFile.close ( );
    }
  }

  const aol::RandomAccessContainer<aol::Vec2<RealType> >& getPeakPositions ( ) const {
    return _peaks;
  }

  const aol::Vector<RealType>& getPeriodicitySpacingsPixels ( ) const {
    return _periodicitySpacingsPixels;
  }
  
  const aol::Vector<RealType>& getPeriodicityAnglesRadians ( ) const {
    return _periodicityAnglesRadians;
  }

  const aol::Vector<RealType>& getPeriodicityAnglesDegrees ( ) const {
    return _periodicityAnglesDegrees;
  }
  
  RealType getPeriodicityAnglesRadians ( const short Axis ) const {
    if ( Axis < 0 || Axis > 2 )
      throw aol::Exception ( "Argument \"Axis\" has to be between 0 and 2!", __FILE__, __LINE__ );
    return _periodicityAnglesRadians[Axis];
  }

  RealType getPeriodicityAnglesDegrees ( const short Axis ) const {
    if ( Axis < 0 || Axis > 2 )
      throw aol::Exception ( "Argument \"Axis\" has to be between 0 and 2!", __FILE__, __LINE__ );
    return _periodicityAnglesDegrees[Axis];
  }

  void saveFourierCoefficients ( const char* Path, const bool LogScale = false, const bool RedPeaks = false, const bool BlueAxes = false ) {
    PictureType u ( _fourierPowerCoefficients );
    RealType uMax = u.getMaxValue ( );
    if ( LogScale ) {
      for ( int k = 0; k<u.size ( ) ; ++k )
        u[k] = log ( 1 + u[k] / uMax * 255 );
      uMax = u.getMaxValue ( );
    }
    
    if ( RedPeaks ) {
      ColoredPictureType v ( u.getNumX ( ), u.getNumY ( ) );
      v[0] = u;
      v[1] = u;
      v[2] = u;
      for ( int k=0; k<_peaks.size ( ) ; ++k ) {
        v[0].set ( _peaks[k][0] + _center[0], _peaks[k][1] + _center[1] - 1, uMax );
        v[1].set ( _peaks[k][0] + _center[0], _peaks[k][1] + _center[1] - 1, 0 );
        v[2].set ( _peaks[k][0] + _center[0], _peaks[k][1] + _center[1] - 1, 0 );
      }
      v.setOverflowHandling ( aol::CLIP_THEN_SCALE, v.getMinValue ( ), v.getMaxValue ( ) );
      v.savePNG ( Path );
    } else if ( BlueAxes ) {
      ColoredPictureType v ( u.getNumX ( ), u.getNumY ( ) );
      v[0] = u;
      v[1] = u;
      v[2] = u;
      aol::Vec2<short> pos;
      for ( short k=0; k<2 ; ++k ) {
        for ( short i=-_data.getNumX ( ); i<_data.getNumX ( ) ; ++i ) {
          pos.set ( _center[0] + cos ( _periodicityAnglesRadians[k] ) * i, _center[1] + sin ( _periodicityAnglesRadians[k] ) * i );
          if ( pos[0] >= 0 && pos[0] < _data.getNumX ( ) && pos[1] >= 0 && pos[1] < _data.getNumY ( ) ) {
            v[0].set ( pos, 0 );
            v[1].set ( pos, 0 );
            v[2].set ( pos, uMax );
          }
        }
      }
      v.setOverflowHandling ( aol::CLIP_THEN_SCALE, v.getMinValue ( ), v.getMaxValue ( ) );
      v.savePNG ( Path );
    } else
      u.save ( Path, qc::PGM_DOUBLE_BINARY );
  }

  void saveFourierPeaks ( const char* Path ) {
    _fourierPowerPeaks.save ( Path, qc::PGM_DOUBLE_BINARY );
  }

  void saveDataPeriodicityAxesImg ( const char* Path, const aol::Vec2<short> &Origin = aol::Vec2<short> ( -1, -1 ), const short Axis = -1 ) {
    if ( Axis < -1 || Axis > 1 )
      throw aol::Exception ( "Argument \"Axis\" has to be between -1 and 1!", __FILE__, __LINE__ );
    
    aol::Vec2<short> origin ( Origin );
    if ( origin[0] < 0 || origin[0] >= _data.getNumX ( ) || origin[1] < 0 || origin[1] >= _data.getNumY ( ) )
      origin.set ( _mapper.splitGlobalIndex ( _data.getMaxIndexAndValue ( ).first )[0], _mapper.splitGlobalIndex ( _data.getMaxIndexAndValue ( ).first )[1] );
    
    ColoredPictureType periodicityAxesImg ( _data.getNumX ( ), _data.getNumY ( ) );
    periodicityAxesImg[0] = _data;
    periodicityAxesImg[1] = _data;
    periodicityAxesImg[2] = _data;
    aol::Vec2<short> pos;
    aol::Vector<short> axes;
    if ( Axis == -1 ) {
      axes.pushBack ( 0 );
      axes.pushBack ( 1 );
    } else axes.pushBack ( Axis );
    for ( short k=0; k<axes.size ( ) ; ++k ) {
      for ( short i=-_data.getNumX ( ); i<_data.getNumX ( ) ; ++i ) {
        pos.set ( origin[0] + cos ( _periodicityAnglesRadians[axes[k]] ) * i, origin[1] + sin ( _periodicityAnglesRadians[axes[k]] ) * i );
        if ( pos[0] >= 0 && pos[0] < _data.getNumX ( ) && pos[1] >= 0 && pos[1] < _data.getNumY ( ) ) {
          periodicityAxesImg[0].set ( pos, 0 );
          periodicityAxesImg[1].set ( pos, 0 );
          periodicityAxesImg[2].set ( pos, _data.getMaxValue ( ) );
        }
      }
    }
    periodicityAxesImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, _data.getMinValue ( ), _data.getMaxValue ( ) );
    periodicityAxesImg.savePNG ( Path );
  }
  
  void saveDataPeriodicPatternImg ( const char* Path, const aol::Vec2<short> &Origin = aol::Vec2<short> ( -1, -1 ), const short SearchWindowSize = 1 ) {
    aol::Vec2<short> origin ( Origin );
    if ( origin[0] < 0 || origin[0] >= _data.getNumX ( ) || origin[1] < 0 || origin[1] >= _data.getNumY ( ) )
      origin.set ( _mapper.splitGlobalIndex ( _data.getMaxIndexAndValue ( ).first )[0], _mapper.splitGlobalIndex ( _data.getMaxIndexAndValue ( ).first )[1] );
    
    const short searchWindowOffset = ( SearchWindowSize - 1 ) / 2;

    ColoredPictureType periodicityAxesImg ( _data.getNumX ( ), _data.getNumY ( ) );
    for ( short i=0; i<_data.getNumX ( ) ; ++i ) {
      for ( short j=0; j<_data.getNumY ( ) ; ++j ) {
        periodicityAxesImg[0].set ( i, j, _data.get ( i, j ) );
        periodicityAxesImg[1].set ( i, j, _data.get ( i, j ) );
        periodicityAxesImg[2].set ( i, j, _data.get ( i, j ) );
      }
    }
    
    aol::Vec2<short> pos1, pos2, pos3;
    pos1.set ( origin );
    for ( short i=-_data.getNumX ( ); i<_data.getNumX ( ) ; ++i ) {
      pos1.set ( origin[0] + cos ( _periodicityAnglesRadians[0] ) * i * _periodicitySpacingsPixels[0],
                 origin[1] + sin ( _periodicityAnglesRadians[0] ) * i * _periodicitySpacingsPixels[0] );
      for ( short j=-_data.getNumX ( ); j<_data.getNumX ( ) ; ++j ) {
        pos2.set ( pos1[0] + cos ( _periodicityAnglesRadians[1] ) * j * _periodicitySpacingsPixels[1],
                   pos1[1] + sin ( _periodicityAnglesRadians[1] ) * j * _periodicitySpacingsPixels[1] );
        for ( short dx=-searchWindowOffset; dx<=searchWindowOffset ; ++dx ) {
          for ( short dy=-searchWindowOffset; dy<=searchWindowOffset ; ++dy ) {
            pos3.set ( pos2[0] + dx, pos2[1] + dy );
            if ( pos3[0] >= 0 && pos3[0] < _data.getNumX ( ) && pos3[1] >= 0 && pos3[1] < _data.getNumY ( ) ) {
              periodicityAxesImg[0].set ( pos3, _data.getMaxValue ( ) );
              periodicityAxesImg[1].set ( pos3, 0 );
              periodicityAxesImg[2].set ( pos3, 0 );
            }
          }
        }
      }
    }
    periodicityAxesImg[0].set ( origin, 0 );
    periodicityAxesImg[1].set ( origin, _data.getMaxValue ( ) );
    periodicityAxesImg.setOverflowHandling ( aol::SCALE, 0, 255 );
    periodicityAxesImg.savePNG ( Path );
  }
  
  void saveIntensityPlotAlongPeriodicAxis ( const char* Path, const short Axis, const aol::Vec2<short> &Origin = aol::Vec2<short> ( -1, -1 ) ) {
    if ( Axis < 0 || Axis > 1 )
      throw aol::Exception ( "Argument \"Axis\" has to be between 0 and 1!", __FILE__, __LINE__ );
    
    aol::Vec2<short> origin ( Origin );
    if ( origin[0] < 0 || origin[0] >= _data.getNumX ( ) || origin[1] < 0 || origin[1] >= _data.getNumY ( ) )
      origin.set ( _mapper.splitGlobalIndex ( _data.getMaxIndexAndValue ( ).first )[0], _mapper.splitGlobalIndex ( _data.getMaxIndexAndValue ( ).first )[1] );
    
    std::vector<std::pair<RealType, RealType> > intensities;
    aol::Vec2<short> pos;
    for ( short sign=-1; sign<=1 ; sign+=2 ) {
      int i = ( sign == -1 ) ? 1 : 0;
      pos.set ( origin[0], origin[1] );
      while ( pos[0] >= 0 && pos[0] < _data.getNumX ( ) && pos[1] >= 0 && pos[1] < _data.getNumY ( ) ) {
        intensities.insert ( ( sign == -1 ) ? intensities.begin ( ) : intensities.end ( ),
                             std::pair<RealType, RealType> ( sign * i * aol::Vec2<RealType> ( cos ( _periodicityAnglesRadians[Axis] ),
                                                                                              sin ( _periodicityAnglesRadians[Axis] ) ).norm ( ),
                                                                                              _data.get ( pos ) ) );
        ++i;
        pos.set ( origin[0] + sign * cos ( _periodicityAnglesRadians[Axis] ) * i, origin[1] + sign * sin ( _periodicityAnglesRadians[Axis] ) * i );
      }
    }
  
    aol::Plotter<RealType> plotter;
    plotter.set_outfile_base_name ( Path );
    aol::PlotDataFileHandler<RealType> plotHandler;
    plotHandler.generateFunctionPlot ( intensities );
    plotter.addPlotCommandsFromHandler ( plotHandler );
    plotter.genPlot( aol::GNUPLOT_PNG );
  }
  
protected:
  void setIntensitiesAlongAxis ( std::vector<std::pair<RealType, RealType> > &Intensities,
                                 const PictureType &Data, const RealType AngleRadians,
                                 const aol::Vec2<short> &Origin = aol::Vec2<short> ( -1, -1 ) ) {
    Intensities.clear ( );
    const qc::FastILexMapper<qc::QC_2D> mapper ( Data.getNumX ( ), Data.getNumY ( ) );
    
    aol::Vec2<short> origin ( Origin );
    if ( origin[0] < 0 || origin[0] >= Data.getNumX ( ) || origin[1] < 0 || origin[1] >= Data.getNumY ( ) )
      origin.set ( mapper.splitGlobalIndex ( Data.getMaxIndexAndValue ( ).first )[0], mapper.splitGlobalIndex ( Data.getMaxIndexAndValue ( ).first )[1] );
    
    aol::Vec2<short> pos;
    for ( short sign=-1; sign<=1 ; sign+=2 ) {
      int i = ( sign == -1 ) ? 1 : 0;
      pos.set ( origin[0], origin[1] );
      while ( pos[0] >= 0 && pos[0] < Data.getNumX ( ) && pos[1] >= 0 && pos[1] < Data.getNumY ( ) ) {
        Intensities.insert ( ( sign == -1 ) ? Intensities.begin ( ) : Intensities.end ( ),
                             std::pair<RealType, RealType> ( sign * i * aol::Vec2<RealType> ( cos ( AngleRadians ), sin ( AngleRadians ) ).norm ( ),
                                                             Data.get ( pos ) ) );
        ++i;
        pos.set ( origin[0] + sign * cos ( AngleRadians ) * i, origin[1] + sign * sin ( AngleRadians ) * i );
      }
    }
  }
  
  RealType getPeriodicitySpacing ( const std::vector<std::pair<RealType, RealType> > &Intensities, RealType &Energy, const int OutputNr = 0 ) {
    // Calculate mean value and create intensities vector with zero mean
    aol::Vector<RealType> intensities;
    for ( int i=0; i<Intensities.size ( ) ; ++i )
      intensities.pushBack ( Intensities[i].second );
    const RealType mean = intensities.getMeanValue ( ), maxZeroMean = intensities.getMaxValue ( ) - mean;
    std::vector<std::pair<RealType, RealType> > intensitiesZeroMean;
    for ( int i=0; i<Intensities.size ( ) ; ++i )
      intensitiesZeroMean.push_back ( std::pair<RealType, RealType> ( Intensities[i].first, Intensities[i].second - mean ) );
    
    // Define energy, derivative and second derivative of a non-linear least squares sum of sines fit target functional
    const int numTerms = 1;
    SumOfSinesTargetFunctional<RealType> F ( intensitiesZeroMean, numTerms );
    SumOfSinesTargetJacobian<RealType, aol::FullMatrix<RealType> > DF ( intensitiesZeroMean, numTerms );
    
    // Define a trust region search method to find optimal parameters (especially the frequency of the sin functions)
    LevenbergMarquardtAlgorithm<RealType, aol::FullMatrix<RealType>, LinearRegressionQR<RealType> > levenbergMarquardtAlg ( intensitiesZeroMean.size ( ), F, DF, 50, 1, 0.2, 0.8, 1e-6, 1e-6, 1e-6, false );
    aol::Vector<RealType> arg ( 3 * numTerms ), dest ( 3 * numTerms );
    arg[0] = maxZeroMean;
    
    // Try different initial values (trust region method does not seem to converge globally)
    aol::Vector<RealType> energies, frequencies, optDest ( dest.size ( ) );
    aol::Vector<RealType> diffs ( intensitiesZeroMean.size ( ) );
    RealType energy = 0;

    setCtrlCHandler ( );
    std::stringstream ss;
    ss << "PatternAnalysis: computing spacing along axis #" << OutputNr << " (step " << OutputNr << "/2)";
    if ( _progressBar != NULL ) _progressBar->setText ( ss.str( ).c_str() );
    if ( _progressBar != NULL ) _progressBar->start ( intensitiesZeroMean.size ( ) - 1 );
    for ( int dx=1; dx<intensitiesZeroMean.size ( ) && !wantsInterrupt ( ) ; ++dx ) {
      arg[1] = 2 * aol::NumberTrait<RealType>::pi / dx;
      levenbergMarquardtAlg.apply ( arg, dest );
      F.apply ( dest, diffs );
      energy = diffs.normSqr ( );
      energies.pushBack ( energy );
      frequencies.pushBack ( dest[1] );
      if ( energy == energies.getMinValue ( ) )
        optDest = dest;
      if ( _progressBar != NULL ) (*_progressBar)++;
    }
    if ( _progressBar != NULL ) _progressBar->finish ( );
    unsetCtrlCHandler ( );
    
    if ( _verbose ) {
      std::cerr << "Optimal initial frequency=" << 2 * aol::NumberTrait<RealType>::pi / ( energies.getMinIndexAndValue ( ).first + 1)
                << "; optimal frequency=" << frequencies[energies.getMinIndexAndValue ( ).first]
                << "; minimal energy=" << energies[energies.getMinIndexAndValue ( ).first] << std::endl;
    }
    
    RealType periodicitySpacing = aol::Abs<RealType> ( 2 * aol::NumberTrait<RealType>::pi / frequencies[energies.getMinIndexAndValue ( ).first] );
    
    if ( _verbose && _outputDir.size ( ) > 0 ) {
      std::vector<std::pair<RealType, RealType> > sumOfSines;
      for ( int i=0; i<intensitiesZeroMean.size ( ) ; ++i )
        sumOfSines.push_back ( std::pair<RealType, RealType> ( intensitiesZeroMean[i].first,
                                                               SumOfSinesLeastSquaresEnergy<RealType>::sumOfSines ( intensitiesZeroMean[i].first, optDest, numTerms ) ) );
      
      std::stringstream ss;
      ss << _outputDir << "/intensitiesZeroMean_vs_sumOfSines_" << OutputNr;
      aol::Plotter<RealType> plotter;
      plotter.set_outfile_base_name ( ss.str ( ).c_str ( ) );
      aol::PlotDataFileHandler<RealType> plotHandler;
      plotHandler.generateFunctionPlot ( intensitiesZeroMean );
      plotHandler.generateFunctionPlot ( sumOfSines );
      plotter.addPlotCommandsFromHandler ( plotHandler );
      plotter.genPlot( aol::GNUPLOT_PNG );
      
      ss.str ( std::string ( ) ); // clear stringstream
      ss << _outputDir << "/intensities_vs_sumOfSines_" << OutputNr << ".csv";
      std::ofstream txtFile ( ss.str ( ).c_str ( ) );
      for ( int i=0; i<intensitiesZeroMean.size ( ) ; ++i )
        txtFile << Intensities[i].first << ", " << Intensities[i].second << ", " << ( sumOfSines[i].second + mean ) << std::endl;
      txtFile.close ( );
    }
    
    if ( _verbose )
      std::cerr << "Periodicity spacing: " << periodicitySpacing << " pixels." << std::endl;
    
    Energy = energies.getMinValue ( );
    return periodicitySpacing;
  }
  
  void setCtrlCHandler () const {
    _previousCtrlCHandler = signal ( InterruptSignal, aol::ctrlCHandler );
  }
  
  void unsetCtrlCHandler () const {
    signal ( InterruptSignal, _previousCtrlCHandler );
  }
  
  bool wantsInterrupt() const {
    if (!aol::getCtrlCState())
      return false;
    else
      return true;
  }
};


template <typename _RealType, typename _PictureType>
class PatternAnalyzerBergmann {
  typedef _RealType RealType;
  typedef aol::Matrix22<int> MatrixType;
  typedef _PictureType PictureType;
  typedef qc::MultiArray<RealType, qc::QC_2D, 3> ColoredPictureType;
protected:
  const PictureType &_data;
  const std::string &_outputDir;
  const bool _verbose;
  const qc::FastILexMapper<qc::QC_2D> _mapper;
  PictureType _fourierPowerCoefficients, _fourierPowerPeaks;
  MatrixType _patternNormalForm;
  aol::RandomAccessContainer<aol::Vec2<short> > _pattern;
public:
  PatternAnalyzerBergmann ( const PictureType &Data, const std::string &OutputDir = "", const bool Verbose = false )
  : _data ( Data ), _outputDir ( OutputDir ), _verbose ( Verbose ), _mapper ( Data.getNumX ( ), Data.getNumY ( ) ),
  _fourierPowerCoefficients ( Data.getNumX ( ), Data.getNumY ( ) ), _fourierPowerPeaks ( Data.getNumX ( ), Data.getNumY ( ) ),
  _patternNormalForm ( ) {
    setFourierPowerCoefficients ( );
    setFourierPowerPeaks ( );
    setPatternNormalForm ( );
    setPattern ( );
  }
  
  static MatrixType getPatternNormalForm ( const MatrixType &Matrix ) {
    if ( Matrix.det ( ) == 0 )
      throw aol::Exception ( "The specified matrix is not regular!", __FILE__, __LINE__ );
    
    MatrixType pnf ( Matrix );
    // Form upper triangular matrix
    gcdOnRows ( pnf, 1, 0 );
    // Make diagonal positive
    for ( short row=0; row<2 ; ++row ) {
      if ( pnf.get ( row, row ) < 0 ) {
        for ( short col=0; col<2 ; ++col )
          pnf.set ( row, col, -pnf.get ( row, col ) );
      }
    }
    // Make upper non-zero values of a column lie between 0 and
    int f = 0;
    if ( pnf.get ( 0, 1 ) < 0 || pnf.get ( 0, 1 ) >= pnf.get ( 1, 1 ) )
      f = -floor ( pnf.get ( 0, 1 ) / pnf.get ( 1, 1 ) );
    for ( short col=0; col<2 ; ++col )
      pnf.set ( 0, col, pnf.get ( 0, col ) + f * pnf.get ( 1, col ) );
    
    return pnf;
  }
  
  static void setPattern ( aol::RandomAccessContainer<aol::Vec2<short> > &Positions, const MatrixType &M, const aol::Vec2<short> &TorusSize = aol::Vec2<short> ( 0, 0 ) ) {
    Positions.clear ( );
    const RealType stepSize1 = 1.0 / abs ( M.get ( 1, 1 ) );
    aol::Vector<RealType> steps1;
    RealType curStep = -0.5;
    do {
      steps1.pushBack ( curStep );
      curStep += stepSize1;
    } while ( curStep < 0.5-stepSize1 );
    const RealType stepSize2 = 1.0 / abs ( M.get ( 0, 0 ) );
    aol::Vector<RealType> steps2;
    curStep = -0.5;
    do {
      steps2.pushBack ( curStep );
      curStep += stepSize2;
    } while ( curStep < 0.5-stepSize2 );
    aol::Vector<RealType> tSums ( steps1 );
    tSums *= M.get ( 0, 1 );
    aol::RandomAccessContainer<aol::Vec2<RealType> > positions;
    for ( int i=0; i<steps1.size ( ) ; ++i ) {
      aol::Vec2<RealType> newPos;
      for ( int j=abs ( M.get ( 0, 0 ) ) * i; j<abs ( M.get ( 0, 0 ) ) * ( i + 1 ) ; ++j ) {
        newPos[0] = steps2[j - abs ( M.get ( 0, 0 ) ) * i] + stepSize2 * ( ceil ( tSums[i] ) - tSums[i] );
        newPos[1] = steps1[i];
        positions.pushBack ( newPos );
      }
    }
    aol::Vec2<short> torusSize ( TorusSize );
    if ( torusSize.norm ( ) == 0 )
      torusSize.set ( M.det ( ), M.det ( ) );
    for ( int i=0; i<positions.size ( ) ; ++i ) {
      int x = round ( positions[i][0] * torusSize[0] ), y = round ( positions[i][1] * torusSize[1] );
      x %= static_cast<int> ( torusSize[0] );
      y %= static_cast<int> ( torusSize[1] );
      x = ( x < 0 ) ? x + torusSize[0] : x;
      y = ( y < 0 ) ? y + torusSize[1] : y;
      Positions.pushBack ( aol::Vec2<short> ( x , y ) );
    }
  }
  
  void savePatternDataImage ( const aol::Vec2<short> &Point = aol::Vec2<short> ( -1, -1 ), const std::string &OutputDir = "" ) {
    if ( _outputDir.size ( ) == 0 && OutputDir.size ( ) == 0 )
      throw aol::Exception ( "No output directory specified!", __FILE__, __LINE__ );
    const std::string outputDir = ( _outputDir.size ( ) > 0 ) ? _outputDir : OutputDir;
    
    aol::Vec2<short> point ( Point );
    if ( point[0] < 0 || point[0] >= _data.getNumX ( ) || point[1] < 0 || point[1] >= _data.getNumY ( ) )
      point.set ( _mapper.splitGlobalIndex ( _data.getMaxIndexAndValue ( ).first )[0], _mapper.splitGlobalIndex ( _data.getMaxIndexAndValue ( ).first )[1] );
    translatePatternTo ( point );
    
    ColoredPictureType patternDataImg ( _data.getNumX ( ), _data.getNumY ( ) );
    for ( short i=0; i<3 ; ++i )
      patternDataImg[i] = _data;
    for ( short i=0; i<_pattern.size ( ) ; ++i ) {
      patternDataImg[0].set ( _pattern[i], _data.getMaxValue ( ) );
      patternDataImg[1].set ( _pattern[i], 0 );
      patternDataImg[2].set ( _pattern[i], 0 );
    }
    patternDataImg[0].set ( _pattern[0], 0 );
    patternDataImg[1].set ( _pattern[0], _data.getMaxValue ( ) );
    
    std::stringstream ss;
    ss << outputDir << "/patternDataImg.png";
    patternDataImg.savePNG ( ss.str ( ).c_str ( ) );
  }
  
private:
  void setFourierPowerCoefficients ( ) {
    qc::MultiArray<RealType, 2, 2> function ( _data.getNumX ( ), _data.getNumY ( ) ), transform ( _data.getNumX ( ), _data.getNumY ( ) );
    function[0] = _data;
    qc::FourierTransform ( function, transform );
    for ( short x=0; x<_data.getNumX ( ) ; ++x )
      for ( short y=0; y<_data.getNumY ( ) ; ++y )
        _fourierPowerCoefficients.set ( x, y, aol::Vec2<RealType> ( transform[0].get ( x, y ), transform[1].get ( x, y ) ).norm ( ) );
  }
  
  void setFourierPowerPeaks ( ) {
    eraseFourierPowerPeak ( 0, 0 );
    if ( _verbose && _outputDir.size ( ) > 0 ) {
      std::stringstream ss;
      ss << _outputDir << "/fourierPowerCoefficients_meanErased.pgm";
      _fourierPowerCoefficients.setOverflowHandlingToCurrentValueRange ( );
      _fourierPowerCoefficients.save ( ss.str ( ).c_str ( ), qc::PGM_UNSIGNED_CHAR_BINARY );
    }
    
    for ( int i=0; i<6 ; ++i ) {
      const std::pair<int, RealType> maxIndexAndValue = _fourierPowerCoefficients.getMaxIndexAndValue ( );
      const aol::Vec2<short> maxIndex ( _mapper.splitGlobalIndex ( maxIndexAndValue.first )[0], _mapper.splitGlobalIndex ( maxIndexAndValue.first )[1] );
      _fourierPowerPeaks.set ( maxIndex, 1 );
      eraseFourierPowerPeak ( maxIndex );
    }
    if ( _verbose && _outputDir.size ( ) > 0 ) {
      std::stringstream ss;
      ss << _outputDir << "/fourierPowerPeaks.pgm";
      _fourierPowerPeaks.setOverflowHandlingToCurrentValueRange ( );
      _fourierPowerPeaks.save ( ss.str ( ).c_str ( ), qc::PGM_UNSIGNED_CHAR_BINARY );
    }
  }
  
  void setPatternNormalForm ( ) {
    // 1. Search for two peaks that are closest to the origin
    MatrixType patternMatrix;
    PictureType fourierPowerPeakDistances ( _fourierPowerPeaks );
    RealType maxDist = aol::Vec2<short> ( _fourierPowerPeaks.getNumX ( ) + 1, _fourierPowerPeaks.getNumY ( ) + 1 ).norm ( );
    for ( int x=0; x<_fourierPowerPeaks.getNumX ( ) ; ++x ) {
      for ( int y=0; y<_fourierPowerPeaks.getNumY ( ) ; ++y ) {
        RealType val = _fourierPowerPeaks.get ( x, y );
        fourierPowerPeakDistances.set ( x, y, ( val > 0 ) ? aol::Vec2<short> ( x, y ).norm ( ) : maxDist );
      }
    }
    std::pair<int, RealType> minIndexAndValue = fourierPowerPeakDistances.getMinIndexAndValue ( );
    patternMatrix.set ( 0, 0, _mapper.splitGlobalIndex ( minIndexAndValue.first )[1] );
    patternMatrix.set ( 1, 0, _mapper.splitGlobalIndex ( minIndexAndValue.first )[0] );
    fourierPowerPeakDistances.set ( minIndexAndValue.first, maxDist );
    minIndexAndValue = fourierPowerPeakDistances.getMinIndexAndValue ( );
    patternMatrix.set ( 0, 1, _mapper.splitGlobalIndex ( minIndexAndValue.first )[1] );
    patternMatrix.set ( 1, 1, _mapper.splitGlobalIndex ( minIndexAndValue.first )[0] );
    _patternNormalForm = getPatternNormalForm ( patternMatrix );
  }
  
  void translatePatternTo ( const aol::Vec2<short> &Point ) {
    aol::Vec2<short> translation ( Point );
    translation -= _pattern[0];
    for ( int i=0; i<_pattern.size ( ) ; ++i ) {
      _pattern[i] += translation;
      _pattern[i][0] %= _data.getNumX ( ) - 1;
      _pattern[i][1] %= _data.getNumY ( ) - 1;
      _pattern[i][0] = ( _pattern[i][0] < 0 ) ? _pattern[i][0] + _data.getNumX ( ) - 1 : _pattern[i][0];
      _pattern[i][1] = ( _pattern[i][1] < 0 ) ? _pattern[i][1] + _data.getNumY ( ) - 1 : _pattern[i][1];
    }
  }
  
  void eraseFourierPowerPeak ( const short X, const short Y, const short RubberSize = 3 ) {
    const short offset = ( RubberSize - 1 ) / 2;
    for ( short dx=-offset; dx<=offset ; ++dx ) {
      for ( short dy=-offset; dy<=offset ; ++dy ) {
        if ( X+dx >= 0 && X+dx < _data.getNumX ( ) && Y+dy >= 0 && Y+dy < _data.getNumY ( ) )
          _fourierPowerCoefficients.set ( X+dx, Y+dy, 0 );
      }
    }
  }
  
  void eraseFourierPowerPeak ( const aol::Vec2<short> &Pos, const short RubberSize = 3 ) {
    eraseFourierPowerPeak ( Pos[0], Pos[1], RubberSize );
  }
  
  static void gcdOnRows ( MatrixType &M, const short Ri, const short Ci ) {
    if ( M.get ( Ri, Ci ) != 0 ) {
      // Modify by row addition, such that M[Ci,Ci] is non-negative
      if ( M.get ( Ci, Ci ) < 0 ) {
        for ( short col=0; col<2 ; ++col )
          M.set ( Ci, col, M.get ( Ci, col ) - floor ( M.get ( Ci, Ci ) / M.get ( Ri, Ci ) ) * M.get ( Ri, col ) );
      }
      // Make M[Ci,Ci] positive
      if ( M.get ( Ci, Ci ) == 0 ) {
        for ( short col=0; col<2 ; ++col )
          M.set ( Ci, col, M.get ( Ci, col ) + sign ( M.get ( Ri, Ci ) ) * M.get ( Ri, col ) );
      }
      // Make second entry in that column positive as well
      if ( M.get ( Ri, Ci ) < 0 ) {
        int f = ceil ( M.get ( Ri, Ci ) / M.get ( Ci, Ci ) );
        if ( f == 0 )
          ++f;
        for ( short col=0; col<2 ; ++col )
          M.set ( Ri, col, M.get ( Ri, col ) - f * M.get ( Ci, col ) );
      }
      // Euclidian algorithm on rows in order to get M[Ri,Ci] to zero
      while ( M.get ( Ri, Ci ) != 0 ) {
        if ( abs ( M.get ( Ci, Ci ) ) > abs ( M.get ( Ri, Ci ) ) ) {
          int f = floor ( M.get ( Ci, Ci ) / M.get ( Ri, Ci ) );
          if ( M.get ( Ci, Ci ) % M.get ( Ri, Ci ) == 0 )
            f = f - sign ( M.get ( Ri, Ci ) ) * sign ( M.get ( Ci, Ci ) );
          for ( short col=0; col<2 ; ++col )
            M.set ( Ci, col, M.get ( Ci, col ) - f * M.get ( Ri, col ) );
        } else {
          int f = floor ( M.get ( Ri, Ci ) / M.get ( Ci, Ci ) );
          for ( short col=0; col<2 ; ++col )
            M.set ( Ri, col, M.get ( Ri, col ) - f * M.get ( Ci, col ) );
        }
      }
    }
  }
  
  void setPattern ( ) {
    setPattern ( _pattern, _patternNormalForm, aol::Vec2<short> ( _data.getNumX ( ) - 1, _data.getNumY (  ) - 1 ) );
  }
  
  static int sign (int val) {
    return ( val < 0 ) ? -1 : ( ( val > 0 ) ? 1 : 0 );
  }
};



#endif /* PATTERNANALYSIS_H_ */
