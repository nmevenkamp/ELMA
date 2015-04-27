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
#include <cellCenteredGrid.h>
#include <linearSmoothOp.h>
#include <clustering.h>


template<typename _RealType, typename _PictureType>
int getNumNonNaNs ( const _PictureType &Data ) {
  int res = 0;
  for ( int i=0; i<Data.size ( ); ++i ) {
    if ( !aol::isNaN<_RealType> ( Data[i] ) ) ++res;
  }
  return res;
}


template <typename _RealType, typename _PictureType>
class PeriodicityTargetFunctional : public aol::Op<aol::Vector<_RealType> > {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  const PictureType &_data;
  const int _numNonNaNs;
public:
  PeriodicityTargetFunctional ( const PictureType &Data )
  : _data ( Data ), _numNonNaNs ( getNumNonNaNs<RealType, PictureType> ( Data ) ) { }
  
  void applyAdd ( const aol::Vector<RealType> &Arg, aol::Vector<RealType> &Dest ) const {
    if ( Arg.size ( ) != 4 )
      throw aol::Exception ( "Arguments size must be equal to 4!", __FILE__, __LINE__ );
    
    if ( Dest.size ( ) != 2 * _numNonNaNs )
      throw aol::Exception ( "Destination vector does not match 2 * #non-NaNs!", __FILE__, __LINE__ );
    
    aol::Vec2<RealType> xpv;
    int i = 0;
    for ( int y=0; y<_data.getNumY ( ) ; ++y ) {
      for ( int x=0; x<_data.getNumX ( ) ; ++x ) {
        if ( !aol::isNaN<RealType> ( _data.get ( x, y ) ) ) {
          for ( int j=0; j<2 ; ++j ) {
            xpv.set ( x + Arg[2*j], y + Arg[2*j+1] );
            if ( xpv[0] >= 0 && xpv[0] < _data.getNumX ( ) && xpv[1] >= 0 && xpv[1] < _data.getNumY ( ) && !aol::isNaN<RealType> ( _data.get ( xpv[0], xpv[1] ) ) )
              Dest[i] += _data.get ( x, y ) - _data.interpolate ( xpv );
            ++i;
          }
        }
      }
    }
  }
};

template <typename _RealType, typename _MatrixType, typename _PictureType>
class PeriodicityTargetJacobian : public aol::Op<aol::Vector<_RealType>, _MatrixType> {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  const PictureType &_data;
  const int _numNonNaNs;
public:
  PeriodicityTargetJacobian ( const PictureType &Data )
    : _data ( Data ), _numNonNaNs ( getNumNonNaNs<RealType, PictureType> ( Data ) ) { }
  
  void applyAdd ( const aol::Vector<RealType> &/*Arg*/, aol::FullMatrix<RealType> &/*Dest*/ ) const {
    throw aol::UnimplementedCodeException( "Not implemented", __FILE__, __LINE__ );
  }
  
  void apply ( const aol::Vector<RealType> &Arg, aol::FullMatrix<RealType> &Dest ) const {
    if ( Arg.size ( ) != 4 )
      throw aol::Exception ( "Arguments size must be equal to 4!", __FILE__, __LINE__ );
    
    if ( Dest.getNumRows ( ) != 2 * _numNonNaNs || Dest.getNumCols ( ) != 4 )
      throw aol::Exception ( "Destination dimensions do not fit 2 * $non-NaNs and parameters!", __FILE__, __LINE__ );
    
    Dest.setAll ( 0.0 );
    aol::Vec2<RealType> xpv;
    int i = 0;
    for ( int y=0; y<_data.getNumY ( ) ; ++y ) {
      for ( int x=0; x<_data.getNumX ( ) ; ++x ) {
        if ( !aol::isNaN<RealType> ( _data.get ( x, y ) ) ) {
          for ( int j=0; j<2 ; ++j ) {
            xpv.set ( x + Arg[2*j], y + Arg[2*j+1] );
            if ( xpv[0] >= 0 && xpv[0] < _data.getNumX ( ) && xpv[1] >= 0 && xpv[1] < _data.getNumY ( ) && !aol::isNaN<RealType> ( _data.get ( xpv[0], xpv[1] ) ) ) {
              Dest.set ( i, 2*j, -_data.dxFD ( xpv[0], xpv[1] ) );
              Dest.set ( i, 2*j+1, -_data.dyFD ( xpv[0], xpv[1] ) );
            }
            ++i;
          }
        }
      }
    }
  }
};


enum PatternAnalysisAlgorithmType {
  FourierPeaksAndSineFit,
  ProjectiveStdDevAndPeriodicityEnergyMinimization,
  PeriodicityFunctionalMinimization
};


template <typename _RealType, typename _PictureType>
class PatternAnalyzer {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
  typedef qc::MultiArray<RealType, qc::QC_2D, 3> ColoredPictureType;
  typedef qc::RectangularGridConfigurator<RealType, qc::QC_1D, aol::GaussQuadrature<RealType,qc::QC_1D,3> > ConfType;
protected:
  const std::string &_outputDir;
  const bool _verbose;
  aol::ProgressBar<> *_progressBar;
  mutable sigfunc _previousCtrlCHandler;
public:
  PatternAnalyzer ( const std::string &OutputDir = "", const bool Verbose = false,
                    aol::ProgressBar<> *ProgressBar = NULL )
    : _outputDir ( OutputDir ), _verbose ( Verbose ), _progressBar ( ProgressBar ) { }
  
  void getLatticeVectors ( aol::MultiVector<RealType> &LatticeVectors, const PictureType &Data,
                           const PatternAnalysisAlgorithmType &AlgorithmType = FourierPeaksAndSineFit,
                           const aol::Vector<RealType> &Params = aol::Vector<RealType> ( ) ) {
    LatticeVectors.resize ( 2, 2 );
    if ( AlgorithmType == FourierPeaksAndSineFit ) getLatticeVectorsByFourierPeakExtractionAndSineFit ( LatticeVectors, Data );
    else if ( AlgorithmType == ProjectiveStdDevAndPeriodicityEnergyMinimization ) getLatticeVectorsByProjectiveStdDevAndPeriodicityEnergyMinimization ( LatticeVectors, Data, Params );
    else if ( AlgorithmType == PeriodicityFunctionalMinimization ) getLatticeVectorsByPeriodicityFunctionalMinimization ( LatticeVectors, Data, Params );
    else throw aol::Exception ( "Could not recognize specified pattern analysis algorithm type!" );
  }
  
  void getLatticeAnglesAndPeriods ( aol::Vec2<RealType> &LatticeAngles, aol::Vec2<RealType> &LatticePeriods, const PictureType &Data,
                                    const PatternAnalysisAlgorithmType &AlgorithmType = FourierPeaksAndSineFit,
                                    const aol::Vector<RealType> &Params = aol::Vector<RealType> ( ) ) {
    aol::MultiVector<RealType> latticeVectors;
    getLatticeVectors ( latticeVectors, Data, AlgorithmType, Params );
    for ( int i=0; i<2 ; ++i ) {
      LatticeAngles[i] = atan ( latticeVectors[i][1] / latticeVectors[i][0] );
      LatticePeriods[i] = latticeVectors[i].norm ( );
    }
  }
  
  void refineLatticeVectorsByPeriodicityFunctionalMinimization ( aol::MultiVector<RealType> &LatticeVectors, const PictureType &Data ) {
    PeriodicityTargetFunctional<RealType, PictureType> F ( Data );
    PeriodicityTargetJacobian<RealType, aol::FullMatrix<RealType>, PictureType> DF ( Data );
    LevenbergMarquardtAlgorithm<RealType, aol::FullMatrix<RealType>, LinearRegressionQR<RealType> > levenbergMarquardtAlg ( 2 * getNumNonNaNs<RealType, PictureType> ( Data ),
                                                                                                                           F, DF, 50, 1, 0.2, 0.8, 1e-6, 1e-6, 1e-6, true );
    aol::Vector<RealType> arg ( 4 ), dest ( 4 );
    for ( int i=0; i<2 ; ++i ) {
      arg[2*i] = LatticeVectors[i][0];
      arg[2*i+1] = LatticeVectors[i][1];
    }
    levenbergMarquardtAlg.apply ( arg, dest );
    for ( int i=0; i<2 ; ++i ) {
      LatticeVectors[i][0] = dest[2*i];
      LatticeVectors[i][1] = dest[2*i+1];
    }
    
    // Create some output for debugging/analysis (if requested)
    if ( _verbose && _outputDir.size ( ) > 0 ) {
      std::stringstream ss;
      ss << _outputDir << "/refined_periodicityAxes.png";
      aol::Vec2<RealType> latticeAngles;
      for ( int i=0; i<2 ; ++i ) latticeAngles[i] = atan ( LatticeVectors[i][1] / LatticeVectors[i][0] );
      saveDataPeriodicityAxesImg ( ss.str ( ).c_str ( ), Data, latticeAngles );
      
      ss.str ( std::string ( ) ); // clear stringstream
      ss << _outputDir << "/refined_periodicityPattern.png";
      saveDataPeriodicPatternImg ( ss.str ( ).c_str ( ), Data, LatticeVectors );
    }
  }
  
protected:
  void getLatticeVectorsByPeriodicityFunctionalMinimization ( aol::MultiVector<RealType> &LatticeVectors, const PictureType &Data, const aol::Vector<RealType> &Params ) {
    getLatticeVectorsByProjectiveStdDevAndPeriodicityEnergyMinimization ( LatticeVectors, Data, Params );
    refineLatticeVectorsByPeriodicityFunctionalMinimization ( LatticeVectors, Data );
  }
  
  
  void getLatticeVectorsByProjectiveStdDevAndPeriodicityEnergyMinimization ( aol::MultiVector<RealType> &LatticeVectors, const PictureType &Data, aol::Vec2<int> NumClusters ) {
    /*
     * BEGIN: RevSTEM like angle optimization
     */
    aol::Vec2<RealType> latticeAnglesDegrees, latticeAnglesRadians;
    
    // Calculate projective standard deviations for each angle with 1 degree increments
    const int l0 = 0.1 * sqrt ( Data.size ( ) );
    
    aol::Vector<RealType> projectiveStandardDeviations ( 180 ), projectedAverageIntensities;
    for ( int angleDegrees=0; angleDegrees<180 ; ++angleDegrees ) {
      getProjectedAverageIntensities ( projectedAverageIntensities, Data, angleDegrees, l0 ) ;
      projectiveStandardDeviations[angleDegrees] = projectedAverageIntensities.getStdDev ( );
    }
    
    if ( _verbose && _outputDir.size ( ) > 0 ) {
      std::stringstream ss;
      ss << _outputDir << "/projectiveStandardDeviations";
      aol::Plotter<RealType> plotter;
      plotter.set_outfile_base_name ( ss.str ( ).c_str ( ) );
      aol::PlotDataFileHandler<RealType> plotHandler;
      std::vector<std::pair<RealType, RealType> > data;
      for ( int i=0; i<projectiveStandardDeviations.size ( ) ; ++i ) data.push_back ( std::pair<RealType, RealType> ( i, projectiveStandardDeviations[i] ) );
      plotHandler.generateFunctionPlot ( data );
      plotter.addPlotCommandsFromHandler ( plotHandler );
      plotter.genPlot( aol::GNUPLOT_PNG );
    }
    
    // Find and refine two largest peaks of projective standard deviation
    
    // Find 1st peak
    aol::Vec2<RealType> optParams = fitGaussianToProjectiveStandardDeviationPeak ( projectiveStandardDeviations );
    latticeAnglesDegrees[0] = to360Degrees ( optParams[0] + 90 );
    if ( _verbose ) std::cerr << "Peak #1: " << optParams << std::endl;
    
    // Remove 2nd peak
    for ( int angleDegrees=aol::Max<int> ( optParams[0]-sqrt(optParams[1]), 0 ); angleDegrees<=aol::Min<int> ( optParams[0]+sqrt(optParams[1]), 179 ) ; ++angleDegrees )
      projectiveStandardDeviations[angleDegrees] = 0;
    
    if ( _verbose && _outputDir.size ( ) > 0 ) {
      std::stringstream ss;
      ss << _outputDir << "/projectiveStandardDeviations_peak0Removed";
      aol::Plotter<RealType> plotter;
      plotter.set_outfile_base_name ( ss.str ( ).c_str ( ) );
      aol::PlotDataFileHandler<RealType> plotHandler;
      std::vector<std::pair<RealType, RealType> > data;
      for ( int i=0; i<projectiveStandardDeviations.size ( ) ; ++i ) data.push_back ( std::pair<RealType, RealType> ( i, projectiveStandardDeviations[i] ) );
      plotHandler.generateFunctionPlot ( data );
      plotter.addPlotCommandsFromHandler ( plotHandler );
      plotter.genPlot( aol::GNUPLOT_PNG );
    }
    
    // Find 2nd peak
    optParams = fitGaussianToProjectiveStandardDeviationPeak ( projectiveStandardDeviations );
    latticeAnglesDegrees[1] = to360Degrees ( optParams[0] + 90 );
    if ( _verbose ) std::cerr << "Peak #2: " << optParams << std::endl;
    
    // Remove 2nd peak
    for ( int angleDegrees=aol::Max<int> ( optParams[0]-sqrt(optParams[1]), 0 ); angleDegrees<=aol::Min<int> ( optParams[0]+sqrt(optParams[1]), 179 ) ; ++angleDegrees )
      projectiveStandardDeviations[angleDegrees] = 0;

    if ( _verbose && _outputDir.size ( ) > 0 ) {
      std::stringstream ss;
      ss << _outputDir << "/projectiveStandardDeviations_peak1Removed";
      aol::Plotter<RealType> plotter;
      plotter.set_outfile_base_name ( ss.str ( ).c_str ( ) );
      aol::PlotDataFileHandler<RealType> plotHandler;
      std::vector<std::pair<RealType, RealType> > data;
      for ( int i=0; i<projectiveStandardDeviations.size ( ) ; ++i ) data.push_back ( std::pair<RealType, RealType> ( i, projectiveStandardDeviations[i] ) );
      plotHandler.generateFunctionPlot ( data );
      plotter.addPlotCommandsFromHandler ( plotHandler );
      plotter.genPlot( aol::GNUPLOT_PNG );
    }
    
    if ( _verbose ) {
      for ( int i=0; i<2 ; ++i ) std::cerr << "Angle #" << i << ": " << latticeAnglesDegrees[i] << std::endl;
    }
    
    for ( int i=0; i<2 ; ++i )
      latticeAnglesRadians[i] = latticeAnglesDegrees[i] * aol::NumberTrait<RealType>::pi / 180.0;
    /*
     * END: RevSTEM like angle optimization
     */
    
    
    /*
     * BEGIN: Periodicity optimization based on Energy minimization
     */
    aol::Vec2<RealType> latticePeriods;
    
    // Preprocess data (to make peak finding and sine fitting more robust)
    const short filterSize = ( Data.getMaxValue ( ) <= 11 ) ? 7 : 5, filterOffset = ( filterSize - 1 ) / 2;
    PictureType preprocessedData ( Data.getNumX ( ) - filterSize + 1, Data.getNumY ( ) - filterSize + 1 );
    PictureType block ( filterSize, filterSize );
    for ( int x=filterOffset; x<Data.getNumX ( )-filterOffset ; ++x ) {
      for ( int y=filterOffset; y<Data.getNumY ( )-filterOffset ; ++ y) {
        Data.copyBlockTo ( x-filterOffset, y-filterOffset, block );
        preprocessedData.set ( x-filterOffset, y-filterOffset, block.getMeanValue ( ) );
      }
    }
    
    if ( _verbose && _outputDir.size ( ) > 0 ) {
      std::stringstream ss;
      ss << _outputDir << "/preprocessedData" << qc::getDefaultArraySuffix ( qc::QC_2D );
      preprocessedData.save ( ss.str ( ).c_str ( ), qc::PGM_DOUBLE_BINARY );
    }
    
    for ( int i=0; i<2 ; ++i ) {
      // Compute energies of difference between image and shifted image for a range of possible periods
      aol::Vector<RealType> energies;
      for ( int period=0; period<0.5*aol::Min<int> ( Data.getNumX ( ), Data.getNumY ( ) ) ; ++period ) {
        RealType energy = 0;
        int numPoints = 0;
        aol::Vec2<RealType> shiftedPos;
        for ( int y=0; y<preprocessedData.getNumY ( ) ; ++y ) {
          for ( int x=0; x<preprocessedData.getNumX ( ) ; ++x ) {
            shiftedPos.set ( x + period * cos ( latticeAnglesRadians[i] ), y + period * sin ( latticeAnglesRadians[i] ) );
            if ( shiftedPos[0] >= 0 && shiftedPos[0] < preprocessedData.getNumX ( ) && shiftedPos[1] >= 0 && shiftedPos[1] < preprocessedData.getNumY ( ) ) {
              energy += aol::Sqr<RealType> ( preprocessedData.get ( x, y ) - preprocessedData.interpolate ( shiftedPos ) );
              ++numPoints;
            }
          }
        }
        energies.pushBack ( energy / static_cast<RealType> ( numPoints ) );
      }
      
      if ( _verbose && _outputDir.size ( ) > 0 ) {
        std::stringstream ss;
        ss << _outputDir << "/energies_" << i;
        aol::Plotter<RealType> plotter;
        plotter.set_outfile_base_name ( ss.str ( ).c_str ( ) );
        aol::PlotDataFileHandler<RealType> plotHandler;
        std::vector<std::pair<RealType, RealType> > data;
        for ( int i=0; i<energies.size ( ) ; ++i ) data.push_back ( std::pair<RealType, RealType> ( i, energies[i] ) );
        plotHandler.generateFunctionPlot ( data );
        plotter.addPlotCommandsFromHandler ( plotHandler );
        plotter.genPlot( aol::GNUPLOT_PNG );
        
        ss << ".csv";
        std::ofstream txtFile ( ss.str ( ).c_str ( ) );
        for ( int i=0; i<energies.size ( ) ; ++i )
          txtFile << i << ", " << energies[i] << std::endl;
        txtFile.close ( );
      }
      
      // Perform periodicity analysis on energies
      // Step 1: Smooth energies
      const short filterSizeEnergies = 3, filterOffsetEnergies = ( filterSizeEnergies - 1 ) / 2;
      aol::Vector<RealType> preprocessedEnergies ( energies.size ( ) );
      for ( int x=filterOffsetEnergies; x<energies.size ( )-filterOffsetEnergies ; ++x ) {
        for ( int dx=-filterOffsetEnergies; dx<=filterOffsetEnergies ; ++dx )
          preprocessedEnergies[x] += energies[x+dx];
        preprocessedEnergies[x] /= static_cast<RealType> ( filterSizeEnergies );
      }
      
      // Step 2: Find all local minima
      aol::Vector<RealType> localMinima;
      aol::Vector<int> localMinimaSpacings;
      for ( int x=1; x<preprocessedEnergies.size ( )-1 ; ++x ) {
        if ( preprocessedEnergies[x] < preprocessedEnergies[x-1] && preprocessedEnergies[x] < preprocessedEnergies[x+1] ) {
          localMinima.pushBack ( preprocessedEnergies[x] );
          localMinimaSpacings.pushBack ( x );
        }
      }
      
      // Step 3: Determine local minimum that corresponds to smallest spacing and which belongs to the same class of local minima as the global minimum
      int minSpacing = energies.size ( );
      
      if ( NumClusters[i] > 1 ) {
        // Step 3.1: Cluster local Minima (currently using k-means and k hast to be set by the user)
        KMeansClusterer<RealType> kMeansClusterer;
        aol::Vector<RealType> clusters;
        aol::Vector<int> clusterLabels;
        kMeansClusterer.apply ( localMinima, clusters, NumClusters[i], clusterLabels );
        
        // Step 3.2: In the cluster with smallest mean (i.e. the one corresponding to the global minimum), identify element that corresponds to smallest spacing
        std::pair<int, RealType> minClusterIndVal = clusters.getMinIndexAndValue ( );
        for ( int j=0; j<localMinima.size ( ) ; ++j ) {
          if ( clusterLabels[j] == minClusterIndVal.first && localMinimaSpacings[j] < minSpacing )
            minSpacing = localMinimaSpacings[j];
        }
      } else {
        // Step 3.1: Identify local minimum that corresponds to smallest spacing
        for ( int j=0; j<localMinima.size ( ) ; ++j ) {
          if ( localMinimaSpacings[j] < minSpacing )
            minSpacing = localMinimaSpacings[j];
        }
      }
      
      latticePeriods[i] = minSpacing;
      if ( _verbose ) std::cerr << "Spacing #" << i+1 << ": " << latticePeriods[i] << std::endl;
    }
    
    /*
     * END: Periodicty optimization based on Energy minimization
     */
    
    for ( int i=0; i<2 ; ++i ) {
      LatticeVectors[i][0] = latticePeriods[i] * cos ( latticeAnglesRadians[i] );
      LatticeVectors[i][1] = latticePeriods[i] * sin ( latticeAnglesRadians[i] );
    }
    
    // Create some output for debugging/analysis (if requested)
    if ( _verbose && _outputDir.size ( ) > 0 ) {
      std::stringstream ss;
      ss << _outputDir << "/periodicityAxes.png";
      saveDataPeriodicityAxesImg ( ss.str ( ).c_str ( ), Data, latticeAnglesRadians );
      
      ss.str ( std::string ( ) ); // clear stringstream
      ss << _outputDir << "/periodicityPattern.png";
      saveDataPeriodicPatternImg ( ss.str ( ).c_str ( ), Data, LatticeVectors );
    }
  }
  
  void getLatticeVectorsByProjectiveStdDevAndPeriodicityEnergyMinimization ( aol::MultiVector<RealType> &LatticeVectors, const PictureType &Data, const aol::Vector<RealType> &Params ) {
    if ( Params.size ( ) == 2 ) {
      aol::Vec2<int> numClusters ( 2 );
      for ( int i=0; i<2 ; ++i ) numClusters[i] = static_cast<int> ( Params[i] );
      getLatticeVectorsByProjectiveStdDevAndPeriodicityEnergyMinimization ( LatticeVectors, Data, numClusters );
    } else throw aol::Exception ( "Periodicity energy minimization requires the expected number of local minima clusters along each periodicity axis!" );
  }
  
  RealType to360Degrees ( const RealType Degrees ) const {
    RealType res = Degrees;
    if ( res < 0 ) res += 360;
    if ( res >= 360 ) res -= 360;
    return res;
  }
  
  const aol::Vec2<RealType> fitGaussianToProjectiveStandardDeviationPeak ( const aol::Vector<RealType> &Data ) const {
    aol::Vector<RealType> data ( Data );
    std::pair<int, RealType> indVal = data.getMaxIndexAndValue ( );
    std::vector<std::pair<RealType, RealType> > dataPairs;
    aol::Vec2<int> angleDegreesMinMax;
    for ( int sign=-1; sign<=1 ; sign+=2 ) {
      int angleDegrees = indVal.first;
      while ( data[angleDegrees] > 0.75 * indVal.second ) angleDegrees += sign;
      angleDegreesMinMax[( sign + 1 ) / 2] = angleDegrees;
    }
    int angleDegreeAbsThreshold = aol::Min<int> ( indVal.first - angleDegreesMinMax[0], angleDegreesMinMax[1] - indVal.first );
    for ( int angleDegrees=indVal.first-angleDegreeAbsThreshold; angleDegrees<=indVal.first+angleDegreeAbsThreshold ; ++angleDegrees )
      dataPairs.push_back ( std::pair<RealType, RealType> ( angleDegrees, data[angleDegrees] ) );

    Gaussian1DTargetFunctional<RealType> F ( dataPairs );
    Gaussian1DTargetJacobian<RealType, aol::FullMatrix<RealType> > DF ( dataPairs );
    LevenbergMarquardtAlgorithm<RealType, aol::FullMatrix<RealType>, LinearRegressionQR<RealType> > levenbergMarquardtAlg ( dataPairs.size ( ),
                                                                                                                            F, DF, 50, 1, 0.2, 0.8, 1e-6, 1e-6, 1e-6, false );
    aol::Vector<RealType> arg ( 3 ), dest ( 3 );
    arg[0] = indVal.first;
    arg[1] = aol::Sqr<RealType> ( ( angleDegreesMinMax[1] - angleDegreesMinMax[0] ) / 6 );
    arg[2] = indVal.second;
    levenbergMarquardtAlg.apply ( arg, dest );
    
    if ( _verbose && _outputDir.size ( ) > 0 ) {
      std::stringstream ss;
      ss << _outputDir << "/peakData_" << indVal.first;
      aol::Plotter<RealType> plotter;
      plotter.set_outfile_base_name ( ss.str ( ).c_str ( ) );
      aol::PlotDataFileHandler<RealType> plotHandler;
      plotHandler.generateFunctionPlot ( dataPairs );
      plotter.addPlotCommandsFromHandler ( plotHandler );
      plotter.genPlot( aol::GNUPLOT_PNG );
    }
    
    return aol::Vec2<RealType> ( dest[0], dest[1] );
  }
  
  void getProjectedAverageIntensities ( aol::Vector<RealType> &ProjectedAverageIntensities, const PictureType &Data, const RealType AngleDegrees, const int L0 ) const {
    RealType angleRadians = AngleDegrees * aol::NumberTrait<RealType>::pi / 180.0;
    std::map<int, RealType> aDelta, nDelta;
    
    // Project image onto line with origin (0,0) and angle angleRadians
    for ( int y=0; y<Data.getNumY ( ) ; ++y ) {
      for ( int x=0; x<Data.getNumX ( ) ; ++x ) {
        const int p = floor ( getProjectedPosition ( x, y, angleRadians ) );
        nDelta[p] = nDelta[p] + 1;
        aDelta[p] = aDelta[p] + Data.get ( x, y );
      }
    }
    
    // Normalize and threshold bins, then convert projected average intensities to vector and compute standard deviation
    ProjectedAverageIntensities.resize ( 0 );
    for ( typename std::map<int, RealType>::iterator it=aDelta.begin ( ); it != aDelta.end ( ); ++it ) {
      if ( nDelta[it->first] >= L0 ) aDelta[it->first] = aDelta[it->first] / nDelta[it->first];
      else aDelta[it->first] = 0;
      ProjectedAverageIntensities.pushBack ( aDelta[it->first] );
    }
    RealType meanVal = ProjectedAverageIntensities.getMeanValue ( );
    for ( int i=0; i<ProjectedAverageIntensities.size ( ) ; ++i ) {
      if ( ProjectedAverageIntensities[i] == 0 ) ProjectedAverageIntensities[i] = meanVal;
    }
  }
  
  RealType getProjectedPosition ( const int X, const int Y, const RealType AngleRadians ) const {
    return X * cos ( AngleRadians ) + Y * sin ( AngleRadians );
  }
  
  
  
  void getLatticeVectorsByFourierPeakExtractionAndSineFit ( aol::MultiVector<RealType> &LatticeVectors, const PictureType &Data ) {
    // Calculate fourier power coefficients
    const RealType mean = Data.getMeanValue ( );
    qc::MultiArray<RealType, 2, 2> dataEliminatedMean ( Data.getNumX ( ), Data.getNumY ( ) );
    for ( short i=0; i<Data.getNumX ( ) ; ++i )
      for ( short j=0; j<Data.getNumY ( ) ; ++j )
        dataEliminatedMean[0].set ( i, j, Data.get ( i, j ) - mean );
    qc::ScalarArray<RealType, qc::QC_2D> modulus ( Data.getNumX ( ), Data.getNumY ( ) );
    qc::computeLogFFTModulus<RealType> ( dataEliminatedMean[0], modulus, 0, false );
    PictureType fourierPowerCoefficients;
    fourierPowerCoefficients.rotate90From ( modulus );
    fourierPowerCoefficients.scaleValuesTo01 ( );

    // Find fourier power peaks
    aol::MultiVector<RealType> peaks ( 2, 2 );
    PictureType fourierPowerPeaks ( fourierPowerCoefficients );
    qc::FastILexMapper<qc::QC_2D> mapper ( fourierPowerPeaks.getNumX ( ), fourierPowerPeaks.getNumY ( ) );
    std::pair<int, RealType> maxIndVal;
    aol::Vec2<RealType> peakPos, center ( Data.getNumX ( ) / 2, Data.getNumY ( ) / 2 );
    int peakIdx = 0;
    while ( peakIdx < 2 ) {
      maxIndVal = fourierPowerPeaks.getMaxIndexAndValue ( );
      for ( short dx=-1; dx<=1 ; ++dx )
        for ( short dy=-1; dy<=1 ; ++dy )
          fourierPowerPeaks.set ( mapper.splitGlobalIndex ( maxIndVal.first )[0] + dx, mapper.splitGlobalIndex ( maxIndVal.first )[1] + dy, 0 );
      peakPos.set ( mapper.splitGlobalIndex ( maxIndVal.first )[0], mapper.splitGlobalIndex ( maxIndVal.first )[1] + 1 ); // TODO: why is the center incorrect?
      peakPos -= center;
      bool angleTooSmall = false;
      for ( short k=0; k<peakIdx ; ++k ) {
        const RealType p1dotp2 = peakPos.dotProduct ( aol::Vec2<RealType> ( peaks[k][0], peaks[k][1] ) ) / ( peakPos.norm ( ) * peaks[k].norm ( ) );
        if ( p1dotp2 < -0.9 || p1dotp2 > 0.9 )
          angleTooSmall = true;
      }
      if ( !angleTooSmall ) {
        peaks[peakIdx][0] = peakPos[0]; peaks[peakIdx][1] = peakPos[1];
        fourierPowerPeaks[maxIndVal.first] = -1;
        ++peakIdx;
      }
    }
    fourierPowerPeaks.clamp ( -1, 0 );
    fourierPowerPeaks *= -1;

    // Calculate angles of main axes of periodicity from the peak positions
    aol::Vec2<RealType> latticeAngles;
    for ( short i=0; i<2 ; ++i )
      latticeAngles[i] = atan ( peaks[i][1] / peaks[i][0] );
    
    /*
     * BEGIN: Calculate periodicity spacings
     */
    aol::Vec2<RealType> latticePeriods;
    
    // Preprocess data (to make peak finding and sine fitting more robust)
    const short filterSize = ( Data.getMaxValue ( ) <= 11 ) ? 7 : 5, filterOffset = ( filterSize - 1 ) / 2;
    PictureType preprocessedData ( Data.getNumX ( ) - filterSize + 1, Data.getNumY ( ) - filterSize + 1 );
    PictureType block ( filterSize, filterSize );
    for ( int x=filterOffset; x<Data.getNumX ( )-filterOffset ; ++x ) {
      for ( int y=filterOffset; y<Data.getNumY ( )-filterOffset ; ++ y) {
        Data.copyBlockTo ( x-filterOffset, y-filterOffset, block );
        preprocessedData.set ( x-filterOffset, y-filterOffset, block.getMeanValue ( ) );
      }
    }
    
    if ( _verbose && _outputDir.size ( ) > 0 ) {
      std::stringstream ss;
      ss << _outputDir << "/preprocessedData.png";
      preprocessedData.setOverflowHandling ( aol::CLIP_THEN_SCALE, preprocessedData.getMinValue ( ), preprocessedData.getMaxValue ( ) );
      preprocessedData.savePNG ( ss.str ( ).c_str ( ) );
    }
    
    // Find brightest peak and see which of the two axes offers the larger intersection with the image
    const qc::FastILexMapper<qc::QC_2D> mapperReduced ( preprocessedData.getNumX ( ), preprocessedData.getNumY ( ) );
    aol::Vec2<short> peak ( mapperReduced.splitGlobalIndex ( preprocessedData.getMaxIndexAndValue ( ).first )[0],
                            mapperReduced.splitGlobalIndex ( preprocessedData.getMaxIndexAndValue ( ).first )[1] );
    aol::RandomAccessContainer<std::vector<std::pair<RealType, RealType> > > intensitiesContainer ( 2 );
    for ( int k=0; k<2 ; ++k ) setIntensitiesAlongAxis ( intensitiesContainer[k], preprocessedData, latticeAngles[k], peak );
    const short kFirst = ( intensitiesContainer[0].size ( ) > intensitiesContainer[1].size ( ) ) ? 0 : 1;

    // Get periodicity spacing along primary axis kFirst
    RealType energy;
    latticePeriods[kFirst] = getPeriodicitySpacing ( intensitiesContainer[kFirst], energy, 1 );
    
    // Move origin of secondary axis along primary axis from brightest peak in periodicity steps,
    // and extract intensities from where the intersection of secondary axis with the image is largest
    std::vector<std::pair<RealType, RealType> > intensities;
    setIntensitiesAlongAxis ( intensities, preprocessedData, latticeAngles[1-kFirst], peak );
    short direction = 1;
    aol::Vec2<RealType> pos ( peak[0], peak[1] );
    aol::Vec2<RealType> stepVector ( direction * latticePeriods[kFirst] * cos ( latticeAngles[kFirst] ),
                                     direction * latticePeriods[kFirst] * sin ( latticeAngles[kFirst] ) );
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
      setIntensitiesAlongAxis ( intensities, preprocessedData, latticeAngles[1-kFirst], localMaxPos );
      if ( intensities.size ( ) > maxIntensitiesSize ) {
        maxIntensitiesSize = intensities.size ( );
        maxIntersectionOrigin.set ( localMaxPos );
      }
    }
    setIntensitiesAlongAxis ( intensities, preprocessedData, latticeAngles[1-kFirst], maxIntersectionOrigin );
    
    // Get periodicity spacing along secondary axis 1-kFirst
    latticePeriods[1-kFirst] = getPeriodicitySpacing ( intensities, energy, 2 );
    
    /*
     * END: Calculate periodicity spacings
     */
    
    for ( int i=0; i<2 ; ++i ) {
      LatticeVectors[i][0] = latticePeriods[i] * cos ( latticeAngles[i] );
      LatticeVectors[i][1] = latticePeriods[i] * sin ( latticeAngles[i] );
    }
    
    
    // Create some output for debugging/analysis (if requested)
    if ( _verbose && _outputDir.size ( ) > 0 ) {
      std::stringstream ss;
      ss << _outputDir << "/periodicityFourierPeaks.png";
      saveFourierCoefficients ( ss.str ( ).c_str ( ), fourierPowerCoefficients, true,
                                true, center, peaks,
                                true, latticeAngles );
      
      ss.str ( std::string ( ) ); // clear stringstream
      ss << _outputDir << "/periodicityAxes.png";
      saveDataPeriodicityAxesImg ( ss.str ( ).c_str ( ), Data, latticeAngles );
      
      ss.str ( std::string ( ) ); // clear stringstream
      ss << _outputDir << "/periodicityPattern.png";
      saveDataPeriodicPatternImg ( ss.str ( ).c_str ( ), Data, LatticeVectors );
    }
    
    if ( _verbose )
      std::cerr << latticePeriods << std::endl;
    
    if ( _verbose && _outputDir.size ( ) > 0 ) {
      std::stringstream ss;
      ss << _outputDir << "/periodicityAnalysis.txt";
      ofstream txtFile;
      txtFile.open ( ss.str ( ).c_str ( ) );
      txtFile << "Estimated grid parameters" << std::endl;
      txtFile << std::endl;
      txtFile << "Delta x_1 = " << latticePeriods[0] << " pixels" << std::endl;
      txtFile << "Delta x_2 = " << latticePeriods[1] << " pixels" << std::endl;
      txtFile << std::endl;
      txtFile << "alpha_1 = " << latticeAngles[0] * 180.0 / aol::NumberTrait<RealType>::pi << " degrees" << std::endl;
      txtFile << "alpha_2 = " << latticeAngles[1] * 180.0 / aol::NumberTrait<RealType>::pi << " degrees" << std::endl;
      txtFile << "alpha_1 = " << latticeAngles[0] << " radians" << std::endl;
      txtFile << "alpha_2 = " << latticeAngles[1] << " radians" << std::endl;
      txtFile.close ( );
    }
  }
  
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
  
  
  
  void saveFourierCoefficients ( const char* Path, const PictureType &Coefficients, const bool LogScale = false,
                                 const bool RedPeaks = false, const aol::Vec2<RealType> &Center = aol::Vec2<RealType> ( -1, -1 ), const aol::MultiVector<RealType> &Peaks = aol::MultiVector<RealType> ( ),
                                 const bool BlueAxes = false, const aol::Vec2<RealType> &LatticeAngles = aol::Vec2<RealType> ( -1, -1 ) ) const {
    PictureType u ( Coefficients );
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
      for ( int k=0; k<Peaks.numComponents ( ) ; ++k ) {
        v[0].set ( Peaks[k][0] + Center[0], Peaks[k][1] + Center[1] - 1, uMax );
        v[1].set ( Peaks[k][0] + Center[0], Peaks[k][1] + Center[1] - 1, 0 );
        v[2].set ( Peaks[k][0] + Center[0], Peaks[k][1] + Center[1] - 1, 0 );
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
        for ( short i=-Coefficients.getNumX ( ); i<Coefficients.getNumX ( ) ; ++i ) {
          pos.set ( Center[0] + cos ( LatticeAngles[k] ) * i, Center[1] + sin ( LatticeAngles[k] ) * i );
          if ( pos[0] >= 0 && pos[0] < Coefficients.getNumX ( ) && pos[1] >= 0 && pos[1] < Coefficients.getNumY ( ) ) {
            v[0].set ( pos, 0 );
            v[1].set ( pos, 0 );
            v[2].set ( pos, uMax );
          }
        }
      }
      v.setOverflowHandling ( aol::CLIP_THEN_SCALE, v.getMinValue ( ), v.getMaxValue ( ) );
      v.savePNG ( Path );
    } else u.save ( Path, qc::PGM_DOUBLE_BINARY );
  }
  
  void saveFourierCoefficientsRedPeaks ( const char* Path, const PictureType Coefficients,
                                         const aol::Vec2<RealType> &Center, const aol::MultiVector<RealType> &Peaks,
                                         const bool LogScale = false ) const {
    aol::Vec2<RealType> latticeAngles;
    saveFourierCoefficients ( Path, Coefficients, LogScale, true, Center, Peaks, false, latticeAngles );
  }
  
  void saveFourierCoefficientsBlueAxes ( const char* Path, const PictureType Coefficients,
                                         const aol::Vec2<RealType> &LatticeAngles,
                                         const bool LogScale = false ) const {
    aol::Vec2<RealType> center;
    aol::MultiVector<RealType> peaks;
    saveFourierCoefficients ( Path, Coefficients, LogScale, false, center, peaks, true, LatticeAngles );
  }

  void saveDataPeriodicityAxesImg ( const char* Path, const PictureType &Data, const aol::Vec2<RealType> &LatticeAnglesRadians,
                                    const aol::Vec2<short> &Origin = aol::Vec2<short> ( -1, -1 ), const short Axis = -1 ) const {
    if ( Axis < -1 || Axis > 1 )
      throw aol::Exception ( "Argument \"Axis\" has to be between -1 and 1!", __FILE__, __LINE__ );
    
    aol::Vec2<short> origin ( Origin );
    if ( origin[0] < 0 || origin[0] >= Data.getNumX ( ) || origin[1] < 0 || origin[1] >= Data.getNumY ( ) ) {
      qc::FastILexMapper<qc::QC_2D> mapper ( Data.getNumX ( ), Data.getNumY ( ) );
      origin.set ( mapper.splitGlobalIndex ( Data.getMaxIndexAndValue ( ).first )[0], mapper.splitGlobalIndex ( Data.getMaxIndexAndValue ( ).first )[1] );
    }
    
    ColoredPictureType periodicityAxesImg ( Data.getNumX ( ), Data.getNumY ( ) );
    periodicityAxesImg[0] = Data;
    periodicityAxesImg[1] = Data;
    periodicityAxesImg[2] = Data;
    aol::Vec2<short> pos;
    aol::Vector<short> axes;
    if ( Axis == -1 ) {
      axes.pushBack ( 0 );
      axes.pushBack ( 1 );
    } else axes.pushBack ( Axis );
    for ( short k=0; k<axes.size ( ) ; ++k ) {
      for ( short i=-Data.getNumX ( ); i<Data.getNumX ( ) ; ++i ) {
        pos.set ( origin[0] + cos ( LatticeAnglesRadians[axes[k]] ) * i, origin[1] + sin ( LatticeAnglesRadians[axes[k]] ) * i );
        if ( pos[0] >= 0 && pos[0] < Data.getNumX ( ) && pos[1] >= 0 && pos[1] < Data.getNumY ( ) ) {
          periodicityAxesImg[0].set ( pos, 0 );
          periodicityAxesImg[1].set ( pos, 0 );
          periodicityAxesImg[2].set ( pos, Data.getMaxValue ( ) );
        }
      }
    }
    periodicityAxesImg.setOverflowHandling ( aol::CLIP_THEN_SCALE, Data.getMinValue ( ), Data.getMaxValue ( ) );
    periodicityAxesImg.savePNG ( Path );
  }
  
  void saveDataPeriodicPatternImg ( const char* Path, const PictureType &Data, const aol::MultiVector<RealType> &LatticeVectors,
                                    const aol::Vec2<short> &Origin = aol::Vec2<short> ( -1, -1 ), const short SearchWindowSize = 1 ) const {
    aol::Vec2<short> origin ( Origin );
    if ( origin[0] < 0 || origin[0] >= Data.getNumX ( ) || origin[1] < 0 || origin[1] >= Data.getNumY ( ) ) {
      qc::FastILexMapper<qc::QC_2D> mapper ( Data.getNumX ( ), Data.getNumY ( ) );
      origin.set ( mapper.splitGlobalIndex ( Data.getMaxIndexAndValue ( ).first )[0], mapper.splitGlobalIndex ( Data.getMaxIndexAndValue ( ).first )[1] );
    }
    
    const short searchWindowOffset = ( SearchWindowSize - 1 ) / 2;

    ColoredPictureType periodicityAxesImg ( Data.getNumX ( ), Data.getNumY ( ) );
    for ( short i=0; i<Data.getNumX ( ) ; ++i ) {
      for ( short j=0; j<Data.getNumY ( ) ; ++j ) {
        periodicityAxesImg[0].set ( i, j, Data.get ( i, j ) );
        periodicityAxesImg[1].set ( i, j, Data.get ( i, j ) );
        periodicityAxesImg[2].set ( i, j, Data.get ( i, j ) );
      }
    }
    
    aol::Vec2<short> pos1, pos2, pos3;
    pos1.set ( origin );
    for ( short i=-Data.getNumX ( ); i<Data.getNumX ( ) ; ++i ) {
      pos1.set ( origin[0] + i * LatticeVectors[0][0], origin[1] + i * LatticeVectors[0][1] );
      for ( short j=-Data.getNumX ( ); j<Data.getNumX ( ) ; ++j ) {
        pos2.set ( pos1[0] + j * LatticeVectors[1][0], pos1[1] + j * LatticeVectors[1][1] );
        for ( short dx=-searchWindowOffset; dx<=searchWindowOffset ; ++dx ) {
          for ( short dy=-searchWindowOffset; dy<=searchWindowOffset ; ++dy ) {
            pos3.set ( pos2[0] + dx, pos2[1] + dy );
            if ( pos3[0] >= 0 && pos3[0] < Data.getNumX ( ) && pos3[1] >= 0 && pos3[1] < Data.getNumY ( ) ) {
              periodicityAxesImg[0].set ( pos3, Data.getMaxValue ( ) );
              periodicityAxesImg[1].set ( pos3, 0 );
              periodicityAxesImg[2].set ( pos3, 0 );
            }
          }
        }
      }
    }
    periodicityAxesImg[0].set ( origin, 0 );
    periodicityAxesImg[1].set ( origin, Data.getMaxValue ( ) );
    periodicityAxesImg.setOverflowHandling ( aol::SCALE, 0, 255 );
    periodicityAxesImg.savePNG ( Path );
  }
  
  void saveIntensityPlotAlongPeriodicAxis ( const char* Path, const PictureType &Data, const aol::Vec2<RealType> &LatticeAnglesRadians,
                                            const short Axis, const aol::Vec2<short> &Origin = aol::Vec2<short> ( -1, -1 ) ) const {
    if ( Axis < 0 || Axis > 1 )
      throw aol::Exception ( "Argument \"Axis\" has to be between 0 and 1!", __FILE__, __LINE__ );
    
    aol::Vec2<short> origin ( Origin );
    if ( origin[0] < 0 || origin[0] >= Data.getNumX ( ) || origin[1] < 0 || origin[1] >= Data.getNumY ( ) ) {
      qc::FastILexMapper<qc::QC_2D> mapper ( Data.getNumX ( ), Data.getNumY ( ) );
      origin.set ( mapper.splitGlobalIndex ( Data.getMaxIndexAndValue ( ).first )[0], mapper.splitGlobalIndex ( Data.getMaxIndexAndValue ( ).first )[1] );
    }
    
    std::vector<std::pair<RealType, RealType> > intensities;
    aol::Vec2<short> pos;
    for ( short sign=-1; sign<=1 ; sign+=2 ) {
      int i = ( sign == -1 ) ? 1 : 0;
      pos.set ( origin[0], origin[1] );
      while ( pos[0] >= 0 && pos[0] < Data.getNumX ( ) && pos[1] >= 0 && pos[1] < Data.getNumY ( ) ) {
        intensities.insert ( ( sign == -1 ) ? intensities.begin ( ) : intensities.end ( ),
                             std::pair<RealType, RealType> ( sign * i * aol::Vec2<RealType> ( cos ( LatticeAnglesRadians[Axis] ),
                                                                                              sin ( LatticeAnglesRadians[Axis] ) ).norm ( ),
                                                                                              Data.get ( pos ) ) );
        ++i;
        pos.set ( origin[0] + sign * cos ( LatticeAnglesRadians[Axis] ) * i, origin[1] + sign * sin ( LatticeAnglesRadians[Axis] ) * i );
      }
    }
  
    aol::Plotter<RealType> plotter;
    plotter.set_outfile_base_name ( Path );
    aol::PlotDataFileHandler<RealType> plotHandler;
    plotHandler.generateFunctionPlot ( intensities );
    plotter.addPlotCommandsFromHandler ( plotHandler );
    plotter.genPlot( aol::GNUPLOT_PNG );
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
