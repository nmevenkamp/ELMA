#ifndef EMIMGQUALITYQUANTIFIER_H_
#define EMIMGQUALITYQUANTIFIER_H_

// standard
#include <iostream>
#include <fstream>
#include <cmath>

// quocmesh
#include "atomFinder.h"

template <typename _RealType,
          typename _MatrixType = aol::FullMatrix<_RealType>,
          typename _LinearRegressionType = LinearRegressionQR<_RealType>,
          typename _ScalarPictureType = qc::ScalarArray<_RealType, qc::QC_2D>,
          typename _ColoredPictureType = qc::MultiArray<_RealType, qc::QC_2D, 3> >
class EMImgQualityQuantifier {
  typedef _RealType RealType;
  typedef _MatrixType MatrixType;
  typedef _LinearRegressionType LinearRegressionType;
  typedef _ScalarPictureType PictureType;
  typedef _ColoredPictureType ColoredPictureType;
protected:
  std::string _outputDir, _outputDirGT, _outputDirEstimate;
  bool _quietMode, _diskOutput;
public:
  EMImgQualityQuantifier ( const std::string &OutputDir = "", const bool Quiet = true )
    : _outputDir ( OutputDir ), _quietMode ( Quiet ), _diskOutput ( OutputDir != "" ) {
    if ( _diskOutput ) {
      std::stringstream ss;
      ss << OutputDir << "/gt";
      _outputDirGT = ss.str ( );
      aol::makeDirectory ( _outputDirGT.c_str ( ) );
    
      ss.str ( std::string ( ) ); // clear stringstream
      ss << OutputDir << "/estimate";
      _outputDirEstimate = ss.str ( );
      aol::makeDirectory ( _outputDirEstimate.c_str ( ) );
    }
  }
  
  const aol::Vec2<RealType> getPrecisions ( const aol::MultiVector<RealType> &Centers,
                                            const RealType PeriodX, const RealType PeriodY, const RealType PeriodDelta,
                                            const RealType AngleX, const RealType AngleY, const RealType AngleDelta,
                                            const char *GnuplotXMatches = NULL, const char *GnuplotYMatches = NULL ) const {
    aol::MultiVector<RealType> distances;
    getInterAtomicDistances ( distances, Centers, PeriodX, PeriodY, PeriodDelta, AngleX, AngleY, AngleDelta, GnuplotXMatches, GnuplotYMatches );
    return getPrecisions ( distances );
  }
  
  const aol::Vec2<RealType> getPrecisions ( const std::string &CentersCSVSrcPath,
                                            const RealType PeriodX, const RealType PeriodY, const RealType PeriodDelta,
                                            const RealType AngleX, const RealType AngleY, const RealType AngleDelta ) const {
    aol::MultiVector<RealType> centers;
    AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDir, _quietMode );
    atomFinder.readCentersFromCSV ( centers, CentersCSVSrcPath );
    return getPrecisions ( centers, PeriodX, PeriodY, PeriodDelta, AngleX, AngleY, AngleDelta );
  }
  
  const aol::Vec2<RealType> getPrecisions ( const aol::MultiVector<RealType> &Centers, const aol::ParameterParser &Parser,
                                            const char *GnuplotXMatches = NULL, const char *GnuplotYMatches = NULL ) const {
    RealType periodX, periodY, periodDelta, angleX, angleY, angleDelta;
    readPrecisionAnalysisParameters ( Parser, periodX, periodY, periodDelta, angleX, angleY, angleDelta );
    return getPrecisions ( Centers, periodX, periodY, periodDelta, angleX, angleY, angleDelta, GnuplotXMatches, GnuplotYMatches );
  }
  
  const aol::Vec2<RealType> getPrecisions ( const std::string &CentersCSVSrcPath, const aol::ParameterParser &Parser ) const {
    aol::MultiVector<RealType> centers;
    AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDir, _quietMode );
    atomFinder.readCentersFromCSV ( centers, CentersCSVSrcPath );
    return getPrecisions ( centers, Parser );
  }
  
  const aol::Vec2<RealType> getPrecisions ( const PictureType &Data, const aol::ParameterParser &Parser ) const {
    aol::MultiVector<RealType> centers;
    AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDir, _quietMode );
    atomFinder.getAtomPositions ( centers, Data, Parser );
    return getPrecisions ( centers, Parser );
  }
  
  RealType getFidelity ( const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate ) const;
  
  RealType getFidelity ( const std::string &CentersRefCSVSrcPath, const std::string &CentersEstimateCSVSrcPath ) {
    aol::MultiVector<RealType> centersRef, centersEstimate;
    AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDirGT, _quietMode );
    atomFinder.readCentersFromCSV ( centersRef, CentersRefCSVSrcPath );
    atomFinder.setOutputDir ( _outputDirEstimate );
    atomFinder.readCentersFromCSV ( centersEstimate, CentersEstimateCSVSrcPath );
    return getFidelity ( centersRef, centersEstimate );
  }
  
  RealType getFidelity ( const PictureType &Ref, const PictureType &Estimate, const aol::ParameterParser &Parser ) {
    aol::MultiVector<RealType> centersRef, centersEstimate;
    AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDirGT, _quietMode );
    atomFinder.getAtomPositions ( centersRef, Ref, Parser );
    atomFinder.setOutputDir ( _outputDirEstimate );
    atomFinder.getAtomPositions ( centersEstimate, Estimate, Parser );
    return getFidelity ( centersRef, centersEstimate );
  }
  
  int getNumCorrespondences ( const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate ) const;
  
  RealType getDetectionFraction ( const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate ) const {
    return static_cast<RealType> ( getNumCorrespondences ( CentersRef, CentersEstimate ) ) / CentersRef.numComponents ( );
  }
  
  RealType getDetectionFraction ( const std::string &CentersRefCSVSrcPath, const std::string &CentersEstimateCSVSrcPath ) {
    aol::MultiVector<RealType> centersRef, centersEstimate;
    AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDirGT, _quietMode );
    atomFinder.readCentersFromCSV ( centersRef, CentersRefCSVSrcPath );
    atomFinder.setOutputDir ( _outputDirEstimate );
    atomFinder.readCentersFromCSV ( centersEstimate, CentersEstimateCSVSrcPath );
    return getDetectionFraction ( centersRef, centersEstimate );
  }
  
  RealType getDetectionFraction ( const PictureType &Ref, const PictureType &Estimate, const aol::ParameterParser &Parser ) {
    aol::MultiVector<RealType> centersRef, centersEstimate;
    AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDirGT, _quietMode );
    atomFinder.getAtomPositions ( centersRef, Ref, Parser );
    atomFinder.setOutputDir ( _outputDirEstimate );
    atomFinder.getAtomPositions ( centersEstimate, Estimate, Parser );
    return getDetectionFraction ( centersRef, centersEstimate );
  }
  
  RealType getMisdetectionFraction ( const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate ) const {
    return 1 - static_cast<RealType> ( getNumCorrespondences ( CentersRef, CentersEstimate ) ) / CentersEstimate.numComponents ( );
  }
  
  RealType getMisdetectionFraction ( const std::string &CentersRefCSVSrcPath, const std::string &CentersEstimateCSVSrcPath ) {
    aol::MultiVector<RealType> centersRef, centersEstimate;
    AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDirGT, _quietMode );
    atomFinder.readCentersFromCSV ( centersRef, CentersRefCSVSrcPath );
    atomFinder.setOutputDir ( _outputDirEstimate );
    atomFinder.readCentersFromCSV ( centersEstimate, CentersEstimateCSVSrcPath );
    return getMisdetectionFraction ( centersRef, centersEstimate );
  }
  
  RealType getMisdetectionFraction ( const PictureType &Reference, const PictureType &Estimate, const aol::ParameterParser &Parser ) {
    aol::MultiVector<RealType> centersRef, centersEstimate;
    AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDirGT, _quietMode );
    atomFinder.getAtomPositions ( centersRef, Reference, Parser );
    atomFinder.setOutputDir ( _outputDirEstimate );
    atomFinder.getAtomPositions ( centersEstimate, Estimate, Parser );
    return getMisdetectionFraction ( centersRef, centersEstimate );
  }
  
  void saveStatistics ( const std::string &Path, const PictureType &Data, const aol::ParameterParser &Parser ) const {
    aol::MultiVector<RealType> centers;
    AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDir, _quietMode );
    atomFinder.getAtomPositions ( centers, Data, Parser );
    saveStatistics ( Path, aol::MultiVector<RealType> ( 0, 0 ), centers, Parser );
  }
  
  const std::string getStatistics ( const PictureType &Data, const aol::ParameterParser &Parser ) const {
    aol::MultiVector<RealType> centers;
    AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDir, _quietMode );
    atomFinder.getAtomPositions ( centers, Data, Parser );
    return getStatistics ( aol::MultiVector<RealType> ( 0, 0 ), centers, Parser );
  }
  
  void saveStatistics ( const std::string &Path, const aol::MultiVector<RealType> &Centers, const aol::ParameterParser &Parser ) const {
    saveStatistics ( Path, aol::MultiVector<RealType> ( 0 , 0 ), Centers, Parser );
  }
  
  const std::string getStatistics ( const aol::MultiVector<RealType> &Centers, const aol::ParameterParser &Parser ) const {
    return getStatistics ( aol::MultiVector<RealType> ( 0 , 0 ), Centers, Parser );
  }
  
  void saveStatistics ( const std::string &Path, const PictureType &Reference, const PictureType &Estimate, const aol::ParameterParser &Parser ) const {
    aol::MultiVector<RealType> centersRef, centersEstimate;
    AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDirGT, _quietMode );
    atomFinder.getAtomPositions ( centersRef, Reference, Parser );
    
    atomFinder.setOutputDir ( _outputDirEstimate );
    atomFinder.getAtomPositions ( centersEstimate, Estimate, Parser );
    
    saveStatistics ( Path, centersRef, centersEstimate, Parser );
  }
  
  const std::string getStatistics ( const PictureType &Reference, const PictureType &Estimate, const aol::ParameterParser &Parser ) const {
    aol::MultiVector<RealType> centersRef, centersEstimate;
    AtomFinder<RealType, MatrixType, LinearRegressionType, PictureType, ColoredPictureType> atomFinder ( _outputDirGT, _quietMode );
    atomFinder.getAtomPositions ( centersRef, Reference, Parser );
    atomFinder.setOutputDir ( _outputDirEstimate );
    atomFinder.getAtomPositions ( centersEstimate, Estimate, Parser );
    return getStatistics ( centersRef, centersEstimate, Parser );
  }
  
  void saveStatistics ( const std::string &Path, const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate, const aol::ParameterParser &Parser ) const {
    ofstream txtFile;
    txtFile.open ( Path.c_str ( ) );
    txtFile << getStatistics ( CentersRef, CentersEstimate, Parser );
    txtFile.close ( );
  }
  
  const std::string getStatistics ( const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate, const aol::ParameterParser &Parser ) const {
    RealType periodX, periodY, periodDelta, angleX, angleY, angleDelta;
    readPrecisionAnalysisParameters ( Parser, periodX, periodY, periodDelta, angleX, angleY, angleDelta );
    return getStatistics ( CentersRef, CentersEstimate, periodX, periodY, periodDelta, angleX, angleY, angleDelta );
  }
  
  void saveStatistics ( const std::string &Path, const aol::MultiVector<RealType> &Centers,
                        const RealType PeriodX = 0, const RealType PeriodY = 0, const RealType PeriodDelta = 0, const RealType AngleX = 0, const RealType AngleY = 0, const RealType AngleDelta = 0 ) const {
    saveStatistics ( Path, aol::MultiVector<RealType> ( 0, 0 ), Centers, PeriodX, PeriodY, PeriodDelta, AngleX, AngleY, AngleDelta );
  }
  
  const std::string getStatistics ( const aol::MultiVector<RealType> &Centers,
                                    const RealType PeriodX = 0, const RealType PeriodY = 0, const RealType PeriodDelta = 0, const RealType AngleX = 0, const RealType AngleY = 0, const RealType AngleDelta = 0 ) const {
    return getStatistics ( aol::MultiVector<RealType> ( 0, 0 ), Centers, PeriodX, PeriodY, PeriodDelta, AngleX, AngleY, AngleDelta );
  }
  
  void saveStatistics ( const std::string &Path, const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate,
                        const RealType PeriodX = 0, const RealType PeriodY = 0, const RealType PeriodDelta = 0, const RealType AngleX = 0, const RealType AngleY = 0, const RealType AngleDelta = 0 ) const {
    ofstream txtFile;
    txtFile.open ( Path.c_str ( ) );
    txtFile << getStatistics ( CentersRef, CentersEstimate, PeriodX, PeriodY, PeriodDelta, AngleX, AngleY, AngleDelta );
    txtFile.close ( );
  }
  
  const std::string getStatistics ( const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate,
                                    const RealType PeriodX = 0, const RealType PeriodY = 0, const RealType PeriodDelta = 0, const RealType AngleX = 0, const RealType AngleY = 0, const RealType AngleDelta = 0 ) const {
    aol::MultiVector<RealType> distances;
    aol::Vec2<RealType> precisions, meanDistances;
    aol::Vec2<int> numNeighbors;
    
    if ( PeriodX > 0 && PeriodY > 0 ) {
      getInterAtomicDistances ( distances, CentersEstimate, PeriodX, PeriodY, PeriodDelta, AngleX, AngleY, AngleDelta );
      precisions.set ( getPrecisions ( distances ) );
      meanDistances.set ( distances[0].getMeanValue ( ), distances[1].getMeanValue ( ) );
      numNeighbors.set ( distances[0].size ( ), distances[1].size ( ) );
    }
    
    std::stringstream ss;
    ss << "---------- General ------------" << std::endl;
    ss << "Number of atoms:                   " << CentersEstimate.numComponents ( ) << std::endl;
    if ( CentersRef.numComponents ( ) > 0 )
      ss << "Number of atoms (ref.):            " << CentersRef.numComponents ( ) << std::endl;
    
    if ( PeriodX > 0 && PeriodY > 0 ) {
      ss << std::endl << std::endl;
      ss << "---------- Precision ----------" << std::endl;
      ss << "Number of horizontal neighbors:    " << numNeighbors[0] << std::endl;
      ss << "Number of vertical neighbors:      " << numNeighbors[1] << std::endl;
      ss << "Mean horizontal atom separation:   " << meanDistances[0] << " px" << std::endl;
      ss << "Mean vertical atom separation:     " << meanDistances[1] << " px" << std::endl;
      ss << "Horizontal precision:              " << precisions[0] << " px" << std::endl;
      ss << "Vertical precision:                " << precisions[1] << " px" << std::endl;
      ss << "Total precision:                   " << precisions.norm ( ) << " px" << std::endl;
    }
    
    if ( CentersRef.numComponents ( ) > 0 ) {
      ss << std::endl << std::endl;
      ss << "---------- Fidelity  ----------" << std::endl;
      ss << "Detection fraction:                " << getDetectionFraction ( CentersRef, CentersEstimate ) << std::endl;
      ss << "Misdetection fraction:             " << getMisdetectionFraction ( CentersRef, CentersEstimate ) << std::endl;
      ss << "Fidelity:                          " << getFidelity ( CentersRef, CentersEstimate ) << " px" << std::endl;
    }
    
    return ss.str ( );
  }
  
  void setOutputDir ( const std::string &OutputDir ) {
    _outputDir = OutputDir;
    _diskOutput = OutputDir != "";
    
    if ( _diskOutput ) {
      std::stringstream ss;
      ss << OutputDir << "/gt";
      _outputDirGT = ss.str ( );
      aol::makeDirectory ( _outputDirGT.c_str ( ) );
      
      ss.str ( std::string ( ) ); // clear stringstream
      ss << OutputDir << "/estimate";
      _outputDirEstimate = ss.str ( );
      aol::makeDirectory ( _outputDirEstimate.c_str ( ) );
    }
  }
  
  void setQuietMode ( const bool Quiet = true ) {
    _quietMode = Quiet;
  }
private:
  void readPrecisionAnalysisParameters ( const aol::ParameterParser &Parser,
                                         RealType &PeriodX, RealType &PeriodY, RealType &PeriodDelta,
                                         RealType &AngleX, RealType &AngleY, RealType &AngleDelta ) const {
    PeriodX = Parser.getDoubleOrDefault ( "periodX", 0 );
    PeriodY = Parser.getDoubleOrDefault ( "periodY", 0 );
    PeriodDelta = Parser.getDoubleOrDefault ( "periodDelta", 0 );
    AngleX = Parser.getDoubleOrDefault ( "angleX", 0 );
    AngleY = Parser.getDoubleOrDefault ( "angleY", 0 );
    AngleDelta = Parser.getDoubleOrDefault ( "angleDelta", 0 );
  }
  
  const aol::Vec2<RealType> getPrecisions ( const aol::MultiVector<RealType> &Distances ) const {
    return aol::Vec2<RealType> ( Distances[0].getStdDev ( ), Distances[1].getStdDev ( ) );
  }
  
  void getInterAtomicDistances ( aol::MultiVector<RealType> &Distances, const aol::MultiVector<RealType> &Centers,
                                const RealType PeriodX, const RealType PeriodY, const RealType PeriodDelta,
                                const RealType AngleX, const RealType AngleY, const RealType AngleDelta,
                                const char *GnuplotXMatches = NULL, const char *GnuplotYMatches = NULL ) const;
  
  void getCorrespondences ( aol::Vector<short> &Correspondences, const aol::MultiVector<RealType> &CentersRef, const aol::MultiVector<RealType> &CentersEstimate ) const;
};

#endif /* EMIMGQUALITYQUANTIFIER_H_ */
