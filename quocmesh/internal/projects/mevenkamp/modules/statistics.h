#ifndef STATISTICS_H_
#define STATISTICS_H_

#include <quoc.h>
#include <gnuplotter.h>
#include <stats.h>

#ifdef USE_BOOST
#ifndef Q_MOC_RUN // See: https://bugreports.qt-project.org/browse/QTBUG-22829
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/special_functions/gamma.hpp>
#endif
#endif

static const double INV_SQRT_2 = 1 / sqrt ( 2 );
template <typename _RealType>
class NormalDistribution {
  typedef _RealType RealType;
protected:
  aol::RandomGenerator _randomGenerator;
public:
  NormalDistribution ( )
    : _randomGenerator ( ) {
    _randomGenerator.randomize ( );
  }
  
  static RealType PDF ( const RealType X, const RealType Mean = 0, const RealType Variance = 1 ) {
    const RealType inv_sigma_sqrt_2 = 1 / sqrt ( Variance ) * INV_SQRT_2;
    return inv_sigma_sqrt_2 * exp ( -pow ( ( X - Mean ) * inv_sigma_sqrt_2 , 2 ) );
  }
  
  static RealType CDF ( const RealType Z, const RealType Mean = 0, const RealType Variance = 1 ) {
    return 0.5 * ( 1 + aol::Erf ( ( Z - Mean ) / sqrt ( Variance ) * INV_SQRT_2 ) );
  }
  
  static RealType InverseCDF ( const RealType P, const RealType Mean, const RealType Variance ) {
#ifdef USE_BOOST
    return Mean + sqrt ( 2 * Variance ) * boost::math::erf_inv<RealType> ( 2 * P - 1 );
#else
    throw aol::Exception ( "Boost required! Compile with -DUSE_BOOST=1", __FILE__, __LINE__ );
#endif
  }
  
  void addNoise ( aol::Vector<RealType> &Vals, const RealType Variance ) {
    for ( int k=0; k<Vals.size ( ) ; ++ k )
      Vals[k] += getRandomObservation ( 0, Variance );
  }
  
  RealType getRandomObservation ( const RealType Mean, const RealType Variance ) {
    return _randomGenerator.normalrReal<RealType> ( Mean, sqrt ( Variance ) );
  }
};


template <typename _RealType>
class PoissonDistribution {
  typedef _RealType RealType;
protected:
  aol::RandomGenerator _randomGenerator;
public:
  PoissonDistribution ( )
    : _randomGenerator ( ) {
    _randomGenerator.randomize ( );
  }
  
  static RealType PDF ( const int K, const RealType Lambda ) {
#ifdef USE_BOOST
    boost::math::poisson_distribution<RealType> poissonDistr ( Lambda );
    return boost::math::pdf<RealType> ( poissonDistr, K );
#else
    throw aol::Exception ( "Boost required! Compile with -DUSE_BOOST=1", __FILE__, __LINE__ );
#endif
  }
  
  static RealType CDF ( const int K, const RealType Lambda ) {
#ifdef USE_BOOST
    boost::math::poisson_distribution<RealType> poissonDistr ( Lambda );
    return boost::math::cdf<RealType> ( poissonDistr, K );
#else
    throw aol::Exception ( "Boost required! Compile with -DUSE_BOOST=1", __FILE__, __LINE__ );
#endif
  }
  
  static int InverseCDF ( const RealType P, const RealType Lambda ) {
#ifdef USE_BOOST
    boost::math::poisson_distribution<RealType> poissonDistr ( Lambda );
    return boost::math::quantile<RealType> ( poissonDistr, P );
#else
    throw aol::Exception ( "Boost required! Compile with -DUSE_BOOST=1", __FILE__, __LINE__ );
#endif
  }
  
  void addNoise ( aol::Vector<RealType> &Vals ) {
    for ( int k=0; k<Vals.size ( ) ; ++ k )
      Vals[k] = getRandomObservation ( Vals[k] );
  }
  
  int getRandomObservation ( const RealType Lambda ) {
    if ( Lambda <= 0.0 ) throw aol::Exception ( "Mean value must be greater than zero!", __FILE__, __LINE__ );
    
    const int res = _randomGenerator.poissonrInt<RealType> ( Lambda );
    if ( res < 0 ) throw aol::Exception ( "Implemented Poisson RNG returned negative observation!", __FILE__, __LINE__ );
    
    return res;
  }
};


template <typename _RealType>
class CoPoissonDistribution {
  typedef _RealType RealType;
protected:
  aol::RandomGenerator _randomGenerator;
public:
  CoPoissonDistribution ( )
  : _randomGenerator ( ) {
    _randomGenerator.randomize ( );
  }
  
  static RealType PDF ( const int K, const RealType Lambda ) {
#ifdef USE_BOOST
    boost::math::poisson_distribution<RealType> poissonDistr ( Lambda );
    return boost::math::pdf<RealType> ( poissonDistr, K );
#else
    throw aol::Exception ( "Boost required! Compile with -DUSE_BOOST=1", __FILE__, __LINE__ );
#endif
  }
  
  static RealType CDF ( const int K, const RealType Lambda ) {
#ifdef USE_BOOST
    return boost::math::gamma_p ( K+1, Lambda );
#else
    throw aol::Exception ( "Boost required! Compile with -DUSE_BOOST=1", __FILE__, __LINE__ );
#endif
  }
  
  static int InverseCDF ( const RealType P, const int K ) {
#ifdef USE_BOOST
    return boost::math::gamma_p_inv ( K, P );
#else
    throw aol::Exception ( "Boost required! Compile with -DUSE_BOOST=1", __FILE__, __LINE__ );
#endif
  }
  
  void addNoise ( aol::Vector<RealType> &Vals ) {
    for ( int k=0; k<Vals.size ( ) ; ++ k )
      Vals[k] = getRandomObservation ( Vals[k] );
  }
  
  int getRandomObservation ( const int K ) {
    return InverseCDF ( K, _randomGenerator.rReal<RealType> ( ) );
  }
};


template <typename _RealType>
class DistributionTester {
  typedef _RealType RealType;
protected:
  const aol::Vector<RealType> &_observations;
  const std::string &_outputDir, _statisticsPath;
  const bool _verbose, _diskOutput;
public:
  DistributionTester ( const aol::Vector<RealType> &Observations, const std::string &OutputDir = "", const bool Verbose = false )
    : _observations ( Observations ), _outputDir ( OutputDir ), _statisticsPath ( getStatisticsPath ( OutputDir ) ),
      _verbose ( Verbose ), _diskOutput ( OutputDir.size ( ) > 0 ) { }
  
  virtual ~DistributionTester ( ) { }
  
  void setBinBounds ( aol::Vector<RealType> &BinBounds, const int NumBins ) const {
    RealType min = _observations.getMinValue ( ), max = _observations.getMaxValue ( );
    BinBounds.reallocate ( NumBins + 1 );
    for ( int i=0; i<=NumBins ; ++i )
      BinBounds[i] = min + static_cast<double> ( i ) / static_cast<double> ( NumBins ) * ( max - min );
  }
  
  void saveHistogramPlot ( const short NumBins, const bool PlotCDF = true ) const {
    std::stringstream ss;
    ss << _outputDir << "/histogram_observations";
    std::vector<std::pair<RealType, int> > histo;
    _observations.createHistogramOfValues ( histo, NumBins );
    aol::plotHistogram<RealType> ( histo, ss.str ( ).c_str ( ), PlotCDF );
  }
  
  RealType getChiSquare ( const int NumBins, aol::Vector<RealType> &CDFParams ) const {
    aol::Vector<int> histo;
    _observations.createHistogramOfValues ( histo, NumBins );
    
    aol::Vector<RealType> binBounds ( NumBins );
    setBinBounds ( binBounds, NumBins );
    RealType chiSqr = 0;
    for ( int m=0; m<NumBins ; ++m ) {
      const RealType expected = ( CDF ( binBounds[m+1], CDFParams ) - CDF ( binBounds[m], CDFParams ) ) * _observations.size ( );
      chiSqr += pow ( histo[m] - expected, 2 ) / expected;
    }
    
    if ( _diskOutput ) {
      std::ofstream txtFile;
      txtFile.open ( _statisticsPath.c_str ( ), std::ios_base::app );
      txtFile << "X^2: " << chiSqr << std::endl;
    }
    
    return chiSqr;
  }
  
  static RealType getChiSquareOfUniformDistribution ( const aol::Vector<RealType> &Observations, const int NumBins ) {
    aol::Vector<int> histo;
    Observations.createHistogramOfValues ( histo, NumBins );
    
    RealType chiSqr = 0, expected = Observations.size ( ) / NumBins;
    for ( short m=0; m<NumBins ; ++m )
      chiSqr += pow ( histo[m] - expected, 2 ) / expected;
    
    return chiSqr;
  }
  
  RealType getKolmogorovSmirnovDistance ( aol::Vector<RealType> &CDFParams ) const {
    aol::Vector<RealType> sortedObservations ( _observations );
    sortedObservations.sortValues ( );
    aol::Vector<RealType> diffs;
    for ( int k=0; k<sortedObservations.size ( ) ; ++k ) {
      diffs.pushBack ( abs ( static_cast<double> ( k ) / static_cast<double> ( sortedObservations.size ( ) ) - CDF ( sortedObservations[k], CDFParams ) ) );
      diffs.pushBack ( abs ( static_cast<double> ( k - 1 ) / static_cast<double> ( sortedObservations.size ( ) ) - CDF ( sortedObservations[k], CDFParams ) ) );
    }
    const RealType res = diffs.getMaxValue ( );
    
    if ( _diskOutput ) {
      std::ofstream txtFile;
      txtFile.open ( _statisticsPath.c_str ( ), std::ios_base::app );
      txtFile << "Kolmogorov-Smirnov distance: " << res << std::endl;
    }
    
    return res;
  }
  
  RealType getPValue ( const RealType Observed, const aol::Vector<RealType> &CDFParams ) const {
    return 2 * aol::Min ( CDF ( Observed, CDFParams ), 1 - CDF ( Observed, CDFParams ) );
  }
  
  RealType getPValue ( const int I, const aol::Vector<RealType> &CDFParams ) const {
    return getPValue ( _observations[I], CDFParams );
  }
  
  void getPValues ( aol::Vector<RealType> &PValues, const aol::Vector<RealType> &CDFParams ) const {
    PValues.reallocate ( _observations.size ( ) );
    for ( int i=0; i<_observations.size ( ) ; ++i )
      PValues[i] = getPValue ( i, CDFParams );
  }
  
  void getPValues ( aol::Vector<RealType> &PValues, const aol::MultiVector<RealType> &CDFParams ) const {
    if ( CDFParams.numComponents ( ) != _observations.size ( ) )
      throw aol::Exception ( "Dimensions of specified CDF parameters don't agree with number of observations", __FILE__, __LINE__ );
    
    PValues.reallocate ( _observations.size ( ) );
    for ( int i=0; i<_observations.size ( ) ; ++i )
      PValues[i] = getPValue ( i, CDFParams[i] );
  }
  
  RealType getChiSquareOfUniformDistributionOfPValues ( const int NumBins, const aol::Vector<RealType> &CDFParams ) const {
    aol::Vector<RealType> pValues;
    getPValues ( pValues, CDFParams );
    return getChiSquareOfUniformDistribution ( pValues, NumBins );
  }
  
  RealType getChiSquareOfUniformDistributionOfPValues ( const int NumBins, const aol::MultiVector<RealType> &CDFParams ) const {
    aol::Vector<RealType> pValues;
    getPValues ( pValues, CDFParams );
    return getChiSquareOfUniformDistribution ( pValues, NumBins );
  }
  
  int getPValuesTendency ( const aol::MultiVector<RealType> &CDFParams ) const {
    aol::Vector<RealType> pValues;
    getPValues ( pValues, CDFParams );
    int tendency = 0;
    for ( int k=0; k<pValues.size ( ) ; ++k ) {
      if ( pValues[k] > 0.5 )
        tendency++;
      if ( pValues[k] < 0.5 )
        tendency--;
    }
    return 0.5 * tendency;
  }
  
  static std::pair<int, int> getNumPValuesInSymmetricQuantiles ( const aol::Vector<RealType> &PValues, const RealType QuantileSize ) {
    std::pair<int, int> numPValuesInQuantiles ( 0 , 0 );
    for ( int k=0; k<PValues.size ( ) ; ++k ) {
      if ( PValues[k] < QuantileSize )
        numPValuesInQuantiles.first++;
      if ( PValues[k] > 1-QuantileSize )
        numPValuesInQuantiles.second++;
    }
    return numPValuesInQuantiles;
  }
  
  std::pair<int, int> getNumPValuesInSymmetricQuantiles ( const aol::MultiVector<RealType> &CDFParams, const RealType QuantileSize ) const {
    aol::Vector<RealType> pValues;
    getPValues ( pValues, CDFParams );
    return getNumPValuesSmallerLargerThanSymmetricThreshold ( pValues, QuantileSize );
  }
  
  static int getPValuesQuantileThreshold ( const int N, const RealType P, const RealType Alpha ) {
    const RealType mu = N * P, sigma = sqrt ( mu * ( 1 - P ) );
    int k = 0;
    while ( 1 - NormalDistribution<RealType>::CDF ( ( k - 1 - mu ) / sigma ) > Alpha )
      k++;
    return k;
  }
  
  static int getPValuesTendencyThreshold ( const int NumObserved, const RealType Significancy ) {
    int threshold = 0;
    while ( NormalDistribution<RealType>::CDF ( -threshold / ( 0.5 * sqrt ( NumObserved ) ) ) > Significancy )
      threshold++;
    return threshold;
  }
  
  static RealType getChiSqrOfOneMinusAlphaMMinusOneQuantil ( const RealType Alpha, const int M ) {
    return 0; // TODO
  }
  
  RealType getFisherNumber ( const aol::MultiVector<RealType> &CDFParams ) {
    RealType fisherNumber = 0;
    aol::Vector<RealType> pValues;
    getPValues ( pValues, CDFParams );
    for ( int k=0; k<pValues.size ( ) ; ++k )
      fisherNumber += log ( pValues[k] );
    return fisherNumber * -2;
  }
  
  virtual RealType CDF ( const RealType Z, const aol::Vector<RealType> &CDFParams ) const = 0;
private:
  std::string getStatisticsPath ( const std::string &OutputDir ) const {
    std::stringstream ss;
    ss << OutputDir << "/statistics_Distribution.txt";
    return ss.str ( );
  }
};
      
      
template <typename _RealType>
class NormalDistributionTester : public DistributionTester<_RealType> {
  typedef _RealType RealType;
public:
  NormalDistributionTester ( const aol::Vector<RealType> &Observations, const std::string &OutputDir = "", const bool Verbose = false )
    : DistributionTester<RealType> ( Observations, OutputDir, Verbose ) { }
  
  RealType CDF ( const RealType Z, const aol::Vector<RealType> &CDFParams ) const {
    if ( CDFParams.size ( ) != 0 && CDFParams.size ( ) != 2 )
      throw aol::Exception ( "Between 0 and 2 parameters of the CDF have to be passed (mean and standard deviation)!", __FILE__, __LINE__ );
    
    if ( CDFParams.size ( ) == 0 )
      return NormalDistribution<RealType>::CDF ( Z );
    else
      return NormalDistribution<RealType>::CDF ( Z, CDFParams[0], CDFParams[1] );
  }
};
      

template <typename _RealType>
class PoissonDistributionTester : public DistributionTester<_RealType> {
  typedef _RealType RealType;
public:
  PoissonDistributionTester ( const aol::Vector<RealType> &Observations, const std::string &OutputDir = "", const bool Verbose = false )
    : DistributionTester<RealType> ( Observations, OutputDir, Verbose ) { }
  
  RealType CDF ( const RealType Z, const aol::Vector<RealType> &CDFParams ) const {
    if ( CDFParams.size ( ) != 1 )
      throw aol::Exception ( "Exactly 1 parameter of the CDF has to be passed (mean = standard deviation)!", __FILE__, __LINE__ );
 
    return PoissonDistribution<RealType>::CDF ( Z, CDFParams[0] );
  }
};



enum TESTIMAGE_TYPE {
  CHECKERBOARD, PERIODIC_ATOMS, NORMAL_DISTRIBUTION, UNIFORM_DISTRIBUTION
};

template <typename _RealType, typename _PictureType, typename _DistributionTesterType = NormalDistributionTester<_RealType> >
class NoiseImageAnalyzer {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
protected:
  const PictureType &_observations, &_means;
  const int _numBins;
  const std::string &_outputDir;
  const bool _verbose, _diskOutput;
  const qc::FastILexMapper<qc::QC_2D> _mapper;
  _DistributionTesterType _distributionTester;
  aol::MultiVector<RealType> _CDFParams;
  aol::RandomGenerator _randomGenerator;
  
public:
  NoiseImageAnalyzer ( const PictureType &Observations, const PictureType &Means, const int NumBins,
                       const std::string &OutputDir = "", const bool Verbose = false )
    : _observations ( Observations ), _means ( Means ),
    _numBins ( NumBins ),
    _outputDir ( OutputDir ), _verbose ( Verbose ), _diskOutput ( OutputDir.size ( ) > 0 ),
    _mapper ( Observations.getNumX ( ), Observations.getNumY ( ) ),
    _distributionTester ( _observations, OutputDir, Verbose ),
    _randomGenerator ( ) {
    _randomGenerator.randomize ( );
  }
  
  NoiseImageAnalyzer ( const PictureType &Observations, const PictureType &Means, const aol::RandomAccessContainer<aol::Vec2<short> > &/*Pixels*/,
                       const int NumBins,
                       const std::string &OutputDir = "", const bool Verbose = false )
  : _observations ( Observations ), _means ( Means ),
  _numBins ( NumBins ),
  _outputDir ( OutputDir ), _verbose ( Verbose ), _diskOutput ( OutputDir.size ( ) > 0 ),
  _mapper ( Observations.getNumX ( ), Observations.getNumY ( ) ),
  _distributionTester ( _observations, OutputDir, Verbose ),
  _randomGenerator ( ) {
    _randomGenerator.randomize ( );
  }
  
  virtual ~NoiseImageAnalyzer ( ) { }
  
  virtual void saveStatistics ( ) const {
    std::stringstream ss;
    ss << _outputDir << "/statistics_Noise.txt";
    ofstream txtFile;
    txtFile.open ( ss.str ( ).c_str ( ) );
    txtFile << "Statistics about the observed and mean image" << std::endl;
    txtFile << std::endl;
    txtFile << "X^2 ( H0 = PValues ~ U(0,1) ) = " << getChiSquareOfUniformDistributionOfPValues ( ) << std::endl;
    txtFile.close ( );
  }
  
  void saveHistogramPlots ( const bool PlotCDF = true ) const {
    saveHistogramPlotObservations ( PlotCDF );
    saveHistogramPlotMeans ( PlotCDF );
    saveHistogramPlotResiduals ( PlotCDF );
  }
  
  void saveHistogramPlotObservations ( const bool PlotCDF = true ) const {
    _distributionTester.saveHistogramPlot ( _numBins, PlotCDF );
  }
  
  void saveHistogramPlotMeans ( const bool PlotCDF = true ) const {
    std::stringstream ss;
    ss << _outputDir << "/histogram_means";
    std::vector<std::pair<RealType, int> > histo;
    _means.createHistogramOfValues ( histo, _numBins );
    aol::plotHistogram<RealType> ( histo, ss.str ( ).c_str ( ), PlotCDF );
  }
  
  void saveHistogramPlotResiduals ( const bool PlotCDF = true ) const {
    PictureType residuals ( _observations );
    residuals -= _means;
    
    std::stringstream ss;
    ss << _outputDir << "/histogram_residuals";
    std::vector<std::pair<RealType, int> > histo;
    residuals.createHistogramOfValues ( histo, _numBins );
    aol::plotHistogram<RealType> ( histo, ss.str ( ).c_str ( ), PlotCDF );
  }
  
  RealType getPValue ( const int I ) const {
    return _distributionTester.getPValue ( I, _CDFParams[I] );
  }
  
  RealType getPValue ( const int X, const int Y ) const {
    const int I = _mapper.getGlobalIndex ( X, Y );
    return getPValue ( I );
  }
  
  RealType getPValue ( const aol::Vec2<short> &Position ) const {
    return getPValue ( Position[0], Position[1] );
  }
  
  void getPValues ( aol::Vector<RealType> &PValues ) const {
    _distributionTester.getPValues ( PValues, _CDFParams );
  }
  
  void getPValues ( aol::Vector<RealType> &PValues, const aol::MultiVector<short> &Positions ) const {
    PValues.resize ( Positions.numComponents ( ) );
    for ( int k=0; k<Positions.numComponents ( ) ; ++k )
      PValues[k] = getPValue ( aol::Vec2<short> ( Positions[k][0], Positions[k][1] ) );
  }
  
  void getPValues ( PictureType &PValues ) const {
    for ( int x=0; x<_observations.getNumX ( ) ; ++x )
      for ( int y=0; y<_observations.getNumY ( ) ; ++y )
        PValues.set ( x, y, getPValue ( x, y ) );
  }
  
  RealType getChiSquareOfUniformDistributionOfPValues ( ) const {
    return _distributionTester.getChiSquareOfUniformDistributionOfPValues ( _numBins, _CDFParams );
  }
  
  RealType getChiSquareOfUniformDistributionOfPValues ( const aol::MultiVector<short> &Positions ) const {
    aol::Vector<RealType> pValues;
    getPValues ( pValues, Positions );
    return _distributionTester.getChiSquareOfUniformDistribution ( pValues, _numBins );
  }
  
  int getPValuesTendency ( ) const {
    return _distributionTester.getPValuesTendency ( _CDFParams );
  }
  
  const std::pair<int, int> getNumPValuesInSymmetricQuantiles ( const RealType QuantileSize ) const {
    return _distributionTester.getNumPValuesInSymmetricQuantiles ( _CDFParams, QuantileSize );
  }
  
  const std::pair<int, int> getNumPValuesInSymmetricQuantiles ( const RealType QuantileSize, const aol::MultiVector<short> &Positions ) const {
    aol::Vector<RealType> pValues;
    getPValues ( pValues, Positions );
    return _distributionTester.getNumPValuesInSymmetricQuantiles ( pValues, QuantileSize );
  }
  
  int getPValuesQuantileThreshold ( const int N, const RealType P, const RealType Alpha ) const {
    return _distributionTester.getPValuesQuantileThreshold ( N, P, Alpha );
  }
  
  int getPValuesTendencyThreshold ( const int NumObserved, const RealType Significancy ) const {
    return _distributionTester.getPValuesSignificancyThreshold ( NumObserved, Significancy );
  }
  
  RealType getChiSqrOfOneMinusAlphaMMinusOneQuantil ( const RealType Alpha, const int M ) const {
    return _distributionTester.getChiSqrOfOneMinusAlphaMMinusOneQuantil ( Alpha, M );
  }
  
  RealType getFisherNumber ( ) {
    return _distributionTester.getFisherNumber ( _CDFParams );
  }
  
  void saveHistogramPlotPValues ( ) {
    aol::Vector<RealType> pValues;
    getPValues ( pValues );
    std::vector<std::pair<RealType, int> > histoPValues;
    pValues.createHistogramOfValues ( histoPValues, _numBins );
    std::stringstream ss;
    ss << _outputDir << "/histogram_pValues";
    aol::plotHistogram<RealType> ( histoPValues, ss.str ( ).c_str ( ), false, false, true );
  }
  
  void saveHistogramPlotPValues ( const aol::MultiVector<short> &Positions ) {
    aol::Vector<RealType> pValues;
    getPValues ( pValues, Positions );
    std::vector<std::pair<RealType, int> > histoPValues;
    pValues.createHistogramOfValues ( histoPValues, _numBins );
    std::stringstream ss;
    ss << _outputDir << "/histogram_pValues";
    aol::plotHistogram<RealType> ( histoPValues, ss.str ( ).c_str ( ), false, false, true );
  }
  
  void saveSortedPValuesPlot ( ) {
    aol::Vector<RealType> pValues;
    getPValues ( pValues );
    pValues.sortValues ( );
    
    std::stringstream ss;
    ss << _outputDir << "/pValues";
    aol::Plotter<RealType> plotter;
    plotter.set_outfile_base_name ( ss.str ( ) );
    aol::PlotDataFileHandler<RealType> plotHandler;
    plotHandler.generateFunctionPlot ( pValues );
    plotter.addPlotCommandsFromHandler ( plotHandler );
    plotter.genPlot( aol::GNUPLOT_PNG );
  }
  
  void savePValuesImage ( ) {
    PictureType pValues ( _observations.getNumX ( ), _observations.getNumY ( ) );
    getPValues ( pValues );
    
    std::stringstream ss;
    ss << _outputDir << "/pValues.pgm";
    pValues.setOverflowHandlingToCurrentValueRange ( );
    pValues.save ( ss.str ( ).c_str ( ), qc::PGM_UNSIGNED_CHAR_BINARY );
    
    ss.str ( std::string ( ) ); // clear stringstream
    ss << _outputDir << "/pValues" << getDefaultArraySuffix ( qc::QC_2D );
    pValues.save ( ss.str ( ).c_str ( ), qc::PGM_DOUBLE_BINARY );
  }
  
  void getRandomNoiseImage ( PictureType &Dest, const aol::MultiVector<RealType> &CDFParams ) {
    if ( Dest.size ( ) != CDFParams.numComponents ( ) )
      throw aol::Exception ( "Dimensions of passed image don't agree with size of specified CDF parameters!", __FILE__, __LINE__ );
    
    for ( int k=0; k<Dest.size ( ) ; ++k )
      Dest[k] = getRandomObservation ( CDFParams[k] );
  }
  
  virtual void getMeanAndNoisyTestImageNonStatic ( PictureType &Means, PictureType &Observations,
                                                   const TESTIMAGE_TYPE Example, const aol::Vector<RealType> &ExampleParams,
                                                   const aol::Vector<RealType> &NoiseParams ) {
    aol::RandomGenerator randomGenerator;
    randomGenerator.randomize ( );
    if ( Example == CHECKERBOARD ) {
      if ( ExampleParams.size ( ) != 3 )
        throw aol::Exception ( "Exactly 3 example parameters have to be specified (min observed, max observed, size of the squares)!", __FILE__, __LINE__ );
      
      for ( short i=0; i<Means.getNumX ( ) / ExampleParams[0] ; ++i ) {
        for ( short j=0; j<Means.getNumY ( ) / ExampleParams[0] ; ++j ) {
          RealType mean = randomGenerator.rReal ( ExampleParams[0], ExampleParams[1] );
          for ( short k=0; k<ExampleParams[2] ; ++k ) {
            for ( short l=0; l<ExampleParams[2] ; ++l ) {
              short x = k + i * ExampleParams[2];
              short y = l + j * ExampleParams[2];
              if ( x >= 0 && x < Means.getNumX ( ) && y >= 0 && y < Means.getNumY ( ) )
                Means.set ( x, y, mean );
            }
          }
        }
      }
    } else if ( Example == PERIODIC_ATOMS ) {
      if ( ExampleParams.size ( ) != 3 )
        throw aol::Exception ( "Exactly 3 example parameters have to be specified (min observed, max observed, size of the atoms)!", __FILE__, __LINE__ );
      
      const RealType q = aol::NumberTrait<long double>::pi * sqrt ( 3. );
      PictureType periodic ( Means.getNumX ( ), Means.getNumY ( ) );
      for ( short x=0; x<Means.getNumX ( ) ; ++x ) {
        for ( short y=0; y<Means.getNumY ( ) ; ++y ) {
          const RealType scaledX = ExampleParams[2] * x;
          const RealType scaledY = ExampleParams[2] * y;
          
          periodic.set ( x, y, cos ( q * scaledX ) * cos ( q * scaledY / sqrt ( 3. ) ) - 0.5 * cos ( 2 * q * scaledY / sqrt( 3. ) ) );
        }
      }
      const RealType periodicMin = periodic.getMinValue ( );
      for ( int k=0; k<periodic.size ( ) ; ++k )
        periodic[k] -= periodicMin;
      const RealType periodicMax = periodic.getMaxValue ( );
      for ( int k=0; k<Means.size ( ) ; ++k )
        Means[k] = ExampleParams[0] + ( 1 - periodic[k] / periodicMax ) * ( ExampleParams[1] - ExampleParams[0] );
    } else if ( Example == UNIFORM_DISTRIBUTION ) {
      if ( ExampleParams.size ( ) != 2 )
        throw aol::Exception ( "Exactly 2 example parameters have to be specified (min observed, max observed)!", __FILE__, __LINE__ );
      for ( int k=0; k<Means.size ( ) ; ++k )
        Means[k] = randomGenerator.rReal ( ExampleParams[0], ExampleParams[1] );
    } else if ( Example == NORMAL_DISTRIBUTION ) {
      throw aol::Exception ( "Exactly 2 example parameters have to be specified (mean, variance)!", __FILE__, __LINE__ );
      for ( short x=0; x<Means.getNumX ( ) ; ++x )
        for ( short y=0; y<Means.getNumY ( ) ; ++y )
          Means.set ( x, y, randomGenerator.normalrReal ( ExampleParams[0], sqrt ( ExampleParams[1] ) ) );
    } else
      throw aol::Exception ( "Could not recognize specified example!", __FILE__, __LINE__ );
    
    aol::MultiVector<RealType> CDFParams;
    setCDFParams ( CDFParams, Means, NoiseParams );
    getRandomNoiseImage ( Observations, CDFParams );
  }
private:
  static PictureType getResiduals ( const PictureType &Observations, const PictureType &Means ) {
    PictureType residuals ( Observations );
    residuals -= Means;
    return residuals;
  }
  
  virtual void setCDFParams ( aol::MultiVector<RealType> &CDFParams, const PictureType &Means, const PictureType &Observations,
                              const aol::RandomAccessContainer<aol::Vec2<short> > &Pixels ) const = 0;
  
  virtual void setCDFParams ( aol::MultiVector<RealType> &CDFParams, const PictureType &Means, const PictureType &Observations ) const = 0;
  
  virtual void setCDFParams ( aol::MultiVector<RealType> &CDFParams, const PictureType &Means, const aol::Vector<RealType> &NoiseParams ) const = 0;
  
  virtual RealType getRandomObservation ( const aol::Vector<RealType> &CDFParams ) = 0;
};



template <typename _RealType, typename _PictureType>
class GaussianNoiseImageAnalyzer : public NoiseImageAnalyzer<_RealType, _PictureType, NormalDistributionTester<_RealType> > {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
public:
  GaussianNoiseImageAnalyzer ( const PictureType &Observations, const PictureType &Means, const int NumBins,
                               const std::string &OutputDir = "", const bool Verbose = false )
    : NoiseImageAnalyzer<RealType, PictureType, NormalDistributionTester<RealType> > ( Observations, Means, NumBins, OutputDir, Verbose ) {
    setCDFParams ( this->_CDFParams, Means, Observations );
  }
  
  GaussianNoiseImageAnalyzer ( const PictureType &Observations, const PictureType &Means, const aol::RandomAccessContainer<aol::Vec2<short > > &Pixels,
                               const int NumBins,
                               const std::string &OutputDir = "", const bool Verbose = false )
  : NoiseImageAnalyzer<RealType, PictureType, NormalDistributionTester<RealType> > ( Observations, Means, NumBins, OutputDir, Verbose ) {
    setCDFParams ( this->_CDFParams, Means, Observations, Pixels );
  }
  
  RealType getRandomObservation ( const aol::Vector<RealType> &CDFParams ) {
    return this->_randomGenerator.normalrReal ( CDFParams[0], sqrt ( CDFParams[1] ) );
  }
  
  static void getMeanAndNoisyTestImage ( PictureType &Means, PictureType &Observations, const int NumBins,
                                         const TESTIMAGE_TYPE Example, const aol::Vector<RealType> &ExampleParams,
                                         const aol::Vector<RealType> &NoiseParams ) {
    PictureType dummy;
    GaussianNoiseImageAnalyzer<RealType, PictureType> imageAnalyzer ( dummy, dummy, NumBins );
    imageAnalyzer.getMeanAndNoisyTestImageNonStatic ( Means, Observations, Example, ExampleParams, NoiseParams );
  }
  
  RealType getChiSquareOfNormalDistributionOfMethodNoise ( ) const {
    aol::Vector<RealType> observations ( this->_observations.size ( ) );
    for ( int k=0; k<this->_observations.size ( ) ; ++k )
      observations[k] = this->_observations[k] - this->_means[k];
    aol::Vector<RealType> CDFParams ( 2 );
    CDFParams[0] = 0;
    CDFParams[1] = aol::Sqr<RealType> ( observations.getStdDev ( ) );
    NormalDistributionTester<RealType> normalDistributionTester ( observations );
    return normalDistributionTester.getChiSquare ( 100, CDFParams );
  }
  
  void saveStatistics ( ) const {
    std::stringstream ss;
    ss << this->_outputDir << "/statistics_Noise.txt";
    ofstream txtFile;
    txtFile.open ( ss.str ( ).c_str ( ) );
    txtFile << "Statistics about the observed and mean image" << std::endl;
    txtFile << std::endl;
    txtFile << "X^2 ( H0 = PValues ~ U(0,1) ) = " << this->getChiSquareOfUniformDistributionOfPValues ( ) << std::endl;
    txtFile << "X^2 ( MethodNoise ~ N(0,sigma^2) ) = " << getChiSquareOfNormalDistributionOfMethodNoise ( ) << std::endl;
    PictureType residuals ( this->_observations );
    residuals -= this->_means;
    txtFile << "sigma = " << residuals.getStdDev ( ) << std::endl;
    txtFile.close ( );
  }
private:
  void setCDFParams ( aol::MultiVector<RealType> &CDFParams, const PictureType &Means, const PictureType &Observations ) const {
    aol::Vector<RealType> noiseParams ( 1 );
    PictureType residuals ( Observations );
    residuals -= Means;
    noiseParams[0] = aol::Sqr<RealType> ( residuals.getStdDev ( ) );
//    noiseParams[0] = aol::Sqr<RealType> ( 13 );
    setCDFParams ( CDFParams, Means, noiseParams );
  }
  
  void setCDFParams ( aol::MultiVector<RealType> &CDFParams, const PictureType &Means, const PictureType &Observations,
                      const aol::RandomAccessContainer<aol::Vec2<short > > &Pixels ) const {
    aol::Vector<RealType> noiseParams ( 1 );
    aol::Vector<RealType> residuals ( Pixels.size ( ) );
    for ( int i=0; i<Pixels.size ( ) ; ++i )
      residuals[i] = Observations.get ( Pixels[i] ) - Means.get ( Pixels[i] );
    noiseParams[0] = aol::Sqr<RealType> ( residuals.getStdDev ( ) );
//    noiseParams[0] = aol::Sqr<RealType> ( 13 );
    setCDFParams ( CDFParams, Means, noiseParams );
  }
  
  void setCDFParams ( aol::MultiVector<RealType> &CDFParams, const PictureType &Means, const aol::Vector<RealType> &NoiseParams ) const {
    CDFParams.reallocate ( Means.size ( ), 2 );
    for ( int k=0; k<Means.size ( ) ; ++k ) {
      CDFParams[k][0] = Means[k];
      CDFParams[k][1] = NoiseParams[0];
    }
  }
};


template <typename _RealType, typename _PictureType>
class PoissonNoiseImageAnalyzer : public NoiseImageAnalyzer<_RealType, _PictureType, PoissonDistributionTester<_RealType> > {
  typedef _RealType RealType;
  typedef _PictureType PictureType;
public:
  PoissonNoiseImageAnalyzer ( const PictureType &Observations, const PictureType &Means, const int NumBins,
                              const std::string &OutputDir = "", const bool Verbose = false )
    : NoiseImageAnalyzer<RealType, PictureType, PoissonDistributionTester<RealType> > ( Observations, Means, NumBins, OutputDir, Verbose ) {
    setCDFParams ( this->_CDFParams, Means, Observations );
  }
  
  PoissonNoiseImageAnalyzer ( const PictureType &Observations, const PictureType &Means, const aol::RandomAccessContainer<aol::Vec2<short > > &/*Pixels*/,
                              const int NumBins,
                              const std::string &OutputDir = "", const bool Verbose = false )
    : NoiseImageAnalyzer<RealType, PictureType, PoissonDistributionTester<RealType> > ( Observations, Means, NumBins, OutputDir, Verbose ) {
    setCDFParams ( this->_CDFParams, Means, Observations );
  }
  
  RealType getRandomObservation ( const aol::Vector<RealType> &CDFParams ) {
    return this->_randomGenerator.poissonrInt ( CDFParams[0] );
  }
  
  static void getMeanAndNoisyTestImage ( PictureType &Means, PictureType &Observations, const int NumBins,
                                         const TESTIMAGE_TYPE Example, const aol::Vector<RealType> &ExampleParams,
                                         const aol::Vector<RealType> &NoiseParams ) {
    PictureType dummy;
    PoissonNoiseImageAnalyzer<RealType, PictureType> imageAnalyzer ( dummy, dummy, NumBins );
    imageAnalyzer.getMeanAndNoisyTestImageNonStatic ( Means, Observations, Example, ExampleParams, NoiseParams );
  }
private:
  void setCDFParams ( aol::MultiVector<RealType> &CDFParams, const PictureType &Means, const PictureType &/*Observations*/ ) const {
    setCDFParams ( CDFParams, Means, aol::Vector<RealType> ( ) );
  }
  
  void setCDFParams ( aol::MultiVector<RealType> &CDFParams, const PictureType &Means, const PictureType &/*Observations*/,
                      const aol::RandomAccessContainer<aol::Vec2<short> > &/*Pixels*/ ) const {
    setCDFParams ( CDFParams, Means, aol::Vector<RealType> ( ) );
  }
  
  void setCDFParams ( aol::MultiVector<RealType> &CDFParams, const PictureType &Means, const aol::Vector<RealType> &/*NoiseParams*/ ) const {
    CDFParams.reallocate ( Means.size ( ), 1 );
    for ( int k=0; k<Means.size ( ) ; ++k )
      CDFParams[k][0] = Means[k];
  }
};

#endif /* STATISTICS_H_ */
