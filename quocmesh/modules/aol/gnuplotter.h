#ifndef __GNUPLOTTER_H
#define __GNUPLOTTER_H

#include <aol.h>
#include <levelSetDrawer.h>
#include <vectorExtensions.h>

namespace aol {

enum PlotOutFileType {
  GNUPLOT_PS,
  GNUPLOT_EPS,
  GNUPLOT_PNG,
  GNUPLOT_GIF,
  GNUPLOT_LATEX,
  GNUPLOT_PDF
};

enum PlotStyle {
  GNUPLOT_LINE,
  GNUPLOT_VECTOR,
  GNUPLOT_3DLINE,
  GNUPLOT_BOXES
};

/**
 * Class to automatically generate gnuplot comaptible data files.
 * Can be used in conjunction with aol::Plotter.
 *
 * Deletes all files created upon destruction.
 *
 * \author Berkels
 */
template <typename RealType>
class PlotDataFileHandler {
  std::vector<std::string> _datafileFileNameVector;
  std::vector<PlotStyle> _datafilePlotStyleVector;

  void deleteDataFiles ( ) const {
    for ( unsigned int i = 0; i < _datafileFileNameVector.size(); i++ )
      remove ( _datafileFileNameVector[i].c_str() );
  }

public:
  PlotDataFileHandler() {}

  ~PlotDataFileHandler() {
    deleteDataFiles ();
  }

  void clear ( ) {
    deleteDataFiles();
    _datafileFileNameVector.clear();
    _datafilePlotStyleVector.clear();
  }

  void generateColoredVectorfieldData ( const qc::Array<RealType> &VectorfieldX,
                                        const qc::Array<RealType> &VectorfieldY,
                                        const std::vector<const qc::BitArray<qc::QC_2D>*> &MaskVector,
                                        const RealType Spacing = 0.1) {
    char vectorFieldTempFileName[1024];
    for ( unsigned int i = 0; i < MaskVector.size(); i++ ) {
      // The segment is defined, but empty. Don't try to draw it or gnuplot will complain.
      if ( (MaskVector[i] != NULL) && MaskVector[i]->numTrue() == 0 )
        continue;

      ofstream vectorFieldTempOutFile;
      sprintf ( vectorFieldTempFileName, "vectorField.datXXXXXX" );
      generateTemporaryFile ( vectorFieldTempFileName, vectorFieldTempOutFile );
      qc::WriteVectorFieldAsGnuplotFile<RealType>( vectorFieldTempOutFile, VectorfieldX, VectorfieldY, Spacing, MaskVector[i] );
      vectorFieldTempOutFile.close();
      std::string vectorFieldTempFileNameString = vectorFieldTempFileName;
      _datafileFileNameVector.push_back(vectorFieldTempFileName);
      _datafilePlotStyleVector.push_back(GNUPLOT_VECTOR);
    }
  }
  void generateColoredVectorfieldData ( const aol::MultiVector<RealType> &Vectorfield,
                                        const int SizeX,
                                        const int SizeY,
                                        const std::vector<const qc::BitArray<qc::QC_2D>*> &MaskVector,
                                        const RealType Spacing = 0.1) {
    qc::Array<RealType> vectorfieldX ( Vectorfield[0], SizeX, SizeY );
    qc::Array<RealType> vectorfieldY ( Vectorfield[1], SizeX, SizeY );
    generateColoredVectorfieldData ( vectorfieldX, vectorfieldY, MaskVector, Spacing );
  }
  void generateColoredVectorfieldData ( const aol::MultiVector<RealType> &Vectorfield,
                                        const qc::GridStructure &Grid,
                                        const std::vector<const qc::BitArray<qc::QC_2D>*> &MaskVector,
                                        const RealType Spacing = 0.1) {
    generateColoredVectorfieldData ( Vectorfield, Grid.getNumX(), Grid.getNumY(), MaskVector, Spacing );
  }
  //! Unfortunately the ConfiguratorType dummy argument is necessary as workaround for a GCC bug.
  template<typename ConfiguratorType>
  void generateIsolineData ( const aol::Vector<RealType> &levelsetFunction,
                             const typename ConfiguratorType::InitType &Grid,
                             const ConfiguratorType &,
                             const double Isovalue = 0. ) {
    // The isoline is empty. Don't try to draw it or gnuplot will complain.
    if ( (levelsetFunction.getMaxValue() <= Isovalue) || (levelsetFunction.getMinValue() >= Isovalue) )
      return;

    qc::LevelSetDrawer<ConfiguratorType> drawer ( Grid );
    qc::ScalarArray<RealType, qc::QC_2D> levelsetFunctionArray ( levelsetFunction, Grid );
    char tempFileName[1024];
    sprintf ( tempFileName, "isoline.datXXXXXX" );
    ofstream tempOutFile;
    generateTemporaryFile ( tempFileName, tempOutFile );
    drawer.drawToGnuplot ( levelsetFunctionArray, tempOutFile, Isovalue, true, true );
    tempOutFile.close();
    std::string tempFileNameString = tempFileName;
    _datafileFileNameVector.push_back(tempFileNameString);
    _datafilePlotStyleVector.push_back(GNUPLOT_LINE);
  }
  //! Unfortunately the ConfiguratorType dummy argument is necessary as workaround for a GCC bug.
  template<typename ConfiguratorType>
  void generateIsolineData ( const aol::MultiVector<RealType> &levelsetFunctions,
                             const typename ConfiguratorType::InitType &Grid,
                             const ConfiguratorType &ConfDummy,
                             const double Isovalue = 0. ) {
    for ( int i = 0; i < levelsetFunctions.numComponents(); i++ )
      generateIsolineData ( levelsetFunctions[i], Grid, ConfDummy, Isovalue );
  }
  template <typename PositionDataType, typename ValueDataType>
  void generateFunctionPlot ( const aol::Vector<PositionDataType> &FunctionPositions,
                              const aol::Vector<ValueDataType> &FunctionValues, const bool PlotAsBoxes = false ) {
    QUOC_ASSERT ( FunctionPositions.size() == FunctionValues.size() );
    char tempFileName[1024];
    sprintf ( tempFileName, "function.datXXXXXX" );
    ofstream tempOutFile;
    generateTemporaryFile ( tempFileName, tempOutFile );
    for ( int i = 0; i < FunctionPositions.size(); ++i )
      tempOutFile << FunctionPositions[i] << " " << FunctionValues[i] << endl;
    tempOutFile.close();
    std::string tempFileNameString = tempFileName;
    _datafileFileNameVector.push_back(tempFileNameString);
    if ( PlotAsBoxes )
      _datafilePlotStyleVector.push_back(GNUPLOT_BOXES);
    else
      _datafilePlotStyleVector.push_back(GNUPLOT_LINE);
  }
  void generateFunctionPlot ( const std::vector<std::pair<RealType, RealType> > &FunctionPositionsAndValues, const bool PlotAsBoxes = false ) {
    const int numPositions = FunctionPositionsAndValues.size();
    aol::Vector<RealType> functionPositions ( numPositions );
    aol::Vector<RealType> functionValues ( numPositions );
    for ( int i = 0; i < numPositions; ++i ) {
      functionPositions[i] = FunctionPositionsAndValues[i].first;
      functionValues[i] = FunctionPositionsAndValues[i].second;
    }
    generateFunctionPlot ( functionPositions, functionValues, PlotAsBoxes );
  }
  template <typename InputDataType>
  void generateFunctionPlot ( const aol::Vector<InputDataType> &FunctionValues, const RealType A = 0, const RealType B = 1, const bool PlotAsBoxes = false ) {
    const int numberOfPoints = FunctionValues.size();
    aol::Vector<RealType> positions ( numberOfPoints );
    for ( int i = 0; i < numberOfPoints; ++i ) {
      positions[i] = A + static_cast<RealType>(i)/(numberOfPoints-1) * ( B - A );
    }
    generateFunctionPlot(positions, FunctionValues, PlotAsBoxes);
  }
  /**
   * Plots Function and its derivative in the interval [A,B], discretized with
   * NumberOfPoints points. Function needs to have the member functions
   * evaluate and evaluateDerivative. For example useful to plot
   * aol::ZeroOneIntervalPenaltyFunction.
   */
  template<typename FunctionType>
  void generateFunctionAndDerivativePlot ( const FunctionType &Function,
                                           const RealType A,
                                           const RealType B,
                                           const int NumberOfPoints ) {
    QUOC_ASSERT ( NumberOfPoints > 0 );
    aol::Vector<RealType> functionValues ( NumberOfPoints );
    aol::Vector<RealType> derivativeValues ( NumberOfPoints );
    aol::Vector<RealType> positions ( NumberOfPoints );
    for ( int i = 0; i < NumberOfPoints; ++i ) {
      const RealType position = A + static_cast<RealType>(i)/(NumberOfPoints-1) * ( B - A );
      functionValues[i] = Function.evaluate ( position );
      derivativeValues[i] = Function.evaluateDerivative ( position );
      positions[i] = position;
    }
    generateFunctionPlot(positions, functionValues);
    generateFunctionPlot(positions, derivativeValues);
  }

  void generateCurvePlot ( const aol::RandomAccessContainer<aol::Vec<2, RealType> > &Points ) {
    char tempFileName[1024];
    sprintf ( tempFileName, "curve.datXXXXXX" );
    ofstream tempOutFile;
    generateTemporaryFile ( tempFileName, tempOutFile );
    for ( int i = 0; i < Points.size()-1; ++i ) {
      tempOutFile << Points[i][0] << " " << Points[i][1] << endl;
      tempOutFile << Points[i+1][0] << " " << Points[i+1][1] << endl << endl;
    }
    tempOutFile.close();
    std::string tempFileNameString = tempFileName;
    _datafileFileNameVector.push_back(tempFileNameString);
    _datafilePlotStyleVector.push_back(GNUPLOT_LINE);
  }

  void generateClosedCurvePlot ( const aol::MultiVector<RealType> &Points ) {
    char tempFileName[1024];
    sprintf ( tempFileName, "curve.datXXXXXX" );
    ofstream tempOutFile;
    generateTemporaryFile ( tempFileName, tempOutFile );
    const int numPoints = Points.getEqualComponentSize();
    for ( int i = 0; i < numPoints; ++i ) {
      tempOutFile << Points[0][i] << " " << Points[1][i] << endl;
      const int next = (i+1)%numPoints;
      tempOutFile << Points[0][next] << " " << Points[1][next] << endl << endl;
    }
    tempOutFile.close();
    std::string tempFileNameString = tempFileName;
    _datafileFileNameVector.push_back(tempFileNameString);
    _datafilePlotStyleVector.push_back(GNUPLOT_LINE);
  }

  void generateCurvePlot ( const aol::RandomAccessContainer<aol::Vec<3, RealType> > & ) {
    throw aol::Exception ( "To plot 3D curves use generate3DCurvePlot and aol::CurvePlotter3D!", __FILE__, __LINE__ );
  }

  //! \note Not compatible with aol::Plotter since that class doesn't work with 3D data.
  void generate3DCurvePlot ( const aol::RandomAccessContainer<aol::Vec<3, RealType> > &Points ) {
    char tempFileName[1024];
    sprintf ( tempFileName, "curve.datXXXXXX" );
    ofstream tempOutFile;
    generateTemporaryFile ( tempFileName, tempOutFile );
    for ( int i = 0; i < Points.size(); ++i ) {
      tempOutFile << Points[i][0] << " " << Points[i][1] << " " << Points[i][2] << endl;
    }
    tempOutFile.close();
    std::string tempFileNameString = tempFileName;
    _datafileFileNameVector.push_back(tempFileNameString);
    _datafilePlotStyleVector.push_back(GNUPLOT_3DLINE);
  }

  void generateCenterLinePlot ( const qc::ScalarArray<RealType, qc::QC_2D> &Image ) {
    const int y = static_cast<int> ( Image.getNumY() / 2 );

    aol::Vector<RealType> functionValues ( Image.getNumX() );

    for ( int x = 0; x < Image.getNumX(); ++x )
      functionValues[x] = Image.get ( x, y );

    generateFunctionPlot ( functionValues, 0, Image.getNumX() - 1 );
  }
  
  void generateCrossSectionPlot ( const qc::ScalarArray<RealType, qc::QC_2D> &Image, aol::Vec2<short> &P1, aol::Vec2<short> &P2 ) {
    int numPoints = ceil ( aol::Vec2<RealType> ( P2[0] - P1[0], P2[1] - P1[1] ).norm ( ) ) + 1;
    aol::Vector<RealType> functionValues ( numPoints );
    RealType lambda;
    aol::Vec2<short> pos;
    for ( int i = 0; i < numPoints ; ++i ) {
      lambda = static_cast<RealType> ( i ) / ( numPoints - 1 );
      pos.set ( aol::Rint ( lambda * P1[0] + ( 1 - lambda ) * P2[0] ), aol::Rint ( lambda * P1[1] + ( 1 - lambda ) * P2[1] ) );
      functionValues[i] = Image.get ( pos );
    }
    
    generateFunctionPlot ( functionValues, 0, numPoints );
  }

  const std::vector<std::string> &getDataFileNames ( ) const {
    return _datafileFileNameVector;
  }
  const std::vector<PlotStyle> &getDataFileStyles ( ) const {
    return _datafilePlotStyleVector;
  }

};

//! Returns true if the system call to gnuplot returned EXIT_SUCCESS, false otherwise.
bool runGnuplot ( const char *GnuplotCommandFileName );

/**
 * Class to create various types of plots with gnuplot.
 *
 * Supports output to PS, EPS, PNG and GIF. Due to licensing
 * issues GIF output is not supported by the Linux version of
 * gnuplot.
 *
 * \author Berkels
 */
template <typename RealType>
class Plotter {
protected:
  char _noutput[1024], _plotfile[1024], _title[1024], _special[1024], _grid[1024], _xlabel[1024], _ylabel[1024], _plotcommand[1024], _outputDir[1024], _terminal[1024];
  string _xRange, _yRange;
  std::vector<std::string> _additionalPlotCommands;
public:
  Plotter ( ) {
    setDefaults();
  }
  //! if DataIsScatterplot = true the value for DataIsVectorField is ignored
  //! the combination DataIsScatterplot = true and DataIsVectorField = true is not defined
  Plotter ( const char *ndatafile, const bool DataIsVectorField = false, const bool DataIsScatterplot = false ) {
    setDefaults();
    sprintf ( _noutput, "%s", ndatafile );
    if ( !DataIsScatterplot ) {
      if( DataIsVectorField )
        sprintf ( _plotcommand, "plot \"%s\" w vec", ndatafile );
      else
        sprintf ( _plotcommand, "plot \"%s\" w l", ndatafile );
    } else {
      sprintf ( _plotcommand, "plot \"%s", ndatafile);
    }
  }
  Plotter ( const std::vector<std::string> &VectorDataFileNames ) {
    setDefaults();
    sprintf ( _noutput, "%s", VectorDataFileNames[0].c_str() );
    string plotcommand = "plot ";
    for( unsigned int i = 0; i < VectorDataFileNames.size(); i++ ){
      plotcommand += " \"";
      plotcommand += VectorDataFileNames[i];
      plotcommand += "\" w vec notitle";
      if( i < VectorDataFileNames.size()-1 )
        plotcommand += ",";
    }
    sprintf ( _plotcommand, "%s", plotcommand.c_str() );
  }
  Plotter ( const char *ndatafile1, const char *ndatafile2 ) {
    setDefaults();
    sprintf ( _noutput, "%s_%s", ndatafile1, ndatafile2 );
    sprintf ( _plotcommand, "plot \"%s\" w l, \"%s\" w l", ndatafile1, ndatafile2 );
  }
  Plotter ( const char *ndatafile1, const char *ndatafile2, const char *nameGraph1, const char *nameGraph2, const char *styleGraph1 = "l", const char *styleGraph2 = "l" ) {
    setDefaults();
    sprintf ( _noutput, "%s_%s", ndatafile1, ndatafile2 );
    sprintf ( _plotcommand, "plot \"%s\" title \"%s\" w %s, \"%s\" title \"%s\" w %s", ndatafile1, nameGraph1, styleGraph1, ndatafile2, nameGraph2, styleGraph2 );
  }
  ~Plotter () {}
protected:
  void setDefaults() {
    _title[0] = '\0';
    _special[0] = '\0';
    _grid[0] = '\0';
    _xlabel[0] = '\0';
    _ylabel[0] = '\0';
    _outputDir[0] = '\0';
    _plotcommand[0] = '\0';
    sprintf ( _noutput, "RENAMEME" );
  }
private:
  void gen() const {
    char gnuplotDatTempFileName[1024];
    sprintf ( gnuplotDatTempFileName, "gnuplot.datXXXXXX" );
    ofstream gnuplotdat;
    generateTemporaryFile ( gnuplotDatTempFileName, gnuplotdat );
    gnuplotdat << _title;
    gnuplotdat << _terminal;
    gnuplotdat << _special;
    gnuplotdat << _grid;
    gnuplotdat << "set xlabel \"" << _xlabel << "\"\n";
    ;
    gnuplotdat << "set ylabel \"" << _ylabel << "\"\n";
    ;
    gnuplotdat << _xRange;
    gnuplotdat << _yRange;
    gnuplotdat << _plotfile;
    // The different plot commands have to be separated by commas,
    // but there must not be a comma before the first one.
    if ( _plotcommand[0] != '\0' ){
      gnuplotdat << _plotcommand;
      if ( _additionalPlotCommands.size() > 0 )
        gnuplotdat << "," << _additionalPlotCommands[0];
    }
    else{
      QUOC_ASSERT ( _additionalPlotCommands.size() > 0 );
      gnuplotdat << "plot " << _additionalPlotCommands[0];
    }
    for ( unsigned int i = 1; i < _additionalPlotCommands.size(); i++ )
      gnuplotdat << "," << _additionalPlotCommands[i];

    gnuplotdat.close();
    runGnuplot ( gnuplotDatTempFileName );
    remove ( gnuplotDatTempFileName );
  }
public:
  void clearAdditionalPlotCommands ( ) {
    _additionalPlotCommands.clear();
  }

  void addPlotCommand( const char *DataFileName, const PlotStyle Style, const int LineStyleNumber = -2, const char *Title = NULL ){
    string additionalPlotCommand = "\"";
    additionalPlotCommand += DataFileName;
    switch( Style ) {
    case GNUPLOT_LINE:
      {
        additionalPlotCommand += "\" w l";
      }
      break;
    case GNUPLOT_VECTOR:
      {
        additionalPlotCommand += "\" w vec";
      }
      break;
    case GNUPLOT_BOXES:
      {
        additionalPlotCommand += "\" w boxes";
      }
      break;
    default:
      {
        // Don't add an empty plot command.
        cerr << "Warning: Ignoring unsupported PlotStyle.\n";
        return;
      }
    }

    if ( Title ) {
      additionalPlotCommand+= " title \"" + string ( Title ) + "\"";
    }
    else
      additionalPlotCommand +=  " notitle";

    if ( LineStyleNumber != -2 ){
      additionalPlotCommand += " ls ";
      additionalPlotCommand += aol::intFormat(LineStyleNumber);
    }
    _additionalPlotCommands.push_back( additionalPlotCommand );
  }
  void addPlotCommandsFromHandler( const PlotDataFileHandler<RealType> &Handler ){
    const std::vector<std::string> &datafileFileNameVector = Handler.getDataFileNames();
    const std::vector<PlotStyle> &datafilePlotStyleVector = Handler.getDataFileStyles();
    for ( unsigned int i = 0; i < datafileFileNameVector.size(); i++ )
      addPlotCommand( datafileFileNameVector[i].c_str(), datafilePlotStyleVector[i] );
  }
  //! \todo Remove code duplication with addPointsPlot.
  void addLinePlot( const char *DataFileName, const char *nameGraph = NULL ){
    string additionalPlotCommand = "\"";
    additionalPlotCommand += DataFileName;
    if ( nameGraph != NULL )
      additionalPlotCommand += "\" w l title \"" + string(nameGraph)+"\"";
    else
      additionalPlotCommand += "\" w l notitle";
    _additionalPlotCommands.push_back( additionalPlotCommand );
  }
  //! \todo Remove code duplication with addLinePlot.
  void addPointsPlot( const char *DataFileName, const char *nameGraph = NULL ){
    string additionalPlotCommand = "\"";
    additionalPlotCommand += DataFileName;
    if ( nameGraph != NULL )
      additionalPlotCommand += "\" w p title \"" + string(nameGraph)+"\"";
    else
      additionalPlotCommand += "\" w p notitle";
    _additionalPlotCommands.push_back( additionalPlotCommand );
  }
  void addGrid( const RealType xGridSize = -1, const RealType yGridSize = -1, const unsigned short int LineType = 0, const unsigned short int LineWidth = 1 ) {
    string gridInfo;
    if ( xGridSize > aol::ZTrait<RealType>::zero ) {
      gridInfo += "set xtics " + aol::to_string(xGridSize) + "\n";
    }
    if ( yGridSize > aol::ZTrait<RealType>::zero ) {
      gridInfo += "set ytics " + aol::to_string(yGridSize) + "\n";
    }
    gridInfo += "set grid xtics ytics lw " + aol::to_string(LineWidth) + " lt " + aol::to_string(LineType) + "\n\0";
    sprintf ( _grid, "%s", gridInfo.c_str() );
  }
  void set_outfile_base_name ( const char outfile_base_name[1024] ) {
    sprintf ( _noutput, "%s", outfile_base_name );
  }
  void set_outfile_base_name ( const string outfile_base_name ) {
    sprintf ( _noutput, "%s", outfile_base_name.c_str() );
  }
  void set_title ( const char title[1024] ) {
    sprintf ( _title, "set title \"%s\"\n", title );
  }
  void set_title ( const char title_line1[1024], const char title_line2[1024] ) {
    sprintf ( _title, "set title \"%s\\n%s\"\n", title_line1, title_line2 );
  }
  void setXLabel ( const char XLabel[1024] ){
    sprintf( _xlabel, "%s", XLabel );
  }
  void setYLabel ( const char YLabel[1024] ){
    sprintf( _ylabel, "%s", YLabel );
  }
  void setLabels ( const char XLabel[1024], const char YLabel[1024] ){
    setXLabel ( XLabel );
    setYLabel ( YLabel );
  }
  void setXRange ( const double XMin, const double XMax ) {
    _xRange = aol::strprintf ( "set xrange [%f:%f]\n", XMin, XMax );
  }
  void setYRange ( const double YMin, const double YMax ) {
    _yRange = aol::strprintf ( "set yrange [%f:%f]\n", YMin, YMax );
  }
  void setSpecial ( const char *Special ){
    sprintf( _special, "%s", Special );
  }
  void genPS() {
    sprintf ( _terminal, "set terminal postscript color solid\n" );
    sprintf ( _plotfile, "set output \"%s%s.ps\"\n", _outputDir, _noutput );
    gen();
  }
  void genEPS() {
    sprintf ( _terminal, "set terminal postscript eps color solid\n" );
    sprintf ( _plotfile, "set output \"%s%s.eps\"\n", _outputDir, _noutput );
    gen();
  }
  void genGrayEPS() {
    sprintf ( _terminal, "set terminal postscript eps\n" );
    sprintf ( _plotfile, "set output \"%s%s.eps\"\n", _outputDir, _noutput );
    gen();
  }
  void genPNG() {
    sprintf ( _terminal, "set terminal png\n" );
    sprintf ( _plotfile, "set output \"%s%s.png\"\n", _outputDir, _noutput );
    gen();
  }
  void genGIF() {
#ifdef WIN32
    sprintf ( _terminal, "set terminal gif small size 640,480\n" );
    sprintf ( _plotfile, "set output \"%s%s.gif\"\n", _outputDir, _noutput );
    gen();
#else
    cerr << "Plot to .gif not implemented under Linux\n";
#endif
  }
  void genLATEX() {
    sprintf ( _terminal, "set terminal latex\n" );
    sprintf ( _plotfile, "set output \"%s%s.tex\"\n", _outputDir, _noutput );
    gen();
  }
  void genPDF() {
    sprintf ( _terminal, "set terminal pdf color solid\n" );
    sprintf ( _plotfile, "set output \"%s%s.pdf\"\n", _outputDir, _noutput );
    gen();
  }
  void plotToScreen() {
    // Under Windows (possibly also under Linux, not tested yet) gnuplot closes itself
    // and the plot window immediately, if it was started with "gnuplot COMMANDFILENMAE".
    // To keep the plot window we tell the terminal to persist. Unfortunately this way
    // the plot window and the gnuplot console window have to be closed manually.
#ifdef WIN32
    sprintf ( _terminal, "set terminal wxt persist\n" );
#else
    // The default terminal on all OSes should render the plot on the screen.
    _terminal[0] = '\0';
#endif
    _plotfile[0] = '\0';
    gen();
  }
  void genPlot( const PlotOutFileType OutType, string *BaseFileName = NULL) {
    switch( OutType ) {
    case GNUPLOT_PS:
      {
        genPS();
        if( BaseFileName != NULL )
          *BaseFileName += ".ps";
      }
      break;
    case GNUPLOT_EPS:
      {
        genEPS();
        if( BaseFileName != NULL )
          *BaseFileName += ".eps";
      }
      break;
    case GNUPLOT_PNG:
      {
        genPNG();
        if( BaseFileName != NULL )
          *BaseFileName += ".png";
      }
      break;
    case GNUPLOT_GIF:
      {
        genGIF();
        if( BaseFileName != NULL )
          *BaseFileName += ".gif";
      }
      break;
    case GNUPLOT_LATEX:
      {
        genLATEX();
        if( BaseFileName != NULL )
          *BaseFileName += ".tex";
      }
      break;
    case GNUPLOT_PDF:
      {
        genPDF();
        if( BaseFileName != NULL )
          *BaseFileName += ".pdf";
      }
      break;
    default:
      throw aol::Exception ( "aol::Plotter::genPlot: unknown PlotOutFileType", __FILE__, __LINE__ );
    }
  }
};

/**
 * Plots a vector of 3D points as curve in 3D.
 *
 * \author Berkels
 */
 template <typename RealType>
class CurvePlotter3D : protected Plotter<RealType> {
  aol::PlotDataFileHandler<RealType> _plotDataFileHandler;

public:
  CurvePlotter3D ( aol::RandomAccessContainer<aol::Vec<3, RealType> > &Points ) {
    _plotDataFileHandler.generate3DCurvePlot ( Points );
    sprintf ( this->_plotcommand, "splot \"%s\" u 1:2:3 with lines notitle", _plotDataFileHandler.getDataFileNames()[0].c_str() );
  }

  using Plotter<RealType>::set_outfile_base_name;
  using Plotter<RealType>::genPlot;
};

/**
 * This class accompanies classes derivied from GradientDescentBase
 * and plots the energy descent in the course of the descent iteration.
 * Changes of the filter width are also plotted, dots on the energy
 * graph mark the iterations where the filter width is adapted.
 *
 * \author Berkels
 */
template<typename DescentType>
class EnergyDescentPlotter{
  std::ofstream _outEnergy, _outSigma;
  string _energyFilename, _sigmaFilename;
  string _saveDirectory;
public:
  EnergyDescentPlotter( const char* EnergyFileName,
                        const char* FilterWidthFileName,
                        const char *SaveDirectory,
                        DescentType &Descent ){
    _saveDirectory = SaveDirectory;
    _energyFilename = SaveDirectory;
    _energyFilename += EnergyFileName;
    _outEnergy.open( _energyFilename.c_str() );
    Descent.setOutStream( _outEnergy );
    _sigmaFilename = SaveDirectory;
    _sigmaFilename += FilterWidthFileName;
    _outSigma.open( _sigmaFilename.c_str() );
    Descent.setFilterWidthOutStream( _outSigma );
  }
  ~EnergyDescentPlotter (){
    _outEnergy.close();
    _outSigma.close();
  }
  void plot( const char *PlotBaseFileName ){
    aol::Plotter<double> plotter( _energyFilename.c_str(), _sigmaFilename.c_str(), "", "", "l", "p" );
    string plotFilename = _saveDirectory;
    plotFilename += PlotBaseFileName;
    plotter.set_outfile_base_name( plotFilename );
    plotter.setLabels( "iterations", "Energy" );
    plotter.genPlot( aol::GNUPLOT_PNG );
    plotter.genPlot( aol::GNUPLOT_PS );
    plotter.genPlot( aol::GNUPLOT_EPS );
  }
};

/**
 * This function plots a histogram (given as vector of integer values) as box plot.
 * Optionally it also includes the cumulative histogram (scaled to match the range
 * of the histogram) as line plot.
 *
 * \author Berkels
 */
template <typename RealType>
void plotHistogram ( const aol::Vector<int> &Histogram, const string &BaseOutName, const bool PlotCumulativeHisto = false, const bool PlotToScreen = false,
                     const bool SetYRangeFromZeroToMaxBinCount = false ) {
  aol::Plotter<RealType> plotter;
  plotter.set_outfile_base_name( BaseOutName );
  aol::PlotDataFileHandler<RealType> plotHandler;
  const int numVals = Histogram.size();
  plotHandler.generateFunctionPlot( Histogram, 0, numVals - 1, true );

  if ( PlotCumulativeHisto ) {
    aol::Vector<RealType> cumulativeHisto ( numVals );
    cumulativeHisto[0] = Histogram[0];
    for ( int i = 1; i < numVals; ++i )
      cumulativeHisto[i] = cumulativeHisto[i-1] + Histogram[i];
    cumulativeHisto *= Histogram.getMaxValue() / cumulativeHisto.getMaxValue();
    plotHandler.generateFunctionPlot( cumulativeHisto, 0, numVals - 1 );
  }

  plotter.addPlotCommandsFromHandler ( plotHandler );
  plotter.setXRange ( 0, numVals-1 );
  if ( SetYRangeFromZeroToMaxBinCount )
    plotter.setYRange ( 0, Histogram.getMaxValue ( ) * 1.15 );
  const string boxplotOptions = "set grid\nset boxwidth 0.95 relative\nset style fill transparent solid 0.5 noborder\n";
  plotter.setSpecial ( boxplotOptions.c_str() );
  if ( PlotToScreen )
    plotter.plotToScreen();
  else
    plotter.genPlot( aol::GNUPLOT_PDF );
}
  
/**
 * This function plots a histogram, given as vector of pairs of bin values (real)
 * and bin counts (integer) as box plot.
 * Optionally it also includes the cumulative histogram (scaled to match the range
 * of the histogram) as line plot.
 *
 * \author mevenkamp
 */
template <typename RealType>
  void plotHistogram ( const std::vector<std::pair<RealType, int> > &Histogram, const string &BaseOutName,
                       const bool PlotCumulativeHisto = false, const bool PlotToScreen = false,
                       const bool SetYRangeFromZeroToMaxBinCount = false ) {
  aol::Plotter<RealType> plotter;
  plotter.set_outfile_base_name( BaseOutName );
  aol::PlotDataFileHandler<RealType> plotHandler;
  const int numVals = Histogram.size();
  std::vector<std::pair<RealType, RealType> > histo;
  for ( int i=0; i < Histogram.size ( ) ; ++i )
    histo.push_back ( std::pair<RealType, RealType> ( Histogram[i].first, static_cast<RealType> ( Histogram[i].second ) ) );
  plotHandler.generateFunctionPlot ( histo, true );
    
    RealType histoMaxCounts = 0, histoMinBinValue = aol::NumberTrait<RealType>::Inf, histoMaxBinValue = 0;
  for ( int i=0; i<histo.size ( ) ; ++i ) {
    if ( histo[i].second > histoMaxCounts )
      histoMaxCounts = histo[i].second;
    if ( histo[i].first < histoMinBinValue )
      histoMinBinValue = histo[i].first;
    if ( histo[i].first > histoMaxBinValue )
      histoMaxBinValue = histo[i].first;
  }
  
  if ( PlotCumulativeHisto ) {
    std::vector<std::pair<RealType, RealType> > cumulativeHisto ( numVals );
    cumulativeHisto[0].first = Histogram[0].first;
    cumulativeHisto[0].second = Histogram[0].second;
    RealType cumulativeHistoMax = cumulativeHisto[0].second;
    for ( int i = 1; i < numVals; ++i ) {
      cumulativeHisto[i].first = Histogram[i].first;
      cumulativeHisto[i].second = cumulativeHisto[i-1].second + Histogram[i].second;
      if ( cumulativeHisto[i].second > cumulativeHistoMax )
        cumulativeHistoMax = cumulativeHisto[i].second;
    }
    for ( int i = 0; i < numVals; ++i )
      cumulativeHisto[i].second *= histoMaxCounts / cumulativeHistoMax;
    plotHandler.generateFunctionPlot( cumulativeHisto );
  }
  
  plotter.addPlotCommandsFromHandler ( plotHandler );
  plotter.setXRange ( histoMinBinValue, histoMaxBinValue );
  if ( SetYRangeFromZeroToMaxBinCount )
    plotter.setYRange ( 0, histoMaxCounts * 1.15 );
  const string boxplotOptions = "set grid\nset boxwidth 0.95 relative\nset style fill transparent solid 0.5 noborder\n";
  plotter.setSpecial ( boxplotOptions.c_str() );
  if ( PlotToScreen )
    plotter.plotToScreen();
  else
    plotter.genPlot( aol::GNUPLOT_PDF );
}

} // namespace aol
#endif //__GNUPLOTTER_H
