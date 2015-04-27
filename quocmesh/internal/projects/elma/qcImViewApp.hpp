#ifndef __QCIMVIEWAPP_H
#define __QCIMVIEWAPP_H

#include <qtIncludes.h>
#include <previousFolderFileDialog.hpp>
#include <scalarArray.h>
#include <multiArray.h>
#include <gnuplotter.h>
#include <convolution.h>

#include "ui_qcImViewApp.h"

/**
 * \brief Simple QT example class that opens, displays and manipulates 2D Quoc array files, as well as PNGs and PGMs
 *
 * \author Berkels, Mevenkamp
 */
class QuocImageViewApp : public QMainWindow, protected Ui::QuocImageViewer {
  Q_OBJECT
protected:
  std::string _dirName;
  std::vector<std::string> _imageList;
  unsigned int _currentImageIndex;
  double _minVal;
  double _maxVal;
  qc::MultiArray<double, qc::QC_2D, 3> _image;
  int _scaleExponent;
  std::vector<QPoint> _selectedFrameCorners;
  std::string _defaultSaveFormat;
  
  const QString DEFAULT_DIR_KEY;
  
  QLabel *_sbarFileName;
  QLabel *_sbarFramePosSize;
  QLabel *_sbarMousePos;
  QLabel *_sbarFileNum;

  void loadImage ( const char *Filename ) {
    imageLabel->setText ( "" );
    if ( aol::fileNameEndsWith ( Filename, ".png" ) ) {
      // Be carefuly when loading the file. Possibly it can't be read or is not even an image.
      try {
        _image.loadPNG ( Filename );
      }
      catch ( aol::Exception &el ) {
        el.dump();
        imageLabel->setText( el.getMessage().c_str() );
        _image.reallocate ( 1, 1 );
      }
      _minVal = 0;
      _maxVal = 255;
    } else {
      // Be careful when loading the file. Possibly it can't be read or is not even an image.
      try {
        qc::ScalarArray<double, qc::QC_2D> scalarImage;
        scalarImage.load ( Filename );
        _image.reallocate ( scalarImage.getNumX ( ), scalarImage.getNumY ( ), 3 );
        for ( unsigned short c=0; c<3 ; ++c )
          _image[c] = scalarImage;
      } catch ( aol::Exception &el ) {
        el.dump();
        imageLabel->setText( el.getMessage().c_str() );
        _image.reallocate ( 1, 1 );
      }
      bool rangeIs8Bit = ( aol::fileNameEndsWith ( Filename, ".pgm" ) );
      _minVal = rangeIs8Bit ? 0 : _image.getMinValue();
      _maxVal = rangeIs8Bit ? 255 : _image.getMaxValue();
    }
    
    if ( imageLabel->text().isEmpty() ) {
      refreshPixmap ( );
      if ( actionAutoContrast->isChecked() )
        on_actionEnhanceContrast_triggered();
      emit imageLoaded ( );
    }
  }
  
  virtual void refreshPixmap ( ) {
    _selectedFrameCorners[0].setX ( 0 );
    _selectedFrameCorners[0].setY ( 0 );
    _selectedFrameCorners[1].setX ( _image[0].getNumX ( ) - 1 );
    _selectedFrameCorners[1].setY ( _image[0].getNumY ( ) - 1 );
    
    QImage qimage;
    quocMultiArrayToQImage ( qimage, _minVal, _maxVal );
    imageLabel->setPixmap ( QPixmap::fromImage(qimage) );
  }
  
  bool imageIsGrayScale ( ) {
    for ( int i=0; i<_image[0].size ( ) ; ++i )
      if ( _image[1][i] != _image[0][i] || _image[2][i] != _image[0][i] )
        return false;
    return true;
  }
  
  void quocMultiArrayToQImage ( QImage &Qimage, const double MinVal, const double MaxVal ) const {
    const double valScaleFactor = 255 / ( MaxVal - MinVal );
    Qimage = QImage ( _image[0].getNumX ( ), _image[0].getNumY ( ), QImage::Format_ARGB32_Premultiplied );
    for ( int j = 0; j < _image[0].getNumY ( ) ; ++j )
      for ( int i = 0; i < _image[0].getNumX ( ) ; ++i )
        Qimage.setPixel ( i, j, qRgb ( static_cast<unsigned char> ( aol::Clamp<double> ( valScaleFactor * ( _image[0].get ( i, j ) - MinVal ), 0, 255 ) ),
                                       static_cast<unsigned char> ( aol::Clamp<double> ( valScaleFactor * ( _image[1].get ( i, j ) - MinVal ), 0, 255 ) ),
                                       static_cast<unsigned char> ( aol::Clamp<double> ( valScaleFactor * ( _image[2].get ( i, j ) - MinVal ), 0, 255 ) ) ) );
  }
  
  enum RGBTOGRAY_METHOD {
    LUMINANCE,  // Gray = Y = 0.2126 R + 0.7152 G + 0.0722 B (this ) http://www.w3.org/Graphics/Color/sRGB
    AVERAGE,    // Gray = ( R + G + B ) / 3
    LIGHTNESS   // Gray = ( max(R, G, B) + min(R, G, B) ) / 2
  };
  
  void rgbToGray ( const qc::MultiArray<double, qc::QC_2D, 3> &ColorImage, qc::ScalarArray<double, qc::QC_2D> &GrayImage,
                   const RGBTOGRAY_METHOD RGBToGrayMethod = LUMINANCE ) {
    GrayImage.reallocate ( ColorImage[0].getNumX ( ), ColorImage[0].getNumY ( ) );
    for ( int i = 0; i < ColorImage[0].size ( ) ; ++i ) {
      aol::Vec3<double> color;
      for ( int c = 0; c < 3; ++c )
        color[c] = _image[c][i];
      if ( RGBToGrayMethod == LUMINANCE ) {
        GrayImage[i] = 0.2126 * color[0] + 0.7152 * color[1] + 0.0722 * color[2];
      } else if ( RGBToGrayMethod == AVERAGE ) {
        GrayImage[i] = color.getMeanValue ( );
      } else if ( RGBToGrayMethod == LIGHTNESS ) {
        GrayImage[i] = 0.5 * ( color.getMinValue ( ) + color.getMaxValue ( ) );
      }
    }
  }
public:
  QuocImageViewApp ( const char * InputFilename = NULL )
    : _scaleExponent ( 0 ),
      _selectedFrameCorners ( 2 ),
      _defaultSaveFormat ( ".q2bz" ),
      DEFAULT_DIR_KEY ( "default_dir" ),
      _sbarFileName ( new QLabel ),
      _sbarFramePosSize ( new QLabel ),
      _sbarMousePos ( new QLabel ),
      _sbarFileNum ( new QLabel( "0/0" ) ) {
    setupUi ( this );
    statusbar->addPermanentWidget( _sbarFileName, 1 );
    statusbar->addPermanentWidget ( _sbarFramePosSize );
    statusbar->addPermanentWidget( _sbarMousePos );
    statusbar->addPermanentWidget( _sbarFileNum );
    imageLabel->setAlignment ( Qt::AlignCenter );
    scrollArea->setBackgroundRole ( QPalette::Dark );
    scrollArea->setAlignment ( Qt::AlignCenter );
    scrollArea->setWidget ( imageLabel );
    scrollArea->setWidgetResizable ( true );
    connect ( actionAboutQt, SIGNAL ( triggered() ), qApp, SLOT ( aboutQt() ) );
    loadAndShowImage ( InputFilename );

// Workaround for a Qt5 bug under OS X: https://bugreports.qt-project.org/browse/QTBUG-38256
#if defined ( __APPLE__ ) && ( QT_VERSION >= 0x050000 )
    foreach ( QAction* a, menuFile->actions())
      QObject::connect ( new QShortcut(a->shortcut(), a->parentWidget()), SIGNAL(activated()), a, SLOT(trigger()) );
    foreach ( QAction* a, menuView->actions())
      QObject::connect ( new QShortcut(a->shortcut(), a->parentWidget()), SIGNAL(activated()), a, SLOT(trigger()) );
#endif
  }

  void loadAndShowImage ( const char * InputFilename ) {
    if ( InputFilename != NULL ) {
      loadImage ( InputFilename );
      const QString path = InputFilename;
      initDirectory ( path );
    }
  }

  ~QuocImageViewApp ( ) {
    delete _sbarFileName;
    delete _sbarFramePosSize;
    delete _sbarMousePos;
    delete _sbarFileNum;
  }

  bool hasImageLoaded ( ) const {
    return ( _image.getSize()[0] != 0 );
  }
protected:
  std::map<std::string,std::string> getSupportedFileFormats ( ) {
    std::map<std::string,std::string> supportedFileFormats;
    
    supportedFileFormats[".png"] = "Portable Network Graphics";
    supportedFileFormats[".dat.bz2"] = "Quocmesh 2D Array (deprecated)";
    supportedFileFormats[".q2bz"] = "Quocmesh 2D Array";
    supportedFileFormats[".pgm"] = "Portable Gray Map";
    supportedFileFormats[".tif"] = "Tagged Image File Format";
    supportedFileFormats[".tiff"] = "Tagged Image File Format";
    supportedFileFormats[".dm3"] = "Digital Micrograph 3";
    supportedFileFormats[".dm4"] = "Digital Micrograph 4";
    
    return supportedFileFormats;
  }
  
  std::map<std::string,std::string> getSupportedLoadFormats ( ) {
    std::map<std::string,std::string> supportedLoadFormats = getSupportedFileFormats ( );
    supportedLoadFormats.erase ( ".dat.bz2" );
    return supportedLoadFormats;
  }
  
  std::map<std::string,std::string> getSupportedSaveFormats ( ) {
    std::map<std::string,std::string> supportedSaveFormats = getSupportedFileFormats ( );
    supportedSaveFormats.erase ( ".dat.bz2" );
    supportedSaveFormats.erase ( ".tiff" );
    supportedSaveFormats.erase ( ".dm3" );
    supportedSaveFormats.erase ( ".dm4" );
    return supportedSaveFormats;
  }
  
  void initDirectory ( const QString Filename ) {
    _dirName = QFileInfo ( Filename ).absolutePath().toLatin1().constData();
    _dirName += "/";
    _imageList.clear();

    std::vector<std::string> dirList;
    aol::createDirectoryListing ( _dirName.c_str(), dirList );

    std::map<std::string,std::string> supportedFileFormats = getSupportedFileFormats ( );

    for ( unsigned int i = 0; i < dirList.size(); ++i ) {
      // Ignore OS X meta data files.
      if ( ( dirList[i].size() > 1 ) && ( dirList[i].compare ( 0, 2, "._" ) == 0 ) )
        continue;

      for ( std::map<std::string, std::string>::iterator it=supportedFileFormats.begin ( ); it != supportedFileFormats.end ( ) ; ++it ) {
        if ( aol::fileNameEndsWith ( dirList[i].c_str(), it->first.c_str() ) ) {
          _imageList.push_back ( dirList[i] );
          break;
        }
      }
    }

    _currentImageIndex = 0;

    // Make sure that the filename uses '/' as dir seperator.
    std::string filename = QDir::fromNativeSeparators ( Filename ).toLatin1().constData ();

    for ( unsigned int i = 0; i < _imageList.size(); ++i ) {
      if ( ( _dirName + _imageList[i] ).compare ( filename ) == 0 )
        _currentImageIndex = i;
    }
    updateStatusBar();
  }

  void updateStatusBar ( ) {
    _sbarFileName->setText( ( _dirName + ( ( _imageList.size() > 0 ) ? _imageList [ _currentImageIndex ] : "unknown" ) ).c_str() );
    _sbarFileNum->setText( aol::strprintf ( "%d/%d", _currentImageIndex + 1, _imageList.size() ).c_str() );
  }
  
  QRect getSelectedFrame ( ) {
    return QRect ( QPoint ( min ( _selectedFrameCorners[0].x ( ), _selectedFrameCorners[1].x ( ) ),
                            min ( _selectedFrameCorners[0].y ( ), _selectedFrameCorners[1].y ( ) ) ),
                   QPoint ( max ( _selectedFrameCorners[0].x ( ), _selectedFrameCorners[1].x ( ) ),
                            max ( _selectedFrameCorners[0].y ( ), _selectedFrameCorners[1].y ( ) ) ) );
  }
  
  virtual void highlightSelectedFrame ( ) {
    QRect frame = getSelectedFrame ( );
    _sbarFramePosSize->setText ( aol::strprintf ( "[(%d,%d) %dx%d]", frame.left ( ), frame.top ( ), frame.width ( ), frame.height ( ) ).c_str ( ) );
    
    QImage qimage ( imageLabel->getPixmap ( ).size ( ), QImage::Format_ARGB32_Premultiplied );
    qimage.fill ( qRgba(255,255,255,255) );
    QPainter painter ( &qimage );
    painter.setCompositionMode ( QPainter::CompositionMode_Clear );
    painter.fillRect ( frame, qRgba(0,0,0,0) );
    imageLabel->setOverlay ( qimage, 0.5 );
  }
  
  void setSelectedFrameCorner ( const short Index, const QPoint &MousePos, Qt::KeyboardModifiers &Modifiers ) {
    if ( Modifiers & Qt::ControlModifier ) {
      const aol::Vec2<double> dists ( MousePos.x ( ) - _selectedFrameCorners[1-Index].x ( ), MousePos.y ( ) - _selectedFrameCorners[1-Index].y ( ) );
      _selectedFrameCorners[Index].setX ( _selectedFrameCorners[1-Index].x ( ) + aol::signum<double> ( dists[0] ) * dists.getMinAbsValue ( ) );
      _selectedFrameCorners[Index].setY ( _selectedFrameCorners[1-Index].y ( ) + aol::signum<double> ( dists[1] ) * dists.getMinAbsValue ( ) );
    } else if ( Modifiers & Qt::ShiftModifier ) {
      _selectedFrameCorners[Index] = MousePos;
      const aol::Vec2<double> dists ( MousePos.x ( ) - _selectedFrameCorners[1-Index].x ( ), MousePos.y ( ) - _selectedFrameCorners[1-Index].y ( ) );
      if ( aol::Abs<double> ( dists[0] ) < aol::Abs<double> ( dists[1] ) ) _selectedFrameCorners[Index].setX ( _selectedFrameCorners[1-Index].x ( ) );
      else _selectedFrameCorners[Index].setY ( _selectedFrameCorners[1-Index].y ( ) );
    } else _selectedFrameCorners[Index] = MousePos;
  }

  void keyPressEvent ( QKeyEvent * event ) {
    if ( event->key() == Qt::Key_F ) {
      if ( this->windowState() & Qt::WindowMaximized )
        this->showNormal();
      else
        this->showMaximized();
    } else if ( event->key() == Qt::Key_Enter ) {
      if ( this->windowState() & Qt::WindowFullScreen )
        this->showNormal();
      else
        this->showFullScreen();
    }
  }

  void scaleImage ( const int ExponentChange ) {
    if ( actionFitToWindow->isChecked() )
      actionFitToWindow->trigger ( );
    _scaleExponent += ExponentChange;
    const double scaleFactor = std::pow ( 2., _scaleExponent );
    imageLabel->resize ( scaleFactor * _image[0].getNumX ( ), scaleFactor * _image[0].getNumY ( ) );
  }
  
  virtual void mouseMoveEvent ( QMouseEvent *event ) {
    QPoint mousePos = imageLabel->mapToPixmapFromGlobal ( QPoint ( event->x ( ), event->y ( ) ) );
    mousePos.setX ( aol::Clamp ( mousePos.x ( ), 0, imageLabel->getPixmap ( ).width ( ) - 1 ) );
    mousePos.setY ( aol::Clamp ( mousePos.y ( ), 0, imageLabel->getPixmap ( ).height ( ) - 1 ) );
    Qt::KeyboardModifiers modifiers = event->modifiers ( );
    if ( event->buttons ( ) == Qt::LeftButton ) setSelectedFrameCorner ( 0, mousePos, modifiers );
    if ( event->buttons ( ) == Qt::RightButton )  setSelectedFrameCorner ( 1, mousePos, modifiers );
    if ( event->buttons ( ) == Qt::LeftButton || event->buttons ( ) == Qt::RightButton ) highlightSelectedFrame ( );
    _sbarMousePos->setText ( aol::strprintf ( "(%d,%d)", mousePos.x ( ), mousePos.y ( ) ).c_str ( ) );
  }
  
  void setDefaultSaveFormat ( const std::string &Format ) {
    _defaultSaveFormat = Format;
  }

protected slots:
  void on_actionShowInfo_triggered () {
    std::stringstream message;
    message << "Size:          " << _image[0].getNumX ( ) << "x" << _image[0].getNumY ( ) << endl;
    message << "Value range:   " << _image.getMinValue ( ) << " to " << _image.getMaxValue ( ) << endl;
    if ( !imageIsGrayScale ( ) ) {
      message << "Channel value ranges:" << std::endl;
      message << "Red:           " << _image[0].getMinValue ( ) << " to " << _image[0].getMaxValue ( ) << endl;
      message << "Green:         " << _image[1].getMinValue ( ) << " to " << _image[1].getMaxValue ( ) << endl;
      message << "Blue:          " << _image[2].getMinValue ( ) << " to " << _image[2].getMaxValue ( ) << endl;
    }
    QMessageBox::information(NULL, "Info", message.str ().c_str());
  }

  void on_actionEnhanceContrast_triggered () {
    if ( hasImageLoaded() ) {
      aol::Vec2<double> minMax = _image.getSaturatedMinMaxValue ( 0.15 );
      _image.clamp ( minMax[0], minMax[1] );
      _minVal = minMax[0];
      _maxVal = minMax[1];
      refreshPixmap ( );
    }
  }

  void on_actionGrayValuesToHSV_triggered () {
    if ( hasImageLoaded() && imageIsGrayScale ( ) ) {
      const aol::RGBColorMap<double> hsvMap ( 0., 1., aol::HSV_BLUE_TO_RED );
      aol::Vec2<double> minMax = _image[0].getMinMaxValue ( );

      for ( int j = 0; j < _image[0].getNumY ( ); ++j ) {
        for ( int i = 0; i < _image[0].getNumX ( ); ++i ) {
          aol::Vec3<double> color;
          hsvMap.scalarToColor ( ( _image[0].get ( i, j ) - minMax[0] ) / ( minMax[1] - minMax[0] ), color );
          color *= 255;
          for ( unsigned short c=0; c<3 ; ++ c )
            _image[c].set ( i, j, color[c] );
        }
      }
      _minVal = 0;
      _maxVal = 255;

      refreshPixmap ( );
    } else QMessageBox::information(NULL, "Not applicable", "Conversion to HSV colormap is only applicable to grayscale images." );
  }

  void on_actionNextFile_triggered () {
    _currentImageIndex = ( _currentImageIndex + 1 ) % _imageList.size();
    loadImage ( ( _dirName + _imageList[_currentImageIndex] ).c_str() );
    updateStatusBar();
  }

  void on_actionPrevFile_triggered () {
    _currentImageIndex = ( _currentImageIndex + _imageList.size() - 1 ) % _imageList.size();
    loadImage ( ( _dirName + _imageList[_currentImageIndex] ).c_str() );
    updateStatusBar();
  }

  void on_actionLoad_triggered () {
    std::map<std::string, std::string> supportedLoadFormats = getSupportedLoadFormats ( );
    std::map<std::string, std::string>::iterator it = supportedLoadFormats.begin ( );
    std::stringstream ss;
    ss << "Images (*" << it->first; ++it;
    for ( ; it != supportedLoadFormats.end ( ) ; ++it )
      ss << " *" << it->first;
    ss << ");;All files (*.*)";
    const QString path = QPreviousFolderFileDialog::getOpenFileName ( this, "Select a 2D image file",
                                                                      QString ( ),
                                                                      tr ( ss.str ( ).c_str ( ) ) );
    if ( path.isEmpty() == false )
      loadAndShowImage ( path.toLatin1().constData ( ) );
  };
  
  void on_actionReload_triggered () {
    if ( _currentImageIndex < _imageList.size ( ) )
      loadImage ( ( _dirName + _imageList[_currentImageIndex] ).c_str() );
    imageLabel->clearPermanentOverlay ( );
  }
  
  void on_actionSaveAs_triggered () {
    std::map<std::string, std::string> supportedSaveFormats = getSupportedSaveFormats ( );
    std::stringstream ss;
    ss << supportedSaveFormats[_defaultSaveFormat] << " (*" << _defaultSaveFormat << ")";
    
    supportedSaveFormats.erase ( _defaultSaveFormat );
    for ( std::map<std::string, std::string>::iterator it = supportedSaveFormats.begin ( ); it != supportedSaveFormats.end ( ) ; ++it )
      ss << ";;" << it->second << " (*" << it->first << ")";
    const QString path = QPreviousFolderFileDialog::getSaveFileName ( this, "Select destination",
                                                                      QString ( ),
                                                                      tr ( ss.str ( ).c_str ( ) ) );    
    const char* dest = path.toStdString ( ).c_str ( );
    
    if ( path.isEmpty() == false ) {
      try {
        if ( aol::fileNameEndsWith ( dest, ".png" ) ) {
          _image.setOverflowHandling ( aol::SCALE, _image.getMinValue ( ), _image.getMaxValue ( ) );
          if ( imageIsGrayScale ( ) ) _image[0].savePNG ( dest );
          else _image.savePNG ( dest );
        } else {
          qc::ScalarArray<double, qc::QC_2D> gray;
          rgbToGray ( _image, gray );
          
          if ( aol::fileNameEndsWith ( dest, ".q2bz" ) ) {
            gray.save ( dest, qc::PGM_DOUBLE_BINARY );
          } else if ( aol::fileNameEndsWith ( dest, ".pgm" ) ) {
            gray.setOverflowHandlingToCurrentValueRange ( );
            gray.save ( dest, qc::PGM_UNSIGNED_CHAR_BINARY );
          } else if ( aol::fileNameEndsWith ( dest, ".tif" ) ) {
            gray.saveTIFF ( dest );
          } else QMessageBox::information(NULL, "Output not supported", "Saving is not supported for the selected image format yet." );
        }
      } catch ( aol::Exception &el ) {
        el.dump ( );
        imageLabel->setText( el.getMessage ( ).c_str ( ) );
        _image.reallocate ( 1, 1 );
      }
    }
  }

  void on_actionPlotHistogram_triggered () {
    if ( aol::runGnuplot ( "--version" ) == false )
      QMessageBox::information(NULL, "Error", "Calling gnuplot failed! Make sure gnuplot is in your search path.");
    else {
      aol::Vector<int> histo;
      qc::ScalarArray<double, qc::QC_2D> grayImage;
      rgbToGray ( _image, grayImage );
      grayImage.createHistogramOfValues ( histo, 256 );
      aol::plotHistogram<double> ( histo, "", true, true );
    }
  };

  void on_actionPlotCrossSection_triggered () {
    if ( aol::runGnuplot ( "--version" ) == false )
      QMessageBox::information(NULL, "Error", "Calling gnuplot failed! Make sure gnuplot is in your search path.");
    else if ( imageIsGrayScale() ) {
      aol::Plotter<double> plotter;
      aol::PlotDataFileHandler<double> plotHandler;
      qc::ScalarArray<double, qc::QC_2D> grayImage;
      rgbToGray ( _image, grayImage );
      aol::Vec2<short> p1 ( _selectedFrameCorners[0].x ( ), _selectedFrameCorners[0].y ( ) );
      aol::Vec2<short> p2 ( _selectedFrameCorners[1].x ( ), _selectedFrameCorners[1].y ( ) );
      plotHandler.generateCrossSectionPlot ( grayImage, p1, p2 );
      plotter.addPlotCommandsFromHandler ( plotHandler );
      plotter.plotToScreen();
    }
  };

  void on_actionFitToWindow_triggered () {
    if ( actionFitToWindow->isChecked() )
      scrollArea->setWidgetResizable ( true );
    else {
      scrollArea->setWidgetResizable ( false );
      imageLabel->resize ( _image[0].getNumX ( ), _image[0].getNumY ( ) );
    }
  };

  void on_actionAutoContrast_triggered () {
    if ( actionAutoContrast->isChecked() )
      on_actionEnhanceContrast_triggered();
  };

  void on_actionZoomIn_triggered () {
    scaleImage ( 1 );
  };

  void on_actionZoomOut_triggered () {
    scaleImage ( -1 );
  };
  
  void on_actionCrop_triggered () {
    QRect frame = getSelectedFrame ( );
    for ( int c=0; c<3 ; ++c )
      _image[c].crop ( aol::Vec2<int> ( frame.left ( ), frame.top ( ) ), aol::Vec2<int> ( frame.width ( ), frame.height ( ) ) );
    imageLabel->cropPixmap ( frame );
    refreshPixmap ( );
  }
  
  void on_actionRGBToGrayLuminance_triggered () {
    qc::ScalarArray<double, qc::QC_2D> gray;
    rgbToGray ( _image, gray, LUMINANCE );
    for ( unsigned short c=0; c<3 ; ++c )
      _image[c] = gray;
    refreshPixmap ( );
  }
  
  void on_actionRGBToGrayAverage_triggered () {
    qc::ScalarArray<double, qc::QC_2D> gray;
    rgbToGray ( _image, gray, AVERAGE );
    for ( unsigned short c=0; c<3 ; ++c )
      _image[c] = gray;
    refreshPixmap ( );
  }
  
  void on_actionRGBToGrayLightness_triggered () {
    qc::ScalarArray<double, qc::QC_2D> gray;
    rgbToGray ( _image, gray, LIGHTNESS );
    for ( unsigned short c=0; c<3 ; ++c )
      _image[c] = gray;
    refreshPixmap ( );
  }
  
  void on_actionFourierPowerSpectrum_triggered () {
    if ( hasImageLoaded() && imageIsGrayScale ( ) ) {
      qc::ScalarArray<double, qc::QC_2D> modulus ( _image[0].getNumX ( ), _image[0].getNumY ( ) );
      qc::computeLogFFTModulus<double> ( _image[0], modulus, 0, false );
      _image[0] = modulus;
      _image[1] = _image[0];
      _image[2] = _image[0];
      refreshPixmap ( );
    }
  }
  
  void on_actionLogScale_triggered () {
    double min = _image.getMinValue ( );
    for ( int i=0; i<_image[0].size ( ) ; ++i )
      for ( unsigned short c=0; c<3 ; ++c )
        _image[c][i] = log ( min + _image[c][i] + 1 );
    on_actionEnhanceContrast_triggered ( );
  }
  
signals:
  void imageLoaded ( );
};

class QuocViewerApp : public QApplication {
  Q_OBJECT

  std::vector<QuocImageViewApp *> _dialogs;
public:
  QuocViewerApp ( int & argc, char **argv ) : QApplication ( argc, argv ) {
    QuocImageViewApp *dialog = new QuocImageViewApp ( ( argc > 1 ) ? argv[1] : NULL );
    dialog->show();
    _dialogs.push_back ( dialog );
  }
  virtual ~QuocViewerApp() {
    for ( unsigned int i = 0; i < _dialogs.size(); ++i )
      delete _dialogs[i];
  }

  bool notify ( QObject *Receiver, QEvent *Event ) {
    try {
      return QApplication::notify( Receiver, Event );
    }
    catch ( aol::Exception &el ) {
      QMessageBox::information(NULL, "Info", el.getMessage().c_str());
      return false;
    }
  }
protected:
  // This is necessary because OS X doesn't tell an app which file to open with command line arguments, but with an event.
  bool event ( QEvent *Event ) {
    switch ( Event->type() ) {
      case QEvent::FileOpen: {
        QString filename = static_cast<QFileOpenEvent *>(Event)->file();
        if ( ( _dialogs.size() == 1 ) && ( _dialogs[0]->hasImageLoaded() == false ) ) {
          _dialogs[0]->loadAndShowImage ( filename.toLatin1().constData () );
        }
        else {
          QuocImageViewApp *dialog = new QuocImageViewApp( filename.toLatin1().constData () );
          dialog->show();
          _dialogs.push_back ( dialog );
        }
        return true;
      }
      default:
        return QApplication::event ( Event );
    }
  }
};

#endif // __QCIMVIEWAPP_H
