#ifndef __ELMAADDON_H
#define __ELMAADDON_H

#include <QDialog>

#include <elmaReadme.h>
#include <qcImViewApp.hpp>
#include "qtWorkers.hpp"
#include "qtInputDialogs.hpp"
#include "previousFolderFileDialog.hpp"
#include <atomFinder.h>
#include <emImgQualityQuantifier.h>
#include <displayLicenses.hpp>
#include <QProgressDialog>

/**
 * \brief ELectron Micrograph Analysis tool
 *
 * \author Mevenkamp, Berkels
 */
class ELMAImageViewApp : public QuocImageViewApp {
  Q_OBJECT
protected:
  qc::ScalarArray<double, qc::QC_2D> _electronMicrograph;
  aol::MultiVector<double> _atomPositions;
  aol::MultiVector<double> _2DGaussianParameters;
  qc::BitArray<qc::QC_2D> _segmented;
  qc::ScalarArray<double, qc::QC_2D> _bumpFunctionImg;
public:
  QMenu *menuEMTools;
  QAction *actionDenoise;
  
  QMenu *menuEMAnalysis;
  QAction *actionDetectAtoms;
  QActionGroup *viewAtomDetectionResultsGroup;
  QAction *actionViewAtomPositions;
  QActionGroup *viewAtomDetectionImagesGroup;
  QAction *actionViewElectronMicrograph;
  QAction *actionViewSegmented;
  QAction *actionView2DGaussians;
  QAction *actionSaveAtomPositions;
  QAction *actionSave2DGaussianParameters;
  QAction *actionViewQuantitativeMeasures;
  
  QMenu *menuAboutExternalLibraries;
  QAction *actionAboutELMA;
  
  ELMAImageViewApp ( const char * InputFilename = NULL ) : QuocImageViewApp ( InputFilename ) {
    this->setWindowTitle(QApplication::translate("QuocImageViewer", "ELMA - ELectron Micrograph Analysis", 0));
    setupMenus ( );
    this->setDefaultSaveFormat ( ".tif" );
  }
protected:
  bool runWorker ( ProgressWorker *Worker ) {
    QThread *thread = new QThread;
    QProgressDialog *progressDialog = new QProgressDialog ( Worker->getText ( ), "Abort", 0, 0, this );
    QProgressDialog *abortingDialog = new QProgressDialog ( "Aborting current computation...", "", 0, 0, this );
    abortingDialog->setCancelButton ( 0 );
    progressDialog->setWindowModality ( Qt::WindowModal );
    progressDialog->setAutoClose ( false );
    progressDialog->show ( );
    connect ( progressDialog, SIGNAL ( canceled ( ) ), this, SLOT ( stop ( ) ) );
    connect ( progressDialog, SIGNAL ( canceled ( ) ), abortingDialog, SLOT ( show ( ) ) );
    connect ( Worker, SIGNAL ( finished ( ) ), thread, SLOT ( quit ( ) ) );
    connect ( Worker, SIGNAL ( progressValueChanged ( int ) ), progressDialog, SLOT ( setValue ( int ) ) );
    connect ( Worker, SIGNAL ( maximumValueChanged ( int ) ), progressDialog, SLOT ( setMaximum ( int ) ) );
    connect ( Worker, SIGNAL ( labelTextChanged ( const QString & ) ), progressDialog, SLOT ( setLabelText ( const QString & ) ) );
    Worker->moveToThread ( thread );
    thread->start ( );
    QMetaObject::invokeMethod ( Worker, "doWork", Qt::QueuedConnection );
    while ( thread->isRunning ( ) ) { QCoreApplication::processEvents ( ); }
    progressDialog->close ( );
    abortingDialog->close ( );
    aol::resetCtrlCState ( );
    
    return Worker->hasFinished ( );
  }
  
protected slots:
  void stop ( ) {
    aol::ctrlCHandler ( 1 );
  }
  
  void clearAnalysis ( ) {
    this->imageLabel->clearPermanentOverlay ( );
    _atomPositions.clear ( );
    _2DGaussianParameters.clear ( );
    _segmented.reallocate ( 0, 0 );
    _bumpFunctionImg.reallocate ( 0, 0 );
    setAtomDetectionResultsEnabled ( false );
  }
  
  void loadElectronMicrographAndClearAnalysis ( ) {
    if ( _electronMicrograph.size ( ) > 0 ) {
      for ( int c=0; c<3 ; ++c )
        _image[c] = _electronMicrograph;
      this->refreshPixmap ( );
    }
    clearAnalysis ( );
  }
  
  void setAtomDetectionResultsEnabled ( bool Enabled ) {
    viewAtomDetectionResultsGroup->setEnabled ( Enabled );
    viewAtomDetectionImagesGroup->setEnabled ( Enabled );
    if ( !Enabled ) {
      uncheckQActionGroup ( viewAtomDetectionResultsGroup );
      uncheckQActionGroup ( viewAtomDetectionImagesGroup );
    }
  }
  
  void enableImageDependentActions ( ) {
    this->actionReload->setEnabled ( true );
    this->actionSaveAs->setEnabled ( true );
    this->menuEMTools->menuAction ( )->setEnabled ( true );
    this->menuEMAnalysis->menuAction ( )->setEnabled ( true );
  }
  
  void on_actionDenoise_triggered () {
    if ( !this->imageIsGrayScale ( ) || this->_image.getMinValue ( ) < 0 ) QMessageBox::information(NULL, "Not applicable", "Denoising is only applicable to images given in #detected electrons per pixel." );
    else {
      QDenoisingDialog dialog ( this );
      dialog.exec ( );
      if ( dialog.result ( ) == QDialog::Accepted ) {
        loadElectronMicrographAndClearAnalysis ( );
        
        DenoisingWorker *worker;
        if ( dialog.filterType ( ) == 0 ) {
          EMBM3DOptions<double, qc::ScalarArray<double, qc::QC_2D> > options ( this->_image[0] );
          dialog.setEMBM3DOptions ( options );
          worker = new BM3DWorker ( this->imageLabel, this->_image, options );
        } else {
          EMNLMOptions<double, qc::ScalarArray<double, qc::QC_2D> > options ( this->_image[0] );
          dialog.setEMNLMOptions ( options );
          worker = new NLMWorker ( this->imageLabel, this->_image, options );
        }
        
        if ( runWorker ( worker ) ) {
          for ( unsigned int c=0; c<3 ; ++c )
            this->_image[c] = worker->getEstimate ( );
          this->refreshPixmap ( );
        }
      }
    }
  }
  
  void on_actionDetectAtoms_triggered () {
    if ( !this->imageIsGrayScale ( ) ) QMessageBox::information(NULL, "Not applicable", "Atom detection is only applicable to gray scale electron micrographs." );
    else {
      QDetectAtomsDialog dialog ( this );
      dialog.exec ( );
      if ( dialog.result ( ) == QDialog::Accepted ) {
        loadElectronMicrographAndClearAnalysis ( );
        
        AtomDetectionWorker *worker = new AtomDetectionWorker ( this->imageLabel, this->_image, this->_image[0], dialog.atomType ( ),
                                                                dialog.gamma ( ), MAXIT, dialog.epsilon ( ) );
        if ( runWorker ( worker ) ) {
          _electronMicrograph.reallocate ( this->_image[0].getNumX ( ), this->_image[0].getNumY ( ) );
          _electronMicrograph = worker->getElectronMicrograph ( );
          actionViewElectronMicrograph->setChecked ( true );
          on_actionViewElectronMicrograph_triggered ();
          _atomPositions.reallocate ( worker->getAtomPositions ( ) );
          _atomPositions = worker->getAtomPositions ( );
          _2DGaussianParameters.reallocate ( worker->get2DGaussianParameters ( ) );
          _2DGaussianParameters = worker->get2DGaussianParameters ( );
          _segmented.reallocate ( this->_image[0].getNumX ( ), this->_image[0].getNumY ( ) );
          _segmented = worker->getSegmented ( );
          actionViewAtomPositions->setChecked ( true );
          on_actionViewAtomPositions_triggered ();
          _bumpFunctionImg.reallocate ( _electronMicrograph.getNumX ( ), _electronMicrograph.getNumY ( ) );
          AtomFinder<double, aol::FullMatrix<double>, LinearRegressionQR<double>,
                     qc::ScalarArray<double, qc::QC_2D>, qc::MultiArray<double, qc::QC_2D, 3> >::getBumpFunctionImage ( _2DGaussianParameters, _bumpFunctionImg );
          _bumpFunctionImg.scaleValuesTo01 ( );
          _bumpFunctionImg *= 255;
          setAtomDetectionResultsEnabled ( true );
        }
      }
    }
  }
  
  void on_actionViewAtomPositions_triggered () {
    if ( _atomPositions.numComponents ( ) > 0 && actionViewAtomPositions->isChecked ( ) ) {
      QImage qimage ( imageLabel->getPixmap ( ).size ( ), QImage::Format_ARGB32_Premultiplied );
      qimage.fill ( qRgba(0,0,0,0) );
      for ( int i=0; i<_atomPositions.numComponents ( ) ; ++i )
        qimage.setPixel ( round ( _atomPositions[i][0] ), round ( _atomPositions[i][1] ), qRgba(255,0,0,255) );
      this->imageLabel->setPermanentOverlay ( qimage, 1.0 );
    } else this->imageLabel->clearPermanentOverlay ( );
  }
  
  void on_actionSaveAtomPositions_triggered () {
    std::stringstream ss;
    ss << "Comma Separated Values (*.csv);;aol::MultiVector (*.dat)";

    const QString path = QPreviousFolderFileDialog::getSaveFileName ( this, "Select destination",
                                                                      QString ( ),
                                                                      tr ( ss.str ( ).c_str ( ) ) );
    const char* dest = path.toStdString ( ).c_str ( );
    
    if ( path.isEmpty() == false ) {
      try {
        if ( aol::fileNameEndsWith ( dest, ".csv" ) ) {
          std::ofstream txtFile ( dest );
          txtFile << "x, y" << std::endl;
          for ( int i=0; i<_atomPositions.numComponents ( ) ; ++i ) {
            for ( short j=0; j<_atomPositions[i].size ( ) ; ++j ) {
              if ( j > 0 ) txtFile << ", ";
              txtFile << _atomPositions[i][j];
            }
            txtFile << std::endl;
          }
          txtFile.close ( );
        } else if ( aol::fileNameEndsWith ( dest, ".dat" ) ) {
          _atomPositions.save ( dest );
        } else QMessageBox::information(NULL, "Output not supported", "Saving is not supported for the selected MultiVector format." );
      } catch ( aol::Exception &el ) {
        el.dump ( );
        imageLabel->setText( el.getMessage ( ).c_str ( ) );
        _image.reallocate ( 1, 1 );
      }
    }
  }
  
  void on_actionSave2DGaussianParameters_triggered () {
    if ( _2DGaussianParameters.numComponents ( ) > 0 ) {
      std::stringstream ss;
      ss << "Comma Separated Values (*.csv);;aol::MultiVector (*.dat)";
      const QString path = QPreviousFolderFileDialog::getSaveFileName ( this, "Select destination",
                                                                        QString ( ),
                                                                        tr ( ss.str ( ).c_str ( ) ) );
      const char* dest = path.toStdString ( ).c_str ( );
      
      if ( path.isEmpty() == false ) {
        try {
          if ( aol::fileNameEndsWith ( dest, ".csv" ) ) {
            int k = 0;
            // Save parameters of fitted single atoms (if any)
            if ( _2DGaussianParameters[0].size ( ) == AsymmetricGaussianBumpFunction<double>::NumberOfParameters ) {
              std::ofstream txtFile ( dest );
              txtFile << "x, y, amplitude, width, height, rotation, background" << std::endl;
              while ( k < _2DGaussianParameters.numComponents ( ) && _2DGaussianParameters[k].size ( ) == AsymmetricGaussianBumpFunction<double>::NumberOfParameters ) {
                for ( short j=0; j<_2DGaussianParameters[k].size ( ) ; ++j ) {
                  if ( j > 0 ) txtFile << ", ";
                  txtFile << _2DGaussianParameters[k][j];
                }
                txtFile << std::endl;
                ++k;
              }
              txtFile.close ( );
            }
            // Save parameters of fitted dumbbells (if any)
            if ( k < _2DGaussianParameters.numComponents ( ) ) {
              const char* dumbbellDest;
              if ( k > 0 ) {
                QFileInfo fi ( path );
                dumbbellDest = ( fi.path ( ) + "/" + fi.baseName ( ) + "_dumbbell.csv" ).toStdString ( ).c_str ( );
              } else dumbbellDest = dest;
              std::ofstream txtFile ( dumbbellDest );
              txtFile << "x1, y1, x2, y2, amplitude1, amplitude2, width1, height1, width2, height2, rotation1, rotation2, background" << std::endl;
              while ( k < _2DGaussianParameters.numComponents ( ) ) {
                for ( short j=0; j<_2DGaussianParameters[k].size ( ) ; ++j ) {
                  if ( j > 0 ) txtFile << ", ";
                  txtFile << _2DGaussianParameters[k][j];
                }
                txtFile << std::endl;
                ++k;
              }
              txtFile.close ( );
            }
          } else if ( aol::fileNameEndsWith ( dest, ".dat" ) ) {
            _2DGaussianParameters.save ( dest );
          } else QMessageBox::information(NULL, "Output not supported", "Saving is not supported for the selected MultiVector format." );
        } catch ( aol::Exception &el ) {
          el.dump ( );
          imageLabel->setText( el.getMessage ( ).c_str ( ) );
          _image.reallocate ( 1, 1 );
        }
      }
    } else QMessageBox::information(NULL, "No parameters available", "There are no parameters available from 2D Gaussian fitting." );
  }
  
  void on_actionViewElectronMicrograph_triggered () {
    if ( _electronMicrograph.size ( ) > 0 ) {
      for ( int c=0; c<3 ; ++c )
        _image[c] = _electronMicrograph;
      this->refreshPixmap ( );
    }
  }
  
  void on_actionViewSegmented_triggered () {
    if ( _segmented.getNumX ( ) == _electronMicrograph.getNumX (  ) && _segmented.getNumY ( ) == _electronMicrograph.getNumY ( ) ) {
      for ( int c=0; c<3 ; ++c )
        for ( int k=0; k<_image[c].size ( ) ; ++k )
          _image[c][k] = 255 * _segmented[k];
      this->refreshPixmap ( );
    }
  }
  
  void on_actionView2DGaussians_triggered () {
    if ( _bumpFunctionImg.getNumX ( ) == _electronMicrograph.getNumX (  ) && _bumpFunctionImg.getNumY ( ) == _electronMicrograph.getNumY ( ) ) {
      for ( int c=0; c<3 ; ++c )
        _image[c] = _bumpFunctionImg;
      this->refreshPixmap ( );
    }
  }
  
  void on_actionViewQuantitativeMeasures_triggered () {
    QQuantitativeMeasuresDialog dialog ( this );
    dialog.exec ( );
    if ( dialog.result ( ) == QDialog::Accepted ) {
      EMImgQualityQuantifier<double, aol::FullMatrix<double>, LinearRegressionQR<double>,
      qc::ScalarArray<double, qc::QC_2D>, qc::MultiArray<double, qc::QC_2D, 3> > emImgQualityQuantifier;
      
      double periodX, periodY, angleX, angleY;
      dialog.setPrecisionParameters ( periodX, periodY, angleX, angleY );
      
      std::string statistics = emImgQualityQuantifier.getStatistics ( _atomPositions, periodX, periodY, 5, angleX, angleY, 5 );
      QMessageBox *msg = new QMessageBox(QMessageBox::Information, "Quantitative measures", statistics.c_str ( )) ;
      QFont font("Monospace");
      font.setStyleHint(QFont::TypeWriter);
      msg->setFont(font);
      msg->exec();
    }
  }
  
  void on_actionAboutELMA_triggered () {
    QMessageBox msgBox(this);
    msgBox.setWindowTitle("ELMA (ELectron Micrograph Analysis)");
    msgBox.setTextFormat(Qt::RichText);
    msgBox.setText ( QString ( getELMAReadme ( ).c_str ( ) ) );
    QSpacerItem* horizontalSpacer = new QSpacerItem(600, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
    QGridLayout* layout = static_cast<QGridLayout*>(msgBox.layout());
    layout->addItem(horizontalSpacer, layout->rowCount(), 0, 1, layout->columnCount());
    msgBox.exec();
  }
protected:
#ifdef ELMA_DEPLOY
  void refreshPixmap ( ) {
    aol::Vec2<double> minMax = _image.getSaturatedMinMaxValue ( 0.15 );
    _image.clamp ( minMax[0], minMax[1] );
    _minVal = minMax[0];
    _maxVal = minMax[1];
    QuocImageViewApp::refreshPixmap ( );
  }
  
  void mouseMoveEvent ( QMouseEvent *event ) {
    QPoint mousePos = imageLabel->mapToPixmapFromGlobal ( QPoint ( event->x ( ), event->y ( ) ) );
    mousePos.setX ( aol::Clamp ( mousePos.x ( ), 0, imageLabel->getPixmap ( ).width ( ) - 1 ) );
    mousePos.setY ( aol::Clamp ( mousePos.y ( ), 0, imageLabel->getPixmap ( ).height ( ) - 1 ) );
    _sbarMousePos->setText ( aol::strprintf ( "(%d,%d)", mousePos.x ( ), mousePos.y ( ) ).c_str ( ) );
  }
  
  void highlightSelectedFrame ( ) {
  }
#endif
  
  void setupMenus ( ) {
#ifdef ELMA_DEPLOY
    const char* emToolsTxt = "Tools";
    const char* emAnalysisTxt = "Analysis";
    
    // Remove unnecessary QuocViewer menus
    this->menubar->removeAction ( this->menuTools->menuAction ( ) );
    this->menubar->removeAction ( this->menuView->menuAction ( ) );
    
    this->actionAutoContrast->setChecked ( true );
#else
    const char* emToolsTxt = "EM Tools";
    const char* emAnalysisTxt = "EM Analysis";
#endif
    QObject::connect(this, SIGNAL(imageLoaded()), this, SLOT(enableImageDependentActions()));
    this->actionReload->setEnabled ( false );
    this->actionSaveAs->setEnabled ( false );
    
    /*
     * BEGIN: EM Tools
     */
    menuEMTools = new QMenu(this->menubar);
    menuEMTools->setObjectName(QString::fromUtf8("menuEMTools"));
    menuEMTools->setTitle(QApplication::translate("QuocImageViewer", emToolsTxt, 0));
    this->menubar->addMenu ( menuEMTools );
    menuEMTools->menuAction ( )->setEnabled ( false );
    
    // Denoise
    actionDenoise = new QAction(this);
    actionDenoise->setObjectName(QString::fromUtf8("actionDenoise"));
    actionDenoise->setText(QApplication::translate("QuocImageViewer", "Denoise...", 0));
    menuEMTools->addAction(actionDenoise);
    QObject::connect(actionDenoise, SIGNAL(triggered()), this, SLOT(on_actionDenoise_triggered()));
    /*
     * END: EM Tools
     */
    
    
    /*
     * BEGIN: EM Analysis
     */
    menuEMAnalysis = new QMenu(this->menubar);
    menuEMAnalysis->setObjectName(QString::fromUtf8("menuEMAnalysis"));
    menuEMAnalysis->setTitle(QApplication::translate("QuocImageViewer", emAnalysisTxt, 0));
    this->menubar->addMenu ( menuEMAnalysis );
    menuEMAnalysis->menuAction ( )->setEnabled ( false );
    
    // Detect atoms
    actionDetectAtoms = new QAction(this);
    actionDetectAtoms->setObjectName(QString::fromUtf8("actionDetectAtoms"));
    actionDetectAtoms->setText(QApplication::translate("QuocImageViewer", "Detect atoms...", 0));
    menuEMAnalysis->addAction(actionDetectAtoms);
    QObject::connect(actionDetectAtoms, SIGNAL(triggered()), this, SLOT(on_actionDetectAtoms_triggered()));
    
    /*
     * BEGIN: Atom detection results
     */
    menuEMAnalysis->addSeparator();
    viewAtomDetectionResultsGroup = new QActionGroup(this);
    viewAtomDetectionResultsGroup->setExclusive(false);
    connect(actionReload, SIGNAL(triggered()), this, SLOT(clearAnalysis()));
    connect(actionCrop, SIGNAL(triggered()), this, SLOT(clearAnalysis()));
    
    // View atom positions
    actionViewAtomPositions = new QAction(viewAtomDetectionResultsGroup);
    actionViewAtomPositions->setObjectName(QString::fromUtf8("actionViewAtomPositions"));
    actionViewAtomPositions->setText(QApplication::translate("QuocImageViewer", "View atom centers", 0));
    actionViewAtomPositions->setCheckable(true);
    menuEMAnalysis->addAction(actionViewAtomPositions);
    QObject::connect(actionViewAtomPositions, SIGNAL(triggered()), this, SLOT(on_actionViewAtomPositions_triggered()));
    
    // Save atom positions
    actionSaveAtomPositions = new QAction(viewAtomDetectionResultsGroup);
    actionSaveAtomPositions->setObjectName(QString::fromUtf8("actionSaveAtomPositions"));
    actionSaveAtomPositions->setText(QApplication::translate("QuocImageViewer", "Save atom centers...", 0));
    menuEMAnalysis->addAction(actionSaveAtomPositions);
    QObject::connect(actionSaveAtomPositions, SIGNAL(triggered()), this, SLOT(on_actionSaveAtomPositions_triggered()));
    
    // Save 2D Gaussian parameters
    actionSave2DGaussianParameters = new QAction(viewAtomDetectionResultsGroup);
    actionSave2DGaussianParameters->setObjectName(QString::fromUtf8("actionSave2DGaussianParameters"));
    actionSave2DGaussianParameters->setText(QApplication::translate("QuocImageViewer", "Save Gaussian fit parameters...", 0));
    menuEMAnalysis->addAction(actionSave2DGaussianParameters);
    QObject::connect(actionSave2DGaussianParameters, SIGNAL(triggered()), this, SLOT(on_actionSave2DGaussianParameters_triggered()));
    
    /*
     * BEGIN: View atom detection images
     */
    menuEMAnalysis->addSeparator();
    viewAtomDetectionImagesGroup = new QActionGroup(this);
    viewAtomDetectionImagesGroup->setExclusive(true);
    
    // View Electron micrograph
    actionViewElectronMicrograph = new QAction(viewAtomDetectionImagesGroup);
    actionViewElectronMicrograph->setObjectName(QString::fromUtf8("actionViewElectronMicrograph"));
    actionViewElectronMicrograph->setText(QApplication::translate("QuocImageViewer", "View electron micrograph", 0));
    actionViewElectronMicrograph->setCheckable(true);
    menuEMAnalysis->addAction(actionViewElectronMicrograph);
    QObject::connect(actionViewElectronMicrograph, SIGNAL(triggered()), this, SLOT(on_actionViewElectronMicrograph_triggered()));
    
    // View Segmented
    actionViewSegmented = new QAction(viewAtomDetectionImagesGroup);
    actionViewSegmented->setObjectName(QString::fromUtf8("actionViewSegmented"));
    actionViewSegmented->setText(QApplication::translate("QuocImageViewer", "View segmented image", 0));
    actionViewSegmented->setCheckable(true);
    menuEMAnalysis->addAction(actionViewSegmented);
    QObject::connect(actionViewSegmented, SIGNAL(triggered()), this, SLOT(on_actionViewSegmented_triggered()));
    
    // View 2D Gaussians
    actionView2DGaussians = new QAction(viewAtomDetectionImagesGroup);
    actionView2DGaussians->setObjectName(QString::fromUtf8("actionView2DGaussians"));
    actionView2DGaussians->setText(QApplication::translate("QuocImageViewer", "View 2D Gaussians", 0));
    actionView2DGaussians->setCheckable(true);
    menuEMAnalysis->addAction(actionView2DGaussians);
    QObject::connect(actionView2DGaussians, SIGNAL(triggered()), this, SLOT(on_actionView2DGaussians_triggered()));
    
    menuEMAnalysis->addSeparator();
    /*
     * END: View atom detection images
     */
    
    // View quantitative measures
    actionViewQuantitativeMeasures = new QAction(viewAtomDetectionResultsGroup);
    actionViewQuantitativeMeasures->setObjectName(QString::fromUtf8("actionViewQuantitativeMeasures"));
    actionViewQuantitativeMeasures->setText(QApplication::translate("QuocImageViewer", "View quantitative measures", 0));
    menuEMAnalysis->addAction(actionViewQuantitativeMeasures);
    QObject::connect(actionViewQuantitativeMeasures, SIGNAL(triggered()), this, SLOT(on_actionViewQuantitativeMeasures_triggered()));
    
    setAtomDetectionResultsEnabled ( false );
    /*
     * END: Atom detection results
     */
    /*
     * END: EM Analysis
     */
    
    /*
     * BEGIN: About
     */
    // About ELMA
    actionAboutELMA = new QAction(this);
    actionAboutELMA->setObjectName(QString::fromUtf8("actionAboutELMA"));
    actionAboutELMA->setText(QApplication::translate("QuocImageViewer", "About ELMA", 0));
    menuAbout->addAction(actionAboutELMA);
    QObject::connect(actionAboutELMA, SIGNAL(triggered()), this, SLOT(on_actionAboutELMA_triggered()));
    
    // About External Libraries
    menuAboutExternalLibraries = new QAboutExternalLibrariesMenu(this->menuAbout);
    menuAboutExternalLibraries->setObjectName(QString::fromUtf8("menuAboutExternalLibraries"));
    menuAboutExternalLibraries->setTitle(QApplication::translate("QuocImageViewer", "About External Libraries", 0));
    menuAbout->addMenu ( menuAboutExternalLibraries );
    
    // Remove about menu and add again, to make it appear in the right-most menu slot
    this->menubar->removeAction ( this->menuAbout->menuAction ( ) );
    this->menubar->addMenu ( this->menuAbout );
    /*
     * END: About
     */
    
    if ( this->hasImageLoaded ( ) ) enableImageDependentActions ( );
  }
  
  void uncheckQActionGroup ( QActionGroup *actionGroup ) {
    QList<QAction *> list = actionGroup->actions ( );
    QList<QAction *>::iterator i;
    for ( i = list.begin(); i != list.end(); ++i)
      (*i)->setChecked ( false );
  }
};


class ELMAViewerApp : public QApplication {
  Q_OBJECT
  
  std::vector<ELMAImageViewApp *> _dialogs;
public:
  ELMAViewerApp ( int & argc, char **argv ) : QApplication ( argc, argv ) {
    ELMAImageViewApp *dialog = new ELMAImageViewApp ( ( argc > 1 ) ? argv[1] : NULL );
    dialog->show();
    _dialogs.push_back ( dialog );
  }
  virtual ~ELMAViewerApp() {
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
          ELMAImageViewApp *dialog = new ELMAImageViewApp( filename.toLatin1().constData () );
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

#endif // __ELMAADDON_H
