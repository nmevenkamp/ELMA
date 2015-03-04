#ifndef __QTINPUTDIALOGS_H
#define __QTINPUTDIALOGS_H

#ifdef __GNUC__
#pragma GCC system_header
#endif

#include <aol.h>
#include <emBM3DFilter.h>
#include <emNonLocalMeansFilter.h>
#include <atomFinder.h>
#include "ui_denoisingDialog.h"
#include "ui_denoisingDialogHelp.h"
#include "ui_detectAtomsDialog.h"
#include "ui_detectAtomsDialogHelp.h"
#include "ui_quantitativeMeasuresDialog.h"
#include "blocksizespinbox.hpp"


class QDenoisingDialogHelp : public QDialog, protected Ui::DenoisingDialogHelp {
  Q_OBJECT
public:
  QDenoisingDialogHelp ( QWidget * parent = 0, Qt::WindowFlags f = 0 ) : QDialog ( parent, f ) {
    setupUi ( this );
  }
};

class QDenoisingDialog : public QDialog, protected Ui::DenoisingDialog {
  Q_OBJECT
protected:
  QDenoisingDialogHelp *denoisingDialogHelp;
public:
  QDenoisingDialog ( QWidget * parent = 0, Qt::WindowFlags f = 0 ) : QDialog ( parent, f ) {
    setupUi ( this );
    selectBM3DFilterOptions ( );
    
    // Enable help
    helpButton->setIcon ( QIcon ( "://help.png" ) );
    denoisingDialogHelp = new QDenoisingDialogHelp ( this, f );
    connect ( helpButton, SIGNAL ( clicked ( ) ), denoisingDialogHelp, SLOT ( show ( ) ) );
    
    // Enable swapping between EM-BM3D filter options and EM-NLM filter options
    connect ( bm3dFilterRadioButton, SIGNAL ( clicked ( ) ), this, SLOT ( selectBM3DFilterOptions ( ) ) );
    connect ( nlmFilterRadioButton, SIGNAL ( clicked ( ) ), this, SLOT ( selectNLMFilterOptions ( ) ) );
    
    for ( int i=0; i<POISSON_NOISE_ADAPTATION::NUM ; ++i )
      poissonNoiseAdaptationComboBox->addItem ( POISSON_NOISE_ADAPTATION::toString ( i ).c_str ( ) );
    populateBM3DSimilaritySearchMethodComboBox ( );
    for ( int i=0; i<EM_SCANDISTORTIONCORRECTION_METHOD::NUM ; ++i )
      scanDistortionCorrectionMethodComboBox->addItem ( aol::strprintf ( "%s: %s", EM_SCANDISTORTIONCORRECTION_METHOD::getIdentifier ( i ).c_str ( ),
                                                                                   EM_SCANDISTORTIONCORRECTION_METHOD::toString ( i ).c_str ( ) ).c_str ( ) );
    for ( int i=0; i<EMBM3D_PROFILE::NUM ; ++i )
      profileComboBox->addItem ( EMBM3D_PROFILE::toString ( i ).c_str ( ) );
  }
  
  int filterType ( ) const {
    if ( bm3dFilterRadioButton->isChecked ( ) ) return 0;
    else return 1;
  }
  
  void setEMBM3DOptions ( EMBM3DOptions<double, qc::ScalarArray<double, qc::QC_2D> > &Options ) {
    Options.blockSize = blockSizeSpinBox->value ( );
    Options.poissonNoiseAdaptation = poissonNoiseAdaptationComboBox->currentIndex ( );
    Options.similaritySearchMethod = similaritySearchMethodComboBox->currentIndex ( );
    Options.scanDistortionCorrectionMethod = scanDistortionCorrectionMethodComboBox->currentIndex ( );
    Options.profile = profileComboBox->currentIndex ( );
  }
  
  void setEMNLMOptions ( EMNLMOptions<double, qc::ScalarArray<double, qc::QC_2D> > &Options ) {
    Options.blockSize = blockSizeSpinBox->value ( );
    Options.poissonNoiseAdaptation = poissonNoiseAdaptationComboBox->currentIndex ( );
    Options.similaritySearchMethod = similaritySearchMethodComboBox->currentIndex ( );
    Options.scanDistortionCorrectionMethod = scanDistortionCorrectionMethodComboBox->currentIndex ( );
    Options.filterParameter = filterParameterDoubleSpinBox->value ( );
  }
protected:
  void populateBM3DSimilaritySearchMethodComboBox ( ) {
    similaritySearchMethodComboBox->clear ( );
    for ( int i=0; i<EMBM3D_SIMILARITYSEARCH_METHOD::NUM ; ++i )
      similaritySearchMethodComboBox->addItem ( EMBM3D_SIMILARITYSEARCH_METHOD::toString ( i ).c_str ( ) );
    similaritySearchMethodComboBox->setCurrentIndex ( 2 );
  }
  
  void populateNLMSimilaritySearchMethodComboBox ( ) {
    similaritySearchMethodComboBox->clear ( );
    for ( int i=0; i<EMNLM_SIMILARITYSEARCH_METHOD::NUM ; ++i )
      similaritySearchMethodComboBox->addItem ( EMNLM_SIMILARITYSEARCH_METHOD::toString ( i ).c_str ( ) );
    similaritySearchMethodComboBox->setCurrentIndex ( 2 );
  }
protected slots:
  void selectBM3DFilterOptions ( ) {
    blockSizeSpinBox->setValue ( 16 );
    blockSizeSpinBox->setFilterType ( 0 );
    blockSizeLabel->setText ( "Block size (powers of two: 4, 8, 16, 32):" );
    populateBM3DSimilaritySearchMethodComboBox ( );
    additionalFilterOptionsStackedWidget->setCurrentWidget ( emBM3DFilterOptionsPage );
  }
  
  void selectNLMFilterOptions ( ) {
    blockSizeSpinBox->setValue ( 5 );
    blockSizeSpinBox->setFilterType ( 1 );
    blockSizeLabel->setText ( "Block size (odd: 3, 5, 7, ..., 29, 31):" );
    populateNLMSimilaritySearchMethodComboBox ( );
    additionalFilterOptionsStackedWidget->setCurrentWidget ( emNLMFilterOptionsPage );
  }
};


class QDetectAtomsDialogHelp : public QDialog, protected Ui::DetectAtomsDialogHelp {
  Q_OBJECT
public:
  QDetectAtomsDialogHelp ( QWidget * parent = 0, Qt::WindowFlags f = 0 ) : QDialog ( parent, f ) {
    setupUi ( this );
  }
};

class QDetectAtomsDialog : public QDialog, protected Ui::DetectAtomsDialog {
  Q_OBJECT
protected:
  QDetectAtomsDialogHelp *detectAtomsDialogHelp;
public:
  QDetectAtomsDialog ( QWidget * parent = 0, Qt::WindowFlags f = 0 ) : QDialog ( parent, f ) {
    setupUi ( this );
    
    // Enable help
    helpButton->setIcon ( QIcon ( "://help.png" ) );
    detectAtomsDialogHelp = new QDetectAtomsDialogHelp ( this, f );
    connect ( helpButton, SIGNAL ( clicked ( ) ), detectAtomsDialogHelp, SLOT ( show ( ) ) );
    
    // Enable switching default options on/off
    connect ( useDefaultSegmentationParametersCheckBox, SIGNAL ( clicked ( ) ), this, SLOT ( useDefaultSegmentationParametersCheckboxClicked ( ) ) );
    
    // Set gamma and epsilon to default parameters specified in atomFinder.h
    gammaDoubleSpinBox->setValue ( GAMMA );
    epsilonDoubleSpinBox->setValue ( EPSILON );
  }
  
  int atomType ( ) const {
    if ( singleAtomsRadioButton->isChecked ( ) ) return 0;
    else if ( dumbbellsRadioButton->isChecked ( ) ) return 1;
    else return 2;
  }
  
  double gamma ( ) const {
    return gammaDoubleSpinBox->value ( );
  }
  
  double epsilon ( ) const {
    return epsilonDoubleSpinBox->value ( );
  }
  protected slots:
  void useDefaultSegmentationParametersCheckboxClicked ( ) {
    gammaDoubleSpinBox->setEnabled ( !useDefaultSegmentationParametersCheckBox->isChecked ( ) );
    epsilonDoubleSpinBox->setEnabled ( !useDefaultSegmentationParametersCheckBox->isChecked ( ) );
  }
};


class QQuantitativeMeasuresDialog : public QDialog, protected Ui::QuantitativeMeasuresDialog {
  Q_OBJECT
public:
  QQuantitativeMeasuresDialog ( QWidget * parent = 0, Qt::WindowFlags f = 0 ) : QDialog ( parent, f ) {
    setupUi ( this );
  }
  
  void setPrecisionParameters ( double &PeriodX, double &PeriodY, double &AngleX, double &AngleY ) const {
    PeriodX = periodXDoubleSpinBox->value ( );
    PeriodY = periodYDoubleSpinBox->value ( );
    AngleX = angleXDoubleSpinBox->value ( );
    AngleY = angleYDoubleSpinBox->value ( );
  }
};

#endif // __QTINPUTDIALOGS_H
