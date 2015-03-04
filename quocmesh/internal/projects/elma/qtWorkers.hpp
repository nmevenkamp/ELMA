#ifndef __QTWORKERS_H
#define __QTWORKERS_H

#ifdef __GNUC__
#pragma GCC system_header
#endif

#include "qtProgressBar.h"

#include <aol.h>
#include <ctrlCCatcher.h>
#include <emNonLocalMeansFilter.h>
#include <emBM3DFilter.h>
#include <atomFinder.h>
#include <aspectratiopixmaplabel.hpp>


class QtProgressBar;

class ProgressWorker : public QObject {
  Q_OBJECT
protected:
  QtProgressBar *_progressBar;
  AspectRatioPixmapLabel *imageLabel;
  qc::MultiArray<double, qc::QC_2D, 3> &_image;
  bool _hasFinished;
public:
  ProgressWorker ( AspectRatioPixmapLabel *ImageLabel, qc::MultiArray<double, qc::QC_2D, 3> &Image );
  
  bool hasFinished ( ) const {
    return _hasFinished;
  }
  
  virtual const char* getText ( ) = 0;
protected:
  virtual void performBlockingOperation ( ) = 0;
  
  void emitProgressValueChanged ( int val ) {
    emit progressValueChanged ( val );
  }
  
  void emitMaximumValueChanged ( int maxVal ) {
    emit maximumValueChanged ( maxVal );
  }
  
  void emitLabelTextChanged ( const char* text ) {
    emit labelTextChanged ( QString ( text ) );
  }
public slots:
  virtual void doWork ( ) {
    _hasFinished = false;
    try {
      performBlockingOperation ( );
    } catch ( aol::Exception &el ) {
      el.dump ( );
      imageLabel->setText ( el.getMessage ( ).c_str ( ) );
      _image.reallocate ( 1, 1 );
    }
    _hasFinished = !aol::getCtrlCState ( );
    emit this->finished ( );
  }
signals:
  void finished ( );
  
  void progressValueChanged ( int );
  
  void maximumValueChanged ( int );
  
  void labelTextChanged ( const QString & );
  
  friend class QtProgressBar;
};


class AtomDetectionWorker : public ProgressWorker {
  Q_OBJECT
protected:
  qc::ScalarArray<double, qc::QC_2D> _electronMicrograph;
  aol::MultiVector<double> _atomPositions;
  aol::MultiVector<double> _2DGaussianParameters;
  qc::BitArray<qc::QC_2D> _segmented;
  const int _atomType;
  const double _gamma, _epsilon;
  const int _maxIt;
public:
  AtomDetectionWorker ( AspectRatioPixmapLabel *ImageLabel, qc::MultiArray<double, qc::QC_2D, 3> &Image,
                        qc::ScalarArray<double, qc::QC_2D> &ElectronMicrograph,
                        const int AtomType,
                        const double Gamma, const double MaxIt, const double Epsilon )
    : ProgressWorker ( ImageLabel, Image ), _electronMicrograph ( ElectronMicrograph ), _atomType ( AtomType ), _gamma ( GAMMA ), _maxIt ( MaxIt ), _epsilon ( Epsilon ) { }
  
  const char* getText ( ) {
    return "Detecting atoms";
  }
  
  const qc::ScalarArray<double, qc::QC_2D>& getElectronMicrograph ( ) {
    return _electronMicrograph;
  }
  
  const aol::MultiVector<double>& getAtomPositions ( ) {
    return _atomPositions;
  }
  
  const aol::MultiVector<double>& get2DGaussianParameters ( ) {
    return _2DGaussianParameters;
  }
  
  const qc::BitArray<qc::QC_2D>& getSegmented ( ) {
    return _segmented;
  }
protected:
  void performBlockingOperation ( );
};


class DenoisingWorker : public ProgressWorker {
  Q_OBJECT
protected:
  qc::ScalarArray<double, qc::QC_2D> _estimate;
public:
  DenoisingWorker ( AspectRatioPixmapLabel *ImageLabel, qc::MultiArray<double, qc::QC_2D, 3> &Image )
    : ProgressWorker ( ImageLabel, Image ) { }
  
  const char* getText ( ) {
    return "Denoising";
  }
  
  const qc::ScalarArray<double, qc::QC_2D>& getEstimate ( ) {
    return _estimate;
  }
};

class NLMWorker : public DenoisingWorker {
  Q_OBJECT
protected:
  EMNLMOptions<double, qc::ScalarArray<double, qc::QC_2D> > _options;
public:
  NLMWorker ( AspectRatioPixmapLabel *ImageLabel, qc::MultiArray<double, qc::QC_2D, 3> &Image,
              EMNLMOptions<double, qc::ScalarArray<double, qc::QC_2D> > &Options )
    : DenoisingWorker ( ImageLabel, Image ), _options ( Options ) { }
  
  const char* getText ( ) {
    return "Denoising with NLM";
  }
protected:
  void performBlockingOperation ( );
};

class BM3DWorker : public DenoisingWorker {
  Q_OBJECT
protected:
  EMBM3DOptions<double, qc::ScalarArray<double, qc::QC_2D> > _options;
public:
  BM3DWorker ( AspectRatioPixmapLabel *ImageLabel, qc::MultiArray<double, qc::QC_2D, 3> &Image,
               EMBM3DOptions<double, qc::ScalarArray<double, qc::QC_2D> > &Options )
    : DenoisingWorker ( ImageLabel, Image ), _options ( Options ) { }
  
  const char* getText ( ) {
    return "Denoising with BM3D";
  }
protected:
  void performBlockingOperation ( );
};


#endif // __QTWORKERS_H
