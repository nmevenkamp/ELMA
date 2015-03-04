#include "qtWorkers.hpp"

ProgressWorker::ProgressWorker ( AspectRatioPixmapLabel *ImageLabel, qc::MultiArray<double, qc::QC_2D, 3> &Image )
  : _progressBar ( new QtProgressBar ( this ) ), imageLabel ( ImageLabel ), _image ( Image ), _hasFinished ( false ) { }


void AtomDetectionWorker::performBlockingOperation ( ) {
  AtomFinder<double, aol::FullMatrix<double>, LinearRegressionQR<double>, qc::ScalarArray<double, qc::QC_2D>, qc::MultiArray<double, qc::QC_2D, 3> > atomFinder;
  atomFinder.setCatchCtrlC ( true );
  atomFinder.setProgressBar ( this->_progressBar );
  if ( _atomType == 0 ) atomFinder.getSingleAtomPositions ( _atomPositions, _2DGaussianParameters, _segmented, _electronMicrograph, _gamma, _maxIt, _epsilon );
  else if ( _atomType == 1 ) atomFinder.getDumbbellAtomPositions ( _atomPositions, _2DGaussianParameters, _segmented, _electronMicrograph, _gamma, _maxIt, _epsilon );
  else atomFinder.getAtomPositions ( _atomPositions, _2DGaussianParameters, _segmented, _electronMicrograph, _gamma, _maxIt, _epsilon );
}

void NLMWorker::performBlockingOperation ( ) {
  _options.progressBar = this->_progressBar;
  EMNonLocalMeansFilter<double, qc::ScalarArray<double, qc::QC_2D> > filter;
  filter.setCatchCtrlC ( true );
  filter.apply ( _options, _estimate );
}

void BM3DWorker::performBlockingOperation ( ) {
  _options.progressBar = this->_progressBar;
  EMBM3DFilter<double, qc::ScalarArray<double, qc::QC_2D> > filter;
  filter.setCatchCtrlC ( true );
  filter.apply ( _options, _estimate );
}