#include "qtProgressBar.h"

void QtProgressBar::start ( const int maxValue, const int delta, const int /*incr*/ ) {
  aol::ProgressBar<>::start ( maxValue, delta );
  _pWorker->emitMaximumValueChanged ( maxValue );
}

void QtProgressBar::operator++ ( int I ) {
  aol::ProgressBar<>::operator++ ( I );
  if ( !aol::getCtrlCState ( ) )
    _pWorker->emitProgressValueChanged ( this->_value );
}

void QtProgressBar::setText ( const char *text ) {
  aol::ProgressBar<>::setText ( text );
  _pWorker->emitLabelTextChanged ( text );
}