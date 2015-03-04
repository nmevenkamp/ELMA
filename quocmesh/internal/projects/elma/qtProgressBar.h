#ifndef __QTPROGRESSBAR_H
#define __QTPROGRESSBAR_H

#ifdef __GNUC__
#pragma GCC system_header
#endif

#include <platformDependent.h>
#include <ctrlCCatcher.h>
#include <progressBar.h>

#include "qtWorkers.hpp"

class ProgressWorker;

class QtProgressBar : public aol::ProgressBar<> {
protected:
  ProgressWorker *_pWorker;
public:
  QtProgressBar ( ProgressWorker *worker, const char *text = NULL ) : aol::ProgressBar<> ( text ), _pWorker ( worker ) { }
  
  void start ( const int maxValue, const int delta = 1, const int /*incr*/ = 1 );
  
  /** postfix-increment: increment the internal counter and print progress bar afterwards
   */
  void operator++ ( int I );
  
  void setText ( const char *text );
  
  void display ( std::ostream &out ) const {
  }
};

#endif // __QTPROGRESSBAR_H
