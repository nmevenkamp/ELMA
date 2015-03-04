#include <gnuplotter.h>

namespace aol {

bool runGnuplot ( const char *GnuplotCommandFileName ) {
  string systemCommand;
  // Since the Windows version of gnuplot finally also has an exe called "gnuplot" the ifdef
  // shouldn't be necessary anymore, it won't hurt to keep it here for now though.
#ifdef WIN32
  systemCommand = "gnuplot.exe ";
#else
  systemCommand = "gnuplot ";
#endif
  systemCommand += GnuplotCommandFileName;
  const bool failed = ( system ( systemCommand.c_str() ) != EXIT_SUCCESS );
  if ( failed )
    cerr << "aol::runGnuplot: Calling gnuplot returned an error.\n";
  return !failed;
}

}
