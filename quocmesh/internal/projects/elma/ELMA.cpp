#include "elmaAddon.hpp"

int main ( int argc, char *argv[] ) {
  // Set some basic information about our application. This influences how Qt stores our settings.
  QCoreApplication::setOrganizationName ( "AICES" );
  QCoreApplication::setOrganizationDomain ( "aices.rwth-aachen.de" );
  QCoreApplication::setApplicationName ( "ELectron Micrograph Analysis" );

  // Check if a path is stored in the application settings. If so, add it to the search path.
  // Such a setting can be added manually with
  // settings.setValue( "path", "/usr/local/macports108/bin" );
  // and is stored persistently in the settings file.
  QSettings settings;
  QString path = settings.value( "path" ).toString();
  if ( path.isEmpty() == false )
    aol::appendSearchPath ( path.toLatin1().constData() );

  ELMAViewerApp app ( argc, argv );
  return app.exec();
}