#ifndef __PREVIOUSFOLDERFILEDIALOG_H
#define __PREVIOUSFOLDERFILEDIALOG_H

#include <QFileDialog>
#include <QDir>
#include <QFileInfo>
#include <QSettings>

class QPreviousFolderFileDialog : public QFileDialog {
  Q_OBJECT
public:
  QPreviousFolderFileDialog ( ) : QFileDialog ( ) { }

  static QString getOpenFileName ( QWidget * parent = 0, const QString & caption = QString ( ), const QString & dir = QString ( ),
                                   const QString & filter = QString ( ), QString * selectedFilter = 0, Options options = 0 ) {
    QSettings settings;
    QString customDir ( dir );
    if ( customDir.isEmpty ( ) ) {
      if ( QDir ( settings.value ( defaultOpenDirKey ( ) ).toString ( ) ).exists ( ) )
        customDir = settings.value ( defaultOpenDirKey ( ) ).toString ( );
    }

    QString openFileName = QFileDialog::getOpenFileName ( parent, caption, customDir, filter, selectedFilter, options );

    if ( !openFileName.isNull ( ) ) {
      QFileInfo fileInfo ( openFileName );
      QDir currentDir ( fileInfo.absoluteDir ( ) );
      if ( currentDir.exists ( ) ) settings.setValue ( defaultOpenDirKey ( ), currentDir.absolutePath ( ) );
    }

    return openFileName;
  }

  static QString getSaveFileName ( QWidget * parent = 0, const QString & caption = QString ( ), const QString & dir = QString ( ),
                                   const QString & filter = QString ( ), QString * selectedFilter = 0, Options options = 0 ) {
    QSettings settings;
    QString customDir ( dir );
    if ( customDir.isEmpty ( ) ) {
      if ( QDir ( settings.value ( defaultSaveDirKey ( ) ).toString ( ) ).exists ( ) )
        customDir = settings.value ( defaultSaveDirKey ( ) ).toString ( );
    }

    QString saveFileName = QFileDialog::getSaveFileName ( parent, caption, customDir, filter, selectedFilter, options );

    if ( !saveFileName.isNull ( ) ) {
      QFileInfo fileInfo ( saveFileName );
      QDir currentDir ( fileInfo.absoluteDir ( ) );
      if ( currentDir.exists ( ) ) settings.setValue ( defaultSaveDirKey ( ), currentDir.absolutePath ( ) );
    }

    return saveFileName;
  }
protected:
  static QString defaultOpenDirKey ( ) {
    return "default_open_dir";
  }

  static QString defaultSaveDirKey ( ) {
    return "default_save_dir";
  }
};


#endif // __PREVIOUSFOLDERFILEDIALOG_H
