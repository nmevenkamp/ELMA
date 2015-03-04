#ifndef __ASPECTRATIOPIXMAPLABEL_HPP
#define __ASPECTRATIOPIXMAPLABEL_HPP

#include <qtIncludes.h>
#include <aol.h>

// Based on AspectRatioPixmapLabel from http://stackoverflow.com/questions/8211982/qt-resizing-a-qlabel-containing-a-qpixmap-while-keeping-its-aspect-ratio
// author http://stackoverflow.com/users/999943/phyatt
class AspectRatioPixmapLabel : public QLabel {
  Q_OBJECT
public:
  explicit AspectRatioPixmapLabel ( QWidget *parent = 0 );
  virtual int heightForWidth ( int width ) const;
  virtual QSize sizeHint() const;
  const QPixmap &getPixmap () const;
  const QImage &getOverlay () const;
signals:
  
  public slots:
  void setPixmap ( const QPixmap & );
  void cropPixmap ( const QRect & );
  void setOverlay ( const QImage &, const double );
  void setPermanentOverlay ( const QImage &, const double );
  void clearPermanentOverlay ( );
  void resizeEvent ( QResizeEvent * );
  const QPoint mapToPixmapFromGlobal ( const QPoint &Pos );
private:
  void displayPixmap ( );
  void resetOverlay ( );
  QPixmap _pix, _canvas;
  QImage _overlay, _permanentOverlay;
  double _overlayOpacity, _permanentOverlayOpacity;
};

#endif // __ASPECTRATIOPIXMAPLABEL_HPP