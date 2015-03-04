#include "aspectratiopixmaplabel.hpp"
#include <aol.h>

// Based on AspectRatioPixmapLabel from http://stackoverflow.com/questions/8211982/qt-resizing-a-qlabel-containing-a-qpixmap-while-keeping-its-aspect-ratio
// author http://stackoverflow.com/users/999943/phyatt
AspectRatioPixmapLabel::AspectRatioPixmapLabel ( QWidget *parent ) : QLabel ( parent ) {
  this->setMinimumSize ( 1, 1 );
}

void AspectRatioPixmapLabel::setPixmap ( const QPixmap & p ) {
  _pix = p;
  resetOverlay ( );
  displayPixmap ( );
}

void AspectRatioPixmapLabel::cropPixmap ( const QRect & r ) {
  if ( _pix.rect ( ).contains ( r ) ) {
    _pix = _pix.copy ( r );
    _permanentOverlay = _permanentOverlay.copy ( r );
  }
  displayPixmap ( );
}
  
void AspectRatioPixmapLabel::displayPixmap ( ) {
  _canvas = QPixmap ( _pix.size ( ) );
  QPoint origin ( 0, 0 );
  QPainter painter ( &_canvas );
  painter.drawPixmap ( origin, _pix );
  painter.setOpacity ( _permanentOverlayOpacity );
  painter.drawImage ( origin, _permanentOverlay );
  painter.setOpacity ( _overlayOpacity );
  painter.drawImage ( origin, _overlay );
  painter.end ( );
  
  // This takes care of rescaling the new pixmap to display it.
  resizeEvent ( NULL );
}

void AspectRatioPixmapLabel::setOverlay ( const QImage & p, const double Opacity ) {
  _overlay = p;
  _overlayOpacity = Opacity;
  displayPixmap ( );
}

void AspectRatioPixmapLabel::setPermanentOverlay ( const QImage & p, const double Opacity ) {
  _permanentOverlay = p;
  _permanentOverlayOpacity = Opacity;
  displayPixmap ( );
}

void AspectRatioPixmapLabel::resetOverlay ( ) {
  _overlay = QImage ( _pix.size ( ), QImage::Format_ARGB32_Premultiplied );
  _overlay.fill ( qRgba(0,0,0,0) );
  _overlayOpacity = 0;
}

void AspectRatioPixmapLabel::clearPermanentOverlay ( ) {
  _permanentOverlay = QImage ( _pix.size ( ), QImage::Format_ARGB32_Premultiplied );
  _permanentOverlay.fill ( qRgba(0,0,0,0) );
  _permanentOverlayOpacity = 0;
  displayPixmap ( );
}

int AspectRatioPixmapLabel::heightForWidth ( int width ) const {
  return static_cast<qreal> ( _pix.height() * width ) / _pix.width();
}

QSize AspectRatioPixmapLabel::sizeHint() const {
  const int w = this->width();
  return QSize ( w, heightForWidth ( w ) );
}

void AspectRatioPixmapLabel::resizeEvent ( QResizeEvent * /*e*/ ) {
  if ( text().isEmpty() ) {
    QSize size = this->size();
#if QT_VERSION >= 0x050000
    // Make use of the full resolution the display offers.
    size *= this->devicePixelRatio();
#endif
    QLabel::setPixmap ( _canvas.scaled ( size, Qt::KeepAspectRatio, Qt::FastTransformation ) );
  }
}

const QPixmap &AspectRatioPixmapLabel::getPixmap () const {
  return _pix;
}

const QImage &AspectRatioPixmapLabel::getOverlay () const {
  return _overlay;
}

const QPoint AspectRatioPixmapLabel::mapToPixmapFromGlobal ( const QPoint &Pos ) {
  double scalingFactor;
  QPoint scaledPixMapCorner ( 0, 0 );
  if ( static_cast<double> ( this->rect ( ).width ( ) ) / _pix.width ( ) < static_cast<double> ( this->rect ( ).height ( ) ) / _pix.height ( ) ) {
    scalingFactor = static_cast<double> ( this->rect ( ).width ( ) ) / _pix.width ( );
    scaledPixMapCorner.ry ( ) += ( this->rect ( ).height ( ) - _pix.height ( ) * scalingFactor ) / 2;
  } else {
    scalingFactor = static_cast<double> ( this->rect ( ).height ( ) ) / _pix.height ( );
    scaledPixMapCorner.rx ( ) += ( this->rect ( ).width ( ) - _pix.width ( ) * scalingFactor ) / 2;
  }
  QPoint pixelPos = mapFromParent ( Pos );
  pixelPos -= scaledPixMapCorner;
  pixelPos.setX ( floor ( pixelPos.x ( ) / scalingFactor ) );
  pixelPos.setY ( floor ( pixelPos.y ( ) / scalingFactor ) );
  return pixelPos;
}