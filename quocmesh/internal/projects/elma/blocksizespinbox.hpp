#ifndef __BLOCKSIZESPINBOX_HPP
#define __BLOCKSIZESPINBOX_HPP

#include <QSpinBox>

class BlockSizeSpinBox : public QSpinBox {
  Q_OBJECT
protected:
  unsigned int _filterType;
public:
  BlockSizeSpinBox ( QWidget * parent = 0 ) : QSpinBox ( parent ), _filterType ( 0 ) { }
  
  void setFilterType ( const int FilterType = 0 ) {
    _filterType =  FilterType;
  }
protected:
  virtual QValidator::State	validate ( QString & text, int & pos ) const {
    if ( _filterType == 0 ) {
      int x = text.toInt ( );
      if ( x >= 4 && x <= 32 && ( x & ( x - 1 ) ) == 0 ) return QValidator::Acceptable;
      else if ( x <=3 && x != 2 ) return QValidator::Intermediate;
      else return QValidator::Invalid;
    } else if ( _filterType == 1 ) {
      int x = text.toInt ( );
      if ( x >= 3 && x <= 31 && x % 2 == 1 ) return QValidator::Acceptable;
      else if ( x <= 3 ) return QValidator::Intermediate;
      else return QValidator::Invalid;
    } else return QSpinBox::validate ( text, pos );
  }
};

#endif // __BLOCKSIZESPINBOX_HPP
