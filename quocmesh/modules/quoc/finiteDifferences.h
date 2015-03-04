#ifndef __FINITEDIFFERENCES_H
#define __FINITEDIFFERENCES_H

#include <multiArray.h>
#include <ctrlCCatcher.h>
#include <ChanVese.h>
#include <quocTimestepSaver.h>

namespace qc {

/**
 * This struct + static function construction is a workaround to C++ limitations.
 *
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
struct doCalculateBackwardFDDivergence {
  static void apply ( const qc::MultiArray<RealType, Dim> &MArg, qc::ScalarArray<RealType, Dim> &Divergence );
};

template <typename RealType>
struct doCalculateBackwardFDDivergence<RealType, qc::QC_2D> {
  static void apply ( const qc::MultiArray<RealType, qc::QC_2D> &MArg, qc::ScalarArray<RealType, qc::QC_2D> &Divergence ) {
    const int numX = Divergence.getNumX();
    const int numY = Divergence.getNumY();
    // inner nodes
    for ( int y = 1; y < numY - 1; ++y ) {
      for ( int x = 1; x < numX - 1; ++x ) {
        Divergence.set ( x, y, MArg[0].get(x,y) -MArg[0].get(x-1,y) + MArg[1].get(x,y) -MArg[1].get(x,y-1) );
      }
    }

    // top and bottom nodes (without corners)
    for ( int y = 1; y < numY - 1; ++y ) {
      Divergence.set ( 0, y, MArg[0].get(0,y) + MArg[1].get(0,y)  -MArg[1].get(0,y-1) );
      Divergence.set ( numX-1, y, -MArg[0].get(numX-2,y) + MArg[1].get(numX-1,y) -MArg[1].get(numX-1,y-1) );
    }

    // left and right nodes (without corners)
    for ( int x = 1; x < numX-1; ++x ) {
      Divergence.set ( x, 0, MArg[0].get(x,0) -MArg[0].get(x-1,0) + MArg[1].get(x,0) );
      Divergence.set ( x, numY-1, MArg[0].get(x,numY-1) -MArg[0].get(x-1,numY-1) -MArg[1].get(x,numY-2) );
    }

    // corner (x,y) = (0,0)
    Divergence.set ( 0, 0, MArg[0].get(0,0) + MArg[1].get(0,0) );

    // corner (x,y) = (0,numY-1)
    Divergence.set ( 0, numY-1, MArg[0].get(0,numY-1) - MArg[1].get(0,numY-2) );

    // corner (x,y) = (numX-1,0)
    Divergence.set ( numX-1, 0, -MArg[0].get(numX-2,0) + MArg[1].get(numX-1,0) );

    // corner (x,y) = (numX-1,numY-1)
    Divergence.set ( numX-1, numY-1, -MArg[0].get(numX-2,numY-1) - MArg[1].get(numX-1,numY-2) );

    // Straightforward and readable, but less efficient implementation:
    /*
    Divergence.setZero();
    for ( int x = 0; x < numX; ++x ) {
      for ( int y = 0; y < numY; ++y ) {
        if ( x < numX - 1 )
          Divergence.add ( x, y, MArg[0].get(x,y) );
        if ( x > 0 )
          Divergence.add ( x, y, -MArg[0].get(x-1,y) );
        if ( y < numY - 1 )
          Divergence.add ( x, y, MArg[1].get(x,y) );
        if ( y > 0 )
          Divergence.add ( x, y, -MArg[1].get(x,y-1) );
      }
    }
    */
  }
};

//! \todo Rewrite to get rid of all the if statements inside the for loops.
template <typename RealType>
struct doCalculateBackwardFDDivergence<RealType, qc::QC_3D> {
  static void apply ( const qc::MultiArray<RealType, qc::QC_3D> &MArg, qc::ScalarArray<RealType, qc::QC_3D> &Divergence ) {
    const int numX = Divergence.getNumX();
    const int numY = Divergence.getNumY();
    const int numZ = Divergence.getNumZ();
    // Straightforward and readable, but not most efficient implementation:
    Divergence.setZero();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int z = 0; z < numZ; ++z ) {
      for ( int y = 0; y < numY; ++y ) {
        for ( int x = 0; x < numX; ++x ) {
          if ( x < numX - 1 )
            Divergence.add ( x, y, z, MArg[0].get(x,y,z) );
          if ( x > 0 )
            Divergence.add ( x, y, z, -MArg[0].get(x-1,y,z) );
          if ( y < numY - 1 )
            Divergence.add ( x, y, z, MArg[1].get(x,y,z) );
          if ( y > 0 )
            Divergence.add ( x, y, z, -MArg[1].get(x,y-1,z) );
          if ( z < numZ - 1 )
            Divergence.add ( x, y, z, MArg[2].get(x,y,z) );
          if ( z > 0 )
            Divergence.add ( x, y, z, -MArg[2].get(x,y,z-1) );
        }
      }
    }
  }
};

/**
 * Calulates the divergence with backward finite differences (not scaled by the grid width),
 * as used in {An Algorithm for Total Variation Minimization and Applications} by Antonin Chambolle.
 *
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
void calculateBackwardFDDivergence ( const qc::MultiArray<RealType, Dim> &MArg, qc::ScalarArray<RealType, Dim> &Divergence ) {
  doCalculateBackwardFDDivergence<RealType, Dim>::apply ( MArg, Divergence );
}

template <typename RealType, qc::Dimension Dim>
struct doCalculateCentralFDDivergencePBC {
  static void apply ( const qc::MultiArray<RealType, Dim> &MArg, qc::ScalarArray<RealType, Dim> &Divergence );
};

//! \todo Rewrite to get rid of all the if statements inside the for loops.
template <typename RealType>
struct doCalculateCentralFDDivergencePBC<RealType, qc::QC_3D> {
  static void apply ( const qc::MultiArray<RealType, qc::QC_3D> &MArg, qc::ScalarArray<RealType, qc::QC_3D> &Divergence ) {
    const int numX = Divergence.getNumX();
    const int numY = Divergence.getNumY();
    const int numZ = Divergence.getNumZ();
    // Straightforward and readable, but not most efficient implementation:
    Divergence.setZero();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int z = 0; z < numZ; ++z ) {
      for ( int y = 0; y < numY; ++y ) {
        for ( int x = 0; x < numX; ++x ) {
          Divergence.add ( x, y, z, MArg[0].get((x+1) %numX,y,z) );
          Divergence.add ( x, y, z, -MArg[0].get((x-1 + numX)%numX,y,z) );
          Divergence.add ( x, y, z, MArg[1].get(x,(y+1) %numY,z) );
          Divergence.add ( x, y, z, -MArg[1].get(x,(y-1 + numY) % numY,z) );
          Divergence.add ( x, y, z, MArg[2].get(x,y,(z+1) %numZ) );
          Divergence.add ( x, y, z, -MArg[2].get(x,y,(z-1 + numZ) % numZ) );
          Divergence.set ( x, y, z, Divergence.get(x,y,z) * 0.5 );
        }
      }
    }
  }
};

/**
 * Calulates the divergence with central finite differences (not scaled by the grid width),
 * and periodic boundary conditions
 *
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
void calculateCentralFDDivergencePBC ( const qc::MultiArray<RealType, Dim> &MArg, qc::ScalarArray<RealType, Dim> &Divergence ) {
  doCalculateCentralFDDivergencePBC<RealType, Dim>::apply ( MArg, Divergence );
}

/**
 * This struct + static function construction is a workaround to C++ limitations.
 *
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
struct doCalculateForwardFDGradient {
  static void apply ( const qc::ScalarArray<RealType, qc::QC_2D> &Arg, qc::MultiArray<RealType, qc::QC_2D> &Gradient );
};

template <typename RealType>
struct doCalculateForwardFDGradient<RealType, qc::QC_2D> {
  static void apply ( const qc::ScalarArray<RealType, qc::QC_2D> &Arg, qc::MultiArray<RealType, qc::QC_2D> &Gradient ) {
    const int numX = Arg.getNumX();
    const int numY = Arg.getNumY();
    for ( int y = 0; y < numY-1; ++y ) {
      for ( int x = 0; x < numX-1; ++x ) {
        Gradient[0].set ( x, y, Arg.get(x+1,y) - Arg.get(x,y) );
        Gradient[1].set ( x, y, Arg.get(x,y+1) - Arg.get(x,y) );
      }
    }

    for ( int x = 0; x < numX-1; ++x ) {
      Gradient[0].set ( x, numY-1, Arg.get(x+1,numY-1) - Arg.get(x, numY-1) );
      Gradient[1].set ( x, numY-1, 0 );
    }

    for ( int y = 0; y < numY-1; ++y ) {
      Gradient[0].set ( numX-1, y, 0 );
      Gradient[1].set ( numX-1, y, Arg.get(numX-1,y+1) - Arg.get(numX-1,y) );
    }

    Gradient[0].set ( numX-1, numY-1, 0 );
    Gradient[1].set ( numX-1, numY-1, 0 );

    // Straightforward and readable, but less efficient implementation:
    /*
    for ( int x = 0; x < numX; ++x ) {
      for ( int y = 0; y < numY; ++y ) {
        if ( x < numX - 1 )
          Gradient[0].set ( x, y, Arg.get(x+1,y) - Arg.get(x,y) );
        else
          Gradient[0].set ( x, y, 0 );
        if ( y < numY - 1 )
          Gradient[1].set ( x, y, Arg.get(x,y+1) - Arg.get(x,y) );
        else
          Gradient[1].set ( x, y, 0 );
      }
    }
    */
  }
};

//! \todo Rewrite to get rid of all the if statements inside the for loops.
template <typename RealType>
struct doCalculateForwardFDGradient<RealType, qc::QC_3D> {
  static void apply ( const qc::ScalarArray<RealType, qc::QC_3D> &Arg, qc::MultiArray<RealType, qc::QC_3D> &Gradient ) {
    const int numX = Arg.getNumX();
    const int numY = Arg.getNumY();
    const int numZ = Arg.getNumZ();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    // Straightforward and readable, but not most efficient implementation:
    for ( int z = 0; z < numZ; ++z ) {
      for ( int y = 0; y < numY; ++y ) {
        for ( int x = 0; x < numX; ++x ) {
          if ( x < numX - 1 )
            Gradient[0].set ( x, y, z, Arg.get(x+1,y,z) - Arg.get(x,y,z) );
          else
            Gradient[0].set ( x, y, z, 0 );
          if ( y < numY - 1 )
            Gradient[1].set ( x, y, z, Arg.get(x,y+1,z) - Arg.get(x,y,z) );
          else
            Gradient[1].set ( x, y, z, 0 );
          if ( z < numZ - 1 )
            Gradient[2].set ( x, y, z, Arg.get(x,y,z+1) - Arg.get(x,y,z) );
          else
            Gradient[2].set ( x, y, z, 0 );
        }
      }
    }
  }
};

/**
 * Calulates the gradient with forward finite differences (not scaled by the grid width),
 * as used in {An Algorithm for Total Variation Minimization and Applications} by Antonin Chambolle.
 *
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
void calculateForwardFDGradient ( const qc::ScalarArray<RealType, Dim> &Arg, qc::MultiArray<RealType, Dim> &Gradient ) {
  doCalculateForwardFDGradient<RealType, Dim>::apply ( Arg, Gradient );
}

template <typename RealType, qc::Dimension Dim>
struct doCalculateCentralFDGradientPBC{
  static void apply ( const qc::ScalarArray<RealType, qc::QC_2D> &Arg, qc::MultiArray<RealType, qc::QC_2D> &Gradient );
};

//! \todo Rewrite to get rid of all the if statements inside the for loops.
template <typename RealType>
struct doCalculateCentralFDGradientPBC<RealType, qc::QC_3D> {
  static void apply ( const qc::ScalarArray<RealType, qc::QC_3D> &Arg, qc::MultiArray<RealType, qc::QC_3D> &Gradient ) {
    const int numX = Arg.getNumX();
    const int numY = Arg.getNumY();
    const int numZ = Arg.getNumZ();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    // Straightforward and readable, but not most efficient implementation:
    for ( int z = 0; z < numZ; ++z ) {
      for ( int y = 0; y < numY; ++y ) {
        for ( int x = 0; x < numX; ++x ) {
            Gradient[0].set ( x, y, z, 0.5 * (Arg.get((x+1) % numX,y,z) - Arg.get((x-1 + numX) % numX,y,z) ));
            Gradient[1].set ( x, y, z, 0.5 * (Arg.get(x,(y+1) %numY,z) - Arg.get(x,(y-1 + numY) % numY,z) ));
            Gradient[2].set ( x, y, z, 0.5 * (Arg.get(x,y,(z+1) % numZ) - Arg.get(x,y,(z-1 + numZ) % numZ) ));
        }
      }
    }
  }
};

/**
 * Calulates the gradient with central finite differences (not scaled by the grid width),
 * and periodic boundary conditions.
 *
 * \author Berkels
 */
template <typename RealType, qc::Dimension Dim>
void calculateCentralFDGradientPBC ( const qc::ScalarArray<RealType, Dim> &Arg, qc::MultiArray<RealType, Dim> &Gradient ) {
  doCalculateCentralFDGradientPBC<RealType, Dim>::apply ( Arg, Gradient );
}

/**
 * Collects very basic stuff for iterative algorithms.
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class TVAlgorithmBase {
protected:
  typedef typename ConfiguratorType::RealType RealType;
  const typename ConfiguratorType::InitType &_grid;
  const RealType _gamma;
  int _maxIterations;
  RealType _stopEpsilon;
  const DefaultArraySaver<RealType, ConfiguratorType::Dim>* _pStepSaver;
  bool _quietMode;
public:

  TVAlgorithmBase ( const typename ConfiguratorType::InitType &Initializer,
                    const RealType Gamma,
                    const int MaxIterations,
                    const RealType StopEpsilon )
    : _grid ( Initializer ),
      _gamma ( Gamma ),
      _maxIterations ( MaxIterations ),
      _stopEpsilon ( StopEpsilon ),
      _pStepSaver ( NULL ),
      _quietMode ( false ) {}

  void setMaxIterations ( const int MaxIterations ) {
    _maxIterations = MaxIterations;
  }

  int getMaxIterations ( ) const {
    return _maxIterations;
  }

  void setStopEpsilon ( const RealType StopEpsilon ) {
    _stopEpsilon = StopEpsilon;
  }

  RealType getStopEpsilon ( ) const {
    return _stopEpsilon;
  }

  void setStepSaverReference ( const DefaultArraySaver<RealType, ConfiguratorType::Dim> &StepSaver ) {
    _pStepSaver = &StepSaver;
  }

  const DefaultArraySaver<RealType, ConfiguratorType::Dim> *getStepSaverPointer ( ) const {
    return _pStepSaver;
  }
  
  void setQuietMode ( bool qmode ) {
    _quietMode = qmode;
  }
};

/**
 * \brief Abstract base class for two phase Mumford Shah segmentation.
 *
 * Minimizes the quadratic Esedoglu like model
 * \f[ \min_{u} \gamma\int_\Omega g|\nabla u|dx + \int_\Omega f_1u^2+f_2(1-u)^2 dx, \f]
 * where \f$ f_1,f_2 \f$ are indicator functions for the two different regions to be segmented.
 *
 * Uses a finite difference scheme closely based on {An Algorithm
 * for Total Variation Minimization and Applications} by Antonin Chambolle:
 * \f[ p^{k+1} = \left(p^k+\tau h^2\nabla\frac{\mathrm{div}p^k-2f_2/\gamma}{2(f_1+f_2)}\right)
 *    /\left(1+\frac{\tau h^2}{g}\left|\nabla\frac{\mathrm{div}p^k-2f_2/\gamma}{2(f_1+f_2)}\right|\right) \f]
 * where \f$ u^k=\frac{2f_2-\gamma\mathrm{div}p^k}{2(f_1+f_2)} \f$,
 * \f$ h \f$ is the grid size, and \f$ \tau \f$ can be taken as \f$ \frac18 \f$
 * as shown by Chambolle.
 *
 * The interface function generateIndicatorFunction, in which \f$ f_1 \f$ and \f$ f_2 \f$ are
 * defined, needs to be implemented in the derived class. Likewise \f$ g \f$ can be defined
 * in the derived class by implementing the interface function generateEdgeWeight.
 *
 * \author Berkels
 */
template <typename ConfiguratorType>
class TwoPhaseMSSegmentor : public TVAlgorithmBase<ConfiguratorType> {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
protected:
  RealType _tau;
  const ArrayType* _pIndicator[2];
  virtual void prepareIndicatorFunctionGeneration ( ) const {}
  virtual void generateIndicatorFunction ( const int IndicatorNumber, ArrayType &IndicatorFunction ) const {
    if ( _pIndicator[IndicatorNumber] == NULL )
      throw ( aol::Exception ( aol::strprintf ( "Indicator %d not set and generateIndicatorFunction() not overloaded.", IndicatorNumber ), __FILE__, __LINE__ ) );

    IndicatorFunction = *( _pIndicator[IndicatorNumber] );
  }
  virtual void generateEdgeWeight ( ArrayType & edgeWeight ) const {
    edgeWeight.setAll ( aol::ZOTrait<RealType>::one );
  }

public:
  TwoPhaseMSSegmentor ( const typename ConfiguratorType::InitType &Initializer,
                        const RealType Gamma )
    : TVAlgorithmBase<ConfiguratorType> ( Initializer, Gamma, 1000, 0.01 ),
      _tau ( 0.25 ) {
    _pIndicator[0] = NULL;
    _pIndicator[1] = NULL;
  }

  virtual ~TwoPhaseMSSegmentor () {}

  void calcPrimalFromDual ( const qc::MultiArray<RealType, ConfiguratorType::Dim> &Dual,
                            ArrayType &Primal,
                            const ArrayType &Indicator1Plus2,
                            const ArrayType &Indicator2 ) const {
    qc::calculateBackwardFDDivergence<RealType, ConfiguratorType::Dim> ( Dual, Primal );
    Primal *= -0.5 * this->_gamma / this->_grid.H();
    Primal += Indicator2;
    Primal /= Indicator1Plus2;
  }

  void segment ( ArrayType &Segmentation, qc::MultiArray<RealType, ConfiguratorType::Dim> * PDual = NULL ) const {
    prepareIndicatorFunctionGeneration ( );
    doSegment ( Segmentation, PDual );
  }

private:
  //! Can assume that prepareIndicatorFunctionGeneration ( ) has just been called.
  virtual void doSegment ( ArrayType &Segmentation, qc::MultiArray<RealType, ConfiguratorType::Dim> *PDual ) const {
    ArrayType indicator2 ( this->_grid );
    ArrayType indicator1Plus2 ( this->_grid );
    ArrayType edgeWeight ( this->_grid );

    generateIndicatorFunction ( 1, indicator2 );
    generateIndicatorFunction ( 0, indicator1Plus2 );
    ArrayType indicator2OverGamma ( indicator2 );
    generateEdgeWeight ( edgeWeight );
    indicator1Plus2 += indicator2;
    indicator2OverGamma *= static_cast<RealType> ( this->_grid.H() ) / this->_gamma;

    qc::MultiArray<RealType, ConfiguratorType::Dim> pOld ( this->_grid );
    if ( PDual != NULL )
      pOld = *PDual;
    qc::MultiArray<RealType, ConfiguratorType::Dim> pNew ( this->_grid );
    ArrayType temp ( this->_grid );
    const int numPrimalDofs = temp.size();
    qc::MultiArray<RealType, ConfiguratorType::Dim> gradient ( this->_grid );

    aol::ProgressBar<> progressBar ( "Segmenting" );
    if ( !this->_quietMode ) {
      progressBar.start ( this->_maxIterations );
      progressBar.display ( cerr );
    }

    for ( int fixpointIterations = 0; fixpointIterations < this->_maxIterations; fixpointIterations++ ) {
      qc::calculateBackwardFDDivergence<RealType, ConfiguratorType::Dim> ( pOld, temp );
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for ( int j = 0; j < numPrimalDofs; ++j ) {
        temp[j] = ( temp[j] * 0.5 - indicator2OverGamma[j] ) / indicator1Plus2[j];
      }
      qc::calculateForwardFDGradient<RealType, ConfiguratorType::Dim> ( temp, gradient );

      typename ConfiguratorType::VecType gradientVec;
#ifdef _OPENMP
#pragma omp parallel for firstprivate ( gradientVec )
#endif
      for ( int j = 0; j < numPrimalDofs; ++j ) {
        for ( int i = 0; i < ConfiguratorType::Dim; ++i )
          gradientVec[i] = gradient[i][j];
        const RealType gradientVecNorm = gradientVec.norm();
        for ( int i = 0; i < ConfiguratorType::Dim; ++i ) {
          pNew[i][j] = ( pOld[i][j] + _tau * gradientVec[i] ) / ( 1 + _tau * gradientVecNorm / edgeWeight[j] );
        }
      }
      pOld -= pNew;
      const RealType change = pOld.norm();
      pOld = pNew;

      if ( this->_pStepSaver ) {
        calcPrimalFromDual ( pNew, temp, indicator1Plus2, indicator2 );
        this->_pStepSaver->saveStep ( temp, fixpointIterations );
      }

      // If the change is small enough, we consider the gradient descent to have converged.
      if ( change < this->_stopEpsilon )
        break;

      if ( !this->_quietMode ) progressBar++;
    }
    if ( !this->_quietMode ) progressBar.finish();
    calcPrimalFromDual ( pOld,  Segmentation, indicator1Plus2, indicator2 );
    if ( PDual != NULL )
      *PDual = pOld;
  }

public:
  void setTau( const RealType Tau ) {
    _tau = Tau;
  }

  void setIndicatorReference ( const int IndicatorNumber, const ArrayType &IndicatorFunction ) {
    _pIndicator[IndicatorNumber] = &IndicatorFunction;
  }
};

template <typename ConfiguratorType>
class WeightedTwoPhaseMSSegmentor : public TwoPhaseMSSegmentor<ConfiguratorType> {

public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;

  WeightedTwoPhaseMSSegmentor ( const typename ConfiguratorType::InitType &grid,
                                const RealType gamma,
                                const ArrayType & edgeWeight )
    : TwoPhaseMSSegmentor<ConfiguratorType> ( grid, gamma ),
      _edgeWeight ( edgeWeight )
  {}

  virtual void generateEdgeWeight ( ArrayType & edgeWeight ) const {
    edgeWeight = _edgeWeight;
  }

protected:
  const ArrayType & _edgeWeight;
};

/**
 * Piecewise constant Mumford Shah image segmentation with support for scalar
 * and vector values images, e.g. color images or vector fields.
 *
 * \author Berkels
 */
template <typename ConfiguratorType, int ImageDimension, typename BaseClass = TwoPhaseMSSegmentor<ConfiguratorType> >
class PiecewiseConstantTwoPhaseMSSegmentor : public BaseClass {
public:
  typedef typename ConfiguratorType::RealType RealType;
  typedef typename ConfiguratorType::ArrayType ArrayType;
protected:
  const aol::MultiVector<RealType> &_imageMVec;
  aol::MultiVector<RealType> _meanValues;
private:
  int _outerIterations;
  aol::ProgressBar<> *_progressBar;
  bool _catchCtrlC;
  mutable sigfunc _previousCtrlCHandler;
public:

  PiecewiseConstantTwoPhaseMSSegmentor ( const typename ConfiguratorType::InitType &Initializer,
                                         const RealType Gamma,
                                         aol::MultiVector<RealType> &ImageMVec )
    : BaseClass ( Initializer, Gamma ),
      _imageMVec ( ImageMVec ),
      _meanValues ( 2, ImageMVec.numComponents ( ) ),
      _outerIterations ( 5 ),
      _catchCtrlC ( false ) {
    _meanValues[0][0] = 0.;
    _meanValues[1][0] = 1.;
  }

  virtual ~PiecewiseConstantTwoPhaseMSSegmentor() {}

protected:
  virtual void generateIndicatorFunction ( const int IndicatorNumber, ArrayType &IndicatorFunction ) const {
    if ( ( IndicatorNumber >= 2 ) || ( IndicatorNumber < 0 ) )
      throw ( aol::OutOfBoundsException ( "PiecewiseConstantTwoPhaseMSSegmentor only defines two indicator functions.", __FILE__, __LINE__ ) );

    const int imageDim = _imageMVec.numComponents ( );
    const RealType shift = 0.5;
    for ( int i = 0; i < IndicatorFunction.size(); ++i) {
      RealType indicator = 0.;
      for ( int j = 0; j < imageDim; j++ )
        // The mean values are assigned to the segments consistently to
        // aol::ClassicalChanVeseEnergyMulti (when used for binary segmentation)
        indicator += aol::Sqr( _imageMVec[j][i] - _meanValues[1-IndicatorNumber][j] );
      IndicatorFunction[i] = indicator + shift;
    }
  }
  virtual void updateGrayValues ( const aol::Vector<RealType> &CurrentSegmentation ) {
    aol::IdentityFunction<RealType>  heavisideFunction;
    aol::ClassicalChanVeseMeanValueUpdater<ConfiguratorType, aol::IdentityFunction<RealType>, 1, ImageDimension> grayValueUpdater
      ( this->_grid, heavisideFunction, _imageMVec );

    aol::MultiVector<RealType> levelsetFunctions ( 1, CurrentSegmentation.size() );
    levelsetFunctions[0] = CurrentSegmentation;

    // The "shift" in generateIndicatorFunction causes a severe loss of contrast
    // in CurrentSegmentation (according to the theory the 0.5 levelset that we are
    // interested in should be unaffected though). To properly calculate the gray
    // values we need to account for this here: The values start in [0,1].
    // First, shift to [-0.5,0.5]
    levelsetFunctions.addToAll ( -0.5 );
    // Rescale such that -0.5 and/or 0.5 is attained (increasing contrast but leaving 0 unaffected).
    levelsetFunctions /= 2*levelsetFunctions.getMaxAbsValue();
    // Finally, shift back to [0,1].
    levelsetFunctions.addToAll ( 0.5 );

    // This is an alternative to the code above (neither always better nor always worse).
    //levelsetFunctions.threshold( 0.5, 0., 1. );

    aol::MultiVector<RealType> newGrayValuesMVec ( 2, ImageDimension );

    grayValueUpdater.update( levelsetFunctions, newGrayValuesMVec );
    
    for ( int j = 0; j < ImageDimension; j++ ) {
      _meanValues[0][j] = newGrayValuesMVec[0][j];
      _meanValues[1][j] = newGrayValuesMVec[1][j];
    }
    if ( !this->_quietMode ) cerr << _meanValues << endl;
  }
public:
  void segmentAndAdjustGrayValues ( ArrayType &Segmentation, qc::MultiArray<RealType, ConfiguratorType::Dim> * PDual = NULL ) {
    setCtrlCHandler ( );
    for ( int i = 0; i < _outerIterations && !wantsInterrupt ( ) ; ++i ) {
      this->segment ( Segmentation, PDual );
      this->updateGrayValues ( Segmentation );
    }
    unsetCtrlCHandler ( );
  }

  const aol::MultiVector<RealType>& getMeanValuesReference () const {
    return _meanValues;
  }

  aol::MultiVector<RealType>& getMeanValuesReference () {
    return _meanValues;
  }

  void setOuterIterations ( const int OuterIterations ) {
    _outerIterations = OuterIterations;
  }
  
  void setCatchCtrlC ( bool catchCtrlC ) {
    _catchCtrlC = catchCtrlC;
  }
protected:
  void setCtrlCHandler () const {
    if (_catchCtrlC)
      _previousCtrlCHandler = signal ( InterruptSignal, aol::ctrlCHandler );
  }
  
  void unsetCtrlCHandler () const {
    if (_catchCtrlC)
      signal ( InterruptSignal, _previousCtrlCHandler );
  }
  
  bool wantsInterrupt() const {
    if (!_catchCtrlC || !aol::getCtrlCState())
      return false;
    else
      return true;
  }
};

} // namespace qc

#endif // __FINITEDIFFERENCES_H
