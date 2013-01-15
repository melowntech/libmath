/**
 * @file filters.hpp
 * @author Ondrej Prochazka <ondrej.prochazka@citationtech.net>
 *
 * FIR filters
 */

#ifndef MATH_FILTERS_HPP
#define MATH_FILTERS_HPP

#include <iomanip>
#include <sstream>

#include <cmath>

namespace math {

/**
 * Abstract class providing horizontal image filtering with a static convolution
 * kernel.
 */


class FIRFilter_t {

public :

    FIRFilter_t( const uint halforder ) : halforder( halforder ) {
        order = 2 * halforder + 1;
        kernel = new float[ order ];
        std::fill( kernel, kernel + order, 0 );
    }

    ~FIRFilter_t() { delete[] kernel; };

    std::string dump() {
        std::ostringstream str;
        for ( uint i = 0; i < order; i++ )
            str << "\t" << i << "\t"
            << std::setprecision( 10 ) << std::fixed << kernel[ i ] << "\n";
        return str.str();
    }

    /**
     * Obtain a single filter response value for a given input. The input is represented
     * by an iterator.
     * This is in fact a convolution operation, with filter used as a convolution kernel.
     * xpos and rowsize parameters signal the position within input sequence and the length
     * of the sequence.
     */ 
    template <typename Iterator_t>
    float convolute( const Iterator_t & pos,
                        uint xpos, uint rowSize ) const {

        float weightSum = 0.0;
        float valueSum = 0.0;

        Iterator_t begin = pos - std::min( xpos, halforder );
        Iterator_t end = pos + std::min( rowSize - 1 - xpos, halforder );


        int index = std::max( (int) ( halforder - xpos ), (int) 0 );

        for ( Iterator_t it = begin; it <= end; ++it ) {
            valueSum += it[0] * kernel[ index ];
            weightSum += kernel[ index ];
            index++;
        }

        return valueSum / weightSum;
    }


    /**
     * boost::gil view iterators need a special treatment, even if theya are grayscale views.
     * I don't know why.
     */
    template <typename Iterator_t>
    float gilConvolute( const Iterator_t & pos,
                        uint xpos, uint rowSize ) const {

        float weightSum = 0.0;
        float valueSum = 0.0;

        Iterator_t begin = pos - std::min( xpos, halforder );
        Iterator_t end = pos + std::min( rowSize - 1 - xpos, halforder );


        int index = std::max( (int) ( halforder - xpos ), (int) 0 );

        for ( Iterator_t it = begin; it <= end; ++it ) {
            valueSum += it[0][0] * kernel[ index ];
            weightSum += kernel[ index ];
            index++;
        }

        return valueSum / weightSum;
    }
     

protected :

    uint halforder, order;
    float * kernel;
};

/**
 * Low-pass filter constructed from sinc function by applying a hamming window.
 */

class LowPassFilter_t : public FIRFilter_t {

public :
    LowPassFilter_t( float cutoffPeriod, uint halfwindow )
    : FIRFilter_t( halfwindow ) {

        for ( uint i = 0; i <= halforder; i++ ) {

            float sinci = i == 0 ? 2 / cutoffPeriod :
                1.0 / ( M_PI * i ) * sin( 2 * M_PI * i / cutoffPeriod );
            float hamming = ( order == 1 ) ? 1 :
                0.54 + 0.46 * cos( 2 * M_PI * i / ( order - 1 ) );
            kernel[ halforder - i ]
                = kernel[ halforder + i ] = hamming * sinci;

        }
    }
};

/*
 * Box filter, used for simple average.
 */

class BoxFilter : public FIRFilter_t {

public:
    BoxFilter( uint halfwindow ) : FIRFilter_t( halfwindow ) {

        for ( uint i = 0; i <= halforder; i++ )
            kernel[ halforder - i ] = 1.0 /  ( 2 * halforder + 1 );
    }
};



/**
 * Filter2 is a 2D analytic window filter, used in image processing.
 * It models the following concepts:
 *
 * class Filter2Concept {
 * public:
 *      double operator() ( const double x, const double y ) const;
 *      double halfwinx() const;
 *      double halfwiny() const;
 * };
 *
 *
 * class LowPassFilter2Concept {
 * public :
 *      LowPassFilter2Concept( double cutoffX, double cutoffY )
 *      double operator() ( const double x, const double y ) const;
 *      double halfwinx() const;
 *      double halfwiny() const;
 * }
 */
 
class SincHamming2 {

public:
    SincHamming2( double cutoffX, double cutoffY  )
        : hwinx( cutoffX ), hwiny( cutoffY ),
          cutoffX( cutoffX ), cutoffY( cutoffY ) {}
          
    SincHamming2( double cutoffX, double hwinx, double cutoffY, double hwiny )
        : hwinx( hwinx ), hwiny( hwiny ),
          cutoffX( cutoffX ), cutoffY( cutoffY ) {}

    double operator() ( const double x, const double y ) const {

        return sinc( x, cutoffX ) * hamming( x, hwinx )
            * sinc( y, cutoffY ) * hamming( y, hwiny );
    }

    double halfwinx() const { return hwinx; }
    double halfwiny() const { return hwiny; }
    

private :

    static double sinc( double x, double cutoff )  {
        if ( fabs( x ) < 1E-15 )
            return 2.0 / cutoff;
        else
            return 1.0 / ( M_PI * x ) * sin( 2 * M_PI * x / cutoff );
    }

    static double hamming( double x, double hwin ) {
        return fabs( x ) > hwin ? 0.0 : 0.54 + 0.46 * cos( M_PI * x / hwin );
    }

    double hwinx, hwiny, cutoffX, cutoffY;
};


class Lanczos2 {

public:

    Lanczos2( double cutoffX, double cutoffY )
        : hwinx( cutoffX ), hwiny( cutoffY ),
          cutoffX( cutoffX ), cutoffY( cutoffY ) {}

    Lanczos2( double cutoffX, double hwinx, double cutoffY, double hwiny )
        : hwinx( hwinx ), hwiny( hwiny ),
          cutoffX( cutoffX ), cutoffY( cutoffY ) {}

    double operator() ( const double x, const double y ) const {

        if ( fabs( x ) <= hwinx && fabs( y ) <= hwiny )
            return
                sinc( x, cutoffX ) * sinc( x, 2.0 * hwinx )
                * sinc( y, cutoffY ) * sinc( y, 2.0 * hwiny );
            else
                return 0;
        
    }

    double halfwinx() const { return hwinx; }
    double halfwiny() const { return hwiny; }


private :

    static double sinc( double x, double cutoff )  {
        if ( fabs( x ) < 1E-15 )
            return 2.0 / cutoff;
        else
            return 1.0 / ( M_PI * x ) * sin( 2 * M_PI * x / cutoff );
    }

    double hwinx, hwiny, cutoffX, cutoffY;    
};


class Box2 {

public:
    Box2( double hwinx, double hwiny )
        : hwinx_( hwinx ), hwiny_( hwiny ) {
        value_ = 1.0 / ( hwinx_ * hwiny_ );
    }

    double operator() ( double x, double y ) const {
        if ( fabs( x ) <= hwinx_ && fabs( y ) <= hwiny_ )
            return value_;
        else
            return 0.0;
    }

    double halfwinx() const { return hwinx_; }
    double halfwiny() const { return hwiny_; }
    
private:
    double hwinx_, hwiny_;
    double value_;
    
};

class Linear2 {
public:
    Linear2(double cutoffX, double cutoffY)
        : hwinx(cutoffX / 2.), hwiny(cutoffY / 2.)
        , x2_T(2. / cutoffX), x4_T2(x2_T * x2_T)
        , y2_T(2. / cutoffY), y4_T2(y2_T * y2_T)
    {}

    double operator()(const double x, const double y) const {
        auto ax(std::abs(x));
        auto ay(std::abs(y));
        return (((ax >= hwinx) ? 0. : (x2_T - x4_T2 * ax))
                * ((ay >= hwiny) ? 0. : (y2_T - y4_T2 * ay)));
    }

    double halfwinx() const { return hwinx; }
    double halfwiny() const { return hwiny; }

private:
    double hwinx, hwiny;
    double x2_T, x4_T2; // x coefficients: 2/T and (2/T)^2
    double y2_T, y4_T2; // y coefficients: 2/T and (2/T)^2
};


/** OBSOLETE
 * Abstract class for 2D analytic window filter
 */

/*class Filter2_t {
public :
    float halfwindowX() const { return hwinx; }
    float halfwindowY() const { return hwiny; }

    virtual float operator() ( float x, float y ) const = 0;

    Filter2_t( float hwinx, float hwiny ) : hwinx( hwinx ), hwiny( hwiny ) {};

protected :
    float hwinx, hwiny;
};*/

/**
 * 2D sinc function with a hamming window
 */

/*class LowPassFilter2_t : public Filter2_t {

public:
    LowPassFilter2_t( float cutoffX, float hwinx, float cutoffY, float hwiny )
        : Filter2_t( hwinx, hwiny ), cutoffX( cutoffX ),
            cutoffY( cutoffY ) {}

    virtual float operator() ( float x, float y ) const {

        return sinc( x, cutoffX ) * hamming( x, hwinx )
            * sinc( y, cutoffY ) * hamming( y, hwiny );
    }

private :

    static float sinc( float x, float cutoff )  {
        if ( fabs( x ) < 1E-15 )
            return 2.0 / cutoff;
        else
            return 1.0 / ( M_PI * x ) * sin( 2 * M_PI * x / cutoff );
    }

    static float hamming( float x, float hwin ) {
        return fabs( x ) > hwin ? 0.0 : 0.54 + 0.46 * cos( M_PI * x / hwin );
    }

    float cutoffX, cutoffY;
};*/

} // namespace math

#endif // MATH_FILTERS_HPP

