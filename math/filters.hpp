/**
 * @file filters.hpp
 * @author Ondrej Prochazka <ondrej.prochazka@citationtech.net>
 *
 * FIR filters
 */

#ifndef MATH_FILTERS_HPP
#define MATH_FILTERS_HPP

#include <cmath>

namespace math {

/**
 * Abstract class providing horizontal image filtering with a static convolution
 * kernel.
 */

#include "volume.hpp"


class FIRFilter_t {

public :

    FIRFilter_t( const uint halforder ) : halforder( halforder ) {
        order = 2 * halforder + 1;
        kernel = new float[ order ];
        std::fill( kernel, kernel + order, 0 );
    }

    ~FIRFilter_t() { delete[] kernel; };

    void dump() {
        for ( uint i = 0; i < order; i++ )
            std::cout << "\t" << i << "\t"
            << std::setprecision( 10 ) << std::fixed << kernel[ i ] << "\n";
    }

    /**
     * Filter a grayscale image (in x direction). GIL image views should be used.
     */
    template <typename SrcView_t, typename DstView_t>
    void filterImage( const SrcView_t & srcView, const DstView_t & dstView ) {

        for ( int y = 0; y < srcView.height(); y++ ) {
            typename SrcView_t::x_iterator srcxit = srcView.row_begin( y );
            typename DstView_t::x_iterator dstxit = dstView.row_begin( y );

            for ( int x = 0; x < srcView.width(); x++ ) {

                *dstxit = gilConvolute( srcxit, x, srcView.width() );
                srcxit++; dstxit++;
            }

        }
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




/**
 * Abstract class for 2D analytic window filter
 */

class Filter2_t {
public :
    /** window sizes */
    float halfwindowX() const { return hwinx; }
    float halfwindowY() const { return hwiny; }

    /** kernel value at given pos */
    virtual float operator() ( float x, float y ) const = 0;

    Filter2_t( float hwinx, float hwiny ) : hwinx( hwinx ), hwiny( hwiny ) {};

protected :
    float hwinx, hwiny;
};

/**
 * 2D sinc function with a hamming window
 */

class LowPassFilter2_t : public Filter2_t {

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
};

} // namespace math

#endif // MATH_FILTERS_HPP

