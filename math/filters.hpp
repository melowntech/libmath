/**
 * @file filters.hpp
 * @author Ondrej Prochazka <ondrej.prochazka@citationtech.net>
 *
 * FIR filters
 */

#ifndef FILTERS_HPP
#define FILTERS_HPP

namespace math {

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
        if ( fabs( x ) < VERY_SMALL_NUMBER )
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

#endif // FILTERS_HPP

