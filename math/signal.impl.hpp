/**
 * @file signal.hpp
 * @author Ondrej Prochazka <ondrej.prochazka@citationtech.net>
 *
 * 2d signal with forward signal reconstruction capability.
 */

#ifndef MATH_SIGNAL_INMPL_HPP
#define MATH_SIGNAL_INMPL_HPP

#include <cmath>

#include "dbglog/dbglog.hpp"

#include "signal.hpp"
#include "math.hpp"

namespace math {

template <typename FilterType>
void Signal2<FilterType>::sample( double x, double y, double value,
    double periodX, double periodY ) {

    double px = ( x - llx ) / pixelx;
    double py = ( y - lly ) / pixely;
    
    double pcutoffX = 2 * std::max( periodX, pixelx ) / pixelx;
    double pcutoffY = 2 * std::max( periodY, pixely ) / pixely;

    /*LOG( debug ) << "Adding sample at [" << x << "," << y << "] ("
        << px << "," << py << "), periods ["
        << periodX << "," << periodY <<"].";*/
    
    FilterType filter( pcutoffX, pcutoffY );

    double phwinx = filter.halfwinx();
    double phwiny = filter.halfwiny();

    
    int minpx = int( floor( px - phwinx ) );
    int maxpx = int( ceil( px + phwinx ) );
    int minpy = std::max( int( floor( py - phwiny ) ), 0 );
    int maxpy = std::min( int( ceil( py + phwiny ) ), sizeY - 1 );

    if ( type == Nonperiodic ) {
        minpx = std::max( minpx, 0 );
        maxpx = std::min( maxpx, sizeX - 1 );
    }

    for ( int j = minpy; j <= maxpy; j++ )
        for ( int i = minpx; i <= maxpx; i++ ) {
            Cell_s & cell = at( ( i + 10 * sizeX ) % sizeX, j ) ;

            unsigned char quadrant = Cell_s::QUADRANT_NONE;

            bool innerWindow =
                math::ccinterval( px - 0.55 * pcutoffX,
                            px + 0.55 * pcutoffX, (double) i ) &&
                math::ccinterval( py - 0.55 * pcutoffY,
                            py + 0.55 * pcutoffY, (double) j );

            if ( innerWindow ) {

                if ( i <= px  &&  j <= py ) quadrant |= Cell_s::QUADRANT_LL;
                if ( i >= px  &&  j <= py ) quadrant |= Cell_s::QUADRANT_LR;
                if ( i >= px  &&  j >= py ) quadrant |= Cell_s::QUADRANT_UR;
                if ( i <= px  &&  j >= py ) quadrant |= Cell_s::QUADRANT_UL;
            }

            cell.add( value, filter( i - px, j - py ), quadrant );
        }
}

template <typename FilterType>
double Signal2<FilterType>::operator() ( const double x, const double y ) const
{
    // transform coordinates to pixel space
    double px = ( x - llx ) / pixelx;
    double py = ( y - lly ) / pixely;

    if ( type & Xperiodic ) {
        int hedge = int( floor( px / sizeX ) );
        px -= hedge * sizeX;
    }

    if ( fabs( px - 0.0 ) < 1E-15 ) px = 0.0;
    if ( fabs( px - sizeX + 1 ) < 1E-15 ) px = sizeX - 1;
    if ( fabs( py - 0.0 ) < 1E-15 ) py = 0.0;
    if ( fabs( py - sizeY + 1 ) < 1E-15 ) py = sizeY - 1;

    if ( ! ( type & Xperiodic ) )
        if ( ! math::ccinterval( 0.0, sizeX - 1.0, px ) ) {

            LOG( err2 ) << "Out of domain value (" << x << "," << y
                        << ") requested from signal.";
            throw std::runtime_error( "Domain error in Signal2" );
            
        }
            
    if ( ! math::ccinterval( 0.0, sizeY - 1.0, py ) ) {
            LOG( err2 ) << "Out of domain value (" << x << "," << y
                        << ") requested from signal.";
            throw std::runtime_error( "Domain error in Signal2" );
    }

    // reconstruct value
    FilterType filter( 2, 2 );

    int minpx, maxpx, minpy, maxpy;
    
    
    if ( type & Xperiodic ) {
        minpx = int( floor( px - filter.halfwinx() ) );
        maxpx = int( ceil( px + filter.halfwinx() ) );
        
    } else {
        minpx = std::max( int( floor( px - filter.halfwinx() ) ), 0 );
        maxpx = std::min( int( ceil( px + filter.halfwinx() ) ), sizeX - 1 );
    }
    
    minpy = std::max( int( floor( py - filter.halfwiny() ) ), 0 );
    maxpy = std::min( int( ceil( py ) + filter.halfwiny() ), sizeY - 1 );
    
    double valueSum( 0.0 ), weightSum( 0.0 );
    
    for ( int i = minpx; i <= maxpx; i++ )
        for ( int j = minpy; j <= maxpy; j++ ) {

            const Cell_s & cell( at( ( i + 10 * sizeX ) % sizeX, j ) );
            
            if ( cell.defined() ) {
                double weight = filter( i - px, j - py );
                
                valueSum += weight * cell.value();
                weightSum += weight;
            }
        }


    if ( weightSum < 1E-15 ) {
        LOG( err1 ) << "Signal value undefined at (" << x << "," << y << ").";
        throw std::runtime_error( "Value error in Signal2" );
    }

    return valueSum / weightSum;
}

} // namespace math

#endif // MATH_SIGNAL_IMPL_HPP
