/*
 * signal.cpp
 */

#include "signal.hpp"

#include "math.hpp"
#include "filters.hpp"
#include <cmath>

#include <dbglog/dbglog.hpp>

namespace math {

Signal2::Signal2( int sizeX, int sizeY,
    const double llx, const double lly, const double urx,
    const double ury, const Type_t type )
    : sizeX( sizeX ), sizeY( sizeY ),
      llx( llx ), lly( lly ), urx( urx ), ury( ury ),
      type( type ) {

    LOG( info2 )
        << "Created signal with extents ["
        << llx << "," << lly << "]->[" << urx << "," << ury << "]";
          
    if ( type & Xperiodic )
        pixelx = ( urx - llx ) / sizeX;
    else
        pixelx = ( urx - llx ) / ( sizeX - 1 );
    
    pixely = ( ury - lly ) / ( sizeY - 1 );

    cells = new Cell_s[ sizeX * sizeY ];
}

Signal2::~Signal2() {
    delete[] cells;
}

void Signal2::sample( double x, double y, double value,
    double periodX, double periodY ) {

    double px = ( x - llx ) / pixelx;
    double py = ( y - lly ) / pixely;
    
    double pcutoffX = 2 * std::max( periodX, pixelx ) / pixelx;
    double pcutoffY = 2 * std::max( periodY, pixely ) / pixely;

    /* LOG( debug ) << "Adding sample at [" << x << "," << y << "] ("
        << px << "," << py << "), periods ["
        << periodX << "," << periodY <<"].";*/
    
    math::LowPassFilter2_t filter( pcutoffX, pcutoffX, pcutoffY, pcutoffY );

    double phwinx = filter.halfwindowX();
    double phwiny = filter.halfwindowY();

    
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

void Signal2::visualize( gil::rgba8_image_t & output ) const {

    gil::rgba8_view_t ov = gil::view( output );
    gil::fill_pixels( ov, gil::rgba8_pixel_t( 0x0, 0x0, 0x0, 0x0 ) );

    // find value extents
    double minValue( 0.0 ), maxValue( 0.0 );
    bool isdef( false );

    for ( const_iterator it = begin(); it < end(); ++it )
        if ( it->defined() ) {
            if ( isdef ) {
                minValue  = std::min( minValue, it->value() );
                maxValue  = std::max( maxValue, it->value() );
            } else {
                minValue = maxValue = it->value();
                isdef = true;
            }
        }

    if ( ! isdef  ) return;
    if ( fabs( maxValue - minValue ) < 1E-15 )
        minValue = maxValue - 1.0;


    // iterate through pixels
    for ( int j = 0; j < sizeY; j++ ) {
        const_iterator sit = row_begin( sizeY - 1 - j );
        gil::rgba8_view_t::x_iterator dit = ov.row_begin( j );

        for ( int i = 0; i < sizeX; i++ ) {

            gil::rgba8_pixel_t & pixel = *dit;

            if ( sit->defined() ) {
                pixel[0] = pixel[1] = pixel[2]
                    = int( round( ( sit->value() - minValue )
                        / ( maxValue - minValue ) * 0xff ) );
                pixel[3] = 0xff;
            } else {
                pixel[0] = pixel[1] = pixel[2] = pixel[3] = 0x00;
            }


            /*if ( sit->defined() )
                std::cout << i << "\t" << j << "\t"
                    << sit->pwsum << "\t" << sit->nwsum << "\t" << sit->value()
                    << "\t" << int( pixel[0] ) << std::endl;*/

            dit++; sit++;
        }
    }
}

double Signal2::operator() ( const double x, const double y ) const {

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
    LowPassFilter2_t filter( 2, 2, 2, 2 );

    int minpx, maxpx, minpy, maxpy;
    
    
    if ( type & Xperiodic ) {
        minpx = int( floor( px - filter.halfwindowX() ) );
        maxpx = int( ceil( px + filter.halfwindowX() ) );
        
    } else {
        minpx = std::max( int( floor( px - filter.halfwindowX() ) ), 0 );
        maxpx = std::min( int( ceil( px + filter.halfwindowX() ) ), sizeX - 1 );
    }
    
    minpy = std::max( int( floor( py - filter.halfwindowY() ) ), 0 );
    maxpy = std::min( int( ceil( py ) + filter.halfwindowY() ), sizeY - 1 );
    
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

void Signal2::getQuads( QuadList_t & quads, const Transform_t & trafo  ) {

    int xupper = ( type & Xperiodic ) ? sizeX : sizeX - 1;

    for ( int i = 0; i < xupper; i++ ) {

        const_col_iterator lit = col_begin( i );
        const_col_iterator rit = col_begin( ( i + 1 ) % sizeX );

        for ( int j = 0; j < sizeY - 1; j++ ) {

            bool isdef = ( lit->defined() && rit->defined()
                && ( lit + 1 )->defined() && ( rit + 1 )->defined() );

            if ( isdef ) {

                ublas::vector<double> ll(3), lr(3), ul(3), ur(3);

                ll[0] = llx + i * pixelx;
                lr[0] = llx + ( i + 1 ) * pixelx;
                ur[0] = llx + ( i + 1 ) * pixelx;
                ul[0] = llx + i * pixelx;

                ll[1] = lly + j * pixely;
                lr[1] = lly + j * pixely;
                ur[1] = lly + ( j + 1 ) * pixely;
                ul[1] = lly + ( j + 1 ) * pixely;

                ll[2] = lit->value();
                lr[2] = rit->value();
                ur[2] = ( rit + 1 )->value();
                ul[2] = ( lit + 1 )->value();

                quads.push_back( trafo( ll ) ); quads.push_back( trafo( lr ) );
                quads.push_back( trafo( ur ) ); quads.push_back( trafo( ul ) );
            }

            ++lit; ++rit;
        }
    }
}

} // namespace math
