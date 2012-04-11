/*
 * signal.cpp
 */

#include "signal.hpp"

#include <cmath>


Signal_t::Signal_t( int sizeX, int sizeY,
    const double llx, const double lly, const double urx,
    const double ury, const Type_t type )
    : sizeX( sizeX ), sizeY( sizeY ),
      llx( llx ), lly( lly ), urx( urx ), ury( ury ),
      type( type ) {

    pixelx = ( urx - llx ) / sizeX;
    pixely = ( ury - lly ) / sizeY;

    cells = new Cell_s[ sizeX * sizeY ];
}

Signal_t::~Signal_t() {
    delete[] cells;
}

void Signal_t::sample( double x, double y, double value,
    double periodX, double periodY ) {

    double px = ( x - llx ) / pixelx;
    double py = ( y - lly ) / pixely;

    double pcutoffX = 2 * std::max( periodX, pixelx ) / pixelx;
    double pcutoffY = 2 * std::max( periodY, pixely ) / pixely;
    double phwinx = pcutoffX;
    double phwiny = pcutoffY;

    LowPassFilter2_t filter( pcutoffX, phwinx, pcutoffY, phwiny );

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
            Cell_s & cell = at( ( i + sizeX ) % sizeX, j ) ;

            unsigned char quadrant = Cell_s::QUADRANT_NONE;

            bool innerWindow =
                ccinterval( px - 0.55 * pcutoffX,
                            px + 0.55 * pcutoffX, (double) i ) &&
                ccinterval( py - 0.55 * pcutoffY,
                            py + 0.55 * pcutoffY, (double) j );

            if ( innerWindow ) {

                if ( i < px  &&  j < py ) quadrant = Cell_s::QUADRANT_LL;
                if ( i >= px  &&  j < py ) quadrant = Cell_s::QUADRANT_LR;
                if ( i >= px  &&  j >= py ) quadrant = Cell_s::QUADRANT_UR;
                if ( i < px  &&  j >= py ) quadrant = Cell_s::QUADRANT_UL;
            }

            /*double weight = 1.0 / ( 15.0 * (
                sqr( ( i - px ) / std::min( periodX, periodY )  )
                + sqr( ( j - py ) / std::min( periodX, periodY ) ) ) + 1.0 );*/

            cell.add( value, filter( i - px, j - py ), quadrant );
            //cell.add( value, weight, quadrant );
            /*if ( ( ( i + sizeX ) % sizeX ) == 95 && j == 220 )
            std::cout << px << "\t" << py << "\t" << filter( i - px, j - py )
                << "\t" << value << "\t" << periodY << std::endl; */
        }
}

void Signal_t::visualize( gil::rgba8_image_t & output ) const {

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
    if ( fabs( maxValue - minValue ) < VERY_SMALL_NUMBER )
        minValue = maxValue - 1.0;


    //std::cout << "min = " << minValue << ", max = " << maxValue << "\n";

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

void Signal_t::getQuads( QuadList_t & quads, const Transform_t & trafo  ) {

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
