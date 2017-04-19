/**
 * Copyright (c) 2017 Melown Technologies SE
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * *  Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * *  Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
/*
 * signal.cpp
 */

#include "signal.hpp"

#include "math.hpp"
#include <cmath>

#include "dbglog/dbglog.hpp"

namespace math {

Signal2Base::Signal2Base( int sizeX, int sizeY,
    const double llx, const double lly, const double urx,
    const double ury, const Type_t type )
    : sizeX( sizeX ), sizeY( sizeY ),
      llx( llx ), lly( lly ), urx( urx ), ury( ury ),
      type( type ), cells(new Cell_s[ sizeX * sizeY ])
{

    LOG( info1 )
        << "Created signal with extents ["
        << llx << "," << lly << "]->[" << urx << "," << ury << "]";
          
    if ( type & Xperiodic )
        pixelx = ( urx - llx ) / sizeX;
    else
        pixelx = ( urx - llx ) / ( sizeX - 1 );
    
    pixely = ( ury - lly ) / ( sizeY - 1 );
}

void Signal2Base::visualize( gil::rgba8_image_t & output ) const {

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

void Signal2Base::dump(std::ostream &f) const
{
    f << "[";
    for (int j = 0; j < sizeY; j++)
    {
        for (int i = 0; i < sizeX; i++)
        {
            const Cell_s &c(at(i, j));
            if (c.defined())
                f << c.value() << " ";
            else
                f << "nan ";
        }
        f << "\n";
    }
    f << "];\n";
}

void Signal2Base::getQuads( QuadList_t & quads, const Transform_t & trafo  )
{

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
