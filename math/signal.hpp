/**
 * @file signal.hpp
 * @author Ondrej Prochazka <ondrej.prochazka@citationtech.net>
 *
 * 2d signal with forward signal reconstruction capability.
 */

#ifndef MATH_SIGNAL_HPP
#define MATH_SIGNAL_HPP

#include <ostream>

#include <boost/gil/gil_all.hpp>
#include <boost/scoped_array.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "filters.hpp"

namespace gil = boost::gil;
namespace ublas = boost::numeric::ublas;

namespace math {

class Signal2Base {

public :

    enum Type_t {
        Nonperiodic = 0x00,
        Xperiodic = 0x01
    };

    Signal2Base( int sizeX, int sizeY,
        const double llx, const double lly, const double urx,
        const double ury, const Type_t type = Nonperiodic );

    int width() const { return sizeX; }
    int height() const { return sizeY; }

    double pixelX() { return pixelx; }
    double pixelY() { return pixely; }

    struct Cell_s {

        enum {

            QUADRANT_LL = 0x01,
            QUADRANT_LR = 0x02,
            QUADRANT_UL = 0x04,
            QUADRANT_UR = 0x08,

            QUADRANT_NONE = 0x00,
            QUADRANT_ALL = 0x0f

        };

        double valueSum, pwsum, nwsum;
        unsigned char quadrants;

        Cell_s() : valueSum( 0.0 ), pwsum( 0.0 ), nwsum( 0.0 ),
            quadrants( QUADRANT_NONE ) {};

        bool defined() const {

            if ( pwsum - nwsum < 1E-15 ) return false;
            if ( pwsum <= 3 * nwsum ) return false;
            if ( quadrants != QUADRANT_ALL ) return false;

            return true;
        }


        double value() const {
            assert( defined() );
            return valueSum / ( pwsum - nwsum );
        }

        void add( double value, double weight, unsigned char & quadrant ) {
            valueSum += weight * value;
            quadrants |= quadrant;
            if ( weight > 0 ) pwsum += weight; else nwsum -= weight;
        };
    };


    typedef Cell_s * iterator;
    typedef const Cell_s * const_iterator;
    typedef const Cell_s * const_row_iterator;

    class const_col_iterator {

    public:

        const Cell_s & operator[] ( int i ) const {
            return *( cell + step * i ); }
        const_col_iterator & operator++()  { cell += step; return *this; }
        const_col_iterator operator + ( int i ) const {
            return const_col_iterator( cell + step * i, step ); }
        const Cell_s & operator* () const { return *cell; }
        const Cell_s * operator -> () const { return cell; }
        bool operator < ( const const_col_iterator & s ) {
            return cell < s.cell; }
        const_col_iterator( const Cell_s * cell, int step )
            : cell( cell ), step( step ) {}

    private:
        const Cell_s * cell;
        int step;
    };

    Cell_s & at( int i, int j ) { return cells[sizeX * j + i]; }
    const Cell_s & at( int i, int j ) const { return cells[ sizeX * j + i]; }

    iterator begin() { return iterator( cells.get() ); }
    const_iterator begin() const { return iterator( cells.get() ); }
    iterator end() { return iterator( cells.get() + sizeX * sizeY ); }
    const_iterator end() const {
        return const_iterator( cells.get() + sizeX * sizeY ); }
    const_iterator row_begin( int j ) const {
        return iterator( cells.get() + j * sizeX ); }
    const_iterator row_end( int j ) const {
        return const_iterator( cells.get() + ( j + 1 ) * sizeX );  };
    const_col_iterator col_begin( int i ) const {
        return const_col_iterator( cells.get() + i, sizeX ); }
    const_col_iterator col_end( int i ) const {
        return const_col_iterator( cells.get() + sizeX * sizeY + i, sizeX ); };

     /**
      * Provide a signal visualisation in form of an image.
      * White pixels correspond to the maximum value, black pixels to the
      * minimum.
      * The alpha channel is used as a binary indication of the part of the
      * signal that has been reconstructed.
      * The input image needs to have the same dimensions as the sampled signal.
      * This function is intended for diagnostics. */
     void visualize( gil::rgba8_image_t & output ) const;

     /** Debug dump in the form of a Matlab matrix. Undefined cells are stored
      *  as "nan"s.
      */
     void dump(std::ostream &f) const;

     class Transform_t {
     public :
         virtual ublas::vector<double> operator () (
            const ublas::vector<double> & v ) const { return v; }

         virtual ~Transform_t() {}
     };

     typedef std::vector<ublas::vector<double> > QuadList_t;

     void getQuads( QuadList_t & quads, const Transform_t & trafo =
         Transform_t() );

protected:
    int sizeX, sizeY;
    double llx, lly, urx, ury, pixelx, pixely;
    Type_t type;
    boost::scoped_array<Cell_s> cells;
};

template <typename FilterType = math::SincHamming2>
class Signal2 : public Signal2Base {
public :
    Signal2( int sizeX, int sizeY,
        const double llx, const double lly, const double urx,
        const double ury, const Type_t type = Nonperiodic )
        : Signal2Base(sizeX, sizeY, llx, lly,  urx, ury, type)
    {}

    /**
     * provide a signal sample at given position, of given value, and with
     * given sampling period hint. */
    void sample( double x, double y, double value,
        double periodX, double periodY );

    /**
     * return interpolated value of the signal at given position.
     */
    double operator () ( const double x, const double y ) const;
};

} // namespace math

// templated functions implementation
#include "signal.impl.hpp"

#endif // MATH_SIGNAL_HPP
