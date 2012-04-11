/**
 * @file geometry.hpp
 * @author Ondrej Prochazka <ondrej.prochazka@citationtech.net>
 *
 * Analytical geometry functions
 */

#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace math {


namespace ublas = boost::numeric::ublas;

/**
 * Find the point where two lines get closest in 3D space.
 * The lines are defined parametrically in euclidian coordinates.
 * Return value is the euclidian metrics of the line distance, plus
 * the two parameter values for each of the lines.
 */

float lineDistance(
    const ublas::vector<float> & p1, const ublas::vector<float> & u1,
    const ublas::vector<float> & p2, const ublas::vector<float> & u2,
    float & r1, float & r2 );


/** normalize vector */
template <class T>
ublas::vector<T> normalize( const ublas::vector<T> & v ) {
    return v / ublas::norm_2( v );
}

/** homogenous coordinates (from euclidian) */
template <class T>
ublas::vector<typename T::value_type> homogenous( const T & src ) {
    
    ublas::vector<typename T::value_type> ret(
        ublas::unit_vector<typename T::value_type>(
            src.size() + 1, src.size() ) );
    ublas::subrange( ret, 0, src.size() ) = src;
    return ret;
}

/** euclidian coordinates (from homogenous) */
template <class T>
ublas::vector<typename T::value_type> euclidian( const T & src ) {
    return ublas::subrange( src, 0, src.size() - 1 ) / src( src.size() - 1 );
}

/** Parametric line, in euclidian 2D
 */
struct Line2_t {
    ublas::vector<double> p, u;

    Line2_t( const ublas::vector<double> p, const ublas::vector<double> u )
        : p( p ), u( u ) {};
};

/**
 * Parametric line, in euclidian 3D
 */

struct Line3_t {
    ublas::vector<double> p, u;

    Line3_t(
        const ublas::vector<double> p = ublas::zero_vector<double>( 3 ),
        const ublas::vector<double> u = ublas::zero_vector<double>( 3 ) )
        : p( p ), u( u ) {}
};

template <typename E, typename T>
static std::basic_ostream<E, T> & operator << (
        std::basic_ostream<E,T> & os,
        const Line3_t & line ) {

    os << line.p << " + t * " << line.u;
    return os;
}


/**
 * Parametric plane, in euclidian 3D
 */

struct Plane3_t {

    ublas::vector<double> p, u, v;

    Plane3_t(
        const ublas::vector<double> p = ublas::zero_vector<double>( 3 ),
        const ublas::vector<double> u = ublas::zero_vector<double>( 3 ),
        const ublas::vector<double> v = ublas::zero_vector<double>( 3 ) )
        : p( p ), u( u ), v( v ) {}
};


template <typename E, typename T>
static std::basic_ostream<E, T> & operator << (
        std::basic_ostream<E,T> & os,
        const Plane3_t & plane ) {

    os << plane.p << " + t * " << plane.u << " + s * " << plane.v;
    return os;
}


/** line and plane intersection */

ublas::vector<double> intersection(
    const Line3_t & line, const Plane3_t & plane );



} // namespace math

#endif //GEOMETRY_HPP
