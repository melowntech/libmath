/**
 * @file geometry.hpp
 * @author Ondrej Prochazka <ondrej.prochazka@citationtech.net>
 *
 * Analytical geometry functions
 */

#ifndef MATH_GEOMETRY_HPP
#define MATH_GEOMETRY_HPP

#include "geometry_core.hpp"

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
ublas::vector<typename T::value_type> normalize( const T & v ) {
    return v / ublas::norm_2( v );
}

/** normalize vector */
template <class T>
Point2_<T> normalize( const Point2_<T> & v ) {
    return v / ublas::norm_2( v );
}

/** normalize vector */
template <class T>
Point3_<T> normalize( const Point3_<T> & v ) {
    return v / ublas::norm_2( v );
}

/** homogeneous coordinates (from euclidian) */
template <class T>
ublas::vector<typename T::value_type> homogeneous( const T & src ) {
    
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

/** cross product, in euclidian 3D */

template<typename T, typename U>
ublas::vector<typename T::value_type> crossProduct(
    const T & u,
    const U & v ) {

    assert( u.size() == 3 && v.size() == 3 );
    ublas::vector<typename T::value_type> retval( 3 );
    retval(0) = u(1) * v(2) - v(1) * u(2);
    retval(1) = -u(0) * v(2) + v(0) * u(2);
    retval(2) = u(0) * v(1) - v(0) * u(1);
    return retval;
}

template <typename T>
Point3_<T> crossProduct(const Point3_<T> & u, const Point3_<T> & v )
{
    return Point3_<T>(u(1) * v(2) - v(1) * u(2)
                     , -u(0) * v(2) + v(0) * u(2)
                     , u(0) * v(1) - v(0) * u(1));
}

template<typename T>
Point3_<T> crossProduct(const Point3_<T> & u, const ublas::vector<T> & v )
{
    assert(v.size() == 3);
    return Point3_<T>(u(1) * v(2) - v(1) * u(2)
                     , -u(0) * v(2) + v(0) * u(2)
                     , u(0) * v(1) - v(0) * u(1));
}

template<typename T>
Point3_<T> crossProduct(const ublas::vector<T> & u, const Point3_<T> & v )
{
    assert(u.size() == 3);
    return Point3_<T>(u(1) * v(2) - v(1) * u(2)
                     , -u(0) * v(2) + v(0) * u(2)
                     , u(0) * v(1) - v(0) * u(1));
}


/** Parametric line, in euclidian 2D
 */
struct Line2 {
    Point2 p, u;

    Line2( const ublas::vector<double> p, const ublas::vector<double> u )
        : p( p ), u( u ) {};
};

/**
 * Parametric line, in euclidian 3D
 */

struct Line3 {
    Point3 p, u;

    Line3(
        const ublas::vector<double> p = ublas::zero_vector<double>( 3 ),
        const ublas::vector<double> u = ublas::zero_vector<double>( 3 ) )
        : p( p ), u( u ) {}
};

template <typename E, typename T>
static std::basic_ostream<E, T> & operator << (
        std::basic_ostream<E,T> & os,
        const Line3 & line ) {

    os << line.p << " + t * " << line.u;
    return os;
}


/**
 * Parametric plane, in euclidian 3D
 */

struct Plane3 {

    Point3 p, u, v;

    Plane3(
        const ublas::vector<double> p = ublas::zero_vector<double>( 3 ),
        const ublas::vector<double> u = ublas::zero_vector<double>( 3 ),
        const ublas::vector<double> v = ublas::zero_vector<double>( 3 ) )
        : p( p ), u( u ), v( v ) {}


    template <class Matrix>
    Plane3 transform( const Matrix & trafo ) const {

        Plane3 tplane;
        
        tplane.p = euclidian( prod( trafo, homogeneous( p ) ) );
        tplane.u = euclidian( prod( trafo, homogeneous( p + u ) ) ) - tplane.p;
        tplane.v = euclidian( prod( trafo, homogeneous( p + v ) ) ) - tplane.p;

        return tplane;
        
    }
};




template <typename E, typename T>
static std::basic_ostream<E, T> & operator << (
        std::basic_ostream<E,T> & os,
        const Plane3 & plane ) {

    os << plane.p << " + t * " << plane.u << " + s * " << plane.v;
    return os;
}


/** line and plane intersection */

ublas::vector<double> intersection(
    const Line3 & line, const Plane3 & plane );



} // namespace math

#endif // MATH_GEOMETRY_HPP
