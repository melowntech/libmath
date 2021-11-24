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
/**
 * @file geometry.hpp
 * @author Ondrej Prochazka <ondrej.prochazka@citationtech.net>
 *
 * Analytical geometry functions
 */

#ifndef MATH_GEOMETRY_HPP
#define MATH_GEOMETRY_HPP

#include "geometry_core.hpp"

#include <algorithm>
#include <array>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace math {

namespace ublas = boost::numeric::ublas;


typedef std::array<math::Point2d, 3> Triangle2d;
typedef std::array<math::Point3d, 3> Triangle3d;

typedef std::vector<Triangle2d> Triangles2d;
typedef std::vector<Triangle3d> Triangles3d;

typedef Points2d Polygon; // single CCW ring, not closed
typedef std::vector<Polygon> MultiPolygon; // multiple rings, holes CW

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
inline ublas::vector<typename T::value_type> normalize( const T & v ) {
    return v / ublas::norm_2( v );
}

/** normalize vector */
template <class T>
inline Point2_<T> normalize( const Point2_<T> & v ) {
    return v / ublas::norm_2( v );
}

/** normalize vector */
template <class T>
inline Point3_<T> normalize( const Point3_<T> & v ) {
    return v / ublas::norm_2( v );
}

/** length of the vector */
template <typename T>
inline double length( const T & v ) {
    return ublas::norm_2(v);
}

/** angle between two 3D vectors (in radians)
 * See Kahan W.: How Futile are Mindless Assessments of Roundoff in
 * Floating-Point Computation?, page 46
 * (https://people.eecs.berkeley.edu/~wkahan/Mindless.pdf)
 *  **/
template <class T>
inline T angle(const Point3_<T>& v1, const Point3_<T>& v2)
{
    return 2
           * std::atan2(length(v1 * length(v2) - length(v1) * v2),
                        length(v1 * length(v2) + length(v1) * v2));
}

/** homogeneous coordinates (from euclidian) */
template <class T>
inline ublas::vector<T, ublas::bounded_array<T, 4> >
homogeneous( const Point3_<T> & src )
{
    ublas::vector<T, ublas::bounded_array<T, 4> > dst(4);
    dst(0) = src(0); dst(1) = src(1); dst(2) = src(2);
    dst(3) = T(1);
    return dst;
}

/** homogeneous coordinates (from euclidian) */
template <class T>
inline ublas::vector<typename T::value_type> homogeneous( const T & src ) {

    ublas::vector<typename T::value_type> ret(
        ublas::unit_vector<typename T::value_type>(
            src.size() + 1, src.size() ) );
    ublas::subrange( ret, 0, src.size() ) = src;
    return ret;
}

/** euclidian coordinates (from homogenous) */
template <class T>
inline Point3_<T>
euclidian( const ublas::vector<T, ublas::bounded_array<T, 4> > & src )
{
    auto div(src(3));
    return Point3_<T>(src(0) / div, src(1) / div, src(2) / div);
}

/** euclidian coordinates (from homogenous) */
template <class T>
inline ublas::vector<typename T::value_type> euclidian( const T & src ) {
    return ublas::subrange( src, 0, src.size() - 1 ) / src( src.size() - 1 );
}

/** cross product, in euclidian 3D */

template<typename T, typename U>
inline ublas::vector<typename T::value_type> crossProduct(
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
inline Point3_<T> crossProduct(const Point3_<T> & u, const Point3_<T> & v )
{
    return Point3_<T>(u(1) * v(2) - v(1) * u(2)
                     , -u(0) * v(2) + v(0) * u(2)
                     , u(0) * v(1) - v(0) * u(1));
}

template<typename T>
inline Point3_<T> crossProduct(const Point3_<T> & u, const ublas::vector<T> & v )
{
    assert(v.size() == 3);
    return Point3_<T>(u(1) * v(2) - v(1) * u(2)
                     , -u(0) * v(2) + v(0) * u(2)
                     , u(0) * v(1) - v(0) * u(1));
}

template<typename T>
inline Point3_<T> crossProduct(const ublas::vector<T> & u, const Point3_<T> & v )
{
    assert(u.size() == 3);
    return Point3_<T>(u(1) * v(2) - v(1) * u(2)
                     , -u(0) * v(2) + v(0) * u(2)
                     , u(0) * v(1) - v(0) * u(1));
}

/** 2D version of cross product; result is not an vector but crossproduct's
 *  z-component.
 */
template <typename T>
inline T crossProduct(const Point2_<T> & u, const Point2_<T> & v )
{
    return u(0) * v(1) - v(0) * u(1);
}

//! return area of a triangle in 2D plane
template <typename T>
inline double triangleArea(const Point2_<T>& a, const Point2_<T>& b,
                           const Point2_<T>& c)
{
    return std::abs(crossProduct(Point2(b - a), Point2(c - a))) * 0.5;
}

//! return area of a triangle in 3D space
template <typename T>
inline double triangleArea(const Point3_<T>& a, const Point3_<T>& b,
                           const Point3_<T>& c)
{
    return norm_2(crossProduct(b - a, c - a)) * 0.5;
}

#ifdef MATH_HAS_OPENCV
    //! return area of a triangle in 2D plane
    template <typename T>
    inline double triangleArea(const cv::Point_<T>& a, const cv::Point_<T>& b,
                               const cv::Point_<T>& c)
    {
        return triangleArea(Point2_<T>(a), Point2_<T>(b), Point2_<T>(c));
    }

    //! return area of a triangle in 3D space
    template <typename T>
    inline double triangleArea(const cv::Point3_<T>& a, const cv::Point3_<T>& b,
                               const cv::Point3_<T>& c)
    {
        return triangleArea(Point3_<T>(a), Point3_<T>(b), Point3_<T>(c));
    }
#endif // MATH_HAS_OPENCV

/** Returns a positive number if the sequence of points {a, b, c} turns
 *  counter-clockwise in the XY plane (negative number otherwise).
 *  Zero (within numerical tolerance) means collinear.
 */
template <typename T>
inline T ccw(const Point2_<T> &a, const Point2_<T> &b, const Point2_<T> &c)
{
    return (b(0) - a(0))*(c(1) - a(1)) - (b(1) - a(1))*(c(0) - a(0));
}


/** Parametric line, in euclidian 2D
 */

template <typename T>
struct Line2_
{
    T p, u;

    Line2_(const T p, const T u) : p(p), u(u) {};
};

typedef Line2_<Point2> Line2;

/** Returns a point where two lines `l1` and `l2` intersect.
 */
template <typename T>
inline T lineIntersection(const Line2_<T>& l1, const Line2_<T>& l2)
{
    T pd = l2.p - l1.p; // diff of origins
    auto num = -l2.u(1) * pd(0) + l2.u(0) * pd(1);
    auto den = l2.u(0) * l1.u(1) - l1.u(0) * l2.u(1);
    return l1.p + (num / den) * l1.u;
}

/**
 * Line segment represented by start and end points
 */
template <typename T>
struct Segment2_
{
    T p1, p2;

    Segment2_(const T p1, const T p2) : p1(p1), p2(p2) {};
};

typedef Segment2_<Point2> Segment2;

/**
 * Parametric line, in euclidian 3D
 */

struct Line3 {
    Point3 p, u;

    Line3(const Point3 p = Point3(), const Point3 u = Point3() )
        : p( p ), u( u ) {}

    /** Returns line's point at given parameter (t)
     */
    Point3 point(double t) const { return p + u * t; }
};

template <typename E, typename T>
inline std::basic_ostream<E, T> & operator << (
        std::basic_ostream<E,T> & os,
        const Line3 & line ) {

    os << line.p << " + t * " << line.u;
    return os;
}

Point3 midpoint( const Line3 & line1, const Line3 & line2
               , double minAngleCos = 0.9962 );

double pointLineDistance(const Point3 &p, const Line3 &line);


/**
 * Parametric plane, in euclidian 3D
 * Legacy representation of plane by 3 points.
 */
namespace legacy {
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
inline std::basic_ostream<E, T> & operator << (
        std::basic_ostream<E,T> & os,
        const Plane3 & plane ) {

    os << plane.p << " + t * " << plane.u << " + s * " << plane.v;
    return os;
}
}  // namespace legacy


/** line and plane intersection */

Point3 intersection( const Line3 & line, const legacy::Plane3 & plane );

/** line and plane intersection
 *  instead of point returns 3 coefficients:
 *      * lines t-parameter
 *      * planes t-parameter
 *      * planes s-parameter
 */
Point3 intersectionParams(const Line3 &line, const legacy::Plane3 &plane);


/**
* Plane represented in general form
* Representation as (a,b,c,d) in equation: ax+by+cz+d=0
*/
struct Plane3
{
    math::Point3 n_;  // vector with (a,b,c) i.e. normal vector
    double d_;

    // Plane from parameters
    Plane3(double a, double b, double c, double d) : n_(a, b, c), d_(d) { }

    // Plane from point and normal
    Plane3(const math::Point3d& pt, const math::Point3d& n)
        : n_(n),
          d_(-ublas::inner_prod(pt, n)) {};

    // Plane from three points
    Plane3(const math::Point3d& p1,
           const math::Point3d& p2,
           const math::Point3d& p3)
        : Plane3(p1, crossProduct(p2 - p1, p3 - p1)) {};

    // Plane from legacy representation (point-normal packed as Line3)
    Plane3(const math::Line3& l) : Plane3(l.p, l.u) {};

    // Returns a plane with opposite orientation
    Plane3 opposite() const {
        Plane3 pl2(*this);
        pl2.n_ *= -1;
        pl2.d_ *= -1;
        return pl2;
    }
};

/**
 * Returns a distance of point to a plane (perpendicular)
*/
double pointPlaneDistance(const Point3 &p, const Plane3 &plane);

/**
 * Returns an orthogonal projection of a point to a plane
 */
Point3 pointPlaneProjection(const Point3 &p, const Plane3 &plane);

/**
 * Returns a point of intersection between plane and line
 */
Point3 linePlaneIntersection(const Line3 &l, const Plane3 &plane);

/**
 * Returns a line of intersection between two planes
 */
Line3 planeIntersection(const Plane3 &p1, const Plane3 &p2);

/**
 * Returns a point of intersection between three planes
 */
Point3 planeIntersection(const Plane3 &p1, const Plane3 &p2, const Plane3 &p3);

/**
 * Returns a measure of triangular polyface regularity. Value of 1.0
 * indicates an equilateral triangle, value of 0.0 indicates a triangle
 * with at least one degenerate edge.
 */
double polygonRegularity(
    const Point3 & v0, const Point3 & v1, const Point3 & v2  );

/**
 * Returns a measure of quad polyface regularity. Value of 1.0
 * indicates a square, value of 0.0 indicates a triangle
 * with at least one degenerate edge.
 */
double polygonRegularity(
    const Point3 & v0, const Point3 & v1, const Point3 & v2, const Point3 & v3 );

/**
 * Returns whether the triangle and rectagle collide.
 */
bool triangleRectangleCollision( math::Point2 triangle[3]
                               , math::Point2 ll, math::Point2 ur);

/**
 * Convert cartesian coordinates (r) to barycentric with respect to the (a,b,c)
 * triangle. See http://en.wikipedia.org/wiki/Barycentric_coordinate_system
 */
math::Point3 barycentricCoords(const math::Point2 &r, const math::Point2 &a
                               , const math::Point2 &b, const math::Point2 &c);
inline
math::Point3 barycentricCoords(const math::Point2 &r, const Triangle2d &tri) {
    return barycentricCoords(r, tri[0], tri[1], tri[2]);
}

namespace detail {
    template <typename T, typename Q> struct ExtentsTypeTraits;

    template <typename T> struct ExtentsTypeTraits<T, Point2_<T>> {
        typedef Extents2_<T> type;
    };

    template <typename T> struct ExtentsTypeTraits<T, Point3_<T>> {
        typedef Extents3_<T> type;
    };

#ifdef MATH_HAS_OPENCV
    template <typename T> struct ExtentsTypeTraits<T, cv::Point_<T>> {
        typedef Extents2_<T> type;
    };

    template <typename T> struct ExtentsTypeTraits<T, cv::Point3_<T>> {
        typedef Extents3_<T> type;
    };
#endif
} // namespace detail

template <typename Iterator>
inline auto computeExtents(Iterator begin, Iterator end)
    -> typename detail::ExtentsTypeTraits
    <typename std::iterator_traits<Iterator>::value_type::value_type
    , typename std::iterator_traits<Iterator>::value_type>::type
{
    typedef typename detail::ExtentsTypeTraits
        <typename std::iterator_traits<Iterator>::value_type::value_type
         , typename std::iterator_traits<Iterator>::value_type>::type Extents;
    typedef typename Extents::point_type point_type;

    if (begin == end) {
        return Extents(InvalidExtents{});
    }

    Extents extents(*begin++);
    std::for_each(begin, end
                  , [&extents](const point_type &p) {
                      update(extents, p);
                  });
    return extents;
}

/** Simplified version for STL container.
 */
template <typename Container>
inline auto computeExtents(const Container &c)
    -> typename detail::ExtentsTypeTraits
    <typename Container::value_type::value_type
     , typename Container::value_type>::type
{
    return computeExtents(c.begin(), c.end());
}

// prefered natural aliases for above

template <typename Iterator>
inline auto extents(Iterator begin, Iterator end)
    -> typename detail::ExtentsTypeTraits
    <typename std::iterator_traits<Iterator>::value_type::value_type
    , typename std::iterator_traits<Iterator>::value_type>::type
{
    return computeExtents<Iterator>(begin,end);
}

template <typename Container>
inline auto extents(const Container &c)
    -> typename detail::ExtentsTypeTraits
    <typename Container::value_type::value_type
     , typename Container::value_type>::type
{
    return computeExtents<Container>(c);
}

} // namespace math

#endif // MATH_GEOMETRY_HPP
