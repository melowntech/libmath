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
 * @file geometry_core.hpp
 * @author Ondrej Prochazka <ondrej.prochazka@citationtech.net>
 * @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 * Basic structures for geometric modeling
 */

#ifndef MATH_GEOMETRY_CORE_HPP
#define MATH_GEOMETRY_CORE_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_match.hpp>
#include <boost/spirit/include/qi_match_auto.hpp>
#include <boost/spirit/include/qi_alternative.hpp>

#include <boost/rational.hpp>

#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include <limits>

#ifdef MATH_HAS_OPENCV
#include <opencv2/core/core.hpp>
#endif

#include "utility/streams.hpp"

namespace math {

namespace ublas = boost::numeric::ublas;

struct GeometryError : std::runtime_error {
    GeometryError(const std::string &msg) : std::runtime_error(msg) {}
};

struct NoIntersectError : GeometryError {
    NoIntersectError(const std::string &msg) : GeometryError(msg) {}
};

/* sizes */

template <class T>
struct Size2_ {
    typedef T value_type;

    Size2_() : width( 0 ), height( 0 ) {};

    Size2_( const T & width, const T & height )
        : width( width ), height( height ) {};

    T width, height;

    template <typename U>
    explicit Size2_(const Size2_<U> &s)
        : width(s.width), height(s.height)
    {}

    template <typename U>
    explicit Size2_(const Size2_<boost::rational<U> > &s)
        : width(boost::rational_cast<T>(s.width))
        , height(boost::rational_cast<T>(s.height))
    {}

    bool operator== (const Size2_<T>& s) const {
        return width == s.width && height == s.height;
    }
    bool operator!= (const Size2_<T>& s) const {
        return !operator==(s);
    }
};


typedef Size2_<int> Size2i;
typedef Size2_<long long> Size2ll;
typedef Size2_<double> Size2f;
typedef Size2i Size2;
typedef Size2_<boost::rational<long long> > Size2r;

namespace detail {

/** Area type. Must be wide enough to hold result of (width * height).
 *
 *  Defaults to long long (should hold any sane area of integral sizes.
 *
 *  Floating types have identical area type.
 */
template <typename T> struct AreaType { typedef long long type; };
template<> struct AreaType<float> { typedef float type; };
template<> struct AreaType<double> { typedef double type; };
template<> struct AreaType<long double> { typedef long double type; };

} // namespace detail

template <typename T>
typename detail::AreaType<T>::type area(const Size2_<T> &size)
{
    return ((typename detail::AreaType<T>::type)(size.width) * size.height);
}

template <typename T>
bool empty(const Size2_<T> &size)
{
    return (size.width <= 0) || (size.height <= 0);
}

template <typename T> struct Size2SimpleReader { math::Size2_<T> *value; };

/** Used to safely read Size2_<T> in format WxH or A (equivalent to AxA) from
 *  input stream.
 */
template <typename T>
Size2SimpleReader<T> size2SimpleReader(math::Size2_<T> &value);

template <class T>
struct Size3_ {
    typedef T value_type;

    Size3_() : width( 0 ), height( 0 ), depth( 0 ) {};

    Size3_( const T & width, const T & height, const T & depth)
        : width( width ), height( height ), depth(depth) {};

    T width, height, depth;

    template <typename U>
    explicit Size3_(const Size3_<U> &s)
        : width(s.width), height(s.height), depth(s.depth)
    {}

    template <typename U>
    explicit Size3_(const Size3_<boost::rational<U> > &s)
        : width(boost::rational_cast<T>(s.width))
        , height(boost::rational_cast<T>(s.height))
        , depth(boost::rational_cast<T>(s.depth))
    {}

    bool operator== (const Size3_<T>& s) const {
        return width == s.width && height == s.height && depth == s.depth;
    }
    bool operator!= (const Size3_<T>& s) const {
        return !operator==(s);
    }

    T & operator()(int idx) {
        switch(idx) {
            case 0 : return width;
            case 1 : return height;
            case 2 : return depth;
            default :
                throw GeometryError("Bad index to Size3.");
        }
    }

    const T & operator()(int idx) const {
        switch(idx) {
            case 0 : return width;
            case 1 : return height;
            case 2 : return depth;
            default :
                throw GeometryError("Bad index to Size3.");
        }
    }
};

typedef Size3_<int> Size3i;
typedef Size3_<long long> Size3ll;
typedef Size3_<double> Size3f;
typedef Size3i Size3;
typedef Size3_<boost::rational<long long> > Size3r;

template <typename T> struct Size3SimpleReader { math::Size3_<T> *value; };

/** Used to safely read Size3_<T> in format WxHxD or WxA (equivalent WxAxA) or A
 *  (equivalent to AxAxA) from input stream.
 */
template <typename T>
Size3SimpleReader<T> size3SimpleReader(math::Size3_<T> &value);

/* viewports */

template <typename T>
struct Viewport2_ {
    typedef T value_type;
    typedef Size2_<T> size_type;

    value_type width;
    value_type height;
    value_type x;
    value_type y;

    Viewport2_() : width(1000), height(1000), x(0), y(0) {}

    Viewport2_(value_type width, value_type height
               , value_type x = 0, value_type y = 0)
        : width(width), height(height), x(x), y(y) {};

    Viewport2_(const size_type &size, value_type x = 0, value_type y = 0)
        : width( size.width ), height( size.height ), x( x ), y( y ) {};

    size_type size() const { return math::Size2i(width, height); }

    template <typename U>
    explicit Viewport2_(const Viewport2_<U> &v)
        : width(v.width), height(v.height), x(v.x), y(v.y)
    {}
};

typedef Viewport2_<int> Viewport2i;
typedef Viewport2_<double> Viewport2f;
typedef Viewport2i Viewport2;

template <class T> class Point2_;

template <typename T1, typename T2>
inline bool inside(const Viewport2_<T1> &v, const Point2_<T2> &p)
{
    return ((p(0) >= v.x) && (p(0) <= (v.x + v.width))
            && (p(1) >= v.y) && (p(1) <= (v.y + v.height)));
}

template <typename T1, typename T2>
inline bool inside(const Viewport2_<T1> &v, T2 x, T2 y)
{
    return ((x >= v.x) && (x <= (v.x + v.width))
            && (y >= v.y) && (y <= (v.y + v.height)));
}

/* points and point vectors */

template <class T>
class Point2_ : public ublas::vector<T, ublas::bounded_array<T, 2> > {
public:
    Point2_( const T & x = 0.0, const T & y = 0.0 )
        : ublas::vector<T, ublas::bounded_array<T, 2> >(2) {
        (*this)(0) = x; (*this)(1) = y;
    }

    template <class AE>
    Point2_( const ublas::vector_expression<AE> & op )
        : ublas::vector<T, ublas::bounded_array<T, 2> >(2) {
        ublas::vector_assign<ublas::scalar_assign>(*this, op );
    }

    Point2_( const ublas::vector<T> & op )
        : ublas::vector<T, ublas::bounded_array<T, 2> >(2) {
        ublas::vector_assign<ublas::scalar_assign>( *this, op );
    }

#ifdef MATH_HAS_OPENCV
    Point2_( const cv::Point_<T> & op )
        : ublas::vector<T, ublas::bounded_array<T, 2> >(2) {
        (*this)(0) = op.x; (*this)(1) = op.y;
    }
#endif

#if 0
    // quite dangerous, imho, should be phased out, use euclidian instead
    Point2_( const cv::Point3_<T> & op ) {
        (*this)(0) = op.x / op.z; (*this)(1) = op.y / op.z;
    }
#endif

    template <typename U>
    explicit Point2_(const Point2_<U> &op)
        : ublas::vector<T, ublas::bounded_array<T, 2> >(2)
    {
        (*this)(0) = op(0); (*this)(1) = op(1);
    }

    bool operator== (const Point2_<T>& p) const {
        return (*this)(0) == p(0) && (*this)(1) == p(1);
    }
    bool operator!= (const Point2_<T>& p) const {
        return !operator==(p);
    }
};

template <class T>
class Point3_ : public ublas::vector<T, ublas::bounded_array<T, 3> > {
public:
    Point3_( const T & x = 0.0, const T & y = 0.0, const T & z = 0.0 )
        : ublas::vector<T, ublas::bounded_array<T, 3> >(3) {
        (*this)(0) = x; (*this)(1) = y; (*this)(2) = z;
    }

    template <class AE>
    Point3_( const ublas::vector_expression<AE> & op )
        : ublas::vector<T, ublas::bounded_array<T, 3> >(3) {
        ublas::vector_assign<ublas::scalar_assign>( *this, op );
    }

    Point3_( const ublas::vector<T> & op )
        : ublas::vector<T, ublas::bounded_array<T, 3> >(3) {
        ublas::vector_assign<ublas::scalar_assign>( *this, op );
    }

    Point3_(const Point3_<T> &) = default;
    Point3_ & operator=(const Point3_<T> &) = default;

    // make point movable
    Point3_( Point3_<T> && op )
        : ublas::vector<T, ublas::bounded_array<T, 3> >(3) {
        ublas::vector_assign<ublas::scalar_assign>( *this, op );
    }

    Point3_ & operator=(Point3_<T> &&op) {
        ublas::vector_assign<ublas::scalar_assign>( *this, op );
        return *this;
    }

    template <typename U>
    explicit Point3_(const Point3_<U> &op)
        : ublas::vector<T, ublas::bounded_array<T, 3> >(3)
    {
        (*this)(0) = op(0); (*this)(1) = op(1); (*this)(2) = op(2);
    }

    bool operator== (const Point3_<T>& p) const {
        return (*this)(0) == p(0) && (*this)(1) == p(1) && (*this)(2) == p(2);
    }
    bool operator!= (const Point3_<T>& p) const {
        return !operator==(p);
    }
};

template <class T>
class Point4_ : public ublas::vector<T, ublas::bounded_array<T, 4> > {
public:
    Point4_( const T & x = 0.0, const T & y = 0.0
           , const T & z = 0.0, const T & w = 0.0 )
        : ublas::vector<T, ublas::bounded_array<T, 4> >(4) {
        (*this)(0) = x; (*this)(1) = y; (*this)(2) = z; (*this)(3) = w;
    }

    template <class AE>
    Point4_( const ublas::vector_expression<AE> & op )
        : ublas::vector<T, ublas::bounded_array<T, 4> >(4) {
        ublas::vector_assign<ublas::scalar_assign>( *this, op );
    }

    Point4_( const ublas::vector<T> & op )
        : ublas::vector<T, ublas::bounded_array<T, 4> >(4) {
        ublas::vector_assign<ublas::scalar_assign>( *this, op );
    }

    bool operator== (const Point4_<T>& p) const {
        return (*this)(0) == p(0) && (*this)(1) == p(1) &&
               (*this)(2) == p(2) && (*this)(3) == p(3);
    }
    bool operator!= (const Point4_<T>& p) const {
        return !operator==(p);
    }
};

typedef Point2_<int> Point2i;
typedef Point2_<long long> Point2ll;
typedef Point2_<float> Point2f;
typedef Point2_<double> Point2d;
typedef Point2d Point2;

typedef Point3_<int> Point3i;
typedef Point3_<long long> Point3ll;
typedef Point3_<float> Point3f;
typedef Point3_<double> Point3d;
typedef Point3d Point3;

typedef Point4_<int> Point4i;
typedef Point4_<long long> Point4ll;
typedef Point4_<float> Point4f;
typedef Point4_<double> Point4d;
typedef Point4d Point4;

typedef std::vector<Point2_<int> > Points2i;
typedef std::vector<Point2_<long long> > Points2ll;
typedef std::vector<Point2_<float> > Points2f;
typedef std::vector<Point2_<double> > Points2d;
typedef std::vector<Point2> Points2;

typedef std::vector<Point3_<int> > Points3i;
typedef std::vector<Point3_<long long> > Points3ll;
typedef std::vector<Point3_<float> > Points3f;
typedef std::vector<Point3_<double> > Points3d;
typedef std::vector<Point3> Points3;

typedef ublas::matrix<double,ublas::row_major,
                      ublas::bounded_array<double, 4> > Matrix2;
typedef ublas::matrix<double,ublas::row_major,
                      ublas::bounded_array<double, 9> > Matrix3;
typedef ublas::matrix<double,ublas::row_major,
                      ublas::bounded_array<double, 16> > Matrix4;


// handy point->size conversion
template<typename T>
inline Size2_<T> size(const Point2_<T> &p) {
    return { p(0), p(1) };
}

// handy size->point conversion
template<typename T>
inline Point2_<T> point(const Size2_<T> &s) {
    return { s.width, s.height };
}

// handy point->size conversion for point3
template<typename T>
inline Size3_<T> size(const Point3_<T> &p) {
    return { p(0), p(1), p(2) };
}

// handy size->point conversion
template<typename T>
inline Point3_<T> point(const Size3_<T> &s) {
    return { s.width, s.height, s.depth };
}

/** Helper tag struct for invalid value extents instationation.
 */
struct InvalidExtents {};

/** Helper tag struct for inclusive size computation:
 *
 *  size(Extents<integral>(x, y, x, y)) is (0, 0)
 *  size(Extents<integral>(x, y, x, y), Inclusive{}) is (1, 1)
 */
struct Inclusive {};

template <typename T>
struct Extents2_ {
    typedef T value_type;
    typedef Point2_<T> point_type;

    point_type ll;
    point_type ur;

    Extents2_() : ll(), ur() {}

    explicit Extents2_(InvalidExtents)
        : ll(std::numeric_limits<T>::max()
             , std::numeric_limits<T>::max())
        , ur(std::numeric_limits<T>::lowest()
             , std::numeric_limits<T>::lowest())
    {}

    explicit Extents2_(const point_type &p)
        : ll(p), ur(p)
    {}

    Extents2_(const point_type &ll, const point_type &ur)
        : ll(ll), ur(ur)
    {}

    Extents2_(const value_type &xll, const value_type &yll
              , const value_type &xur, const value_type &yur)
        : ll(xll, yll), ur(xur, yur)
    {}

    template <typename U>
    explicit Extents2_(const Extents2_<U> &e)
        : ll(e.ll), ur(e.ur)
    {}

    T area() const {
        if ( ur[1] < ll[1] || ur[0] < ll[0] ) return 0;
        return ( ur[1] - ll[1] ) * ( ur[0] - ll[0] ); }

    T size() const {
        return std::max(ur[0] - ll[0], ur[1] - ll[1]);
    }
};

typedef Extents2_<int> Extents2i;
typedef Extents2_<long long> Extents2ll;
typedef Extents2_<double> Extents2f;
typedef Extents2f Extents2;

template <typename T1, typename T2>
inline bool inside(const Extents2_<T1> &e, const T2 &p0, const T2 &p1)
{
    return ((p0 >= e.ll(0)) && (p0 <= e.ur(0))
            && (p1 >= e.ll(1)) && (p1 <= e.ur(1)));
}

template <typename T1, typename T2>
inline bool inside(const Extents2_<T1> &e, const Point2_<T2> &p)
{
    return ((p(0) >= e.ll(0)) && (p(0) <= e.ur(0))
            && (p(1) >= e.ll(1)) && (p(1) <= e.ur(1)));
}

/** Checks only for x and y components of Point3
 */
template <typename T1, typename T2>
inline bool inside(const Extents2_<T1> &e, const Point3_<T2> &p)
{
    return ((p(0) >= e.ll(0)) && (p(0) <= e.ur(0))
            && (p(1) >= e.ll(1)) && (p(1) <= e.ur(1)));
}

template <typename T1, typename T2>
inline Extents2_<T1> operator+(const Extents2_<T1> &e, const T2 &diff)
{
    return { T1(e.ll(0) - diff), T1(e.ll(1) - diff)
            , T1(e.ur(0) + diff), T1(e.ur(1) + diff) };
}

template <typename T>
Viewport2_<T> viewport(const Extents2_<T> &e)
{
    return Viewport2_<T>(e.ur(0) - e.ll(0), e.ur(1) - e.ll(1)
                         , e.ll(0), e.ll(1));
}

template <typename T>
Extents2_<T> extents(const Viewport2_<T> &v)
{
    return Extents2_<T>(v.x, v.y, v.x + v.width, v.y + v.height);
}

template <typename T>
const typename Extents2_<T>::point_type& ll(const Extents2_<T> &e) {
    return e.ll;
}

template <typename T>
const typename Extents2_<T>::point_type& ur(const Extents2_<T> &e) {
    return e.ur;
}

template <typename T>
typename Extents2_<T>::point_type ul(const Extents2_<T> &e) {
    return typename Extents2_<T>::point_type(e.ll(0), e.ur(1));
}

template <typename T>
typename Extents2_<T>::point_type lr(const Extents2_<T> &e) {
    return typename Extents2_<T>::point_type(e.ur(0), e.ll(1));
}

template<typename T>
inline Point2_<T> clip(const Extents2_<T> &e, const Point2_<T> &p) {
    return Point2_<T>(std::max(e.ll(0), std::min(e.ur(0), p(0)))
                      , std::max(e.ll(1), std::min(e.ur(1), p(1))));
}

template<typename T>
inline Point2_<T> center(const Extents2_<T> &e) {
    return (e.ur + e.ll) / 2;
}

template<typename T>
inline Size2_<T> size(const Extents2_<T> &e) {
    return Size2_<T>(e.ur(0) - e.ll(0), e.ur(1) - e.ll(1));
}

template<typename T>
inline Size2_<T> size(const Extents2_<T> &e, const Inclusive&) {
    return Size2_<T>(e.ur(0) - e.ll(0) + T(1), e.ur(1) - e.ll(1) + T(1));
}

template<typename T>
inline T area(const Extents2_<T> &e) {
    return ((e.ur[1] < e.ll[1]) || (e.ur[0] < e.ll[0]))
        ? 0
        : (e.ur[1] - e.ll[1]) * (e.ur[0] - e.ll[0]);
}

template<typename T>
inline bool empty(const Extents2_<T> &e) {
    return !area(e);
}

template<typename T>
inline bool valid(const Extents2_<T> &e) {
    return (e.ll(0) <= e.ur(0)) && (e.ll(1) <= e.ur(1));
}

template <typename T>
inline Extents2_<T> unite( const Extents2_<T> &a, const Extents2_<T> &b ) {

    return Extents2_<T>(
        typename Extents2_<T>::point_type (
            std::min( a.ll[0], b.ll[0] ),
            std::min( a.ll[1], b.ll[1] ) ),
        typename Extents2_<T>::point_type (
            std::max( a.ur[0], b.ur[0] ),
            std::max( a.ur[1], b.ur[1] ) ) );
}

template <typename T1, typename T2>
inline bool overlaps(const Extents2_<T1> &a, const Extents2_<T2> &b)
{
    return ((a.ll(0) < b.ur(0)) && (b.ll(0) < a.ur(0))
             && (a.ll(1) < b.ur(1)) && (b.ll(1) < a.ur(1)));
}

template <typename T>
inline Extents2_<T> intersect( const Extents2_<T> &a, const Extents2_<T> &b ) {
    if (!overlaps(a, b)) {
        throw NoIntersectError
            ("Extents do not overlap, cannot compute intersection");
    }

    return Extents2_<T>(
        typename Extents2_<T>::point_type (
            std::max( a.ll[0], b.ll[0] ),
            std::max( a.ll[1], b.ll[1] ) ),
        typename Extents2_<T>::point_type (
            std::min( a.ur[0], b.ur[0] ),
            std::min( a.ur[1], b.ur[1] ) ) );
}


template <typename T>
inline double overlap( const Extents2_<T> & a, const Extents2_<T> & b ) {
    try {
        return ( intersect( a, b ).area() / unite( a, b ).area() );
    } catch (const NoIntersectError &) {
        return .0;
    }
}

template <typename T>
inline void update(Extents2_<T> &e, const T &x, const T &y) {
    e.ll(0) = std::min(e.ll(0), x);
    e.ll(1) = std::min(e.ll(1), y);

    e.ur(0) = std::max(e.ur(0), x);
    e.ur(1) = std::max(e.ur(1), y);
}

template <typename T>
inline void update(Extents2_<T> &e, const Point2_<T> &p) {
    e.ll(0) = std::min(e.ll(0), p(0));
    e.ll(1) = std::min(e.ll(1), p(1));

    e.ur(0) = std::max(e.ur(0), p(0));
    e.ur(1) = std::max(e.ur(1), p(1));
}

template <typename T>
inline void update(Extents2_<T> &e, const Point3_<T> &p) {
    e.ll(0) = std::min(e.ll(0), p(0));
    e.ll(1) = std::min(e.ll(1), p(1));

    e.ur(0) = std::max(e.ur(0), p(0));
    e.ur(1) = std::max(e.ur(1), p(1));
}

template <typename T>
inline void update(Extents2_<T> &e, const Extents2_<T> &other) {
    update(e, other.ll);
    update(e, other.ur);
}

template <typename P, typename R, typename T>
inline P snapToGrid(const P & point, const P & origin, T step, R roundFcn) {
    P res;

    for (std::size_t i(0); i < res.size(); ++i) {
        res(i) = roundFcn((point(i) - origin(i))/step) * step + origin(i);
    }

    return res;
}

template <typename E, typename T>
inline E snapToGrid( const E &ext, const typename E::point_type &origin
                   , T step, bool inscribe = false)
{
    typedef typename E::value_type VT;
    if (inscribe) {
        return { snapToGrid(ext.ll, origin, step
                            , [](VT value) { return std::ceil(value); })
                , snapToGrid(ext.ur, origin, step
                             , [](VT value) { return std::floor(value); }) };
    }

    return { snapToGrid(ext.ll, origin, step
                        , [](VT value) { return std::floor(value); })
            , snapToGrid(ext.ur, origin, step
                         , [](VT value) { return std::ceil(value); }) };
}

template <typename T>
struct Extents3_ {
    typedef T value_type;
    typedef Point3_<T> point_type;

    point_type ll;
    point_type ur;

    Extents3_() : ll(), ur() {}

    explicit Extents3_(InvalidExtents)
        : ll(std::numeric_limits<T>::max(), std::numeric_limits<T>::max()
             , std::numeric_limits<T>::max())
        , ur(std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest()
             , std::numeric_limits<T>::lowest())
    {}

    explicit Extents3_(const point_type &p)
        : ll(p), ur(p)
    {}

    Extents3_(const point_type &ll, const point_type &ur)
        : ll(ll), ur(ur)
    {}

    Extents3_(const value_type &xll, const value_type &yll
              , const value_type &zll
              , const value_type &xur, const value_type yur
              , const value_type &zur)
        : ll(xll, yll, zll), ur(xur, yur, zur)
    {}

    template <typename U>
    explicit Extents3_(const Extents3_<U> &e)
        : ll(e.ll), ur(e.ur)
    {}
};

typedef Extents3_<int> Extents3i;
typedef Extents3_<long long> Extents3ll;
typedef Extents3_<double> Extents3f;
typedef Extents3f Extents3;

template<typename T>
inline Point3_<T> center(const Extents3_<T> &e) {
    return (e.ur + e.ll) / 2;
}

template<typename T>
inline Point3_<T> clip(const Extents3_<T> &e, const Point3_<T> &p) {
    return Point3_<T>(std::max(e.ll(0), std::min(e.ur(0), p(0)))
                      , std::max(e.ll(1), std::min(e.ur(1), p(1)))
                      , std::max(e.ll(2), std::min(e.ur(2), p(2))));
}

template<typename T>
inline Size3_<T> size(const Extents3_<T> &e) {
    return Size3_<T>(e.ur(0) - e.ll(0), e.ur(1) - e.ll(1), e.ur(2) - e.ll(2));
}

template<typename T>
inline Size3_<T> size(const Extents3_<T> &e, const Inclusive&) {
    return Size3_<T>(e.ur(0) - e.ll(0) + T(1)
                     , e.ur(1) - e.ll(1) + T(1)
                     , e.ur(2) - e.ll(2) + T(1));
}

template<typename T>
inline T volume(const Extents3_<T> &e) {
    return ((e.ur[2] < e.ll[2]) || (e.ur[1] < e.ll[1]) || (e.ur[0] < e.ll[0]))
        ? 0
        : (e.ur[2] - e.ll[2]) * (e.ur[1] - e.ll[1]) * (e.ur[0] - e.ll[0]);
}

template<typename T>
inline bool empty(const Extents3_<T> &e) {
    return !volume(e);
}

template<typename T>
inline bool valid(const Extents3_<T> &e) {
    return ((e.ll(0) <= e.ur(0))
            && (e.ll(1) <= e.ur(1))
            && (e.ll(2) <= e.ur(2)));
}

template <typename T>
inline void update(Extents3_<T> &e, const Point3_<T> &p) {
    e.ll(0) = std::min(e.ll(0), p(0));
    e.ll(1) = std::min(e.ll(1), p(1));
    e.ll(2) = std::min(e.ll(2), p(2));

    e.ur(0) = std::max(e.ur(0), p(0));
    e.ur(1) = std::max(e.ur(1), p(1));
    e.ur(2) = std::max(e.ur(2), p(2));
}


template <typename T>
inline void update(Extents3_<T> &e, const Extents3_<T> &other) {
    update(e, other.ll);
    update(e, other.ur);
}

template <typename T>
inline Extents3_<T> unite(const Extents3_<T> &a, const Extents3_<T> &b ) {
    Extents3_<T> res(a);
    update(res, b.ll);
    update(res, b.ur);
    return res;
}

template <typename T1, typename T3>
inline bool overlaps(const Extents3_<T1> &a, const Extents3_<T3> &b)
{
    return ((a.ll(0) < b.ur(0)) && (b.ll(0) < a.ur(0))
             && (a.ll(1) < b.ur(1)) && (b.ll(1) < a.ur(1))
             && (a.ll(2) < b.ur(2)) && (b.ll(2) < a.ur(2)));
}

template <typename T>
inline Extents3_<T> intersect( const Extents3_<T> &a, const Extents3_<T> &b ) {
    if (!overlaps(a, b)) {
        throw NoIntersectError
            ("Extents do not overlap, cannot compute intersection");
    }

    typedef Extents3_<T> E3;
    typedef typename E3::point_type point_type;
    return Extents3(point_type(std::max(a.ll[0], b.ll[0])
                               , std::max(a.ll[1], b.ll[1])
                               , std::max(a.ll[2], b.ll[2]))
                    , point_type(std::min(a.ll[0], b.ll[0])
                                 , std::min(a.ll[1], b.ll[1])
                                 , std::min(a.ll[2], b.ll[2])));
}

template <typename T>
inline double overlap(const Extents3_<T> & a, const Extents3_<T> & b) {
    try {
        return (volume(intersect(a, b)) / volume(unite(a, b)));
    } catch (const NoIntersectError&) {
        return .0;
    }
}

template<typename CharT, typename Traits, typename T>
inline std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits> &os, const Viewport2_<T> &v)
{
    std::ios::fmtflags flags(os.flags());
    os << v.width << "x" << v.height << std::showpos << v.x << v.y;
    os.flags(flags);
    return os;
}

template<typename CharT, typename Traits, typename T>
inline std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits> &is, Viewport2_<T> &v)
{
    using boost::spirit::qi::auto_;
    using boost::spirit::qi::char_;
    using boost::spirit::qi::omit;
    using boost::spirit::qi::match;

    char sign1, sign2;

    is >> match((auto_ >> omit['x'] >> auto_
                 >> (char_('+') | char_('-')) >> auto_
                 >> (char_('+') | char_('-')) >> auto_)
                , v.width, v.height, sign1, v.x, sign2, v.y);

    if (sign1 == '-') { v.x = -v.x; }
    if (sign2 == '-') { v.y = -v.y; }
    return is;
}

#if 0
// original unsafe code
template <class UblasContainer>
inline bool operator < (
    const UblasContainer & op1, const UblasContainer & op2 ) {

    return std::lexicographical_compare(
                op1.data().begin(), op1.data().end(),
                op2.data().begin(), op2.data().end() );
}

#else

inline bool operator < (
    const Matrix4 & op1, const Matrix4 & op2 ) {
    return std::lexicographical_compare(
                op1.data().begin(), op1.data().end(),
                op2.data().begin(), op2.data().end() );
}

template <typename T1, typename T2>
inline bool operator < (
    const Point3_<T1> & op1, const Point3_<T2> & op2 ) {

    return std::lexicographical_compare(
                op1.data().begin(), op1.data().end(),
                op2.data().begin(), op2.data().end() );
}

template <typename T1, typename T2>
inline bool operator < (
    const Point2_<T1> & op1, const Point2_<T2> & op2 ) {

    return std::lexicographical_compare(
                op1.data().begin(), op1.data().end(),
                op2.data().begin(), op2.data().end() );
}

#endif

template<typename CharT, typename Traits, typename T>
inline std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits> &os, const Extents2_<T> &e)
{
    os << e.ll(0) << ',' << e.ll(1) << ':' << e.ur(0) << ',' << e.ur(1);
    return os;
}

template<typename CharT, typename Traits, typename T>
inline std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits> &is, Extents2_<T> &e)
{
    using boost::spirit::qi::auto_;
    using boost::spirit::qi::char_;
    using boost::spirit::qi::omit;
    using boost::spirit::qi::match;

    is >> match((auto_ >> omit[','] >> auto_ >> omit[':']
                 >> auto_ >> omit[','] >> auto_)
                , e.ll(0), e.ll(1), e.ur(0), e.ur(1));

    return is;
}

template<typename CharT, typename Traits, typename T>
inline std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits> &os, const Extents3_<T> &e)
{
    os << e.ll(0) << ',' << e.ll(1) << ',' << e.ll(2) << ':'
       << e.ur(0) << ',' << e.ur(1) << ',' << e.ur(2);
    return os;
}

template<typename CharT, typename Traits, typename T>
inline std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits> &is, Extents3_<T> &e)
{
    using boost::spirit::qi::auto_;
    using boost::spirit::qi::char_;
    using boost::spirit::qi::omit;
    using boost::spirit::qi::match;

    is >> match((auto_ >> omit[','] >> auto_ >> omit[','] >> auto_ >> omit[':']
                 >> auto_ >> omit[','] >> auto_ >> omit[','] >> auto_)
                , e.ll(0), e.ll(1), e.ll(2), e.ur(0), e.ur(1), e.ur(2));

    return is;
}

template <typename T1, typename T2>
inline bool inside(const Extents3_<T1> &e, const T2 &p0, const T2 &p1
                   , const T2 &p2)
{
    return ((p0 >= e.ll(0)) && (p0 <= e.ur(0))
            && (p1 >= e.ll(1)) && (p1 <= e.ur(1))
            && (p2 >= e.ll(2)) && (p2 <= e.ur(2)));
}

template <typename T1, typename T2>
inline bool inside(const Extents3_<T1> &e, const Point3_<T2> &p)
{
    return ((p(0) >= e.ll(0)) && (p(0) <= e.ur(0))
            && (p(1) >= e.ll(1)) && (p(1) <= e.ur(1))
            && (p(2) >= e.ll(2)) && (p(2) <= e.ur(2)));
}

template <typename T1, typename T2>
inline Extents3_<T1> operator+(const Extents3_<T1> &e, const T2 &diff)
{
    return { e.ll(0) - diff, e.ll(1) - diff, e.ll(2) - diff
            , e.ur(0) + diff, e.ur(1) + diff, e.ur(2) + diff };
}

template <typename T>
const typename Extents3_<T>::point_type& bll(const Extents3_<T> &e) {
    return e.ll;
}

template <typename T>
typename Extents3_<T>::point_type bul(const Extents3_<T> &e) {
    return typename Extents3_<T>::point_type(e.ll(0), e.ur(1), e.ll(2));
}

template <typename T>
typename Extents3_<T>::point_type bur(const Extents3_<T> &e) {
    return typename Extents3_<T>::point_type(e.ur(0), e.ur(1), e.ll(2));
}

template <typename T>
typename Extents3_<T>::point_type blr(const Extents3_<T> &e) {
    return typename Extents3_<T>::point_type(e.ur(0), e.ll(1), e.ll(2));
}

template <typename T>
const typename Extents3_<T>::point_type tll(const Extents3_<T> &e) {
    return typename Extents3_<T>::point_type(e.ll(0), e.ll(1), e.ur(2));
}

template <typename T>
typename Extents3_<T>::point_type tul(const Extents3_<T> &e) {
    return typename Extents3_<T>::point_type(e.ll(0), e.ur(1), e.ur(2));
}

template <typename T>
const typename Extents3_<T>::point_type& tur(const Extents3_<T> &e) {
    return e.ur;
}

template <typename T>
typename Extents3_<T>::point_type tlr(const Extents3_<T> &e) {
    return typename Extents3_<T>::point_type(e.ur(0), e.ll(1), e.ur(2));
}

/** Get all 4 vertices from Extents2 as a vector of points.
 */
template <typename T>
const std::vector<typename Extents2_<T>::point_type>
vertices(const Extents2_<T> &e)
{
    return { ll(e), ul(e), ur(e), lr(e) };
}

/** Get all 8 vertices from Extents3 as a vector of points.
 */
template <typename T>
const std::vector<typename Extents3_<T>::point_type>
vertices(const Extents3_<T> &e)
{
    return { bll(e), bul(e), bur(e), blr(e), tll(e), tul(e), tur(e), tlr(e) };
}

/** Helper functions to convert between 2d and 3d entities
 */

// Size2_ <-> Size3_

template <typename T>
Size2_<T> size2(const Size3_<T> &s)
{
    return Size2_<T>(s.width, s.height);
}

template <typename T>
const Size2_<T>& size2(const Size2_<T> &s) { return s; }

template <typename T>
Size3_<T> size3(const Size2_<T> &s)
{
    return Size3_<T>(s.width, s.height, T(0));
}

template <typename T>
const Size3_<T>& size3(const Size3_<T> &s) { return s; }

// Extents2_ <-> Extents3_

template <typename T>
Extents2_<T> extents2(const Extents3_<T> &e)
{
    return Extents2_<T>(e.ll(0), e.ll(1), e.ur(0), e.ur(1));
}

template <typename T>
const Extents2_<T>& extents2(const Extents2_<T> &e) { return e; }

template <typename T>
Extents3_<T> extents3(const Extents2_<T> &e)
{
    return Extents3_<T>(e.ll(0), e.ll(1), T(0), e.ur(0), e.ur(1), T(0));
}

template <typename T>
const Extents3_<T>& extents3(const Extents3_<T> &e) { return e; }

// Point2_ <-> Point3_

template <typename T>
Point2_<T> point2(const Point3_<T> &p)
{
    return Point2_<T>(p(0), p(1));
}

template <typename T>
const Point2_<T>& point2(const Point2_<T> &p) { return p; }

template <typename T>
Point3_<T> point3(const Point2_<T> &p)
{
    return Point3_<T>(p(0), p(1), T(0));
}

template <typename T>
const Point3_<T>& point3(const Point3_<T> &p) { return p; }


#define MATH_GEOMETRY_CORE_HPP_INLINES_
#include "detail/size2.inline.hpp"
#include "detail/size3.inline.hpp"
#include "detail/extents2.inline.hpp"
#undef MATH_GEOMETRY_CORE_HPP_INLINES_

} // namespace math

#endif // MATH_GEOMETRY_CORE_HPP

