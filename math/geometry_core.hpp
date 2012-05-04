/**
 * @file geometry_core.hpp
 * @author Ondrej Prochazka <ondrej.prochazka@citationtech.net>
 *
 * Basic structures for geometric modeling
 */

#ifndef MATH_GEOMETRY_CORE_HPP
#define MATH_GEOMETRY_CORE_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_match.hpp>
#include <boost/spirit/include/qi_match_auto.hpp>
#include <boost/spirit/include/qi_alternative.hpp>

#include <vector>
#include <iostream>
#include <iomanip>

namespace ublas = boost::numeric::ublas;

namespace math {

/* sizes */

template <class T>
struct Size2_ {
    Size2_() : width( 0 ), height( 0 ) {};

    Size2_( const T & width, const T & height )
        : width( width ), height( height ) {};

    T width, height;
};


typedef Size2_<int> Size2i;
typedef Size2_<double> Size2f;
typedef Size2i Size2;

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

    Viewport2_(value_type width, value_type height)
        : width(width), height(height), x(0), y(0) {};

    Viewport2_(const size_type &size)
        : width( size.width ), height( size.height ), x( 0 ), y( 0 ) {};

    size_type size() const { return math::Size2i(width, height); }
};

typedef Viewport2_<int> Viewport2i;
typedef Viewport2_<double> Viewport2f;
typedef Viewport2i Viewport2;

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

    bool operator== (const Point3_<T>& p) const {
        return (*this)(0) == p(0) && (*this)(1) == p(1) && (*this)(2) == p(2);
    }
    bool operator!= (const Point3_<T>& p) const {
        return !operator==(p);
    }
};

typedef Point2_<int> Point2i;
typedef Point2_<float> Point2f;
typedef Point2_<double> Point2d;
typedef Point2d Point2;

typedef Point3_<int> Point3i;
typedef Point3_<float> Point3f;
typedef Point3_<double> Point3d;
typedef Point3d Point3;

typedef std::vector<Point2_<int> > Points2i;
typedef std::vector<Point2_<float> > Points2f;
typedef std::vector<Point2_<double> > Points2d;
typedef std::vector<Point2> Points2;

typedef std::vector<Point3_<int> > Points3i;
typedef std::vector<Point3_<float> > Points3f;
typedef std::vector<Point3_<double> > Points3d;
typedef std::vector<Point3> Points3;

typedef ublas::matrix<double,ublas::row_major,
                      ublas::bounded_array<double, 16> > Matrix4;


template<typename CharT, typename Traits, typename T>
inline std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits> &os, const Size2_<T> &s)
{
    return os << s.width << "x" << s.height;
}

template<typename CharT, typename Traits, typename T>
inline std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits> &is, math::Size2_<T> &s)
{
    using boost::spirit::qi::auto_;
    using boost::spirit::qi::char_;
    using boost::spirit::qi::omit;
    using boost::spirit::qi::match;

    return is >> match(auto_ >> omit['x'] >> auto_, s.width, s.height);
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

} // namespace math

#endif
