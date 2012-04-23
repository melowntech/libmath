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

#include <vector>

namespace ublas = boost::numeric::ublas;

namespace math {

/* sizes */

template <class T>
struct Size2_ {

    Size2_( const T & width, const T & height )
        : width( width ), height( height ) {};
        
    T width, height; 
};

typedef Size2_<int> Size2i;
typedef Size2_<double> Size2f;
typedef Size2i Size2;

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
};

typedef Point2_<int> Point2i;
typedef Point2_<double> Point2f;
typedef Point2f Point2;

typedef Point3_<int> Point3i;
typedef Point3_<double> Point3f;
typedef Point3f Point3;

typedef std::vector<Point2_<int> > Points2i;
typedef std::vector<Point2_<double> > Points2f;
typedef std::vector<Point2> Points2;

typedef std::vector<Point3_<int> > Points3i;
typedef std::vector<Point3_<double> > Points3f;
typedef std::vector<Point3> Points3;


} // namespace math

#endif
