/*
 * geometry.cpp
 */
 
#include "geometry.hpp"

#include "math.hpp"

namespace math {

float lineDistance(
    const ublas::vector<float> & p1, const ublas::vector<float> & u1,
    const ublas::vector<float> & p2, const ublas::vector<float> & u2,
    float & r1, float & r2 ) {

    ublas::matrix<float> a(2,2);
    ublas::vector<float> b(2);

    a(0,0) = - inner_prod( u1, u2 );
    a(0,1) = inner_prod( u2, u2 );
    a(1,0) = inner_prod( u1, u1 );
    a(1,1) = - inner_prod( u1, u2 );
    b(0) = inner_prod( u2, p1 - p2 );
    b(1) = inner_prod( u1, p2 - p1 );

    // test for parallel vectors
    if ( fabs( a(0,0) * a(1,1) - a(1,0) * a(0,1) ) < 1E-15 ) {

        // well, not quite right - but this is a singular case
        r1 = r2 = 0.0;
        return ublas::norm_2( p1 - p2 );
    }


    ublas::matrix<float> ai = matrixInvert( a );
    ublas::vector<float> r = ublas::prod( ai, b );
    r1 = r[0]; r2 = r[1];


    return ublas::norm_2( p1 + r1 * u1 - ( p2 + r2 * u2 ) );
}


ublas::vector<double> intersection(
    const Line3_t & line, const Plane3_t & plane ) {

    ublas::matrix<double> a( 3, 3 ), ai( 3, 3 );
    ublas::vector<double> op( 3 ), pars( 3 );

    ublas::column( a, 0 ) =  line.u;
    ublas::column( a, 1 ) = plane.u;
    ublas::column( a, 2 ) = plane.v;
    op = plane.p - line.p;

    ai = matrixInvert( a );
    pars = ublas::prod( ai, op );
    return line.p + line.u * pars( 0 );
}


} // namespace math

