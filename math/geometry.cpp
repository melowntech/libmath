/*
 * geometry.cpp
 */
 
#include "geometry.hpp"

#include "math.hpp"

#include <dbglog/dbglog.hpp>

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


Point3 midpoint( const Line3 & line1, const Line3 & line2 ) {

    Matrix2 a(2, 2);
    Point2 b(2);
    double r1, r2;
    
    a(0,0) = - inner_prod( line1.u, line2.u );
    a(0,1) = inner_prod( line2.u, line2.u );
    a(1,0) = inner_prod( line1.u, line1.u );
    a(1,1) = - inner_prod( line1.u, line2.u );
    b(0) = inner_prod( line2.u, line1.p - line2.p );
    b(1) = inner_prod( line1.u, line2.p - line1.p );

    if ( fabs( inner_prod( line1.u, line2.u ) / ( norm_2( line1.u ) * norm_2( line2.u ) ) ) > 0.9962 )
      LOGTHROW( err1, std::runtime_error )
          << "Lines close to paralel, midpoint not reliable";

    // test for parallel vectors
    double detA = a(0,0) * a(1,1) - a(1,0) * a(0,1);
          
    if ( fabs( detA ) < 1E-15 ) {

        LOGTHROW( err1, std::runtime_error )
            << "Paralel lines detected, no midpoint.";        
    }


    /*
    ublas::matrix<float> ai = matrixInvert( a );
    ublas::vector<float> r = ublas::prod( ai, b );*/

    r1 = ( b(0) * a(1,1) - b(1) * a(0,1) ) / detA;
    r2 = ( a(0,0) * b(1) - a(1,0) * b(0) ) / detA;

    return ( line1.p + r1 * line1.u + line2.p + r2 * line2.u ) * 0.5;
}


ublas::vector<double> intersection(
    const Line3 & line, const Plane3 & plane ) {

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


double polygonRegularity(
    const Point3 & v0, const Point3 & v1, const Point3 & v2  ) {

    double areaToCircumverence =
        ublas::norm_2( crossProduct( v1 - v0, v2 - v1 ) * 0.5 ) /
        sqr( ublas::norm_2( v1 - v0 ) + ublas::norm_2( v2 - v1 ) +
          ublas::norm_2( v0 - v2 ) ) ;

    //LOG( debug ) << areaToCircumverence / ( sqrt( 3 ) / 36.0 );
        
    return areaToCircumverence / ( sqrt( 3 ) / 36.0 );
}

double polygonRegularity(
    const Point3 & v0, const Point3 & v1, const Point3 & v2,
    const Point3 & v3  ) {

    double areaToCircumverence =
        ( ublas::norm_2( crossProduct( v1 - v0, v2 - v1 ) * 0.5 )
          + ublas::norm_2( crossProduct( v3 - v2, v0 - v3 ) * 0.5 ) ) /
        sqr( ublas::norm_2( v1 - v0 ) + ublas::norm_2( v2 - v1 ) +
          ublas::norm_2( v3 - v2 ) + ublas::norm_2 ( v0 - v3 ) );

    //LOG( debug ) << areaToCircumverence / ( 1 / 16.0 );
        
    return areaToCircumverence / ( 1 / 16.0 );
}


} // namespace math

