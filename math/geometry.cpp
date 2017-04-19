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


Point3 midpoint( const Line3 & line1, const Line3 & line2, double minAngleCos ) 
{
    Matrix2 a(2, 2);
    Point2 b(2);
    double r1, r2;
    
    a(0,0) = - inner_prod( line1.u, line2.u );
    a(0,1) = inner_prod( line2.u, line2.u );
    a(1,0) = inner_prod( line1.u, line1.u );
    a(1,1) = - inner_prod( line1.u, line2.u );
    b(0) = inner_prod( line2.u, line1.p - line2.p );
    b(1) = inner_prod( line1.u, line2.p - line1.p );

    if ( fabs( inner_prod( line1.u, line2.u ) / ( norm_2( line1.u ) * norm_2( line2.u ) ) ) > minAngleCos )
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


double pointLineDistance(const Point3 &p, const Line3 &line)
{
    return norm_2(crossProduct(line.u, line.p - p)) / norm_2(line.u);
}


Point3 intersectionParams(const Line3 &line, const Plane3 &plane)
{
    ublas::matrix<double> a( 3, 3 ), ai( 3, 3 );
    ublas::vector<double> op( 3 ), pars( 3 );

    ublas::column( a, 0 ) =  line.u;
    ublas::column( a, 1 ) = plane.u;
    ublas::column( a, 2 ) = plane.v;
    op = plane.p - line.p;

    ai = matrixInvert( a );
    pars = ublas::prod( ai, op );

    return pars;
}

Point3 intersection( const Line3 & line, const Plane3 & plane )
{
    return line.point(intersectionParams(line, plane)(0));
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

inline float signedDistance(math::Point2 &point, math::Point2 &normal ,double d){
        return point(0)*normal(0) + point(1)*normal(1) + d;
}

bool triangleRectangleCollision( math::Point2 triangle[3]
                               , math::Point2 ll, math::Point2 ur){
    //find collision using SAT, triangles are presumed to have points in CCW order
    //first try all half spaces of the rectangle
    //on the left side of rectangle
    bool allout = true;
    for(uint i=0;i<3;++i){
        if(triangle[i](0)> ll(0)){
            allout = false;
        }
    }
    if(allout){
        return false;
    }
    //on the right side of rectangle
    allout = true;
    for(uint i=0;i<3;++i){
        if(triangle[i](0)< ur(0)){
            allout = false;
        }
    }
    if(allout){
        return false;
    }
    //on the bottom side of rectangle
    allout = true;
    for(uint i=0;i<3;++i){
        if(triangle[i](1)> ll(1)){
            allout = false;
        }
    }
    if(allout){
        return false;
    }
    //on the top side of rectangle
    allout = true;
    for(uint i=0;i<3;++i){
        if(triangle[i](1)< ur(1)){
            allout = false;
        }
    }
    if(allout){
        return false;
    }

    //now try all half spaces of the triangle
    math::Point2 corners[4];
    corners[0] = math::Point2(ll(0),ll(1));
    corners[1] = math::Point2(ll(0),ur(1));
    corners[2] = math::Point2(ur(0),ll(1));
    corners[3] = math::Point2(ur(0),ur(1));

    for(uint i=0;i<3;++i){
        bool allout = true;
        math::Point2 lvec = triangle[i]-triangle[(i+1)%3];
        math::Point2 lnormal = normalize(math::Point2(lvec(1), -lvec(0))); 
        double ld = -lnormal(0)*triangle[i](0)-lnormal(1)*triangle[i](1);
        for(uint c=0;c<4;++c){
            allout = allout && signedDistance(corners[c], lnormal, ld) < 0; 
        }
        if(allout){
            return false;
        }
    }
    return true;
}

} // namespace math

