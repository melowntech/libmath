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


Point3 intersectionParams(const Line3 &line, const legacy::Plane3 &plane)
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

Point3 intersection( const Line3 & line, const legacy::Plane3 & plane )
{
    return line.point(intersectionParams(line, plane)(0));
}

double pointPlaneDistance(const Point3& p, const Plane3& plane)
{
    return std::abs(ublas::inner_prod(plane.n_, p) + plane.d_)
           / ublas::norm_2(plane.n_);
}

Point3 pointPlaneProjection(const Point3& p, const Plane3& plane)
{
    double k = (ublas::inner_prod(plane.n_, p) + plane.d_)
               / ublas::inner_prod(plane.n_, plane.n_);

    return Point3(p(0) - k * plane.n_(0),
                  p(1) - k * plane.n_(1),
                  p(2) - k * plane.n_(2));
}

Point3 linePlaneIntersection(const Line3& l, const Plane3& plane)
{
    double t = (-ublas::inner_prod(plane.n_, l.p) - plane.d_)
               / (ublas::inner_prod(plane.n_, l.u));
    return l.p + t * l.u;
}

Line3 planeIntersection(const Plane3& p1, const Plane3& p2)
{
    // normalize plane representation
    double norm1 = ublas::norm_2(p1.n_);
    Point3 n1 = p1.n_ / norm1;
    double d1 = p1.d_ / norm1;

    double norm2 = ublas::norm_2(p2.n_);
    Point3 n2 = p2.n_ / norm2;
    double d2 = p2.d_ / norm2;

    // get one point on line
    double n1n2 = ublas::inner_prod(n1, n2);   
    double den = 1 - std::pow(n1n2, 2);
    double c1 = ((d2 * n1n2) - d1) / den;
    double c2 = ((d1 * n1n2) - d2) / den;
    Point3 pt = c1 * n1 + c2 * n2;

    return Line3(pt, crossProduct(n1, n2));
}

Point3 planeIntersection(const Plane3& p1, const Plane3& p2, const Plane3& p3)
{
    Matrix3 D = ublas::zero_matrix<double>(3, 3);
    ublas::row(D, 0) = p1.n_;
    ublas::row(D, 1) = p2.n_;
    ublas::row(D, 2) = p3.n_;

    Point3 p(-p1.d_, -p2.d_, -p3.d_);

    Matrix3 Dx = D;
    ublas::column(Dx, 0) = p;
    Matrix3 Dy = D;
    ublas::column(Dy, 1) = p;
    Matrix3 Dz = D;
    ublas::column(Dz, 2) = p;

    double den = determinant(D);

    return Point3(determinant(Dx) / den,
                  determinant(Dy) / den,
                  determinant(Dz) / den);
}

math::Matrix4 createPlaneCrs(const math::Plane3& plane)
{
    // origin of the crs
    math::Point3 og = math::pointPlaneProjection(math::Point3(0, 0, 0), plane);

    // get some second point on the plane
    math::Point3 p1 = math::pointPlaneProjection(math::Point3(1, 0, 0), plane);
    math::Point3 p2 = math::pointPlaneProjection(math::Point3(0, 1, 0), plane);

    // choose the more distant to avoid singularity
    math::Point3 p = math::length(og - p1) > math::length(og - p2) ? p1 : p2;

    // two base vectors
    math::Point3 n1 = math::normalize(p - og);
    math::Point3 n3 = math::normalize(plane.n_);

    // get the third to form a right-hand crs
    math::Point3 n2 = math::normalize(math::crossProduct(n3, n1));

    math::Matrix4 tf = boost::numeric::ublas::identity_matrix<double>(4, 4);

    auto col1 = ublas::column(tf, 0);
    auto col2 = ublas::column(tf, 1);
    auto col3 = ublas::column(tf, 2);
    auto col4 = ublas::column(tf, 3);
    ublas::subrange(col1, 0, 3) = n1;
    ublas::subrange(col2, 0, 3) = n2;
    ublas::subrange(col3, 0, 3) = n3;
    ublas::subrange(col4, 0, 3) = og;
    return tf;
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
    return float(point(0)*normal(0) + point(1)*normal(1) + d);
}

bool triangleRectangleCollision( math::Point2 triangle[3]
                               , math::Point2 ll, math::Point2 ur){
    //find collision using SAT, triangles are presumed to have points in CCW order
    //first try all half spaces of the rectangle
    //on the left side of rectangle
    bool allout = true;
    for(unsigned int i=0;i<3;++i){
        if(triangle[i](0)> ll(0)){
            allout = false;
        }
    }
    if(allout){
        return false;
    }
    //on the right side of rectangle
    allout = true;
    for(unsigned int i=0;i<3;++i){
        if(triangle[i](0)< ur(0)){
            allout = false;
        }
    }
    if(allout){
        return false;
    }
    //on the bottom side of rectangle
    allout = true;
    for(unsigned int i=0;i<3;++i){
        if(triangle[i](1)> ll(1)){
            allout = false;
        }
    }
    if(allout){
        return false;
    }
    //on the top side of rectangle
    allout = true;
    for(unsigned int i=0;i<3;++i){
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

    for(unsigned int i=0;i<3;++i){
        bool allout = true;
        math::Point2 lvec = triangle[i]-triangle[(i+1)%3];
        math::Point2 lnormal = normalize(math::Point2(lvec(1), -lvec(0))); 
        double ld = -lnormal(0)*triangle[i](0)-lnormal(1)*triangle[i](1);
        for(unsigned int c=0;c<4;++c){
            allout = allout && signedDistance(corners[c], lnormal, ld) < 0; 
        }
        if(allout){
            return false;
        }
    }
    return true;
}


math::Point3 barycentricCoords(const math::Point2 &r, const math::Point2 &a
                               , const math::Point2 &b, const math::Point2 &c)
{
    math::Matrix2 A(2, 2);
    math::Point2 l12, rhs;

    for (int i = 0; i < 2; i++) {
        A(i,0) = a(i) - c(i);
        A(i,1) = b(i) - c(i);
        rhs(i) = r(i) - c(i);
    }

    math::solve2x2(A, rhs, l12);

    return {l12(0), l12(1), 1.0 - l12(0) - l12(1)};
}


math::Point3 centroid(const math::Points3& points,
                      const std::vector<double>& weights,
                      const double eps)
{
    // calculate the (weighted) centroid
    bool weigthed = !weights.empty();
    math::Point3 centroid(0, 0, 0);
    double wSum = 0.0;
    for (std::size_t i = 0; i < points.size(); i++) {
        math::Point3 pt(points[i]);
        if (weigthed) {
            pt *= weights[i];
            wSum += weights[i];
        }
        centroid += pt;
    }
    if (weigthed && (std::abs(wSum) < eps)) {
        LOGTHROW(err4, std::runtime_error) << "Point weights are too small.";
    }

    if (weigthed) { centroid /= wSum; }
    else { centroid /= points.size(); }

    return centroid;
}

math::Point3 centroid(const math::Points3& points)
{
    return centroid(points, {}, 0);
}

} // namespace math

