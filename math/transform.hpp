/**
 * @file transform.hpp
 * @author Jakub Cerveny <jakub.cerveny@ext.citationtech.net>
 * @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 * Point transformations
 */

#ifndef MATH_TRANSFORM_HPP
#define MATH_TRANSFORM_HPP

#include "geometry_core.hpp"

namespace math {

inline Point3 transform(const Matrix4& tr, const Point3 &pt)
{
    return Point3(
        tr(0,0)*pt(0) + tr(0,1)*pt(1) + tr(0,2)*pt(2) + tr(0,3),
        tr(1,0)*pt(0) + tr(1,1)*pt(1) + tr(1,2)*pt(2) + tr(1,3),
        tr(2,0)*pt(0) + tr(2,1)*pt(1) + tr(2,2)*pt(2) + tr(2,3));
}

inline Point2 transform(const Matrix4& tr, const Point2 &pt)
{
    return Point2(
        tr(0,0)*pt(0) + tr(0,1)*pt(1) + tr(0,3),
        tr(1,0)*pt(0) + tr(1,1)*pt(1) + tr(1,3));
}

inline void transform(const Matrix4& tr, Points3& points)
{
    for (Point3& pt : points)
        pt = transform(tr, pt);
}

inline void transform(const Matrix4& tr, Points2& points)
{
    for (Point2& pt : points)
        pt = transform(tr, pt);
}

inline void transform(const Matrix4& tr, Extents2& extents)
{
    transform(tr, extents.ll);
    transform(tr, extents.ur);
}

} // namespace math

#endif // MATH_TRANSFORM_HPP
