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

// matrix4

/** 4x4 identity matrix.
 */
inline Matrix4 identity4() {
    return boost::numeric::ublas::identity_matrix<double>(4, 4);
}

/** 4x4 zero matrix.
 */
inline Matrix4 zero4() {
    return boost::numeric::ublas::zero_matrix<double>(4, 4);
}

/** Makes 4x4 transformation matrix that shifts points by \param delta.
 */
template <typename T> Matrix4 translate(const math::Point3_<T> &delta);

/** Makes 4x4 transformation matrix that scales points by \param scale.
 */
template <typename T> Matrix4 scale(const math::Point3_<T> &scale);

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

inline Extents2 transform(const Matrix4& tr, const Extents2& extents)
{
    return Extents2(
        transform(tr, extents.ll),
        transform(tr, extents.ur));
}

inline Extents3 transform(const Matrix4& tr, const Extents3 &extents)
{
    return Extents3(transform(tr, extents.ll)
                    , transform(tr, extents.ur));
}

#ifdef MATH_HAS_OPENCV
template<typename T>
inline cv::Point3_<T> transform(const Matrix4& tr, const cv::Point3_<T>& pt)
{
    return cv::Point3_<T>(
        tr(0,0)*pt.x + tr(0,1)*pt.y + tr(0,2)*pt.z + tr(0,3),
        tr(1,0)*pt.x + tr(1,1)*pt.y + tr(1,2)*pt.z + tr(1,3),
        tr(2,0)*pt.x + tr(2,1)*pt.y + tr(2,2)*pt.z + tr(2,3));
}
#endif

// in-place transformations for vectors:

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

#ifdef MATH_HAS_OPENCV
template<typename T>
inline void transform(const Matrix4& tr, std::vector<cv::Point3_<T> >& points)
{
    for (cv::Point3_<T>& pt : points)
        pt = transform(tr, pt);
}
#endif


// matrix3

/** 3x3 identity matrix.
 */
inline Matrix3 identity3() {
    return boost::numeric::ublas::identity_matrix<double>(3, 3);
}

/** 3x3 zero matrix.
 */
inline Matrix3 zero3() {
    return boost::numeric::ublas::zero_matrix<double>(3, 3);
}

/** Makes 3x3 transformation matrix that shifts points by \param delta.
 */
template <typename T> Matrix3 translate(const math::Point2_<T> &delta);

/** Makes 3x3 transformation matrix that scales points by \param scale.
 */
template <typename T> Matrix3 scale(const math::Point2_<T> &scale);

inline Point3 transform(const Matrix3 &tr, const Point3 &pt)
{
    return Point3(
        tr(0,0)*pt(0) + tr(0,1)*pt(1) + tr(0,2),
        tr(1,0)*pt(0) + tr(1,1)*pt(1) + tr(1,2),
        tr(2,0)*pt(0) + tr(2,1)*pt(1) + tr(2,2));
}

inline Point2 transform(const Matrix3 &tr, const Point2 &pt)
{
    return Point2(
        tr(0,0)*pt(0) + tr(0,1)*pt(1) + tr(0,2),
        tr(1,0)*pt(0) + tr(1,1)*pt(1) + tr(1,2));
}

inline Extents2 transform(const Matrix3 &tr, const Extents2 &extents)
{
    return Extents2(
        transform(tr, extents.ll),
        transform(tr, extents.ur));
}

#ifdef MATH_HAS_OPENCV
template<typename T>
inline cv::Point3_<T> transform(const Matrix3 &tr, const cv::Point_<T> &pt)
{
    return cv::Point3_<T>(
        tr(0,0)*pt.x + tr(0,1)*pt.y + tr(0,2),
        tr(1,0)*pt.x + tr(1,1)*pt.y + tr(1,2),
        tr(2,0)*pt.x + tr(2,1)*pt.y + tr(2,2));
}
#endif

// in-place transformations for vectors:

inline void transform(const Matrix3 &tr, Points2 &points)
{
    for (auto &pt : points) {
        pt = transform(tr, pt);
    }
}

#ifdef MATH_HAS_OPENCV
template<typename T>
inline void transform(const Matrix3 &tr, std::vector<cv::Point_<T> > &points)
{
    for (auto &pt : points) {
        pt = transform(tr, pt);
    }
}
#endif

// template inlines

template <typename T> Matrix4 translate(const math::Point3_<T> &delta)
{
    auto m(identity4());
    m(0, 3) = delta(0);
    m(1, 3) = delta(1);
    m(2, 3) = delta(2);
    return m;
}

template <typename T> Matrix4 scale(const math::Point3_<T> &scale)
{
    auto m(identity4());
    m(0, 0) = scale(0);
    m(1, 1) = scale(1);
    m(2, 2) = scale(2);
    return m;
}

template <typename T> Matrix3 translate(const math::Point2_<T> &delta)
{
    auto m(identity3());
    m(0, 2) = delta(0);
    m(1, 2) = delta(1);
    return m;
}

template <typename T> Matrix3 scale(const math::Point2_<T> &scale)
{
    auto m(identity4());
    m(0, 0) = scale(0);
    m(1, 1) = scale(1);
    return m;
}

} // namespace math

#endif // MATH_TRANSFORM_HPP
