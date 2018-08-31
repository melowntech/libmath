/**
 * Copyright (c) 2018 Melown Technologies SE
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

#include <sstream>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/raw_function.hpp>
#include <boost/python/slice.hpp>
#include <boost/python/call.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include <stdint.h>

#include "dbglog/dbglog.hpp"

#include "pysupport/package.hpp"

#include "../geometry_core.hpp"
#include "../math.hpp"
#include "../geometry.hpp"
#include "../transform.hpp"

namespace bp = boost::python;

namespace math { namespace py {

// point2

template <typename T>
bp::object repr_from_ostream(const T &t)
{
    std::ostringstream os;
    os << std::fixed << t;
    return bp::str(os.str());
}

template <typename T>
T Point2_getItem(const math::Point2_<T> &p, int index)
{
    return p[index];
}

template <typename T>
void Point2_setItem(math::Point2_<T> &p, int index, T value)
{
    p[index] = value;
}

template <typename T>
math::Extents2_<T> computeExtents2(const std::vector<Point2_<T>> &points)
{
    return math::computeExtents(points.begin(), points.end());
}

template <typename T>
bp::class_<Point2_<T>> point2(const char *name, const char *listName)
{
    using namespace bp;

    typedef math::Point2_<T> Point;

    auto cls = class_<Point>
        (name, init<const Point&>())
        .def(init<>())
        .def(init<T, T>())

        .def("__repr__", &py::repr_from_ostream<Point>)
        .def("__getitem__", &Point2_getItem<T>)
        .def("__setitem__", &Point2_setItem<T>)
        ;

        // wrap vector of point2
    class_<std::vector<Point>>(listName)
        .def(vector_indexing_suite<std::vector<Point>>())
        ;

    def("computeExtents", &computeExtents2<T>);

    return cls;
}

// point3

template <typename T>
T Point3_getItem(const math::Point3_<T> &p, int index)
{
    return p[index];
}

template <typename T>
void Point3_setItem(math::Point3_<T> &p, int index, T value)
{
    p[index] = value;
}

template <typename T>
math::Extents3_<T> computeExtents3(const std::vector<Point3_<T>> &points)
{
    return math::computeExtents(points.begin(), points.end());
}

template <typename T>
bp::class_<Point3_<T>> point3(const char *name, const char *listName)
{
    using namespace bp;

    typedef math::Point3_<T> Point;

    auto cls = class_<Point>
        (name, init<const Point&>())
        .def(init<>())
        .def(init<T, T, T>())

        .def("__repr__", &py::repr_from_ostream<Point>)
        .def("__getitem__", &Point3_getItem<T>)
        .def("__setitem__", &Point3_setItem<T>)
        ;

    // wrap vector of point3
    class_<std::vector<Point>>(listName)
        .def(vector_indexing_suite<std::vector<Point>>())
        ;

    def("computeExtents", &computeExtents3<T>);

    return cls;
}

// size2

template <typename T>
bp::class_<Size2_<T>> size2(const char *name)
{
    using namespace bp;

    typedef math::Size2_<T> Size;

    auto cls = class_<Size>
        (name, init<const Size&>())
        .def(init<>())
        .def(init<T, T>())

        .def("__repr__", &py::repr_from_ostream<Size>)
        .def_readonly("width", &Size::width)
        .def_readonly("height", &Size::height)
        .template def<bool (*)(const Size&)
                      >("empty", &math::empty<T>)
        ;
    return cls;
}

// size3

template <typename T>
bp::class_<Size3_<T>> size3(const char *name)
{
    using namespace bp;

    typedef math::Size3_<T> Size;

    auto cls = class_<Size>
        (name, init<const Size&>())
        .def(init<>())
        .def(init<T, T, T>())

        .def("__repr__", &py::repr_from_ostream<Size>)
        .def_readonly("width", &Size::width)
        .def_readonly("height", &Size::height)
        .def_readonly("depth", &Size::depth)
        ;
    return cls;
}

// extents2

template <typename Extents>
std::shared_ptr<Extents> invalidExtents(void*)
{
    return std::make_shared<Extents>(math::InvalidExtents{});
}

template <typename Extents>
auto Extents_center(const Extents &extents)
    -> decltype(math::center(extents))
{
    return math::center(extents);
}

template <typename Extents>
auto Extents_size(const Extents &extents)
    -> decltype(math::size(extents))
{
    return math::size(extents);
}

template <typename T>
bp::class_<Extents2_<T>> extents2(const char *name)
{
    using namespace bp;

    typedef math::Extents2_<T> Extents;
    typedef math::Point2_<T> Point;

    auto cls = class_<Extents>
        (name, init<const Extents&>())
        .def(init<>())
        .def(init<const Point&>())
        .def(init<const Point&, const Point&>())
        .def(init<T, T, T, T>())
        .def("__init__", make_constructor(&py::invalidExtents<Extents>))

        .def("__repr__", &py::repr_from_ostream<Extents>)
        .def_readonly("ll", &Extents::ll)
        .def_readonly("ur", &Extents::ur)
        .template def<math::Size2_<T> (*)(const Extents&)
                      >("size", &math::size<T>)
        .template def<bool (*)(const Extents&)
                      >("empty", &math::empty<T>)
        .def("center", &py::Extents_center<Extents>)
        .def("size", &py::Extents_size<Extents>)
        ;
    return cls;
}

// extents3

template <typename T>
bp::class_<Extents3_<T>> extents3(const char *name)
{
    using namespace bp;

    typedef math::Extents3_<T> Extents;
    typedef math::Point3_<T> Point;

    auto cls = class_<Extents>
        (name, init<const Extents&>())
        .def(init<>())
        .def(init<const Point&>())
        .def(init<const Point&, const Point&>())
        .def(init<T, T, T, T, T, T>())
        .def("__init__", make_constructor(&py::invalidExtents<Extents>))

        .def("__repr__", &py::repr_from_ostream<Extents>)
        .def_readonly("ll", &Extents::ll)
        .def_readonly("ur", &Extents::ur)
        .template def<math::Size3_<T> (*)(const Extents&)
                      >("size", &math::size<T>)
        .template def<bool (*)(const math::Extents2_<T>&)
                      >("empty", &math::empty<T>)

        .def("center", &py::Extents_center<Extents>)
        .def("size", &py::Extents_size<Extents>)
        ;
    return cls;
}

template <typename Matrix, int size>
std::shared_ptr<Matrix> Matrix_create()
{
    return std::make_shared<Matrix>
        (boost::numeric::ublas::zero_matrix<double>(size));
}

template <typename Matrix, int size>
Matrix Matrix_eye()
{
    return Matrix(boost::numeric::ublas::identity_matrix<double>(size));
}

template <typename Matrix>
double Matrix_get(const Matrix &m, int j, int i)
{
    return m(j, i);
}

template <typename Matrix>
double Matrix_set(Matrix &m, int j, int i, double value)
{
    return m(j, i) = value;
}

template <typename Matrix>
Matrix Matrix_mul(const Matrix &m1, const Matrix &m2)
{
    return ublas::prod(m1, m2);
}

template <typename Matrix>
void Matrix_imul(Matrix &m1, const Matrix &m2)
{
    m1 = ublas::prod(m1, m2);
}

template <typename Matrix>
Matrix Matrix_invert(const Matrix &m)
{
    return math::matrixInvert(m);
}

template <typename Matrix, int size>
bp::class_<Matrix> matrix(const char *name)
{
    using namespace bp;

    auto cls = class_<Matrix>
        (name, init<const Matrix&>())
        .def("__init__", make_constructor(&py::Matrix_create<Matrix, size>))
        .def("eye", &py::Matrix_eye<Matrix, size>)
        .staticmethod("eye")
        .def("__repr__", &py::repr_from_ostream<Matrix>)
        .def("__call__", &py::Matrix_get<Matrix>)
        .def("__call__", &py::Matrix_set<Matrix>)
        .def("__mul__",  &py::Matrix_mul<Matrix>)
        .def("__imul__",  &py::Matrix_imul<Matrix>)
        .def("invert",  &py::Matrix_invert<Matrix>)
        ;

    return cls;
}

} } // namespace math::py

BOOST_PYTHON_MODULE(melown_math)
{
    using namespace bp;
    namespace py = math::py;

    // geometry core: 2D stuff
    py::point2<int>("Point2i", "Points2i");
    py::point2<float>("Point2f", "Points2f");
    py::point2<double>("Point2", "Points2");

    py::size2<int>("Size2");
    py::size2<double>("Size2f");

    py::extents2<int>("Extents2i");
    py::extents2<double>("Extents2");

    // geometry core: 3D stuff
    py::point3<int>("Point3i", "Points3i");
    py::point3<float>("Point3f", "Points3f");
    py::point3<double>("Point3", "Points3");

    py::size3<int>("Size3");
    py::size3<double>("Size3f");

    py::extents3<int>("Extents3i");
    py::extents3<double>("Extents3");

    // matrices (only doubles)
    py::matrix<math::Matrix2, 2>("Matrix2");
    py::matrix<math::Matrix3, 3>("Matrix3");
    py::matrix<math::Matrix4, 4>("Matrix4");

    // transform.hpp
    def<math::Point3 (*)(const math::Matrix4&, const math::Point3&)>
        ("transform", &math::transform);
    def<math::Point2 (*)(const math::Matrix4&, const math::Point2&)>
        ("transform", &math::transform);
    def<math::Extents2 (*)(const math::Matrix4&, const math::Extents2&)>
        ("transform", &math::transform);
    def<void (*)(const math::Matrix4&, math::Points3&)>
        ("transform", &math::transform);
    def<void (*)(const math::Matrix4&, math::Points2&)>
        ("transform", &math::transform);
}

namespace math { namespace py {
PYSUPPORT_MODULE_IMPORT(math)
} } // namespace math::py
