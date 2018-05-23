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
bp::class_<Point2_<T>> point2(const char *name)
{
    using namespace bp;

    typedef math::Point2_<T> Point;

    auto cls = class_<Point>
        (name, init<const Point&>())
        .def(init<>())
        .def(init<T, T>())

        .def("__repr__", &py::repr_from_ostream<Point>)
        .def("__getitem__", &Point2_getItem<T>)
        ;
    return cls;
}

// point3

template <typename T>
T Point3_getItem(const math::Point3_<T> &p, int index)
{
    return p[index];
}

template <typename T>
bp::class_<Point3_<T>> point3(const char *name)
{
    using namespace bp;

    typedef math::Point3_<T> Point;

    auto cls = class_<Point>
        (name, init<const Point&>())
        .def(init<>())
        .def(init<T, T, T>())

        .def("__repr__", &py::repr_from_ostream<Point>)
        .def("__getitem__", &Point3_getItem<T>)
        ;
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
        ;
    return cls;
}

} } // namespace math::py

BOOST_PYTHON_MODULE(melown_math)
{
    using namespace bp;
    namespace py = math::py;

    // geometry core: 2D stuff
    py::point2<int>("Point2i");
    py::point2<float>("Point2f");
    py::point2<double>("Point2");

    py::size2<int>("Size2");
    py::size2<double>("Size2f");

    py::extents2<int>("Extents2i");
    py::extents2<double>("Extents2");

    // geometry core: 3D stuff
    py::point3<int>("Point3i");
    py::point3<float>("Point3f");
    py::point3<double>("Point3");

    py::size3<int>("Size3");
    py::size3<double>("Size3f");

    py::extents3<int>("Extents3i");
    py::extents3<double>("Extents3");
}

namespace math { namespace py {
PYSUPPORT_MODULE_IMPORT(math)
} } // namespace math::py
