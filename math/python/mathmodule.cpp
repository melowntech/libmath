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
#include <mutex>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/raw_function.hpp>
#include <boost/python/slice.hpp>
#include <boost/python/call.hpp>

#include <stdint.h>

#include "dbglog/dbglog.hpp"

#include "pysupport/package.hpp"

#include "math/geometry_core.hpp"

namespace bp = boost::python;

namespace math { namespace py {

// point2

template <typename T>
bp::object Point2_repr(const math::Point2_<T> &p)
{
    std::ostringstream os;
    os << std::fixed << p;
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

    auto cls = class_<math::Point2_<T>>(name, init<math::Point2_<T>>())
        .def("__repr__", &py::Point2_repr<T>)
        .def("__getitem__", &Point2_getItem<T>)
        ;
    return cls;
}

// point3

template <typename T>
bp::object Point3_repr(const math::Point3_<T> &p)
{
    std::ostringstream os;
    os << std::fixed << p;
    return bp::str(os.str());
}

template <typename T>
T Point3_getItem(const math::Point3_<T> &p, int index)
{
    return p[index];
}

template <typename T>
bp::class_<Point3_<T>> point3(const char *name)
{
    using namespace bp;

    auto cls = class_<math::Point3_<T>>(name, init<math::Point3_<T>>())
        .def("__repr__", &py::Point3_repr<T>)
        .def("__getitem__", &Point3_getItem<T>)
        ;
    return cls;
}

// size2

template <typename T>
bp::object Size2_repr(const math::Size2_<T> &e)
{
    std::ostringstream os;
    os << std::fixed << e;
    return bp::str(os.str());
}

template <typename T>
bp::class_<Size2_<T>> size2(const char *name)
{
    using namespace bp;

    auto cls = class_<math::Size2_<T>>(name, init<math::Size2_<T>>())
        .def("__repr__", &py::Size2_repr<T>)
        .def_readonly("width", &math::Size2_<T>::width)
        .def_readonly("height", &math::Size2_<T>::height)
        .template def<bool (*)(const math::Size2_<T>&)
                      >("empty", &math::empty<T>)
        ;
    return cls;
}

// size3

template <typename T>
bp::object Size3_repr(const math::Size3_<T> &e)
{
    std::ostringstream os;
    os << std::fixed << e;
    return bp::str(os.str());
}

template <typename T>
bp::class_<Size3_<T>> size3(const char *name)
{
    using namespace bp;

    auto cls = class_<math::Size3_<T>>(name, init<math::Size3_<T>>())
        .def("__repr__", &py::Size3_repr<T>)
        .def_readonly("width", &math::Size3_<T>::width)
        .def_readonly("height", &math::Size3_<T>::height)
        .def_readonly("depth", &math::Size3_<T>::depth)
        ;
    return cls;
}

// extents2

template <typename T>
bp::object Extents2_repr(const math::Extents2_<T> &e)
{
    std::ostringstream os;
    os << std::fixed << e;
    return bp::str(os.str());
}

template <typename T>
bp::class_<Extents2_<T>> extents2(const char *name)
{
    using namespace bp;

    auto cls = class_<math::Extents2_<T>>(name, init<math::Extents2_<T>>())
        .def("__repr__", &py::Extents2_repr<T>)
        .def_readonly("ll", &math::Extents2_<T>::ll)
        .def_readonly("ur", &math::Extents2_<T>::ur)
        .template def<math::Size2_<T> (*)(const math::Extents2_<T>&)
                      >("size", &math::size<T>)
        .template def<bool (*)(const math::Extents2_<T>&)
                      >("empty", &math::empty<T>)
        ;
    return cls;
}

// extents3

template <typename T>
bp::object Extents3_repr(const math::Extents3_<T> &e)
{
    std::ostringstream os;
    os << std::fixed << e;
    return bp::str(os.str());
}

template <typename T>
bp::class_<Extents3_<T>> extents3(const char *name)
{
    using namespace bp;

    auto cls = class_<math::Extents3_<T>>(name, init<math::Extents3_<T>>())
        .def("__repr__", &py::Extents3_repr<T>)
        .def_readonly("ll", &math::Extents3_<T>::ll)
        .def_readonly("ur", &math::Extents3_<T>::ur)
        .template def<math::Size3_<T> (*)(const math::Extents3_<T>&)
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

namespace {
std::once_flag onceFlag;
} // namespace

boost::python::object import()
{
    std::call_once(onceFlag, [&]()
    {
        typedef bp::handle< ::PyObject> Handle;
        Handle module(PyInit_melown_math());

        auto package(pysupport::package());

        if (::PyModule_AddObject(package.ptr(), "math"
                                 , bp::incref(module.get())) == -1)
        {
            LOG(err2) << "PyModule_AddObject failed";
        }

        auto sys(bp::import("sys"));
        sys.attr("modules")["melown.math"] = module;
    });

    return bp::import("melown.math");
}

} } // namespace dbglog::py
