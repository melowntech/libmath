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

template <typename T>
bp::object Extents2_repr(const math::Extents2_<T> &e)
{
    std::ostringstream os;
    os << std::fixed << e;
    return bp::str(os.str());
}

} } // namespace math::py

BOOST_PYTHON_MODULE(melown_math)
{
    using namespace bp;
    namespace py = math::py;

    auto Point2 = class_<math::Point2>("Point2", init<math::Point2>())
        ;

    auto Extents2 = class_<math::Extents2>("Extents2", init<math::Extents2>())
        .def("__repr__", &py::Extents2_repr<double>)
        .def_readonly("ll", &math::Extents2::ll)
        .def_readonly("ur", &math::Extents2::ur)
        ;
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
