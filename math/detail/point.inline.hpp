/**
 * Copyright (c) 2020 Melown Technologies SE
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
 * @file detail/point.inline.hpp
 * @author Vaclav Blazek <vaclav.blazek@melowntech.com>
 *
 * Point I/O.
 * Do not include directly!
 */

#ifndef MATH_GEOMETRY_CORE_HPP_INLINES_
#    error Do not include this implementation header directly!
#endif

template<class E, class T, class PT>
std::basic_istream<E, T>&
operator>>(std::basic_istream<E, T> &is, Point2_<PT> &v)
{
    ublas::vector<PT, ublas::bounded_array<PT, 2>> tmp;
    is >> tmp;
    v = tmp;
    return is;
}

template<class E, class T, class PT>
std::basic_istream<E, T>&
operator>>(std::basic_istream<E, T> &is, Point3_<PT> &v)
{
    ublas::vector<PT, ublas::bounded_array<PT, 3>> tmp;
    is >> tmp;
    v = tmp;
    return is;
}

template<class E, class T, class PT>
std::basic_istream<E, T>&
operator>>(std::basic_istream<E, T> &is, Point4_<PT> &v)
{
    ublas::vector<PT, ublas::bounded_array<PT, 4>> tmp;
    is >> tmp;
    v = tmp;
    return is;
}
