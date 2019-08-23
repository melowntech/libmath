/**
 * Copyright (c) 2019 Melown Technologies SE
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

#ifndef math_extent_hpp_included_
#define math_extent_hpp_included_

#include "geometry_core.hpp"

namespace math {

template <typename T>
struct Extent_ {
    typedef T value_type;
    typedef T point_type;

    point_type l;
    point_type r;

    Extent_() : l(), r() {}

    explicit Extent_(InvalidExtents)
        : l(std::numeric_limits<T>::max())
        , r(std::numeric_limits<T>::lowest())
    {}

    explicit Extent_(const point_type &p)
        : l(p), r(p)
    {}

    Extent_(const point_type &l, const point_type &r)
        : l(l), r(r)
    {}

    template <typename U>
    explicit Extent_(const Extent_<U> &e)
        : l(e.l), r(e.r)
    {}

    T size() const {
        return r[0] - l[0];
    }
};

typedef Extent_<int> Extenti;
typedef Extent_<double> Extentf;
typedef Extentf Extent;

template <typename T1, typename T2>
inline bool inside(const Extent_<T1> &e, const T2 &p)
{
    return ((p >= e.l) && (p <= e.r));
}

template <typename T1, typename T2>
inline Extent_<T1> operator+(const Extent_<T1> &e, const T2 &diff)
{
    return { T1(e.l - diff), T1(e.r + diff) };
}

template <typename T>
const typename Extent_<T>::point_type& l(const Extent_<T> &e) {
    return e.l;
}

template <typename T>
const typename Extent_<T>::point_type& r(const Extent_<T> &e) {
    return e.r;
}

template<typename T>
inline T center(const Extent_<T> &e) {
    return (e.r + e.l) / 2;
}

template<typename T>
inline T size(const Extent_<T> &e) {
    return (e.r - e.l);
}

template<typename T>
inline Size2_<T> size(const Extent_<T> &e, const Inclusive&) {
    return (e.r - e.l + T(1));
}

template<typename T>
inline bool empty(const Extent_<T> &e) {
    return (e.l == e.r);
}

template<typename T>
inline bool valid(const Extent_<T> &e) {
    return (e.l <= e.r);
}

template <typename T>
inline Extent_<T> unite(const Extent_<T> &a, const Extent_<T> &b) {
    return Extent_<T>(std::min(a.l, b.l), std::max(a.r, b.r));
}

template <typename T1, typename T2>
inline bool overlaps(const Extent_<T1> &a, const Extent_<T2> &b)
{
    return ((a.r > b.l) && (b.r > a.l));
}

template <typename T>
inline Extent_<T> intersect(const Extent_<T> &a, const Extent_<T> &b) {
    if (!overlaps(a, b)) {
        throw NoIntersectError
            ("Extents do not overlap, cannot compute intersection");
    }
    return Extent_<T>(std::max(a.l, b.l), std::min(a.r, b.r));
}


template <typename T>
inline void update(Extent_<T> &e, const T &x) {
    e.l = std::min(e.l, x);
    e.r = std::max(e.r, x);
}

template <typename T>
inline void update(Extent_<T> &e, const Extent_<T> &other) {
    update(e, other.l);
    update(e, other.r);
}

template <typename T, typename U>
inline bool operator==(const math::Extent_<T> &l
                       , const math::Extent_<U> &r)
{
    return (l.l == r.l) && (l.r == r.r);
}

template <typename T, typename U>
inline bool operator!=(const math::Extent_<T> &l
                       , const math::Extent_<U> &r)
{
    return !operator==(l, r);
}

template<typename CharT, typename Traits, typename T>
inline std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits> &os, const Extent_<T> &e)
{
    os << e.l << ':' << e.r(0);
    return os;
}

template<typename CharT, typename Traits, typename T>
inline std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits> &is, Extent_<T> &e)
{
    using boost::spirit::qi::auto_;
    using boost::spirit::qi::char_;
    using boost::spirit::qi::omit;
    using boost::spirit::qi::match;

    is >> match((auto_ >> omit[':'] >> auto_)
                , e.l, e.r);

    return is;
}

} // namespace math

#endif // math_extent_hpp_included_
