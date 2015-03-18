/**
 * @file math/io.hpp
 * @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 * I/O functions for math.
 */

#include "utility/streams.hpp"

#include "./geometry_core.hpp"

namespace math {

/** IO wrapper around Point2.
 */
template <typename T> struct Point2IOWrapper;

/** Convenient generator for Point2IOWrapper.
 */
template <typename T> Point2IOWrapper<T> point2IOWrapper(math::Point2_<T> &p);

// definitions

template <typename T>
struct Point2IOWrapper {
    Point2IOWrapper(math::Point2_<T> &p) : p(&p) {}
    math::Point2_<T> *p;
};

template<typename CharT, typename Traits, typename T>
inline std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits> &is
           , const Point2IOWrapper<T> &r)
{
    return is >> (*r.p)(0) >> utility::expect(',') >> (*r.p)(1);
}

template<typename CharT, typename Traits, typename T>
inline std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits> &os
           , const Point2IOWrapper<T> &r)
{
    return os << (*r.p)(0) << ',' << (*r.p)(1);
}

template <typename T>
inline Point2IOWrapper<T> point2IOWrapper(math::Point2_<T> &p)
{
    return Point2IOWrapper<T>(p);
}

} // namespace math
