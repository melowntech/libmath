/**
 * @file detail/extents2.inline.hpp
 * @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 * Inline stuff for math::Extents2_
 * Do not include directly!
 */

#ifndef MATH_GEOMETRY_CORE_HPP_INLINES_
#    error Do not include this implementation header directly!
#endif

template <typename T, typename U>
inline bool operator==(const math::Extents2_<T> &l
                       , const math::Extents2_<U> &r)
{
    return (l.ll == r.ll) && (l.ur == r.ur);
}

template <typename T, typename U>
inline bool operator!=(const math::Extents2_<T> &l
                       , const math::Extents2_<U> &r)
{
    return !operator==(l, r);
}
