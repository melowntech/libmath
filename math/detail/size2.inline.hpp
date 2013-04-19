/**
 * @file detail/size2.inline.hpp
 * @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 * Inline stuff for math::Size2_
 * Do not include directly!
 */

#ifndef MATH_GEOMETRY_CORE_HPP_INLINES_
#    error Do not include this implementation header directly!
#endif

template <typename T, typename U>
auto operator*(const math::Size2_<T> &l, const math::Size2_<U> &r)
    -> math::Size2_<decltype(l.width * r.width)>
{
    return { l.width * r.width, l.height * r.height };
}

template <typename T, typename U>
auto operator/(const math::Size2_<T> &l, const math::Size2_<U> &r)
    -> math::Size2_<decltype(l.width / r.width)>
{
    return { l.width / r.width, l.height / r.height };
}

template <typename T, typename U>
auto operator+(const math::Size2_<T> &l, const math::Size2_<U> &r)
    -> math::Size2_<decltype(l.width + r.width)>
{
    return { l.width + r.width, l.height + r.height };
}

template <typename T, typename U>
auto operator-(const math::Size2_<T> &l, const math::Size2_<U> &r)
    -> math::Size2_<decltype(l.width - r.width)>
{
    return { l.width - r.width, l.height - r.height };
}
