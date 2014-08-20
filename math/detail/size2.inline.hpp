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

template <typename T>
inline Size2SimpleReader<T> size2SimpleReader(math::Size2_<T> &value)
{
    return { &value };
}

template<typename CharT, typename Traits, typename T>
inline std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits> &is
           , const Size2SimpleReader<T> &s)
{
    if (!(is >> s.value->width)) { return is; }

    auto x(utility::match<CharT>('x'));

    if (!(is >> x)) { return is; }
    if (x.matched) {
        return is >> s.value->height;
    }

    // use width value as height
    s.value->height = s.value->width;

    return is;
}

template<typename CharT, typename Traits, typename T>
inline std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits> &os, const Size2_<T> &s)
{
    return os << s.width << "x" << s.height;
}

template<typename CharT, typename Traits, typename T>
inline std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits> &is, math::Size2_<T> &s)
{
    return is >> s.width >> utility::expect('x') >> s.height;
}

template <typename T, typename U>
inline auto operator*(const math::Size2_<T> &l, const math::Size2_<U> &r)
    -> math::Size2_<decltype(l.width * r.width)>
{
    return { l.width * r.width, l.height * r.height };
}

template <typename T, typename U>
inline auto operator*(const math::Size2_<T> &l, const U &r)
    -> math::Size2_<decltype(l.width * r)>
{
    return { l.width * r, l.height * r };
}

template <typename T, typename U>
inline auto operator/(const math::Size2_<T> &l, const math::Size2_<U> &r)
    -> math::Size2_<decltype(l.width / r.width)>
{
    return { l.width / r.width, l.height / r.height };
}

template <typename T, typename U>
inline auto operator/(const math::Size2_<T> &l, U r)
    -> math::Size2_<decltype(l.width / r)>
{
    return { l.width / r, l.height / r };
}

template <typename T, typename U>
inline auto operator+(const math::Size2_<T> &l, const math::Size2_<U> &r)
    -> math::Size2_<decltype(l.width + r.width)>
{
    return { l.width + r.width, l.height + r.height };
}

template <typename T, typename U>
inline auto operator+(const math::Size2_<T> &l, const U &r)
    -> math::Size2_<decltype(l.width + r)>
{
    return { l.width + r, l.height + r };
}

template <typename T, typename U>
inline auto operator-(const math::Size2_<T> &l, const math::Size2_<U> &r)
    -> math::Size2_<decltype(l.width - r.width)>
{
    return { l.width - r.width, l.height - r.height };
}

template <typename T, typename U>
inline auto operator-(const math::Size2_<T> &l, const U &r)
    -> math::Size2_<decltype(l.width - r)>
{
    return { l.width - r, l.height - r };
}

template <typename T, typename U>
inline bool operator<(const math::Size2_<T> &l, const math::Size2_<U> &r)
{
    if (l.width < r.width) { return true; }
    if (r.width < l.width) { return false; }
    return l.height < r.height;
}
