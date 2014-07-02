/**
 * @file detail/size3.inline.hpp
 * @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 * Inline stuff for math::Size3_
 * Do not include directly!
 */

#ifndef MATH_GEOMETRY_CORE_HPP_INLINES_
#    error Do not include this implementation header directly!
#endif

template <typename T>
inline Size3SimpleReader<T> size3SimpleReader(math::Size3_<T> &value)
{
    return { &value };
}

template<typename CharT, typename Traits, typename T>
inline std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits> &is
           , const Size3SimpleReader<T> &s)
{
    if (!(is >> s.value->width)) { return is; }

    auto x(utility::match<CharT>('x'));

    if (!(is >> x)) { return is; }
    if (x.matched) {
        if (!(is >> s.value->height)) { return is; };
        if (!(is >> x)) { return is; }

        if (x.matched) {
            return is >> s.value->depth;
        }
        s.value->depth = s.value->height;
        return is;
    }

    // use width value as height and depth
    s.value->height = s.value->depth = s.value->width;

    return is;
}

template<typename CharT, typename Traits, typename T>
inline std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits> &os, const Size3_<T> &s)
{
    return os << s.width << "x" << s.height << "x" << s.depth;
}

template<typename CharT, typename Traits, typename T>
inline std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits> &is, math::Size3_<T> &s)
{
    return is >> s.width >> utility::expect('x') >> s.height
              >> utility::expect('x') >> s.depth;
}

template <typename T, typename U>
inline auto operator*(const math::Size3_<T> &l, const math::Size3_<U> &r)
    -> math::Size3_<decltype(l.width * r.width)>
{
    return { l.width * r.width, l.height * r.height, l.depth * r.depth };
}

template <typename T, typename U>
inline auto operator*(const math::Size3_<T> &l, U r)
    -> math::Size3_<decltype(l.width * r)>
{
    return { l.width * r, l.height * r, l.depth * r };
}

template <typename T, typename U>
inline auto operator/(const math::Size3_<T> &l, const math::Size3_<U> &r)
    -> math::Size3_<decltype(l.width / r.width)>
{
    return { l.width / r.width, l.height / r.height, l.depth / r.depth };
}

template <typename T, typename U>
inline auto operator/(const math::Size3_<T> &l, U r)
    -> math::Size3_<decltype(l.width / r)>
{
    return { l.width / r, l.height / r, l.depth / r };
}

template <typename T, typename U>
inline auto operator+(const math::Size3_<T> &l, const math::Size3_<U> &r)
    -> math::Size3_<decltype(l.width + r.width)>
{
    return { l.width + r.width, l.height + r.height, l.depth + r.depth };
}

template <typename T, typename U>
inline auto operator+(const math::Size3_<T> &l, U r)
    -> math::Size3_<decltype(l.width + r)>
{
    return { l.width + r, l.height + r, l.depth + r };
}

template <typename T, typename U>
inline auto operator-(const math::Size3_<T> &l, const math::Size3_<U> &r)
    -> math::Size3_<decltype(l.width - r.width)>
{
    return { l.width - r.width, l.height - r.height, l.depth - r.depth };
}

template <typename T, typename U>
inline auto operator-(const math::Size3_<T> &l, U &r)
    -> math::Size3_<decltype(l.width - r)>
{
    return { l.width - r, l.height - r, l.depth - r };
}

template <typename T, typename U>
inline bool operator<(const math::Size3_<T> &l, const math::Size3_<U> &r)
{
    if (l.width < r.width) { return true; }
    if (r.width < l.width) { return false; }
    if (l.height < r.height) { return true; }
    if (r.height < l.height) { return false; }
    return l.depth < r.depth;
}
