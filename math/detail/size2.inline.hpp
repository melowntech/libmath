/**
 * Copyright (c) 2017 Melown Technologies SE
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

#ifdef MATH_CAN_USE_BOOST_SPIRIT
template<typename CharT, typename Traits, typename T>
inline std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits> &is, math::Size2_<T> &s)
{
    using boost::spirit::qi::auto_;
    using boost::spirit::qi::omit;
    using boost::spirit::qi::phrase_match;
    using boost::spirit::ascii::blank;
    using boost::spirit::lexeme;
    using boost::spirit::qi::skip_flag;

    typename std::basic_istream<CharT, Traits>::sentry sentry(is);

    boost::io::ios_flags_saver ifs(is);
    is.unsetf(std::ios::skipws);

    is >> phrase_match(lexeme[auto_ >> omit['x'] >> auto_]
                       , blank, skip_flag::dont_postskip
                       , s.width, s.height);

    return is;
}

template<typename CharT, typename Traits, typename T>
inline std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits> &is
           , math::Size2_<boost::rational<T>> &s)
{
    using boost::spirit::qi::auto_;
    using boost::spirit::qi::omit;
    using boost::spirit::qi::phrase_match;
    using boost::spirit::ascii::blank;
    using boost::spirit::lexeme;
    using boost::spirit::qi::skip_flag;

    typename std::basic_istream<CharT, Traits>::sentry sentry(is);

    boost::io::ios_flags_saver ifs(is);
    is.unsetf(std::ios::skipws);

    std::pair<T, T> width, height;
    is >> phrase_match(lexeme[auto_ >> omit['/'] >> auto_
                              >> omit['x']
                              >> auto_ >> omit['/'] >> auto_]
                       , blank, skip_flag::dont_postskip
                       , width.first, width.second
                       , height.first, height.second);
    if (is) {
        s.width.assign(width.first, width.second);
        s.height.assign(height.first, height.second);
    }

    return is;
}
#endif

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
