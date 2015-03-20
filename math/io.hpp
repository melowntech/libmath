/**
 * @file math/io.hpp
 * @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 * I/O functions for math.
 */

#include "utility/streams.hpp"

#include "dbglog/dbglog.hpp"

#include "./geometry_core.hpp"

namespace math {

/** IO wrappers around Point2.
 */

/** Non-const Point2_ wrapper.
 *
 *  Input/Output format: %d,%d
 *
 *  If there is no bound point all changes are made on local member.
 *
 *  Assignment to wrapper changes:
 *      * bound point if wrapper is bound to valid point
 *      * local point if wrapper is not bound
 */
template <typename T> struct Point2IOWrapper;

/** Const Point2_ wrapper.
 *
 *  Input/Output format: %d,%d
 *
 *  Allows only output.
 */
template <typename T> struct ConstPoint2IOWrapper;

/** Convenient generator for Point2IOWrapper.
 *  Returns proper type of Point2IOWrapper based on type of Point2_
 *  Binds wrapper to this point.
 */
template <typename T>
Point2IOWrapper<T> point2IOWrapper(Point2_<T> &p);

/** Convenient generator for ConstPoint2IOWrapper.
 *  Returns proper type of Point2IOWrapper based on type of Point2_
 */
template <typename T>
ConstPoint2IOWrapper<T> point2IOWrapper(const Point2_<T> &p);

// definitions

template <typename T>
struct Point2IOWrapper {
    Point2IOWrapper() : p_() {}
    Point2IOWrapper(Point2_<T> &p) : p_(&p) {}
    Point2IOWrapper(const Point2IOWrapper &other);
    Point2IOWrapper& operator=(const Point2IOWrapper &other);

    Point2_<T> &get() const;

private:
    Point2_<T> *p_;
    mutable Point2_<T> self_;
};

template <typename T>
struct ConstPoint2IOWrapper {
    ConstPoint2IOWrapper(const Point2_<T> &p) : p_(&p) {}
    const Point2_<T> &get() const { return *p_; }

private:
    const Point2_<T> *p_;
};

template <typename T>
inline Point2IOWrapper<T>::Point2IOWrapper(const Point2IOWrapper &other)
    : p_(other.p_), self_(other.self_) {}

template <typename T>
inline Point2IOWrapper<T>&
Point2IOWrapper<T>::operator=(const Point2IOWrapper &other) {
    if (&other != this) {
        if (p_) {
            *p_ = other.get();
        } else {
            self_ = other.get();
        }
    }
    return *this;
}

template <typename T>
Point2_<T>& Point2IOWrapper<T>::get() const { return p_ ? *p_ : self_; }

template<typename CharT, typename Traits, typename T>
inline std::basic_istream<CharT, Traits>&
operator>>(std::basic_istream<CharT, Traits> &is
           , const Point2IOWrapper<T> &r)
{
    return is >> r.get()(0) >> utility::expect(',') >> r.get()(1);
}

template<typename CharT, typename Traits, typename T>
inline std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits> &os
           , const Point2IOWrapper<T> &r)
{
    return os << r.get()(0) << ',' << r.get()(1);
}

template<typename CharT, typename Traits, typename T>
inline std::basic_ostream<CharT, Traits>&
operator<<(std::basic_ostream<CharT, Traits> &os
           , const ConstPoint2IOWrapper<T> &r)
{
    return os << r.get()(0) << ',' << r.get()(1);
}

template <typename T>
inline Point2IOWrapper<T> point2IOWrapper(Point2_<T> &p)
{
    return Point2IOWrapper<T>(p);
}

template <typename T>
inline ConstPoint2IOWrapper<T> point2IOWrapper(const Point2_<T> &p)
{
    return ConstPoint2IOWrapper<T>(p);
}

} // namespace math
