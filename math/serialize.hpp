/**
 * @file serialization.hpp
 * @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 * Implementation of all serialization methods found in math library.
 * You have to link your binary with Boost.Serialization library.
 */


#ifndef MATH_SERIALIZATION_HPP
#define MATH_SERIALIZATION_HPP

#include <math/geometry_core.hpp>

namespace math {

template<typename Archive, typename T>
inline void serialize(Archive &ar, Size2_<T> &t
                      , const unsigned int version)
{
    (void) version;

    ar & t.width;
    ar & t.height;
}

template<class Archive, typename T>
inline void serialize(Archive &ar, Point2_<T> &t
                      , const unsigned int version)
{
    (void) version;

    ar & t(0);
    ar & t(1);
}

template<class Archive, typename T>
inline void serialize(Archive &ar, Point3_<T> &t
                      , const unsigned int version)
{
    (void) version;

    ar & t(0);
    ar & t(1);
    ar & t(2);
}

} // namespace optics

#endif // MATH_SERIALIZATION_HPP
