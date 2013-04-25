/**
 *  @file rotation.hpp
 *  @author Jakub Cerveny <jakub.cerveny@ext.citationtech.net>
 *
 *  Math related to rotations.
 */

#ifndef MATH_ROTATION_HPP
#define MATH_ROTATION_HPP

#include "geometry_core.hpp"

namespace math {


/** Converts a rotation vector (see below) to a rotation matrix.
 */
inline Matrix4 rotationMatrix(const Point3 & rvec)
{
    Matrix4 m( ublas::identity_matrix<double>(4) );

    double angle = ublas::norm_2(rvec);
    if (std::abs(angle) < 1e-10) { return m; }

    Point3 vec(rvec * (1.0 / angle));

    double s = sin(angle),
           c = cos(angle);

    double xx = vec(0) * vec(0),
           yy = vec(1) * vec(1),
           zz = vec(2) * vec(2),
           xy = vec(0) * vec(1),
           yz = vec(1) * vec(2),
           zx = vec(2) * vec(0);

    double xs = vec(0) * s,
           ys = vec(1) * s,
           zs = vec(2) * s;

    double one_c = 1.0 - c;

    m(0, 0) = (one_c * xx) + c;
    m(0, 1) = (one_c * xy) - zs;
    m(0, 2) = (one_c * zx) + ys;

    m(1, 0) = (one_c * xy) + zs;
    m(1, 1) = (one_c * yy) + c;
    m(1, 2) = (one_c * yz) - xs;

    m(2, 0) = (one_c * zx) - ys;
    m(2, 1) = (one_c * yz) + xs;
    m(2, 2) = (one_c * zz) + c;

    return m;
}


/** Converts a rotation matrix to a rotation vector, whose direction represents
 *  the axis of rotation and its magnitude represents the angle of rotation
 *  in radians.
 */
inline Point3 rotationVector(const Matrix4 & mat)
{
    double angle = acos(0.5 * (mat(0,0) + mat(1,1) + mat(2,2) - 1.0));

    Point3 axis(mat(2,1) - mat(1,2),
                mat(0,2) - mat(2,0),
                mat(1,0) - mat(0,1));

    return angle / (2*sin(angle)) * axis;
}



/** Quaternion -- represents a rotation.
 */
class Quaternion
{
public:

    Quaternion()
        : x(0.0), y(0.0), z(0.0), w(1.0) {}

    Quaternion(double x, double y, double z, double w)
        : x(x), y(y), z(z), w(w) {}

    Quaternion(const Quaternion& q)
        { x = q.x;  y = q.y;  z = q.z;  w = q.w; }

    Quaternion& operator=(const Quaternion& q)
        { x = q.x;  y = q.y;  z = q.z;  w = q.w;  return *this; }

    Quaternion(const Matrix4& mat)
        { fromMatrix(mat); }

    /// Converts a matrix to a quaternion.
    void fromMatrix(const Matrix4& m);

    /// Converts this quaternion to a matrix.
    Matrix4 toMatrix() const;

    Quaternion& normalize();

    Quaternion operator+(const Quaternion& other) const;
    Quaternion operator*(const Quaternion& other) const;
    Quaternion operator*(double scalar) const;
    Quaternion& operator*=(double scalar);
    Quaternion& operator*=(const Quaternion& other);

    double x, y, z; // imaginary part
    double w; // real part
};


inline Quaternion Quaternion::operator*(double s) const
{
    return Quaternion(s*x, s*y, s*z, s*w);
}

inline Quaternion& Quaternion::operator*=(double s)
{
    x*=s;  y*=s;  z*=s;  w*=s;
    return *this;
}

inline Quaternion& Quaternion::operator*=(const Quaternion& other)
{
    return (*this = other * (*this));
}

inline Quaternion Quaternion::operator+(const Quaternion& b) const
{
    return Quaternion(x+b.x, y+b.y, z+b.z, w+b.w);
}

inline Quaternion Quaternion::operator*(const Quaternion& other) const
{
    Quaternion tmp;
    tmp.w = (other.w * w) - (other.x * x) - (other.y * y) - (other.z * z);
    tmp.x = (other.w * x) + (other.x * w) + (other.y * z) - (other.z * y);
    tmp.y = (other.w * y) + (other.y * w) + (other.z * x) - (other.x * z);
    tmp.z = (other.w * z) + (other.z * w) + (other.x * y) - (other.y * x);
    return tmp;
}


inline void Quaternion::fromMatrix(const Matrix4& m)
{
    const double diag = m(0,0) + m(1,1) + m(2,2) + 1;

    if (diag > 0.0)
    {
        const double scale = sqrt(diag) * 2.0; // get scale from diagonal

        x = (m(2,1) - m(1,2)) / scale;
        y = (m(0,2) - m(2,0)) / scale;
        z = (m(1,0) - m(0,1)) / scale;
        w = 0.25 * scale;
    }
    else
    {
        if (m(0,0) > m(1,1) && m(0,0) > m(2,2))
        {
            // 1st element of diag is greatest value
            // find scale according to 1st element, and double it
            const double scale = sqrt( 1.0 + m(0,0) - m(1,1) - m(2,2)) * 2.0;

            x = 0.25 * scale;
            y = (m(0,1) + m(1,0)) / scale;
            z = (m(2,0) + m(0,2)) / scale;
            w = (m(2,1) - m(1,2)) / scale;
        }
        else if (m(1,1) > m(2,2))
        {
            // 2nd element of diag is greatest value
            // find scale according to 2nd element, and double it
            const double scale = sqrt( 1.0 + m(1,1) - m(0,0) - m(2,2)) * 2.0;

            x = (m(0,1) + m(1,0) ) / scale;
            y = 0.25 * scale;
            z = (m(1,2) + m(2,1) ) / scale;
            w = (m(0,2) - m(2,0) ) / scale;
        }
        else
        {
            // 3rd element of diag is greatest value
            // find scale according to 3rd element, and double it
            const double scale = sqrt( 1.0 + m(2,2) - m(0,0) - m(1,1)) * 2.0;

            x = (m(0,2) + m(2,0)) / scale;
            y = (m(1,2) + m(2,1)) / scale;
            z = 0.25 * scale;
            w = (m(1,0) - m(0,1)) / scale;
        }
    }

    normalize();
}


inline Quaternion& Quaternion::normalize()
{
    double n = x*x + y*y + z*z + w*w;
    *this *= 1.0 / sqrt(n);
    return *this;
}


inline Matrix4 Quaternion::toMatrix() const
{
    Matrix4 m = ublas::identity_matrix<double>(4);

    m(0,0) = 1.0 - 2.0*y*y - 2.0*z*z;
    m(1,0) = 2.0*x*y + 2.0*z*w;
    m(2,0) = 2.0*x*z - 2.0*y*w;
    m(3,0) = 0.0;

    m(0,1) = 2.0*x*y - 2.0*z*w;
    m(1,1) = 1.0 - 2.0*x*x - 2.0*z*z;
    m(2,1) = 2.0*z*y + 2.0*x*w;
    m(3,1) = 0.0;

    m(0,2) = 2.0*x*z + 2.0*y*w;
    m(1,2) = 2.0*z*y - 2.0*x*w;
    m(2,2) = 1.0 - 2.0*x*x - 2.0*y*y;
    m(3,2) = 0.0;

    m(0,3) = 0.0;
    m(1,3) = 0.0;
    m(2,3) = 0.0;
    m(3,3) = 1.0;

    return m;
}

} // namespace math

#endif
