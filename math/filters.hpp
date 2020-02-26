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
 * @file filters.hpp
 * @author Ondrej Prochazka <ondrej.prochazka@citationtech.net>
 *
 * FIR filters
 */

#ifndef MATH_FILTERS_HPP
#define MATH_FILTERS_HPP

#include <iomanip>
#include <sstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm> // std::min

namespace math {

template <typename T> struct FilterTraits{};

/**
 * Abstract class providing horizontal image filtering with a static convolution
 * kernel.
 */

class FIRFilter_t {

public :

    FIRFilter_t( const unsigned int halforder ) : halforder( halforder ) {
        order = 2 * halforder + 1;
        kernel = new double[ order ];
        std::fill( kernel, kernel + order, 0 );
    }

    template<typename Filter1>
    FIRFilter_t(FilterTraits<Filter1>,const double cutOffX){
        Filter1 filter(cutOffX);

        halforder = std::floor(filter.halfwinx());
        order = 2 * halforder + 1;
        kernel = new double[ order ];

        for(unsigned int i=0; i<order; ++i){
            kernel[i]=filter((int)i - (int)halforder);
        }
    }

    ~FIRFilter_t() { delete[] kernel; };

    std::string dump() {
        std::ostringstream str;
        for ( unsigned int i = 0; i < order; i++ )
            str << "\t" << i << "\t"
            << std::setprecision( 10 ) << std::fixed << kernel[ i ] << "\n";
        return str.str();
    }

    /**
     * Obtain a single filter response value for a given input. The input is represented
     * by an iterator.
     * This is in fact a convolution operation, with filter used as a convolution kernel.
     * xpos and rowsize parameters signal the position within input sequence and the length
     * of the sequence.
     */
    template <typename Iterator_t>
    double convolute( const Iterator_t & pos,
                        unsigned int xpos, unsigned int rowSize ) const {

        double weightSum = 0.0;
        double valueSum = 0.0;

        Iterator_t begin = pos - std::min( xpos, halforder );
        Iterator_t end = pos + std::min( rowSize - 1 - xpos, halforder );


        int index = std::max( (int) ( halforder - xpos ), (int) 0 );

        for ( Iterator_t it = begin; it <= end; ++it ) {
            valueSum += it[0] * kernel[ index ];
            weightSum += kernel[ index ];
            index++;
        }

        return valueSum / weightSum;
    }


    /**
     * boost::gil view iterators need a special treatment, even if they are grayscale views.
     * I don't know why.
     */
    template <typename Iterator_t>
    float gilConvolute( const Iterator_t & pos,
                        unsigned int xpos, unsigned int rowSize ) const {

        float weightSum = 0.0;
        float valueSum = 0.0;

        Iterator_t begin = pos - std::min( xpos, halforder );
        Iterator_t end = pos + std::min( rowSize - 1 - xpos, halforder );


        int index = std::max( (int) ( halforder - xpos ), (int) 0 );

        for ( Iterator_t it = begin; it <= end; ++it ) {
            valueSum += it[0][0] * kernel[ index ];
            weightSum += kernel[ index ];
            index++;
        }

        return valueSum / weightSum;
    }

    /** Fetches kernel as a vector.
     */
    std::vector<double> getKernel() const {
        return std::vector<double>(kernel, kernel + order);
    }

protected :

    unsigned int halforder, order;
    double * kernel;
};

/**
 * Low-pass filter constructed from sinc function by applying a hamming window.
 */

class LowPassFilter_t : public FIRFilter_t {

public :
    LowPassFilter_t( float cutoffPeriod, unsigned int halfwindow )
    : FIRFilter_t( halfwindow ) {

        for ( unsigned int i = 0; i <= halforder; i++ ) {

            auto sinci = i == 0 ? 2 / cutoffPeriod :
                1.0 / ( M_PI * i ) * sin( 2 * M_PI * i / cutoffPeriod );
            auto hamming = ( order == 1 ) ? 1 :
                0.54 + 0.46 * cos( 2 * M_PI * i / ( order - 1 ) );
            kernel[ halforder - i ]
                = kernel[ halforder + i ] = hamming * sinci;

        }
    }
};

/*
 * Box filter, used for simple average.
 */

class BoxFilter : public FIRFilter_t {

public:
    BoxFilter( unsigned int halfwindow ) : FIRFilter_t( halfwindow ) {

        for ( unsigned int i = 0; i <= halforder; i++ )
            kernel[ halforder - i ] = 1.0 /  ( 2 * halforder + 1 );
    }
};



namespace detail {
    /**
     * Def is set of parameters defining a shape of a spline.
     * It models the following concepts:
     *
     * struct DefConcept{
     *     static constexpr double B;
     *     static constexpr double C;
     * };
     */

    struct CatmullRom2Spline {
        static constexpr double B = .0;
        static constexpr double C = .5;
    };
} // namespace detail


/**
 * Filter1 is a 1D analytic window filter, used in image processing.
 * It models the following concepts:
 *
 * class Filter1Concept {
 * public:
 *      double operator() ( const double x ) const;
 *      double halfwinx() const;
 * };
 *
 *
 * class LowPassFilter1Concept {
 * public :
 *      LowPassFilter1Concept( double cutoffX )
 *      double operator() ( const double x ) const;
 *      double halfwinx() const;
 * }
 */

class Box1 {
public:
    Box1( double cutOffX)
        : hwinx_( cutOffX / 4 ) {
        value_ = 1.0 / ( hwinx_ );
    }

    double operator() ( double x) const {
        if ( fabs( x ) <= hwinx_ )
            return value_;
        else
            return 0.0;
    }
    double halfwinx() const { return hwinx_; }

    Box1 scaled(double scaleX) const { return { 4.0 * hwinx_ * scaleX }; }

private:
    double hwinx_;
    double value_;
};

class SincHamming1 {

public:
    SincHamming1( double cutoffX)
        : hwinx( cutoffX ), cutoffX( cutoffX ){}

    SincHamming1( double cutoffX, double hwinx)
        : hwinx( hwinx ), cutoffX( cutoffX ) {}

    double operator() ( const double x ) const {

        return sinc( x, cutoffX ) * hamming( x, hwinx );
    }

    double halfwinx() const { return hwinx; }

    SincHamming1 scaled(double scaleX) const { return { cutoffX * scaleX }; }

private :

    static double sinc( double x, double cutoff )  {
        if ( fabs( x ) < 1E-15 )
            return 2.0 / cutoff;
        else
            return 1.0 / ( M_PI * x ) * sin( 2 * M_PI * x / cutoff );
    }

    static double hamming( double x, double hwin ) {
        return fabs( x ) > hwin ? 0.0 : 0.54 + 0.46 * cos( M_PI * x / hwin );
    }

    double hwinx, cutoffX;
};

class Lanczos1 {

public:

    Lanczos1( double cutoffX)
        : hwinx( cutoffX ), cutoffX( cutoffX ){}

    Lanczos1( double cutoffX, double hwinx)
        : hwinx( hwinx ), cutoffX( cutoffX ) {}

    double operator() ( const double x) const {
        if ( fabs( x ) <= hwinx )
            return
                sinc( x, cutoffX ) * sinc( x, 2.0 * hwinx );
            else
                return 0;
    }

    double halfwinx() const { return hwinx; }

    Lanczos1 scaled(double scaleX) const { return { cutoffX * scaleX }; }

private :

    static double sinc( double x, double cutoff )  {
        if ( fabs( x ) < 1E-15 )
            return 2.0 / cutoff;
        else
            return 1.0 / ( M_PI * x ) * sin( 2 * M_PI * x / cutoff );
    }

    double hwinx, cutoffX;
};

class Linear1 {
public:
    Linear1(double cutoffX)
        : hwinx(cutoffX / 2.), x2_T(2. / cutoffX), x4_T2(x2_T * x2_T)
    {}

    double operator()(const double x) const {
        auto ax(std::abs(x));
        return (((ax >= hwinx) ? 0. : (x2_T - x4_T2 * ax)) );
    }

    double halfwinx() const { return hwinx; }

    Linear1 scaled(double scaleX) const { return { hwinx * 2.0 * scaleX }; }

private:
    double hwinx;
    double x2_T, x4_T2; // x coefficients: 2/T and (2/T)^2
};

template <typename Def>
class Cubic1 {
public:
    Cubic1(double cutoffX)
        : cutoffX(cutoffX)
        , x2_T(2. / cutoffX)
    {}

    double operator()(const double x) const {
        return value(std::abs(x * x2_T));
    }

    double halfwinx() const { return cutoffX; }

    Cubic1<Def> scaled(double scaleX) const { return { cutoffX * scaleX }; }

protected:
    static inline double value(const double x) {
        if (x < 1.) {
            return (x * (x * x * (12. - 9. * Def::B - 6. * Def::C)
                         + x * (-18. + 12. * Def::B + 6. * Def::C))
                    + 6. - 2. * Def::B) / 6.;
        } else if (x < 2.) {
            return (x * (x * (x * (-Def::B - 6. * Def::C)
                              + (6. * Def::B + 30. * Def::C))
                         + (-12. * Def::B - 48. * Def::C))
                    + (8. * Def::B + 24. * Def::C)) / 6.;
        }
        return 0.;
    }

    double cutoffX;
    double x2_T;
};

class CatmullRom1 : public Cubic1<detail::CatmullRom2Spline> {
public:
    CatmullRom1(double cutoffX)
        : Cubic1<detail::CatmullRom2Spline>(cutoffX)
    {}

    CatmullRom1 scaled(double scaleX) const { return { cutoffX * scaleX }; }
};


/**
 * Filter2 is a 2D analytic window filter, used in image processing.
 * It models the following concepts:
 *
 * class Filter2Concept {
 * public:
 *      double operator() ( const double x, const double y ) const;
 *      double halfwinx() const;
 *      double halfwiny() const;
 * };
 *
 *
 * class LowPassFilter2Concept {
 * public :
 *      LowPassFilter2Concept( double cutoffX, double cutoffY )
 *      double operator() ( const double x, const double y ) const;
 *      double halfwinx() const;
 *      double halfwiny() const;
 * }
 */

class SincHamming2 {

public:
    SincHamming2( double cutoffX, double cutoffY  )
        : hwinx( cutoffX ), hwiny( cutoffY ),
          cutoffX( cutoffX ), cutoffY( cutoffY ) {}

    SincHamming2( double cutoffX, double hwinx, double cutoffY, double hwiny )
        : hwinx( hwinx ), hwiny( hwiny ),
          cutoffX( cutoffX ), cutoffY( cutoffY ) {}

    double operator() ( const double x, const double y ) const {

        return sinc( x, cutoffX ) * hamming( x, hwinx )
            * sinc( y, cutoffY ) * hamming( y, hwiny );
    }

    double halfwinx() const { return hwinx; }
    double halfwiny() const { return hwiny; }


private :

    static double sinc( double x, double cutoff )  {
        if ( fabs( x ) < 1E-15 )
            return 2.0 / cutoff;
        else
            return 1.0 / ( M_PI * x ) * sin( 2 * M_PI * x / cutoff );
    }

    static double hamming( double x, double hwin ) {
        return fabs( x ) > hwin ? 0.0 : 0.54 + 0.46 * cos( M_PI * x / hwin );
    }

    double hwinx, hwiny, cutoffX, cutoffY;
};


class Lanczos2 {

public:

    Lanczos2( double cutoffX, double cutoffY )
        : hwinx( cutoffX ), hwiny( cutoffY ),
          cutoffX( cutoffX ), cutoffY( cutoffY ) {}

    Lanczos2( double cutoffX, double hwinx, double cutoffY, double hwiny )
        : hwinx( hwinx ), hwiny( hwiny ),
          cutoffX( cutoffX ), cutoffY( cutoffY ) {}

    double operator() ( const double x, const double y ) const {

        if ( fabs( x ) <= hwinx && fabs( y ) <= hwiny )
            return
                sinc( x, cutoffX ) * sinc( x, 2.0 * hwinx )
                * sinc( y, cutoffY ) * sinc( y, 2.0 * hwiny );
            else
                return 0;

    }

    double halfwinx() const { return hwinx; }
    double halfwiny() const { return hwiny; }

    Lanczos2 scaled(double scaleX, double scaleY) const {
        return { cutoffX * scaleX, cutoffY * scaleY };
    }

private :

    static double sinc( double x, double cutoff )  {
        if ( fabs( x ) < 1E-15 )
            return 2.0 / cutoff;
        else
            return 1.0 / ( M_PI * x ) * sin( 2 * M_PI * x / cutoff );
    }

    double hwinx, hwiny, cutoffX, cutoffY;
};


class Box2 {

public:
    Box2( double hwinx, double hwiny )
        : hwinx_( hwinx ), hwiny_( hwiny ) {
        value_ = 1.0 / ( hwinx_ * hwiny_ );
    }

    double operator() ( double x, double y ) const {
        if ( fabs( x ) <= hwinx_ && fabs( y ) <= hwiny_ )
            return value_;
        else
            return 0.0;
    }

    double halfwinx() const { return hwinx_; }
    double halfwiny() const { return hwiny_; }

    Box2 scaled(double scaleX, double scaleY) const {
        return { hwinx_ * scaleX, hwiny_ * scaleY };
    }

private:
    double hwinx_, hwiny_;
    double value_;

};

class Linear2 {
public:
    Linear2(double cutoffX, double cutoffY)
        : hwinx(cutoffX / 2.), hwiny(cutoffY / 2.)
        , x2_T(2. / cutoffX), x4_T2(x2_T * x2_T)
        , y2_T(2. / cutoffY), y4_T2(y2_T * y2_T)
    {}

    double operator()(const double x, const double y) const {
        auto ax(std::abs(x));
        auto ay(std::abs(y));
        return (((ax >= hwinx) ? 0. : (x2_T - x4_T2 * ax))
                * ((ay >= hwiny) ? 0. : (y2_T - y4_T2 * ay)));
    }

    double halfwinx() const { return hwinx; }
    double halfwiny() const { return hwiny; }

    Linear2 scaled(double scaleX, double scaleY) const {
        return { 2.0 * hwinx * scaleX, 2.0 * hwiny * scaleY };
    }

private:
    double hwinx, hwiny;
    double x2_T, x4_T2; // x coefficients: 2/T and (2/T)^2
    double y2_T, y4_T2; // y coefficients: 2/T and (2/T)^2
};

template <typename Def>
class Cubic2 {
public:
    Cubic2(double cutoffX, double cutoffY)
        : cutoffX(cutoffX), cutoffY(cutoffY)
        , x2_T(2. / cutoffX), y2_T(2. / cutoffY)
    {}

    double operator()(const double x, const double y) const {
        return value(std::abs(x * x2_T)) * value(std::abs(y * y2_T));
    }

    double halfwinx() const { return cutoffX; }
    double halfwiny() const { return cutoffY; }

    Cubic2<Def> scaled(double scaleX, double scaleY) const {
        return { cutoffX * scaleX, cutoffY * scaleY };
    }

protected:
    static inline double value(const double x) {
        if (x < 1.) {
            return (x * (x * x * (12. - 9. * Def::B - 6. * Def::C)
                         + x * (-18. + 12. * Def::B + 6. * Def::C))
                    + 6. - 2. * Def::B) / 6.;
        } else if (x < 2.) {
            return (x * (x * (x * (-Def::B - 6. * Def::C)
                              + (6. * Def::B + 30. * Def::C))
                         + (-12. * Def::B - 48. * Def::C))
                    + (8. * Def::B + 24. * Def::C)) / 6.;
        }
        return 0.;
    }

    double cutoffX, cutoffY;
    double x2_T, y2_T;
};

class CatmullRom2
    : public Cubic2<detail::CatmullRom2Spline>
{
public:
    CatmullRom2(double cutoffX, double cutoffY)
        : Cubic2<detail::CatmullRom2Spline>(cutoffX, cutoffY)
    {}

    CatmullRom2 scaled(double scaleX, double scaleY) const {
        return { cutoffX * scaleX, cutoffY * scaleY };
    }
};

/**
 * Filter3 is a 3D analytic window filter, used in volume processing.
 * It models the following concepts:
 *
 * class Filter3Concept {
 * public:
 *      double operator() ( const double x, const double y, const double z ) const;
 *      double halfwinx() const;
 *      double halfwiny() const;
 *      double halfwinz() const;
 * };
 *
 *
 * class LowPassFilter3Concept {
 * public :
 *      LowPassFilter3Concept( double cutoffX, double cutoffY, double cutoffZ )
 *      double operator() ( const double x, const double y, const double z ) const;
 *      double halfwinx() const;
 *      double halfwiny() const;
 *      double halfwinz() const;
 * }
 */

template <typename Def>
class Cubic3 {
public:
    Cubic3(double cutoffX, double cutoffY, double cutoffZ)
        : cutoffX(cutoffX), cutoffY(cutoffY), cutoffZ(cutoffZ)
        , x2_T(2. / cutoffX), y2_T(2. / cutoffY), z2_T(2. / cutoffZ)
    {}

    double operator()(const double x, const double y, const double z) const {
        return value(std::abs(x * x2_T)) * value(std::abs(y * y2_T))
                * value(std::abs(z * z2_T));
    }

    double halfwinx() const { return cutoffX; }
    double halfwiny() const { return cutoffY; }
    double halfwinz() const { return cutoffZ; }

    Cubic3<Def> scaled(double scaleX, double scaleY, double scaleZ) const {
        return { cutoffX * scaleX, cutoffY * scaleY, cutoffZ * scaleZ };
    }

protected:
    static inline double value(const double x) {
        if (x < 1.) {
            return (x * (x * x * (12. - 9. * Def::B - 6. * Def::C)
                         + x * (-18. + 12. * Def::B + 6. * Def::C))
                    + 6. - 2. * Def::B) / 6.;
        } else if (x < 2.) {
            return (x * (x * (x * (-Def::B - 6. * Def::C)
                              + (6. * Def::B + 30. * Def::C))
                         + (-12. * Def::B - 48. * Def::C))
                    + (8. * Def::B + 24. * Def::C)) / 6.;
        }
        return 0.;
    }

    double cutoffX, cutoffY, cutoffZ;
    double x2_T, y2_T, z2_T;
};

class CatmullRom3 : public Cubic3<detail::CatmullRom2Spline> {
public:
    CatmullRom3(double cutoffX, double cutoffY, double cutoffZ)
        : Cubic3<detail::CatmullRom2Spline>(cutoffX, cutoffY, cutoffZ)
    {}

    CatmullRom3 scaled(double scaleX, double scaleY, double scaleZ) const {
        return { cutoffX * scaleX, cutoffY * scaleY, cutoffZ * scaleZ };
    }
};


} // namespace math

#endif // MATH_FILTERS_HPP

