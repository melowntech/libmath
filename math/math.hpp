/**
 * @file math.hpp
 * @author Ondrej Prochazka <ondrej.prochazka@citationtech.net>
 *
 * Low level math helper functions
 */
      
#ifndef MATH_MATH_HPP
#define MATH_MATH_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>

namespace ublas = boost::numeric::ublas; 

namespace math {

/* square function */

template <typename T>
T sqr( T val ) { return val * val; }

/** even and odd */

template <typename T>
bool even( T val ) { return ( val >> 1 << 1 == val ); }

template <typename T>
bool odd( T val ) { return ( val >> 1 << 1 != val ); }


/** interval check */

template<typename T>
bool ccinterval( const T & lb, const T  & ub, const T & value ) {
    return ( lb <= value && value <= ub );
}
    

/* cross product */

template<typename T>
ublas::vector<T> crossProduct( 
    const ublas::vector<T> & u,   
    const ublas::vector<T> & v ) {
    
    assert( u.size() == 3 && v.size() == 3 );
    ublas::vector<T> retval( 3 );
    retval(0) = u(1) * v(2) - v(1) * u(2); 
    retval(1) = -u(0) * v(2) + v(0) * u(2);
    retval(2) = u(0) * v(1) - v(0) * u(1);
    return retval;
}
                                

/**
  * Signum function
  */
  
template <typename Value_t>
int sgn( const Value_t & value ) {
   if ( value > 0 ) return 1;
   if ( value < 0 ) return -1;
   return 0;
}
                  
/**
  * Matrix inversion
  */
                       
template<typename T, typename L, typename C>
ublas::matrix<T> matrixInvert( const ublas::matrix<T,L,C> & input ) {
                       
    using namespace boost::numeric::ublas;
    typedef permutation_matrix<std::size_t> pmatrix;
                               
    // create a working copy of the input
    matrix<T,L,C> A(input);
                                       
    // create a permutation matrix for the LU-factorization
    pmatrix pm(A.size1());
                                               
    // perform LU-factorization
    int res = lu_factorize(A,pm);
    if( res != 0 ) abort();
                                                           
    // create identity matrix of "inverse"
    matrix<T> inverse = identity_matrix<T>(A.size1());
                                                         
    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);
    return inverse;
}
                                                                               


} // namespace math

#endif // MATH_MATH_HPP
      