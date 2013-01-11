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


namespace math {

namespace ublas = boost::numeric::ublas;

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
                       
template <typename T, typename L, typename C>
ublas::matrix<T,L,C> matrixInvert( const ublas::matrix<T,L,C> & input ) {
                       
    typedef ublas::permutation_matrix<std::size_t> pmatrix;
                               
    // create a working copy of the input
    ublas::matrix<T,L,C> A(input);
                                       
    // create a permutation matrix for the LU-factorization
    pmatrix pm(A.size1());
                                               
    // perform LU-factorization
    int res = lu_factorize(A,pm);
    if( res != 0 )
    {
        std::cerr << "Singular matrix in math::matrixInvert. Aborting.";
        abort();
    }
                                                           
    // create identity matrix of "inverse"
    ublas::matrix<T,L,C> inverse = ublas::identity_matrix<T>(A.size1());
                                                         
    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);
    return inverse;
}
                                                                               
/**
 * Matrix determinant
 */

template <class matrix_T>
double determinant(ublas::matrix_expression<matrix_T> const& mat_r)
{
  double det = 1.0;

  matrix_T mLu(mat_r() );
  ublas::permutation_matrix<std::size_t> pivots(mat_r().size1() );

  int is_singular = lu_factorize(mLu, pivots);

  if (!is_singular)
  {
    for (std::size_t i=0; i < pivots.size(); ++i)
    {
      if (pivots(i) != i)
        det *= -1.0;

      det *= mLu(i,i);
    }
  }
  else
    det = 0.0;

  return det;
} 

                                                                               
} // namespace math

#endif // MATH_MATH_HPP
      
