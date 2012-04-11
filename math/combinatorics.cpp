/*
 * filters.cpp
 */
 
#include "combinatorics.hpp"

namespace math {

/* class deBruijn_t */

DeBruijn_t::DeBruijn_t( int k, int n ) {

    std::vector<int> a = std::vector<int>( k * n, 0 );
    db( 1, 1, n, k, a );
}

void DeBruijn_t::db( const int t, const int p, const int n, const int k,
                std::vector<int> & a ) {
    if ( t > n ) {

        if ( n % p == 0 )
            for ( int j = 1; j < p + 1; j++ ) push_back( a[j] );
    } else {

        a[t] = a[t - p];
        db( t + 1, p, n, k, a );
        for ( int j = a[ t - p ] + 1; j < k; j++ ) {
            a[t] = j;
            db( t + 1, t, n, k, a );
        }
    }
}


} // namespace math

