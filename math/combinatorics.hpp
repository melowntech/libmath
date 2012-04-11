/**
 * @file combinatorics.hpp
 * @author Ondrej Prochazka <ondrej.prochazka@citationtech.net>
 *
 * Combinatorics functions
 */

#ifndef MATH_COMBINATORICS_HPP
#define MATH_COMBINATORICS_HPP

#include <vector>

#include <sys/types.h>

namespace math {

/**
 * Return de Brujin sequence for alphabet size k and subsequences of size n
 */

class DeBruijn_t : public std::vector<uint> {

public :
    DeBruijn_t( const int k, const int n );

private :
    void db( const int t, const int p,
             const int n, const int k, std::vector<int> & a );

};


}

#endif // MATH_COMBINATORICS_HPP