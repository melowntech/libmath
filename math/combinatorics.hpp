/**
 * @file combinatorics.hpp
 * @author Ondrej Prochazka <ondrej.prochazka@citationtech.net>
 *
 * Combinatorics functions
 */

#ifndef COMBINATORICS_HPP
#define COMBINATORICS_HPP

#include <vector>

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

#endif // COMBINATORICS_HPP