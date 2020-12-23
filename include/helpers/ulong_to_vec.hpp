/*! file ulong_to_vec.hpp
 *
 * @copyright Licensed under MIT license by hyperQ – Ewa Hendzel
 *
 */
#ifndef ULONG_TO_VEC_HPP_
#define ULONG_TO_VEC_HPP_

#include <iostream>
#include <vector>

namespace helpers {
/**
 * @brief Convert a ulong value representing a binary string of up to 32
 * bits into vector of chars containing individual bit values of the
 * first n_bits.
 *
 * @param val A ulong value representing binary string.
 * @param n_bits Number of bits to be store in the vector.
 * @return std::vector<char> Vector of zeros and ones containing binary
 * representation of first n_bits of val.
 */
std::vector<char> ulong_to_vec(ulong val, uint n_bits) {
  if (n_bits > 32) {
    throw std::invalid_argument("state can be up to 32 bit");
  } 
  std::vector<char> state;
  for (auto i = 0; i < n_bits; ++i) {
    state.push_back((val >> i) & 1);
  }
  return state;
}

}; // namespace helpers

#endif