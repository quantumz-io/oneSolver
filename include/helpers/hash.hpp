/*! file hash.hpp
 *
 * @copyright Licensed under MIT license by hyperQ – Ewa Hendzel
 *
 */
#ifndef HASH_PAIR_HPP__
#define HASH_PAIR_HPP__

#include <cstdint>
#include <iostream>
#include <utility>

namespace helpers {

/**
 * @brief A helper structure allowing to hash pairs of objects.
 * 
 */
struct hash_pair {
  /**
   * @brief Returns a hash of hashes stored in the pair.
   * 
   * @tparam T1 Type of the first element in the pair.
   * @tparam T2 Type of the second element in the pair.
   * @param p Reference to pair to hash.
   * @return size_t Hash value.
   */

  template <class T1, class T2>
  size_t operator()(const std::pair<T1, T2> &p) const {
    auto lhs = std::hash<T1>{}(p.first);
    auto rhs = std::hash<T2>{}(p.second);
    return lhs ^ (rhs);
  }
};
} // namespace helpers
#endif