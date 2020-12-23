/*! file insert.hpp
 *
 * @copyright Licensed under MIT license by hyperQ – Ewa Hendzel
 *
 */
#ifndef INSERT_HPP_
#define INSERT_HPP_

#include <unordered_map>

namespace helpers {

/**
 * @brief A helper function that inserts linear or quadratic
 * coefficients int o QUBO model.
 *
 * @tparam NodeType An integer number type used to index nodes of the
 * QUBO model.
 * @tparam CoefType A real type used to store values of the QUBO model.
 * @tparam Hash Type of the pair storing linear or quadratic
 * coefficients of the QUBO model.
 * @param um Linear or quadratic coefficients store.
 * @param key Address of a variable or addresses pair of the variables
 * to add to the model.
 * @param val Coefficient value.
 */
template <class NodeType, class CoefType, class Hash>
void insert_model(std::unordered_map<NodeType, CoefType, Hash> &um,
                  const NodeType &key, const CoefType &val) {
  if (um.count(key) == 0) {

    um.insert({{key, val}});

  } else {

    um[key] = val;
  }
}

}; // namespace helpers

#endif