/*! file qubo_helpers.hpp
 *
 * @copyright Licensed under MIT license by hyperQ – Ewa Hendzel
 *
 */
#include "model/qubo.hpp"
#include <vector>

#ifndef QUBO_HELPERS_HPP_
#define QUBO_HELPERS_HPP_

/*!
 * @brief Helper functions and utilities for use with QUBO models.
 */
namespace helpers {
template <class CoeffType>

/*!
 * @brief Convert QUBO to a flat vector with problem coefficients.
 *
 * The flattened QUBO is stored in dense vector in a row-wise fashion.
 *
 * @param instance a problem to be flattened.
 * @returns Flattened qubo.
 */
std::vector<CoeffType> flatten_qubo(qubo::QUBOModel<int, CoeffType> instance) {
  auto num_variables = instance.get_nodes();

  std::vector<CoeffType> result(num_variables * num_variables);

  for (auto i = 0; i < num_variables; i++) {
    for (auto j = 0; j < num_variables; j++) {
      if (i == j) {
        result.at(i + j * num_variables) = instance.get_variable(i);
      } else {
        auto coefficient = instance.get_connection(std::pair(i, j));
        result.at(i + j * num_variables) += coefficient;
        result.at(j + i * num_variables) += coefficient;
      }
    }
  }

  return result;
}
} // namespace helpers

#endif