/*! file solution.hpp
 *
 * @copyright Licensed under MIT license by hyperQ – Ewa Hendzel
 *
 */
#include <vector>

#ifndef SOLUTION_HPP_
#define SOLUTION_HPP_

/*!
 * @brief Classes representing QUBO Model and its solutions.
 */
namespace qubo {

/*!
 * @brief Solution of QUBO model.
 */
class Solution {
public:
  /*! State, i.e. configuration of the system. i-th element of the
      vector corresponds to i-th variable.
  */
  std::vector<char> state;
  /*! Energy of the state stored in this solution. */
  double energy;

  /*!
   * @brief Create a solution, taking values of variables from
   *        some iterator.
   *
   * An iterator pair provided to this constructor should define
   * a range of N 0-1 values, where i-th value corresponds to
   * i-th variable. This implicitly assumes 0-based indexing
   * of variables.
   *
   * @param begin start of the variables' configuration.
   * @param end end of the variables' configuration.
   * @param energy energy of the configuration.
   */
  template <typename InputIt>
  Solution(InputIt begin, InputIt end, double energy)
      : state(begin, end), energy(energy){};

  /*!
   * @brief serialize this solution to the given output stream.
   *
   * The solution is formatted in two-row CSV format. Both lines
   * are N-items long, where N is the number of variables.
   * The first line serves as a header, and contains comma-separated
   * indices of variables followed by string "energy".
   * The second line contains comma separated values of variables
   * followed by an energy of the configuration.
   *
   * @param stream output stream to which this soulution should be
   *               serialized.
   */
  void save(std::ostream &stream) {
    for (auto var_index = 0; var_index < state.size(); var_index++) {
      stream << var_index << ',';
    }
    stream << "energy" << std::endl;
    for (auto &bit : state) {
      stream << int(bit) << ",";
    }
    stream << energy << std::endl;
  }
};
} // namespace qubo

#endif