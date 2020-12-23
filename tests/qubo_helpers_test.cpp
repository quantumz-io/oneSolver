/*! file qubo_helpers_test.cpp
 *
 * @copyright Licensed under MIT license by hyperQ – Ewa Hendzel
 *
 */
#include "helpers/qubo_helpers.hpp"
#include "model/qubo.hpp"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(qubo_helpers_test)

BOOST_AUTO_TEST_CASE(qubo_flattening) {

  qubo::LinearCoef<int, double> linear_c{
      {0, 0.5}, {1, -2.0}, {2, 1.0}, {4, -1.5}};

  qubo::QuadraticCoef<int, double> quadratic_c{
      {std::pair(0, 1), 1.0}, {std::pair(0, 3), 7.2},  {std::pair(1, 4), -1.0},
      {std::pair(2, 3), 2.0}, {std::pair(2, 4), -1.5}, {std::pair(3, 4), -3.5}};
  qubo::QUBOModel<int, double> instance(linear_c, quadratic_c);
  instance.set_nodes(5);

  auto flat_qubo = helpers::flatten_qubo(instance);

  std::vector<double> expected_result(
      {0.5, 1.0,  0.0, 7.2, 0.0, 1.0, -2.0, 0.0, 0.0,  -1.0, 0.0,  0.0, 1.0,
       2.0, -1.5, 7.2, 0.0, 2.0, 0.0, -3.5, 0.0, -1.0, -1.5, -3.5, -1.5});

  BOOST_TEST(flat_qubo == expected_result);
}

BOOST_AUTO_TEST_SUITE_END()