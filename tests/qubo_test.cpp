/*! file qubo_test.cpp
 *
 * @copyright Licensed under MIT license by hyperQ – Ewa Hendzel
 *
 */

#define BOOST_TEST_MODULE qubo_test_module

#include "model/qubo.hpp"
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(qubo_test)

BOOST_AUTO_TEST_CASE(qubo_linear_variable) {

  qubo::LinearCoef<int, double> linear_c{
      {1, -0.4}, {2, 1.1}, {3, -1.0}, {4, -1.2}};

  qubo::QuadraticCoef<int, double> quadratic_c{
      {std::make_pair(1, 2), 1.0}, {std::make_pair(1, 3), 7.2},
      {std::make_pair(1, 4), 1.0}, {std::make_pair(2, 3), 2.0},
      {std::make_pair(2, 4), 1.9}, {std::make_pair(3, 4), 3.0}};

  qubo::QUBOModel<int, double> qubos(linear_c, quadratic_c);
  qubos.add_variable(1, 20);

  BOOST_TEST_REQUIRE(qubos.get_variable(1) == 20);
}

BOOST_AUTO_TEST_CASE(qubo_quadratic_variable) {

  qubo::LinearCoef<int, double> linear_c{};

  qubo::QuadraticCoef<int, double> quadratic_c{};

  qubo::QUBOModel<int, double> qubos(linear_c, quadratic_c);
  qubos.add_connection(std::make_pair(1, 2), 20);

  BOOST_TEST_REQUIRE(qubos.get_connection(std::make_pair(1, 2)) == 20);
}

BOOST_AUTO_TEST_CASE(qubo_print) {

  qubo::LinearCoef<int, int> linear_c{};

  qubo::QuadraticCoef<int, int> quadratic_c{};

  qubo::QUBOModel<int, int> qubos(linear_c, quadratic_c);
  qubos.add_variable(1, 10);
  qubos.add_connection(std::make_pair(1, 2), 20);

  auto s = qubos.str();

  BOOST_TEST_REQUIRE(s.compare("QUBO model 1--1:10 1--2:20") == 0);
}

BOOST_AUTO_TEST_SUITE_END()
