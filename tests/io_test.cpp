/*! file io_test.cpp
 *
 * @copyright Licensed under MIT license by hyperQ – Ewa Hendzel
 *
 */
#include "model/qubo.hpp"
#include "model/solution.hpp"
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <sstream>

const std::string
    qubo_file_contents("c qubo Target MaxNodes NumNodes NumLinks\n"
                       "p qubo 0 100 2 2\n"
                       "c comment\n"
                       "0 0 -0.5\n"
                       "0 1 2.0\n"
                       "1 2 4\n"
                       "2 2 -0.7\n");

const std::string grammatically_incorrect_file_contents[]{
    std::string("p qubo 0 1 100 3 492\n"
                "0 0 -0.5\n"
                "0 1 2.0\n"),
    std::string("p qubo 0 100 3 492\n"
                "0 0 1 -0.5\n"
                "0 1 2.0\n"),
    std::string("p qubo 0 100 3 492\n"
                "0 0 -0.5\n"
                "0 2.0\n"),
    std::string("p qubo 0 100 492\n"
                "0 0 -0.5\n"
                "0 1 2.0\n"),
    std::string("p qubo 0 100 3 492\n"
                "0 0 -0.5\n"
                "unexpected string\n"
                "0 1 2.0\n")};

BOOST_AUTO_TEST_SUITE(io_test)

BOOST_AUTO_TEST_CASE(number_of_variable_agrees_with_given_in_input_stream) {
  std::stringstream stream(qubo_file_contents);
  stream.unsetf(std::ios::skipws);
  qubo::QUBOModel<int, double> problem =
      qubo::QUBOModel<int, double>::load(stream);

  BOOST_TEST_REQUIRE(problem.get_nodes() == 3);
}

BOOST_AUTO_TEST_CASE(all_linear_coefficients_are_present_in_the_model) {
  std::stringstream stream(qubo_file_contents);
  stream.unsetf(std::ios::skipws);
  qubo::QUBOModel<int, double> problem =
      qubo::QUBOModel<int, double>::load(stream);

  BOOST_TEST_REQUIRE(problem.get_variable(0) == -0.5);
  BOOST_TEST_REQUIRE(problem.get_variable(2) == -0.7);
}

BOOST_AUTO_TEST_CASE(all_quadratic_coefficients_are_present_in_the_model) {
  std::stringstream stream(qubo_file_contents);
  stream.unsetf(std::ios::skipws);
  qubo::QUBOModel<int, double> problem =
      qubo::QUBOModel<int, double>::load(stream);

  BOOST_TEST_REQUIRE(problem.get_connection(std::pair(0, 1)) == 2.0);
  BOOST_TEST_REQUIRE(problem.get_connection(std::pair(1, 2)) == 4.0);
}

BOOST_AUTO_TEST_CASE(
    read_qubo_throws_if_interaction_is_placed_in_reversed_order) {
  std::string invalid_qubo_file_contents("p qubo 0 100 2 2\n"
                                         "0 0 -0.5\n"
                                         "1 0 2.0\n"
                                         "1 2 4\n"
                                         "2 2 -0.7\n");

  std::stringstream stream(invalid_qubo_file_contents);
  stream.unsetf(std::ios::skipws);

  BOOST_CHECK_THROW((qubo::QUBOModel<int, double>::load(stream)),
                    std::invalid_argument);
}

BOOST_DATA_TEST_CASE(read_qubo_throws_if_stream_contents_dont_adhere_to_grammar,
                     grammatically_incorrect_file_contents) {
  std::stringstream stream(sample);
  stream.unsetf(std::ios::skipws);

  BOOST_CHECK_THROW((qubo::QUBOModel<int, double>::load(stream)),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(number_of_quadratic_terms_is_equal_to_the_declared_one) {
  std::string invalid_qubo_file_contents("p qubo 0 100 3 3\n"
                                         "0 0 -0.5\n"
                                         "0 1 2.0\n"
                                         "1 2 4\n"
                                         "2 2 -0.7\n");

  std::stringstream stream(invalid_qubo_file_contents);
  stream.unsetf(std::ios::skipws);

  BOOST_CHECK_THROW((qubo::QUBOModel<int, double>::load(stream)),
                    std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(saved_solution_has_correct_structure) {
  char state[] = {0, 1, 1, 0, 1};
  qubo::Solution solution(state, state + 5, -12.5);

  std::string expected_output = "0,1,2,3,4,energy\n0,1,1,0,1,-12.5\n";

  std::stringstream stream("");
  solution.save(stream);

  BOOST_TEST_REQUIRE(stream.str() == expected_output);
}

BOOST_AUTO_TEST_SUITE_END()
