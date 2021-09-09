/*! file exhaustive_test.cpp
 *
 * @copyright Licensed under MIT license by hyperQ – Ewa Hendzel
 *
 */
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "exhaustive/exhaustive.hpp"
#include "model/qubo.hpp"

const std::pair<std::string, double> qubo_file_contents[]{
    std::make_pair<std::string, double>(
        "c qubo Target MaxNodes NumNodes NumLinks\n"
        "p qubo 0 4 4 6\n"
        "c comment\n"
        "0 0 -5\n"
        "0 1 2\n"
        "0 2 4\n"
        "0 3 0\n"
        "1 1 -3\n"
        "1 2 1\n"
        "1 3 0\n"
        "2 2 -8\n"
        "2 3 5\n"
        "3 3 -6\n",
        -12.0),
    std::make_pair<std::string, double>(
        "c qubo Target MaxNodes NumNodes NumLinks\n"
        "p qubo 0 3 3 1\n"
        "c comment\n"
        "0 0 -0.5\n"
        "1 1 -0.5\n"
        "1 2 4\n"
        "2 2 -0.7\n",
        -1.2),
    std::make_pair<std::string, double>(
        "c qubo Target MaxNodes NumNodes NumLinks\n"
        "p qubo 0 5 5 10\n"
        "0 0 -6.000000\n"
        "0 1 -4.000000\n"
        "0 2 12.000000\n"
        "0 3 16.000000\n"
        "0 4 -12.000000\n"
        "1 1 12.000000\n"
        "1 2 -8.000000\n"
        "1 3 -20.000000\n"
        "1 4 8.000000\n"
        "2 2 -8.000000\n"
        "2 3 4.000000\n"
        "2 4 8.000000\n"
        "3 3 -2.000000\n"
        "3 4 4.000000\n"
        "4 4 -4.000000\n",
        -22.0),
    std::make_pair<std::string, double>(
        "c qubo Target MaxNodes NumNodes NumLinks\n"
        "p qubo 0 7 5 20\n"
        "0 0 6.000000\n"
        "0 1 -16.000000\n"
        "0 2 8.000000\n"
        "0 3 0.000000\n"
        "0 4 -4.000000\n"
        "0 5 -12.000000\n"
        "0 6 12.000000\n"
        "1 2 4.000000\n"
        "1 3 0.000000\n"
        "1 4 8.000000\n"
        "1 5 0.000000\n"
        "1 6 4.000000\n"
        "2 2 -4.000000\n"
        "2 4 4.000000\n"
        "2 5 -16.000000\n"
        "2 6 8.000000\n"
        "3 3 14.000000\n"
        "3 4 -16.000000\n"
        "3 5 -4.000000\n"
        "3 6 -8.000000\n"
        "4 5 4.000000\n"
        "4 6 4.000000\n"
        "5 5 16.000000\n"
        "5 6 -4.000000\n"
        "6 6 -8.000000\n",
        -14.0),
    std::make_pair<std::string, double>(
        "c qubo Target MaxNodes NumNodes NumLinks\n"
        "p qubo 0 13 11 42\n"
        "0 0 -2.000000\n"
        "0 1 -4.000000\n"
        "0 2 -4.000000\n"
        "0 4 -4.000000\n"
        "0 5 8.000000\n"
        "0 6 4.000000\n"
        "0 9 0.000000\n"
        "0 10 4.000000\n"
        "1 1 -6.000000\n"
        "1 2 4.000000\n"
        "1 3 8.000000\n"
        "1 5 -4.000000\n"
        "1 6 -8.000000\n"
        "1 7 0.000000\n"
        "1 9 12.000000\n"
        "1 12 4.000000\n"
        "2 2 -12.000000\n"
        "2 3 4.000000\n"
        "2 5 4.000000\n"
        "2 7 4.000000\n"
        "2 9 8.000000\n"
        "2 10 4.000000\n"
        "2 12 0.000000\n"
        "3 3 -8.000000\n"
        "3 5 4.000000\n"
        "3 6 -8.000000\n"
        "3 7 4.000000\n"
        "3 9 4.000000\n"
        "3 10 8.000000\n"
        "3 12 -8.000000\n"
        "4 4 4.000000\n"
        "4 5 -4.000000\n"
        "4 6 0.000000\n"
        "4 7 -4.000000\n"
        "4 12 4.000000\n"
        "5 5 -12.000000\n"
        "5 6 4.000000\n"
        "5 7 4.000000\n"
        "5 9 4.000000\n"
        "5 10 8.000000\n"
        "5 12 -4.000000\n"
        "6 6 8.000000\n"
        "6 7 -8.000000\n"
        "6 9 -4.000000\n"
        "6 12 4.000000\n"
        "7 7 2.000000\n"
        "7 12 -4.000000\n"
        "9 9 -16.000000\n"
        "9 10 4.000000\n"
        "9 12 4.000000\n"
        "10 10 -10.000000\n"
        "10 12 -8.000000\n"
        "12 12 4.000000\n",
        -32.0)
        };

BOOST_AUTO_TEST_SUITE(calculate_energy)

BOOST_AUTO_TEST_CASE(minimal_energy_properly_calculated) {
  for (auto &qubo_data : qubo_file_contents) {
    std::stringstream stream(qubo_data.first);
    stream.unsetf(std::ios::skipws);
    qubo::QUBOModel<int, double> problem =
        qubo::QUBOModel<int, double>::load(stream);

    sycl::queue q(sycl::cpu_selector{});
    auto solution = exhaustive::solve(q, problem, 0, 1 <<problem.get_nodes() );  //changed signature

    BOOST_TEST_REQUIRE(std::abs(solution.energy - qubo_data.second) < 1e-13);
  }
}

BOOST_AUTO_TEST_SUITE_END()