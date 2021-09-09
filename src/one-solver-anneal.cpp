/*! file one-solver-anneal.cpp
 *
 * @copyright Licensed under MIT license by hyperQ – Ewa Hendzel
 *
 */

#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <ctime>

#include "CL/sycl.hpp"
#include "helpers/devices.hpp"
#include "helpers/qubo_helpers.hpp"
#include "model/qubo.hpp"
#include "model/solution.hpp"
#include "simulated_annealing/annealing.hpp"

namespace po = boost::program_options;
namespace sycl = cl::sycl;

using queue_ptr = std::unique_ptr<sycl::queue>;

void construct_linear_beta_schedule(std::vector<double> &schedule,
                                    double beta_min, double beta_max,
                                    uint num_iter) {
  for (auto i = 0; i < num_iter; i++) {
    schedule[i] = beta_min + beta_max * i / static_cast<double>(num_iter - 1);
  }
}

void construct_geometric_beta_schedule(std::vector<double> &schedule,
                                       double beta_min, double beta_max,
                                       uint num_iter) {
  schedule[0] = beta_min;
  auto alpha = pow(beta_max / beta_min, 1.0 / (num_iter - 1));
  for (auto i = 1; i < num_iter; i++) {
    schedule[i] = (schedule[i - 1]) * alpha;
  }
}

int main(int argc, char *argv[]) {
  const clock_t begin_time = clock();

  std::string input_file;
  std::string output_file;
  std::string schedule_type;
  std::string device_type;

  uint num_iter;
  uint num_tries;

  double beta_min, beta_max;

  try {

    po::options_description options("Allowed options");
    options.add_options()("help", "produce help message")(
        "input", po::value<std::string>(&input_file), "input file")(
        "output", po::value<std::string>(&output_file), "output file")(
        "num-iter", po::value<uint>(&num_iter)->default_value(100),
        "number of iterations of the algorithm")(
        "num-tries", po::value<uint>(&num_tries)->default_value(100),
        "number of trajectories to try")(
        "schedule-type",
        po::value<std::string>(&schedule_type)->default_value("geometric"),
        "type of beta schedule tu use, either linear or geometric")(
        "beta-min", po::value<double>(&beta_min)->default_value(0.1),
        "minimum value of beta in the annealing schedule (default 0.1)")(
        "beta-max", po::value<double>(&beta_max)->default_value(1.0),
        "maximum value of beta in the annealing schedule (default 1.0)")(
        "device-type",
        po::value<std::string>(&device_type)->default_value("host"),
        "device type to use (cpu, gpu or host)");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, options), vm);
    po::notify(vm);

    if (vm.count("help")) {
      std::cout << options << std::endl;
      return 0;
    }

    if (!vm.count("input")) {
      std::cerr << "No input file provided." << std::endl;
      return -1;
    }

    if (!vm.count("output")) {
      std::cerr << "No output file provided." << std::endl;
      return -1;
    }

    if (device_type != "cpu" && device_type != "gpu" && device_type != "host") {
      std::cerr << "Unknown device type: " << device_type << std::endl;
      return -1;
    }

    if (schedule_type != "linear" && schedule_type != "geometric") {
      std::cerr << "Unknown beta schedule: " << schedule_type << std::endl;
      return -1;
    }

    if (beta_max < 0 || beta_min < 0) {
      std::cerr
          << "Invalid schedule, both ends of beta range need to be positive"
          << std::endl;
      return -1;
    }

    if (beta_min >= beta_max) {
      std::cerr
          << "Invalid schedule, initial beta is not lesser than final beta"
          << std::endl;
      return -1;
    }

    std::cout << "Reading input from: " << input_file << std::endl;

    std::cout << "Output will be saved to: " << output_file << std::endl;

    std::cout << "Schedule type: " << schedule_type << std::endl;

    std::cout << "Beta range: [" << beta_min << ", " << beta_max << "]"
              << std::endl;

    std::cout << "Number of iterations: " << num_iter << std::endl;
    std::cout << "Number of tries: " << num_tries << std::endl;

    std::ifstream qubo_file(input_file);

    if (!qubo_file) {
      std::cerr << "can not open input file: " << input_file << std::endl;
      return -1;
    }

    qubo_file.unsetf(std::ios::skipws);
    auto instance = qubo::QUBOModel<int, double>::load(qubo_file);

    queue_ptr q_ptr;

    try {
      q_ptr = queue_ptr(
          new sycl::queue(*devices::construct_device_selector(device_type)));
    } catch (std::runtime_error) {
      std::cerr << "No devices of given type could be initialized."
                << std::endl;
    }

    std::cout << "Using device: "
              << q_ptr->get_device().get_info<sycl::info::device::name>()
              << std::endl;

    std::vector<double> beta_schedule(num_iter);

    if (schedule_type == "linear") {
      construct_linear_beta_schedule(beta_schedule, beta_min, beta_max,
                                     num_iter);
    } else {
      construct_geometric_beta_schedule(beta_schedule, beta_min, beta_max,
                                        num_iter);
    }

    auto solution =
        sa::anneal(instance, *q_ptr, beta_schedule, num_iter, num_tries);

    std::ofstream results_file(output_file);

    solution.save(results_file);

    results_file.close();
  } catch (std::exception &e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  } catch (...) {
    std::cerr << "Exception of unknown type!\n";
  }

  std::cout << "Calculation time [s]: "<< float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\n";
  return 0;
}