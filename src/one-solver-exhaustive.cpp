/*! file one-solver-exhaustive.cpp
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
#include "exhaustive/exhaustive.hpp"
#include "helpers/qubo_helpers.hpp"
#include "model/qubo.hpp"
#include "model/solution.hpp"

namespace po = boost::program_options;
namespace sycl = cl::sycl;

using queue_ptr = std::unique_ptr<sycl::queue>;

int main(int argc, char *argv[]) {
  const clock_t begin_time = clock();

  std::string input_file;
  std::string output_file;
  std::string device_type;

  try {

    po::options_description options("Allowed options");
    options.add_options()("help", "produce help message")(
        "input", po::value<std::string>(&input_file), "input file")(
        "output", po::value<std::string>(&output_file), "output file")(
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

    std::cout << "Reading input from: " << input_file << std::endl;

    std::cout << "Output will be saved to: " << output_file << std::endl;

    std::ifstream qubo_file(input_file);

    if (!qubo_file) {
      std::cerr << "can not open input file: " << input_file << std::endl;
      return -1;
    }
    
    // Start MPI.
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
      std::cout << "Failed to initialize MPI\n";
      exit(-1);
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

    auto solution = exhaustive::solve(*q_ptr, instance);

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
  MPI_Finalize();
  return 0;
}