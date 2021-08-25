/*! file one-solver-exhaustive.cpp
 *
 * @copyright Licensed under MIT license by hyperQ – Ewa Hendzel
 *
 */

#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/binary_object.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/pool/object_pool.hpp>
#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <ctime>
#include <vector>

#include <mpi.h>
#include "CL/sycl.hpp"
#include "helpers/devices.hpp"
#include "exhaustive/exhaustive.hpp"
#include "helpers/qubo_helpers.hpp"
#include "model/qubo.hpp"
#include "model/solution.hpp"

#define MAX_BUFFER_SIZE 12000

namespace po = boost::program_options;
namespace sycl = cl::sycl;

using queue_ptr = std::unique_ptr<sycl::queue>;


qubo::Solution sycl_native(qubo::QUBOModel<int, double> instance, std::string device_type)
{
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

    return solution;

}

int main(int argc, char *argv[]) {
  const clock_t begin_time = clock();

  char machine_name[MPI_MAX_PROCESSOR_NAME];
  int name_len=0;
  int rank=0;
  int num_procs=0;
  int size=0;

  std::string input_file;
  std::string output_file;
  std::string device_type;

  try {

    //MPI::Init(argc, argv);
    // Start MPI.
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
      std::cout << "Failed to initialize MPI\n";
      exit(-1);
    }

    // Create the communicator, and retrieve the number of MPI ranks.
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // Determine the rank number.
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Get the machine name.
    MPI_Get_processor_name(machine_name, &name_len);

    std::cout << "Rank #" << rank << " runs on: " << machine_name;
              // << ", uses device: "
              // << myQueue.get_device().get_info<info::device::name>() << "\n";

    if(rank == 0){
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
        MPI_Finalize();
        return 0;
      }

      if (!vm.count("input")) {
        std::cerr << "No input file provided." << std::endl;
        MPI_Finalize();
        return -1;
      }

      if (!vm.count("output")) {
        std::cerr << "No output file provided." << std::endl;
        MPI_Finalize();
        return -1;
      }

      if (device_type != "cpu" && device_type != "gpu" && device_type != "host") {
        std::cerr << "Unknown device type: " << device_type << std::endl;
        MPI_Finalize();
        return -1;
      }

      std::cout << "Reading input from: " << input_file << std::endl;

      std::cout << "Output will be saved to: " << output_file << std::endl;

      std::ifstream qubo_file(input_file);

      if (!qubo_file) {
        std::cerr << "can not open input file: " << input_file << std::endl;
        MPI_Finalize();
        return -1;
      }

      qubo_file.unsetf(std::ios::skipws);
      auto instance = qubo::QUBOModel<int, double>::load(qubo_file);

      std::ostringstream outstream;
      //boost::archive::binary_oarchive oa(stream);
      boost::archive::binary_oarchive ar(outstream, boost::archive::no_header);
      ar << boost::serialization::make_binary_object(&instance, sizeof(instance));

      std::cout << "Size: " << outstream.str().size() << "\n";

      //MPI_Send(stream.str().c_str(), stream.str().size(), MPI_BYTE, 0, 1, MPI_COMM_WORLD);

      MPI::COMM_WORLD.Send(outstream.str().c_str(),
                           outstream.str().size(),
                           MPI::Datatype(MPI_BYTE),
                           1,
                           1);

      auto solution = sycl_native(instance, device_type);

      std::ofstream results_file(output_file);

      solution.save(results_file);

      results_file.close();

    }
    else{

      std::vector<char> incomingBuffer(MAX_BUFFER_SIZE);
      MPI::Status msgStatus; 

      MPI::COMM_WORLD.Recv(&incomingBuffer[0], incomingBuffer.size(),
                           MPI::Datatype(MPI_BYTE),
                           0,
                           1, msgStatus);

      incomingBuffer.resize(msgStatus.Get_count(MPI::Datatype(MPI_BYTE)));

      std::cout << "msg size: " << incomingBuffer.size() << "\n";
      //std::cout << "error: " << msgStatus.Get_error() << "\n";

      if (incomingBuffer.size() > 0) {
      std::istringstream instream(std::string(&incomingBuffer[0], incomingBuffer.size()));
      //std::istringstream instream(std::string("hi"));
      boost::archive::binary_iarchive ar(instream);

      //boost::object_pool<qubo::QUBOModel<int, double>> pool;

      // qubo::QUBOModel<int, double> *instance;
      // {
      //   //boost::archive::binary_iarchive ia(instream);
      // }

      qubo::QUBOModel<int, double> *instance;
      ar >> boost::serialization::make_binary_object(instance, sizeof(instance));

      //auto solution = sycl_native(*instance, device_type);

      }

      


    }
    

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