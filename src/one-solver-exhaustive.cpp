/*! file one-solver-exhaustive.cpp
 *
 * @copyright Licensed under MIT license by hyperQ – Ewa Hendzel
 *
 */

// include headers that implement a archive in simple text format
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
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
  int count =0;
  long long int msg_size =0;

  std::string input_file;
  std::string output_file;
  std::string device_type;
  std::vector<char> msg_buff;

  MPI::Status msg_status;

  qubo::QUBOModel<int, double> instance;
  

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

    //std::cout << "Rank #" << rank << " runs on: " << machine_name;
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
      instance = qubo::QUBOModel<int, double>::load(qubo_file);

      std::ostringstream oss;
      // save data to archive
      {
         boost::archive::text_oarchive ar(oss);
         // write class instance to archive
         ar << instance;
         // archive and stream closed when destructors are called
      }

      //std::cout << "sSize: " << oss.str().size() << std::endl;
      //std::cout << "sData: " << oss.str().c_str() << std::endl;

      msg_buff.assign(oss.str().c_str(), oss.str().c_str() + oss.str().size() +1);
      msg_size = oss.str().size();
    }

    MPI::COMM_WORLD.Bcast(&msg_size,
                         1,
                         MPI::LONG_LONG_INT,
                         0);
  
    if(rank != 0){
      msg_buff.resize(msg_size);
    }

    MPI::COMM_WORLD.Bcast(msg_buff.data(),
                         msg_size,
                         MPI::Datatype(MPI::CHAR),
                         0);
    
    if(rank != 0){
      if (msg_buff.size() > 0) {
      std::istringstream iss(std::string(msg_buff.data(), msg_buff.size()));
      //std::cout << "rSize: " << iss.str().size() << std::endl;
      //std::cout << "rData: " << iss.str().c_str() << std::endl;

       {
         boost::archive::text_iarchive ar(iss);   
         ar >> instance;
       }
        //auto solution = sycl_native(instance, device_type);
     }
    }

    if(rank == 0){
      msg_size = device_type.size();
      msg_buff.assign(device_type.c_str(),device_type.c_str() + device_type.size() +1);
    }
  
    MPI::COMM_WORLD.Bcast(&msg_size,
                         1,
                         MPI::LONG_LONG_INT,
                         0);

    if(rank != 0){
      msg_buff.resize(msg_size);
    }
  
    MPI::COMM_WORLD.Bcast(msg_buff.data(),
                         msg_size,
                         MPI::Datatype(MPI::CHAR),
                         0);

    if(rank != 0){
      device_type = std::string(&msg_buff[0], msg_size);
      //std::cout << "device_type: " << device_type << std::endl;
    }

    auto solution = sycl_native(instance, device_type);

    std::cout << "Solution: " << solution.energy << std::endl;

    if(rank == 0){
      std::ofstream results_file(output_file);
      solution.save(results_file);
      results_file.close();
    }

    
  } catch (std::exception &e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  } catch (...) {
    std::cerr << "Exception of unknown type!\n";
  }
  if(rank == 0){
    std::cout << "Calculation time [s]: "<< float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\n";
  }
  
  MPI_Finalize();
  return 0;
}