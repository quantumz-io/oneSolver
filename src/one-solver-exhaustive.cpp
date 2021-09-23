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

namespace po = boost::program_options;
namespace sycl = cl::sycl;

using queue_ptr = std::unique_ptr<sycl::queue>;

/*!
 * @brief Compute QUBO solution in a give solution space on a given device.
 * @param Instance of the proble to solve.
 * @param Accelerator device type.
 * @param Start state in the solution space.
 * @param End state in the solution space.
 *
 * @returns solution to the given QUBO problem and acclerator device name.
 */
std::tuple<qubo::Solution, std::string> sycl_native(qubo::QUBOModel<int, double> instance, std::string device_type, ulong start_state, ulong end_state)
{
  queue_ptr q_ptr;
  std::string accl_device_name;
  
  try {
    q_ptr = queue_ptr(
        new sycl::queue(*devices::construct_device_selector(device_type)));
  } catch (std::runtime_error) {
    std::cerr << "No devices of given type could be initialized."
              << std::endl;
  }
 
  accl_device_name = q_ptr->get_device().get_info<sycl::info::device::name>();
  // std::cout << "Using device: "
  //           << q_ptr->get_device().get_info<sycl::info::device::name>()
  //           << std::endl;

  auto solution = exhaustive::solve(*q_ptr, instance, start_state, end_state);

  return {solution, accl_device_name} ;
}

/*!
 * @brief Parse command line.
 * @param argument count.
 * @param arguments array.
 *
 * @returns returns input file name, output file name and device type.
 */
std::tuple<std::string, std::string, std::string> parse_cmd_line(int argc, char *argv[])
{
  std::string input_file;
  std::string output_file;
  std::string device_type;

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
    MPI_Abort(MPI_COMM_WORLD, MPI_SUCCESS);
  }

  if (!vm.count("input")) {
    throw std::runtime_error("No input file provided.");
  }

  if (!vm.count("output")) {
    throw std::runtime_error("No output file provided.");
  }
  
  if (device_type != "cpu" && device_type != "gpu" && device_type != "host") {
    throw std::runtime_error("Unknown device type: " + device_type);
  }

  return {input_file, output_file, device_type};
}

int main(int argc, char *argv[]) {
  const clock_t begin_time = clock();

  char machine_name[MPI_MAX_PROCESSOR_NAME];
  int name_len = 0;
  int rank = 0;
  const int root_rank = 0; // Rank zero process.
  int min_energy_process_rank = 0; // Rank of process with minimum energy. 
  int num_procs = 0;

  long long int msg_size = 0;

  ulong start_state = 0;
  ulong end_state = 0;

  ulong *ranges_buf = NULL;
  ulong recv_arr[2];

  std::string input_file;
  std::string output_file;
  std::string device_type;
  std::vector<char> msg_buff;
  
  qubo::QUBOModel<int, double> instance;
  
  try {
    // Start MPI.
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
      throw std::runtime_error("Failed to initialize MPI\n");
    }

    // Create the communicator, and retrieve the number of MPI ranks.
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // Determine the rank number.
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Get the machine name.
    MPI_Get_processor_name(machine_name, &name_len);

    if (rank == root_rank) {

      std::tie(input_file, output_file, device_type) = parse_cmd_line(argc, argv);

      std::cout << "Reading input from: " << input_file << std::endl;

      std::cout << "Output will be saved to: " << output_file << std::endl;

      std::ifstream qubo_file(input_file);
      
      if (!qubo_file) {
        throw std::runtime_error("can not open input file:" + input_file);
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
#ifdef DEBUG
      std::cout << "sSize: " << oss.str().size() << std::endl;
      std::cout << "sData: " << oss.str().c_str() << std::endl;
#endif
      msg_buff.assign(oss.str().c_str(), oss.str().c_str() + oss.str().size() +1);
      msg_size = oss.str().size();
    }

    // Send QUBOModel instance to all the proecesses.
    if (MPI_Bcast(&msg_size, 1, MPI_LONG_LONG_INT, root_rank, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("MPI broadcast failed\n");
    }

    if (rank != root_rank) {
      msg_buff.resize(msg_size);
    }

    if (MPI_Bcast(msg_buff.data(), msg_size, MPI_CHAR, root_rank, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("MPI broadcast failed\n");
    }
    
    if (rank != root_rank) {
      if (msg_buff.size() > 0) {
      std::istringstream iss(std::string(msg_buff.data(), msg_buff.size()));
#ifdef DEBUG
      std::cout << "rSize: " << iss.str().size() << std::endl;
      std::cout << "rData: " << iss.str().c_str() << std::endl;
#endif
       {
         boost::archive::text_iarchive ar(iss);   
         ar >> instance;
       }
     }
    }

    // Send device_type to all the processes.
    if (rank == root_rank) {
      msg_size = device_type.size();
      msg_buff.assign(device_type.c_str(),device_type.c_str() + device_type.size() +1);
    }
  
    if (MPI_Bcast(&msg_size, 1, MPI_LONG_LONG_INT, root_rank, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("MPI broadcast failed\n");
    }

    if (rank != root_rank) {
      msg_buff.resize(msg_size);
    }
  
    if (MPI_Bcast(msg_buff.data(), msg_size, MPI_CHAR, root_rank, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("MPI broadcast failed\n");
    }

    if (rank != root_rank) {
      device_type = std::string(&msg_buff[0], msg_size);
#ifdef DEBUG
      std::cout << "device_type: " << device_type << std::endl;
#endif
    }

    if (rank == root_rank) {
      // Divide task to each node and let sycl decide the distribution of task on each node.
      auto n_bits = instance.get_nodes();
      ranges_buf = new ulong[num_procs*2];
      // note that all the variables related to states have to be unsigned 64 bit integers
      ulong n_states = 1 << n_bits; 
      ulong states_per_node = n_states / num_procs;
      ulong states_per_node_remainder = n_states % num_procs;

      ulong reminder_count = 0;
      ulong states_count = 0;

      for (auto i = 0; i < num_procs; ++i) {
        start_state = states_count;
        states_count += states_per_node;
        if (reminder_count < states_per_node_remainder) {
          states_count += 1;
          reminder_count++;
        }
        end_state = states_count;

        ranges_buf[i*2] = start_state;
        ranges_buf[i*2+1] = end_state;
#ifdef DEBUG
        std::cout << ranges_buf[i*2];
        std::cout << ", ";
        std::cout << ranges_buf[i*2+1] << std::endl;
#endif
      }
    }

    if (MPI_Scatter(ranges_buf, 2, MPI_UNSIGNED_LONG, recv_arr, 2, MPI_UNSIGNED_LONG, root_rank, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("MPI scatter failed\n");
    }

    start_state = recv_arr[0];
    end_state = recv_arr[1];

    auto [solution, accl_device_name] = sycl_native(instance, device_type, start_state, end_state);

    if(rank == root_rank){
      //std::unique_ptr<char> log_buff_ptr(new char[100]);
      MPI_Status status;
      int msg_length;
      
      //char log_buff_ptr[100];
      std::cout << "Node ID [" << machine_name << "], " << "Rank [" << rank << "], " << "Using device: " << accl_device_name << std::endl;
      for(auto rank_indx = 1; rank_indx < num_procs; ++rank_indx){
        MPI_Probe(rank_indx, 1, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_CHAR, &msg_length);

        std::unique_ptr<char> log_buff_ptr(new char[msg_length]{});
        std::cout << log_buff_ptr.get();
        //std::cout << "rank_indx: " << rank_indx << std::endl;
        MPI_Recv(log_buff_ptr.get(), msg_length, MPI_CHAR, rank_indx, 1 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::string log_msg(log_buff_ptr.get());
        //std::cout << log_buff_ptr.get() << std::endl;
        std::cout << log_msg;
      }
    }else{
      std::ostringstream log_stream;
      log_stream << "Node ID [" << machine_name << "], " << "Rank [" << rank << "], " << "Using device: " << accl_device_name << std::endl;
      MPI_Send(log_stream.str().c_str(), log_stream.str().size(), MPI_CHAR, root_rank, 1, MPI_COMM_WORLD);
    }
    

#ifdef DEBUG
    std::cout << "Solution: " << solution.energy << std::endl;
#endif
    
    std::unique_ptr<double> energy_buff_ptr(new double[num_procs]);

    if (MPI_Gather(&solution.energy, 1, MPI_DOUBLE, energy_buff_ptr.get(), 1, MPI_DOUBLE, root_rank, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("MPI gather failed\n");
    }

    if (rank == root_rank) {
       std::vector<double> energies(&energy_buff_ptr.get()[0], &energy_buff_ptr.get()[num_procs]);
       auto min_energy = std::min_element(energies.begin(), energies.end());
       auto min_idx_energy = std::distance(energies.begin(), min_energy);
#ifdef DEBUG
       std::cout << "min_energy: " << energies[min_idx_energy] << std::endl;
#endif
       solution.energy = energies[min_idx_energy];

       min_energy_process_rank = min_idx_energy;
    }

    if (MPI_Bcast(&min_energy_process_rank, 1, MPI_INT, root_rank, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("MPI broadcast failed\n");
    }
    
    if ((min_energy_process_rank != root_rank) && (rank == min_energy_process_rank)) {
      auto state = solution.state; 
      if (MPI_Send(state.data(), state.size(), MPI_CHAR, root_rank, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
        throw std::runtime_error("MPI send failed\n");
      }
    }

    if (rank == root_rank) {
      char buff[32];
      if (min_energy_process_rank != root_rank) {
        if (MPI_Recv(&buff, 32, MPI_CHAR, min_energy_process_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
          throw std::runtime_error("MPI receive failed\n");
        }
#ifdef DEBUG
        std::cout << "Min Energy State: "       
        for (auto i = 0; i < solution.state.size(); ++i) {
        std::cout << (int)buff[i] << " ";
        }
        std::cout << std::endl;
#endif
        for (auto i = 0; i < solution.state.size(); ++i ) {
          solution.state[i] = buff[i];
        } 
      }

#ifdef DEBUG
      std::cout << "Min Energy State: "
      for (auto i = 0; i < solution.state.size(); ++i) {
        std::cout << (int)solution.state[i] << " ";
      }
      std::cout << std::endl;
#endif
      std::ofstream results_file(output_file);
      solution.save(results_file);
      results_file.close();
    }
      
  } catch (std::exception &e) {
    std::cerr << "error: " << e.what() << "\n";
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  } catch (...) {
    std::cerr << "Exception of unknown type!\n";
  }

  if (rank == root_rank) {
    delete [] ranges_buf;
    std::cout << "Calculation time [s]: "<< float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\n";
  }

  MPI_Finalize();
  return 0;
}