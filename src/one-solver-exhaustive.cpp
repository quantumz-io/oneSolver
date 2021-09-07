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


qubo::Solution sycl_native(qubo::QUBOModel<int, double> instance, std::string device_type, ulong start_state, ulong end_state)
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

    auto solution = exhaustive::solve(*q_ptr, instance, start_state, end_state);

    return solution;

}

int main(int argc, char *argv[]) {
  const clock_t begin_time = clock();

  char machine_name[MPI_MAX_PROCESSOR_NAME];
  int name_len=0;
  int rank=0;
  int process_rank = 0;
  int num_procs=0;
  int size=0;
  int count =0;
  long long int msg_size =0;

  ulong start_state = 0;
  ulong end_state = 0;

  ulong *ranges_buf = NULL;
  ulong recv_arr[2];

  std::string input_file;
  std::string output_file;
  std::string device_type;
  std::vector<char> msg_buff;
  

  

  //MPI::Status msg_status;

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

    if(rank == 0){
      // Divide task to each node and let sycl decide the distribution of task on each node.
      auto n_bits = instance.get_nodes();
      //std::vector<std::vector<ulong>> ranges_buf(num_procs, std::vector<ulong> (2,0));
      //ulong ** ranges_buf = new ulong[num_procs][2];
      //auto ranges_buf = new ulong [num_procs][2];
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
      // ulong* p =  &ranges_buf[0];
      // for (auto i = 0; i < num_procs*2; ++i) {
       
      //  std::cout << p << ":" << *p << std::endl;
      //  ++p;
      //  //std::cout << p << ":" << *p << std::endl;
      // }
    }

    MPI::COMM_WORLD.Scatter(ranges_buf, 2, MPI::UNSIGNED_LONG,recv_arr, 2, MPI::UNSIGNED_LONG, 0);
    //std::cout << "Rank #" << rank << " range: " << recv_arr[0] << "-" << recv_arr[1] << std::endl;
    start_state = recv_arr[0];
    end_state = recv_arr[1];


    auto solution = sycl_native(instance, device_type, start_state, end_state);

    //std::cout << "Solution: " << solution.energy << std::endl;

    
    double *energy_buff = new double[num_procs];

    MPI::COMM_WORLD.Gather(&solution.energy, 1, MPI::DOUBLE, energy_buff, 1, MPI::DOUBLE, 0);
    
    if(rank == 0){
       std::vector<double> energies(&energy_buff[0], &energy_buff[num_procs]);
       auto min_energy = std::min_element(energies.begin(), energies.end());
       auto min_idx_energy = std::distance(energies.begin(), min_energy);
       //std::cout << "min_idx: " << energies[min_idx_energy] << std::endl;
       solution.energy = energies[min_idx_energy];

       process_rank = min_idx_energy;
       //std::cout << min_energy << std::endl;
       
      //  for (auto i = energies.begin(); i != energies.end(); ++i)
      //     std::cout << *i << ' ';

      //  for (auto i = 0; i < num_procs; ++i) {
      //  std::cout << "energy: " << ":" <<  energy_buff[i] << std::endl;
      //  }
    }

    MPI::COMM_WORLD.Bcast(&process_rank, 1, MPI::INT, 0);

    
    if ((process_rank !=0) && (rank == process_rank))
    {
       auto state = solution.state; 
       MPI::COMM_WORLD.Send(state.data(), state.size(), MPI_CHAR, 0, 0);
       //std::cout << "rank: " << rank << std::endl;
    }

    if(rank == 0){
      char buff[32];
      if(process_rank != 0){
        MPI::COMM_WORLD.Recv(&buff, 32, MPI_CHAR, process_rank, 0 );
        // for (auto i = 0; i < solution.state.size(); ++i) {
        // std::cout << (int)buff[i] << " ";
        // }
        // std::cout << std::endl;
        for(auto i=0; i < solution.state.size(); ++i ){
          solution.state[i] = buff[i];
        } 
      }
      // for (auto i = 0; i < solution.state.size(); ++i) {
      //   std::cout << (int)solution.state[i] << " ";
      // }
      std::cout << std::endl;
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
    delete ranges_buf;
    std::cout << "Calculation time [s]: "<< float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\n";
  }
  
  MPI_Finalize();
  return 0;
}