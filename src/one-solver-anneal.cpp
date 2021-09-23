/*! file one-solver-anneal.cpp
 *
 * @copyright Licensed under MIT license by hyperQ – Ewa Hendzel
 *
 */

#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <ctime>
#include <random>

#include <mpi.h>
#include "CL/sycl.hpp"
#include "helpers/devices.hpp"
#include "helpers/qubo_helpers.hpp"
#include "model/qubo.hpp"
#include "model/solution.hpp"
#include "simulated_annealing/annealing.hpp"

#define SEED 1234

namespace po = boost::program_options;
namespace sycl = cl::sycl;

using queue_ptr = std::unique_ptr<sycl::queue>;

/*!
 * @brief Constructs a linear beta schedule.
 *
 * @param a reference to a vector for storing the schedule
 * @param minimum value of beta in the annealing schedule.
 * @param maximum value of beta in the annealing schedule.
 * @param number of iterations of the algorithm.
 *
 * @returns returns schedule.
 */
void construct_linear_beta_schedule(std::vector<double> &schedule,
                                    double beta_min, double beta_max,
                                    uint num_iter) {
  for (auto i = 0; i < num_iter; i++) {
    schedule[i] = beta_min + beta_max * i / static_cast<double>(num_iter - 1);
  }
}

/*!
 * @brief Constructs a geometric beta schedule.
 *
 * @param a reference to a vector for storing the schedule
 * @param minimum value of beta in the annealing schedule.
 * @param maximum value of beta in the annealing schedule.
 * @param number of iterations of the algorithm.
 *
 */
void construct_geometric_beta_schedule(std::vector<double> &schedule,
                                       double beta_min, double beta_max,
                                       uint num_iter) {
  schedule[0] = beta_min;
  auto alpha = pow(beta_max / beta_min, 1.0 / (num_iter - 1));
  for (auto i = 1; i < num_iter; i++) {
    schedule[i] = (schedule[i - 1]) * alpha;
  }
}

/*!
 * @brief Compute QUBO solution in a give solution space on a given device.
 *
 * @param Instance of the proble to solve.
 * @param Accelerator device type.
 * @param beta schedule type linear or geometric.
 * @param seed for random number generator.
 * @param minimum value of beta in the annealing schedule.
 * @param maximum value of beta in the annealing schedule.
 * @param number of iterations of the algorithm.
 * @param number of trajectories to try.
 *
 * @returns solution to the given QUBO problem and acclerator device name.
 */
std::tuple<qubo::Solution, std::string> sycl_native_anneal(qubo::QUBOModel<int, double> instance, std::string device_type, std::string schedule_type, std::uint64_t seed, double beta_min, double beta_max, uint num_iter, uint num_tries)
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
  // std::cout << "Node ID [" << machine_name << "], " << "Rank [" << rank << "], " << "Using device: "
  //           << q_ptr->get_device().get_info<sycl::info::device::name>()
  //           << std::endl;

  std::vector<double> beta_schedule(num_iter);

  if (schedule_type == "linear") {
    construct_linear_beta_schedule(beta_schedule, beta_min, beta_max,
                                   num_iter);
  } else {
    construct_geometric_beta_schedule(beta_schedule, beta_min, beta_max,
                                      num_iter);
  }

  auto solution =
      sa::anneal(instance, *q_ptr, beta_schedule, seed, num_iter, num_tries);

  return {solution, accl_device_name};
}

/*!
 * @brief Parse command line.
 * @param argument count.
 * @param arguments array.
 *
 * @returns returns input file name, output file name, device type, beta_min, beta_max, 
 *                  number of iteration and number of tries.
 */
std::tuple<std::string, std::string, std::string, std::string, double, double, uint, uint> parse_cmd_line_anneal(int argc, char *argv[])
{
  std::string input_file;
  std::string output_file;
  std::string schedule_type;
  std::string device_type;

  uint num_iter;
  uint num_tries;

  double beta_min, beta_max;

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
     "type of beta schedule to use, either linear or geometric")(
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

  if (schedule_type != "linear" && schedule_type != "geometric") {
    throw std::runtime_error("Unknown beta schedule: " + schedule_type);
  }

  if (beta_max < 0 || beta_min < 0) {
    throw std::runtime_error("Invalid schedule, both ends of beta range need to be positive.");
  }

  if (beta_min >= beta_max) {
    throw std::runtime_error("Invalid schedule, initial beta is not lesser than final beta.");
  }
  
  return {input_file, output_file, device_type, schedule_type, beta_min, beta_max, num_iter, num_tries};
}

int main(int argc, char *argv[]) {
  const clock_t begin_time = clock();

  char machine_name[MPI_MAX_PROCESSOR_NAME];
  int name_len = 0;
  int rank = 0;
  const int root_rank = 0; // Rank zero process.
  int min_energy_process_rank = 0; // Rank of process with minimum energy.
  int num_procs = 0;
  int size = 0;
  
  long long int msg_size = 0;

  std::string input_file;
  std::string output_file;
  std::string schedule_type;
  std::string device_type;
  std::vector<char> msg_buff;
  std::uint64_t *seed_buff = NULL;
  std::uint64_t seed;
  
  uint num_iter;
  uint num_tries;
  uint param_buff[2];

  double beta_min, beta_max;
  double beta_buff[2];

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
      // Parse command line to get the parameters.
      std::tie(input_file, output_file, device_type, schedule_type, beta_min, beta_max, num_iter, num_tries) = parse_cmd_line_anneal(argc, argv);

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
      instance = qubo::QUBOModel<int, double>::load(qubo_file);

      std::ostringstream oss;
      // save data to archive
      {
         boost::archive::text_oarchive ar(oss);
         // write class instance to archive
         ar << instance;
         // archive and stream closed when destructors are called.
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

    // Send schedule_type to all the processes.
    if (rank == root_rank) {
      msg_size = schedule_type.size();
      msg_buff.assign(schedule_type.c_str(),schedule_type.c_str() + schedule_type.size() +1);
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
      schedule_type = std::string(&msg_buff[0], msg_size);
#ifdef DEBUG
      std::cout << "schedule_type: " << schedule_type << std::endl;
#endif
    }

    // Send beta_min and beta_max to all the processes.
    if(rank == root_rank){
      beta_buff[0] = beta_min;
      beta_buff[1] = beta_max;
    }

    if (MPI_Bcast(&beta_buff, 2, MPI::DOUBLE, root_rank, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("MPI broadcast failed\n");
    }

    if (rank != root_rank) {
      beta_min = beta_buff[0];
      beta_max = beta_buff[1];
#ifdef DEBUG
      std::cout << "Beta range: [" << beta_min << ", " << beta_max << "]"
                << std::endl;
#endif
    }

    // Send num_iter and num_tries to all the processes.
    if (rank == root_rank) {
      param_buff[0] = num_iter;
      param_buff[1] = num_tries;
    }

    if (MPI_Bcast(&param_buff, 2, MPI_UNSIGNED, root_rank, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("MPI broadcast failed\n");
    }

    if (rank != root_rank) {
      num_iter = param_buff[0];
      num_tries = param_buff[1];
#ifdef DEBUG
      std::cout << "Number of iterations: " << num_iter << std::endl;
      std::cout << "Number of tries: " << num_tries << std::endl;
#endif
    }

    // Send the seed for random number generator to all the processes.
    if (rank == root_rank) {
      seed_buff = new std::uint64_t[num_procs];
      // Seed generation using mt19937_64 i.e 64 bit Mersenne Twister by Matsumoto, 2000.
      //std::random_device rd; // for non-deterministic generator uncomment this line and
      std::mt19937_64 gen(SEED); // replace SEED with rd() in mersenne twister.
                                 
      for (auto i=0; i < num_procs; ++i) {
        seed_buff[i] = gen();
#ifdef DEBUG
        std::cout << "seed: " << seed_buff[i] << std::endl;
#endif
      }    
    }
    
    if (MPI_Scatter(seed_buff, 1, MPI_INT64_T, &seed, 1, MPI_INT64_T, root_rank, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("MPI scatter failed\n");
    }
    
    if (rank == root_rank) {
      delete [] seed_buff;
    }

#ifdef DEBUG   
    std::cout << "seed: " << seed << std::endl;
#endif

    auto [solution, accl_device_name] = sycl_native_anneal(instance, device_type, schedule_type, seed, beta_min, beta_max, num_iter, num_tries);

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

    double *energy_buff = new double[num_procs];

    if (MPI_Gather(&solution.energy, 1, MPI_DOUBLE, energy_buff, 1, MPI_DOUBLE, root_rank, MPI_COMM_WORLD) != MPI_SUCCESS) {
      throw std::runtime_error("MPI gather failed\n");
    }

    if (rank == root_rank) {
       std::vector<double> energies(&energy_buff[0], &energy_buff[num_procs]);
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
      if (min_energy_process_rank != 0) {
        if (MPI_Recv(&buff, 32, MPI_CHAR, min_energy_process_rank, root_rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
         throw std::runtime_error("MPI receive failed\n");
        }
#ifdef DEBUG
        std::cout << "Min Energy State: "
        for (auto i = 0; i < solution.state.size(); ++i) {
        std::cout << (int)buff[i] << " ";
        }
        std::cout << std::endl;
#endif
        for (auto i=0; i < solution.state.size(); ++i ) {
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

    delete [] energy_buff;

  } catch (std::exception &e) {
    std::cerr << "error: " << e.what() << "\n";
    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  } catch (...) {
    std::cerr << "Exception of unknown type!\n";
  }

  if (rank == root_rank) {
    std::cout << "Calculation time [s]: "<< float( clock () - begin_time ) /  CLOCKS_PER_SEC << "\n";
  }

  MPI_Finalize();
  return 0;
}