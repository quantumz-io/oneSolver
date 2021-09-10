/*! file exhaustive.hpp
 *
 * @copyright Licensed under MIT license by hyperQ – Ewa Hendzel
 *
 */
#ifndef EXHAUSTIVE_CALCULATION_HPP__
#define EXHAUSTIVE_CALCULATION_HPP__

#include "CL/sycl.hpp"
#include "helpers/ulong_to_vec.hpp"
#include "model/qubo.hpp"
#include "model/solution.hpp"



namespace exhaustive {
namespace sycl = cl::sycl;

/**
 * @brief Solve, by finding the minimum energy state, a QUBO model using
 * exhuastive search.
 *
 * @tparam NodeType Integer type of the variables addresses.
 * @tparam CoefType Real type of the variables values.
 * @param q Queue on which the algorithm will be run.
 * @param qubos QUBO model to solve.
 * @param start state in solution search space.
 * @param end state in solution search space.
 * @return qubo::Solution Solution, contaning minimal energy and a
 * minimal energy state.
 */

template <class NodeType, class CoefType>
qubo::Solution solve(sycl::queue &q,
                     qubo::QUBOModel<NodeType, CoefType> &qubos, ulong start_state, ulong end_state) {

  auto wg_size =
      q.get_device().get_info<sycl::info::device::max_work_group_size>();
  auto max_compute_units =
      q.get_device().get_info<sycl::info::device::max_compute_units>();
  auto n_bits = qubos.get_nodes();

  std::vector<double> energies;
  std::vector<ulong> states;

  {

    sycl::buffer<double, 2> qubo_buf(sycl::range<2>(n_bits, n_bits));
    for (auto i = 0; i < n_bits; ++i) {
      for (auto j = i; j < n_bits; ++j) {
        if (i == j) {
#ifdef DEBUG
          std::cout << i << " " << j << " " << qubos.get_variable(i) << "\n";
#endif
          qubo_buf.get_host_access()[i][i] = qubos.get_variable(i);
        } else {
#ifdef DEBUG
          std::cout << i << " " << j << " "
                    << qubos.get_connection(std::make_pair(i, j)) << "\n";
#endif
          qubo_buf.get_host_access()[i][j] =
              qubos.get_connection(std::make_pair(i, j));
        }
      }
    }

    auto num_threads = max_compute_units;
    sycl::buffer<ulong, 2> ranges_buf(sycl::range<2>(num_threads, 2));

    // note that all the variables related to states have to be unsigned 64 bit integers
    //ulong n_states = 1 << n_bits; 
    ulong n_states = end_state - start_state;
    ulong states_per_thread = n_states / num_threads;
    ulong states_per_thread_reminder = n_states % num_threads;

    ulong reminder_count = 0;
    //ulong states_count = 0;
    ulong states_count = start_state;

    for (auto i = 0; i < num_threads; ++i) {
      ulong state_start = states_count;
      states_count += states_per_thread;
      if (reminder_count < states_per_thread_reminder) {
        states_count += 1;
        reminder_count++;
      }
      auto state_end = states_count;

      ranges_buf.get_host_access()[i][0] = state_start;
      ranges_buf.get_host_access()[i][1] = state_end;
#ifdef DEBUG
      std::cout << ranges_buf.get_host_access()[i][0];
      std::cout << ", ";
      std::cout << ranges_buf.get_host_access()[i][1] << std::endl;
#endif
    }
    sycl::buffer<CoefType, 1> energy_buf(num_threads);
    sycl::buffer<ulong, 1> state_buf(num_threads);

    q.submit([&](sycl::handler &h) {
       auto qubo_acc =
           qubo_buf.template get_access<sycl::access::mode::read>(h);
       auto ranges_acc =
           ranges_buf.template get_access<sycl::access::mode::read>(h);
       auto energy_acc =
           energy_buf.template get_access<sycl::access::mode::write>(h);
       auto state_acc =
           state_buf.template get_access<sycl::access::mode::write>(h);

       h.parallel_for<class calc_energy>(
           sycl::range<1>(num_threads), [=](sycl::id<1> item) {
             CoefType e_best = std::numeric_limits<CoefType>::max();
             ulong state_best = 0;
             auto local_start = ranges_acc[item][0];
             auto local_end = ranges_acc[item][1];
             for (ulong state = local_start; state < local_end; ++state) {
               // calculate energy for a given state
               CoefType e = 0.0;

               for (size_t i = 0; i < n_bits; ++i) {
                 for (size_t j = i; j < n_bits; ++j) {
                   if (i == j) {
                     auto i_bit = ((state >> i) & 1);
                     if (i_bit) {
                       e += qubo_acc[i][i];
                     }
                   } else {
                     auto i_bit = ((state >> i) & 1);
                     auto j_bit = ((state >> j) & 1);
                     if (i_bit and j_bit) {
                       e += qubo_acc[i][j];
                     }
                   }
                 }
               }
               if (e < e_best) {
                 e_best = e;
                 state_best = state;
               }
             }
             state_acc[item] = state_best;
             energy_acc[item] = e_best;
           });  // end parallel for
     }).wait(); // end submit

    for (auto i = 0; i < energy_buf.get_count(); ++i) {
#ifdef DEBUG
      std::cout << i << ", " << energy_buf.get_host_access()[i] << std::endl;
#endif
      energies.push_back(energy_buf.get_host_access()[i]);
    }
    for (auto i = 0; i < state_buf.get_count(); ++i) {
#ifdef DEBUG
      std::cout << i << ", ";
      for (auto &bit : helpers::ulong_to_vec(state_buf.get_host_access()[i], n_bits)) {
        std::cout << (int)bit << " ";
      }
      std::cout << std::endl;
#endif
      states.push_back(state_buf.get_host_access()[i]);
    }
  }

  auto min_energy = std::min_element(energies.begin(), energies.end());
  auto min_idx_energy = std::distance(energies.begin(), min_energy);
#ifdef DEBUG
  std::cout << "min element at: " << min_idx_energy << std::endl;
#endif

  auto result = helpers::ulong_to_vec(states[min_idx_energy], n_bits);

  return qubo::Solution(result.begin(), result.end(), energies[min_idx_energy]);
}

}; // namespace exhaustive

#endif
