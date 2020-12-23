/*! file annealing.hpp
 *
 * @copyright Licensed under MIT license by hyperQ – Ewa Hendzel
 *
 */
#include "CL/sycl.hpp"
#include "helpers/qubo_helpers.hpp"
#include "model/solution.hpp"
#include "oneapi/mkl/rng/device.hpp"
#include "simulated_annealing/device_rng.hpp"

namespace sycl = cl::sycl;
namespace mkl = oneapi::mkl;

using mode = sycl::access::mode;

/*!
 * @brief Implementation of simulated annealing.
 */
namespace sa {

/*!
 * @brief Compute QUBO energy for given system configuration.
 * @param flat_qubo array with qubo coefficients. Only the upper
 *                  diagonal is used.
 * @param state system configuration.
 * @param N number of variables.
 *
 * @returns energy of the system.
 */
template <typename QuboArray, typename StateArray>
double energy(const QuboArray flat_qubo, StateArray state, int N) {
  double result = 0.0;
  for (auto i = 0; i < N; i++) {
    for (auto j = i; j < N; j++) {
      result += flat_qubo[i * N + j] * state[i] * state[j];
    }
  }
  return result;
}

/*!
 * @brief Solve given QUBO using simulated annealing.
 * @param instance the problem to solve.
 * @param q queue to use for computations.
 * @param h_beta_schedule array with values of beta. Note that
 *                        those values should follow decreasing order.
 * @param num_iter number of iterations of every pass of the algorithm.
 *                 It should conincide with length of \p h_beta_schedule
 * @param num_tries number of independent passes of the algorithm.
 * @param sweeps_per_beta number of spin flips for each
 *                        annealing step.
 * @return An object containing the lowest energy solution.
 */
template <typename T>
qubo::Solution anneal(qubo::QUBOModel<int, T> instance, sycl::queue q,
                      std::vector<double> &h_beta_schedule, int num_iter,
                      unsigned int num_tries, int sweeps_per_beta = 1) {
  auto N = instance.get_nodes();
  auto qubo_vector = helpers::flatten_qubo(instance);
  std::vector<double> best_energies(num_tries);
  std::vector<char> best_states(num_tries * N);

  try {
    sycl::buffer<double, 1> best_energies_buf(best_energies.data(),
                                              best_energies.size());
    sycl::buffer<double, 1> flat_qubo_buf(qubo_vector.data(),
                                          qubo_vector.size());
    sycl::buffer<double, 1> beta_schedule_buf(h_beta_schedule);
    sycl::buffer<char, 2> best_states_buf(best_states.data(),
                                          sycl::range<2>(num_tries, N));
    sycl::buffer<char, 2> current_states_buf(sycl::range<2>(num_tries, N));

    q.submit([&](sycl::handler &h) {
       auto best_energies_acc =
           best_energies_buf.template get_access<mode::write>(h);
       auto flat_qubo = flat_qubo_buf.template get_access<mode::read>(h);
       auto beta_schedule =
           beta_schedule_buf.template get_access<mode::read>(h);
       auto best_states =
           best_states_buf.template get_access<mode::read_write>(h);
       auto current_states =
           current_states_buf.template get_access<mode::read_write>(h);

       h.parallel_for<class annealing>(
           sycl::range<1>{num_tries}, [=](sycl::id<1> i) {
             mkl::rng::device::philox4x32x10 engine(1234, i);
             sa::RandomGenerator random(engine, N);

             for (auto j = 0; j < N; j++) {
               best_states[i][j] = current_states[i][j] = random.bit();
             }

             double best_energy = energy(flat_qubo, current_states[i], N);
             double current_energy = best_energy;

             for (auto iter = 0; iter < num_iter; iter++) {
               double beta = beta_schedule[iter];

               for (auto sweep = 0; sweep < sweeps_per_beta; sweep++) {
                 auto spin_to_flip = random.bit_index();
                 current_states[i][spin_to_flip] =
                     1 - current_states[i][spin_to_flip];
                 auto new_energy = energy(flat_qubo, current_states[i], N);

                 if ((new_energy < current_energy) ||
                     (sycl::exp((current_energy - new_energy) / beta) >
                      random.uniform())) {
                   current_energy = new_energy;
                 } else {
                   current_states[i][spin_to_flip] =
                       1 - current_states[i][spin_to_flip];
                 }

                 if (current_energy < best_energy) {
                   best_energy = current_energy;

                   for (auto j = 0; j < N; j++) {
                     best_states[i][j] = current_states[i][j];
                   }
                 }
               }
             }

             best_energies_acc[i] = best_energy;
           });
     }).wait();
  } catch (sycl::exception &e) {
    std::cerr
        << "Execution of kernel failed. See below exception info for details";
    throw e;
  }

  auto best_idx = std::min_element(best_energies.begin(), best_energies.end()) -
                  best_energies.begin();

  return qubo::Solution(best_states.begin() + N * best_idx,
                        best_states.begin() + N * (best_idx + 1),
                        best_energies[best_idx]);
}
} // namespace sa