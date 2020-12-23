/*! file device_rng.hpp
 *
 * @copyright Licensed under MIT license by hyperQ – Ewa Hendzel
 *
 */
#include "CL/sycl.hpp"
#include "oneapi/mkl/rng/device.hpp"

namespace mkl = oneapi::mkl;

namespace sa {
/*!
 * @brief Random number generator for use inside kernels.
 *
 * This is a wrapper around oneapi::mkl Engines, suitable for
 * drawing random spin indices for system of given size,
 * generating random bits and random numbers from [0, 1).
 */
template <typename Engine> class RandomGenerator {
protected:
  /*! Random engine to be used by this generator. */
  Engine engine;
  /*! System size this generato was created for. */
  int N;
  /*! Distribution of bits. */
  mkl::rng::device::uniform<int> bit_distr;
  /*! Distribution of bit indices. */
  mkl::rng::device::uniform<int> bit_index_distr;
  /*! Distribution of uniform numbers in [0, 1) interval. */
  mkl::rng::device::uniform<double> uniform_distr;

public:
  /*!
   * @brief Construct new RandomGenerator
   * @param engine random engine to be used.
   * @param N system size.
   */
  RandomGenerator(Engine engine, int N)
      : bit_distr(0, 2), bit_index_distr(0, N), uniform_distr(0.0, 1.0),
        engine(engine) {}

  /*!
   * @brief Chose value from {0, 1} at random.
   * @returns A random bit.
   */
  int bit() { return mkl::rng::device::generate(bit_distr, engine); }

  /*!
   * @brief Chose variable index uniformly at random.
   * @returns An integer from the set {0, ..., N-1},
   *          where N is the system size this generator was
   *          created for.
   */
  int bit_index() {
    return mkl::rng::device::generate(bit_index_distr, engine);
  }

  /*!
   * @brief Chose random number from [0, 1) interval.
   * @returns A random integer in [0, 1) interval.
   */
  double uniform() { return mkl::rng::device::generate(uniform_distr, engine); }
};
} // namespace sa