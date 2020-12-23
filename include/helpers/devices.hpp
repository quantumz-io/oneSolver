/*! file devices.hpp
 *
 * @copyright Licensed under MIT license by hyperQ – Ewa Hendzel
 *
 */
#include "CL/sycl.hpp"

namespace sycl = cl::sycl;

/*!
 * @brief Unique pointer to a device selector.
 */
using device_selector_ptr = std::unique_ptr<sycl::device_selector>;

/*!
 * @brief Helper functions for dealing with devices and selectors.
 */
namespace devices {

/*!
 * @brief Construct device selector based on specified device type.
 * @param device_type string identifying device type. Must be one of
 *                    "cpu", "gpu" or "host".
 * @throw std::invalid_argument if \p device_type is incorrect.
 * @returns Pointer to a selector corresponding to \p device_type.
 */
device_selector_ptr construct_device_selector(std::string device_type) {

  device_selector_ptr result;
  if (device_type == "cpu") {
    result = device_selector_ptr(new sycl::cpu_selector());
  } else if (device_type == "gpu") {
    result = device_selector_ptr(new sycl::gpu_selector());
  } else if (device_type == "host") {
    result = device_selector_ptr(new sycl::host_selector());
  } else {
    std::ostringstream error_stream;
    error_stream << "Unknown device type: " << device_type;
    throw std::invalid_argument(error_stream.str());
  }
  return result;
}
} // namespace devices