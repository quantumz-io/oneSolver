/*! file devices_test.cpp
 *
 * @copyright Licensed under MIT license by hyperQ – Ewa Hendzel
 *
 */
#include "CL/sycl.hpp"
#include "helpers/devices.hpp"
#include <boost/test/unit_test.hpp>

namespace sycl = cl::sycl;

BOOST_AUTO_TEST_SUITE(devices_test)

BOOST_AUTO_TEST_CASE(specifying_host_as_device_type_constructs_host_selector) {
  auto selector_ptr = devices::construct_device_selector("host");
  BOOST_TEST_REQUIRE(dynamic_cast<sycl::host_selector *>(selector_ptr.get()) !=
                     nullptr);
}

BOOST_AUTO_TEST_CASE(specifying_cpu_as_device_type_constructs_cpu_selector) {
  auto selector_ptr = devices::construct_device_selector("cpu");
  BOOST_TEST_REQUIRE(dynamic_cast<sycl::cpu_selector *>(selector_ptr.get()) !=
                     nullptr);
}

BOOST_AUTO_TEST_CASE(specifying_gpu_as_device_type_constructs_gpu_selector) {
  auto selector_ptr = devices::construct_device_selector("gpu");
  BOOST_TEST_REQUIRE(dynamic_cast<sycl::gpu_selector *>(selector_ptr.get()) !=
                     nullptr);
}

BOOST_AUTO_TEST_SUITE_END()