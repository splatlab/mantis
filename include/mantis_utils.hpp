#ifndef __MANTIS_UTILS_HPP__
#define __MANTIS_UTILS_HPP__

#include <chrono>

namespace mantis{
  std::string get_current_time_as_string() {
    // Get the time at the start of the run
    std::time_t result = std::time(nullptr);
    auto time = std::string(std::asctime(std::localtime(&result)));
    time.pop_back(); // remove the newline
    return time;
  }
}

#endif // __MANTIS_UTILS_HPP__
