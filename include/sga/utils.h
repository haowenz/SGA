#ifndef SGA_UTILS_H_
#define SGA_UTILS_H_

#include <sys/resource.h>
#include <sys/time.h>

#include <cstdint>
#include <iostream>
#include <vector>

namespace sga {
// For timing
inline static double GetRealTime() {
  struct timeval tp;
  struct timezone tzp;
  gettimeofday(&tp, &tzp);
  return tp.tv_sec + tp.tv_usec * 1e-6;
}
inline static double GetCPUTime() {
  struct rusage r;
  getrusage(RUSAGE_SELF, &r);
  return r.ru_utime.tv_sec + r.ru_stime.tv_sec +
         1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}

inline bool FileExist(const std::string& file_path) {
  FILE* file = fopen(file_path.c_str(), "r");
  return (file != NULL);
  // ifstream file_stream(file_path.c_str());
  // return file_stream.good();
  // struct stat buffer;
  // return (stat (file_path.c_str(), &buffer) == 0);
}

}  // namespace sga

#endif  // SGA_UTILS_H_
