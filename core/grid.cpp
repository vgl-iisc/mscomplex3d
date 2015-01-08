#include <grid.h>

#include <omp.h>

using namespace std;

namespace grid
{
utl::timer g_timer;

namespace openmp{

std::string get_info()
{
  std::stringstream ss;

  ss << "====================================" << endl
     << "    OpenMP Info                     " << endl
     << "------------------------------------" << endl
     << "MaxThreads     = " << omp_get_max_threads() << endl
     << "NumProcs       = " << omp_get_num_procs() << endl
     << "====================================" << endl;

  return ss.str();
}




}
}
