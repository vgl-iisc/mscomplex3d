#ifndef GRID_DATASET_CL_H
#define GRID_DATASET_CL_H

#include <grid.h>

namespace cl {class Image3D;}

namespace grid
{
namespace opencl
{

class worker
{
protected:
  std::shared_ptr<cl::Image3D> flag_img;
public:
  void assign_gradient(dataset_ptr_t ds, mscomplex_ptr_t msc=mscomplex_ptr_t());
  void owner_extrema(dataset_ptr_t ds);
  worker();
};

bool is_gpu_context();
bool is_cpu_context();

}
}


#endif
