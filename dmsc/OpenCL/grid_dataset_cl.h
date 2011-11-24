#ifndef GRID_DATASET_CL_H
#define GRID_DATASET_CL_H

#include <grid.h>

#define __CL_ENABLE_EXCEPTIONS
#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.hpp>
#else
#include <cl.hpp>
#endif

namespace grid
{
  namespace opencl
  {
    void init();

    class worker
    {
    protected:
      cl::Image3D flag_img;
    public:
      void assign_gradient(dataset_ptr_t ds, mscomplex_ptr_t msc);
      void owner_extrema(dataset_ptr_t ds);
    };

    void assign_gradient_and_owner_extrema(dataset_ptr_t ds);
    void update_to_surv_extrema(dataset_ptr_t ds,mscomplex_ptr_t msc);

//    void check_assign_gradient(dataset_ptr_t ds, int dim,rect_t check_rect,cell_flag_t mask);
  }
}


#endif
