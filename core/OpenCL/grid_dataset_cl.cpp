#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <stdexcept>

#define __CL_ENABLE_EXCEPTIONS

#ifndef __has_cpp_attribute
#define __has_cpp_attribute(x) 0
#endif

#ifndef __has_attribute
#define __has_attribute(x) 0
#endif

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.hpp>
#else
#if defined(__linux__)
#define CL_HPP_TARGET_OPENCL_VERSION 220
#endif
#include<opencl.hpp>
#endif

#include <utl.h>

#include<grid_dataset_cl.h>
#include<grid_dataset.h>
#include<grid_mscomplex.h>

const int WI_SIZE = OPENCL_NUM_WORK_ITEMS_PER_GROUP;
const int WG_NUM  = OPENCL_NUM_WORK_GROUPS;
const int WG_SIZE = WG_NUM*WI_SIZE;

cl::Platform      s_platform;
cl::Device        s_device;
cl::Context       s_context;
cl::CommandQueue  s_queue;

/*
cl::KernelFunctor s_assign_max_facet_edge;
cl::KernelFunctor s_assign_max_facet_face;
cl::KernelFunctor s_assign_max_facet_cube;
cl::KernelFunctor s_assign_pairs;
cl::KernelFunctor s_assign_pairs2;
cl::KernelFunctor s_assign_pairs3;

cl::KernelFunctor s_mark_cps;
cl::KernelFunctor s_mark_boundry_cps;
cl::KernelFunctor s_count_cps;
cl::KernelFunctor s_count_boundry_cps;
cl::KernelFunctor s_save_boundry_cps;
cl::KernelFunctor s_save_boundry_cps;


cl::KernelFunctor s_scan_local_sums;
cl::KernelFunctor s_scan_group_sums;
cl::KernelFunctor s_scan_update_sums;

cl::KernelFunctor s_init_propagate;
cl::KernelFunctor s_propagate;
cl::KernelFunctor s_update_cp_cell_to_cp_no;
cl::KernelFunctor s_init_update_to_surv_cp_no;
cl::KernelFunctor s_update_to_surv_cp_no;
*/

// Declare kernels outside the scope for reuse
cl::Kernel s_assign_max_facet_edge;
cl::Kernel s_assign_max_facet_face;
cl::Kernel s_assign_max_facet_cube;
cl::Kernel s_assign_pairs;
cl::Kernel s_assign_pairs2;
cl::Kernel s_assign_pairs3;

cl::Kernel s_mark_cps;
cl::Kernel s_mark_boundry_cps;
cl::Kernel s_count_cps;
cl::Kernel s_count_boundry_cps;
cl::Kernel s_save_boundry_cps;
cl::Kernel s_save_cps;


cl::Kernel s_scan_local_sums;
cl::Kernel s_scan_group_sums;
cl::Kernel s_scan_update_sums;

cl::Kernel s_init_propagate;
cl::Kernel s_propagate;
cl::Kernel s_update_cp_cell_to_cp_no;
cl::Kernel s_init_update_to_surv_cp_no;
cl::Kernel s_update_to_surv_cp_no;


extern const char *GRID_DATASET_CLH;
extern const char *GRID_DATASET_ASSIGNGRADIENT_CL;
extern const char *GRID_DATASET_MARKANDCOLLECT_CL;
extern const char *GRID_DATASET_OWNEREXTREMA_CL;

using namespace std;

namespace utl
{
  template <> inline std::string to_string (const unsigned char & t)
  {
    char t_str[9];

    t_str[0] = '0'+char((t>>7)&1);
    t_str[1] = '0'+char((t>>6)&1);
    t_str[2] = '0'+char((t>>5)&1);
    t_str[3] = '0'+char((t>>4)&1);
    t_str[4] = '0'+char((t>>3)&1);
    t_str[5] = '0'+char((t>>2)&1);
    t_str[6] = '0'+char((t>>1)&1);
    t_str[7] = '0'+char((t>>0)&1);
    t_str[8] =  0;
    return t_str;
  }
}

namespace grid
{
  namespace opencl
  {
    typedef cl_short4 cell_t;
    typedef cl_short8   cell_pair_t;

    inline cell_t to_cell(const cellid_t & b)
    {
      cell_t a;

      a.x = b[0];
      a.y = b[1];
      a.z = b[2];

      return a;
    }

    inline cellid_t from_cell(const cell_t & b)
    {
      cellid_t a;

      a[0] = b.x;
      a[1] = b.y;
      a[2] = b.z;

      return a;
    }


    inline cell_pair_t to_cell_pair(const rect_t & b)
    {
      cell_pair_t a;

      a.lo = to_cell(b.lc());
      a.hi = to_cell(b.uc());

      return a;
    }

    inline rect_t from_cell_pair(const cell_pair_t & r)
    {
      return rect_t(from_cell(r.lo),from_cell(r.hi));
    }

    inline cell_t to_cell(int x , int y , int z)
    {
      cell_t a;

      a.x = x;
      a.y = y;
      a.z = z;

      return a;
    }

    //inline cl::size_t<3> to_size(const cellid_t & b)
    inline std::array<size_t, 3> to_size(const cellid_t & b)
    {
    	std::array<size_t, 3> a;

      //a.assign(b.begin(),b.end());
        for (size_t i = 0; i < 3; ++i) {
            a[i] = static_cast<size_t>(b[i]);  // Convert from short to size_t
        }


      return a;
    }

    //inline cl::size_t<3> to_size(int x,int y,int z)
    inline std::array<size_t, 3> to_size(int x,int y,int z)
    {
      //cl::size_t<3> a;

        std::array<size_t, 3> a;

      a[0] = x;
      a[1] = y;
      a[2] = z;

      return a;
    }

    //inline cl::size_t<3> get_size(const cell_pair_t &p)
    inline std::array<size_t, 3> get_size(const cell_pair_t &p)
    {
      //cl::size_t<3> s;
      std::array<size_t, 3> s;

      s[0] = p.hi.x -p.lo.x +1;
      s[1] = p.hi.y -p.lo.y +1;
      s[2] = p.hi.z -p.lo.z +1;

      return s;
    }

    inline int num_cells(const cell_pair_t &p)
    {
      return grid::num_cells(from_cell_pair(p));
    }

    template<typename T>
    void log_buffer(T*  buf_cpu,int n,int xrepeat = -1, int yrepeat = -1,std::ostream &os=cout)
    {
      for( int i = 0 ; i < n; ++i)
      {
        if( (xrepeat > 0)  && (i%xrepeat == 0))
        {
          os<<endl;

          if((yrepeat > 0 ) && ((i/xrepeat)%yrepeat == 0))
            os<<"sheet no" << i/(xrepeat*yrepeat)<<endl;
        }

        os << buf_cpu[i]<< " ";
      }

      os<<endl;
    }

    template<typename T>
    void log_buffer(cl::Buffer buf,int n,int xrepeat = -1, int yrepeat = -1,std::ostream &os=cout)
    {
      std::vector<T> buf_cpu(n);

      s_queue.finish();

      s_queue.enqueueReadBuffer(buf,true,0,sizeof(T)*n,buf_cpu.data());

      log_buffer<T>(buf_cpu.data(),n,xrepeat,yrepeat,os);
    }

    std::string get_info()
    {
      std::stringstream ss;

      ss << "====================================" << endl
         << "    OpenCL platform/device Info     " << endl
         << "------------------------------------" << endl
         << "PlatformName   = " << s_platform.getInfo<CL_PLATFORM_NAME>() << endl
         << "PlatformVendor = " << s_platform.getInfo<CL_PLATFORM_VENDOR>() << endl
         << "DeviceName     = " << s_device.getInfo<CL_DEVICE_NAME>() << endl
         << "DeviceVendor   = " << s_device.getInfo<CL_DEVICE_VENDOR>() << endl
         << "DeviceMemory(b)= " << s_device.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>() << "(b)" << endl
         << "MaxMemAlloc(b) = " << s_device.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>() <<"(b)" << endl
         << "====================================" << endl;


      return ss.str();
    }


    void init(int device)
    {
      std::vector<cl::Device> devices;
      cl::Program             program1;
      cl::Program             program2;
      cl::Program             program3;


      string headerCode (GRID_DATASET_CLH);
      string sourceCode1(GRID_DATASET_ASSIGNGRADIENT_CL);
      string sourceCode2(GRID_DATASET_MARKANDCOLLECT_CL);
      string sourceCode3(GRID_DATASET_OWNEREXTREMA_CL);

      const string s1("@OPENCL_NUM_WORK_ITEMS_PER_GROUP@");
      const string s2("@OPENCL_NUM_WORK_GROUPS@");

      headerCode.replace(headerCode.find(s1),s1.size(),
                         utl::to_string(OPENCL_NUM_WORK_ITEMS_PER_GROUP));

      headerCode.replace(headerCode.find(s2),s2.size(),
                         utl::to_string(OPENCL_NUM_WORK_GROUPS));
      try
      {
          std::vector<cl::Platform> platforms;
          cl::Platform::get(&platforms);

          if (platforms.empty())
              throw std::runtime_error("No OpenCL platforms found!");

          std::vector<std::pair<int, int>> availablePlatformDeviceIDs;
          int selectedPlatformDeviceID_GPU = -1;
          int selectedPlatformDeviceID_CPU = -1;

          for (int i = 0; i < platforms.size(); ++i)
          {
              std::vector<cl::Device> devices;
              platforms[i].getDevices(CL_DEVICE_TYPE_ALL, &devices);

              for (int j = 0; j < devices.size(); ++j)
              {
                  availablePlatformDeviceIDs.push_back(std::make_pair(i, j));

                  if (selectedPlatformDeviceID_GPU == -1 && devices[j].getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_GPU)
                      selectedPlatformDeviceID_GPU = availablePlatformDeviceIDs.size() - 1;

                  if (devices[j].getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU)
                      selectedPlatformDeviceID_CPU = availablePlatformDeviceIDs.size() - 1;
              }
          }

          std::pair<int, int> selectedPlatformDevicePair;
          if (device == 0 && selectedPlatformDeviceID_GPU != -1)
          {
              selectedPlatformDevicePair = availablePlatformDeviceIDs[selectedPlatformDeviceID_GPU];
              std::cout << "\nGPU SELECTED\n";
          }
          else if (device == 1 && selectedPlatformDeviceID_CPU != -1)
          {
              selectedPlatformDevicePair = availablePlatformDeviceIDs[selectedPlatformDeviceID_CPU];
              std::cout << "\nCPU SELECTED\n";
          }
          else
          {
              selectedPlatformDevicePair = availablePlatformDeviceIDs[0]; // Default to first available
              std::cout << "\nDEFAULT DEVICE SELECTED\n";
          }

          int platformIndex = selectedPlatformDevicePair.first;
          int deviceIndex = selectedPlatformDevicePair.second;

          if (platformIndex >= platforms.size())
              throw std::runtime_error("Error: Selected platform index is out of bounds.");

          s_platform = platforms[platformIndex];

          std::vector<cl::Device> platformDevices;
          s_platform.getDevices(CL_DEVICE_TYPE_ALL, &platformDevices);

          if (deviceIndex >= platformDevices.size())
              throw std::runtime_error("Error: Selected device index is out of bounds.");

          s_device = platformDevices[deviceIndex];

          cl_context_properties properties[] = { CL_CONTEXT_PLATFORM, (cl_context_properties)s_platform(), 0 };
          s_context = cl::Context({ s_device }, properties);

          s_queue = cl::CommandQueue(s_context, s_device);

          std::cout << "Using OpenCL device: " << s_device.getInfo<CL_DEVICE_NAME>() << std::endl;
          std::cout << "Using platform: " << s_platform.getInfo<CL_PLATFORM_NAME>() << std::endl;
      }
      catch (cl::Error& err)
      {
          std::cerr << "SETUP QUEUE ERROR: " << err.what() << " (" << err.err() << ")" << std::endl;
          throw;
      }
        
      try {
          // Create program sources
          cl::Program::Sources sources = {
              { headerCode.c_str(), headerCode.size() },
              { sourceCode1.c_str(), sourceCode1.size() }
          };

          // Build the program
          cl::Program program1(s_context, sources);
          program1.build(devices);

          // Define NDRange sizes
          cl::NDRange global(1024);  // 1024 work-items in total
          //cl::NDRange local(256);    // 256 work-items per work group
          cl::NDRange local(64);    // 256 work-items per work group

          // Define and initialize kernels
          s_assign_max_facet_edge=cl::Kernel (program1, "assign_max_facet_edge");
          s_assign_max_facet_face=cl::Kernel (program1, "assign_max_facet_face");
          s_assign_max_facet_cube=cl::Kernel (program1, "assign_max_facet_cube");
          s_assign_pairs=cl::Kernel (program1, "assign_pairs");
          s_assign_pairs2=cl::Kernel (program1, "assign_pairs2");
          s_assign_pairs3=cl::Kernel (program1, "assign_pairs3");

          // Enqueue kernels with specific NDRange settings
          //s_queue.enqueueNDRangeKernel(s_assign_max_facet_edge, cl::NullRange, global, local);
          //s_queue.enqueueNDRangeKernel(s_assign_max_facet_face, cl::NullRange, global, local);
          //s_queue.enqueueNDRangeKernel(s_assign_max_facet_cube, cl::NullRange, global, local);
          //s_queue.enqueueNDRangeKernel(s_assign_pairs, cl::NullRange, global, cl::NDRange(local[0] / 2));
          //s_queue.enqueueNDRangeKernel(s_assign_pairs2, cl::NullRange, global, cl::NDRange(local[0] / 2));
          //s_queue.enqueueNDRangeKernel(s_assign_pairs3, cl::NullRange, global, cl::NDRange(local[0] / 2));

          // Wait for queue completion
          //s_queue.finish();

      }
      catch (const cl::Error& e) {
          std::cerr << "PROGRAM1 ERROR: " << e.what() << " (" << e.err() << ")\n";
          std::cerr << program1.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]) << "\n";

          throw;
      }
      
      
      try {
          // Create program sources
          cl::Program::Sources sources = {
              { headerCode.c_str(), headerCode.size() },
              { sourceCode2.c_str(), sourceCode2.size() }
          };

          // Build the program
          cl::Program program2(s_context, sources);
          program2.build(devices);

          // Initialize kernels
          s_count_cps = cl::Kernel(program2, "count_cps");
          s_count_boundry_cps = cl::Kernel(program2, "count_boundry_cps");
          s_scan_local_sums = cl::Kernel(program2, "scan_local_sums");
          s_scan_group_sums = cl::Kernel(program2, "scan_group_sums");
          s_scan_update_sums = cl::Kernel(program2, "scan_update_sums");
          s_mark_cps = cl::Kernel(program2, "mark_cps");
          s_save_cps=cl::Kernel(program2, "save_cps");
          s_mark_boundry_cps= cl::Kernel(program2, "mark_boundry_cps");
          s_save_boundry_cps = cl::Kernel(program2, "save_boundry_cps");



          // Define NDRange sizes
          //cl::NDRange globalWorkSize(WG_SIZE);
          //cl::NDRange localWorkSize(WI_SIZE);
          //cl::NDRange groupWorkSize(WG_NUM);

          // Enqueue kernels with specific NDRange settings
          //s_queue.enqueueNDRangeKernel(s_mark_cps, cl::NullRange, globalWorkSize, localWorkSize);
          //s_queue.enqueueNDRangeKernel(s_count_cps, cl::NullRange, globalWorkSize, localWorkSize);
          //s_queue.enqueueNDRangeKernel(s_save_cps, cl::NullRange, globalWorkSize, localWorkSize);
          //s_queue.enqueueNDRangeKernel(s_mark_boundry_cps, cl::NullRange, globalWorkSize, localWorkSize);
          //s_queue.enqueueNDRangeKernel(s_count_boundry_cps, cl::NullRange, globalWorkSize, localWorkSize);
          //s_queue.enqueueNDRangeKernel(s_save_boundry_cps, cl::NullRange, globalWorkSize, localWorkSize);
          //s_queue.enqueueNDRangeKernel(s_scan_local_sums, cl::NullRange, globalWorkSize, localWorkSize);
          //s_queue.enqueueNDRangeKernel(s_scan_group_sums, cl::NullRange, groupWorkSize, groupWorkSize);
          //s_queue.enqueueNDRangeKernel(s_scan_update_sums, cl::NullRange, globalWorkSize, localWorkSize);

          // Wait for queue completion
          //s_queue.finish();

      }
      catch (const cl::Error& e) {
          std::cerr << "OpenCL error: " << e.what() << " (" << e.err() << ")\n";
          throw;
      }
        /*
      try
      {
        cl::Program::Sources sources;
        //sources.push_back(make_pair(headerCode.c_str(),headerCode.size()));
        sources.push_back({ headerCode.c_str() ,(headerCode.size()) });

        //sources.push_back(make_pair(sourceCode2.c_str(), sourceCode2.size()));
        sources.push_back({ sourceCode2.c_str() ,(sourceCode2.size()) });

        program2 = cl::Program(s_context, sources);
        program2.build(devices);

        s_mark_cps = cl::Kernel(program2, "mark_cps").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_count_cps = cl::Kernel(program2, "count_cps").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_save_cps = cl::Kernel(program2, "save_cps").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_mark_boundry_cps = cl::Kernel(program2, "mark_boundry_cps").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_count_boundry_cps = cl::Kernel(program2, "count_boundry_cps").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_save_boundry_cps = cl::Kernel(program2, "save_boundry_cps").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_scan_local_sums = cl::Kernel(program2, "scan_local_sums").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_scan_group_sums= cl::Kernel(program2, "scan_group_sums").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_NUM),cl::NDRange(WG_NUM));

        s_scan_update_sums = cl::Kernel(program2, "scan_update_sums").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));





      }
      catch (cl::Error err)
      {
       cerr<< "PROGRAM2 ERROR: "<< err.what()<< "("<< err.err()<< ")"<< endl;
       cerr<<program2.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0])<<endl;

       throw;
      }
*/

      try {
          // Create program sources
          cl::Program::Sources sources = {
              { headerCode.c_str(), headerCode.size() },
              { sourceCode3.c_str(), sourceCode3.size() }
          };

          // Build the program
          cl::Program program3(s_context, sources);
          program3.build(devices);

          // Define and initialize kernels
          s_init_propagate=cl::Kernel (program3, "init_propagate");
          s_propagate=cl::Kernel (program3, "propagate");
          s_update_cp_cell_to_cp_no=cl::Kernel (program3, "update_cp_cell_to_cp_no");
          s_init_update_to_surv_cp_no=cl::Kernel (program3, "init_update_to_surv_cp_no");
          s_update_to_surv_cp_no=cl::Kernel (program3, "update_to_surv_cp_no");

          // Define NDRange sizes
          //cl::NDRange globalWorkSize(WG_SIZE);
          //cl::NDRange localWorkSize(WI_SIZE);

          // Enqueue kernels with specific NDRange settings
          //s_queue.enqueueNDRangeKernel(s_init_propagate, cl::NullRange, globalWorkSize, localWorkSize);
          //s_queue.enqueueNDRangeKernel(s_propagate, cl::NullRange, globalWorkSize, localWorkSize);
          //s_queue.enqueueNDRangeKernel(s_update_cp_cell_to_cp_no, cl::NullRange, globalWorkSize, localWorkSize);
          //s_queue.enqueueNDRangeKernel(s_init_update_to_surv_cp_no, cl::NullRange, globalWorkSize, localWorkSize);
          //s_queue.enqueueNDRangeKernel(s_update_to_surv_cp_no, cl::NullRange, globalWorkSize, localWorkSize);

          // Wait for queue completion
          //s_queue.finish();

      }
      catch (const cl::Error& e) {
          std::cerr << "PROGRAM3 ERROR: " << e.what() << " (" << e.err() << ")\n";
          std::cerr << program3.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]) << "\n";
          throw;
      }
      

        /*
      try
      {
        cl::Program::Sources sources;
        sources.push_back(make_pair(headerCode.c_str(),headerCode.size()));
        sources.push_back(make_pair(sourceCode3.c_str(),sourceCode3.size()));

        program3 = cl::Program(s_context, sources);
        program3.build(devices);

        s_init_propagate =  cl::Kernel(program3, "init_propagate").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_propagate =  cl::Kernel(program3, "propagate").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_update_cp_cell_to_cp_no =  cl::Kernel(program3, "update_cp_cell_to_cp_no").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_init_update_to_surv_cp_no = cl::Kernel(program3, "init_update_to_surv_cp_no").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

        s_update_to_surv_cp_no =  cl::Kernel(program3, "update_to_surv_cp_no").
            bind(s_queue,cl::NullRange,cl::NDRange(WG_SIZE),cl::NDRange(WI_SIZE));

      }
      catch (cl::Error err)
      {
       cerr<< "PROGRAM3 ERROR: "<< err.what()<< "("<< err.err()<< ")"<< endl;
       cerr<<program3.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0])<<endl;

       throw;
      }
      */


    }

    bool is_gpu_context()
    {
      return (s_device.getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_GPU);
    }

    bool is_cpu_context()
    {
        return (s_device.getInfo<CL_DEVICE_TYPE>() == CL_DEVICE_TYPE_CPU);
    }

    void __assign_gradient
      ( cell_pair_t rct,
        cell_pair_t ext,
        cell_pair_t dom,
        cl::Image3D &func_img,
        cl::Image3D &flag_img,
        cell_fn_t   *h_func,
        cell_flag_t *h_flag)
    {
      //cl::size_t<3> func_size = to_size(from_cell_pair(ext).span()/2+ 1);
      std::array<size_t,3> func_size = to_size(from_cell_pair(ext).span()/2+ 1);
      //cl::size_t<3> flag_size = get_size(ext);
      std::array<size_t, 3> flag_size = get_size(ext);

      int cell_ct = num_cells(ext);

      cl::NDRange globalSize(1024); // Total number of work items
      cl::NDRange localSize(256);          // Work items per work group

      
      try
      {
        func_img = cl::Image3D(s_context,CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,
                               cl::ImageFormat(CL_R,CL_FLOAT),
                               func_size[0],func_size[1],func_size[2],0,0,h_func);

        // Check device limits
        cl::Device device = s_context.getInfo<CL_CONTEXT_DEVICES>()[0];
        size_t max_width = device.getInfo<CL_DEVICE_IMAGE3D_MAX_WIDTH>();
        size_t max_height = device.getInfo<CL_DEVICE_IMAGE3D_MAX_HEIGHT>();
        size_t max_depth = device.getInfo<CL_DEVICE_IMAGE3D_MAX_DEPTH>();
        
        

        // Check flag size validity
        if (flag_size[0] == 0 || flag_size[1] == 0 || flag_size[2] == 0 ||
            flag_size[0] > max_width || flag_size[1] > max_height || flag_size[2] > max_depth) {
            std::cout << "\nFlag Size: " << flag_size[0] << " x " << flag_size[1] << " x " << flag_size[2];
            std::cout << "\nDevice max 3D image size: "
                << max_width << " x " << max_height << " x " << max_depth << std::endl;
            throw std::runtime_error("Invalid image dimensions for flag_img!");
        }
        
        flag_img = cl::Image3D(s_context, CL_MEM_READ_ONLY,
                               cl::ImageFormat(CL_R,CL_UNSIGNED_INT8),
                               flag_size[0],flag_size[1],flag_size[2],0,0);
                               
        cl::Buffer flag_buf(s_context,CL_MEM_READ_WRITE,cell_ct*sizeof(cell_flag_t));

        s_assign_max_facet_edge.setArg(0, func_img);
        s_assign_max_facet_edge.setArg(1, flag_img);
        s_assign_max_facet_edge.setArg(2, rct.lo);
        s_assign_max_facet_edge.setArg(3, rct.hi);
        s_assign_max_facet_edge.setArg(4, ext.lo);
        s_assign_max_facet_edge.setArg(5, ext.hi);
        s_assign_max_facet_edge.setArg(6, dom.lo);
        s_assign_max_facet_edge.setArg(7, dom.hi);
        s_assign_max_facet_edge.setArg(8, flag_buf);
        s_queue.enqueueNDRangeKernel(
            s_assign_max_facet_edge,
            cl::NullRange,
            globalSize,
            localSize
        );
        s_queue.finish();

        s_queue.enqueueCopyBufferToImage
        (flag_buf, flag_img, 0, to_size(0, 0, 0), flag_size);
        s_queue.finish();


        s_assign_max_facet_face.setArg(0, func_img);
        s_assign_max_facet_face.setArg(1, flag_img);
        s_assign_max_facet_face.setArg(2, rct.lo);
        s_assign_max_facet_face.setArg(3, rct.hi);
        s_assign_max_facet_face.setArg(4, ext.lo);
        s_assign_max_facet_face.setArg(5, ext.hi);
        s_assign_max_facet_face.setArg(6, dom.lo);
        s_assign_max_facet_face.setArg(7, dom.hi);
        s_assign_max_facet_face.setArg(8, flag_buf);
        s_queue.enqueueNDRangeKernel(
            s_assign_max_facet_face,
            cl::NullRange,
            globalSize,
            localSize
        );
        s_queue.finish();
        s_queue.enqueueCopyBufferToImage
        (flag_buf, flag_img, 0, to_size(0, 0, 0), flag_size);
        s_queue.finish();


        s_assign_max_facet_cube.setArg(0, func_img);
        s_assign_max_facet_cube.setArg(1, flag_img);
        s_assign_max_facet_cube.setArg(2, rct.lo);
        s_assign_max_facet_cube.setArg(3, rct.hi);
        s_assign_max_facet_cube.setArg(4, ext.lo);
        s_assign_max_facet_cube.setArg(5, ext.hi);
        s_assign_max_facet_cube.setArg(6, dom.lo);
        s_assign_max_facet_cube.setArg(7, dom.hi);
        s_assign_max_facet_cube.setArg(8, flag_buf);
        s_queue.enqueueNDRangeKernel(
            s_assign_max_facet_cube,
            cl::NullRange,
            globalSize,
            localSize
        );
        s_queue.finish();
        s_queue.enqueueCopyBufferToImage
        (flag_buf, flag_img, 0, to_size(0, 0, 0), flag_size);
        s_queue.finish();

        s_assign_pairs.setArg(0, func_img);
        s_assign_pairs.setArg(1, flag_img);
        s_assign_pairs.setArg(2, rct.lo);
        s_assign_pairs.setArg(3, rct.hi);
        s_assign_pairs.setArg(4, ext.lo);
        s_assign_pairs.setArg(5, ext.hi);
        s_assign_pairs.setArg(6, dom.lo);
        s_assign_pairs.setArg(7, dom.hi);
        s_assign_pairs.setArg(8, flag_buf);
        s_queue.enqueueNDRangeKernel(
            s_assign_pairs,
            cl::NullRange,
            globalSize,
            localSize
        );
        s_queue.finish();
        s_queue.enqueueCopyBufferToImage
        (flag_buf, flag_img, 0, to_size(0, 0, 0), flag_size);
        s_queue.finish();


        s_assign_pairs2.setArg(0, func_img);
        s_assign_pairs2.setArg(1, flag_img);
        s_assign_pairs2.setArg(2, rct.lo);
        s_assign_pairs2.setArg(3, rct.hi);
        s_assign_pairs2.setArg(4, ext.lo);
        s_assign_pairs2.setArg(5, ext.hi);
        s_assign_pairs2.setArg(6, dom.lo);
        s_assign_pairs2.setArg(7, dom.hi);
        s_assign_pairs2.setArg(8, flag_buf);
        s_queue.enqueueNDRangeKernel(
            s_assign_pairs2,
            cl::NullRange,
            globalSize,
            localSize
        );
        s_queue.finish();
        s_queue.enqueueCopyBufferToImage
        (flag_buf, flag_img, 0, to_size(0, 0, 0), flag_size);
        s_queue.finish();

        
        s_mark_cps.setArg(0, rct);
        s_mark_cps.setArg(1, ext);
        s_mark_cps.setArg(2, dom);
        s_mark_cps.setArg(3, flag_buf);

      	s_queue.enqueueNDRangeKernel(
            s_mark_cps,
            cl::NullRange,
            globalSize,
            localSize
        );
        

        //s_mark_cps(rct, ext, dom, flag_buf);
        
        rect_list_t bnds;

        get_boundry_rects(from_cell_pair(rct), from_cell_pair(ext), bnds);

        for (uint i = 0; i < bnds.size(); ++i)
        {
            cell_pair_t bnd = to_cell_pair(bnds[i]);
            cell_t  bnd_dir = to_cell(bnds[i].get_normal());

            s_mark_boundry_cps.setArg(0, rct);
            s_mark_boundry_cps.setArg(1, ext);
            s_mark_boundry_cps.setArg(2, dom);
            s_mark_boundry_cps.setArg(3, bnd);
            s_mark_boundry_cps.setArg(4, bnd_dir);
            s_mark_boundry_cps.setArg(5, flag_buf);

            //s_mark_boundry_cps(rct, ext, dom, bnd, bnd_dir, flag_buf);

            s_queue.enqueueNDRangeKernel(
                s_mark_boundry_cps,
                cl::NullRange,
                globalSize,
                localSize
            );
        }

        s_queue.enqueueCopyBufferToImage(flag_buf, flag_img, 0, to_size(0, 0, 0), flag_size);
        s_queue.enqueueReadBuffer(flag_buf, false, 0, cell_ct * sizeof(cell_flag_t), h_flag);

        s_queue.finish();
        

      	/*
        s_assign_max_facet_edge
            (func_img,flag_img,rct.lo,rct.hi,ext.lo,ext.hi,dom.lo,dom.hi,flag_buf);
        s_queue.finish();

        s_queue.enqueueCopyBufferToImage
            (flag_buf,flag_img,0,to_size(0,0,0),flag_size);
        s_queue.finish();

        s_assign_max_facet_face
            (func_img,flag_img,rct.lo,rct.hi,ext.lo,ext.hi,dom.lo,dom.hi,flag_buf);
        s_queue.finish();

        s_queue.enqueueCopyBufferToImage
            (flag_buf,flag_img,0,to_size(0,0,0),flag_size);
        s_queue.finish();

        s_assign_max_facet_cube
            (func_img,flag_img,rct.lo,rct.hi,ext.lo,ext.hi,dom.lo,dom.hi,flag_buf);
        s_queue.finish();

        s_queue.enqueueCopyBufferToImage
            (flag_buf,flag_img,0,to_size(0,0,0),flag_size);
        s_queue.finish();

        s_assign_pairs
            (func_img,flag_img,rct.lo,rct.hi,ext.lo,ext.hi,dom.lo,dom.hi,flag_buf);
        s_queue.finish();

        s_queue.enqueueCopyBufferToImage
            (flag_buf,flag_img,0,to_size(0,0,0),flag_size);
        s_queue.finish();

        s_assign_pairs2
            (func_img,flag_img,rct.lo,rct.hi,ext.lo,ext.hi,dom.lo,dom.hi,flag_buf);
        s_queue.finish();

        s_queue.enqueueCopyBufferToImage
            (flag_buf,flag_img,0,to_size(0,0,0),flag_size);
        s_queue.finish();

//        s_assign_pairs3
//            (func_img,flag_img,rct.lo,rct.hi,ext.lo,ext.hi,dom.lo,dom.hi,flag_buf);
//        s_queue.finish();

//        s_queue.enqueueCopyBufferToImage
//            (flag_buf,flag_img,0,to_size(0,0,0),flag_size);
//        s_queue.finish();

        s_mark_cps(rct,ext,dom,flag_buf);

        rect_list_t bnds;

        get_boundry_rects(from_cell_pair(rct),from_cell_pair(ext),bnds);

        for( uint i = 0 ; i < bnds.size(); ++i)
        {
          cell_pair_t bnd = to_cell_pair(bnds[i]);
          cell_t  bnd_dir = to_cell(bnds[i].get_normal());

          s_mark_boundry_cps(rct,ext,dom,bnd,bnd_dir,flag_buf);
        }

        s_queue.enqueueCopyBufferToImage(flag_buf,flag_img,0,to_size(0,0,0),flag_size);
        s_queue.enqueueReadBuffer(flag_buf,false,0,cell_ct*sizeof(cell_flag_t),h_flag);

        s_queue.finish();
        */


      }
      catch(cl::Error err)
      {
        ENSURES(false)<<"ERROR: "<< err.what()<< "("<< err.err()<< ")"<< std::endl;
        throw;
      }
      
    }

    void __count_and_scan_cps
    ( cell_pair_t rct,
      cell_pair_t ext,
      cell_pair_t dom,
      cl::Image3D &flag_img,
      cl::Buffer  &cp_count_buf,
      int &num_cps)
    {
        
      try
      {
        cp_count_buf = cl::Buffer(s_context,CL_MEM_READ_WRITE,sizeof(int)*WG_SIZE);
        cl::Buffer group_sums_buf(s_context,CL_MEM_READ_WRITE,sizeof(int)*(WG_NUM));
        
        
        //s_count_cps(rct,ext,dom,flag_img,cp_count_buf);


        s_count_cps.setArg(0, rct);
        s_count_cps.setArg(1, ext);
        s_count_cps.setArg(2, dom);
        s_count_cps.setArg(3, flag_img);
        s_count_cps.setArg(4, cp_count_buf);
        
        
        cl::NDRange globalSize(1024); // Total number of work items
        cl::NDRange localSize(256);          // Work items per work group
        //cl::NDRange localSize(32);          // Work items per work group

      	s_queue.enqueueNDRangeKernel(
            s_count_cps,
            cl::NullRange,          
            globalSize,            
            localSize               
        );

        s_queue.finish();
        
         
        rect_list_t bnds;

        get_boundry_rects(from_cell_pair(rct),from_cell_pair(ext),bnds);
        

    	
        for( uint i = 0 ; i < bnds.size(); ++i)
        {
          cell_pair_t bnd = to_cell_pair(bnds[i]);
          cell_t  bnd_dir = to_cell(bnds[i].get_normal());

          s_count_boundry_cps.setArg(0, rct);
          s_count_boundry_cps.setArg(1, ext);
          s_count_boundry_cps.setArg(2, bnd);
          s_count_boundry_cps.setArg(3, bnd_dir);
          s_count_boundry_cps.setArg(4, dom);
          s_count_boundry_cps.setArg(5, flag_img);
          s_count_boundry_cps.setArg(6, cp_count_buf);
          s_queue.enqueueNDRangeKernel(
              s_count_boundry_cps,    
              cl::NullRange,          
              globalSize,            
              localSize               
          );
          //s_count_boundry_cps(rct,ext,dom,bnd,bnd_dir,flag_img,cp_count_buf);
        }
        s_queue.finish();
        
        s_scan_local_sums.setArg(0, cp_count_buf);
        s_scan_local_sums.setArg(1, group_sums_buf);
        s_queue.enqueueNDRangeKernel(
            s_scan_local_sums,
            cl::NullRange,
            globalSize,
            localSize
        );
        s_queue.finish();

        /*  FIXME: enqueuing this kernel leads to future read/writes from the buffer failing
            exception indicates some sort of lifetime issue, a destructor getting called prematurely*/
  
        // s_scan_group_sums.setArg(0, group_sums_buf);
        // s_queue.enqueueNDRangeKernel(
        //     s_scan_group_sums,
        //     cl::NullRange,
        //     globalSize,
        //     localSize
        // );
        
        
        s_scan_update_sums.setArg(0, cp_count_buf);
        s_scan_update_sums.setArg(1, group_sums_buf);
        s_queue.enqueueNDRangeKernel(
            s_scan_update_sums,
            cl::NullRange,
            globalSize,
            localSize
        );

        //s_scan_local_sums(cp_count_buf,group_sums_buf);
        //s_scan_group_sums(group_sums_buf);
        //s_scan_update_sums(cp_count_buf,group_sums_buf);
        
        s_queue.enqueueReadBuffer(group_sums_buf,false,sizeof(int)*(WG_NUM-1),sizeof(int),&num_cps);
        s_queue.finish();

      }
      catch(cl::Error err)
      {
        ENSURES(false)<< "ERROR: "<< err.what()<< "("<< err.err()<< ")"<< std::endl;
        throw;
      }
      
    }

    void __save_cps
      ( cell_pair_t rct,
        cell_pair_t ext,
        cell_pair_t dom,
        cl::Image3D &func_img,
        cl::Image3D &flag_img,
        cl::Buffer  &cp_offset_buf,
        int num_cps,
        cellid_t   *h_cellid,
        cellid_t   *h_vertid,
        int        *h_pair_idx,
        int8_t     *h_index,
        cell_fn_t  *h_func)
    {
      try
      {
        cl::Buffer cp_cellid_buf(s_context,CL_MEM_READ_WRITE,sizeof(cellid_t)*num_cps);
        cl::Buffer cp_vertid_buf(s_context,CL_MEM_READ_WRITE,sizeof(cellid_t)*num_cps);
        cl::Buffer cp_pair_idx_buf(s_context,CL_MEM_READ_WRITE,sizeof(int)*num_cps);
        cl::Buffer cp_index_buf(s_context,CL_MEM_READ_WRITE,sizeof(int8_t)*num_cps);
        cl::Buffer cp_func_buf(s_context,CL_MEM_READ_WRITE,sizeof(cell_fn_t)*num_cps);

        cl::NDRange globalSize(1024); // Total number of work items
        cl::NDRange localSize(256);          // Work items per work group
        //cl::NDRange localSize(32);          // Work items per work group

        s_save_cps.setArg(0, rct);
        s_save_cps.setArg(1, ext);
        s_save_cps.setArg(2, dom);
        s_save_cps.setArg(3, func_img);
        s_save_cps.setArg(4, flag_img);
        s_save_cps.setArg(5, cp_offset_buf);
        s_save_cps.setArg(6, cp_cellid_buf);
        s_save_cps.setArg(7, cp_index_buf);
        s_save_cps.setArg(8, cp_pair_idx_buf);
        s_save_cps.setArg(9, cp_vertid_buf);
        s_save_cps.setArg(10, cp_func_buf);

        s_queue.enqueueNDRangeKernel(
            s_save_cps,
            cl::NullRange,
            globalSize,
            localSize
        );
        s_queue.finish();


        //s_save_cps(rct,ext,dom,func_img,flag_img,cp_offset_buf,cp_cellid_buf,
        //           cp_index_buf,cp_pair_idx_buf,cp_vertid_buf,cp_func_buf);

        rect_list_t bnds;

        get_boundry_rects(from_cell_pair(rct),from_cell_pair(ext),bnds);

        for( uint i = 0 ; i < bnds.size(); ++i)
        {
          cell_pair_t bnd = to_cell_pair(bnds[i]);
          cell_t  bnd_dir = to_cell(bnds[i].get_normal());



          s_save_boundry_cps.setArg(0, rct);
          s_save_boundry_cps.setArg(1, ext);
          s_save_boundry_cps.setArg(2, dom);
          s_save_boundry_cps.setArg(3, bnd);
          s_save_boundry_cps.setArg(4, bnd_dir);
          s_save_boundry_cps.setArg(5, func_img);
          s_save_boundry_cps.setArg(6, flag_img);
          s_save_boundry_cps.setArg(7, cp_offset_buf);
          s_save_boundry_cps.setArg(8, cp_cellid_buf);
          s_save_boundry_cps.setArg(9, cp_index_buf);
          s_save_boundry_cps.setArg(10, cp_pair_idx_buf);
          s_save_boundry_cps.setArg(11, cp_vertid_buf);
          s_save_boundry_cps.setArg(12, cp_func_buf);

          s_queue.enqueueNDRangeKernel(
              s_save_boundry_cps,
              cl::NullRange,
              globalSize,
              localSize
          );

          //s_save_boundry_cps(rct,ext,dom,bnd,bnd_dir,func_img,flag_img,
          //                   cp_offset_buf,cp_cellid_buf,cp_index_buf,
          //                   cp_pair_idx_buf,cp_vertid_buf,cp_func_buf);
        }

        s_queue.finish();

        s_queue.enqueueReadBuffer(cp_cellid_buf,false,0,sizeof(cellid_t)*num_cps,h_cellid);
        s_queue.enqueueReadBuffer(cp_vertid_buf,false,0,sizeof(cellid_t)*num_cps,h_vertid);
        s_queue.enqueueReadBuffer(cp_pair_idx_buf,false,0,sizeof(int)*num_cps,h_pair_idx);
        s_queue.enqueueReadBuffer(cp_index_buf,false,0,sizeof(int8_t)*num_cps,h_index);
        s_queue.enqueueReadBuffer(cp_func_buf,false,0,sizeof(cell_fn_t)*num_cps,h_func);

        s_queue.finish();
      }
      catch(cl::Error err)
      {
        ENSURES(false)<< "ERROR: "<< err.what()<< "("<< err.err()<< ")"<< std::endl;
        throw;
      }
    }

    worker::worker():flag_img(new cl::Image3D){}

    /**
     * \brief Assigns gradient and identifies critical points in the scalar field
     * \param ds dataset with scalar field data
     * \param msc mscomplex data 
     */
    void worker::assign_gradient(dataset_ptr_t ds, mscomplex_ptr_t msc)
    {
      cl::Image3D  func_img;
      cell_pair_t rct = to_cell_pair(ds->m_rect);
      cell_pair_t ext = to_cell_pair(ds->m_ext_rect);
      cell_pair_t dom = to_cell_pair(ds->m_domain_rect);

      __assign_gradient(rct,ext,dom,func_img,*flag_img,
                        //ds->m_vert_fns.data(),ds->m_cell_flags.data());
                        ds->m_vert_fns.data.data(),ds->m_cell_flags.data.data());
      
      if(msc)
      {
        cl::Buffer cp_offset_buf;
        int          num_cps;
		
       __count_and_scan_cps(rct,ext,dom,*flag_img,cp_offset_buf,num_cps);

        msc->resize(num_cps);
        
        __save_cps(rct,ext,dom,func_img,*flag_img,cp_offset_buf,num_cps,
                   msc->m_cp_cellid.data(),msc->m_cp_vertid.data(),
                   msc->m_cp_pair_idx.data(),msc->m_cp_index.data(),
                   msc->m_cp_fn.data());
                   
        //for (int i = 0; i < std::min(10000, msc->get_num_critpts()); ++i) {
            //if(msc->m_cp_cellid[i][0]==256)
            //std::cout << "Critical Point " << i << ": Cell ID = " << msc->m_cp_cellid[i] << std::endl;
        //}
        
      }
      
    }


    void __owner_extrema
    ( cell_pair_t rct,
      cell_pair_t ext,
      cell_pair_t dom,
      cell_pair_t ex_rect,
      cl::Image3D &flag_img,
      int * h_ex_own)
    {
      int num_ex = num_cells2(from_cell_pair(ex_rect));
      try
      {
        cl::Buffer own_buf1(s_context,CL_MEM_READ_WRITE,num_ex*sizeof(int));
        cl::Buffer own_buf2(s_context,CL_MEM_READ_WRITE,num_ex*sizeof(int));
        cl::Buffer is_updated_buf(s_context,CL_MEM_READ_WRITE,sizeof(int));

        cl::NDRange globalSize(1024); // Total number of work items
        cl::NDRange localSize(256);          // Work items per work group

        s_init_propagate.setArg(0, rct);
        s_init_propagate.setArg(1, ext);
        s_init_propagate.setArg(2, dom);
        s_init_propagate.setArg(3, ex_rect);
        s_init_propagate.setArg(4, flag_img);
        s_init_propagate.setArg(5, own_buf1);

        s_queue.enqueueNDRangeKernel(
            s_init_propagate,
            cl::NullRange,
            globalSize,
            localSize
        );

        //s_init_propagate(rct,ext,dom,ex_rect,flag_img,own_buf1);

        int is_updated;

        do
        {
          is_updated = 0;

          s_queue.enqueueWriteBuffer(is_updated_buf,true,0,sizeof(int),&is_updated);

          s_propagate.setArg(0, rct);
          s_propagate.setArg(1, ext);
          s_propagate.setArg(2, dom);
          s_propagate.setArg(3, ex_rect);
          s_propagate.setArg(4, own_buf1);
          s_propagate.setArg(5, own_buf2);
          s_propagate.setArg(6, is_updated_buf);

          s_queue.enqueueNDRangeKernel(
              s_propagate,
              cl::NullRange,
              globalSize,
              localSize
          );
        	//s_propagate(rct,ext,dom,ex_rect,own_buf1,own_buf2,is_updated_buf);


        	s_queue.finish();

          s_queue.enqueueReadBuffer(is_updated_buf,true,0,sizeof(int),&is_updated);
          s_queue.finish();

          std::swap(own_buf1,own_buf2);
        }
        while(is_updated == 1);

        s_queue.enqueueReadBuffer(own_buf1,true,0,num_ex*sizeof(int),h_ex_own);

        s_queue.finish();
      }
      catch(cl::Error err)
      {
        ENSURES(false)<< "ERROR: "<< err.what()<< "("<< err.err()<< ")"<< std::endl;
        throw;
      }
    }

    void worker::owner_extrema(dataset_ptr_t ds)
    {
      cell_pair_t rct      = to_cell_pair(ds->m_rect);
      cell_pair_t ext      = to_cell_pair(ds->m_ext_rect);
      cell_pair_t dom      = to_cell_pair(ds->m_domain_rect);
      cell_pair_t max_rect = to_cell_pair(ds->get_extrema_rect<GDIR_DES>());
      cell_pair_t min_rect = to_cell_pair(ds->get_extrema_rect<GDIR_ASC>());

      //__owner_extrema(rct,ext,dom,max_rect,*flag_img,ds->m_owner_maxima.data());
      __owner_extrema(rct,ext,dom,max_rect,*flag_img,ds->m_owner_maxima.data.data());
     // __owner_extrema(rct,ext,dom,min_rect,*flag_img,ds->m_owner_minima.data());
      __owner_extrema(rct,ext,dom,min_rect,*flag_img,ds->m_owner_minima.data.data());
    }

//    void assign_gradient_and_owner_extrema(dataset_ptr_t ds)
//    {
//      cl::Image3D  flag_img;

//      cell_pair_t rct      = to_cell_pair(ds->m_rect);
//      cell_pair_t ext      = to_cell_pair(ds->m_ext_rect);
//      cell_pair_t dom      = to_cell_pair(ds->m_domain_rect);
//      cell_pair_t max_rect = to_cell_pair(ds->get_extrema_rect(GDIR_DES));
//      cell_pair_t min_rect = to_cell_pair(ds->get_extrema_rect(GDIR_ASC));

//      {
//        cl::Image3D  func_img;

//        __assign_gradient(rct,ext,dom,func_img,flag_img,
//                          ds->m_vert_fns.data(),ds->m_cell_flags.data());
//      }

//      __owner_extrema(rct,ext,dom,max_rect,flag_img,ds->m_owner_maxima.data());
//      __owner_extrema(rct,ext,dom,min_rect,flag_img,ds->m_owner_minima.data());
//    }

    void __update_to_surv_extrema
      ( cell_pair_t rct,
        cell_pair_t ext,
        cell_pair_t dom,
        cell_pair_t ex_rect,
        cl::Buffer &cp_cellid_buf,
        cl::Buffer &surv_cp_no_buf,
        int        num_cps,
        int *      h_ex_own
       )
    {
      int num_ex  = num_cells2(from_cell_pair(ex_rect));
      cl::NDRange globalSize(1024); // Total number of work items
      cl::NDRange localSize(256);          // Work items per work group
      try
      {
        cl::Buffer own_buf1(s_context,CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,num_ex*sizeof(int),h_ex_own);
        cl::Buffer own_buf2(s_context,CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,num_ex*sizeof(int),h_ex_own);


        s_init_update_to_surv_cp_no.setArg(0, rct);
        s_init_update_to_surv_cp_no.setArg(1,ext);
        s_init_update_to_surv_cp_no.setArg(2,dom);
        s_init_update_to_surv_cp_no.setArg(3,cp_cellid_buf);
        s_init_update_to_surv_cp_no.setArg(4,surv_cp_no_buf);
        s_init_update_to_surv_cp_no.setArg(5,num_cps);
        s_init_update_to_surv_cp_no.setArg(6,own_buf1);

        s_queue.enqueueNDRangeKernel(
            s_init_update_to_surv_cp_no,
            cl::NullRange,
            globalSize,
            localSize
        );

        //s_init_update_to_surv_cp_no(rct,ext,dom,ex_rect,cp_cellid_buf,surv_cp_no_buf,num_cps,own_buf1);


        s_update_to_surv_cp_no.setArg(0, rct);
        s_update_to_surv_cp_no.setArg(1, ext);
        s_update_to_surv_cp_no.setArg(2, dom);
        s_update_to_surv_cp_no.setArg(3, ex_rect);
        s_update_to_surv_cp_no.setArg(4, own_buf1);
        s_update_to_surv_cp_no.setArg(5, own_buf2);

        s_queue.enqueueNDRangeKernel(
            s_update_to_surv_cp_no,
            cl::NullRange,
            globalSize,
            localSize
        );


      	//s_update_to_surv_cp_no(rct,ext,dom,ex_rect,own_buf1,own_buf2);



        s_queue.finish();

        s_queue.enqueueReadBuffer(own_buf2,true,0,num_ex*sizeof(int),h_ex_own);
      }
      catch(cl::Error err)
      {
        ENSURES(false)<< "ERROR: "<< err.what()<< "("<< err.err()<< ")"<< std::endl;
        throw;
      }
    }

    void update_to_surv_extrema(dataset_ptr_t ds,mscomplex_ptr_t msc)
    {
      int num_cps = msc->get_num_critpts();

      int_list_t h_surv_cp_no(num_cps);

      #pragma omp parallel for
      for(int i = 0; i < num_cps;++i)
        h_surv_cp_no[i] = (msc->is_extrema(i))?(msc->surv_extrema(i)):(-1);

      cell_pair_t rct = to_cell_pair(ds->m_rect);
      cell_pair_t ext = to_cell_pair(ds->m_ext_rect);
      cell_pair_t dom = to_cell_pair(ds->m_domain_rect);

      cell_pair_t max_rect = to_cell_pair(ds->get_extrema_rect(GDIR_DES));
      cell_pair_t min_rect = to_cell_pair(ds->get_extrema_rect(GDIR_ASC));

      cl::Buffer cp_cellid_buf,surv_cp_no_buf;

      try
      {
        cp_cellid_buf  = cl::Buffer(s_context,CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
                                    num_cps*sizeof(cellid_t),msc->m_cp_cellid.data());

        surv_cp_no_buf = cl::Buffer(s_context,CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
                                    num_cps*sizeof(int),h_surv_cp_no.data());
      }
      catch(cl::Error err)
      {
        ENSURES(false)<< "ERROR: "<< err.what()<< "("<< err.err()<< ")"<< std::endl;
        throw;
      }

      __update_to_surv_extrema(rct,ext,dom,max_rect,cp_cellid_buf,surv_cp_no_buf,
                               //msc->get_num_critpts(),ds->m_owner_maxima.data());
                               msc->get_num_critpts(),ds->m_owner_maxima.data.data());

      __update_to_surv_extrema(rct,ext,dom,min_rect,cp_cellid_buf,surv_cp_no_buf,
                               //msc->get_num_critpts(),ds->m_owner_minima.data());
                               msc->get_num_critpts(),ds->m_owner_minima.data.data());

    }
  }
}

//    void check_assign_gradient_opencl
//      (dataset_ptr_t ds,int dim,rect_t check_rect,cell_flag_t mask)
//    {
//      rect_size_t   span = ds->m_ext_rect.span() + 1;
//      rect_point_t   bl  = ds->m_ext_rect.lower_corner();

//      dataset_t::cellflag_array_t flag(span,boost::fortran_storage_order());

//      flag.reindex(bl);

////      __assign_gradient(ds->m_rect,ds->m_ext_rect,ds->m_domain_rect,
////                              ds->m_vert_fns.data(),flag.data());

//      for(int d = 0 ; d <= dim; ++d)
//      {
//        cellid_t c,s(0,0,0),stride(2,2,2);

//        for(int i = 0 ;  i < d; ++i)
//          s[i] = 1;

//        while(true)
//        {
//          rect_t rect  = rect_t(check_rect.lc()+s,check_rect.uc()-s);

//          int n = c_to_i(rect.uc(),rect,stride) + 1;

//          for( int i = 0; i < n; i ++)
//          {
//            c = i_to_c(i,rect,stride);

//            if((flag(c)&mask) != (ds->m_cell_flags(c)&mask))
//            {
//              cell_flag_t f_gpu =flag(c)&mask;
//              cell_flag_t f_cpu =(ds->m_cell_flags(c)&mask) &mask;

//              cout<<SVAR(c);
//              if(ds->isCellPaired(c))
//              {
//                cout<<"     "<<SVAR(ds->getCellPairId(c));
//              }
//              cout<<endl;

//              cout<<hex<<SVAR(f_cpu)<<" "<<SVAR(f_gpu);
//              cout<<endl;

//            }
//          }

//          if(!next_permutation(s.rbegin(),s.rend()))
//            break;
//        }
//      }
//    }




