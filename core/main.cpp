/*=========================================================================

  Program:   mscomplex3d

  Copyright (c) Nithin Shivashankar, Vijay Natarajan

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/

#include <grid_dataset.h>
#include <grid_mscomplex.h>
#include <CL/cl.h>

#define __CL_ENABLE_EXCEPTIONS

#include "grid_dataset_cl.h"
#include "opencl.hpp"
using namespace std;
using namespace grid;
#include <iostream>
#include <iterator>

int main(const int argc, const char *argv[])
{
    if (argc < 5) {
        std::cout << "Please specify a dataset (.raw) file and its dimensions (3 integers)";
        return 0;
    }

    // 0 = GPU, 1 = CPU
    int device = 0; 

    if (argc > 5) {
        device = atoi(argv[5]);
    }

    std::vector<cl::Platform> platforms;
    cl::Platform::get(&platforms);


    if (platforms.empty())
    {
        std::cerr << "No OpenCL platforms found!" << std::endl;
        return -1;
    }

    std::cout << "Available OpenCL platforms:" << std::endl;
    for (size_t i = 0; i < platforms.size(); ++i)
    {
        std::string platformName = platforms[i].getInfo<CL_PLATFORM_NAME>();
        std::cout << "Platform " << i << ": " << platformName << std::endl;

        std::vector<cl::Device> devices;
        platforms[i].getDevices(CL_DEVICE_TYPE_ALL, &devices);

        for (size_t j = 0; j < devices.size(); ++j)
        {
            std::string deviceName = devices[j].getInfo<CL_DEVICE_NAME>();
            cl_device_type deviceType = devices[j].getInfo<CL_DEVICE_TYPE>();
            std::cout << "  Device " << j << ": " << deviceName
                << " [" << ((deviceType == CL_DEVICE_TYPE_CPU) ? "CPU" : "GPU") << "]"
                << std::endl;
        }
    }

//    string         filename = "C:\\Users\\sachi\\OneDrive\\Documents\\PYMS3D_EXAMPLES\\Hydrogen_128x128x128.raw";
    // string         filename = "C:\\Users\\sachi\\OneDrive\\Documents\\PYMS3D_EXAMPLES\\grid_data.raw";
    string filename = argv[1];

    cellid_t       size = cellid_t(atoi(argv[2]),atoi(argv[3]),atoi(argv[4]));
    // cellid_t       size = cellid_t(3,4,5);

    //opencl::init(1);
    get_hw_info(device);
    const rect_t dom(cellid_t::zero, (size - cellid_t::one) * 2);

    DLOG << "Entered :"
        << endl << "\t" << SVAR(filename)
        << endl << "\t" << SVAR(size);

    mscomplex_ptr_t msc(new mscomplex_t);
    dataset_ptr_t ds;

    msc->m_rect = dom;
    msc->m_domain_rect = dom;
    msc->m_ext_rect = dom;

    ds.reset(new dataset_t(dom, dom, dom));

    ds->init(filename);

    auto minmax = std::minmax_element(ds->m_vert_fns.data.begin(), ds->m_vert_fns.data.end());
    std::cout << "Data range: " << *minmax.first << " to " << *minmax.second << std::endl;

    try
    {
        std::cout << "\nOpenCL Context is GPU: " << opencl::is_gpu_context() << std::endl;
        std::cout << "\nOpenCL Context is CPU: " << opencl::is_cpu_context() << std::endl;

        ds->computeMsGraph(msc);

    }
    catch (cl::Error& err)
    {
            std::cerr << "SETUP QUEUE ERROR: " << err.what() << " (" << err.err() << ")" << std::endl;

            if (err.err() == CL_INVALID_PLATFORM)
                std::cerr << "Invalid platform! The selected platform may not support CPU execution." << std::endl;
            else if (err.err() == CL_INVALID_DEVICE)
                std::cerr << "Invalid device! The CPU device may not be available." << std::endl;
            else if (err.err() == CL_INVALID_CONTEXT)
                std::cerr << "Invalid context! OpenCL failed to create a CPU context." << std::endl;
            else if (err.err() == CL_OUT_OF_HOST_MEMORY)
                std::cerr << "Out of host memory! Your system may be running low on RAM." << std::endl;

            throw;
        
	}

    std::vector<int> minima, saddle1, saddle2, maxima;
    std::copy_if(msc->cpno_range().begin(), msc->cpno_range().end(), std::back_insert_iterator(minima), [&](int i) { return !msc->is_canceled(i) && msc->is_index_i_cp_(i, 0); });
    std::copy_if(msc->cpno_range().begin(), msc->cpno_range().end(), std::back_insert_iterator(saddle1), [&](int i) { return !msc->is_canceled(i) && msc->is_index_i_cp_(i, 1); });
    std::copy_if(msc->cpno_range().begin(), msc->cpno_range().end(), std::back_insert_iterator(saddle2), [&](int i) { return !msc->is_canceled(i) && msc->is_index_i_cp_(i, 2); });
    std::copy_if(msc->cpno_range().begin(), msc->cpno_range().end(), std::back_insert_iterator(maxima), [&](int i) { return !msc->is_canceled(i) && msc->is_index_i_cp_(i, 3); });

    std::cout << "minima: " << minima.size() << ", saddle(1): " << saddle1.size() << ", saddle(2) " << saddle2.size() << ", maxima: " << maxima.size() << std::endl;

    minima.clear();
    saddle1.clear();
    saddle2.clear();
    maxima.clear();

    const float simplification_thresh = 0.05f;
    msc->simplify_pers(simplification_thresh);
    std::cout << "after simplification (threshold: " << simplification_thresh << ")" << std::endl;

    std::copy_if(msc->cpno_range().begin(), msc->cpno_range().end(), std::back_insert_iterator(minima), [&](int i) { return !msc->is_canceled(i) && msc->is_index_i_cp_(i, 0); });
    std::copy_if(msc->cpno_range().begin(), msc->cpno_range().end(), std::back_insert_iterator(saddle1), [&](int i) { return !msc->is_canceled(i) && msc->is_index_i_cp_(i, 1); });
    std::copy_if(msc->cpno_range().begin(), msc->cpno_range().end(), std::back_insert_iterator(saddle2), [&](int i) { return !msc->is_canceled(i) && msc->is_index_i_cp_(i, 2); });
    std::copy_if(msc->cpno_range().begin(), msc->cpno_range().end(), std::back_insert_iterator(maxima), [&](int i) { return !msc->is_canceled(i) && msc->is_index_i_cp_(i, 3); });

    std::cout << "minima: " << minima.size() << ", saddle(1): " << saddle1.size() << ", saddle(2) " << saddle2.size() << ", maxima: " << maxima.size() << std::endl;
    
    return 0;

}
