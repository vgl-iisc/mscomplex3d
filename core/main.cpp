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

int main()
{
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
    string filename = "Debug\\Hydrogen_128x128x128.raw";

    cellid_t       size = cellid_t(128,128,128);
    // cellid_t       size = cellid_t(3,4,5);

    //opencl::init(1);
    get_hw_info(0);
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
    return 0;

}
