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

//#include <boost/program_options.hpp>

#include <grid_dataset.h>
#include <grid_mscomplex.h>
#include <CL/cl.h>

#define __CL_ENABLE_EXCEPTIONS

#include "grid_dataset_cl.h"
#include "opencl.hpp"
//#include <grid_outcore.h>

//#include <sys/resource.h>


using namespace std;
using namespace grid;
/*
namespace bpo = boost::program_options ;

void compute_mscomplex_basic(std::string filename, cellid_t size, double simp_tresh)
{
  g_timer.restart();

  LOG(info) << get_hw_info();

  LOG(info) <<"===================================="<<endl
            <<"         Starting Processing        "<<endl
            <<"------------------------------------"<<endl;

  rect_t d(cellid_t::zero,(size-cellid_t::one)*2);
  dataset_ptr_t   ds (new dataset_t(d,d,d));
  mscomplex_ptr_t msc(new mscomplex_t(d,d,d));

  string basename(filename);

  int ext_pos = basename.size() -4;

  if(ext_pos >=0 && basename.substr(ext_pos,4) == ".raw")
    basename = basename.substr(0,ext_pos);

  ds->init(filename);
  LOG(info) <<"data read ---------------- "<<g_timer.elapsed()<<endl;

  ds->computeMsGraph(msc);
  LOG(info) <<"msgraph done ------------- "<<g_timer.elapsed()<<endl;

  if(simp_tresh >=0)
  {
    msc->simplify_pers(simp_tresh);
    LOG(info) <<"simplification done ------ "<<g_timer.elapsed()<<endl;
  }
  msc->collect_mfolds(ds);

  LOG(info) <<"mfold collection done ---- "<<g_timer.elapsed()<<endl;

  msc->save(basename+".msc");

  LOG(info) <<"write data done ---------- "<<g_timer.elapsed()<<endl;

  LOG(info) <<"------------------------------------"<<endl
            <<"        Finished Processing         "<<endl
            <<"===================================="<<endl;
}

*/



//#include <CL/cl.hpp>
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
    string         filename = "C:\\Users\\sachi\\OneDrive\\Documents\\PYMS3D_EXAMPLES\\grid_data.raw";

    //cellid_t       size = cellid_t(128,128,128);
    cellid_t       size = cellid_t(3,4,5);

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

/*
int main(int ac , char **av)
{
  //string         filename = "C:\\Users\\sachi\\OneDrive\\Documents\\PYMS3D_EXAMPLES\\Hydrogen_128x128x128.raw";
  //string         filename = "C:\\Users\\sachi\\OneDrive\\Documents\\PYMS3D_EXAMPLES\\SquareCylinderOkuboWeiss_t0816.raw";
  string         filename = "C:\\Users\\sachi\\OneDrive\\Documents\\PYMS3D_EXAMPLES\\neghip_64x64x64_uint8.raw";
  //string         filename = "C:\\Users\\sachi\\OneDrive\\Documents\\PYMS3D_EXAMPLES\\grid_data.raw";
  //cellid_t       size = cellid_t(192,64,48);
  cellid_t       size = cellid_t(64,64,64);
  //cellid_t       size = cellid_t(128,128,128);
  //cellid_t       size = cellid_t(129, 4, 4);
  //double         simp_tresh;
  //cellid_t       levels;

  opencl::init();


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

  ds->computeMsGraph(msc);

  //msc->simplify_pers(0.05);




    /*
  int num_critpts = msc->get_num_critpts();
  std::ranges::iota_view<int, int> a= std::views::iota(0, num_critpts);


  auto rng = a | std::ranges::views::filter([this](const auto& item) {
      return this->is_not_canceled(item);
          });
          */
  

  //std::cout << "\nRNG:\n";
  //for (const auto& element : rng) {
  //    std::cout << element << std::endl;
  //}

  /*
  TLOG << "Computed:";
  int i = 0;

      auto filtered_rng = std::ranges::views::filter(rng, [this, i](int x) {
          return this->is_index_i_cp_(x, i);
          });
      std::vector<int> filtered_vector(filtered_rng.begin(), filtered_rng.end());


    for(int l=0;l<filtered_vector.size();i++)
      std::cout <<"\n" << filtered_vector[i];
  
  */

  //ds->computeMsGraph()

  //msc->simplify_pers(0.05);

  

  


  /*
  {
    struct rlimit rl;

    long lim = 1024*1024;
    lim *= 1024*4;

    rl.rlim_cur = lim;
    rl.rlim_max = lim;
    setrlimit (RLIMIT_AS, &rl);
  }
  */

/*
  bpo::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "produce help message")
      ("file,f",bpo::value<string >(&filename)->required(),
       "grid file name")
      ("dim,d", bpo::value<cellid_t>(&size)->required(),
       "dim of grid entered as [x,y,z]")
//      ("levels,l",bpo::value<cellid_t>(&levels)->default_value(cellid_t(0,0,0)),
//       "number of subdivision levels in each dim .. entered as [x,y,z]")
      ("simp-tresh,t",bpo::value<double>(&simp_tresh)->default_value(0.0),
       "simplification treshold")
      ;

  bpo::variables_map vm;
  bpo::store(bpo::parse_command_line(ac, av, desc), vm);

  if (vm.count("help"))
  {
    cout << desc << endl;
    return 0;
  }
  try
  {
    bpo::notify(vm);
  }
  catch(bpo::required_option e)
  {
    LOG(fatal) <<e.what()<<endl;
    LOG(fatal) <<desc<<endl;
    return 1;
  }

  opencl::init();
*/
//  if(levels == cellid_t::zero)
  //{
    //compute_mscomplex_basic(filename,size,simp_tresh);
  //}
//  else
//  {
//    data_manager_ptr_t gdm(new data_manager_t(filename,size,levels,simp_tresh));

//    gdm->work();
//  }
//}



/*
#include <iostream>
#include <string>
#include <filesystem>
#include <optional>
#include <stdexcept>

//namespace fs = std::filesystem;

void parse_arguments(int argc, char* argv[], std::string& filename, cellid_t& size, double& simp_tresh)
{
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if (arg == "--help" || arg == "-h")
        {
            std::cout << "Usage: " << argv[0] << " [options]" << std::endl
                << "--file,-f  : grid file name (required)" << std::endl
                << "--dim,-d   : dim of grid entered as [x,y,z] (required)" << std::endl
                << "--simp-tresh,-t : simplification threshold (default: 0.0)" << std::endl;
            std::exit(0);
        }
        else if (arg == "--file" || arg == "-f")
        {
            if (i + 1 < argc)
            {
                filename = argv[++i];
            }
            else
            {
                throw std::invalid_argument("Missing filename argument.");
            }
        }
        else if (arg == "--dim" || arg == "-d")
        {
            if (i + 1 < argc)
            {
                // Assuming cellid_t can be constructed from string, adjust as needed
                size = std::stoi(argv[++i]);
            }
            else
            {
                throw std::invalid_argument("Missing grid dimension argument.");
            }
        }
        else if (arg == "--simp-tresh" || arg == "-t")
        {
            if (i + 1 < argc)
            {
                simp_tresh = std::stod(argv[++i]);
            }
            else
            {
                throw std::invalid_argument("Missing simplification threshold argument.");
            }
        }
    }

    // Validate required parameters
    if (filename.empty() || size == cellid_t(0))
    {
        throw std::invalid_argument("Filename and size are required.");
    }
}

void compute_mscomplex_basic(const std::string& filename, cellid_t size, double simp_tresh)
{
    g_timer.restart();

    LOG(info) << get_hw_info();
    LOG(info) << "====================================" << std::endl
        << "         Starting Processing        " << std::endl
        << "------------------------------------" << std::endl;

    rect_t d(cellid_t::zero, (size - cellid_t::one) * 2);
    dataset_ptr_t ds(new dataset_t(d, d, d));
    mscomplex_ptr_t msc(new mscomplex_t(d, d, d));

    std::string basename(filename);
    int ext_pos = basename.size() - 4;
    if (ext_pos >= 0 && basename.substr(ext_pos, 4) == ".raw")
        basename = basename.substr(0, ext_pos);

    ds->init(filename);
    LOG(info) << "data read ---------------- " << g_timer.elapsed() << std::endl;

    ds->computeMsGraph(msc);
    LOG(info) << "msgraph done ------------- " << g_timer.elapsed() << std::endl;

    if (simp_tresh >= 0)
    {
        msc->simplify_pers(simp_tresh);
        LOG(info) << "simplification done ------ " << g_timer.elapsed() << std::endl;
    }
    msc->collect_mfolds(ds);

    LOG(info) << "mfold collection done ---- " << g_timer.elapsed() << std::endl;

    msc->save(basename + ".msc");

    LOG(info) << "write data done ---------- " << g_timer.elapsed() << std::endl;

    LOG(info) << "------------------------------------" << std::endl
        << "        Finished Processing         " << std::endl
        << "====================================" << std::endl;
}

int main(int argc, char* argv[])
{
    std::string filename;
    cellid_t size;
    double simp_tresh = 0.0;

    try
    {
        parse_arguments(argc, argv, filename, size, simp_tresh);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    opencl::init();

    compute_mscomplex_basic(filename, size, simp_tresh);

    return 0;
}
*/