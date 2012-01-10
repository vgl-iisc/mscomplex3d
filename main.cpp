#include <exception>
#include <string>

#include <config.h>

#ifdef BUILD_EXEC_GUI
#include <grid_viewer_mainwindow.h>
#endif

#include <grid_datamanager.h>
#include <cpputils.h>

#include <boost/program_options.hpp>

#include <stdexcept>
#include <iostream>

using namespace std;
using namespace grid;

namespace bpo = boost::program_options ;

#ifdef BUILD_EXEC_CUDA
extern "C" void init_cuda();
#endif

#ifdef BUILD_EXEC_OPENCL
namespace grid
{
  namespace opencl
  {
    void init();
  }
}

#endif


int main(int ac , char **av)
{
  string         filename;
  cellid_t       size;
  double         simp_tresh;
  cellid_t       levels;
  bool           incr_simp;

#ifdef BUILD_EXEC_GUI
  bool           gui;
#endif

  bpo::options_description desc("Allowed options");
  desc.add_options()
      ("help,h", "produce help message")
      ("file,f",bpo::value<string >(&filename)->required(),
       "grid file name")
      ("dim,d", bpo::value<cellid_t>(&size)->required(),
       "dim of grid entered as [x,y,z]")
      ("levels,l",bpo::value<cellid_t>(&levels)->default_value(cellid_t(0,0,0)),
       "number of subdivision levels in each dim .. entered as [x,y,z]")
      ("simp-tresh,t",bpo::value<double>(&simp_tresh)->default_value(0.0),
       "simplification treshold")
      ("incr-simp,i", bpo::value<bool>(&incr_simp)->default_value(false),
       "Incrementally simplify the MS complex\n"\
       "simp-tresh is increased in steps of t (-t option) till 1\n"\
       "Results for each step are stored")
#ifdef BUILD_EXEC_GUI
      ("gui,g",bpo::value<bool>(&gui)->default_value(false),"show gui")
#endif
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
    cout<<e.what()<<endl;
    cout<<desc<<endl;
    return 1;
  }

#ifdef BUILD_EXEC_CUDA
  init_cuda();
#endif

#ifdef BUILD_EXEC_OPENCL
    opencl::init();
#endif

  if(incr_simp)
  {
    compute_mscomplex_incr_simp(filename,size,simp_tresh);
  }
  else if(levels == cellid_t::zero)
  {
    compute_mscomplex_basic(filename,size,simp_tresh);
  }
  else
  {
    data_manager_ptr_t gdm(new data_manager_t(filename,size,levels,simp_tresh));

    gdm->work();
  }

#ifdef BUILD_EXEC_GUI
  if(gui)
  {
    QApplication application(ac,av);

    viewer_mainwindow gvmw(gdm);

    gvmw.setWindowTitle("ms complex vis");

    gvmw.show();

    application.exec();
  }
#endif
}
