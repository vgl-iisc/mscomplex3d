#include <boost/program_options.hpp>

#include <config.h>

#include <grid_dataset.h>
#include <grid_mscomplex.h>
#include <grid_outcore.h>

using namespace std;
using namespace grid;

namespace bpo = boost::program_options ;

namespace grid{namespace opencl{void init();}}

void compute_mscomplex_basic(std::string filename, cellid_t size, double simp_tresh)
{
  g_timer.restart();

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
    msc->simplify(simp_tresh,-1);
    LOG(info) <<"simplification done ------ "<<g_timer.elapsed()<<endl;
  }
  msc->collect_mfolds(ds);

  msc->save(basename+".msc.bin",false);

  LOG(info) <<"write msmfolds done ------ "<<g_timer.elapsed()<<endl;

  LOG(info) <<"------------------------------------"<<endl
            <<"        Finished Processing         "<<endl
            <<"===================================="<<endl;
}


int main(int ac , char **av)
{
  string         filename;
  cellid_t       size;
  double         simp_tresh;
  cellid_t       levels;
  bool           incr_simp;

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

  if(levels == cellid_t::zero)
  {
    compute_mscomplex_basic(filename,size,simp_tresh);
  }
  else
  {
    data_manager_ptr_t gdm(new data_manager_t(filename,size,levels,simp_tresh));

    gdm->work();
  }
}
