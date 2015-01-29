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

#include <boost/program_options.hpp>

#include <grid_dataset.h>
#include <grid_mscomplex.h>
//#include <grid_outcore.h>

#include <sys/resource.h>

using namespace std;
using namespace grid;

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


int main(int ac , char **av)
{
  string         filename;
  cellid_t       size;
  double         simp_tresh;
  cellid_t       levels;

  {
    struct rlimit rl;

    long lim = 1024*1024;
    lim *= 1024*4;

    rl.rlim_cur = lim;
    rl.rlim_max = lim;
    setrlimit (RLIMIT_AS, &rl);
  }


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

//  if(levels == cellid_t::zero)
  {
    compute_mscomplex_basic(filename,size,simp_tresh);
  }
//  else
//  {
//    data_manager_ptr_t gdm(new data_manager_t(filename,size,levels,simp_tresh));

//    gdm->work();
//  }
}
