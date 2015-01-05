#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>

#include <iostream>


#include <grid_dataset.h>
#include <grid_mscomplex.h>

using namespace std;
using namespace boost::python;
using namespace grid;
using namespace utl;

namespace bp = boost::python;


/*****************************************************************************/
namespace pyms3d {


/*****************************************************************************/
/******** Python wrapped Morse-Smale complex class                    ********/
/*****************************************************************************/

// Wrapper to hold on to references to some other objects
class mscomplex_pyms3d_t: public mscomplex_t
{
public:
  dataset_ptr_t     ds;
};

typedef  boost::shared_ptr<mscomplex_pyms3d_t> mscomplex_pyms3d_ptr_t;

mscomplex_pyms3d_ptr_t new_msc()
{
  mscomplex_pyms3d_ptr_t msc(new mscomplex_pyms3d_t);
  return msc;
}

void mscomplex_compute_bin
(mscomplex_pyms3d_ptr_t msc,
 cellid_t s,
 std::string bin_file,
 std::string bin_fmt="float32")
{
  ENSURES(bin_fmt == "float32") <<"Only float32 is supported" <<endl;

  rect_t r(cellid_t::zero,(s-cellid_t::one)*2);

  msc->ds.reset(new dataset_t(r,r,r));
  msc->ds->init(bin_file);
  msc->ds->computeMsGraph(msc);
}


int mscomplex_num_canc(mscomplex_pyms3d_ptr_t msc)
{
  return msc->m_canc_list.size();
}

bp::tuple mscomplex_canc(mscomplex_pyms3d_ptr_t msc,int i)
{
  ASSERT(is_in_range(i,0,msc->m_canc_list.size()));
  return bp::make_tuple(msc->m_canc_list[i][0],msc->m_canc_list[i][1]);
}

bp::tuple mscomplex_frange(mscomplex_pyms3d_ptr_t msc)
{
//  return bp::make_tuple(msc->m_fmin,msc->m_fmax);
  LOG(error) << "Not yet implemented" << endl;
  return bp::make_tuple();
}

template <eGDIR dir>
bp::list mscomplex_conn(mscomplex_pyms3d_ptr_t msc, int i)
{
  ASSERT(is_in_range(i,0,msc->get_num_critpts()));
  bp::list r;
  BOOST_FOREACH(int_int_t c,msc->m_conn[dir][i])
  {
    r.append(c);
  }
  return r;
}

bp::list mscomplex_cps(mscomplex_pyms3d_ptr_t msc)
{
  bp::list r;

  for(int i = 0 ; i < msc->get_num_critpts(); ++i)
    if(msc->is_not_paired(i))
      r.append(i);

  return r;
}

//void mscomplex_gen_pers_pairs(mscomplex_pyms3d_ptr_t msc)
//{
//  int lastVer = msc->get_multires_version();
//  msc->simplify(1.0,true);
//  msc->set_multires_version(lastVer);
//}

void mscomplex_collect_mfolds(mscomplex_pyms3d_ptr_t msc)
{
  ENSURES(msc->ds !=0)
      << "Gradient information unavailable" <<endl
      << "Did you load the mscomplex from a file!!!"<<endl;

  msc->collect_mfolds(msc->ds);
}

template <eGDIR dir>
bp::list mscomplex_geom(mscomplex_pyms3d_ptr_t msc, int i)
{
  ASSERT(is_in_range(i,0,msc->get_num_critpts()));
  bp::list r;
  BOOST_FOREACH(cellid_t c,msc->m_mfolds[dir][i])
  {
    r.append(c);
  }
  return r;
}

//bp::list mscomplex_arc_geom(mscomplex_ptr_t msc, int a, int b)
//{
//  bp::list r;

//  if(msc->m_arc_geom.count(int_pair_t(a,b)) != 0)
//  {
//      BOOST_FOREACH(int c,msc->m_arc_geom[int_pair_t(a,b)])
//      {
//        r.append(c);
//      }
//    }

//  return r;
//}


void wrap_mscomplex_t()
{

  docstring_options local_docstring_options(true, false, false);


  class_<mscomplex_pyms3d_t,mscomplex_pyms3d_ptr_t>
      ("mscomplex","The Morse-Smale complex object",no_init)
      .def("__init__", make_constructor( &new_msc),
           "ctor")
      .def("num_cp",&mscomplex_t::get_num_critpts,
           "Number of Critical Points")
      .def("fn",&mscomplex_t::fn,
           "Function value at critical point i")
      .def("index",&mscomplex_t::index,
           "Morse index od critical point i")
      .def("pair_idx",&mscomplex_t::pair_idx,
           "Index of the cp that is paired with i (-1 if it is not paired)")
//      .def("is_boundry",&mscomplex_t::is_boundry,
//           "If the cp is on the boundary or not")
      .def("vertid",&mscomplex_t::vertid,
           "vertex id of maximal vertex of critical cell")
      .def("cellid",&mscomplex_t::cellid,
           "cell id of critical cell")
      .def("load",&mscomplex_t::load,
           "Load mscomplex from file")
      .def("save",&mscomplex_t::save,
           "Save mscomplex to file")
      .def("num_canc",&mscomplex_num_canc,
           "Number of cancellation pairs")
      .def("canc",&mscomplex_canc,
           "The ith cancellation pair")
      .def("frange",&mscomplex_frange,
           "Range of function values")
      .def("asc",&mscomplex_conn<ASC>,
           "List of ascending cps connected to a given critical point i")
      .def("des",&mscomplex_conn<DES>,
           "List of descending cps connected to a given critical point i")
      .def("asc_geom",&mscomplex_geom<ASC>,
           "Ascending manifold geometry of a given critical point i")
      .def("des_geom",&mscomplex_geom<DES>,
           "Descending manifold geometry of a given critical point i")
      .def("cps",&mscomplex_cps,
           "Returns a list of surviving critical cps")
//      .def("gen_pers_hierarchy",&mscomplex_gen_pers_pairs,
//           "Generates the persistence hierarchy using topo simplification")
      .def("compute_off_bin",&mscomplex_compute_bin,
           "Compute the Mscomplex from a triagulation given in the .off format\n"\
           "and scalar function from the given bin file in the given bin format\n"\
           "\n"\
           "Parameters: \n"\
           "    off_file: the off file containing the triangulation.\n"\
           "    bin_file: the bin file containing the scalar function.\n"\
           "    bin_fmt: binary format .\n"\
           "             Acceptable values = (\"float32\",\"float64\")\n"
           "\n"\
           "Note: This only computes the combinatorial structure\n"\
           "     Call collect_mfold(s) to extract geometry\n"
           )
      .def("collect_geom",&mscomplex_collect_mfolds,
           "Collect the geometry of all survivng critical points\n"\
           "\n"\
           "Note: This must be called only after any of the compute functions are called. \n"\
           )

      .def("simplify_pers",&mscomplex_t::simplify,
           "Simplify the Morse-Smale complex using topological persistence\n"\
           "Parameters:\n"\
           "    tresh: persistence threshold\n"\
//           "    is_nrm: is the threshold normalized to [0,1] or not.\n"\
//           "            if not then thresh is in same scale as input function\n"\
//           "    req_nmax,req_nmin: num maxima/minima that should be retained\n"\
//           "                       set to 0 to ignore"
           )


//      .def("arc_geom",&mscomplex_arc_geom,
//           "seqence of cells of arc connecting 2 cps (empty list if no arc exists)")
      ;
}

struct cellid_to_tup
{
  static PyObject* convert(cellid_t v)
  {return boost::python::incref(bp::make_tuple(v[0],v[1],v[2]).ptr());}
};

struct int_int_to_tup
{
  static PyObject* convert(int_int_t v)
  {return boost::python::incref(bp::make_tuple(v.first,v.second).ptr());}
};


/*****************************************************************************/
/******** Define the pyms3d module                                    ********/
/*****************************************************************************/

BOOST_PYTHON_MODULE(pyms3d)
{
  numeric::array::set_module_and_type("numpy", "ndarray");

  boost::python::to_python_converter<cellid_t,cellid_to_tup>();
  boost::python::to_python_converter<int_int_t,int_int_to_tup>();

  wrap_mscomplex_t();
}

}/****************************************************************************/

