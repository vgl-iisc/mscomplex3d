/*=========================================================================

  Program:   pyms3d

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

#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/numpy.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/iterator_range.hpp>

#include <iostream>


#include <grid_dataset.h>
#include <grid_mscomplex.h>

using namespace std;
using namespace boost::python;
using namespace grid;
using namespace utl;

namespace bp = boost::python;
namespace np = boost::numpy;
namespace br = boost::range;
namespace ba = boost::adaptors;

namespace utl {
template<class rng_t>
inline std::string to_string_rng(rng_t rng,const char * delim = ", ")
{
  std::stringstream ss;

  BOOST_AUTO(b,boost::begin(rng));
  BOOST_AUTO(e,boost::end(rng));

  ss<<"["; for(; b!=e ; b) ss << *b++ << delim; ss<<"]";
  return ss.str();
}
}

template<typename DTYPE,int NCOL,typename T>
np::ndarray vector_to_ndarray(const std::vector<T> & vec)
{
  BOOST_STATIC_ASSERT(sizeof(DTYPE)*NCOL==sizeof(T));

  int           N = vec.size();
  np::dtype    dt = np::dtype::get_builtin<DTYPE>();
  bp::tuple   dim = (NCOL==1)?(bp::make_tuple(N)):(bp::make_tuple(N,NCOL));
  np::ndarray arr = np::zeros(dim,dt);
  br::copy(vec,reinterpret_cast<T*>(arr.get_data()));
  return arr;
}

template<typename DTYPE,int NCOL,typename rng_t>
np::ndarray range_to_ndarray(const rng_t & rng)
{
  typedef typename boost::range_value<rng_t>::type    T;

  BOOST_STATIC_ASSERT(sizeof(DTYPE)*NCOL==sizeof(T));

  size_t        N = boost::distance(rng);
  np::dtype    dt = np::dtype::get_builtin<DTYPE>();
  bp::tuple   dim = (NCOL==1)?(bp::make_tuple(N)):(bp::make_tuple(N,NCOL));
  np::ndarray arr = np::zeros(dim,dt);
  char       *ptr = arr.get_data();
  BOOST_AUTO(iter,boost::begin(rng));

  for(;iter!=boost::end(rng);)
  {
    T v  = *iter++;
    memcpy((void*)ptr,(void*)(&v),sizeof(T));
    ptr += sizeof(T);
  }
  return arr;
}



/*===========================================================================*/
namespace pyms3d {
/// \brief Wrapper for python
class mscomplex_pyms3d_t: public mscomplex_t
{
public:
  dataset_ptr_t     ds;
  typedef mscomplex_t base_t;

  /*-------------------------------------------------------------------------*/

  void save(const string &f) const
  {
    DLOG << "Entered :" << SVAR(f);
    std::fstream fs(f.c_str(),std::ios::out|std::ios::binary);
    ENSUREV(fs.is_open(),"file not found!!",f);
    save_bin(fs);
    ds->save_bin(fs);
    DLOG << "Exited  :";
  }

  /*-------------------------------------------------------------------------*/

  void load(const string &f)
  {
    DLOG << "Entered :" << SVAR(f);
    std::fstream fs(f.c_str(),std::ios::in|std::ios::binary);
    ENSUREV(fs.is_open(),"file not found!!",f);
    load_bin(fs);
    ds->load_bin(fs);
    DLOG << "Exited  :";
  }

  /*-------------------------------------------------------------------------*/

  template <eGDIR dir>
  np::ndarray conn(int cp)
  {
    ENSURES(is_in_range(cp,0,get_num_critpts())) << "out of range "<<SVAR(cp);
    TLOG <<SVAR(cp); return range_to_ndarray<int,2>(m_conn[dir][cp]);
  }

  /*-------------------------------------------------------------------------*/

  template <eCCTYPE ccTYPE,int dim>
  np::ndarray __mfold_to_point_indices(const mfold_t & mfold)
  {
    const int nCellPoints = (ccTYPE == CC_PRIM)?(1<<(dim)):(1<<(gc_grid_dim - dim));
    const int           O = (ccTYPE == CC_PRIM)?(0):(1);

    rect_t prect = (ccTYPE == CC_PRIM)?(m_rect):(rect_t(m_rect.lc()+1,m_rect.uc()-1));

    int_list_t points;
    int_list_t point_idxs(mfold.size(),-1);

    int num_pts;

    #pragma omp parallel for
    for(int i = 0; i < mfold.size(); ++i)
    {
      cellid_t c = mfold[i],j;

      for(    j[2] = -((c[2]+O)&1) ; j[2] <= ((c[2]+O)&1) ;j[2]+=2)
        for(  j[1] = -((c[1]+O)&1) ; j[1] <= ((c[1]+O)&1) ;j[1]+=2)
          for(j[0] = -((c[0]+O)&1) ; j[0] <= ((c[0]+O)&1) ;j[0]+=2)
            if(!prect.contains(c+j))
              point_idxs[i] = __sync_fetch_and_add(&num_pts,1);
    }

    points.resize(num_pts);

    #pragma omp parallel for
    for(int i = 0; i < mfold.size(); ++i)
    {
      cellid_t c = mfold[i],j;
      if(point_idxs[i] != -1)
        for(    j[2] = -((c[2]+O)&1) ; j[2] <= ((c[2]+O)&1) ;j[2]+=2)
          for(  j[1] = -((c[1]+O)&1) ; j[1] <= ((c[1]+O)&1) ;j[1]+=2)
            for(j[0] = -((c[0]+O)&1) ; j[0] <= ((c[0]+O)&1) ;j[0]+=2)
              points[point_idxs[i]] =  c_to_i2(prect,c+j);
    }
    TLOG << "Transformed:";

    np::ndarray arr = vector_to_ndarray<int,1>(points);
    if (nCellPoints != 1)
      arr = arr.reshape(bp::make_tuple(points.size()/nCellPoints,nCellPoints));
    TLOG << "Copied  :";

    return arr;
  }

  /// \brief get the geometry of a critical point cp
  ///
  /// \param cp     id of the critical point
  /// \param hver   hierarchical version of the geometry (-1 indicates current)
  /// \param dir    Asc/Des geometry
  /// \param ToPts  convert the geometry to point indices.
  ///               False  : returns a list of cellids
  ///               True   : returns a list of Point idxs in Primal/Dual grid
  ///
  /// \note : if ToPts is enabled, then the cellids are converted to point
  ///          indices in primal/Dual Grid according to the following table.
  ///
  ///         DIR     index(cp)    pt-Index-Type    ret-ArrayShape
  ///         ASC         0         Primal             [NC]
  ///         ASC         1         Dual               [NC',4]
  ///         ASC         2         Dual               [NC',2]
  ///         DES         1         Primal             [NC,2]
  ///         DES         2         Primal             [NC,4]
  ///         DES         3         Dual               [NC]
  ///
  ///     Here NC  #cells in Asc/Des mfold of cp.
  ///          NC' #cells in Asc/Des mfold of cp whose all dual pts are in Primal Grid

  template <eGDIR dir> np::ndarray geom(int cp,int hver=-1,bool ToPts=true)
  {
    TLOG <<"Entered :" << SVAR(cp) << SVAR(hver);

    if( hver == -1) hver = get_hversion();

    ENSURES(is_in_range(cp,0,get_num_critpts()))
        << "out of range "<<SVAR(cp);
    ENSURES(is_in_range(hver,0,m_canc_list.size()+1))
        << "hversion not in range "<<SVAR(hver);
    BOOST_STATIC_ASSERT(   dir==ASC     ||    dir==DES );

    int_list_t l;

    int dim = index(cp);

    m_merge_dag->get_contrib_cps
        (l,dir,cp,hver,m_geom_hversion[dir][dim]);

    mfold_t mfold;

    for(int j = 0 ; j < l.size(); ++j)
      br::copy(m_mfolds[dir][l[j]],std::back_inserter(mfold));

    TLOG << "Computed:" << SVAR(mfold.size());

    if (ToPts)
    {
      if(dir==ASC && dim==0) return __mfold_to_point_indices<CC_PRIM,0>(mfold);
      if(dir==ASC && dim==1) return __mfold_to_point_indices<CC_DUAL,1>(mfold);
      if(dir==ASC && dim==2) return __mfold_to_point_indices<CC_DUAL,2>(mfold);

      if(dir==DES && dim==1) return __mfold_to_point_indices<CC_PRIM,1>(mfold);
      if(dir==DES && dim==2) return __mfold_to_point_indices<CC_PRIM,2>(mfold);
      if(dir==DES && dim==3) return __mfold_to_point_indices<CC_DUAL,3>(mfold);

      ENSURES(false) <<"Should never reach here";
    }

    return vector_to_ndarray<cell_coord_t,gc_grid_dim>(mfold);
  }

  /*-------------------------------------------------------------------------*/

  template <eCCTYPE pTYPE> np::ndarray points()
  {
    rect_t prect = m_rect;

    if(pTYPE == CC_DUAL)
      prect = rect_t(m_rect.lc()+1,m_rect.uc()-1);

    int           N = 1 + c_to_i2(prect,prect.uc());
    np::dtype    dt = np::dtype::get_builtin<float>();
    np::ndarray arr = np::zeros(bp::make_tuple(N,gc_grid_dim),dt);
    float    *  iter= reinterpret_cast<float*>(arr.get_data());

    BOOST_STATIC_ASSERT(gc_grid_dim == 3 && "defined for 3-manifolds only");

    cellid_t c;

    for(    c[2] = prect.lc()[2] ; c[2] <= prect.uc()[2] ;c[2]+=2)
      for(  c[1] = prect.lc()[1] ; c[1] <= prect.uc()[1] ;c[1]+=2)
        for(c[0] = prect.lc()[0] ; c[0] <= prect.uc()[0] ;c[0]+=2)
        {
          *iter++ = float(c[0])/2;
          *iter++ = float(c[1])/2;
          *iter++ = float(c[2])/2;
        }

    return arr;
  }

  /*-------------------------------------------------------------------------*/

  np::ndarray cps(int i)
  {
    TLOG << "Entered :";

    BOOST_AUTO(rng,cpno_range()|ba::filtered
                    (bind(&mscomplex_pyms3d_t::is_not_canceled,this,_1)));

    TLOG << "Computed:";

    if(i == -1) return range_to_ndarray<int,1>(rng);
    else        return range_to_ndarray<int,1>
        (rng|ba::filtered(bind(&mscomplex_pyms3d_t::is_index_i_cp_,this,_1,i)));

    TLOG << "Exited  :";

  }

  /*-------------------------------------------------------------------------*/

  np::ndarray canc_pairs() {TLOG;return vector_to_ndarray<int,2>         (m_canc_list);}
  np::ndarray cps_func()   {TLOG;return vector_to_ndarray<cell_fn_t,1>   (m_cp_fn);}
  np::ndarray cps_index()  {TLOG;return vector_to_ndarray<int8_t,1>      (m_cp_index);}
  np::ndarray cps_pairid() {TLOG;return vector_to_ndarray<int,1>         (m_cp_pair_idx);}
  np::ndarray cps_cellid() {TLOG;return vector_to_ndarray<cell_coord_t,3>(m_cp_cellid);}
  np::ndarray cps_vertid() {TLOG;return vector_to_ndarray<cell_coord_t,3>(m_cp_vertid);}

  /*-------------------------------------------------------------------------*/

  void collect_mfolds(int dir ,int dim)
  {
    DLOG << "Entered :" << SVAR(dir) << SVAR(dim);
    ENSURES(ds !=0)<< "Gradient information unavailable";
    base_t::collect_mfolds((eGDIR)dir,dim,ds);
    DLOG << "Exited  :";
  }
  /*-------------------------------------------------------------------------*/
};

void ping()
{
  ILOG <<"Pinged  :";
}
}

/*===========================================================================*/

namespace pyms3d {

typedef  boost::shared_ptr<mscomplex_pyms3d_t> mscomplex_pyms3d_ptr_t;

mscomplex_pyms3d_ptr_t new_msc()
{
  mscomplex_pyms3d_ptr_t msc(new mscomplex_pyms3d_t);
  return msc;
}

void mscomplex_compute_bin
(mscomplex_pyms3d_ptr_t msc, 
 std::string            bin_file,
 bp::tuple              tp)
{
//  ENSURES(bin_fmt == "float32") <<"Only float32 is supported" <<endl;
  int x = bp::extract<int>(tp[0]);
  int y = bp::extract<int>(tp[1]);
  int z = bp::extract<int>(tp[2]);
  cellid_t size = cellid_t(x,y,z);

  rect_t dom(cellid_t::zero,(size-cellid_t::one)*2);

  DLOG << "Entered :"
       << endl << "\t" << SVAR(bin_file)
       << endl << "\t" << SVAR(size);

  msc->m_rect        = dom;
  msc->m_domain_rect = dom;
  msc->m_ext_rect    = dom;

  msc->ds.reset(new dataset_t(dom,dom,dom));
  msc->ds->init(bin_file);
  msc->ds->computeMsGraph(msc);

  DLOG << "Exited  :";
}


void wrap_mscomplex_t()
{
  docstring_options local_docstring_options(true, false, false);

  class_<mscomplex_pyms3d_t,mscomplex_pyms3d_ptr_t>
      ("mscomplex","The Morse-Smale complex object",no_init)
      .def("__init__", make_constructor( &new_msc),
           "ctor")

      .def("num_cps",&mscomplex_t::get_num_critpts,
           "Number of Critical Points")
      .def("cps",&mscomplex_pyms3d_t::cps,(bp::arg("dim")=-1),
           "Returns a list of surviving critical cps\n"\
           "\n"
           "Parameters   :\n"
           "          dim: index of the cps. -1 signifies all. default=-1\n"
           )


      .def("cp_func",&mscomplex_t::fn,
           "Function value at critical point i")
      .def("cp_index",&mscomplex_t::index,
           "Morse index od critical point i")
      .def("cp_pairid",&mscomplex_t::pair_idx,
           "Index of the cp that is paired with i (-1 if it is not paired)")
//      .def("is_boundry",&mscomplex_t::is_boundry,
//           "If the cp is on the boundary or not")
      .def("cp_vertid",&mscomplex_t::vertid,
           "vertex id of maximal vertex of critical cell")
      .def("cp_cellid",&mscomplex_t::cellid,
           "cell id of critical cell")

      .def("cps_index",&mscomplex_pyms3d_t::cps_index,
           "get Morse Indices of all critical points\n")
      .def("cps_func",&mscomplex_pyms3d_t::cps_func,
           "get Function values of all critical points\n")
      .def("cps_pairid",&mscomplex_pyms3d_t::cps_pairid,
           "get the cancellation pair ids of all critical points\n")
      .def("cps_vertid",&mscomplex_pyms3d_t::cps_vertid,
           "get maximal vertex id sof all critical points\n")
      .def("cps_cellid",&mscomplex_pyms3d_t::cps_cellid,
           "get cellids of all critical points\n")


      .def("asc",&mscomplex_pyms3d_t::conn<ASC>,
           "List of ascending cps connected to a given critical point i")
      .def("des",&mscomplex_pyms3d_t::conn<DES>,
           "List of descending cps connected to a given critical point i")
      .def("asc_geom",&mscomplex_pyms3d_t::geom<ASC>,
           (bp::arg("cp"),bp::arg("hversion")=-1,bp::arg("ToPts")=true),
           "Ascending manifold geometry of a given critical point i"
           "\n"
           "Parameters  : \n"
           "          cp: the critical point id\n"
           "    hversion: desired hierarchical version.\n"
           "              -1 indicates current version (default)"
           "       ToPts: convert the geometry data to point indices.\n"
           "               default = True     \n"
           "\n"
           "               False  : returns a list of cellids \n"
           "               True   : returns a list of Point idxs in \n"
           "                        Primal/Dual grid according to following\n"
           "                        table\n"
           "\n"
           "                         index(cp) pt-Index-Type  ArrayShape\n"
           "                             0      Primal          [NC]\n"
           "                             1      Dual            [NC',4]\n"
           "                             2      Dual            [NC',2]\n"
           "\n"
           "                       NC  #cells in Asc/Des mfold of cp.\n"
           "                       NC' #cells in Asc/Des mfold of cp whose \n"
           "                             dual pts are inside Primal Grid\n"
           )
      .def("des_geom",&mscomplex_pyms3d_t::geom<DES>,
           (bp::arg("cp"),bp::arg("hversion")=-1,bp::arg("ToPts")=true),
           "Descending manifold geometry of a given critical point i"
           "\n"
           "Parameters  : \n"
           "          cp: the critical point id\n"
           "    hversion: desired hierarchical version.\n"
           "              -1 indicates current version (default)"
           "       ToPts: convert the geometry data to point indices.\n"
           "               default = True     \n"
           "\n"
           "               False  : returns a list of cellids \n"
           "               True   : returns a list of Point idxs in \n"
           "                        Primal/Dual grid according to following\n"
           "                        table\n"
           "\n"
           "                         index(cp) pt-Index-Type  ArrayShape\n"
           "                             1      Primal          [NC,2]\n"
           "                             2      Primal          [NC,4]\n"
           "                             3      Dual            [NC]\n"
           "\n"
           "                       NC  #cells in Asc/Des mfold of cp.\n"
           "                       NC' #cells in Asc/Des mfold of cp whose \n"
           "                             dual pts are inside Primal Grid\n"

           )
      .def("primal_points",&mscomplex_pyms3d_t::points<CC_PRIM>,
           "Get coordinates of primal points of the grid")
      .def("dual_points",&mscomplex_pyms3d_t::points<CC_DUAL>,
           "Get coordinates of dual points of the grid")

      .def("compute_bin",&mscomplex_compute_bin,
           "Compute the Mscomplex from a structured grid with scalars given \n"\
           "in a raw/bin format \n"\
           "\n"\
           "Parameters  : \n"\
           "    bin_file: the bin/raw file containing the scalar function\n"\
           "              in float32 format \n"\
           "    size    : size of each dimension in x,y,z ordering .\n"\
//         "    bin_fmt : binary format .\n"\
//         "              Acceptable values = (\"float32\")\n"
           "\n"\
           "Note: This only computes the combinatorial structure\n"\
           "     Call collect_mfold(s) to extract geometry\n"
           )
      .def("collect_geom",&mscomplex_pyms3d_t::collect_mfolds,
           (bp::arg("dir")=2,bp::arg("dim")=-1),
           "Collect the geometry of all survivng critical points\n"\
           "\n"\
           "Parameters  : \n"\
           "         dir: Geometry type \n"\
           "              dir=0 --> Descending \n"\
           "              dir=1 --> Ascending \n"\
           "              dir=2 --> Both (default) \n"\
           "         dim: Critical point type \n"\
           "              dim=-1      --> All (default)\n"\
           "              dim=0,1,2,3 --> Minima, 1-saddle,2-saddle,Maxima \n"\
           "\n"\
           "Note        : Call only after the compute function is called.\n"\
           )
      .def("simplify_pers",&mscomplex_t::simplify_pers,
           (bp::arg("thresh")=1.0,bp::arg("is_nrm")=true,
            bp::arg("nmax")=0,bp::arg("nmin")=0),
           "Simplify the Morse-Smale complex using topological persistence\n"\
           "\n"
           "Parameters   :\n"\
           "    tresh    : persistence threshold\n"\
           "    is_nrm   : is the threshold normalized to [0,1] or not.\n"\
           "               if not then thresh is in scale of input function\n"\
           "    nmax,nmin: num maxima/minima that should be retained\n"\
           "                       set to 0 to ignore\n"\
           "\n"\
           "Note         : \n"\
           "    Any combination of the above criteria may be set\n"\
           "    Simplification will stop when any of the criteria is reached\n"
           "    Call only after the compute function is called.\n"
           )

      .def("load",&mscomplex_pyms3d_t::load,
           "Load mscomplex from file")
      .def("save",&mscomplex_pyms3d_t::save,
           "Save mscomplex to file")


      .def("get_hversion",&mscomplex_t::get_hversion,
           "Get the current Hierarchical Ms Complex version number")
      .def("set_hversion",&mscomplex_t::set_hversion,
           "Set the current Hierarchical Ms Complex version number")
      .def("get_hversion_pers",&mscomplex_t::get_hversion_pers,
           (bp::arg("thresh"),bp::arg("is_nrm")=true),
           "Get the highest hierarchical version num wherein all pairs \n"
           "with persistence less than the given value are canceled"
           "\n"
           "Parameters   :\n"
           "    tresh    : persistence threshold\n"
           "    is_nrm   : is the threshold normalized to [0,1] or not.\n"
           "               if not then thresh is in scale of input function\n"
           "\n"
           "Note         : \n"
           "    this presumes that the hierarchy generated was monotonic  \n"
           "    in the persistence of each canceled pair.                 \n"
           "\n"
           "    Here, persistence is used interchangably to mean the      \n"
           "    the absolute difference in function values of pairs.      \n"
           )
      .def("get_hversion_nextrema",&mscomplex_t::get_hversion_nextrema,
           (bp::arg("nmax")=0,bp::arg("nmin")=0),
           "Get the highest hierarchical no version which retains atleast \n"
           "nmax/nmin maxima/minima"
           "\n"
           "Parameters   :\n"
           "    nmax,nmin: num maxima/minima that should be retained\n"\
           "\n"
           )
      .def("canc_pairs",&mscomplex_pyms3d_t::canc_pairs,
           "get all cancellation pairs \n")

      ;
}

struct cellid_to_tup
{
  static PyObject* convert(cellid_t v)
  {return boost::python::incref(bp::make_tuple(v[0],v[1],v[2]).ptr());}
};


/*****************************************************************************/
/******** Define the pyms3d module                                    ********/
/*****************************************************************************/

BOOST_PYTHON_MODULE(pyms3d)
{
  boost::python::to_python_converter<cellid_t,cellid_to_tup>();

  np::initialize();

  grid::opencl::init();

  def("get_hw_info",&get_hw_info);

  def("ping",&pyms3d::ping);

  wrap_mscomplex_t();
}

int main(int argc, char **argv)
{
    // This line makes our module available to the embedded Python intepreter.
# if PY_VERSION_HEX >= 0x03000000
    PyImport_AppendInittab("pyms3d", &PyInit_pyms3d);
# else
    PyImport_AppendInittab("pyms3d", &initpyms3d);
# endif
    // Initialize the Python runtime.
    Py_Initialize();

    ENSURES(argc == 2) << "Usage : " << argv[0] << " <filename.py>" << endl;

    std::ifstream t(argv[1]);

    ENSURES(t.is_open()) << "Unable to open File="<<argv[1] << endl;

    std::string s((std::istreambuf_iterator<char>(t)),
                  std::istreambuf_iterator<char>());

    PyRun_SimpleString( s.c_str()   );
    Py_Finalize();

    return 0;
}


}/****************************************************************************/

