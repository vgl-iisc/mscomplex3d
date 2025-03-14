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


#include <iostream>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>


#include <grid_dataset.h>
#include <grid_mscomplex.h>

#include <chrono>

struct timer
{
    int flag = 0;
    std::chrono::duration<double> time;
    std::chrono::time_point<std::chrono::steady_clock> start;

    void timer_start()
    {
        if (flag == 1)
        {
            std::cout << "Timer already started, cannot start timer again. End it to start.\n";
        }
        else
        {
            start = std::chrono::steady_clock::now();  // Use steady_clock for better cross-platform consistency
            flag = 1;
        }
    }

    void timer_end()
    {
        if (flag == 0)
        {
            std::cout << "Timer was not started. Please start the timer before ending it.\n";
            return;
        }
        auto end = std::chrono::steady_clock::now();
        time = end - start;
        flag = 0; // Reset flag so timer can be restarted
    }

    void timer_print(const std::string& process_name)
    {
        if (flag == 1)
        {
            std::cout << "Warning: Timer is still running. Stopping it now.\n";
            timer_end();  // Ensure timer stops before printing
        }
        std::cout << process_name << " took: " << time.count() << " seconds\n";
        time = std::chrono::duration<double>(0);
    }
}time_keeper;







#ifdef _WIN32
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif

using namespace std;
using namespace grid;
using namespace utl;

namespace py = pybind11;

namespace utl {
template<class rng_t>
inline std::string to_string_rng(rng_t rng,const char * delim = ", ")
{
  std::stringstream ss;

  auto b = std::begin(rng);
  auto e = std::end(rng);

  ss<<"["; for(; b!=e ; b) ss << *b++ << delim; ss<<"]";
  return ss.str();
}
}

/**
 * \brief Converts a std::vector<T> to a python array
 * \param vec the std::vector object
 * \return 
 */
template<typename DTYPE, int NCOL, typename T>
py::array_t<DTYPE> vector_to_ndarray(const std::vector<T>& vec) {
    // Ensure the sizes match
    static_assert(sizeof(DTYPE) * NCOL == sizeof(T), "Size mismatch between DTYPE and T");

    // Get the size of the vector
    int N = vec.size();

    // Define the shape of the NumPy array
    py::array_t<DTYPE> arr;
    if (NCOL == 1) {
        arr = py::array_t<DTYPE>({ N });
    }
    else {
        arr = py::array_t<DTYPE>({ N, NCOL });
    }

    // Access raw data buffer for writing
    auto buf = arr.request();
    DTYPE* ptr = static_cast<DTYPE*>(buf.ptr);

    // Copy data from the vector into the NumPy array
    std::memcpy(ptr, vec.data(), vec.size() * sizeof(T));

    return arr;
}

template <typename DTYPE, int NCOL, typename rng_t>
py::array_t<DTYPE> range_to_ndarray(rng_t& rng) {
    // Ensure that rng_t is a valid range
    static_assert(std::ranges::range<rng_t>, "rng_t must be a valid range.");

    // Deduce the value type from the range
    using T = std::iter_value_t<std::ranges::iterator_t<rng_t>>;

    // Ensure size compatibility between the data types
    //static_assert(sizeof(DTYPE) * NCOL == sizeof(T), "Size mismatch between DTYPE * NCOL and T.");
    //static_assert(NCOL > 0, "NCOL must be a positive integer.");

    // Get the number of elements in the range
    ssize_t N = std::ranges::distance(rng);

    // Determine the shape of the array
    std::vector<ssize_t> shape = (NCOL == 1) ? std::vector<ssize_t>{N} : std::vector<ssize_t>{ N, NCOL };

    // Create a zero-initialized NumPy array
    py::array_t<DTYPE> arr(shape);

    //std::cout << "\nARRAY FOUND: \n";
   // py::print(arr);

    // Get a pointer to the data in the array
    auto buf = arr.request(); // Get buffer information
    char* ptr = static_cast<char*>(buf.ptr);

    // Copy data from the range into the array
    for (const auto& v : rng) { // Range-based for loop works for all ranges
        static_assert(std::is_trivially_copyable_v<T>, "T must be trivially copyable.");
        std::memcpy(ptr, &v, sizeof(T));
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

    //ds->load_bin(fs);
    if (!ds) {
        throw std::runtime_error("ds is null");
    }
    else
    ds->m_rect = m_rect;

    DLOG << "Exited  :";
  }

  /*-------------------------------------------------------------------------*/

  template <eGDIR dir>
  py::array_t<int> conn(int cp) {

  	ENSURES(is_in_range(cp, 0, get_num_critpts())) << "out of range " << SVAR(cp);

  	TLOG << SVAR(cp);

  	return range_to_ndarray<int, 2>(m_conn[dir][cp]);
  }

  /*-------------------------------------------------------------------------*/
  
  template <eCCTYPE ccTYPE, int dim>
  py::array_t<int> __mfold_to_point_indices(const mfold_t& mfold) {
      const int nCellPoints = (ccTYPE == CC_PRIM) ? (1 << dim) : (1 << (gc_grid_dim - dim));
      const int O = (ccTYPE == CC_PRIM) ? (0) : (1);

      rect_t prect = (ccTYPE == CC_PRIM) ? (m_rect) : (rect_t(m_rect.lc() + 1, m_rect.uc() - 1));

      std::vector<int> points;

      for (int i = 0; i < mfold.size(); ++i) {
          cellid_t c = mfold[i];
          std::array<int, 3> j;

          bool need_cell = true;

          // Loop over neighbors (3D case, can be extended for different dimensions)
          for (j[2] = -((c[2] + O) & 1); j[2] <= ((c[2] + O) & 1); j[2] += 2) {
              for (j[1] = -((c[1] + O) & 1); j[1] <= ((c[1] + O) & 1); j[1] += 2) {
                  for (j[0] = -((c[0] + O) & 1); j[0] <= ((c[0] + O) & 1); j[0] += 2) {
                      if (!prect.contains(c + j)) {
                          need_cell = false;
                      }
                  }
              }
          }

          if (need_cell) {
              for (j[2] = -((c[2] + O) & 1); j[2] <= ((c[2] + O) & 1); j[2] += 2) {
                  for (j[1] = -((c[1] + O) & 1); j[1] <= ((c[1] + O) & 1); j[1] += 2) {
                      for (j[0] = -((c[0] + O) & 1); j[0] <= ((c[0] + O) & 1); j[0] += 2) {
                          points.push_back(c_to_i2(prect, c + j));
                      }
                  }
              }
          }
      }

      // Create a NumPy array from the vector of points
      py::array_t<int> arr(points.size(), points.data());

      if (nCellPoints != 1) {

          // Create a ShapeContainer using std::vector<std::size_t>
          py::array::ShapeContainer shape = { static_cast<std::size_t>(points.size() / nCellPoints), static_cast<std::size_t>(nCellPoints) };

          // Resize the array with the new shape
          arr.resize(shape);
          
      }

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

  template <eGDIR dir> py::array geom(int cp, int hver = -1, bool ToPts = true)
  {
      TLOG << "Entered :" << SVAR(cp) << SVAR(hver);

      if (hver == -1) hver = get_hversion();

      ENSURES(is_in_range(cp, 0, get_num_critpts()))
          << "out of range " << SVAR(cp);
      ENSURES(is_in_range(hver, 0, m_canc_list.size() + 1))
          << "hversion not in range " << SVAR(hver);
      //BOOST_STATIC_ASSERT(dir == ASC || dir == DES);

      int_list_t l;

      int dim = index(cp);

      m_merge_dag->get_contrib_cps
      (l, dir, cp, hver, m_geom_hversion[dir][dim]);

      mfold_t mfold;

      for (int j = 0; j < l.size(); ++j) {
          const auto& sublist = m_mfolds[dir][l[j]];

      	  std::ranges::copy(sublist, std::back_inserter(mfold));
      }

      TLOG << "Computed:" << SVAR(mfold.size());

      if (ToPts)
      {
          if (dir == ASC && dim == 0) return __mfold_to_point_indices<CC_PRIM, 0>(mfold);
          if (dir == ASC && dim == 1) return __mfold_to_point_indices<CC_DUAL, 1>(mfold);
          if (dir == ASC && dim == 2) return __mfold_to_point_indices<CC_DUAL, 2>(mfold);

          if (dir == DES && dim == 1) return __mfold_to_point_indices<CC_PRIM, 1>(mfold);
          if (dir == DES && dim == 2) return __mfold_to_point_indices<CC_PRIM, 2>(mfold);
          if (dir == DES && dim == 3) return __mfold_to_point_indices<CC_DUAL, 3>(mfold);

          ENSURES(false) << "Should never reach here";
      }

      return vector_to_ndarray<cell_coord_t, gc_grid_dim>(mfold);
  }
 


  /*-------------------------------------------------------------------------*/

  template <eCCTYPE pTYPE>
  py::array_t<float> points() {
      // Assuming m_rect and other data structures are defined appropriately
      rect_t prect = m_rect;

      // Adjust prect if the pTYPE is CC_DUAL
      if (pTYPE == CC_DUAL) {
          prect = rect_t(m_rect.lc() + 1, m_rect.uc() - 1);
      }

      // Calculate N based on prect
      int N = 1 + c_to_i2(prect, prect.uc());

      // Create a numpy array with the dtype of float32
		// Create a NumPy array with the given shape
      py::array_t<float> arr({ N, 3 });

      // Access the raw pointer to the data and initialize it with zeros
      auto r = arr.mutable_unchecked<2>(); // 2D array accessor for modification
      for (ssize_t i = 0; i < r.shape(0); ++i) {
          for (ssize_t j = 0; j < r.shape(1); ++j) {
              r(i, j) = 0.0f; // Set all elements to zero
          }
      }
      // Pointer to the data of the numpy array
      float* iter = arr.mutable_data();

      // Static assert for 3-manifolds
      static_assert(gc_grid_dim == 3, "Defined for 3-manifolds only");

      cellid_t c;

      // Loop through the grid and populate the array with coordinates
      for (c[2] = prect.lc()[2]; c[2] <= prect.uc()[2]; c[2] += 2) {
          for (c[1] = prect.lc()[1]; c[1] <= prect.uc()[1]; c[1] += 2) {
              for (c[0] = prect.lc()[0]; c[0] <= prect.uc()[0]; c[0] += 2) {
                  // Store the coordinates in the array (scaled by 1/2)
                  *iter++ = static_cast<float>(c[0]) / 2;
                  *iter++ = static_cast<float>(c[1]) / 2;
                  *iter++ = static_cast<float>(c[2]) / 2;
              }
          }
      }

      return arr;
  }

  /*-------------------------------------------------------------------------*/

  /**
   * \brief Returns a python array of all critical points of index i
   * \param i Critical point index
   * \return 
   */
  py::array cps(int i)
  {
      TLOG << "Entered :";

      auto rng = cpno_range()
          | std::ranges::views::filter([this](const auto& item) {
          return this->is_not_canceled(item);
              });

      TLOG << "Computed:";

      if (i == -1) return range_to_ndarray<int, 1>(rng);
      else {

          auto filtered_rng = std::ranges::views::filter(rng, [this, i](int x) {
              return this->is_index_i_cp_(x, i);
              });
          std::vector<int> filtered_vector(filtered_rng.begin(), filtered_rng.end());
          
          return range_to_ndarray<int, 1>(filtered_rng);
      }
  	TLOG << "Exited  :";

  }

  /*-------------------------------------------------------------------------*/

  py::array_t<int> canc_pairs() { TLOG; return vector_to_ndarray<int, 2>(m_canc_list); }
  py::array_t<int> cps_func() { TLOG; return vector_to_ndarray<cell_fn_t, 1>(m_cp_fn); }
  py::array_t<int> cps_index() { TLOG; return vector_to_ndarray<int8_t, 1>(m_cp_index); }
  py::array_t<int> cps_pairid() { TLOG; return vector_to_ndarray<int, 1>(m_cp_pair_idx); }
  py::array_t<int> cps_cellid() { TLOG; return vector_to_ndarray<cell_coord_t, 3>(m_cp_cellid); }
  py::array_t<int> cps_vertid() { TLOG; return vector_to_ndarray<cell_coord_t, 3>(m_cp_vertid); }

  /*-------------------------------------------------------------------------*/

  void collect_mfolds(int dir ,int dim)
  {
    DLOG << "Entered :" << SVAR(dir) << SVAR(dim);
    ENSURES(ds !=0)<< "Gradient information unavailable";
    base_t::collect_mfolds((eGDIR)dir,dim,ds);
    DLOG << "Exited  :";
  }
  /*-------------------------------------------------------------------------*/

  void compute_arr(py::array arr)
  {
      ENSURES(arr.ndim() == 3) << " 3D array expected " << std::endl;

      ENSURES(arr.flags())// & np::ndarray::ALIGNED) we dont need to check for aligned with pybind
          << "Contiguous array expected " << std::endl;

      //arr = arr.astype(np::dtype::get_builtin<float>());

  	  py::array float_arr = arr.cast<py::array_t<float>>();

  	  //  ENSURES(bin_fmt == "float32") <<"Only float32 is supported" <<endl;
      int x = arr.shape(0);
      int y = arr.shape(1);
      int z = arr.shape(2);
      cellid_t size = cellid_t(x, y, z);

      rect_t dom(cellid_t::zero, (size - cellid_t::one) * 2);

      DLOG << "Entered :"
          << endl << "\t" << SVAR(size);

      m_rect = dom;
      m_domain_rect = dom;
      m_ext_rect = dom;

      bool is_fortran = (arr.flags() & py::array::f_style);

      ds.reset(new dataset_t(dom, dom, dom));
      const cell_fn_t* data = reinterpret_cast<const cell_fn_t*>(arr.data());
      ds->init(data, is_fortran);
      ds->computeMsGraph(shared_from_this());

      DLOG << "Exited  :";
  }

  /*-------------------------------------------------------------------------*/

  cell_fn_t vert_fn(int x, int y, int z)
  {
    return ds->m_vert_fns(cellid_t(x,y,z));
  }

};

void ping()
{
  ILOG <<"Pinged  :";
}
}

/*===========================================================================*/

namespace pyms3d {

typedef  std::shared_ptr<mscomplex_pyms3d_t> mscomplex_pyms3d_ptr_t;

mscomplex_pyms3d_ptr_t new_msc()
{
  mscomplex_pyms3d_ptr_t msc(new mscomplex_pyms3d_t);
  return msc;
}

	/**
	 * \brief Read 3D scalar grid values from a given .raw file
	 * \param msc the datastructure used to store the grid values and subsequent critical points
	 * \param bin_file the .raw file used to sotre scalar values
	 * \param tp dimensions of the 3D scalar field
	 */
	void mscomplex_compute_bin
(mscomplex_pyms3d_ptr_t &msc,
    std::string            bin_file,
    py::tuple              tp)
{
    //  ENSURES(bin_fmt == "float32") <<"Only float32 is supported" <<endl;
    int x = tp[0].cast<int>();
    int y = tp[1].cast<int>();
    int z = tp[2].cast<int>();
	
    cellid_t size = cellid_t(x, y, z);

    rect_t dom(cellid_t::zero, (size - cellid_t::one) * 2);

    DLOG << "Entered :"
	<< endl << "\t" << SVAR(bin_file)
	<< endl << "\t" << SVAR(size);

    msc->m_rect = dom;
    msc->m_domain_rect = dom;
    msc->m_ext_rect = dom;


    time_keeper.timer_start();
    msc->ds.reset(new dataset_t(dom, dom, dom));
    time_keeper.timer_end();
    time_keeper.timer_print("reset");

    time_keeper.timer_start();
	msc->ds->init(bin_file);
    time_keeper.timer_end();
    time_keeper.timer_print("file read");

    time_keeper.timer_start();
	msc->ds->computeMsGraph(msc);
    time_keeper.timer_end();
    time_keeper.timer_print("MS Complex Processing");

    DLOG << "Exited  :";
}


/**
 * \brief Performs binding of C++ functions to Python functions that are exposed via the MsComplexPyms3D object 
 * \param m 
 */
void init_mscomplex(py::module_& m) {

    
    py::class_<mscomplex_t, std::shared_ptr<mscomplex_t>>(m, "MsComplex")
        .def("num_cps", &mscomplex_t::get_num_critpts, "Number of Critical Points")
        .def("cp_func", &mscomplex_t::fn, "Function value at critical point i")
        .def("cp_index", &mscomplex_t::index, "Morse index of critical point i")
        .def("cp_pairid", &mscomplex_t::pair_idx, "Index of the paired critical point i (-1 if unpaired)")
        .def("cp_vertid", &mscomplex_t::vertid, "Vertex id of the critical cell")
        .def("cp_cellid", &mscomplex_t::cellid, "Cell id of the critical cell")
        //.def("is_boundry", &mscomplex_t::, "If the cp is on the boundary or not")
        .def("simplify_pers", &mscomplex_t::simplify_pers,
            py::arg("thresh") = 1.0,
            py::arg("is_nrm") = true,
            py::arg("nmax") = 0,
            py::arg("nmin") = 0,
            R"doc(
Simplify the Morse-Smale complex using topological persistence.

Parameters:
    thresh: Persistence threshold.
    is_nrm: Indicates if the threshold is normalized to [0,1].
            If false, the threshold is in the scale of the input function.
    nmax, nmin: Number of maxima/minima that should be retained.
                Set to 0 to ignore.

Note:
    Any combination of the above criteria may be set.
    Simplification will stop when any of the criteria is reached.
    Call only after the compute function is called.
)doc");
    //the reason we define this second python class is because this is how we take an inherited function from an existing one in Pybind
    py::class_<mscomplex_pyms3d_t, mscomplex_t, std::shared_ptr<mscomplex_pyms3d_t>>(m, "MsComplexPyms3D")
        .def(py::init<>())
	    .def("cps", &mscomplex_pyms3d_t::cps, py::arg("dim") = -1,
            "Returns a list of surviving critical points based on dimension")
        .def("cps_func", &mscomplex_pyms3d_t::cps_func, "Get function values of all critical points")
        .def("cps_index", &mscomplex_pyms3d_t::cps_index, "Get Morse indices of all critical points")
        .def("cps_pairid", &mscomplex_pyms3d_t::cps_pairid, "Get cancellation pair ids of all critical points")
        .def("cps_vertid", &mscomplex_pyms3d_t::cps_vertid, "Get maximal vertex ids of all critical points")
        .def("cps_cellid", &mscomplex_pyms3d_t::cps_cellid, "Get cell ids of all critical points")
        .def("asc", &mscomplex_pyms3d_t::conn<ASC>, "List of ascending cps connected to a critical point")
        .def("des", &mscomplex_pyms3d_t::conn<DES>, "List of descending cps connected to a critical point")
        .def("primal_points", &mscomplex_pyms3d_t::points<CC_PRIM>, "Get primal grid point coordinates")
        .def("dual_points", &mscomplex_pyms3d_t::points<CC_DUAL>, "Get dual grid point coordinates")
        .def("compute_bin", &mscomplex_compute_bin, "Compute the MsComplex from a raw/bin file")
        .def("compute_arr", &mscomplex_pyms3d_t::compute_arr, "Compute the MsComplex from a 3D numpy array")
        .def("load", &mscomplex_pyms3d_t::load, "Load mscomplex from file")
        .def("save", &mscomplex_pyms3d_t::save, "Save mscomplex to file")
        .def("vert_func", &mscomplex_pyms3d_t::vert_fn,
            "Scalar value at vertex coordinate")
        .def("canc_pairs", &mscomplex_pyms3d_t::canc_pairs,
            "get all cancellation pairs \n")
        .def("collect_geom", &mscomplex_pyms3d_t::collect_mfolds,
            py::arg("dir") = 2, 
            py::arg("dim") = -1,
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
        );
        /*
	.def("asc_geom", &mscomplex_pyms3d_t::geom<ASC>,
            py::arg("cp"), py::arg("hversion") = -1, py::arg("ToPts") = true,
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
        );
        */
    /*
        .def("des_geom", &mscomplex_pyms3d_t::geom<DES>,
            (py::arg("cp"), py::arg("hversion") = -1, py::arg("ToPts") = true),
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
        );
        */
}

struct cellid_to_tup {
    static py::tuple convert(const cellid_t& v) {
        return py::make_tuple(v[0], v[1], v[2]);
    }
};


void debug_print() {
    std::cout << "debug";
}

/// <summary>
/// This is where the PyBind takes place. Make sure the first argument has the same name as the module you want to create.
/// For now, the name of this module is pyms3d_core
/// </summary>
PYBIND11_MODULE(pyms3d_core, m) {
    // Call the function that registers the bindings
    init_mscomplex(m);
    m.doc() = "This is a custom module implemented in C++ using Pybind11.";  // Module docstring
        m.def("debug_print", &debug_print, "A function that prints a simple debug message");
        m.def("get_hw_info",&get_hw_info, py::arg("device") = 0);

}

}/****************************************************************************/

