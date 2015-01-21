#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>

#include <aabb.h>

namespace grid
{
  const uint gc_grid_dim = 3;

  typedef int16_t                                         cell_coord_t;
  typedef u_int8_t                                        cell_flag_t;
  typedef float                                           cell_fn_t;
  typedef boost::shared_ptr<cell_fn_t >                   cell_fn_ptr_t;
  typedef aabb::aabb_t<cell_coord_t,gc_grid_dim>          rect_t;
  typedef aabb::aabb_t<cell_coord_t,gc_grid_dim>::point_t cellid_t;
  typedef aabb::aabb_t<cell_coord_t,gc_grid_dim>::point_t rect_point_t;
  typedef aabb::aabb_t<cell_coord_t,gc_grid_dim>::point_t rect_size_t;
  typedef aabb::aabb_t<cell_coord_t,gc_grid_dim>::range_t rect_range_t;
  typedef std::vector<cellid_t>                           cellid_list_t;
  typedef std::vector<int>                                int_list_t;
  typedef std::vector<int8_t>                             int8_list_t;
  typedef std::vector<cell_fn_t>                          cell_fn_list_t;
  typedef std::vector<int8_t>                             bool_list_t;
  typedef std::vector<rect_t>                             rect_list_t;

  typedef boost::shared_ptr<int_list_t>                   int_list_ptr_t;
  typedef boost::shared_ptr<cellid_list_t>                cellid_list_ptr_t;

  typedef boost::multi_array<int,gc_grid_dim>             int_marray_t;
  typedef std::pair<cellid_t,int>                         cellid_int_pair_t;
  typedef std::vector<cellid_int_pair_t>                  cellid_int_pair_list_t;

  typedef cellid_list_t                                   mfold_t;
  typedef std::vector<mfold_t>                            mfold_list_t;

  typedef n_vector_t<int,2>                               int_pair_t;
  typedef std::vector<int_pair_t>                         int_pair_list_t;

  typedef std::pair<int,int>                              int_int_t;
  typedef std::vector<int_int_t>                          int_int_list_t;
  typedef std::map<int,int>                               int_int_map_t;

  typedef std::map<int,int>                               conn_t;
  typedef std::vector<conn_t>                             conn_list_t;

  enum eGDIR  {DES=0,ASC,GDIR_CT};

  const eGDIR GDIR_DES=DES;
  const eGDIR GDIR_ASC=ASC;

  /// \brief cell complex type
  enum eCCTYPE
  {
    CC_NONE=0,
    CC_PRIM=1, // primal
    CC_DUAL=2, // Dual
    CC_BOTH=3
  };

  class dataset_t;
  class mscomplex_t;
  class data_manager_t;

  typedef boost::shared_ptr<dataset_t>            dataset_ptr_t;
  typedef boost::shared_ptr<mscomplex_t>          mscomplex_ptr_t;
  typedef boost::shared_ptr<data_manager_t>       data_manager_ptr_t;

  inline int c_to_i(const rect_t &r,cellid_t c)
  {
    cellid_t s = r.span()+1;
    c = (c - r.lc());
    return (s[0]*s[1]*c[2] + s[0]*c[1] + c[0]);
  }

  inline cellid_t i_to_c(const rect_t &r,int i)
  {
    cellid_t s = r.span()+1;
    cellid_t c = r.lc() + (cellid_t(i%s[0],(i%(s[0]*s[1]))/s[0],i/(s[0]*s[1])));
    ASSERT(r.contains(c));
    return c;
  }

  inline int num_cells(const rect_t &r)
  {return c_to_i(r,r.uc()) + 1;}

  extern "C"
  utl::timer g_timer;

  inline int get_cell_dim ( cellid_t c )
  {return ( c[0]&0x01 ) + ( c[1]&0x01 ) + ( c[2]&0x01 );}

}

namespace grid{
namespace opencl{

/// \brief Init OpenCL runtime
void init();

/// \brief Return the platform/device info after init
std::string get_info();
}

namespace openmp {
/// \brief Return the openmp info
std::string get_info();
}

inline std::string get_hw_info(){return opencl::get_info()+openmp::get_info();}
}

#endif
