#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include <vector>
#include <timer.h>

#include <boost/shared_ptr.hpp>

#include <aabb.h>

namespace grid
{
  const uint gc_grid_dim = 3;
  const int g_num_threads = 8;

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
  typedef std::vector<char>                               char_list_t;
  typedef std::vector<cell_fn_t>                          cell_fn_list_t;
  typedef std::vector<bool>                               bool_list_t;
  typedef std::vector<rect_t>                             rect_list_t;

  enum eGDIR
  {
    GDIR_DES,
    GDIR_ASC,
    GDIR_CT
  };

  class dataset_t;
  class mscomplex_t;
  class data_manager_t;

  typedef boost::shared_ptr<dataset_t>            dataset_ptr_t;
  typedef boost::shared_ptr<mscomplex_t>          mscomplex_ptr_t;
  typedef boost::shared_ptr<data_manager_t>       data_manager_ptr_t;

  typedef boost::shared_ptr<const dataset_t>      dataset_const_ptr_t;
  typedef boost::shared_ptr<const mscomplex_t>    mscomplex_const_ptr_t;
  typedef boost::shared_ptr<const data_manager_t> data_manager_const_ptr_t;

  class cp_producer_t;
  typedef boost::shared_ptr<cp_producer_t>  cp_producer_ptr_t;

  extern "C"
  Timer g_timer;

  template <typename T>
  inline std::ostream& log_range(T b,T e,std::ostream &os = std::cout,char sep = ' ')
  {
    for (;b!=e;++b) os<<*b<<sep;
    return os;
  }
}

namespace std
{
template<typename _IIter>
  typename iterator_traits<_IIter>::difference_type
  count(_IIter b, _IIter e)
  {
    typename iterator_traits<_IIter>::difference_type val = 0;

    for( ; b != e; ++b)
      ++val;

    return val;
  }
}

#define _FFL            (std::string("\n")+FILEFUNCLINE)

#endif
