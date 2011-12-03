#ifndef GRID_H_INCLUDED
#define GRID_H_INCLUDED

#include <vector>
#include <timer.h>

#include <boost/shared_ptr.hpp>
#include <boost/multi_array.hpp>

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
  typedef std::vector<char>                               bool_list_t;
  typedef std::vector<rect_t>                             rect_list_t;

  typedef boost::shared_ptr<int_list_t>                   int_list_ptr_t;
  typedef boost::shared_ptr<cellid_list_t>                cellid_list_ptr_t;

  typedef boost::multi_array<int,gc_grid_dim>             int_marray_t;

  enum eGDIR  {GDIR_DES,GDIR_ASC,GDIR_CT};

  class dataset_t;
  class mscomplex_t;
  class data_manager_t;

  typedef boost::shared_ptr<dataset_t>            dataset_ptr_t;
  typedef boost::shared_ptr<mscomplex_t>          mscomplex_ptr_t;
  typedef boost::shared_ptr<data_manager_t>       data_manager_ptr_t;

  typedef boost::shared_ptr<const dataset_t>      dataset_const_ptr_t;
  typedef boost::shared_ptr<const mscomplex_t>    mscomplex_const_ptr_t;
  typedef boost::shared_ptr<const data_manager_t> data_manager_const_ptr_t;

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
  Timer g_timer;

  template <typename T>
  inline std::ostream& log_range(T b,T e,std::ostream &os = std::cout,char sep = ' ')
  {
    for (;b!=e;++b) os<<*b<<sep;
    return os;
  }

  inline int get_cell_dim ( cellid_t c )
  {return ( c[0]&0x01 ) + ( c[1]&0x01 ) + ( c[2]&0x01 );}

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
