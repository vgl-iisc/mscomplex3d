#ifndef __GRID_MSCOMPLEX_INL_H_INCLUDED_
#define __GRID_MSCOMPLEX_INL_H_INCLUDED_

#include <fstream>

#include <boost/range/algorithm.hpp>

#include <grid_mscomplex.h>


namespace grid {

/*===========================================================================*/

inline int  mscomplex_t::get_num_critpts() const
{return m_cp_cellid.size();}

/*---------------------------------------------------------------------------*/

inline int8_t mscomplex_t::index(int i) const
{
  ASSERTV(is_in_range(i,0,(int)m_cp_index.size()),i);
  return m_cp_index[i];
}

/*---------------------------------------------------------------------------*/

inline int mscomplex_t::pair_idx(int i) const
{
  ASSERTV(is_in_range(i,0,(int)m_cp_pair_idx.size()),i);
  ASSERTV(is_in_range(m_cp_pair_idx[i],0,(int)m_cp_pair_idx.size()),i);
  return m_cp_pair_idx[i];
}

/*---------------------------------------------------------------------------*/

inline bool mscomplex_t::is_paired(int i) const
{
  ASSERTV(is_in_range(i,0,(int)m_cp_pair_idx.size()),i);
  return (m_cp_pair_idx[i] != -1);
}

/*---------------------------------------------------------------------------*/

inline bool mscomplex_t::is_not_paired(int i) const
{
  ASSERTV(is_in_range(i,0,(int)m_cp_pair_idx.size()),i);

  return (m_cp_pair_idx[i] == -1);
}

/*---------------------------------------------------------------------------*/

inline bool mscomplex_t::is_canceled(int i) const
{
  ASSERTV(is_in_range(i,0,(int)m_cp_is_cancelled.size()),i);
  return m_cp_is_cancelled[i];
}

/*---------------------------------------------------------------------------*/

inline bool mscomplex_t::is_not_canceled(int i) const
{
  return !is_canceled(i);
}

/*---------------------------------------------------------------------------*/

inline cellid_t mscomplex_t::cellid(int i) const
{
  ASSERTV(is_in_range(i,0,(int)m_cp_cellid.size()),i);
  return m_cp_cellid[i];
}

/*---------------------------------------------------------------------------*/

inline cellid_t mscomplex_t::vertid(int i) const
{
  ASSERTV(is_in_range(i,0,(int)m_cp_vertid.size()),i);
  return m_cp_vertid[i];
}

/*---------------------------------------------------------------------------*/

inline cell_fn_t mscomplex_t::fn(int i) const
{
  ASSERTV(is_in_range(i,0,(int)m_cp_fn.size()),i);
  return m_cp_fn[i];
}

/*---------------------------------------------------------------------------*/

inline bool mscomplex_t::is_extrema(int i) const
{
  return (index(i) == 0 || index(i) == 3);
}

/*---------------------------------------------------------------------------*/

inline bool mscomplex_t::is_saddle(int i) const
{
  return (index(i) == 1 || index(i) == 2);
}

/*---------------------------------------------------------------------------*/

template<int i> inline bool mscomplex_t::is_index_i_cp(int cp) const
{return (index(cp) == i);}

/*---------------------------------------------------------------------------*/

inline bool mscomplex_t::is_index_i_cp_(int cp, int i) const
{return (index(cp) == i);}


/*---------------------------------------------------------------------------*/

inline int mscomplex_t::surv_extrema(int i) const
{
  ASSERT(is_extrema(i));

  if(is_canceled(i) == false)
    return i;

  eGDIR dir = (index(i) == 3)?(ASC):(DES);

  ASSERT(m_conn[dir][pair_idx(i)].size() == 1);

  int j = m_conn[dir][pair_idx(i)].begin()->first;

  ASSERT(is_not_canceled(j));

  return j;
}

/*---------------------------------------------------------------------------*/

inline cell_fn_t mscomplex_t::fn_min() const
{return *boost::range::min_element(m_cp_fn);}

/*---------------------------------------------------------------------------*/

inline cell_fn_t mscomplex_t::fn_max() const
{return *boost::range::max_element(m_cp_fn);}

/*---------------------------------------------------------------------------*/

inline std::string mscomplex_t::cp_info (int cp_no) const
{
  std::stringstream ss;

  ss<<std::endl;
  ss<<"cp_no        ::"<<cp_no<<std::endl;
  ss<<"cellid       ::"<<cellid(cp_no)<<std::endl;
//    ss<<"vert cell    ::"<<vertid(cp_no)<<std::endl;
  ss<<"index        ::"<<(int)index(cp_no)<<std::endl;
//      ss<<"fn           ::"<<fn(cp_no)<<std::endl;
//    ss<<"is_cancelled ::"<<is_canceled(cp_no)<<std::endl;
//    ss<<"is_paired    ::"<<is_paired(cp_no)<<std::endl;
  ss<<"pair_idx     ::"<<pair_idx(cp_no)<<std::endl;
  return ss.str();
}

/*---------------------------------------------------------------------------*/

inline boost::iterator_range<mscomplex_t::iterator_t> mscomplex_t::cpno_range() const
{return boost::make_iterator_range(iterator_t(0),iterator_t(get_num_critpts()));}

/*---------------------------------------------------------------------------*/

inline void mscomplex_t::save(const std::string &f)
{
  std::fstream fs(f.c_str(),std::ios::out|std::ios::binary);
  ENSUREV(fs.is_open(),"file not found!!",f);
  save_bin(fs);
}

/*---------------------------------------------------------------------------*/

inline void mscomplex_t::load(const std::string &f)
{
  std::fstream fs(f.c_str(),std::ios::in|std::ios::binary);
  ENSUREV(fs.is_open(),"file not found!!",f);
  load_bin(fs);
}

/*---------------------------------------------------------------------------*/

inline int mscomplex_t::get_hversion() const {return m_hversion;}

/*===========================================================================*/



/*===========================================================================*/

inline void order_pr_by_cp_index(const mscomplex_t &msc,int &p,int &q)
{if(msc.index(p) < msc.index(q))std::swap(p,q);}

/*---------------------------------------------------------------------------*/

template<> inline int_pair_t order_pair<DES>(mscomplex_ptr_t msc,int_pair_t pr)
{if(msc->index(pr[0]) < msc->index(pr[1]))
    std::swap(pr[0],pr[1]);return pr;}

template<> inline int_pair_t order_pair<ASC>(mscomplex_ptr_t msc,int_pair_t pr)
{if(msc->index(pr[0]) > msc->index(pr[1]))
    std::swap(pr[0],pr[1]); return pr;}

/*===========================================================================*/




/*===========================================================================*/

/// \brief Data structure to extract MSC geometry in hierarchical versions
class mscomplex_t::merge_dag_t
{
public:
  /// \brief ctor
  merge_dag_t();

  /// \brief init
  void init(int ncps);

  /// \brief clear
  void clear();



  /// \brief For the given critical point (cp) , get the cps that
  /// that contribute to its geometry at the given hierarchical version.
  ///
  /// \param  hver : the hierarchical version of the geom requested
  /// \param ghver : the hierarchical version at which the fine geom is available
  void get_contrib_cps(std::vector<int> &l, eGDIR dir, int cp, int hver, int gver) const;

  /// \brief update merge_dag after new cancellations
  void update(mscomplex_ptr_t msc);



  /// \brief save into a binary stream
  void save_bin(std::ostream &os) const;

  /// \brief load from a binary stream
  void load_bin(std::istream &is);


//private:

  struct node_t
  {
    int base;
    int other;
    int hversion;
    int cpid;

    node_t():cpid(-1),base(-1),other(-1),hversion(0){}
    node_t(int c):cpid(c),base(-1),other(-1),hversion(0){}
    node_t(int c,int b,int o,int h):cpid(c),base(b),other(o),hversion(h){}
  };

  std::vector<node_t>  m_nodes;
  std::vector<int>     m_cp_geom[2];
  int                  m_last_hversion;

  inline node_t get_node(int i) const;
  inline int  get_ncps() const;
};

/*===========================================================================*/

}

#endif
