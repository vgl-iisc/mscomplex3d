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

inline char mscomplex_t::index(int i) const
{
  ASSERTS(is_in_range(i,0,(int)m_cp_index.size()))<<STRMVAR(i);
  return m_cp_index[i];
}

/*---------------------------------------------------------------------------*/

inline int mscomplex_t::pair_idx(int i) const
{
  ASSERTS(is_in_range(i,0,(int)m_cp_pair_idx.size()))<<STRMVAR(i);
  ASSERTS(is_in_range(m_cp_pair_idx[i],0,(int)m_cp_pair_idx.size()))<<STRMVAR(i);
  return m_cp_pair_idx[i];
}

/*---------------------------------------------------------------------------*/

inline bool mscomplex_t::is_paired(int i) const
{
  ASSERTS(is_in_range(i,0,(int)m_cp_pair_idx.size()))<<STRMVAR(i);
  return (m_cp_pair_idx[i] != -1);
}

/*---------------------------------------------------------------------------*/

inline bool mscomplex_t::is_canceled(int i) const
{
  ASSERTS(is_in_range(i,0,(int)m_cp_is_cancelled.size()))<<STRMVAR(i);
  return m_cp_is_cancelled[i];
}

/*---------------------------------------------------------------------------*/

inline cellid_t mscomplex_t::cellid(int i) const
{
  ASSERTS(is_in_range(i,0,(int)m_cp_cellid.size()))<<STRMVAR(i);
  return m_cp_cellid[i];
}

/*---------------------------------------------------------------------------*/

inline cellid_t mscomplex_t::vertid(int i) const
{
  ASSERTS(is_in_range(i,0,(int)m_cp_vertid.size()))<<STRMVAR(i);
  return m_cp_vertid[i];
}

/*---------------------------------------------------------------------------*/

inline cell_fn_t mscomplex_t::fn(int i) const
{
  ASSERTS(is_in_range(i,0,(int)m_cp_fn.size()))<<STRMVAR(i);
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

inline bool mscomplex_t::is_unpaired_saddle(int i) const
{
  return (is_saddle(i) &&(is_paired(i) == false));
}

/*---------------------------------------------------------------------------*/

inline bool mscomplex_t::is_unpaired_two_saddle(int i) const
{return ((index(i) == 2) &&(is_paired(i) == false));}

/*---------------------------------------------------------------------------*/

inline bool mscomplex_t::is_unpaired_one_saddle(int i) const
{return ((index(i) == 1) &&(is_paired(i) == false));}

/*---------------------------------------------------------------------------*/

inline int mscomplex_t::surv_extrema(int i) const
{
  ASSERT(is_extrema(i));

  if(is_paired(i) == false)
    return i;

  eGDIR dir = (index(i) == 3)?(ASC):(DES);

  ASSERT(m_conn[dir][pair_idx(i)].size() == 1);

  int j = m_conn[dir][pair_idx(i)].begin()->first;

  ASSERT(!is_paired(j));

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


/*===========================================================================*/



/*===========================================================================*/
inline void mscomplex_t::save(const std::string &f,bool purge_data)
{std::fstream fs(f.c_str(),std::ios::out|std::ios::binary);save(fs,purge_data);}

/*---------------------------------------------------------------------------*/

void mscomplex_t::load(const std::string &f)
{std::fstream fs(f.c_str(),std::ios::in|std::ios::binary);load(fs);}

/*---------------------------------------------------------------------------*/

inline int mscomplex_t::load_merge(const std::string &f1,const std::string &f2)
{std::fstream fs1(f1.c_str(),std::ios::in|std::ios::binary);
 std::fstream fs2(f2.c_str(),std::ios::in|std::ios::binary);
 return load_merge(fs1,fs2);}

/*---------------------------------------------------------------------------*/

inline void mscomplex_t::unmerge_save(const std::string &f1,const std::string &f2)
{std::fstream fs1(f1.c_str(),std::ios::in|std::ios::out|std::ios::binary);
 std::fstream fs2(f2.c_str(),std::ios::in|std::ios::out|std::ios::binary);
 unmerge_save(fs1,fs2);}

/*---------------------------------------------------------------------------*/

inline boost::iterator_range<mscomplex_t::iterator_t> mscomplex_t::cpno_range() const
{return boost::make_iterator_range(iterator_t(0),iterator_t(get_num_critpts()));}

/*===========================================================================*/



/*===========================================================================*/

inline void order_pr_by_cp_index(const mscomplex_t &msc,int &p,int &q)
{if(msc.index(p) < msc.index(q))std::swap(p,q);}

/*===========================================================================*/


}

#endif
