/***************************************************************************
 *   Copyright (C) 2009 by nithin,,,   *
 *   nithin@gauss   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef __GRID_MSCOMPLEX_H_INCLUDED_
#define __GRID_MSCOMPLEX_H_INCLUDED_

#include <boost/enable_shared_from_this.hpp>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/range/iterator_range.hpp>

#include <grid.h>

namespace grid
{
class mscomplex_t:public boost::enable_shared_from_this<mscomplex_t>
{
public:

  rect_t        m_rect;
  rect_t        m_ext_rect;
  rect_t        m_domain_rect;
  int           m_canc_pos;

  cellid_list_t   m_cp_cellid;
  cellid_list_t   m_cp_vertid;
  int_list_t      m_cp_pair_idx;
  int8_list_t     m_cp_index;
  bool_list_t     m_cp_is_cancelled;
  cell_fn_list_t  m_cp_fn;

  int_pair_list_t m_canc_list;

  conn_list_t   m_conn[GDIR_CT];
  conn_list_t  &m_des_conn;
  conn_list_t  &m_asc_conn;

  mfold_list_t  m_mfolds[GDIR_CT];
  mfold_list_t &m_des_mfolds;
  mfold_list_t &m_asc_mfolds;


public:

  mscomplex_t();
  mscomplex_t(rect_t r,rect_t e,rect_t d);
  ~mscomplex_t();
  void clear();

  // mscomplex basic query functions
  inline int       get_num_critpts()              const;
  inline int       pair_idx(int i)                const;
  inline cellid_t  cellid(int i)                  const;
  inline cellid_t  vertid(int i)                  const;
  inline int8_t    index(int i)                   const;
  inline cell_fn_t fn(int i)                      const;
  inline int       surv_extrema(int i)            const;
  inline bool      is_paired(int i)               const;
  inline bool      is_not_paired(int i)           const;
  inline bool      is_extrema(int i)              const;
  inline bool      is_saddle(int i)               const;
  inline bool      is_canceled(int i)             const;
  inline bool      is_unpaired_saddle(int i)      const;
  inline bool      is_unpaired_two_saddle(int i)  const;
  inline bool      is_unpaired_one_saddle(int i)  const;  
  inline cell_fn_t fn_min()                       const; // O(#cp) complexity
  inline cell_fn_t fn_max()                       const; // O(#cp) complexity  
  template <int i>
  inline bool      is_index_i_cp(int cp)          const;

  // iterator range to go over the set of critical points
  typedef boost::counting_iterator<int> iterator_t;
  inline boost::iterator_range<iterator_t> cpno_range() const;

  // functions to create a mscomplex from a dataset
  void  resize(int i);
  void  set_critpt(int i,cellid_t c,char idx,cell_fn_t f,cellid_t vert_cell);
  void  connect_cps(int p, int q,int m=1);

  // simplification related stuff
  void cancel_pair();
  void cancel_pair(int p, int q);
  bool persistence_cmp(int_pair_t p0,int_pair_t p1) const;
  void simplify_pers(double thresh=1.0,bool is_nrm=true,int nmax=0,int nmin=0);

  // geometry collection related stuff
  void collect_mfolds(eGDIR dir, int dim, dataset_ptr_t ds);
  void collect_mfolds(dataset_ptr_t ds);

  // misc functions
  inline std::string info() const;
  inline std::string cp_info (int cp_no) const;
  std::string        cp_conn (int cp_no) const;

  // simplification related things used during outcore processing
  void dir_connect_cps(int p , int q,int m=1);
  void un_simplify();
  void invert_for_collection();
  void uncancel_pair( int p, int q);
  void simplify_pers_outcore(double f_tresh, double f_range);

  // functions to enable outcore merging and merge history traversal etc.
  void save(const std::string &f, bool purge_data=true);
  void load(const std::string &f);
  int  load_merge(const std::string &f1,const std::string &f2);
  void unmerge_save(const std::string &f1,const std::string &f2);

};

inline void order_pr_by_cp_index(const mscomplex_t &msc,int &p,int &q);
}

#include <grid_mscomplex_inl.h>

#endif
