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

#include <set>
#include <map>

#include <boost/noncopyable.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/function.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/iterator/counting_iterator.hpp>

#include <boost/bind.hpp>

#include <grid.h>

namespace grid
{
  struct critpt_t;

  typedef std::vector<cell_fn_t>     cp_fn_list_t;
  typedef n_vector_t<int,2>          int_pair_t;
  typedef std::vector<int_pair_t>    int_pair_list_t;
  typedef std::map<cellid_t,uint>    id_cp_map_t;
  typedef std::vector<critpt_t >     critpt_list_t;

  typedef std::multiset<uint>                 conn_t;
  typedef std::multiset<uint>::iterator       conn_iter_t;
  typedef std::multiset<uint>::const_iterator const_conn_iter_t;
  typedef std::vector<conn_t>                 conn_list_t;


  class mscomplex_t:public boost::enable_shared_from_this<mscomplex_t>
  {
  public:

    cellid_list_t   m_cp_cellid;
    cellid_list_t   m_cp_vertid;
    int_list_t      m_cp_pair_idx;
    char_list_t     m_cp_index;
    bool_list_t     m_cp_is_cancelled;
    cell_fn_list_t  m_cp_fn;

    id_cp_map_t   m_id_cp_map;

    conn_list_t   m_conn[GDIR_CT];
    conn_list_t  &m_des_conn;
    conn_list_t  &m_asc_conn;

    rect_t        m_rect;
    rect_t        m_ext_rect;

  public:

    mscomplex_t(rect_t r,rect_t e);
    ~mscomplex_t();

    inline int  get_num_critpts() const;
    int  resize(int i);

    int  add_critpt(cellid_t c,char i,cell_fn_t f,cellid_t vert_cell);
    void set_critpt(int i,cellid_t c,char idx,cell_fn_t f,cellid_t vert_cell);

    void build_id_cp_map();

    void connect_cps(cellid_t c1,cellid_t c2);
    void connect_cps(int p, int q);

    void dir_connect_cps(cellid_t c1,cellid_t c2);
    void dir_connect_cps(int p , int q);

    void pair_cps(cellid_t c1,cellid_t c2);
    void pair_cps(int p , int q);

    inline char& index(int i);
    inline const char& index(int i) const;

    inline int& pair_idx(int i);
    inline const int& pair_idx(int i) const;

    inline bool is_paired(int i) const;

    inline void set_is_canceled(int i,bool b);
    inline bool is_canceled(int i) const;

    inline cellid_t& cellid(int i);
    inline const cellid_t& cellid(int i) const;

    inline cellid_t& vertid(int i);
    inline const cellid_t& vertid(int i) const;

    inline cell_fn_t& fn(int i);
    inline const cell_fn_t& fn(int i) const;

    inline int surv_extrema(int i) const;
    inline bool is_extrema(int i) const;
    inline bool is_saddle(int i) const;

    inline bool is_unpaired_saddle(int i) const
    {
      return (is_saddle(i) &&(is_paired(i) == false));
    }

    inline bool is_twosaddle_des_valid(int i) const
    {
      return (index(i) == 2 &&(!is_paired(i) || index(pair_idx(i) == 1) ));
    }

    inline bool is_onesaddle_asc_valid(int i) const
    {
      return (index(i) == 1 &&(!is_paired(i) || index(pair_idx(i) == 2) ));
    }

  public:

    void simplify(int_pair_list_t &,double simplification_treshold);

    void un_simplify(const int_pair_list_t &);

    void simplify_un_simplify(double simplification_treshold );

    void invert_for_collection();

    void cancel_pair(int p, int q);
    void uncancel_pair( int p, int q);

    void clear();

    void merge_up  (const mscomplex_t& ,const mscomplex_t& ,const rect_t&);
    void merge_down(mscomplex_t& ,mscomplex_t& ,const rect_t&);

    void write_graph(std::ostream & os) const;
    void write_graph(const std::string & fn) const;

    inline std::string cp_info (int cp_no) const;
    inline std::string cp_conn (int cp_no) const;

    typedef boost::counting_iterator<int> iterator_t;
    typedef boost::function<bool (int)> filter_t;
    typedef bool (mscomplex_t::*memb_filter_t)(int) const;
    typedef boost::filter_iterator<filter_t,iterator_t> fiterator_t;

    inline iterator_t begin() const;
    inline iterator_t end() const;

    inline fiterator_t fbegin(filter_t f) const;
    inline fiterator_t fend(filter_t f) const;

    inline fiterator_t fbegin(memb_filter_t f) const;
    inline fiterator_t fend(memb_filter_t f) const;

  };

  inline void order_pr_by_cp_index(mscomplex_t &msc,int &p,int &q)
  {
    if(msc.index(p) < msc.index(q))
      std::swap(p,q);
  }

  template<typename Titer>
  class mt_producer_t:boost::noncopyable
  {
  protected:
    Titer          m_begin;
    Titer          m_end;
    boost::mutex   m_mutex;

    mt_producer_t(){};

  public:

    typedef Titer iterator;

    class no_more_items:public std::exception
    {
    public:
      no_more_items(){}
      virtual ~no_more_items() throw(){}

      virtual const char*
      what() const throw()
      {return "nothing more to produce";}
    };

    mt_producer_t(Titer begin,Titer end)
    {
      boost::mutex::scoped_lock scoped_lock(m_mutex);

      m_begin = begin;
      m_end   = end;
    }

    Titer next()
    {
      boost::mutex::scoped_lock scoped_lock(m_mutex);

      if(m_begin == m_end)
        throw no_more_items();

      return m_begin++;
    }
  };

  class mt_cp_producer_t:public mt_producer_t<mscomplex_t::fiterator_t>
  {
  public:
    mt_cp_producer_t(mscomplex_const_ptr_t msc, mscomplex_t::memb_filter_t f)
    {
      boost::mutex::scoped_lock scoped_lock(m_mutex);

      m_begin = msc->fbegin(f);
      m_end   = msc->fend(f);
    }

    mt_cp_producer_t(mscomplex_const_ptr_t msc, mscomplex_t::filter_t f)
    {
      boost::mutex::scoped_lock scoped_lock(m_mutex);

      m_begin = msc->fbegin(f);
      m_end   = msc->fend(f);
    }


    inline int count()
    {
      return std::count(m_begin,m_end);
    }
  };

  class cp_producer_t: boost::noncopyable
  {
  public:
    typedef boost::function<bool (mscomplex_const_ptr_t,int)> cp_filter_t;

    boost::mutex           m_mutex;
    int                    m_ni;
    mscomplex_const_ptr_t  m_msc;
    cp_filter_t            m_cp_filter;

    cp_producer_t(mscomplex_const_ptr_t msc,cp_filter_t cf = pass_filter);

    bool next(int & i);
    int  count() const;

    static bool pass_filter(mscomplex_const_ptr_t  msc, int i);
    static bool extrema_filter(mscomplex_const_ptr_t  msc, int i);
    static bool twosaddle_filter(mscomplex_const_ptr_t  msc, int i);
    static bool saddle_filter(mscomplex_const_ptr_t  msc, int i);
    static bool unpaired_cp_filter(mscomplex_const_ptr_t  msc, int i);
    static bool unpaired_saddle_filter(mscomplex_const_ptr_t  msc, int i);

    // Intention:
    //    for multiple threads to pick distinct cp's to work with
    //
    //     Usage:
    // ----main----thread----
    // cp_producer_ptr_t p(new cp_producer_t(msc,cp_filter))
    //
    // ----worker--threads---
    // for(int i ; p->next(i) ;)
    // {
    //   -- each thread works on a unique cp --
    //
    // }

  };

  inline mscomplex_t::iterator_t mscomplex_t::begin() const
  {
    return iterator_t(0);
  }

  inline mscomplex_t::iterator_t mscomplex_t::end() const
  {
    return iterator_t(get_num_critpts());
  }

  inline mscomplex_t::fiterator_t mscomplex_t::fbegin(mscomplex_t::filter_t f) const
  {
    return boost::make_filter_iterator(f,begin(),end());
  }

  inline mscomplex_t::fiterator_t mscomplex_t::fend(mscomplex_t::filter_t f) const
  {
    return boost::make_filter_iterator(f,end(),end());
  }

  inline mscomplex_t::fiterator_t mscomplex_t::fbegin(mscomplex_t::memb_filter_t f) const
  {
    return fbegin(boost::bind(f,this,_1));
  }

  inline mscomplex_t::fiterator_t mscomplex_t::fend(mscomplex_t::memb_filter_t f) const
  {
    return fend(boost::bind(f,this,_1));
  }


  inline int  mscomplex_t::get_num_critpts() const
  {
    return m_cp_cellid.size();
  }

  inline char& mscomplex_t::index(int i)
  {
    try{ASSERT(is_in_range(i,0,m_cp_index.size()));}
    catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

    return m_cp_index[i];
  }

  inline const char& mscomplex_t::index(int i) const
  {
    try{ASSERT(is_in_range(i,0,m_cp_index.size()));}
    catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

    return m_cp_index[i];
  }

  inline int& mscomplex_t::pair_idx(int i)
  {
    try{ASSERT(is_in_range(i,0,m_cp_pair_idx.size()));}
    catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

    return m_cp_pair_idx[i];
  }

  inline const int& mscomplex_t::pair_idx(int i) const
  {
    try{ASSERT(is_in_range(i,0,m_cp_pair_idx.size()));}
    catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

    return m_cp_pair_idx[i];
  }

  inline bool mscomplex_t::is_paired(int i) const
  {
    try{ASSERT(is_in_range(i,0,m_cp_pair_idx.size()));}
    catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

    return (m_cp_pair_idx[i] != -1);
  }

  inline void mscomplex_t::set_is_canceled(int i,bool c)
  {
    try{ASSERT(is_in_range(i,0,m_cp_is_cancelled.size()));}
    catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

    m_cp_is_cancelled[i] = c;
  }

  inline bool mscomplex_t::is_canceled(int i) const
  {
    try{ASSERT(is_in_range(i,0,m_cp_is_cancelled.size()));}
    catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

    return m_cp_is_cancelled[i];
  }

  inline cellid_t& mscomplex_t::cellid(int i)
  {
    try{ASSERT(is_in_range(i,0,m_cp_cellid.size()));}
    catch(assertion_error e)
    {e.push(_FFL).push(SVAR(i)).push(SVAR(m_cp_cellid.size()));throw;}

    return m_cp_cellid[i];
  }

  inline const cellid_t& mscomplex_t::cellid(int i) const
  {
    try{ASSERT(is_in_range(i,0,m_cp_cellid.size()));}
    catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

    return m_cp_cellid[i];
  }

  inline cellid_t& mscomplex_t::vertid(int i)
  {
    try{ASSERT(is_in_range(i,0,m_cp_vertid.size()));}
    catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

    return m_cp_vertid[i];
  }

  inline const cellid_t& mscomplex_t::vertid(int i) const
  {
    try{ASSERT(is_in_range(i,0,m_cp_vertid.size()));}
    catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

    return m_cp_vertid[i];
  }

  inline cell_fn_t& mscomplex_t::fn(int i)
  {
    try{ASSERT(is_in_range(i,0,m_cp_fn.size()));}
    catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

    return m_cp_fn[i];
  }

  inline const cell_fn_t& mscomplex_t::fn(int i) const
  {
    try{ASSERT(is_in_range(i,0,m_cp_fn.size()));}
    catch(assertion_error e){e.push(_FFL).push(SVAR(i));throw;}

    return m_cp_fn[i];
  }

  inline bool mscomplex_t::is_extrema(int i) const
  {
    return (index(i) == 0 || index(i) == 3);
  }

  inline bool mscomplex_t::is_saddle(int i) const
  {
    return (index(i) == 1 || index(i) == 2);
  }

  inline int mscomplex_t::surv_extrema(int i) const
  {
    ASSERT(is_extrema(i));

    if(is_paired(i) == false)
      return i;

    eGDIR dir = (index(i) == 3)?(GDIR_ASC):(GDIR_DES);

    ASSERT(m_conn[dir][pair_idx(i)].size() == 1);

    int j = *m_conn[dir][pair_idx(i)].begin();

    ASSERT(!is_paired(j));

    return j;
  }

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

  inline std::string mscomplex_t::cp_conn (int i) const
  {
    std::stringstream ss;

    ss<<std::endl<<"des = ";

    for(const_conn_iter_t it = m_des_conn[i].begin(); it != m_des_conn[i].end(); ++it)
      ss<<cellid(*it);

    ss<<std::endl<<"asc = ";

    for(const_conn_iter_t it = m_asc_conn[i].begin(); it != m_asc_conn[i].end(); ++it)
      ss<<cellid(*it);

    ss<<std::endl;

    return ss.str();
  }
}


//#include <boost/serialization/array.hpp>
//#include <boost/serialization/base_object.hpp>
//
//namespace boost
//{
//  namespace serialization
//  {
//    template<class Archive>
//    void serialize(Archive & ar, grid::mscomplex_t & g, const unsigned int );
//
//  } // namespace serialization
//}
#endif
