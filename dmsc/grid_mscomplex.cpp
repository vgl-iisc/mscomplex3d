#include <cmath>
#include <queue>
#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>

#include <boost/typeof/typeof.hpp>

#include <grid_mscomplex.h>

using namespace std;

namespace grid
{
  inline std::string edge_to_string(mscomplex_t *msc,int_pair_t e)
  {
    std::stringstream ss;

    ss<<utls::to_string(msc->cellid(e[0]))<<"----"<<utls::to_string(msc->cellid(e[0]));

    return ss.str();
  }

  mscomplex_t::mscomplex_t(rect_t r,rect_t e,rect_t d)
    :m_rect(r),m_ext_rect(e),m_domain_rect(d),
      m_des_conn(m_conn[0]),m_asc_conn(m_conn[1]){}

  mscomplex_t::~mscomplex_t(){clear();
  }

  void mscomplex_t::set_critpt(int i,cellid_t c,char idx,cell_fn_t f,cellid_t v)
  {
    cellid(i) = c;
    vertid(i) = v;
    index(i)  = idx;
    fn(i)     = f;
  }

  int mscomplex_t::add_critpt(cellid_t c,char idx,cell_fn_t f,cellid_t v)
  {
    int i = get_num_critpts();

    ASSERT(m_id_cp_map.count(c) == 0);
    m_id_cp_map.insert(std::make_pair(c,i));

    resize(i+1);
    set_critpt(i,c,idx,f,v);

    return (i);
  }

  int  mscomplex_t::resize(int i)
  {
    m_cp_cellid.resize(i,cellid_t(-1,-1,-1));
    m_cp_vertid.resize(i,cellid_t(-1,-1,-1));
    m_cp_index.resize(i,-1);
    m_cp_pair_idx.resize(i,-1);
    m_cp_is_cancelled.resize(i,false);
    m_cp_fn.resize(i);
    m_des_conn.resize(i);
    m_asc_conn.resize(i);
  }

  void mscomplex_t::build_id_cp_map()
  {
    for(int i = 0 ; i < get_num_critpts(); ++i)
    {
      try
      {
        ASSERT(cellid(i) != cellid_t(-1,-1,-1));
        ASSERT(m_id_cp_map.count(cellid(i)) == 0);
        m_id_cp_map.insert(std::make_pair(cellid(i),i));
      }
      catch(assertion_error e)
      {
        e.push(_FFL);
        e.push(SVAR(cp_info(i)));

        throw;
      }
    }
  }

  void mscomplex_t::connect_cps(cellid_t c0,cellid_t c1)
  {
    try
    {
      ASSERT(m_id_cp_map.count(c0) == 1);
      ASSERT(m_id_cp_map.count(c1) == 1);

      connect_cps(m_id_cp_map[c0],m_id_cp_map[c1]);
    }
    catch(assertion_error e)
    {
      e.push(_FFL);
      e.push(SVAR2(c0,m_id_cp_map.count(c0)));
      e.push(SVAR2(c1,m_id_cp_map.count(c1)));
      throw;
    }
  }

  void mscomplex_t::connect_cps(int p, int q)
  {
    try
    {
      order_pr_by_cp_index(*this,p,q);

      ASSERT(index(p) == index(q)+1);

      // if a d-cp hits a d+-1 cp and the d+-1 cp is paired
      // then the connection is useful iff the dimension of the pair is d

      ASSERT(!(is_paired(p) && index(pair_idx(p))!= index(q)));
      ASSERT(!(is_paired(q) && index(pair_idx(q))!= index(p)));
      ASSERT(m_des_conn[p].count(q) == m_asc_conn[q].count(p));

      if(m_des_conn[p].count(q) == 2)
        return;

      m_des_conn[p].insert(q);
      m_asc_conn[q].insert(p);
    }
    catch(assertion_error e)
    {
      e.push(_FFL);

      e.push(SVAR(cp_info(p)));
      if(is_paired(p))
        e.push(SVAR(cp_info(pair_idx(p))));

      e.push(SVAR(cp_info(q)));
      if(is_paired(q))
        e.push(SVAR(cp_info(pair_idx(q))));


      throw;
    }
  }

  void mscomplex_t::dir_connect_cps(cellid_t c1,cellid_t c2)
  {
    ASSERT(m_id_cp_map.count(c1) == 1);
    ASSERT(m_id_cp_map.count(c2) == 1);

    dir_connect_cps(m_id_cp_map[c1],m_id_cp_map[c2]);
  }

  void mscomplex_t::dir_connect_cps(int p, int q)
  {
    try
    {
      ASSERT(is_paired(p) != is_paired(q));
      ASSERT(abs(index(p)-index(q)) == 1);

      if(is_paired(q))
        std::swap(p,q);

      conn_t &conn = (index(p) > index(q))?(m_des_conn[p]):(m_asc_conn[p]);

      if(conn.count(q) == 0)
        conn.insert(q);
    }
    catch(assertion_error e)
    {
      e.push(_FFL);
      e.push(SVAR(cp_info(p)));
      e.push(SVAR(cp_info(q)));
      throw;
    }
  }

  void mscomplex_t::pair_cps(cellid_t c1,cellid_t c2)
  {
    ASSERT(m_id_cp_map.count(c1) == 1);
    ASSERT(m_id_cp_map.count(c2) == 1);

    pair_cps(m_id_cp_map[c1],m_id_cp_map[c2]);
  }

  void mscomplex_t::pair_cps(int p, int q)
  {
    pair_idx(p) = q;
    pair_idx(q) = p;
  }

  void mscomplex_t::cancel_pair ( int p, int q)
  {
    order_pr_by_cp_index(*this,p,q);

    try
    {
      ASSERT(index(p) == index(q)+1);
      ASSERT(pair_idx(p) == q);
      ASSERT(pair_idx(q) == p);
      ASSERT(is_canceled(p) == false);
      ASSERT(is_canceled(q) == false);
      ASSERT(m_des_conn[p].count(q) == 1);
      ASSERT(m_asc_conn[q].count(p) == 1);

      conn_iter_t i,j;

      m_des_conn[p].erase(q);
      m_asc_conn[q].erase(p);

      // cps in lower of u except l
      for(i = m_des_conn[p].begin();i != m_des_conn[p].end();++i)
        for(j = m_asc_conn[q].begin();j != m_asc_conn[q].end();++j)
        {
          ASSERT(is_canceled(*i) == false);
          ASSERT(is_canceled(*j) == false);

          connect_cps(*i,*j);
        }

      for(j = m_des_conn[p].begin();j != m_des_conn[p].end();++j)
        m_asc_conn[*j].erase(p);

      for(j = m_asc_conn[p].begin();j != m_asc_conn[p].end();++j)
        m_des_conn[*j].erase(p);

      for(j = m_des_conn[q].begin();j != m_des_conn[q].end();++j)
        m_asc_conn[*j].erase(q);

      for(j = m_asc_conn[q].begin();j != m_asc_conn[q].end();++j)
        m_des_conn[*j].erase(q);

      set_is_canceled(p,true);
      set_is_canceled(q,true);

      m_asc_conn[p].clear();
      m_des_conn[q].clear();
    }
    catch (assertion_error ex)
    {
      ex.push(_FFL).push(SVAR(cp_info(p))).push(SVAR(cp_info(q)));
      ex.PUSHVAR(m_des_conn[p].count(q));
      ex.PUSHVAR(m_asc_conn[q].count(p));
      throw;
    }
  }

  void mscomplex_t::uncancel_pair(int p, int q)
  {
    order_pr_by_cp_index(*this,p,q);

    try
    {
      ASSERT(is_canceled(p) == true && is_canceled(q) == true);
      ASSERT(index(p) == index(q)+1);
      ASSERT(pair_idx(p) == q && pair_idx(q) == p);

      set_is_canceled(p,false);
      set_is_canceled(q,false);

      conn_iter_t i,j;

      for(int d = 0 ; d <2 ; ++d)
      {
        int ed = (d == 0)?(p):(q);

        conn_t old_conn(m_conn[d][ed].begin(),m_conn[d][ed].end());

        m_conn[d][ed].clear();

        for(i = old_conn.begin();i != old_conn.end() ; ++i)
        {
          if(is_paired(*i) == false)
          {
            dir_connect_cps(ed,*i);
            continue;
          }

          int r = pair_idx(*i);

          if(index(ed) != index(r))
            continue;

          try
          {
            ASSERT(is_canceled(*i) ==false && is_canceled(r) ==false);
            ASSERT(abs(index(*i) - index(r)) == 1);
            ASSERT(pair_idx( r) == *i && pair_idx(*i) ==  r);
          }
          catch (assertion_error ex)
          {
            ex.push(_FFL).push(SVAR(cp_info(r))).push(SVAR(cp_info(*i)));
            ex.push(SVAR(cp_info(ed)));
            throw;
          }

          try
          {
            for(j = m_conn[d][r].begin(); j!= m_conn[d][r].end() ; ++j )
              dir_connect_cps(ed,*j);
          }
          catch(assertion_error ex)
          {
            ex.push(_FFL)
              .push("failed to connect ed to *j via pair (*i,r)")
              .push(SVAR(cp_info(ed)))
              .push(SVAR(cp_info(*i)))
              .push(SVAR(cp_info(r)))
              .push(SVAR(cp_info(*j)));

            if(is_paired(*j))
              ex.push(SVAR(cp_info(pair_idx(*j))));
            throw;
          }
        }
      }
    }
    catch (assertion_error ex)
    {
      ex.push(_FFL).push(SVAR(cp_info(p))).push(SVAR(cp_info(q)));
      throw;
    }
  }

  void mscomplex_t::merge_down
      (mscomplex_t& msc1,
       mscomplex_t& msc2,
       const rect_t& bnd)
  {
    mscomplex_t *msc_arr[] = {&msc1,&msc2};

    try
    {
      ASSERT(bnd.eff_dim() == gc_grid_dim-1);
      ASSERT(m_ext_rect.intersection(bnd) == bnd);

      rect_t ixn = m_rect.intersection(bnd);

      for(cellid_t c = ixn.upper_corner() ; c[2] >= ixn[2][0]; --c[2])
      {
        for(c[1] = ixn[1][1] ; c[1] >= ixn[1][0]; --c[1])
        {
          for(c[0] = ixn[0][1] ; c[0] >= ixn[0][0]; --c[0])
          {
            if (m_id_cp_map.count(c) == 0)
              continue;

            int p = m_id_cp_map[c];

            if(!is_paired(p))
              continue;

            int q =  pair_idx(p);

            if(bnd.contains(cellid(q)))
              continue;

            uncancel_pair(p,q);
          }
        }
      }
    }
    catch(assertion_error e)
    {
      e.push(_FFL);
      e.push(SVAR(bnd));
      e.push(SVAR2(m_rect,m_ext_rect));
      e.push(SVAR2(msc1.m_rect,msc1.m_ext_rect));
      e.push(SVAR2(msc2.m_rect,msc2.m_ext_rect));
      throw;
    }

    int_list_t cpi_to_mcpi(get_num_critpts());

    for(int m = 0 ; m < 2; ++m)
    {
      mscomplex_t &msc = *msc_arr[m];

      for(int cpi = 0 ; cpi < get_num_critpts(); ++cpi)
      {
        if(msc.m_id_cp_map.count(cellid(cpi)) == 0 )
          cpi_to_mcpi[cpi] = -1;
        else
          cpi_to_mcpi[cpi] = msc.m_id_cp_map[cellid(cpi)];
      }

      for(int p = 0 ; p < get_num_critpts(); ++p)
      {
        if(cpi_to_mcpi[p] == -1 )
          continue;

        int mp       = cpi_to_mcpi[p];

        msc.m_des_conn[mp].clear();
        msc.m_asc_conn[mp].clear();

        if(is_paired(p) == false)
          continue;

        int q = pair_idx(p);

        if(cpi_to_mcpi[q] != -1)
          continue;

        cpi_to_mcpi[q] = msc.add_critpt(cellid(q),index(q),fn(q),vertid(q));
        msc.pair_cps(mp,cpi_to_mcpi[q]);
      }

      for(int i = 0 ; i < get_num_critpts(); ++i)
      {
        if(cpi_to_mcpi[i] == -1 )
          continue;

        if(is_paired(i) == false)
          continue;

        int d = (index(pair_idx(i)) - index(i) + 1)/2;

        if(d == 0 && (msc.m_rect.contains(cellid(i))))
          continue;

        if(d == 1 && (msc.m_rect.contains(cellid(pair_idx(i))) == false))
          continue;

        for(conn_iter_t j = m_conn[d][i].begin();j != m_conn[d][i].end();++j)
        {
          try
          {
            ASSERT(is_paired(*j) == false);

            if(cpi_to_mcpi[*j] == -1)
              cpi_to_mcpi[*j] = msc.add_critpt(cellid(*j),index(*j),fn(*j),vertid(*j));

            msc.dir_connect_cps(cpi_to_mcpi[i],cpi_to_mcpi[*j]);
          }
          catch(assertion_error e)
          {
            e.push(_FFL);
            e.push(SVAR(cp_info(i)));
            e.push(SVAR(cp_info(pair_idx(i))));
            e.push(SVAR(cp_info(*j)));
            e.push(SVAR2(msc.m_rect,msc.m_ext_rect));
            e.push(SVAR2(m_rect,m_ext_rect));
            throw;
          }
        }
      }
    }
  }

  void mscomplex_t::clear()
  {
    m_cp_cellid.clear();
    m_cp_vertid.clear();
    m_cp_pair_idx.clear();
    m_cp_index.clear();
    m_cp_is_cancelled.clear();
    m_cp_fn.clear();
    m_id_cp_map.clear();
    m_des_conn.clear();
    m_asc_conn.clear();
  }

  inline bool is_valid_canc_edge(const mscomplex_t &msc,const std::vector<bool> &is_inc_ext, int_pair_t e )
  {
    order_pr_by_cp_index(msc,e[0],e[1]);

    if(msc.is_canceled(e[0])||msc.is_canceled(e[1]))
      return false;

    if(msc.is_paired(e[0]) || msc.is_paired(e[1]))
      return false;

    if(is_inc_ext[e[0]] && is_inc_ext[e[1]])
      return false;

    if(msc.m_domain_rect.isOnBoundry(msc.cellid(e[0])) !=
       msc.m_domain_rect.isOnBoundry(msc.cellid(e[1])))
      return false;

    ASSERT(msc.m_des_conn[e[0]].count(e[1]) == msc.m_asc_conn[e[1]].count(e[0]));

    if(msc.m_des_conn[e[0]].count(e[1]) != 1)
      return false;

    return true;
  }

  inline bool is_epsilon_persistent(const mscomplex_t &msc,int_pair_t e )
  {
    return (msc.vertid(e[0]) == msc.vertid(e[1]));
  }

  inline int get_num_new_edges(const mscomplex_t &msc, int_pair_t pr)
  {
    order_pr_by_cp_index(msc,pr[0],pr[1]);

    int pd = msc.m_conn[GDIR_DES][pr[0]].size();
    int pa = msc.m_conn[GDIR_ASC][pr[0]].size();

    int qd = msc.m_conn[GDIR_DES][pr[1]].size();
    int qa = msc.m_conn[GDIR_ASC][pr[1]].size();

    return (pd - 1)*(qa - 1) - (pd + qd + pa +qa -1);
  }

  inline bool fast_persistence_lt(const mscomplex_t &msc, int_pair_t p0, int_pair_t p1)
  {
    bool eps_p0 = is_epsilon_persistent(msc,p0);
    bool eps_p1 = is_epsilon_persistent(msc,p1);

    if( eps_p0 != eps_p1)
      return eps_p0;

    return (get_num_new_edges(msc,p0) < get_num_new_edges(msc,p1));
  }

  inline cell_fn_t get_persistence(const mscomplex_t & msc,int_pair_t e)
  {
    return std::abs(msc.fn(e[0]) - msc.fn(e[1]));
  }

  inline bool persistence_lt(const mscomplex_t &msc, int_pair_t p0, int_pair_t p1)
  {/*
    bool eps_p0 = is_epsilon_persistent(msc,p0);
    bool eps_p1 = is_epsilon_persistent(msc,p1);

    if( eps_p0 != eps_p1)
      return eps_p0;

    cell_fn_t d0 = get_persistence(msc,p0);
    cell_fn_t d1 = get_persistence(msc,p1);

    if(d0 != d1)
      return d0 < d1;

    return p0 < p1;*/

    order_pr_by_cp_index(msc,p0[0],p0[1]);
    order_pr_by_cp_index(msc,p1[0],p1[1]);

    cellid_t v00 = msc.vertid(p0[0]);
    cellid_t v01 = msc.vertid(p0[1]);
    cellid_t v10 = msc.vertid(p1[0]);
    cellid_t v11 = msc.vertid(p1[1]);

    cellid_t c00 = msc.cellid(p0[0]);
    cellid_t c01 = msc.cellid(p0[1]);
    cellid_t c10 = msc.cellid(p1[0]);
    cellid_t c11 = msc.cellid(p1[1]);

    if( (v00 == v01 ) != (v10 == v11))
      return (v00 == v01 );

    if( (v00 == v01 ) &&(v10 == v11))
    {
      if(v00 == v10)
      {
        if(c00 != c10)
          return c00 < c10;
        else
          return c01 < c11;
      }
      else
      {
        return (v00 < v10);
      }
    }

    cell_fn_t f00 = msc.fn(p0[0]);
    cell_fn_t f01 = msc.fn(p0[1]);
    cell_fn_t f10 = msc.fn(p1[0]);
    cell_fn_t f11 = msc.fn(p1[1]);

    cell_fn_t d1 = std::abs(f01-f00);
    cell_fn_t d2 = std::abs(f11-f10);

    if(d1 != d2)
      return d1 < d2;

    if(c00 != c10)
      return c00 < c10;

    return c01 < c11;
  }

  inline bool is_within_treshold(const mscomplex_t & msc,int_pair_t e,cell_fn_t t)
  {
    return (is_epsilon_persistent(msc,e) || get_persistence(msc,e) < t);
  }

  inline void make_is_inc_ext(const mscomplex_t &msc, vector<bool> &is_inc_on_ext)
  {
    is_inc_on_ext.resize(msc.get_num_critpts(),false);

    for(uint i = 0 ;i < msc.get_num_critpts();++i)
    {
      if(!msc.is_paired(i) || msc.is_canceled(i)) continue;

      conn_iter_t db = msc.m_des_conn[i].begin();
      conn_iter_t de = msc.m_des_conn[i].end();

      conn_iter_t ab = msc.m_asc_conn[i].begin();
      conn_iter_t ae = msc.m_asc_conn[i].end();

      for(;db!=de;++db) is_inc_on_ext[*db] = true;
      for(;ab!=ae;++ab) is_inc_on_ext[*ab] = true;
    }
  }

  void mscomplex_t::simplify(double f_tresh, double f_range)
  {
    BOOST_AUTO(cmp,bind(fast_persistence_lt,boost::cref(*this),_2,_1));

    priority_queue<int_pair_t,int_pair_list_t,typeof(cmp)> pq(cmp);

    if(f_range <= 0)
    {
      f_range = *max_element(m_cp_fn.begin(),m_cp_fn.end()) -
          *min_element(m_cp_fn.begin(),m_cp_fn.end());
    }

    f_tresh *= f_range;

    vector<bool> is_inc_ext;

    make_is_inc_ext(*this,is_inc_ext);

    for(uint i = 0 ;i < get_num_critpts();++i)
    {
      for(conn_iter_t j = m_des_conn[i].begin();j != m_des_conn[i].end() ;++j)
      {
        int_pair_t pr(i,*j);

        if(is_valid_canc_edge(*this,is_inc_ext,pr) && is_within_treshold(*this,pr,f_tresh))
          pq.push(pr);
      }
    }

    uint num_cancellations = 0;

    while (pq.size() !=0)
    {
      int_pair_t pr = pq.top();

      pq.pop();

      if(is_valid_canc_edge(*this,is_inc_ext,pr) == false)
        continue;

      int p = pr[0],q = pr[1];

      order_pr_by_cp_index(*this,p,q);

      pair_cps(p,q);

      cancel_pair(p,q);

      num_cancellations++;

      m_canc_list.push_back(pr);

      for(conn_iter_t i = m_des_conn[p].begin();i != m_des_conn[p].end();++i)
        for(conn_iter_t j = m_asc_conn[q].begin();j != m_asc_conn[q].end();++j)
        {
          int_pair_t pr(*i,*j);

          if(is_valid_canc_edge(*this,is_inc_ext,pr) && is_within_treshold(*this,pr,f_tresh))
            pq.push(pr);
        }
    }
    cout<<"num_cancellations    ::"<<(num_cancellations)<<endl;
  }

  void mscomplex_t::un_simplify()
  {
    typedef int_pair_list_t::const_reverse_iterator revit_t;

    for(revit_t it = m_canc_list.rbegin();it != m_canc_list.rend() ; ++it)
      uncancel_pair((*it)[0],(*it)[1]);

    m_canc_list.clear();
  }

  void mscomplex_t::invert_for_collection()
  {
    for(uint i = 0 ; i < get_num_critpts(); ++i)
    {
      if(is_paired(i)== true)
        continue;

      m_des_conn[i].clear();
      m_asc_conn[i].clear();
      continue;
    }

    for(int i = 0 ; i < get_num_critpts(); ++i)
    {
      if(is_paired(i) == false)
        continue;

      ASSERT(is_paired(i) && is_paired(pair_idx(i))== true);
      ASSERT(abs(index(i)- index(pair_idx(i))) == 1);
      ASSERT(pair_idx(pair_idx(i))  == i);

      int dir = (index(i) > index(pair_idx(i)))?(0):(1);

      for(conn_iter_t j  = m_conn[dir][i].begin(); j != m_conn[dir][i].end(); ++j)
      {
        ASSERT(is_paired(*j) == false);

        m_conn[dir^1][*j].insert(i);
      }

//      m_conn[dir][i].clear();
    }
  }

  void mscomplex_t::write_graph(std::ostream & os) const
  {
    using namespace std;

    os<<"# Num Cps"<<endl;
    os<<get_num_critpts()<<endl;

    os<<"# grid "<<endl;
    os<<m_rect<<endl;

    os<<"# SL.No idx isPaired pairIdx cpCell vertCell fn"<<endl;

    for(int i = 0 ; i < get_num_critpts();++i)
    {
      os<<i<<" ";
      os<<(int)index(i)<<" ";
      os<<(bool)is_paired(i)<<" ";
      os<<(int)pair_idx(i)<<" ";
      os<<cellid(i)<<" ";
      os<<vertid(i)<<" ";
      os<<fn(i)<<" ";
      os<<endl;
    }

    os<<"# SL.No  numDes numAsc connList"<<std::endl;

    for(uint i = 0 ; i < get_num_critpts();++i)
    {
      os<<(int)i<<" ";
      os<<(int)m_des_conn[i].size()<<" ";
      os<<(int)m_asc_conn[i].size()<<" ";

      for(conn_iter_t j = m_des_conn[i].begin(); j != m_des_conn[i].end(); ++j)
        os<<*j<<" ";

      for(conn_iter_t j = m_asc_conn[i].begin(); j != m_asc_conn[i].end(); ++j)
        os<<*j<<" ";

      os<<endl;
    }
  }

  template<typename T>
  inline void bin_write_vec(std::ostream &os, std::vector<T> &v)
  {os.write((const char*)(const void*)v.data(),v.size()*sizeof(T));v.clear();}

  template<typename T>
  inline void bin_write(std::ostream &os, const T &v)
  {os.write((const char*)(const void*)&v,sizeof(T));}

  template<typename T>
  inline void bin_read_vec(std::istream &is, std::vector<T> &v,int n)
  {v.resize(n);is.read((char*)(void*)v.data(),n*sizeof(T));}

  template<typename T>
  inline void bin_read(std::istream &is, const T &v)
  {is.read((char*)(void*)&v,sizeof(T));}

  void mscomplex_t::stow(std::ostream &os)
  {
    int N = get_num_critpts();

    bin_write(os,N);
    bin_write(os,m_rect);
    bin_write(os,m_ext_rect);
    bin_write(os,m_domain_rect);

    bin_write_vec(os,m_cp_cellid);
    bin_write_vec(os,m_cp_vertid);
    bin_write_vec(os,m_cp_pair_idx);
    bin_write_vec(os,m_cp_index);
    bin_write_vec(os,m_cp_is_cancelled);
    bin_write_vec(os,m_cp_fn);

    int_list_t nconn(2*N);
    int_list_t adj;

    for(int i = 0 ; i < N; ++i)
    {
      nconn[2*i]   = m_des_conn[i].size();
      nconn[2*i+1] = m_asc_conn[i].size();

      std::copy(m_des_conn[i].begin(),m_des_conn[i].end(),back_inserter(adj));
      std::copy(m_asc_conn[i].begin(),m_asc_conn[i].end(),back_inserter(adj));
    }

    bin_write(os,(int)adj.size());
    bin_write_vec(os,nconn);
    bin_write_vec(os,adj);

    m_id_cp_map.clear();
    m_conn[0].clear();
    m_conn[1].clear();
  }

  void mscomplex_t::load(std::istream &is)
  {
    clear();

    int N,NC;
    int_list_t nconn,adj;

    bin_read(is,N);
    bin_read(is,m_rect);
    bin_read(is,m_ext_rect);
    bin_read(is,m_domain_rect);

    bin_read_vec(is,m_cp_cellid,N);
    bin_read_vec(is,m_cp_vertid,N);
    bin_read_vec(is,m_cp_pair_idx,N);
    bin_read_vec(is,m_cp_index,N);
    bin_read_vec(is,m_cp_is_cancelled,N);
    bin_read_vec(is,m_cp_fn,N);

    bin_read(is,NC);
    bin_read_vec(is,nconn,2*N);
    bin_read_vec(is,adj,NC);

    int_list_t::iterator a,b,c = adj.begin();

    m_des_conn.resize(N);
    m_asc_conn.resize(N);

    cout<<g_timer.getElapsedTimeInMilliSec()
        <<"\t: read adj and off"<<endl;

    for(int i = 0 ; i < N; ++i)
    {
      a = c;
      b = a + (nconn[2*i]);
      c = b + (nconn[2*i+1]);

      m_des_conn[i].insert(a,b);
      m_asc_conn[i].insert(b,c);
    }

    cout<<g_timer.getElapsedTimeInMilliSec()
        <<"\t: setup conn"<<endl;

    build_id_cp_map();

    cout<<g_timer.getElapsedTimeInMilliSec()
        <<"\t: built Id cp"<<endl;

  }

  inline rect_t get_ixn(rect_t r1,rect_t r2,rect_t e1,rect_t e2)
  {
    rect_t ir = r1.intersection(r2);
    rect_t ie = e1.intersection(e2);

    ASSERT((r1.lc() == ir.lc() && r2.uc() == ir.uc()) ||
           (r2.lc() == ir.lc() && r1.uc() == ir.uc()));

    ASSERT((e1.lc() == ie.lc() && e2.uc() == ie.uc()) ||
           (e2.lc() == ie.lc() && e1.uc() == ie.uc()));

    cellid_t s = ir.span();

    int d = (s[0] == 0)?(0):((s[1] == 0)?(1):(2));

    ASSERT(s[d] == 0);

    ie[d] = ir[d];

    ASSERT(ir.eff_dim() == gc_grid_dim-1);
    ASSERT(ie.eff_dim() == gc_grid_dim-1);

    return ie;
  }

  inline int get_idx_maps(int_list_t & idx_map1,int_list_t & idx_map2,
                           rect_t & ixn,int_marray_t & ixn_idx,
                           cellid_list_t &cl1,cellid_list_t &cl2,
                           bool_list_t &ic1,bool_list_t &ic2)
  {
    ixn_idx.resize(ixn.span()+1);
    ixn_idx.reindex(ixn.lc());
    memset(ixn_idx.data(),-1,num_cells(ixn)*sizeof(int));

    int pos = 0,N1=cl1.size(),N2 = cl2.size();
    idx_map1.resize(N1,-1);idx_map2.resize(N2,-1);

    for(int i =0; i < N1; ++i)
      if(!ic1[i])
      {
        if(ixn.contains(cl1[i]))
          ixn_idx(cl1[i]) = pos;

        idx_map1[i] = pos++;
      }

    for(int i =0; i < N2; ++i)
      if(!ic2[i])
      {
        idx_map2[i] = (ixn.contains(cl2[i]))?(ixn_idx(cl2[i])):(pos++);
      }

    return pos;
  }

  inline void copy_cp_info(mscomplex_t &msc,int_list_t & idx_map,
  cellid_list_t  &cl,cellid_list_t  &vl,int_list_t &pi,char_list_t &ci,cell_fn_list_t &fn)
  {
    int N = idx_map.size();

    for(int i = 0 ; i < N; ++i)
    {
      int j = idx_map[i];

      if(j >=0)
      {
        msc.m_cp_cellid[j]       = cl[i];
        msc.m_cp_vertid[j]       = vl[i];
        msc.m_cp_index[j]        = ci[i];
        msc.m_cp_is_cancelled[j] = false;
        msc.m_cp_fn[j]           = fn[i];

        if(pi[i] != -1)
          msc.m_cp_pair_idx[j] = idx_map[pi[i]];

      }
    }
  }

  inline void copy_adj_info(mscomplex_t & msc,int_list_t & idx_map,
                            int_list_t nconn,int_list_t adj)
  {
    int N1 = nconn.size()/2;

    int_list_t::iterator a,b,c = adj.begin();

    for(int i = 0 ; i < N1; ++i)
    {
      a = c;
      b = a + (nconn[2*i]);
      c = b + (nconn[2*i+1]);

      int j = idx_map[i];

      if( j == -1)
        continue;

      for(;a != b; ++a)
      {
        ASSERT(is_in_range(idx_map[*a],0,msc.get_num_critpts()));
        msc.m_des_conn[j].insert(idx_map[*a]);
      }

      for(;b != c; ++b)
      {
        ASSERT(is_in_range(idx_map[*b],0,msc.get_num_critpts()));
        msc.m_asc_conn[j].insert(idx_map[*b]);
      }
    }
  }

  inline void copy_adj_info(mscomplex_t & msc,int_list_t & idx_map,
  int_list_t nconn,int_list_t adj,rect_t &ixn,cellid_list_t & cl)
  {
    int N1 = nconn.size()/2;

    int_list_t::iterator a,b,c = adj.begin();

    for(int i = 0 ; i < N1; ++i)
    {
      a = c;
      b = a + (nconn[2*i]);
      c = b + (nconn[2*i+1]);

      int j = idx_map[i];

      if( j == -1)
        continue;

      bool in_ixn = ixn.contains(cl[i]);

      for(;a != b; ++a)
        if(!in_ixn || !ixn.contains(cl[*a]))
        {
          ASSERT(is_in_range(idx_map[*a],0,msc.get_num_critpts()));
          msc.m_des_conn[j].insert(idx_map[*a]);
        }

      for(;b != c; ++b)
        if(!in_ixn || !ixn.contains(cl[*b]))
        {
          ASSERT(is_in_range(idx_map[*b],0,msc.get_num_critpts()));
          msc.m_asc_conn[j].insert(idx_map[*b]);
        }
    }
  }

  void mscomplex_t::load_merge(std::istream &is1,std::istream &is2)
  {
    int_list_t    idx_map1,idx_map2;
    rect_t        ixn;
    int_marray_t  ixn_idx;

    {
      int            N1,N2,NC1,NC2;
      rect_t         r1,r2,e1,e2,d1,d2;
      int_list_t     nconn1,nconn2,adj1,adj2;
      cellid_list_t  cl1,cl2;
      cellid_list_t  vl1,vl2;
      int_list_t     pi1,pi2;
      char_list_t    ci1,ci2;
      bool_list_t    ic1,ic2;
      cell_fn_list_t fn1,fn2;

      bin_read(is1,N1);               bin_read(is2,N2);
      bin_read(is1,r1);               bin_read(is2,r2);
      bin_read(is1,e1);               bin_read(is2,e2);
      bin_read(is1,d1);               bin_read(is2,d2);

      ixn = get_ixn(r1,r2,e1,e2);

      bin_read_vec(is1,cl1,N1);       bin_read_vec(is2,cl2,N2);
      bin_read_vec(is1,vl1,N1);       bin_read_vec(is2,vl2,N2);
      bin_read_vec(is1,pi1,N1);       bin_read_vec(is2,pi2,N2);
      bin_read_vec(is1,ci1,N1);       bin_read_vec(is2,ci2,N2);
      bin_read_vec(is1,ic1,N1);       bin_read_vec(is2,ic2,N2);
      bin_read_vec(is1,fn1,N1);       bin_read_vec(is2,fn2,N2);

      bin_read(is1,NC1);              bin_read(is2,NC2);
      bin_read_vec(is1,nconn1,2*N1);  bin_read_vec(is2,nconn2,2*N2);
      bin_read_vec(is1,adj1,NC1);     bin_read_vec(is2,adj2,NC2);

      resize(get_idx_maps(idx_map1,idx_map2,ixn,ixn_idx,cl1,cl2,ic1,ic2));
      copy_cp_info(*this,idx_map1,cl1,vl1,pi1,ci1,fn1);
      copy_cp_info(*this,idx_map2,cl2,vl2,pi2,ci2,fn2);

      copy_adj_info(*this,idx_map1,nconn1,adj1);
      copy_adj_info(*this,idx_map2,nconn2,adj2,ixn,cl2);
    }

    int n_cp = 0;

    for(cellid_t c = ixn.lc() ; c[2] <= ixn[2][1]; ++c[2])
    {
      for(c[1] = ixn[1][0] ; c[1] <= ixn[1][1]; ++c[1])
      {
        for(c[0] = ixn[0][0] ; c[0] <= ixn[0][1]; ++c[0])
        {
          int p = ixn_idx(c);

          if(p < 0) continue;

          int q = m_cp_pair_idx[p];

          if(q < 0) continue;

          ASSERT(pair_idx(pair_idx(p)) == p);

          if(ixn.contains(cellid(q))) continue;

          cancel_pair(p,q);

          ++n_cp;
        }
      }
    }
  }


  void mscomplex_t::write_graph(const std::string &fn) const
  {
    std::fstream os(fn.c_str(),std::ios::out);

    ensure(os.is_open(),"failed to open file");

    write_graph(os);

    os.close();
  }

}


//#include <boost/serialization/vector.hpp>
//#include <boost/serialization/array.hpp>
//#include <boost/serialization/map.hpp>
//#include <boost/serialization/set.hpp>
//#include <boost/serialization/base_object.hpp>
//#include <boost/serialization/binary_object.hpp>
//
//#include <boost/archive/binary_iarchive.hpp>
//#include <boost/archive/binary_oarchive.hpp>
//
//using namespace grid;
//
//namespace boost
//{
//  namespace serialization
//  {
//    template<class Archive>
//    void serialize(Archive & ar, rect_range_t & r, const unsigned int )
//    {
//      typedef boost::array<rect_range_t::value_type,rect_range_t::static_size>
//          rect_range_base_t;
//
//      ar & boost::serialization::base_object<rect_range_base_t>(r);
//    }
//
//    template<class Archive>
//    void serialize(Archive & ar, rect_point_t & p, const unsigned int )
//    {
//      typedef boost::array<rect_point_t::value_type,rect_point_t::static_size>
//          rect_point_base_t;
//
//      ar & boost::serialization::base_object<rect_point_base_t>(p);
//    }
//
//    template<class Archive>
//    void serialize(Archive & ar, rect_t & r, const unsigned int )
//    {
//      typedef boost::array<rect_t::value_type,rect_t::static_size>
//          rect_base_t;
//
//      ar & boost::serialization::base_object<rect_base_t>(r);
//    }
//
//    template<class Archive>
//    void serialize(Archive & ar, critpt_t & c, const unsigned int )
//    {
//      ar & c.cellid;
//      for(uint dir = 0 ;dir <0 ;++dir)
//        ar& c.conn[dir];
//      ar & c.is_paired;
//      ar & c.isCancelled;
//      ar & c.pair_idx;
//      ar & c.fn;
//    }
//
//
//    template<class Archive>
//    void serialize(Archive & ar, mscomplex_t & g, const unsigned int )
//    {
//      ar & g.m_rect;
//      ar & g.m_ext_rect;
//      ar & g.m_id_cp_map;
//      ar & g.m_cps;
//    }
//
//    //    template<class Archive>
//    //    void serialize(Archive & ar, GridDataset & ds, const unsigned int )
//    //    {
//    //       ar & ds.m_rect;
//    //       ar & ds.m_ext_rect;
//    //
//    //       GridDataset::rect_size_t ext_sz = ds.m_ext_rect.size();
//    //       uint num_data_items = (ext_sz[0]+1)*(ext_sz[1]+1);
//    //
//    //       if(Archive::is_loading::value)
//    //         ds.init(NULL);
//    //
//    //       ar & make_binary_object(ds.(*m_cell_flags).data(),num_data_items*sizeof(GridDataset::cell_flag_t));
//    //       ar & make_binary_object(ds.m_cell_pairs.data(),num_data_items*sizeof(GridDataset::cellid_t));
//    //    }
//  }
//}
//
////// without the explicit instantiations below, the program will
////// fail to link for lack of instantiantiation of the above function
////// The impls are visible only in this file to save compilation time..
////
////template void boost::serialization::serialize<boost::archive::text_iarchive>(
////    boost::archive::text_iarchive & ar,
////    GridDataset & g,
////    const unsigned int file_version
////);
////
////template void boost::serialization::serialize<boost::archive::text_oarchive>(
////    boost::archive::text_oarchive & ar,
////    GridDataset & g,
////    const unsigned int file_version
////);
//
//
//// without the explicit instantiations below, the program will
//// fail to link for lack of instantiantiation of the above function
//// The impls are visible only in this file to save compilation time..
//
//template void boost::serialization::serialize<boost::archive::binary_iarchive>(
//    boost::archive::binary_iarchive & ar,
//    mscomplex_t & g,
//    const unsigned int file_version
//    );
//template void boost::serialization::serialize<boost::archive::binary_oarchive>(
//    boost::archive::binary_oarchive & ar,
//    mscomplex_t & g,
//    const unsigned int file_version
//    );
//
//
