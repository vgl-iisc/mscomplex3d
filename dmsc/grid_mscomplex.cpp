#include <cmath>
#include <queue>
#include <iostream>
#include <fstream>
#include <limits>
#include <cstdlib>

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

  mscomplex_t::mscomplex_t(rect_t r,rect_t e)
    :m_rect(r),m_ext_rect(e),m_des_conn(m_conn[0]),m_asc_conn(m_conn[1])
  {
  }

  mscomplex_t::~mscomplex_t()
  {
    clear();
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

  void mscomplex_t::merge_up
      (const mscomplex_t& msc1,
       const mscomplex_t& msc2,
       const rect_t& bnd)
  {
    typedef conn_t::const_iterator cconn_it_t;

    try
    {
      ASSERT(bnd.eff_dim() == gc_grid_dim-1);
      ASSERT(m_ext_rect.intersection(bnd) == bnd);

      const mscomplex_t *msc_arr[] = {&msc1,&msc2};

      for (int m = 0 ; m < 2 ; ++m)
      {
        const mscomplex_t &msc = *msc_arr[m];

        for(int j = 0; j < msc.get_num_critpts();++j)
        {
          try
          {
            if(msc.is_canceled(j))
              continue;

            ASSERT((bnd.contains(msc.cellid(j)) == false) || (msc_arr[m^1]->m_id_cp_map.count(msc.cellid(j)) == 1));
            ASSERT((bnd.contains(msc.cellid(j)) == true)  || (msc_arr[m^1]->m_id_cp_map.count(msc.cellid(j)) == 0));

            if((m == 1) && (bnd.contains(msc.cellid(j))))
              continue;

            add_critpt(msc.cellid(j),msc.index(j),msc.fn(j),msc.vertid(j));
          }
          catch(assertion_error e)
          {
            e.push(_FFL);
            e.push(SVAR(msc.cp_info(j)));
            e.push(SVAR(msc.m_rect));
            e.push(SVAR(msc.m_ext_rect));
            throw;
          }
        }
      }

      for (int m = 0 ; m < 2 ; ++m)
      {
        const mscomplex_t &msc = *msc_arr[m];

        for(int i = 0; i < msc.get_num_critpts();++i)
        {
          if(msc.is_canceled(i))
            continue;

          bool is_cp_in_bnd = bnd.contains(msc.cellid(i));

          for(cconn_it_t j = msc.m_des_conn[i].begin();j != msc.m_des_conn[i].end();++j)
          {
            if(msc.is_canceled(*j))
              continue;

            if((m == 1) && is_cp_in_bnd && bnd.contains(msc.cellid(*j)))
              continue;

            connect_cps(msc.cellid(i),msc.cellid(*j));
          }

          if(!msc.is_paired(i))
            continue;

          pair_cps(msc.cellid(i),msc.cellid(msc.pair_idx(i)));
        }
      }

      rect_t ixn = m_rect.intersection(bnd);

      for(cellid_t c = ixn.lower_corner() ; c[2] <= ixn[2][1]; ++c[2])
      {
        for(c[1] = ixn[1][0] ; c[1] <= ixn[1][1]; ++c[1])
        {
          for(c[0] = ixn[0][0] ; c[0] <= ixn[0][1]; ++c[0])
          {
            if (m_id_cp_map.count(c) == 0)
              continue;

            int p = m_id_cp_map[c];

            if(!is_paired(p))
              continue;

            ASSERT(pair_idx(pair_idx(p)) == p);

            int q = pair_idx(p);

            if(bnd.contains(cellid(q)))
              continue;

            cancel_pair(p,q);
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

  struct persistence_comparator_t
  {
    mscomplex_t *m_msc;

    persistence_comparator_t(mscomplex_t *m):m_msc(m){}

    bool operator()(const int_pair_t & p0, const int_pair_t &p1)
    {
      return cmp_lt(p1,p0);
    }

    bool cmp_lt(int_pair_t p0, int_pair_t p1)
    {
      order_pr_by_cp_index(*m_msc,p0[0],p0[1]);
      order_pr_by_cp_index(*m_msc,p1[0],p1[1]);

      cellid_t v00 = m_msc->vertid(p0[0]);
      cellid_t v01 = m_msc->vertid(p0[1]);
      cellid_t v10 = m_msc->vertid(p1[0]);
      cellid_t v11 = m_msc->vertid(p1[1]);

      cellid_t c00 = m_msc->cellid(p0[0]);
      cellid_t c01 = m_msc->cellid(p0[1]);
      cellid_t c10 = m_msc->cellid(p1[0]);
      cellid_t c11 = m_msc->cellid(p1[1]);

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

      cell_fn_t f00 = m_msc->fn(p0[0]);
      cell_fn_t f01 = m_msc->fn(p0[1]);
      cell_fn_t f10 = m_msc->fn(p1[0]);
      cell_fn_t f11 = m_msc->fn(p1[1]);

      cell_fn_t d1 = std::abs(f01-f00);
      cell_fn_t d2 = std::abs(f11-f10);

      if(d1 != d2)
        return d1 < d2;

      if(c00 != c10)
        return c00 < c10;

      return c01 < c11;
    }
  };

  bool is_valid_canc_edge(mscomplex_t *msc,int_pair_t e )
  {
    order_pr_by_cp_index(*msc,e[0],e[1]);

    if(msc->is_canceled(e[0])||msc->is_canceled(e[1]))
      return false;

    if(msc->m_rect.isOnBoundry(msc->cellid(e[0])) !=
       msc->m_rect.isOnBoundry(msc->cellid(e[1])))
      return false;

    ASSERT(msc->m_des_conn[e[0]].count(e[1]) == msc->m_asc_conn[e[1]].count(e[0]));

    if(msc->m_des_conn[e[0]].count(e[1]) != 1)
      return false;

    return true;
  }

  bool is_epsilon_persistent(mscomplex_t *msc,int_pair_t e )
  {
    return (msc->vertid(e[0]) == msc->vertid(e[1]));
  }

  void mscomplex_t::simplify(int_pair_list_t & canc_pairs_list,
                               double simplification_treshold)
  {
    typedef std::priority_queue
        <int_pair_t,int_pair_list_t,persistence_comparator_t>
        canc_pair_priq_t;

    persistence_comparator_t comp(this);

    canc_pair_priq_t  canc_pair_priq(comp);

    cell_fn_t max_val = std::numeric_limits<cell_fn_t>::min();
    cell_fn_t min_val = std::numeric_limits<cell_fn_t>::max();

    for(uint i = 0 ;i < get_num_critpts();++i)
    {
      max_val = std::max(max_val,fn(i));
      min_val = std::min(min_val,fn(i));

      for(conn_iter_t j = m_des_conn[i].begin();j != m_des_conn[i].end() ;++j)
      {
        if(is_valid_canc_edge(this,int_pair_t(i,*j)))
          canc_pair_priq.push(int_pair_t(i,*j));
      }
    }

    double max_persistence = max_val - min_val;

    uint num_cancellations = 0;

    uint num_cancellations_eps = 0;

    while (canc_pair_priq.size() !=0)
    {
      int_pair_t pr = canc_pair_priq.top();

      canc_pair_priq.pop();

      double persistence = std::abs(fn(pr[0])-fn(pr[1]))/max_persistence;

      if(is_valid_canc_edge(this,pr) == false)
        continue;

      if(is_epsilon_persistent(this,pr) == false)
      {
        if(persistence >= simplification_treshold)
          break;
      }
      else
      {
        num_cancellations_eps++;
      }

//      std::cout
//             <<   "no = "<<num_cancellations<<" "
//             << "pers = "<<persistence<<" "
//             <<"index = ("<<(int)m_cps[pr[0]]->index<<","<<(int)m_cps[pr[1]]->index<<") "
//             << "edge = "<<pr<<" "
//             << "edge = "<<edge_to_string(this,pr)<<" "
//             <<std::endl;

      int p = pr[0],q = pr[1];

      order_pr_by_cp_index(*this,p,q);

      pair_cps(p,q);

      cancel_pair(p,q);

      num_cancellations++;

      canc_pairs_list.push_back(pr);

      for(conn_iter_t i = m_des_conn[p].begin();i != m_des_conn[p].end();++i)
        for(conn_iter_t j = m_asc_conn[q].begin();j != m_asc_conn[q].end();++j)
            canc_pair_priq.push(int_pair_t(*i,*j));

    }
    cout<<"num_cancellations    ::"<<(num_cancellations)<<endl;
    cout<<"num_cancellations_eps::"<<(num_cancellations_eps)<<endl;
  }

  void mscomplex_t::un_simplify(const int_pair_list_t &canc_pairs_list)
  {
    typedef int_pair_list_t::const_reverse_iterator revit_t;

    for(revit_t it = canc_pairs_list.rbegin();it != canc_pairs_list.rend() ; ++it)
      uncancel_pair((*it)[0],(*it)[1]);
  }

  void mscomplex_t::simplify_un_simplify(double simplification_treshold)
  {
    int_pair_list_t canc_pairs_list;

    simplify(canc_pairs_list,simplification_treshold);

    un_simplify(canc_pairs_list);
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

  void mscomplex_t::write_graph(const std::string &fn) const
  {
    std::fstream os(fn.c_str(),std::ios::out);

    ensure(os.is_open(),"failed to open file");

    write_graph(os);

    os.close();
  }

  cp_producer_t::cp_producer_t(mscomplex_const_ptr_t msc, cp_filter_t cf)
  {
    boost::mutex::scoped_lock scoped_lock(m_mutex);

    m_msc       = msc;
    m_cp_filter = cf;
    m_ni        = 0;


    while(m_ni < m_msc->get_num_critpts() && (m_cp_filter(m_msc,m_ni) == false))
      ++m_ni;
  }

  bool cp_producer_t::next(int & i)
  {
    boost::mutex::scoped_lock scoped_lock(m_mutex);

    if(!is_in_range(m_ni,0,m_msc->get_num_critpts()))
      return false;

    i = m_ni++;

    while(m_ni < m_msc->get_num_critpts() && (m_cp_filter(m_msc,m_ni) == false))
      ++m_ni;

    return true;
  }

  int cp_producer_t::count() const
  {
    if (m_cp_filter == pass_filter)
      return m_msc->get_num_critpts();

    int ct = 0;

    for(int i = 0 ; i < m_msc->get_num_critpts(); ++i)
    {
      if (m_cp_filter(m_msc,i))
        ++ct;
    }

    return ct;
  }


  bool cp_producer_t::pass_filter(mscomplex_const_ptr_t msc, int i)
  {
    return true;
  }
  bool cp_producer_t::extrema_filter(mscomplex_const_ptr_t msc, int i)
  {
    return msc->index(i) == 3 || msc->index(i) == 0;
  }
  bool cp_producer_t::twosaddle_filter(mscomplex_const_ptr_t msc, int i)
  {
    return msc->index(i) == 2;
  }
  bool cp_producer_t::saddle_filter(mscomplex_const_ptr_t msc, int i)
  {
    return msc->index(i) == 2 || msc->index(i) == 1;
  }
  bool cp_producer_t::unpaired_cp_filter(mscomplex_const_ptr_t msc, int i)
  {
    return !msc->is_paired(i);
  }
  bool cp_producer_t::unpaired_saddle_filter(mscomplex_const_ptr_t msc, int i)
  {
    return (!msc->is_paired(i)) && (msc->index(i) == 2 || msc->index(i) == 1);
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
