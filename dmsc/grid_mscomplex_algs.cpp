#include <queue>

#include <boost/foreach.hpp>
#include <boost/range/adaptors.hpp>

#include <grid_mscomplex.h>

using namespace std;

namespace br = boost::range;
namespace ba = boost::adaptors;

namespace grid
{


void mscomplex_t::pair_cps(int p, int q)
{
  m_cp_pair_idx[p] = q;
  m_cp_pair_idx[q] = p;
}

void mscomplex_t::cancel_pair ( int p, int q)
{
  order_pr_by_cp_index(*this,p,q);

  ASSERT(index(p) == index(q)+1);
  ASSERT(pair_idx(p) == q);
  ASSERT(pair_idx(q) == p);
  ASSERT(is_canceled(p) == false);
  ASSERT(is_canceled(q) == false);
  ASSERT(m_des_conn[p].count(q) == 1);
  ASSERT(m_asc_conn[q].count(p) == 1);
  ASSERT(m_des_conn[p][q] == 1);
  ASSERT(m_asc_conn[q][p] == 1);

  m_des_conn[p].erase(q);
  m_asc_conn[q].erase(p);

  // cps in lower of u except l
  BOOST_FOREACH(int_int_t i,m_des_conn[p])
      BOOST_FOREACH(int_int_t j,m_asc_conn[q])
  {
    int u = i.first;
    int v = j.first;
    int m = i.second*j.second;

    ASSERT(is_canceled(u) == false);
    ASSERT(is_canceled(v) == false);

    connect_cps(u,v,m);
  }

  BOOST_FOREACH(conn_t::value_type pr,m_des_conn[p]) m_asc_conn[pr.first].erase(p);
  BOOST_FOREACH(conn_t::value_type pr,m_asc_conn[p]) m_des_conn[pr.first].erase(p);
  BOOST_FOREACH(conn_t::value_type pr,m_des_conn[q]) m_asc_conn[pr.first].erase(q);
  BOOST_FOREACH(conn_t::value_type pr,m_asc_conn[q]) m_des_conn[pr.first].erase(q);

  m_cp_is_cancelled[p] =true;
  m_cp_is_cancelled[q] =true;

  m_asc_conn[p].clear();
  m_des_conn[q].clear();
}


inline bool is_valid_canc_edge
(const mscomplex_t &msc,const std::vector<bool> &is_inc_ext, int_pair_t e )
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

  ASSERT(msc.m_des_conn[e[0]].count(e[1]) == 1);
  ASSERT(msc.m_asc_conn[e[1]].count(e[0]) == 1);
  ASSERT(msc.m_des_conn[e[0]][e[1]] == msc.m_asc_conn[e[1]][e[0]]);

  if(msc.m_des_conn[e[0]][e[1]] != 1)
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

inline bool fast_persistence_lt
(const mscomplex_t &msc, int_pair_t p0, int_pair_t p1)
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

inline bool is_in_treshold(const mscomplex_t & msc,int_pair_t e,cell_fn_t t)
{
  return (is_epsilon_persistent(msc,e) || get_persistence(msc,e) < t);
}

template<typename T>
inline void set_vec_value(std::vector<T> & vec, int i,const T& v){vec[i] = v;}


inline void make_is_inc_ext(const mscomplex_t &msc, vector<bool> &inc_on_ext)
{
  inc_on_ext.resize(msc.get_num_critpts(),false);

  using boost::ref;
  using boost::cref;

  BOOST_AUTO(ftor,bind(&set_vec_value<bool>,ref(inc_on_ext),_1,true));

  for(int i = 0 ;i < msc.get_num_critpts();++i)
  {
    if(msc.is_canceled(i)) continue;

    cellid_t c = msc.cellid(i);

    if( !msc.is_paired(i) && msc.m_rect.boundryCount(c) == msc.m_ext_rect.boundryCount(c))
      continue;

    inc_on_ext[i] = true;

    br::for_each(msc.m_des_conn[i]|ba::map_keys,ftor);
    br::for_each(msc.m_asc_conn[i]|ba::map_keys,ftor);
  }
}

inline void update_is_inc_ext(const mscomplex_t &msc, vector<bool> &is_inc_on_ext,int_pair_t pr)
{
  int p = pr[0],q = pr[1];

  cellid_t c_p = msc.cellid(p);
  cellid_t c_q = msc.cellid(q);

  if(msc.is_paired(p) || msc.m_rect.boundryCount(c_p) != msc.m_ext_rect.boundryCount(c_p))
    is_inc_on_ext[q] = true;

  if(msc.is_paired(q) || msc.m_rect.boundryCount(c_q) != msc.m_ext_rect.boundryCount(c_q))
    is_inc_on_ext[p] = true;
}

void mscomplex_t::simplify(double f_tresh, double f_range)
{
  BOOST_AUTO(cmp,bind(persistence_lt,boost::cref(*this),_2,_1));

  priority_queue<int_pair_t,int_pair_list_t,typeof(cmp)> pq(cmp);

  if(f_range <= 0)
    f_range = *br::max_element(m_cp_fn) - *br::min_element(m_cp_fn);

  f_tresh *= f_range;

  vector<bool> is_inc_ext;

  make_is_inc_ext(*this,is_inc_ext);

  for(int i = 0 ;i < get_num_critpts();++i)
  {
    BOOST_FOREACH(int_int_t j,m_des_conn[i])
    {
      int_pair_t pr(i,j.first);

      if(is_valid_canc_edge(*this,is_inc_ext,pr) && is_in_treshold(*this,pr,f_tresh))
        pq.push(pr);
    }
  }

  while (pq.size() !=0)
  {
    int_pair_t pr = pq.top();

    pq.pop();

    if(is_valid_canc_edge(*this,is_inc_ext,pr) == false)
      continue;

    int &p = pr[0],&q = pr[1];

    order_pr_by_cp_index(*this,p,q);

    pair_cps(p,q);

    cancel_pair(p,q);

    m_canc_list.push_back(pr);

    BOOST_FOREACH(int_int_t i,m_des_conn[p])
        BOOST_FOREACH(int_int_t j,m_asc_conn[q])
    {
      int_pair_t npr(i.first,j.first);

      update_is_inc_ext(*this,is_inc_ext,npr);

      if(is_valid_canc_edge(*this,is_inc_ext,npr) && is_in_treshold(*this,npr,f_tresh))
        pq.push(npr);

    }
  }
}
}
