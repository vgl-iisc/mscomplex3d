#include <queue>

#include <boost/foreach.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/adaptors.hpp>

#include <grid_mscomplex.h>

using namespace std;

namespace br = boost::range;
namespace ba = boost::adaptors;

namespace grid
{


void mscomplex_t::dir_connect_cps(int p, int q,int m)
{
  ASSERT(is_paired(p) != is_paired(q));
  ASSERT(abs(index(p)-index(q)) == 1);

  if(is_paired(q))
    std::swap(p,q);

  eGDIR dir = (index(p) > index(q))?(DES):(ASC);

  if(m_conn[dir][p].count(q) == 0)
    m_conn[dir][p][q] =0;

  m_conn[dir][p][q] += m;
}

void mscomplex_t::uncancel_pair(int p, int q)
{
  order_pr_by_cp_index(*this,p,q);

  ASSERT(is_canceled(p) == true && is_canceled(q) == true);
  ASSERT(index(p) == index(q)+1);
  ASSERT(pair_idx(p) == q && pair_idx(q) == p);

  m_cp_is_cancelled[p] =false;
  m_cp_is_cancelled[q] =false;

  for(int d = 0 ; d <2 ; ++d)
  {
    int ed = (d == 0)?(p):(q);

    conn_t old_conn(m_conn[d][ed].begin(),m_conn[d][ed].end());

    m_conn[d][ed].clear();

    BOOST_FOREACH(int_int_t i,old_conn)
    {
      int c = i.first;
      int m = i.second;

      if(is_paired(c) == false)
      {
        dir_connect_cps(ed,c,m);
        continue;
      }

      int r = pair_idx(c);

      if(index(ed) != index(r))
        continue;

      ASSERT(is_canceled(c) ==false && is_canceled(r) ==false);
      ASSERT(abs(index(c) - index(r)) == 1);
      ASSERT(pair_idx(r) == int(c) && pair_idx(c) ==  r);

      BOOST_FOREACH(int_int_t j,m_conn[d][r])
      {
        dir_connect_cps(ed,j.first,j.second*m);
      }
    }
  }
}


void mscomplex_t::un_simplify()
{
  BOOST_FOREACH(int_pair_t p,boost::make_iterator_range(m_canc_list)|ba::reversed)
  {
    if(is_canceled(p[0]) == false) break;
    uncancel_pair(p[0],p[1]);
  }
}


inline bool is_valid_canc_edge
(const mscomplex_t &msc, int_pair_t e,
 const std::vector<bool> &is_inc_ext,cell_fn_t thr )
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

  bool   is_epsilon_persistent = (msc.vertid(e[0]) == msc.vertid(e[1]));
  bool   is_pers_lt_t          = std::abs(msc.fn(e[0]) - msc.fn(e[1])) < thr;

  return (is_epsilon_persistent || is_pers_lt_t);
}



template<typename T>
inline void set_vec_value(std::vector<T> & vec, int i,const T& v){vec[i] = v;}

inline void make_is_inc_ext(const mscomplex_t &msc, vector<bool> &inc_on_ext)
{
  inc_on_ext.resize(msc.get_num_critpts(),false);


  auto ftor = bind(&set_vec_value<bool>,std::ref(inc_on_ext),_1,true);

  for(int i = 0 ;i < msc.get_num_critpts();++i)
  {
    if(msc.is_canceled(i)) continue;

    cellid_t c = msc.cellid(i);

    if( !msc.is_paired(i) && msc.m_rect.boundryCount(c) ==
        msc.m_ext_rect.boundryCount(c))
      continue;

    inc_on_ext[i] = true;

    br::for_each(msc.m_des_conn[i]|ba::map_keys,ftor);
    br::for_each(msc.m_asc_conn[i]|ba::map_keys,ftor);
  }
}

inline void update_is_inc_ext
(const mscomplex_t &msc, vector<bool> &is_inc_on_ext,int_pair_t pr)
{
  int p = pr[0],q = pr[1];

  cellid_t c_p = msc.cellid(p);
  cellid_t c_q = msc.cellid(q);

  if(msc.is_paired(p) || msc.m_rect.boundryCount(c_p) !=
     msc.m_ext_rect.boundryCount(c_p))
    is_inc_on_ext[q] = true;

  if(msc.is_paired(q) || msc.m_rect.boundryCount(c_q) !=
     msc.m_ext_rect.boundryCount(c_q))
    is_inc_on_ext[p] = true;
}

void mscomplex_t::simplify_pers_outcore(double f_tresh, double f_range)
{
  auto cmp = bind(&mscomplex_t::persistence_cmp,this,_2,_1);

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

      if(is_valid_canc_edge(*this,pr,is_inc_ext,f_tresh))
        pq.push(pr);
    }
  }

  while (pq.size() !=0)
  {
    int_pair_t pr = pq.top();

    pq.pop();

    if(is_valid_canc_edge(*this,pr,is_inc_ext,f_tresh) == false)
      continue;

    int &p = pr[0],&q = pr[1];

    order_pr_by_cp_index(*this,p,q);

    cancel_pair(p,q);

    BOOST_FOREACH(int_int_t i,m_des_conn[p])
    BOOST_FOREACH(int_int_t j,m_asc_conn[q])
    {
      int_pair_t npr(i.first,j.first);

      update_is_inc_ext(*this,is_inc_ext,npr);

      if(is_valid_canc_edge(*this,pr,is_inc_ext,f_tresh))
        pq.push(npr);
    }
  }
}


void mscomplex_t::invert_for_collection()
{
  for(int i = 0 ; i < get_num_critpts(); ++i)
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

    BOOST_FOREACH(int_int_t j,m_conn[dir][i])
    {
      ASSERT(is_paired(j.first) == false);

      if(m_conn[dir^1][j.first].count(i) == 0)
        m_conn[dir^1][j.first][i] = 0;

      m_conn[dir^1][j.first][i] += j.second;
    }

//      m_conn[dir][i].clear();
  }
}


inline cellid_t get_null_axes(rect_t r)
{
  cellid_t spn = r.span();

  spn[0] = (spn[0] != 0)?(0):(1);
  spn[1] = (spn[1] != 0)?(0):(1);
  spn[2] = (spn[2] != 0)?(0):(1);

  return spn;
}

inline void get_ixn(rect_t &ixn,cellid_t &ixn_dir, rect_t r1,rect_t r2,rect_t e1,rect_t e2)
{
  rect_t ir = r1.intersection(r2);
  rect_t ie = e1.intersection(e2);

  ASSERT((r1.lc() == ir.lc() && r2.uc() == ir.uc()) ||
         (r2.lc() == ir.lc() && r1.uc() == ir.uc()));

  ASSERT((e1.lc() == ie.lc() && e2.uc() == ie.uc()) ||
         (e2.lc() == ie.lc() && e1.uc() == ie.uc()));

  ixn_dir = get_null_axes(ir);

  ixn = rect_t(ie.lc()+2*ixn_dir,ie.uc()-2*ixn_dir);

  ASSERT(euclid_norm2(ixn_dir) == 1);
  ASSERT(ir.eff_dim() == gc_grid_dim-1);
  ASSERT(ixn.eff_dim() == gc_grid_dim-1);

  ixn_dir *= (r1.contains(ir.lc()+ixn_dir))?(-1):(1);

  ASSERT(r1.contains(ir.lc()-ixn_dir));
  ASSERT(r2.contains(ir.lc()+ixn_dir));
}

enum eBCROSS_TYPE {BCROSS_NOT_BND=0,BCROSS_BND_NO_CROSS,
                   BCROSS_D_DP1,BCROSS_DP1_D};

inline eBCROSS_TYPE get_bcross_type
(rect_t &ixn,cellid_t &ixn_dir,cellid_list_t &cl,int_list_t &pi,int i)
{
  if(pi[i] >= 0 && dot_product(ixn_dir,cl[i]-cl[pi[i]]) != 0)
  {
    if(ixn.contains(cl[i]))     return BCROSS_D_DP1;
    if(ixn.contains(cl[pi[i]])) return BCROSS_DP1_D;

    return BCROSS_NOT_BND;
  }
  else
  {
    return ixn.contains(cl[i])?(BCROSS_BND_NO_CROSS):(BCROSS_NOT_BND);
  }
}

template<bool KEEP_IXN_CPS>
inline void get_idx_map(int_list_t & idx_map,
                       rect_t & ixn,int_marray_t & ixn_idx,cellid_t ixn_dir,
                       cellid_list_t cl,bool_list_t &ic,int_list_t &pi,
                        int & off)
{
  int N = cl.size();
  idx_map.resize(N,-1);

  for(int i =0; i < N; ++i)
  {
    if(ic[i]) continue;

    eBCROSS_TYPE bct = get_bcross_type(ixn,ixn_dir,cl,pi,i);

    switch(bct)
    {
    case BCROSS_NOT_BND:
    {
      ASSERT(idx_map[i] == -1);
      idx_map[i] = off++;
      break;
    }
    case BCROSS_BND_NO_CROSS:
    {
      if(KEEP_IXN_CPS)
      {
        ASSERT(ixn_idx(cl[i]) == -1);
        ixn_idx(cl[i]) = off;
        off++;
      }
      else
      {
        ASSERT(ixn_idx(cl[i]) != -1);
      }
      ASSERT(idx_map[i] == -1);
      idx_map[i] = ixn_idx(cl[i]);
      break;
    }
    case BCROSS_D_DP1:
    {
      if(KEEP_IXN_CPS)
      {
        ASSERT(ixn_idx(cl[i]) == -1);
        ixn_idx(cl[i]) = off;
        off += 2;
      }
      else
      {
        ASSERT(ixn_idx(cl[i]) != -1);
      }

      ASSERT(idx_map[i] == -1);
      ASSERT(idx_map[pi[i]] == -1);

      idx_map[i]     = ixn_idx(cl[i]);
      idx_map[pi[i]] = ixn_idx(cl[i])+1;
      break;
    }
    case BCROSS_DP1_D:break;
    default:
      ASSERT(false&&"incorrect classification");
    };
  }
}

inline void copy_cp_info(mscomplex_t &msc,int_list_t & idx_map,
cellid_list_t  &cl,cellid_list_t  &vl,int_list_t &pi,int8_list_t &ci,cell_fn_list_t &fn)
{
  int N = idx_map.size();

  for(int i = 0 ; i < N; ++i)
  {
    int j = idx_map[i];

    if(j >=0)
    {
      if(pi[i] != -1)
      {
        ASSERT(msc.m_cp_pair_idx[j] == -1 || msc.m_cp_pair_idx[j] == idx_map[pi[i]]);
        msc.m_cp_pair_idx[j] = idx_map[pi[i]];
      }

      msc.m_cp_cellid[j]       = cl[i];
      msc.m_cp_vertid[j]       = vl[i];
      msc.m_cp_index[j]        = ci[i];
      msc.m_cp_is_cancelled[j] = false;
      msc.m_cp_fn[j]           = fn[i];

    }
  }
}

template<bool KEEP_IXN_EDGES>
inline void copy_adj_info(mscomplex_t & msc,int_list_t & idx_map,
                          int_list_t& nconn,int_int_list_t &adj,
                          rect_t &ixn,cellid_list_t & cl)
{
  int N1 = nconn.size()/2;

  int_int_list_t::iterator a,b,c = adj.begin();

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
      if(KEEP_IXN_EDGES || !(in_ixn && ixn.contains(cl[a->first])))
      {
        ASSERT(is_in_range(idx_map[a->first],0,msc.get_num_critpts()));
        msc.m_des_conn[j].insert(make_pair(idx_map[a->first],a->second));
      }

    for(;b != c; ++b)
      if(KEEP_IXN_EDGES || !(in_ixn && ixn.contains(cl[b->first])))
      {
        ASSERT(is_in_range(idx_map[b->first],0,msc.get_num_critpts()));
        msc.m_asc_conn[j].insert(make_pair(idx_map[b->first],b->second));
      }
  }
}

inline bool check_all_cps_in
(const mscomplex_t &msc,const cellid_list_t &cl,const bool_list_t &ic,
 const int_list_t &pi, rect_t r1,rect_t r2)
{
  std::multiset<cellid_t> cset(msc.m_cp_cellid.begin(),msc.m_cp_cellid.end());

  int N = cl.size();

  for( int i = 0 ; i < N; ++i)
  {
    int expect_cct = (ic[i])?(0):(1);
    int actual_cct = cset.count(cl[i]);

    ASSERT(expect_cct == actual_cct);
  }

  return true;
}

inline bool check_boundry_consistency
(cellid_list_t &cl1,cellid_list_t &cl2,
 int_list_t &pi1, int_list_t &pi2,
 bool_list_t &ic1,bool_list_t &ic2,
 rect_t ixn,cellid_t ixn_dir)
{
  int_marray_t ixn_idx;
  ixn_idx.resize(ixn.span()+1);
  ixn_idx.reindex(ixn.lc());

  memset(ixn_idx.data(),-1,(ixn.pt_end()-ixn.pt_begin())*sizeof(int));

  int N1  = cl1.size();
  int N2  = cl2.size();

  for( int i1 = 0 ; i1 < N1; ++i1)
  {
    if(ic1[i1]) continue;

    eBCROSS_TYPE bct = get_bcross_type(ixn,ixn_dir,cl1,pi1,i1);

    if(bct == BCROSS_D_DP1 ||bct == BCROSS_BND_NO_CROSS)
    {
      ASSERT(ixn_idx(cl1[i1]) == -1);
      ixn_idx(cl1[i1]) = i1;
    }
  }

  for( int i2 = 0 ; i2 < N2; ++i2)
  {
    if(ic2[i2]) continue;

    eBCROSS_TYPE bct = get_bcross_type(ixn,ixn_dir,cl2,pi2,i2);

    if(bct == BCROSS_D_DP1)
    {
      int i1 = ixn_idx(cl2[i2]);

      ASSERT(i1>=0);
      ASSERT(cl1[i1] == cl2[i2]);
      ASSERT(pi1[i1] >=0);
      ASSERT(cl1[pi1[i1]] == cl2[pi2[i2]]);

      ixn_idx(cl2[i2]) = -1;
    }

    if( bct == BCROSS_BND_NO_CROSS)
    {
      int i1 = ixn_idx(cl2[i2]);

      ASSERT(i1>=0);
      ASSERT(cl1[i1] == cl2[i2]);

      if(pi2[i2] >= 0 )
      {
        ASSERT(pi1[i1] >=0);
        ASSERT(cl1[pi1[i1]] == cl2[pi2[i2]]);
      }

      ixn_idx(cl2[i2]) = -1;
    }

  }

  for( int i1 = 0 ; i1 < N1; ++i1)
  {
    if(ic1[i1]) continue;

    eBCROSS_TYPE bct = get_bcross_type(ixn,ixn_dir,cl1,pi1,i1);

    if(bct == BCROSS_D_DP1 ||bct == BCROSS_BND_NO_CROSS)
    {
      ASSERT(ixn_idx(cl1[i1]) == -1);
    }
  }
  return true;
}

inline void copy_from_streams
(mscomplex_t &msc, std::istream &is1,std::istream &is2, rect_t &ixn,int_marray_t &ixn_idx,cellid_t &ixn_dir)
{
  rect_t         r1,r2,e1,e2,d1,d2;

  utl::bin_read(is1,r1);utl::bin_read(is2,r2);
  utl::bin_read(is1,e1);utl::bin_read(is2,e2);
  utl::bin_read(is1,d1);utl::bin_read(is2,d2);

  get_ixn(ixn,ixn_dir,r1,r2,e1,e2);

  ixn_idx.resize(ixn.span()+1);
  ixn_idx.reindex(ixn.lc());

  memset(ixn_idx.data(),-1,(ixn.pt_end()-ixn.pt_begin())*sizeof(int));

  int_list_t    idx_map1,idx_map2;

  cellid_list_t  cl1,cl2;
  cellid_list_t  vl1,vl2;
  int_list_t     pi1,pi2;
  int8_list_t    ci1,ci2;
  bool_list_t    ic1,ic2;
  cell_fn_list_t fn1,fn2;

  int_list_t     nconn1,nconn2;
  int_int_list_t adj1,adj2;
  int            NC1,NC2;

  utl::bin_read_vec(is1,cl1);      utl::bin_read_vec(is2,cl2);
  utl::bin_read_vec(is1,vl1);      utl::bin_read_vec(is2,vl2);
  utl::bin_read_vec(is1,pi1);      utl::bin_read_vec(is2,pi2);
  utl::bin_read_vec(is1,ci1);      utl::bin_read_vec(is2,ci2);
  utl::bin_read_vec(is1,ic1);      utl::bin_read_vec(is2,ic2);
  utl::bin_read_vec(is1,fn1);      utl::bin_read_vec(is2,fn2);

  utl::bin_read_vec(is1,nconn1);   utl::bin_read_vec(is2,nconn2);
  utl::bin_read_vec(is1,adj1);     utl::bin_read_vec(is2,adj2);

  int offset = 0;

  get_idx_map<true>(idx_map1,ixn,ixn_idx,ixn_dir,cl1,ic1,pi1,offset);
  get_idx_map<false>(idx_map2,ixn,ixn_idx,ixn_dir,cl2,ic2,pi2,offset);

  msc.resize(offset);

  ASSERT(check_boundry_consistency(cl1,cl2,pi1,pi2,ic1,ic2,ixn,ixn_dir));
  ASSERT(check_boundry_consistency(cl2,cl1,pi2,pi1,ic2,ic1,ixn,ixn_dir));

  copy_cp_info(msc,idx_map1,cl1,vl1,pi1,ci1,fn1);
  copy_cp_info(msc,idx_map2,cl2,vl2,pi2,ci2,fn2);

  ASSERT(check_all_cps_in(msc,cl1,ic1,pi1,r1,r2));
  ASSERT(check_all_cps_in(msc,cl2,ic2,pi2,r2,r1));

  copy_adj_info<true>(msc,idx_map1,nconn1,adj1,ixn,cl1);
  copy_adj_info<false>(msc,idx_map2,nconn2,adj2,ixn,cl2);
}

int mscomplex_t::load_merge(const string &f1, const string &f2)
{
  rect_t        ixn;
  int_marray_t  ixn_idx;
  cellid_t      ixn_dir;

  std::fstream is1(f1.c_str(),std::ios::in|std::ios::binary);
  std::fstream is2(f2.c_str(),std::ios::in|std::ios::binary);

  copy_from_streams(*this,is1,is2,ixn,ixn_idx,ixn_dir);

  int num_c = 0;

  for(rect_t::pt_iterator b= ixn.pt_begin(),e=ixn.pt_end(); b != e; ++b)
  {
    cellid_t c = *b;

    int p = ixn_idx(c);

    if(p == -1 ) continue;

    if(!is_paired(p)) continue;

    int q = pair_idx(p);

    ASSERT(pair_idx(pair_idx(p)) == p);

    if(dot_product(ixn_dir,cellid(p) - cellid(q)) == 0) continue;

    cancel_pair(p,q);

    num_c ++;
  }
  return num_c;
}

inline void fill_ixn_idx(int_marray_t &ixn_idx,rect_t ixn,const cellid_list_t &l)
{
  int N = l.size();

  for(int i = 0 ; i < N; ++i)
    if(ixn.contains(l[i]))
      ixn_idx(l[i]) = i;
}

inline void make_rev_map(const int_list_t & idx_map, int N,int_list_t &ridx_map)
{
  ridx_map.resize(N,-1);

  int n = idx_map.size();

  for(int i = 0; i < n; ++i)
    if(idx_map[i] >= 0)
      ridx_map[idx_map[i]] = i;
}

inline int_int_t adj_converter
  (const int_int_t &p,int_list_t &idx_map,int_list_t &ridx_map)
{
  int j = p.first;

  if(ridx_map[j] == -1)
  {
    ridx_map[j] = idx_map.size();
    idx_map.push_back(j);
  }

  return make_pair(ridx_map[j],p.second);
}

inline bool is_surv(int_int_t i, const mscomplex_t &msc, const int_list_t &idx_map)
{
  ASSERT(idx_map[i.first] >= 0);
  return !(msc.is_paired(idx_map[i.first]));
}

inline void copy_into_new_adj
( int_list_t &idx_map,int_list_t &ridx_map,
  int_int_list_t &adj,int_list_t &nconn,
  int_int_list_t &new_adj,
 const mscomplex_t &msc)
{
  int N = idx_map.size();


  auto ac_ftr = bind(adj_converter,_1,std::ref(idx_map),std::ref(ridx_map));
  auto surv_ftor = bind(is_surv,_1,std::cref(msc),std::cref(idx_map));

  int_int_list_t::iterator a,b,c = adj.begin();

  for(int i = 0 ; i < N; ++i)
  {
    a = c;
    b = a + (nconn[2*i]);
    c = b + (nconn[2*i+1]);

    int j = idx_map[i];

    if(j >= 0)
    {
      if( msc.is_paired(j))
      {
        br::transform(msc.m_des_conn[j],back_inserter(new_adj),ac_ftr);
        br::transform(msc.m_asc_conn[j],back_inserter(new_adj),ac_ftr);

        nconn[2*i]   = msc.m_des_conn[j].size();
        nconn[2*i+1] = msc.m_asc_conn[j].size();
      }
      else
      {
        int s_a = new_adj.size();
        br::copy(make_pair(a,b)|ba::filtered(surv_ftor), back_inserter(new_adj));
        int s_b = new_adj.size();
        br::copy(make_pair(b,c)|ba::filtered(surv_ftor), back_inserter(new_adj));
        int s_c = new_adj.size();

        nconn[2*i]   = s_b - s_a;
        nconn[2*i+1] = s_c - s_b;
      }
    }
    else
    {
      copy(a,c,back_inserter(new_adj));
    }
  }
}

inline void copy_new_cp_info
( cellid_list_t  &cl,cellid_list_t  &vl,int_list_t &pi,int8_list_t &ci,
  bool_list_t    &ic,cell_fn_list_t &fn,
  const mscomplex_t &msc, const int_list_t &idx_map,const int_list_t &ridx_map)
{
  int n = cl.size();
  int N = idx_map.size();

  ASSERT(n <= N);

  cl.resize(N);
  vl.resize(N);
  pi.resize(N);
  ci.resize(N);
  ic.resize(N,false);
  fn.resize(N);

  for(int i = 0 ;i < n; ++i)
  {
    int j = idx_map[i];

    if(j >=0)
    {
      ASSERT(!msc.is_paired(j) || ridx_map[msc.pair_idx(j)] != -1);
      ASSERT(pi[i] == -1 || (msc.is_paired(j) && pi[i] == ridx_map[msc.pair_idx(j)]));

      pi[i] = (msc.is_paired(j))?(ridx_map[msc.pair_idx(j)]):(-1);

      ASSERT(!msc.is_paired(j) ||is_in_range(pi[i],0,N));
    }
  }

  for(int i = n ;i < N; ++i)
  {
    int j = idx_map[i];

    ASSERT(j >=0);
    ASSERT(!msc.is_paired(j) || ridx_map[msc.pair_idx(j)] != -1);

    cl[i] = msc.cellid(j);
    vl[i] = msc.vertid(j);
    pi[i] = (msc.is_paired(j))?(ridx_map[msc.pair_idx(j)]):(-1);
    ASSERT(!msc.is_paired(j) ||is_in_range(pi[i],0,N));
    ci[i] = msc.index(j);
//      ic[i] = msc.is_canceled(j);
    ASSERT(msc.is_canceled(j) == false);
    fn[i] = msc.fn(j);
  }
}

inline void update_maps_for_new_pairs
( const mscomplex_t &msc,int_list_t &idx_map, int_list_t &ridx_map)
{
  int N = idx_map.size();

  for(int i = 0 ; i < N; ++i)
  {
    int j = idx_map[i];

    if(j >= 0)
    {
      if(msc.is_paired(j))
      {
        if(ridx_map[msc.pair_idx(j)] == -1)
        {
          ridx_map[msc.pair_idx(j)] = idx_map.size();
          idx_map.push_back(msc.pair_idx(j));
        }
      }
    }
  }
}


template<bool KEEP_IXN_CPS>
void copy_into_stream
(const mscomplex_t &msc, std::iostream &io,rect_t r,rect_t e,rect_t d,
 int_marray_t &ixn_idx,rect_t ixn,cellid_t ixn_dir,int &off)
{
  int_list_t    idx_map,ridx_map;

  cellid_list_t  cl;
  cellid_list_t  vl;
  int_list_t     pi;
  int8_list_t    ci;
  bool_list_t    ic;
  cell_fn_list_t fn;

  int_list_t      nconn;
  int_int_list_t  adj;
  int_pair_list_t cancl;

  utl::bin_read_vec(io,cl);
  utl::bin_read_vec(io,vl);
  utl::bin_read_vec(io,pi);
  utl::bin_read_vec(io,ci);
  utl::bin_read_vec(io,ic);
  utl::bin_read_vec(io,fn);

  utl::bin_read_vec(io,nconn);
  utl::bin_read_vec(io,adj);

  get_idx_map<KEEP_IXN_CPS>(idx_map,ixn,ixn_idx,ixn_dir,cl,ic,pi,off);

  make_rev_map(idx_map,msc.get_num_critpts(),ridx_map);

  int_int_list_t new_adj;

  update_maps_for_new_pairs(msc,idx_map,ridx_map);
  nconn.resize(idx_map.size()*2,0);
  copy_into_new_adj(idx_map,ridx_map,adj,nconn,new_adj,msc);

  nconn.resize(idx_map.size()*2,0);
  copy_new_cp_info(cl,vl,pi,ci,ic,fn,msc,idx_map,ridx_map);

  idx_map.clear();ridx_map.clear();adj.clear();

  utl::bin_read_vec(io,cancl);

  ASSERT(int(nconn.size()) == cl.size()*2);

  io.seekp(0,ios::beg);

  utl::bin_write(io,r);
  utl::bin_write(io,e);
  utl::bin_write(io,d);

  utl::bin_write_vec(io,cl);
  utl::bin_write_vec(io,vl);
  utl::bin_write_vec(io,pi);
  utl::bin_write_vec(io,ci);
  utl::bin_write_vec(io,ic);
  utl::bin_write_vec(io,fn);

  utl::bin_write_vec(io,nconn);
  utl::bin_write_vec(io,new_adj);

  utl::bin_write_vec(io,cancl);
}

void mscomplex_t::unmerge_save(const string &f1, const string &f2)
{
  std::fstream is1(f1.c_str(),std::ios::in|std::ios::out|std::ios::binary);
  std::fstream is2(f2.c_str(),std::ios::in|std::ios::out|std::ios::binary);

  rect_t        ixn;
  int_marray_t  ixn_idx;
  cellid_t      ixn_dir;

  rect_t         r1,r2,e1,e2,d1,d2;

  utl::bin_read(is1,r1);utl::bin_read(is2,r2);
  utl::bin_read(is1,e1);utl::bin_read(is2,e2);
  utl::bin_read(is1,d1);utl::bin_read(is2,d2);

  get_ixn(ixn,ixn_dir,r1,r2,e1,e2);

  ixn_idx.resize(ixn.span()+1);
  ixn_idx.reindex(ixn.lc());

  memset(ixn_idx.data(),-1,(ixn.pt_end()-ixn.pt_begin())*sizeof(int));

  fill_ixn_idx(ixn_idx,ixn,m_cp_cellid);

  for(rect_t::pt_riterator b= ixn.pt_rbegin(),e=ixn.pt_rend(); b != e; ++b)
  {
    cellid_t c = *b;

    int p = ixn_idx(c);

    if(p == -1 ) continue;

    if(!is_paired(p)) continue;

    int q = pair_idx(p);

    ASSERT(pair_idx(pair_idx(p)) == p);

    if(is_canceled(p) == false) continue;

    if(dot_product(ixn_dir,cellid(p) - cellid(q)) == 0) continue;

    uncancel_pair(p,q);
  }

  memset(ixn_idx.data(),-1,(ixn.pt_end()-ixn.pt_begin())*sizeof(int));

  int off = 0;

  copy_into_stream<true>(*this,is1,r1,e1,d1,ixn_idx,ixn,ixn_dir,off);
  copy_into_stream<false>(*this,is2,r2,e2,d2,ixn_idx,ixn,ixn_dir,off);
}


}
