#include <queue>

#include <boost/foreach.hpp>
#include <boost/range/adaptors.hpp>

#include <grid_dataset.h>
#include <grid_mscomplex.h>

using namespace std;

namespace br = boost::range;
namespace ba = boost::adaptors;

namespace grid
{

/*===========================================================================*/

void mscomplex_t::cancel_pair ( int p, int q)
{
  order_pr_by_cp_index(*this,p,q);

  ENSURE(m_canc_pos == m_canc_list.size(),
         "Cannot cancel pair !! Ms complex resolution is not coarsest.");
  ENSURE(index(p) == index(q)+1,
         "indices do not differ by 1");
  ENSURE(m_cp_pair_idx[p] == -1 && m_cp_pair_idx[q] == -1,
         "p/q has already been paired");
  ENSURE(m_des_conn[p].count(q) == 1 && m_des_conn[p].count(q),
         "p is not connected to q");
  ENSURE(m_des_conn[p][q] == 1 && m_asc_conn[q][p] == 1,
         "p and q are multiply connected");

//  m_cp_cancno[p] = m_canc_list.size();
//  m_cp_cancno[q] = m_canc_list.size();
  m_canc_list.push_back(int_pair_t(p,q));

  cancel_pair();
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::cancel_pair ()
{
  ENSURE(is_in_range(m_canc_pos,0,m_canc_list.size()),
         "invalid cancellation position");

  int p = m_canc_list[m_canc_pos][0];
  int q = m_canc_list[m_canc_pos][1];

  m_canc_pos++;

  ASSERT(index(p) == index(q)+1);
  ASSERT(m_cp_pair_idx[p] == -1);
  ASSERT(m_cp_pair_idx[q] == -1);
  ASSERT(m_des_conn[p].count(q) == 1);
  ASSERT(m_asc_conn[q].count(p) == 1);
  ASSERT(m_des_conn[p][q] == 1);
  ASSERT(m_asc_conn[q][p] == 1);

  m_cp_pair_idx[p] = q;
  m_cp_pair_idx[q] = p;

  m_des_conn[p].erase(q);
  m_asc_conn[q].erase(p);

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

  BOOST_FOREACH(int_int_t pr,m_des_conn[p]) m_asc_conn[pr.first].erase(p);
  BOOST_FOREACH(int_int_t pr,m_asc_conn[p]) m_des_conn[pr.first].erase(p);
  BOOST_FOREACH(int_int_t pr,m_des_conn[q]) m_asc_conn[pr.first].erase(q);
  BOOST_FOREACH(int_int_t pr,m_asc_conn[q]) m_des_conn[pr.first].erase(q);

  m_cp_is_cancelled[p] =true;
  m_cp_is_cancelled[q] =true;

  m_asc_conn[p].clear();
  m_des_conn[q].clear();
}

/*---------------------------------------------------------------------------*/

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

/*---------------------------------------------------------------------------*/

inline bool persistence_lt(const mscomplex_t &msc,int_pair_t p0,int_pair_t p1)
{
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

/*---------------------------------------------------------------------------*/

inline bool is_in_treshold(const mscomplex_t & msc,int_pair_t e,cell_fn_t t)
{
  bool   is_epsilon_persistent = (msc.vertid(e[0]) == msc.vertid(e[1]));
  bool   is_pers_lt_t          = std::abs(msc.fn(e[0]) - msc.fn(e[1])) < t;

  return (is_epsilon_persistent || is_pers_lt_t);
}

/*---------------------------------------------------------------------------*/

template<typename T>
inline void set_vec_value(std::vector<T> & vec, int i,const T& v){vec[i] = v;}

/*---------------------------------------------------------------------------*/

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

    if( !msc.is_paired(i) && msc.m_rect.boundryCount(c) ==
        msc.m_ext_rect.boundryCount(c))
      continue;

    inc_on_ext[i] = true;

    br::for_each(msc.m_des_conn[i]|ba::map_keys,ftor);
    br::for_each(msc.m_asc_conn[i]|ba::map_keys,ftor);
  }
}

/*---------------------------------------------------------------------------*/

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

/*---------------------------------------------------------------------------*/

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

      if(is_valid_canc_edge(*this,is_inc_ext,pr) &&
         is_in_treshold(*this,pr,f_tresh))
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

    cancel_pair(p,q);

    BOOST_FOREACH(int_int_t i,m_des_conn[p])
    BOOST_FOREACH(int_int_t j,m_asc_conn[q])
    {
      int_pair_t npr(i.first,j.first);

      update_is_inc_ext(*this,is_inc_ext,npr);

      if(is_valid_canc_edge(*this,is_inc_ext,npr) &&
         is_in_treshold(*this,npr,f_tresh))
        pq.push(npr);
    }
  }
}

/*---------------------------------------------------------------------------*/

typedef std::vector<int_list_t> contrib_list_t;

template<eGDIR dir>
int_pair_t order_pair(mscomplex_ptr_t msc,int_pair_t pr);

template <int dim,int odim>
inline bool is_dim_pair(mscomplex_ptr_t msc,int_pair_t pr)
{return (msc->index(pr[0]) == dim) && (msc->index(pr[1]) == odim);}


template<eGDIR dir,int dim>
void get_contrib_cps(mscomplex_ptr_t msc,contrib_list_t& scp_contrib)
{
  const eGDIR odir = (dir == ASC)?(DES):(ASC);
  const int   odim = (dir == ASC)?(dim +1):(dim -1);

  contrib_list_t ccp_contrib;
  std::map<int,int> ccp_map,scp_map;

  // Stage 1:
  // I)   Obtain a sequence of cancellations so that
  //      a) each pair is ordered by dir ..
  //         i.e if dir is DES then look at a (2,1) as (2,1) and not (1,2)
  //      b) pairs that have not yet been cancelled are removed
  //      c) pairs whose first element have index == dim
  //
  // II)  Obtain a list of surviving dim-cps
  //
  // III) Construct separte index mappings for the relevant canceled
  //      and surviving dim cps
  //
  // IV)  Allocate memory/datastructures
  //
  // v)   Insert each surviging cp into its own contrib. i.e. it contributes
  //      to itself

  BOOST_AUTO(ccp_rng,msc->m_canc_list
             // |ba::sliced(0,msc->m_multires_version)
             |ba::transformed(bind(order_pair<dir>,msc,_1))
             |ba::filtered(bind(is_dim_pair<dim,odim>,msc,_1))
             );

  BOOST_AUTO(scp_rng, msc->cpno_range()
             |ba::filtered(bind(&mscomplex_t::is_not_paired,msc,_1))
             |ba::filtered(bind(&mscomplex_t::is_index_i_cp<dim>,msc,_1)));

  BOOST_FOREACH(int_pair_t pr,ccp_rng)
      ccp_map.insert(int_int_t(pr[0],ccp_map.size()));

  BOOST_FOREACH(int cp,scp_rng)
      scp_map.insert(int_int_t(cp,scp_map.size()));

  ccp_contrib.resize(ccp_map.size());
  scp_contrib.resize(scp_map.size());

  BOOST_FOREACH(int_int_t scp_i,scp_map)
    scp_contrib[scp_i.second].push_back(scp_i.first);

  // Stage 2:
  // This part computes for each cancelled critical point,
  // the surviving critical points to which it contributes its
  // finest resolution geometry .
  BOOST_FOREACH(int_pair_t pr,ccp_rng|ba::reversed)
  {
    int p = pr[0],q = pr[1];

    int_list_t & pcontrib = ccp_contrib[ccp_map.at(p)];

    // for each qa in the asc conn of q:
    BOOST_FOREACH(int qa, msc->m_conn[odir][q]|ba::map_keys)
    {
      // a) if qa is not paired ..
      if(msc->is_not_paired(qa))
      {
        // .. p contributes to qa.
        pcontrib.push_back(qa);
      }
      // b) if qa is paired and qa's pair and q have same index ..
      else if(msc->index(q) == msc->index(msc->pair_idx(qa)))
      {
        // .. then foreach qaqa that qa contributes to  ..
        BOOST_FOREACH(int qaqa,ccp_contrib[ccp_map.at(qa)])
        {
          // .. p contributes to qaqa.
          pcontrib.push_back(qaqa);

          // pdpd has to be a surviving cp
          ASSERT(msc->is_not_paired(qaqa));
        }
      }
    }
  }

  // Stage 3:
  // This part just reverses info computed earlier.
  // i.e each surviving critical point will have a
  // list of cancelled critical points that contribute
  // their finest res geometry to it.

  BOOST_FOREACH(int_int_t ccp_i, ccp_map)
  {
    int p = ccp_i.first;

    BOOST_FOREACH(int p_, ccp_contrib[ccp_i.second])
    {
      scp_contrib[scp_map.at(p_)].push_back(p);
    }
  }

  // Stage 4:
  // Debug sanity checks
  BOOST_FOREACH(contrib_list_t::value_type pr, scp_contrib)
  {
    BOOST_FOREACH(int ccp, pr)
    {
      ASSERTS(msc->index(ccp) == dim)
          <<"other dim cp is attempting to contribute";

      if(ccp != pr[0])
      {
        ASSERTS(msc->is_paired(ccp))
            <<"an unpaired cp is attempting to contribute";
        ASSERTS(msc->index(msc->pair_idx(ccp)) == odim)
            <<"wrong cancellation type cp is attempting to contribute";
      }
    }
  }
}

/* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

template <eGDIR dir, int dim>
inline void __collect_mfolds(mscomplex_ptr_t msc, dataset_ptr_t ds)
{
  contrib_list_t contrib;
  get_contrib_cps<dir,dim>(msc,contrib);

  #pragma omp parallel for
  for(int i = 0 ; i < contrib.size() ; ++i)
  {
    msc->m_mfolds[dir][contrib[i][0]].clear();

    cellid_list_t contrib_cells;

    br::copy(contrib[i]|ba::transformed(bind(&mscomplex_t::cellid,msc,_1)),
             back_inserter(contrib_cells));

    contrib[i].clear();

    ds->getManifold(msc->m_mfolds[dir][contrib[i][0]],contrib_cells,dim,dir);
  }
}

/* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

void mscomplex_t::collect_mfolds(eGDIR dir, int dim, dataset_ptr_t ds)
{  
  if(dir == ASC && dim == 0 )__collect_mfolds<ASC,0>(shared_from_this(),ds);
  if(dir == ASC && dim == 1 )__collect_mfolds<ASC,1>(shared_from_this(),ds);
  if(dir == ASC && dim == 2 )__collect_mfolds<ASC,2>(shared_from_this(),ds);

  if(dir == DES && dim == 1 )__collect_mfolds<DES,1>(shared_from_this(),ds);
  if(dir == DES && dim == 2 )__collect_mfolds<DES,2>(shared_from_this(),ds);
  if(dir == DES && dim == 3 )__collect_mfolds<DES,3>(shared_from_this(),ds);
}


void mscomplex_t::collect_mfolds(dataset_ptr_t ds)
{
  __collect_mfolds<ASC,0>(shared_from_this(),ds);
  __collect_mfolds<ASC,1>(shared_from_this(),ds);
  __collect_mfolds<ASC,2>(shared_from_this(),ds);

  __collect_mfolds<DES,1>(shared_from_this(),ds);
  __collect_mfolds<DES,2>(shared_from_this(),ds);
  __collect_mfolds<DES,3>(shared_from_this(),ds);
}

/*---------------------------------------------------------------------------*/

template<>
int_pair_t order_pair<DES>(mscomplex_ptr_t msc,int_pair_t pr)
{if(msc->index(pr[0]) < msc->index(pr[1]))
    std::swap(pr[0],pr[1]);return pr;}

template<>
int_pair_t order_pair<ASC>(mscomplex_ptr_t msc,int_pair_t pr)
{if(msc->index(pr[0]) > msc->index(pr[1]))
    std::swap(pr[0],pr[1]); return pr;}


/*===========================================================================*/

}
