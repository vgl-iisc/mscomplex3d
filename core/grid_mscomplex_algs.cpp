#include <queue>

//#include <boost/foreach.hpp>
//#include <boost/range/adaptors.hpp>
//#include <boost/range/algorithm.hpp>
//#include <boost/bind.hpp>

#include<ranges>

#include <grid_dataset.h>
#include <grid_mscomplex.h>
#include <grid_dataset_cl.h>


//#include "cl.hpp"
//#include <opencl.hpp>

#if defined(__linux__)
#define __CL_ENABLE_EXCEPTIONS
#define CL_HPP_TARGET_OPENCL_VERSION 220
#endif
#include <OpenCL/opencl.hpp>

using namespace std;

//namespace br = boost::range;
//namespace ba = boost::adaptors;


#ifdef _MSC_VER
#define NOMINMAX // To prevent Windows.h from defining min and max macros
#define WIN32_LEAN_AND_MEAN
#include <Windows.h> // Include Windows-specific headers
//#include <cstddef> // This is where std::byte is defined in C++17
//#include <winnt.h>

// Define atomic increment macro for MSVC
#define ATOMIC_INCREMENT(ptr) InterlockedIncrement(ptr)

// Define atomic decrement macro for MSVC
#define ATOMIC_DECREMENT(ptr) InterlockedDecrement(ptr)
#else
// Define atomic increment macro for GCC/Clang
#define ATOMIC_INCREMENT(ptr) __sync_add_and_fetch(ptr, 1)

// Define atomic decrement macro for GCC/Clang
#define ATOMIC_DECREMENT(ptr) __sync_sub_and_fetch(ptr, 1)
#endif


namespace grid
{

/*===========================================================================*/

void mscomplex_t::cancel_pair ( int p, int q)
{
  order_pr_by_cp_index(*this,p,q);

  ENSURE(m_hversion == m_canc_list.size(),
         "Cannot cancel pair !! Ms complex resolution is not coarsest.");
  ENSURE(index(p) == index(q)+1,
         "indices do not differ by 1");
  ENSURE(m_cp_pair_idx[p] == -1 && m_cp_pair_idx[q] == -1,
         "p/q has already been paired");
  ENSURE(m_des_conn[p].count(q) == 1 && m_asc_conn[q].count(p) == 1,
         "p is not connected to q");
  ENSURE(m_des_conn[p][q] == 1 && m_asc_conn[q][p] == 1,
         "p and q are multiply connected");

  m_cp_pair_idx[p] = q;
  m_cp_pair_idx[q] = p;
  m_canc_list.push_back(int_pair_t(p,q));

  cancel_pair();
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::cancel_pair ()
{
  ENSURE(is_in_range(m_hversion,0,m_canc_list.size()),
         "invalid cancellation position");

  int p = m_canc_list[m_hversion][0];
  int q = m_canc_list[m_hversion][1];

  m_hversion++;

  ASSERT(index(p) == index(q)+1);
  ASSERT(m_cp_pair_idx[p] == q);
  ASSERT(m_cp_pair_idx[q] == p);
  ASSERT(m_des_conn[p].count(q) == 1);
  ASSERT(m_asc_conn[q].count(p) == 1);
  ASSERT(m_des_conn[p][q] == 1);
  ASSERT(m_asc_conn[q][p] == 1);

  m_des_conn[p].erase(q);
  m_asc_conn[q].erase(p);
  /*
  BOOST_FOREACH(int_int_t i,m_des_conn[p])
      BOOST_FOREACH(int_int_t j,m_asc_conn[q])
  {
    int u = i.first;
    int v = j.first;
    int m = i.second*j.second;

    ASSERT(is_canceled(u) == false);
    ASSERT(is_canceled(v) == false);

    BTRACE_ERROR(connect_cps(u,v,m));
  }
  */
  for (const auto& i : m_des_conn[p]) {
      for (const auto& j : m_asc_conn[q]) {
          int u = i.first;
          int v = j.first;
          int m = i.second * j.second;

          ASSERT(is_canceled(u) == false);
          ASSERT(is_canceled(v) == false);

          BTRACE_ERROR(connect_cps(u, v, m));
      }
  }

  for (const auto& pr : m_des_conn[p]) { m_asc_conn[pr.first].erase(p); }
  for (const auto& pr : m_asc_conn[p]) { m_des_conn[pr.first].erase(p); }
  for (const auto& pr : m_des_conn[q]) { m_asc_conn[pr.first].erase(q); }
  for (const auto& pr : m_asc_conn[q]) { m_des_conn[pr.first].erase(q); }
  //BOOST_FOREACH(int_int_t pr,m_des_conn[p]) m_asc_conn[pr.first].erase(p);
  //BOOST_FOREACH(int_int_t pr,m_asc_conn[p]) m_des_conn[pr.first].erase(p);
  //BOOST_FOREACH(int_int_t pr,m_des_conn[q]) m_asc_conn[pr.first].erase(q);
  //BOOST_FOREACH(int_int_t pr,m_asc_conn[q]) m_des_conn[pr.first].erase(q);

  m_cp_is_cancelled[p] =true;
  m_cp_is_cancelled[q] =true;
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::anticancel_pair()
{
  ENSURE(is_in_range(m_hversion-1,0,m_canc_list.size()),
         "invalid cancellation position");

  m_hversion--;

  int p = m_canc_list[m_hversion][0];
  int q = m_canc_list[m_hversion][1];

  ASSERT(index(p) == index(q)+1);
  ASSERT(m_cp_pair_idx[p] == q);
  ASSERT(m_cp_pair_idx[q] == p);

    /*
  BOOST_FOREACH(int_int_t pr,m_des_conn[p]) m_asc_conn[pr.first][p] = pr.second;
  BOOST_FOREACH(int_int_t pr,m_asc_conn[p]) m_des_conn[pr.first][p] = pr.second;
  BOOST_FOREACH(int_int_t pr,m_des_conn[q]) m_asc_conn[pr.first][q] = pr.second;
  BOOST_FOREACH(int_int_t pr,m_asc_conn[q]) m_des_conn[pr.first][q] = pr.second;

  // cps in lower of u except l
  BOOST_FOREACH(int_int_t i,m_des_conn[p])
      BOOST_FOREACH(int_int_t j,m_asc_conn[q])
  {
    int u = i.first;
    int v = j.first;
    int m = i.second*j.second;

    ASSERT(is_canceled(u) == false);
    ASSERT(is_canceled(v) == false);

    BTRACE_ERROR(connect_cps(u,v,-m));
  }

  */

  // Update m_asc_conn and m_des_conn with values from each other
  for (const auto& pr : m_des_conn[p]) {
      m_asc_conn[pr.first][p] = pr.second;
  }

  for (const auto& pr : m_asc_conn[p]) {
      m_des_conn[pr.first][p] = pr.second;
  }

  for (const auto& pr : m_des_conn[q]) {
      m_asc_conn[pr.first][q] = pr.second;
  }

  for (const auto& pr : m_asc_conn[q]) {
      m_des_conn[pr.first][q] = pr.second;
  }

  // Connect cps in lower of u except l
  for (const auto& i : m_des_conn[p]) {
      for (const auto& j : m_asc_conn[q]) {
          int u = i.first;
          int v = j.first;
          int m = i.second * j.second;

          ASSERT(is_canceled(u) == false);
          ASSERT(is_canceled(v) == false);

          BTRACE_ERROR(connect_cps(u, v, -m));
      }
  }

  BTRACE_ERROR(connect_cps(p,q,1));

  ASSERT(m_des_conn[p].count(q) == 1);
  ASSERT(m_asc_conn[q].count(p) == 1);
  ASSERT(m_des_conn[p][q] == 1);
  ASSERT(m_asc_conn[q][p] == 1);

  m_cp_is_cancelled[p] =false;
  m_cp_is_cancelled[q] =false;
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::set_hversion(int hver)
{
  int hversion_old = m_hversion;
  for(int i = m_hversion; i>hver && i>0; --i)  anticancel_pair();
  for(int i = m_hversion; i<hver && i<m_canc_list.size(); ++i)cancel_pair();
  TLOG << SVAR(hver) <<SVAR(hversion_old) << SVAR(m_hversion);
}

/*---------------------------------------------------------------------------*/

int mscomplex_t::get_hversion_nextrema(int nmax, int nmin) const
{
    /*
  int ns_max = boost::distance(cpno_range()
    |ba::filtered(bind(&mscomplex_t::is_index_i_cp<gc_grid_dim>,this,_1)));

  int ns_min = boost::distance(cpno_range()
    |ba::filtered(bind(&mscomplex_t::is_index_i_cp<0>,this,_1)));
    */
    
    int ns_max = std::ranges::distance(cpno_range()
        | std::views::filter(std::bind(&mscomplex_t::is_index_i_cp<gc_grid_dim>, this, std::placeholders::_1)));

    int ns_min = std::ranges::distance(cpno_range()
        | std::views::filter(std::bind(&mscomplex_t::is_index_i_cp<0>, this, std::placeholders::_1)));
    

  int hver = 0;

  for(hver = 0 ; hver < m_canc_list.size() ; ++hver)
  {
    int p = m_canc_list[hver][0],q = m_canc_list[hver][1];

    order_pr_by_cp_index(*this,p,q);

    if (ns_max <= nmax || ns_min <= nmin)
      break;

    if(index(p) == gc_grid_dim) --ns_max;
    if(index(q) == 0          ) --ns_min;

  }

  TLOG << SVAR(nmax) <<SVAR(nmin) << SVAR(hver);

  return hver;
}

/*---------------------------------------------------------------------------*/

inline bool is_valid_canc_edge
(const mscomplex_t &msc, int_pair_t e,cell_fn_t thr )
{
  order_pr_by_cp_index(msc,e[0],e[1]);

  if(msc.is_canceled(e[0])||msc.is_canceled(e[1]))
    return false;

  if(msc.is_paired(e[0]) || msc.is_paired(e[1]))
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

/*---------------------------------------------------------------------------*/

bool mscomplex_t::persistence_cmp(int_pair_t p0,int_pair_t p1) const
{
  order_pr_by_cp_index(*this,p0[0],p0[1]);
  order_pr_by_cp_index(*this,p1[0],p1[1]);

  cellid_t v00 = vertid(p0[0]);
  cellid_t v01 = vertid(p0[1]);
  cellid_t v10 = vertid(p1[0]);
  cellid_t v11 = vertid(p1[1]);

  cellid_t c00 = cellid(p0[0]);
  cellid_t c01 = cellid(p0[1]);
  cellid_t c10 = cellid(p1[0]);
  cellid_t c11 = cellid(p1[1]);

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

  cell_fn_t f00 = fn(p0[0]);
  cell_fn_t f01 = fn(p0[1]);
  cell_fn_t f10 = fn(p1[0]);
  cell_fn_t f11 = fn(p1[1]);

  cell_fn_t d1 = std::abs(f01-f00);
  cell_fn_t d2 = std::abs(f11-f10);

  if(d1 != d2)
    return d1 < d2;

  if(c00 != c10)
    return c00 < c10;

  return c01 < c11;
}

/*---------------------------------------------------------------------------*/

void mscomplex_t::simplify_pers(double thresh, bool is_nrm, int nmax, int nmin)
{
  DLOG << "Entered :" << SVAR(thresh) <<SVAR(is_nrm) << SVAR(nmax)<<SVAR(nmin);

  //auto cmp = bind(&mscomplex_t::persistence_cmp,this, _2, _1);
  auto cmp = std::bind(&mscomplex_t::persistence_cmp, this, std::placeholders::_2, std::placeholders::_1);
  priority_queue<int_pair_t,int_pair_list_t,decltype (cmp)> pq(cmp);

  //if(is_nrm)
  //  thresh *= (*br::max_element(m_cp_fn) - *br::min_element(m_cp_fn));

  if (is_nrm)
      thresh *= (*std::max_element(m_cp_fn.begin(),m_cp_fn.end()) - *std::min_element(m_cp_fn.begin(),m_cp_fn.end()));
	

  for(int i = 0 ;i < get_num_critpts();++i)
  {
    /*
    BOOST_FOREACH(int_int_t j,m_des_conn[i])
    {
      int_pair_t pr(i,j.first);

      if(is_valid_canc_edge(*this,pr,thresh))
        pq.push(pr);
    }
    */
	for (const auto& j : m_des_conn[i]) {
        int_pair_t pr(i, j.first);

        if (is_valid_canc_edge(*this, pr, thresh)) {
            pq.push(pr);
        }
    }
  }
  /*
  int ns_max = boost::distance(cpno_range()
    |ba::filtered(bind(&mscomplex_t::is_index_i_cp<gc_grid_dim>,this,_1))
    |ba::filtered(bind(&mscomplex_t::is_not_canceled,this,_1)));
  int ns_min = boost::distance(cpno_range()
    |ba::filtered(bind(&mscomplex_t::is_index_i_cp<0>,this,_1))
    |ba::filtered(bind(&mscomplex_t::is_not_canceled,this,_1)));
    */
  
  int ns_max = std::ranges::distance(cpno_range()
      | std::views::filter(std::bind(&mscomplex_t::is_index_i_cp<gc_grid_dim>, this, std::placeholders::_1))
      | std::views::filter(std::bind(&mscomplex_t::is_not_canceled, this, std::placeholders::_1)));

  int ns_min = std::ranges::distance(cpno_range()
      | std::views::filter(std::bind(&mscomplex_t::is_index_i_cp<0>, this, std::placeholders::_1))
      | std::views::filter(std::bind(&mscomplex_t::is_not_canceled, this, std::placeholders::_1)));
  
  while (pq.size() !=0)
  {
    int_pair_t pr = pq.top();

    pq.pop();

    if(!is_valid_canc_edge(*this,pr,thresh))
      continue;

    if (ns_max <= nmax || ns_min <= nmin)
      break;

    ASSERT(order_pair<DES>(shared_from_this(),pr) == pr);

    if(index(pr[0]) == gc_grid_dim) --ns_max;
    if(index(pr[1]) == 0          ) --ns_min;
    
    cancel_pair(pr[0],pr[1]);
    /*
    BOOST_FOREACH(int_int_t i,m_des_conn[pr[0]])
    BOOST_FOREACH(int_int_t j,m_asc_conn[pr[1]])
    {
      int_pair_t e(j.first,i.first);

      if(is_valid_canc_edge(*this,e,thresh))
        pq.push(e);
    }
    */
    
    for (const auto& i : m_des_conn[pr[0]]) {
        for (const auto& j : m_asc_conn[pr[1]]) {
            int_pair_t e(j.first, i.first);

            if (is_valid_canc_edge(*this, e, thresh)) {
                pq.push(e);
            }
        }
    }
    
  }

  m_merge_dag->update(shared_from_this());

  DLOG << "Exited  :" << SVAR(m_hversion) <<SVAR(ns_max) <<SVAR(ns_min);
}

/*---------------------------------------------------------------------------*/

int mscomplex_t::get_hversion_pers(double thresh, bool is_nrm) const
{
  int hver = 0;

  //if(is_nrm)
    //thresh *= (*br::max_element(m_cp_fn) - *br::min_element(m_cp_fn));

  if (is_nrm)
      thresh *= (*std::max_element(m_cp_fn.begin(),m_cp_fn.end()) - *std::min_element(m_cp_fn.begin(),m_cp_fn.end()));

  for(hver = 0 ; hver < m_canc_list.size() ; ++hver)
  {
    int p = m_canc_list[hver][0],q = m_canc_list[hver][1];

    order_pr_by_cp_index(*this,p,q);

    if (!(abs(fn(p)-fn(q)) < thresh))
      break;
  }

  TLOG << SVAR(thresh) << SVAR(is_nrm)<< SVAR(hver);
  return hver;
}

/*---------------------------------------------------------------------------*/

typedef std::vector<int_list_t> contrib_list_t;


template <int dim,int odim>
inline bool is_dim_pair(mscomplex_ptr_t msc,int_pair_t pr)
{return (msc->index(pr[0]) == dim) && (msc->index(pr[1]) == odim);}

/// \brief For each canceled dim-cp get a list of surviving cps it
///        contributes its dir-geometry
///
/// \note The first entry in each list is the id of the canceled cp
template<eGDIR dir,int dim>
void getCanceledCpContrib(mscomplex_ptr_t msc,contrib_list_t& ccp_contrib)
{
  const eGDIR odir = (dir == ASC)?(DES):(ASC);
  const int   odim = (dir == ASC)?(dim +1):(dim -1);

  // Stage 1:
  // I)   Obtain a sequence of cancellations so that
  //      a) each pair is ordered by dir ..
  //         i.e if dir is DES then look at a (2,1) as (2,1) and not (1,2)
  //      b) pairs that have not yet been cancelled are removed
  //      c) pairs whose first element have index == dim
  //
  // II) Construct an index mapping for the relevant canceled cps
  //
  // III)  Allocate data and ddd in the id of the canceled cp at the last

  std::map<int,int> ccp_map;
  /*
  auto ccp_rng = msc->m_canc_list
             |ba::sliced(0,msc->m_hversion)
             |ba::transformed(bind(order_pair<dir>,msc,_1))
             |ba::filtered(bind(is_dim_pair<dim,odim>,msc,_1));
  */

    
  auto ccp_rng = msc->m_canc_list
      | std::views::drop(0) // Start from index 0
      | std::views::take(msc->m_hversion) // Take `msc->m_hversion` elements
      | std::views::transform(std::bind(order_pair<dir>, msc, std::placeholders::_1))
      | std::views::filter(std::bind(is_dim_pair<dim, odim>, msc, std::placeholders::_1));
      

    /*
  BOOST_FOREACH(int_pair_t pr,ccp_rng)
      ccp_map.insert(int_int_t(pr[0],ccp_map.size()));
      */
  for (const auto& pr : ccp_rng) {
      ccp_map.insert(int_int_t(pr[0], ccp_map.size()));
  }


  ccp_contrib.resize(ccp_map.size());
  /*
  BOOST_FOREACH(int_int_t pr,ccp_map)
    ccp_contrib[pr.second].push_back(pr.first);
    */
	for (const auto& pr : ccp_map) {
      ccp_contrib[pr.second].push_back(pr.first);
  }


  // Stage 2:
  // This part computes for each cancelled critical point,
  // the surviving critical points to which it contributes its
  // finest resolution geometry .
    /*
  BOOST_FOREACH(int_pair_t pr,ccp_rng|ba::reversed)
  {
    int p = pr[0],q = pr[1];

    int_list_t & pcontrib = ccp_contrib[ccp_map.at(p)];

    // for each qa in the asc conn of q:
    BOOST_FOREACH(int qa, msc->m_conn[odir][q]|ba::map_keys)
    {
      // a) if qa is not canceled ..
      if(msc->is_not_canceled(qa))
      {
        // .. p contributes to qa.
        pcontrib.push_back(qa);
      }
      // b) if qa is paired and qa's pair and q have same index ..
      else if(msc->index(q) == msc->index(msc->pair_idx(qa)))
      {
        int_list_t &qa_contrib = ccp_contrib[ccp_map.at(qa)];

        // .. then foreach qaqa that qa contributes to  ..
        BOOST_FOREACH(int qaqa,qa_contrib|ba::sliced(1,qa_contrib.size()))
        {
          // .. p contributes to qaqa.
          pcontrib.push_back(qaqa);

          // pdpd has to be a surviving cp
          ASSERT(msc->is_not_canceled(qaqa));
        }
      }
    }
  }

  // Stage 3:
  // Debug sanity checks
  BOOST_FOREACH(int_list_t & ccp_l, ccp_contrib)
  {
    int ccp = ccp_l.front();

    ASSERTS (msc->index(ccp) == dim) << "incorrect dim";
    ASSERTS (msc->is_canceled(ccp))    << ccp << " should be cancelled";

    BOOST_FOREACH(int scp, ccp_l|ba::sliced(1,ccp_l.size()))
    {
      ASSERTS(msc->index(scp) == dim)  << "incorrect dim";
      ASSERTS(msc->is_not_canceled(scp)) << scp << " should be calcelled";
    }
  }
  */

  //for (const auto& pr : ccp_rng | ba::reversed)
  for (const auto& pr : ccp_rng | std::views::reverse) 
  {
	int p = pr[0], q = pr[1];
  	int_list_t& pcontrib = ccp_contrib[ccp_map.at(p)];
  	// For each qa in the asc conn of q:
    //for (const auto& qa : msc->m_conn[odir][q] | ba::map_keys) {
  	for (const auto& [qa, _] : msc->m_conn[odir][q]) {
        // a) If qa is not canceled ...
        if (msc->is_not_canceled(qa)) {
            // .. p contributes to qa.
            pcontrib.push_back(qa);

        }
        // b) If qa is paired and qa's pair and q have the same index ...
        else if (msc->index(q) == msc->index(msc->pair_idx(qa)))
        {
            int_list_t& qa_contrib = ccp_contrib[ccp_map.at(qa)];
            // .. then foreach qaqa that qa contributes to ..
            for (size_t idx = 1; idx < qa_contrib.size(); ++idx) {
                const int qaqa = qa_contrib[idx];
                // .. p contributes to qaqa.
                pcontrib.push_back(qaqa);
                // pdpd has to be a surviving cp
                ASSERT(msc->is_not_canceled(qaqa));
            }
        }
    }
  }


    // Stage 3: Debug sanity checks
    for (auto& ccp_l : ccp_contrib) {
        int ccp = ccp_l.front();

        ASSERTS(msc->index(ccp) == dim) << "incorrect dim";
        ASSERTS(msc->is_canceled(ccp)) << ccp << " should be canceled";

        for (size_t idx = 1; idx < ccp_l.size(); ++idx) {
            const int scp = ccp_l[idx];

            ASSERTS(msc->index(scp) == dim) << "incorrect dim";
            ASSERTS(msc->is_not_canceled(scp)) << scp << " should be canceled";
        }
    }

  
}

/* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

/// \brief For each survigin dim-cp get a list of canceled cps it
///        contributes its dir-geometry to it
///
/// \note The last entry in each list is the id of the surviving cp
template<eGDIR dir,int dim>
void getSurvivingCpContrib(mscomplex_ptr_t msc,contrib_list_t& scp_contrib)
{
  // Stage 1:
  // Construct an index mapping for each surviving dim-cp and allocate data
  std::map<int,int> scp_map;
  /*
  auto scp_rng = msc->cpno_range()
             |ba::filtered(bind(&mscomplex_t::is_not_canceled,msc,_1))
             |ba::filtered(bind(&mscomplex_t::is_index_i_cp<dim>,msc,_1));
  */
	auto scp_rng = msc->cpno_range()
	| std::views::filter([msc](auto cp) { return msc->is_not_canceled(cp); })
	| std::views::filter([msc](auto cp) { return msc->is_index_i_cp<dim>(cp); });
    /*
  BOOST_FOREACH(int cp,scp_rng)
      scp_map.insert(int_int_t(cp,scp_map.size()));

  scp_contrib.resize(scp_map.size());

  BOOST_FOREACH(int_int_t pr,scp_map)
    scp_contrib[pr.second].push_back(pr.first);
    */

  for (int cp : scp_rng) {
      scp_map.insert(int_int_t(cp, scp_map.size()));
  }

  scp_contrib.resize(scp_map.size());

  for (const auto& pr : scp_map) {
      scp_contrib[pr.second].push_back(pr.first);
  }


  // Stage 2:
  // For each canceld dim cp, get a list of surviging cps it contributes to
  contrib_list_t ccp_contrib;
  getCanceledCpContrib<dir,dim>(msc,ccp_contrib);

  // Stage 3:
  // Invert the above data i.e.  for each surviving cp get a list of
  // canceled cps that contribute to it.
  /*BOOST_FOREACH(int_list_t & ccp_l, ccp_contrib)
  {
    int ccp = ccp_l.front();

    BOOST_FOREACH(int scp, ccp_l|ba::sliced(1,ccp_l.size()))
    {
      scp_contrib[scp_map.at(scp)].push_back(ccp);
    }
  }*/
  for (auto& ccp_l : ccp_contrib)
  {
      int ccp = ccp_l.front();

      /*for (int scp : ccp_l | ba::sliced(1, ccp_l.size())) {
          scp_contrib[scp_map.at(scp)].push_back(ccp);
      }*/
      for (int scp : std::ranges::subrange(ccp_l.begin() + 1, ccp_l.end())) {
          scp_contrib[scp_map.at(scp)].push_back(ccp);
      }
  }


  // Stage 4:
  // Debug sanity checks
    /*
  BOOST_FOREACH(int_list_t & scp_l, scp_contrib)
  {
    int scp = scp_l.front();

    ASSERTS (msc->index(scp) == dim)   << "incorrect dim";
    ASSERTS (msc->is_not_canceled(scp))  << scp << " should be cancelled";

    BOOST_FOREACH(int ccp, scp_l|ba::sliced(1,scp_l.size()))
    {
      ASSERTS(msc->index(ccp) == dim)     << "incorrect dim";
      ASSERTS(msc->is_canceled(ccp)) << ccp << " should be calcelled";
    }
  }
  */

  for (auto& scp_l : scp_contrib)
  {
      int scp = scp_l.front();

      ASSERTS(msc->index(scp) == dim) << "incorrect dim";
      ASSERTS(msc->is_not_canceled(scp)) << scp << " should be cancelled";

      /*
      for (int ccp : scp_l | ba::sliced(1, scp_l.size()))
      {
          ASSERTS(msc->index(ccp) == dim) << "incorrect dim";
          ASSERTS(msc->is_canceled(ccp)) << ccp << " should be cancelled";
      }
      */
      for (int ccp : std::ranges::subrange(scp_l.begin() + 1, scp_l.end())) {
          ASSERTS(msc->index(ccp) == dim) << "incorrect dim";
          ASSERTS(msc->is_canceled(ccp)) << ccp << " should be cancelled";
      }
  }

}


/* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

template <eGDIR dir, int dim>
inline void __collect_mfolds(mscomplex_ptr_t msc, dataset_ptr_t ds)
{
  contrib_list_t contrib;
  getSurvivingCpContrib<dir,dim>(msc,contrib);

  #pragma omp parallel for
  for(int i = 0 ; i < contrib.size() ; ++i)
  {
    msc->m_mfolds[dir][contrib[i][0]].clear();

    cellid_list_t contrib_cells;


    //br::copy(contrib[i]|ba::transformed(bind(&mscomplex_t::cellid,msc,_1)),
      //  back_inserter(contrib_cells));

    // Transform each element using the `cellid` function
    std::transform(
        contrib[i].begin(),
        contrib[i].end(),
        std::back_inserter(contrib_cells),
        [&msc](const auto& elem) {
            return msc->cellid(elem);
        });


    contrib[i].clear();

    ds->getManifold(msc->m_mfolds[dir][contrib[i][0]],contrib_cells,dim,dir);
  }

  msc->m_geom_hversion[dir][dim] = msc->get_hversion();
}

/* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

template <eGDIR dir>
inline void __collect_extrema_mfolds(mscomplex_ptr_t msc, dataset_ptr_t ds)
{
  const int dim = (dir ==ASC)?(0):(gc_grid_dim);

  int_list_t scp_id(msc->get_num_critpts());
  int_list_t scp_ncells(msc->get_num_critpts(),0);

  rect_t    r        = ds->get_extrema_rect<dir>();
  int_marray_t &oarr = ds->owner_extrema<dir>();
  cellid_t lc        = r.lc();
  cellid_t uc        = r.uc();

  {
    contrib_list_t ccp_contrib;
    getCanceledCpContrib<dir,dim>(msc,ccp_contrib);

    #pragma omp parallel for
    for(int i = 0 ; i< scp_id.size(); ++i)
      scp_id[i] = i;

    #pragma omp parallel for
    for(int i = 0 ; i < ccp_contrib.size(); ++i)
      scp_id[ccp_contrib[i][0]] =  ccp_contrib[i][1];
  }

  #pragma omp parallel for
  for(int i = 0 ; i< msc->get_num_critpts(); ++i)
    if(msc->index(i)==dim)
      oarr(msc->cellid(i)/2) = i;




  #pragma omp parallel for
  for(    int z = lc[2]; z <= uc[2]; z +=2 )
    for(  int y = lc[1]; y <= uc[1]; y +=2 )
      for(int x = lc[0]; x <= uc[0]; x +=2 )
      {
        cellid_t c(x,y,z);
        int oe   = oarr(c/2 );
        int oe_i = (ds->isCellCritical(c))?(oe):(oarr(i_to_c2(r,oe)/2));
        int i    = scp_id[oe_i];
        //InterlockedIncrement(&scp_ncells[i]);
        //InterlockedIncrement(reinterpret_cast<LONG*>(&scp_ncells[i]));
        //ATOMIC_INCREMENT(reinterpret_cast<LONG*>(&scp_ncells[i]));
        #ifdef _MSC_VER
        ATOMIC_INCREMENT(reinterpret_cast<LONG*>(&scp_ncells[i]));
        #endif
      }

  #pragma omp parallel for
  for(int i = 0 ; i< msc->get_num_critpts(); ++i)
    if(msc->index(i)==dim && scp_ncells[i] != 0)
      msc->m_mfolds[dir][i].resize(scp_ncells[i]);


  #pragma omp parallel for
  for(    int z = lc[2]; z <= uc[2]; z +=2 )
    for(  int y = lc[1]; y <= uc[1]; y +=2 )
      for(int x = lc[0]; x <= uc[0]; x +=2 )
      {
        cellid_t c(x,y,z);
        int oe   = oarr(c/2 );
        int oe_i = (ds->isCellCritical(c))?(oe):(oarr(i_to_c2(r,oe)/2));
        int i    = scp_id[oe_i];
        //int p    = InterlockedDecrement(&scp_ncells[i]);
        //int p    = InterlockedDecrement(reinterpret_cast<LONG*>(&scp_ncells[i]));
        //int p = atomic_fetch_add(&scp_ncells[i], 1);
        int p=0;
        #ifdef _MSC_VER
        p=ATOMIC_DECREMENT(reinterpret_cast<LONG*>(&scp_ncells[i]));
        #else
        p=ATOMIC_DECREMENT(&scp_ncells[i]);
        #endif
        msc->m_mfolds[dir][i][p] = c;
      }

  #pragma omp parallel for
  for(int i = 0 ; i< msc->get_num_critpts(); ++i)
    if(msc->index(i)==dim)
      oarr(msc->cellid(i)/2) = c_to_i2(r,msc->cellid(i));

  msc->m_geom_hversion[dir][dim] = msc->get_hversion();
}


/* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

void mscomplex_t::collect_mfolds(eGDIR dir, int dim, dataset_ptr_t ds)
{
  ENSURES(dir==DES || dir==ASC || dir==GDIR_CT) <<SVAR(dir);
  ENSURES(dim == -1 || (0<=dim && dim <= gc_grid_dim) ) <<SVAR(dim);

  if (dir == ASC || dir == GDIR_CT)
  {
    if(dim==0 || dim==-1)
    {
      if(opencl::is_gpu_context())
        __collect_extrema_mfolds<ASC>(shared_from_this(),ds);
      else
        __collect_mfolds<ASC,0>(shared_from_this(),ds);
    }
    if(dim==1 || dim==-1)__collect_mfolds<ASC,1>(shared_from_this(),ds);
    if(dim==2 || dim==-1)__collect_mfolds<ASC,2>(shared_from_this(),ds);
  }

  if (dir == DES || dir == GDIR_CT)
  {
    if(dim==1 || dim==-1)__collect_mfolds<DES,1>(shared_from_this(),ds);
    if(dim==2 || dim==-1)__collect_mfolds<DES,2>(shared_from_this(),ds);
    if(dim==3 || dim==-1)
    {
      if(opencl::is_gpu_context())
        __collect_extrema_mfolds<DES>(shared_from_this(),ds);
      else
        __collect_mfolds<DES,3>(shared_from_this(),ds);
    }
  }
}

/* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - */

void mscomplex_t::collect_mfolds(dataset_ptr_t ds)
{
  if(opencl::is_gpu_context())
  {
    __collect_extrema_mfolds<ASC>(shared_from_this(),ds);
    __collect_extrema_mfolds<DES>(shared_from_this(),ds);
  }
  else
  {
    __collect_mfolds<ASC,0>(shared_from_this(),ds);
    __collect_mfolds<DES,3>(shared_from_this(),ds);
  }

  __collect_mfolds<ASC,1>(shared_from_this(),ds);
  __collect_mfolds<ASC,2>(shared_from_this(),ds);

  __collect_mfolds<DES,1>(shared_from_this(),ds);
  __collect_mfolds<DES,2>(shared_from_this(),ds);
}

/*===========================================================================*/

}
