#include <boost/range/adaptors.hpp>
#include <boost/foreach.hpp>

#include <grid_dataset.h>
#include <grid_mscomplex.h>

#include <grid_dataset_cl.h>

#define static_assert BOOST_STATIC_ASSERT

using namespace std;
namespace br = boost::range;
namespace ba = boost::adaptors;

namespace grid
{

/*===========================================================================*/

/// \brief A thread safe class to connect two critical cells in the mscomplex
class mscomplex_connector_t
{
  std::map<cellid_t,int> id_cp_map;
  mscomplex_ptr_t        msc;

public:
  void init()
  {
    #pragma omp critical(mscomplex_connector_critical_section)
    {
      for(int i = 0 ; i < msc->get_num_critpts(); ++i)
        id_cp_map.insert(make_pair(msc->cellid(i),i));
    }
  }

public:
  mscomplex_connector_t(mscomplex_ptr_t msc):msc(msc){}

  void connect_cells(cellid_t c1,cellid_t c2, int m)
  {
    int i = id_cp_map.at(c1);
    int j = id_cp_map.at(c2);

    // Exceptions thrown inside critical sections cannot pass outside
    // So this tomfoolery has to be resorted to.
    std::string e_what;

    #pragma omp critical(mscomplex_connector_critical_section)
    {
      //        if (id_cp_map.size() == 0)
      //          init();
      try{msc->connect_cps(i,j,m);}
      catch(const std::exception &e){e_what = e.what();}
    }

    ENSURES(e_what.empty()) << SVAR(c1) <<SVAR(c2);
  }
};

/*===========================================================================*/




/*===========================================================================*/

template<int dim,eGDIR dir,typename range_t>
inline void mark_reachable(const range_t &rng,dataset_ptr_t ds)
{
  cellid_t cets[40];

  cellid_list_t stk(boost::begin(rng),boost::end(rng));

  while(stk.size() != 0 )
  {
    cellid_t c = stk.back(); stk.pop_back();

    ASSERT(ds->getCellDim(c) == dim);

    ds->visitCell(c);

    for(cellid_t * b = cets, *e = cets + ds->get_cets<dir>(c,cets);b!=e;++b)
    {
      if(!ds->isCellCritical(*b))
      {
        cellid_t p = ds->getCellPairId(*b);
        if(ds->getCellDim(p) == dim && !ds->isCellVisited(p))
        {
          stk.push_back(p);
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/

template <int dim,eGDIR dir>
inline bool cellid_int_pair_cmp
(dataset_ptr_t ds, cellid_int_pair_t c1,cellid_int_pair_t c2)
{return ds->compare_cells_pp_<dir,dim>(c1.first,c2.first);}

template <int dim,eGDIR dir>
inline void compute_inc_pairs_pq
(cellid_t s, dataset_ptr_t ds, mscomplex_connector_t &msc_connector)
{
  const int pdim = (dir == DES)?(dim - 1):(dim + 1);
  const bool no_vcheck =
      ((dim==1) && (dir==DES)) || (dim==0) ||
      ((dim==2) && (dir==ASC)) || (dim==3);

  auto cmp_dim = bind(cellid_int_pair_cmp< dim, dir>,ds,_1,_2);
  auto cmp_pdim = bind(cellid_int_pair_cmp<pdim, dir>,ds,_1,_2);

  priority_queue<cellid_int_pair_t,cellid_int_pair_list_t,typeof(cmp_dim) >
      pq(cmp_dim );
  priority_queue<cellid_int_pair_t,cellid_int_pair_list_t,typeof(cmp_pdim)>
      inc_pq(cmp_pdim);

  cellid_t f[40];

  pq.push(make_pair(s,1));

  while(pq.size() != 0 )
  {
    cellid_t c = pq.top().first;

    ASSERT(ds->getCellDim(c) == dim);

    int n = 0 ;

    do {n += pq.top().second; pq.pop();}
    while(pq.size() != 0 && pq.top().first == c);

    for(cellid_t *b = f,*e = f + ds->get_cets<dir>(c,f);b != e; ++b)
    {
      if(ds->isCellCritical(*b))
      {
        ASSERT(ds->getCellDim(*b) == pdim);
        inc_pq.push(make_pair(*b,n));
      }
      else
      {
        cellid_t p = ds->getCellPairId(*b);

        if ((p != c) &&
            (no_vcheck||ds->isCellVisited(*b)) &&
            (ds->getCellDim(p) == dim ))
        {
          pq.push(make_pair(p,n));
        }
      }
    }
  }

  while(inc_pq.size() != 0 )
  {
    cellid_t p = inc_pq.top().first;

    ASSERT(ds->getCellDim(p) == pdim);

    int n = 0 ;

    do {n += inc_pq.top().second; inc_pq.pop();}
    while(inc_pq.size() != 0 && inc_pq.top().first == p);

    msc_connector.connect_cells(s,p,n);
  }
}

/*---------------------------------------------------------------------------*/

template <int dim,eGDIR dir,bool no_vcheck,typename range_t >
inline void collect_reachable_saddle
(mfold_t &mfold,const range_t &rng, dataset_ptr_t ds_ptr)
{
  const dataset_t & ds = *ds_ptr;

  const int pdim = (dir == DES)?(dim - 1):(dim + 1);

  auto cmp_dim = bind(cellid_int_pair_cmp< dim, dir>,ds_ptr,_1,_2);

  priority_queue<cellid_int_pair_t,cellid_int_pair_list_t,typeof(cmp_dim) >
      pq(cmp_dim);

  BOOST_FOREACH(cellid_t c, rng) {pq.push(cellid_int_pair_t(c,1));}

  cellid_t f[40];

  while(pq.size() != 0 )
  {
    cellid_t c = pq.top().first;

    ASSERT(ds.getCellDim(c) == dim);

    int n = 0 ;

    do {n += pq.top().second; pq.pop();}
    while(pq.size() != 0 && pq.top().first == c);

    mfold.push_back(c);

    for(cellid_t *b = f,*e = f + ds.get_cets<dir>(c,f);b != e; ++b)
    {
      if(!ds.isCellCritical(*b))
      {
        cellid_t p = ds.getCellPairId(*b);

        if (p != c && (no_vcheck||ds.isCellVisited(*b)) && ds.getCellDim(p) == dim )
        {
          pq.push(make_pair(p,n));
        }
      }
    }
  }
}

/*---------------------------------------------------------------------------*/

template <eGDIR dir,typename range_t>
inline void collect_reachable_extrema
(mfold_t &mfold,const range_t &rng, dataset_ptr_t ds)
{
  const int dim = (dir == DES)?(gc_grid_dim):(0);

  cellid_t cets[40];

  cellid_list_t stk(boost::begin(rng),boost::end(rng));

  while(stk.size() != 0 )
  {
    cellid_t c = stk.back(); stk.pop_back();

    ASSERT(ds->getCellDim(c) == dim);

    mfold.push_back(c);

    for(cellid_t * b = cets, *e = cets + ds->get_cets<dir>(c,cets);b!=e;++b)
    {
      if(!ds->isCellCritical(*b))
      {
        cellid_t p = ds->getCellPairId(*b);
        if(ds->getCellDim(p) == dim && p!= c)
        {
          stk.push_back(p);
        }
      }
    }
  }
}

/*===========================================================================*/



/*===========================================================================*/

template<int dim,eGDIR dir>
inline bool is_required_cp(const mscomplex_t& msc,int i)
{
  static const int pdim = (dir == GDIR_DES)? (dim - 1):(dim +1);

  return msc.index(i) == dim && msc.m_rect.contains(msc.cellid(i))
      && (!msc.is_paired(i) || msc.index(msc.pair_idx(i)) == pdim);
}

template <int dim, eGDIR dir>
void computeConnections(mscomplex_ptr_t msc,dataset_ptr_t ds,
                        mscomplex_connector_t &msc_connector)
{
  ENSURES(!(dim == 1 && dir == ASC));

  if(dim == 2 && dir == DES)
  {
    cellid_list_t cps_1asc;

    br::copy(msc->cpno_range()|
             ba::filtered(bind(is_required_cp<1,GDIR_ASC>,boost::cref(*msc),_1))|
             ba::transformed(bind(&mscomplex_t::cellid,msc,_1)),
             std::back_inserter(cps_1asc));

    mark_reachable<1,GDIR_ASC,typeof(cps_1asc)>(cps_1asc,ds);
  }

  cellid_list_t cps;

  br::copy(msc->cpno_range()|
           ba::filtered(bind(is_required_cp<dim,dir>,boost::cref(*msc),_1))|
           ba::transformed(bind(&mscomplex_t::cellid,boost::cref(msc),_1)),
           std::back_inserter(cps));

  #pragma omp parallel for
  for(int i = 0 ; i < cps.size(); ++i)
  {
    compute_inc_pairs_pq<dim,dir>(cps[i],ds,msc_connector);
  }
}

/*---------------------------------------------------------------------------*/

template<eGDIR ex_dir>
void computeExtremaConnections(mscomplex_ptr_t msc,dataset_ptr_t ds,
                               mscomplex_connector_t &mscomplex_connector)
{
  const int   sad_dim = (ex_dir == DES)?(2):(1);
  const int   ex_dim  = (ex_dir == DES)?(3):(0);
  const eGDIR sad_dir = (ex_dir == DES)?(ASC):(DES);

  cellid_list_t sad_cps;

  br::copy(msc->cpno_range()|
           ba::filtered(bind(is_required_cp<sad_dim,sad_dir>,boost::cref(*msc),_1))|
           ba::transformed(bind(&mscomplex_t::cellid,boost::cref(msc),_1)),
           std::back_inserter(sad_cps));

  #pragma omp parallel for
  for(int i = 0 ; i < sad_cps.size(); ++i)
  {
    cellid_t c = sad_cps[i],e1,e2;

    get_adj_extrema<sad_dim>(c,e1,e2);

    if(ds->m_rect.contains(e1))
    {
      e1 = i_to_c2(ds->get_extrema_rect<ex_dir>(),ds->owner_extrema<ex_dir>()(e1/2));
      mscomplex_connector.connect_cells(c,e1,1);
    }

    if(ds->m_rect.contains(e2))
    {
      e2 = i_to_c2(ds->get_extrema_rect<ex_dir>(),ds->owner_extrema<ex_dir>()(e2/2));
      mscomplex_connector.connect_cells(c,e2,1);
    }
  }
}

/*---------------------------------------------------------------------------*/

void  dataset_t::computeMsGraph(mscomplex_ptr_t msc)
{
  opencl::worker w;
  w.assign_gradient(shared_from_this(),msc);

  mscomplex_connector_t msc_connector(msc);
  msc_connector.init();

  #pragma omp sections
  {
    #pragma omp section
    {computeConnections<2,DES>(msc,shared_from_this(),msc_connector);}

    #pragma omp section
    {
      if(opencl::is_gpu_context())
      {
        w.owner_extrema(shared_from_this());
        computeExtremaConnections<DES>(msc,shared_from_this(),msc_connector);
        computeExtremaConnections<ASC>(msc,shared_from_this(),msc_connector);
      }
      else
      {
        computeConnections<2,ASC>(msc,shared_from_this(),msc_connector);
        computeConnections<1,DES>(msc,shared_from_this(),msc_connector);
      }
    }
  }
}

/*---------------------------------------------------------------------------*/

void  dataset_t::getManifold
(mfold_t &mfold,const cellid_list_t &rng,int dim,eGDIR dir)
{
  if(dim == 0 && dir == ASC) collect_reachable_extrema<  ASC>      (mfold,rng,shared_from_this());
  if(dim == 1 && dir == ASC) collect_reachable_saddle <1,ASC,true> (mfold,rng,shared_from_this());
  if(dim == 2 && dir == ASC) collect_reachable_saddle <2,ASC,true> (mfold,rng,shared_from_this());

  if(dim == 1 && dir == DES) collect_reachable_saddle <1,DES,true> (mfold,rng,shared_from_this());
  if(dim == 2 && dir == DES) collect_reachable_saddle <2,DES,true> (mfold,rng,shared_from_this());
  if(dim == 3 && dir == DES) collect_reachable_extrema<  DES>      (mfold,rng,shared_from_this());
}

/*---------------------------------------------------------------------------*/

void dataset_t::compute_owner_grad()
{
  opencl::worker w;
  w.assign_gradient(shared_from_this());

  if(opencl::is_gpu_context())
    w.owner_extrema(shared_from_this());
}

/*===========================================================================*/

}
