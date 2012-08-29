#include <grid_dataset.h>
#include <map>
#include <queue>
#include <stack>
#include <set>

#include <tr1/tuple>

#include <boost/bind/bind.hpp>
#include <boost/range.hpp>

using namespace std;

namespace grid
{

template<int dim,eGDIR dir,typename range_t>
inline void mark_reachable(range_t rng,dataset_t &ds)
{
 cellid_t cets[40];

 cellid_list_t stk(boost::begin(rng),boost::end(rng));

 while(stk.size() != 0 )
 {
   cellid_t c = stk.back(); stk.pop_back();

   ASSERT(ds.getCellDim(c) == dim);

   ds.visitCell(c);

   for(cellid_t * b = cets, *e = cets + ds.get_cets<dir>(c,cets);b!=e;++b)
   {
     if(!ds.isCellCritical(*b) && ds.m_rect.contains(*b))
     {
       cellid_t p = ds.getCellPairId(*b);
       if(ds.m_rect.contains(p) && ds.getCellDim(p) == dim && !ds.isCellVisited(p) )
       {
         stk.push_back(p);
       }
     }
   }
 }
}

typedef tr1::tuple<cellid_t,int> pq_item_t;
typedef std::vector<pq_item_t>   pq_list_t;

template <int dim>
inline bool compare_pq_items(const dataset_t &ds,const pq_item_t &c1,const pq_item_t &c2)
{
  return ds.compare_cells<dim>(tr1::get<0>(c1),tr1::get<0>(c2));
}

template <int dim>
inline cellid_list_ptr_t compute_inc_pairs
(cellid_t s, const dataset_t &ds)
{
  const int pdim = dim - 1;

  BOOST_AUTO(cmp_dim,bind(compare_pq_items<dim>,boost::cref(ds),_1,_2));
  BOOST_AUTO(cmp_pdim,bind(compare_pq_items<pdim>,boost::cref(ds),_1,_2));

  priority_queue<pq_item_t,pq_list_t,typeof(cmp_dim)>  pq(cmp_dim);
  priority_queue<pq_item_t,pq_list_t,typeof(cmp_pdim)> inc_pq(cmp_pdim);

  cellid_t f[40];

  pq.push(tr1::make_tuple(s,1));

  while(pq.size() != 0 )
  {
    cellid_t c = tr1::get<0>(pq.top());

    ASSERT(ds.getCellDim(c) == dim);

    int n = 0 ;

    do {n += tr1::get<1>(pq.top()); pq.pop();}
    while(pq.size() != 0 && tr1::get<0>(pq.top()) == c);


    for(cellid_t *b = f,*e = f + ds.get_cets<GDIR_DES>(c,f);b != e; ++b)
    {
      if(ds.isCellCritical(*b))
      {
        if(ds.isCellVisited(*b))
        {
          ASSERT(ds.getCellDim(*b) == pdim);
          inc_pq.push(tr1::make_tuple(*b,n));
        }
      }
      else
      {
        cellid_t p = ds.getCellPairId(*b);

        if (p != c && ds.isCellVisited(*b) )
        {
          pq.push(tr1::make_tuple(p,n));
        }
      }
    }
  }

  cellid_list_ptr_t cp_inc_pairs(new cellid_list_t);

  while(inc_pq.size() != 0 )
  {
    cellid_t p = tr1::get<0>(inc_pq.top());

    ASSERT(ds.getCellDim(p) == pdim);

    int n = 0 ;

    do {n += tr1::get<1>(inc_pq.top()); inc_pq.pop();}
    while((inc_pq.size() != 0)  && (tr1::get<0>(inc_pq.top()) == p));

    cp_inc_pairs->push_back(s);
    cp_inc_pairs->push_back(p);

    if(n > 1)
    {
      cp_inc_pairs->push_back(s);
      cp_inc_pairs->push_back(p);
    }
  }

  return cp_inc_pairs;
}

template<typename T,typename Titer>
int check_unique_elems(Titer b,Titer e)
{
  std::multiset<T> eset(b,e);

  for(;b!=e;++b)
  {
    ensure(eset.count(*b) == 1,"duplicate elems");
  }

  return true;
}


template <int dim,eGDIR dir,bool thruReachedPairs,typename Titer,typename cmp_t>
inline void compute_mfold
(Titer b, Titer e,
 const dataset_t &ds,
 cellid_list_t &mfold,
 cmp_t cmp)
{
  priority_queue<cellid_t,cellid_list_t,cmp_t> pq(b,e,cmp);

  cellid_t f[40];

  cellid_t lc(-1,-1,-1);

  while(pq.size() != 0 )
  {
    cellid_t c = pq.top(); pq.pop();

    ASSERT(ds.getCellDim(c) == dim);

    if(lc == c) continue; lc = c;

    mfold.push_back(c);

    for(cellid_t *b = f,*e = f + ds.get_cets<dir>(c,f);b != e; ++b)
    {
      if(!ds.isCellCritical(*b) && ds.m_rect.contains(*b))
      {
        cellid_t p = ds.getCellPairId(*b);

        if (p != c && ds.getCellDim(p) == dim &&
            (!thruReachedPairs|| ds.isCellVisited(*b)) )
          pq.push(p);
      }
    }
  }

  ASSERT(check_unique_elems<cellid_t>(mfold.begin(),mfold.end()));
}

}

