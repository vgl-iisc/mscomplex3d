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

typedef tr1::tuple<cellid_t,int> pq_item_t;
typedef std::vector<pq_item_t>   pq_list_t;

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

//  ASSERT(check_unique_elems<cellid_t>(mfold.begin(),mfold.end()));
}

}

