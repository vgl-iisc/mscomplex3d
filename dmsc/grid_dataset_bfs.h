#include <grid_dataset.h>
#include <map>

using namespace std;

#define MUTEX_CALL(__c) \
  {static boost::mutex __mutex;\
   boost::mutex::scoped_lock scoped_lock(__mutex);\
   (__c);}

namespace grid
{

static uint ( dataset_t::*getcets[2] ) ( cellid_t,cellid_t * ) const =
{
  &dataset_t::getCellFacets,
  &dataset_t::getCellCofacets
};

namespace bfs
{

template <typename frontier_iter_t>
inline cellid_t read_frontier(frontier_iter_t it);
template <typename frontier_t>
inline void add_to_frontier(frontier_t &f, cellid_t c);
template <typename frontier_t,typename frontier_iter_t>
inline void add_to_frontier(frontier_t &f, cellid_t c,frontier_iter_t p);




typedef map<cellid_t,int>   frontier_map_t;
template <>
inline cellid_t read_frontier(frontier_map_t::iterator it)
{
  return it->first;
}
template <>
inline void add_to_frontier(frontier_map_t &f, cellid_t c)
{
  if(f.count(c) == 0)
    f[c] = 1;
}
template <>
inline void add_to_frontier(frontier_map_t &f, cellid_t c,frontier_map_t::iterator p)
{
  if(f.count(c) == 0)
    f[c] = 0;

  f[c] += p->second;
}




typedef cellid_list_t       frontier_list_t;
template <>
inline cellid_t read_frontier(frontier_list_t::iterator it)
{
  return *it;
}
template <>
inline void add_to_frontier(frontier_list_t &f, cellid_t c)
{
  f.push_back(c);
}
template <>
inline void add_to_frontier(frontier_list_t &f, cellid_t c,frontier_list_t::iterator p)
{
  f.push_back(c);
}

}


template <> inline std::ostream& log_range
  (bfs::frontier_map_t::iterator b,bfs::frontier_map_t::iterator e,
   std::ostream &os,char sep)
{
  for (;b!=e;++b) os<<b->first<<sep;
  return os;
}

namespace bfs
{

template <eGDIR dir> inline bool compare_cells(dataset_const_ptr_t ds,cellid_t p,cellid_t q);
template <> inline bool compare_cells<GDIR_DES>(dataset_const_ptr_t ds,cellid_t p,cellid_t q)
{
  return ds->compareCells(p,q);
}
template <> inline bool compare_cells<GDIR_ASC>(dataset_const_ptr_t ds,cellid_t p,cellid_t q)
{
  return ds->compareCells(q,p);
}

template <typename frontier_t,typename visit_ftor_t,typename cp_visit_ftor_t,eGDIR dir>
void do_bfs
  ( dataset_const_ptr_t ds,
    cellid_t start_cell,
    visit_ftor_t visit_ftor,
    cp_visit_ftor_t cp_visit_ftor)
{
  typedef boost::shared_ptr<frontier_t> frontier_ptr_t;

  uint dim = ds->getCellDim(start_cell);

  frontier_ptr_t frontier0(new frontier_t);
  frontier_ptr_t frontier1(new frontier_t);

  ASSERT(ds->m_rect.contains(start_cell));

  add_to_frontier(*frontier0,start_cell);

  for(int level = 0 ;;++level)
  {
    frontier_t &frontier     = ((level&1) == 0)?(*frontier0):(*frontier1);
    frontier_t &new_frontier = ((level&1) == 1)?(*frontier0):(*frontier1);

    if(frontier.size() == 0)
      break;

    for(typename frontier_t::iterator it = frontier.begin() ; it != frontier.end() ;++it)
    {
      cellid_t c = read_frontier(it);

      ASSERT(ds->m_rect.contains(c));

      cellid_t      cets[20];

      uint cet_ct = ( ds.get()->*getcets[dir] ) ( c,cets );

      for ( uint i = 0 ; i < cet_ct ; i++ )
      {
        if ( ds->isCellExterior ( cets[i] ) )
          continue;

        if ( ds->isCellCritical ( cets[i] ) )
        {
          cp_visit_ftor(cets[i]);
          continue;
        }

        cellid_t next_cell = ds->getCellPairId ( cets[i] );

        ASSERT(ds->m_rect.contains(next_cell));

        bool is_dim   = (dim  == ds->getCellDim ( next_cell ));
        bool is_vnext = compare_cells<dir>(ds,next_cell,c);

        if (is_dim && is_vnext && visit_ftor(next_cell))
        {
          add_to_frontier(new_frontier,next_cell,it);
        }
      }
    }// end for i in [0,frontier.size)

    frontier.clear();
  }// end while(frontier.size != 0 )
}

template <typename frontier_t,typename visit_ftor_t,typename cp_visit_ftor_t>
inline void do_bfs(dataset_const_ptr_t ds,cellid_t start_cell,
                   eGDIR dir,visit_ftor_t visit_ftor,
                   cp_visit_ftor_t cp_visit_ftor)
{
  if(dir == GDIR_DES)
  {
    do_bfs<frontier_t,visit_ftor_t,cp_visit_ftor_t,GDIR_DES>
    (ds,start_cell,visit_ftor,cp_visit_ftor);
  }
  else if (dir == GDIR_ASC)
  {
    do_bfs<frontier_t,visit_ftor_t,cp_visit_ftor_t,GDIR_ASC>
    (ds,start_cell,visit_ftor,cp_visit_ftor);
  }
  else
  {
    ASSERT(false&&"what the hell man!!!");
  }
}

inline bool pass_visit(cellid_t)
{
  return true;
}

// connecting cp's

inline void make_connection(dataset_ptr_t ds, mscomplex_ptr_t msc,cellid_t p,cellid_t q)
{
  if(ds->isCellPaired(p) && ds->getCellDim(ds->getCellPairId(p)) != ds->getCellDim(q))
    return;

  if(ds->isCellPaired(q) && ds->getCellDim(ds->getCellPairId(q)) != ds->getCellDim(p))
    return;

  MUTEX_CALL(msc->connect_cps(p,q));
}

void connect_cps(dataset_ptr_t ds,mscomplex_ptr_t msc,cellid_t c,eGDIR dir)
{
  try
  {
    do_bfs<frontier_list_t>(ds,c,dir,pass_visit,bind(make_connection,ds,msc,c,_1));
  }
  catch(assertion_error e)
  {
    e.push(_FFL).push(SVAR(c)).push(SVAR(ds->m_rect)).push(SVAR(msc->m_rect));
    throw;
  }
}

// visiting cells
inline bool visit_cell(dataset_ptr_t ds,cellid_t c)
{
  if (ds->isCellVisited(c))
    return false;

  ds->visitCell(c);
  return true;
}
void mark_visits(dataset_ptr_t ds,cellid_t c,eGDIR dir)
{
  do_bfs<frontier_list_t>(ds,c,dir,bind(visit_cell,ds,_1),pass_visit);
}

// connecting cps thru visited pairs.

inline bool visit_if_pair_visited(dataset_ptr_t ds,cellid_t c)
{
  return (ds->isCellVisited(c) && ds->isCellVisited(ds->getCellPairId(c)));
}

void connect_thru_visted_pairs(dataset_ptr_t ds,mscomplex_ptr_t msc,cellid_t c,eGDIR dir)
{
  do_bfs<frontier_list_t>
      (ds,c,dir,bind(visit_if_pair_visited,ds,_1),bind(make_connection,ds,msc,c,_1));
}

// collection of manifolds

typedef dataset_t::mfold_t mfold_t;

inline bool add_to_mfold(mfold_t *mfold,cellid_t c)
{
  mfold->push_back(c);

  return true;
}

void collect_manifolds(dataset_const_ptr_t ds,mfold_t *mfold,cellid_t c,eGDIR dir)
{
  do_bfs<frontier_map_t>(ds,c,dir,bind(add_to_mfold,mfold,_1),pass_visit);
}

// marking owner extrema

inline bool write_owner_extrema(cellid_t c,dataset_t::owner_array_t * own_arr,int i)
{
  (*own_arr)(c/2) = i;
  return true;
}

void mark_owner_extrema(dataset_ptr_t ds,cellid_t c,int i)
{
  ASSERT(get_cell_dim(c) == 0 || get_cell_dim(c) == 3);

  eGDIR dir = (get_cell_dim(c) == 3)?(GDIR_DES):(GDIR_ASC);

  dataset_t::owner_array_t *own_arr = (dir == GDIR_DES)?(&ds->m_owner_maxima):(&ds->m_owner_minima);

  do_bfs<frontier_list_t>(ds,c,dir,bind(write_owner_extrema,_1,own_arr,i),pass_visit);
}

}// end namsspace bfs
}// end namespace grid
