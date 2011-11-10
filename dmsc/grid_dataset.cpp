#include <stack>
#include <queue>
#include <list>
#include <set>

#include <fstream>

#include <boost/typeof/typeof.hpp>
#include <boost/function.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/static_assert.hpp>
#include <boost/iterator/filter_iterator.hpp>

#define static_assert BOOST_STATIC_ASSERT

#include <config.h>

#include <grid_dataset.h>
#include <grid_mscomplex.h>
#include <grid_dataset_bfs.h>

#ifdef BUILD_EXEC_OPENCL
#include <grid_dataset_cl.h>
#endif

using namespace std;

namespace grid
{
  inline uint   dataset_t::getCellIncCells( cellid_t c,cellid_t * inc) const
  {
    for(uint i = 0; i < gc_grid_dim; ++i)
    {
      inc[i*2+0] = c;
      inc[i*2+1] = c;

      inc[i*2+0][i] -= 1;
      inc[i*2+0][i] += 1;
    }
    return gc_grid_dim*2;
  }

  dataset_t::dataset_t (const rect_t &r,const rect_t &e,const rect_t &d) :
      m_rect (r),
      m_ext_rect (e),
      m_domain_rect(d),
      m_cell_flags(cellid_t::zero,boost::fortran_storage_order()),
      m_vert_fns(cellid_t::zero,boost::fortran_storage_order()),
      m_owner_maxima(cellid_t::zero,boost::fortran_storage_order()),
      m_owner_minima(cellid_t::zero,boost::fortran_storage_order())
  {
    // TODO: assert that the given rect is of even size..
    //       since each vertex is in the even positions


    cmp_ftor     = bind(&dataset_t::compareCells,this,_1,_2);
    cmp_ftors[0] = bind(&dataset_t::compareCells,this,_1,_2);
    cmp_ftors[1] = bind(&dataset_t::compareCells,this,_2,_1);
  }

  dataset_t::~dataset_t ()
  {
    clear();
  }

  void dataset_t::init(const std::string &filename)
  {
    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    rect_size_t   span   = m_ext_rect.span() + 1;
    rect_size_t  pt_span = (m_ext_rect.span()/2)+1;
    rect_point_t bl = m_ext_rect.lower_corner();

    m_cell_flags.resize(span);
    m_vert_fns.resize(pt_span);

    uint num_cells = span[0]*span[1]*span[2];
    uint num_pts   = pt_span[0]*pt_span[1]*pt_span[2];

    m_cell_flags.reindex(bl);
    m_vert_fns.reindex(bl/2);

    std::fill_n(m_cell_flags.data(),num_cells,0);

    ifstream ifs(filename.c_str(),ios::in|ios::binary);
    ensure(ifs.is_open(),"unable to open file");

    ifs.read((char*)(void*)m_vert_fns.data(),sizeof(cell_fn_t)*num_pts);
    ensure(ifs.fail()==false,"failed to read some data");

    ifs.seekg(0,ios::end);
    ensure(ifs.tellg()==num_pts*sizeof(cell_fn_t),"file/piece size mismatch");

    ifs.close();

    m_owner_maxima.resize(m_rect.span()/2);
    m_owner_minima.resize((m_rect.span()/2)+1);

    m_owner_maxima.reindex(m_rect.lc()/2);
    m_owner_minima.reindex(m_rect.lc()/2);


  }

  void  dataset_t::clear()
  {
    m_cell_flags.resize(cellid_t::zero);

    m_vert_fns.resize(cellid_t::zero);

    m_owner_maxima.resize(cellid_t::zero);
    m_owner_minima.resize(cellid_t::zero);
  }

  inline cellid_t flag_to_mxfct(cellid_t c,cell_flag_t f)
  {
    cell_flag_t d = f&0x07;
    ASSERT(is_in_range(d,1,7));
    c[(d-1)>>1] += (d&1)?(-1):(+1);
    return c;
  }

  inline cell_flag_t mxfct_to_flag(cellid_t c,cellid_t fct)
  {
    ASSERT(euclid_norm2(c-fct) == 1);

    int d = 0;

    if(c[1] != fct[1])
      d = 1;
    else if(c[2] != fct[2])
      d = 2;

    if(c[d] > fct[d])
      return (1 + d*2 + 0);
    else
      return (1 + d*2 + 1);
  }

  inline cellid_t flag_to_pair(cellid_t c,cell_flag_t f)
  {
    return flag_to_mxfct(c,(f>>3)&0x07);
  }

  inline cell_flag_t pair_to_flag(cellid_t c,cellid_t p)
  {
    return (mxfct_to_flag(c,p)<<3);
  }

  bool dataset_t::areCellsIncident(cellid_t c1,cellid_t c2) const
  {
    return ( euclid_norm2(c1-c2) == 1);
  }

  cellid_t dataset_t::getCellPairId (cellid_t c) const
  {
    ASSERT(isCellPaired(c));

    return flag_to_pair(c,m_cell_flags(c));
  }

  cellid_t dataset_t::getCellMaxFacetId (cellid_t c) const
  {
    return flag_to_mxfct(c,m_cell_flags(c));
  }

  cellid_t dataset_t::getCellSecondMaxFacetId (cellid_t c) const
  {
    return (2*c - flag_to_mxfct(c,m_cell_flags(c)));
  }

  inline bool compareCells_dim
      (const dataset_t * ds,
       const cellid_t &c1,
       const cellid_t &c2,
       const int & dim)
  {
    if(dim == 0)
      return ds->ptLt(c1,c2);

    cellid_t f1 = ds->getCellMaxFacetId(c1);
    cellid_t f2 = ds->getCellMaxFacetId(c2);

    if(f1 != f2)
      return compareCells_dim(ds,f1,f2,dim-1);

    f1 = ds->getCellSecondMaxFacetId(c1);
    f2 = ds->getCellSecondMaxFacetId(c2);

    int boundry_ct1 = ds->m_domain_rect.boundryCount(f1);
    int boundry_ct2 = ds->m_domain_rect.boundryCount(f2);

    if(boundry_ct1 != boundry_ct2)
      return (boundry_ct2 < boundry_ct1);

    return compareCells_dim(ds,f1,f2,dim-1);
  }

  bool dataset_t::compareCells( cellid_t c1,cellid_t  c2 ) const
  {
    int c1_dim = getCellDim(c1);
    int c2_dim = getCellDim(c2);

    cellid_t &p = (c1_dim > c2_dim)?(c1):(c2);

    for(int i = 0 ; i < std::abs(c1_dim-c2_dim);++i)
      p = getCellMaxFacetId(p);

    if(c1 == c2)
      return (c1_dim < c2_dim);

    return compareCells_dim(this,c1,c2,std::min(c1_dim,c2_dim));
  }

  uint dataset_t::getCellPoints (cellid_t c,cellid_t  *p) const
  {
    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    static_assert(((cell_coord_t)-1) < 0 && "coord_t needs to support -1 ");

    uint pos = 0;

    cellid_t i;

    for(i[2] = -(c[2]&1) ; i[2] <= (c[2]&1) ;i[2]+=2)
    {
      for(i[1] = -(c[1]&1) ; i[1] <= (c[1]&1) ;i[1]+=2)
      {
        for(i[0] = -(c[0]&1) ; i[0] <= (c[0]&1) ;i[0]+=2)
        {
          p[pos++] = c+i;
        }
      }
    }

    return (1<<getCellDim (c));
  }

  uint dataset_t::getCellFacets (cellid_t c,cellid_t *f) const
  {
    uint pos = 0;

    for(uint d = 0; d< gc_grid_dim; ++d)
    {
      for(uint i = 0 ; i < (c[d]&1);++i)
      {
        f[pos] = c; f[pos++][d] += 1;
        f[pos] = c; f[pos++][d] -= 1;
      }
    }
    return getCellDim (c)*2;
  }

  uint dataset_t::getCellCofacets (cellid_t c,cellid_t *cf) const
  {
    uint cf_ct = (gc_grid_dim - getCellDim (c))*2 ;

    uint pos = 0;

    for(uint d = 0; d< gc_grid_dim; ++d)
    {
      for(uint i = 0 ; i < ((c[d]+1)&1);++i)
      {

        cf[pos] = c; cf[pos++][d] += 1;
        cf[pos] = c; cf[pos++][d] -= 1;
      }
    }

    uint cf_nv_pos = 0;

    for (uint i = 0 ;i < cf_ct;++i)
      if (m_ext_rect.contains (cf[i]))
        cf[cf_nv_pos++] = cf[i];

    return cf_nv_pos;

  }

  uint dataset_t::getCellCofaces (cellid_t c,cellid_t *cf) const
  {
    uint cf_ct = std::pow(3,(gc_grid_dim - getCellDim (c))) ;

    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    static_assert(((cell_coord_t)-1) < 0 && "coord_t needs to support -1 ");

    uint pos = 0;

    cellid_t i,l = (c+cellid_t::one)&(cellid_t::one);

    for(i[2] = -(l[2]) ; i[2] <= (l[2]) ;i[2]+=1)
    {
      for(i[1] = -(l[1]) ; i[1] <= (l[1]) ;i[1]+=1)
      {
        for(i[0] = -(l[0]) ; i[0] <= (l[0]) ;i[0]+=1)
        {
          cf[pos++] = c+i;
        }
      }
    }

    uint cf_nv_pos = 0;

    for (uint i = 0 ;i < cf_ct;++i)
      if (m_ext_rect.contains (cf[i]))
        cf[cf_nv_pos++] = cf[i];

    return cf_nv_pos;
  }

  uint dataset_t::getCellEst (cellid_t c,cellid_t* est)  const
  {
    cellid_t cfs[40];

    int cfs_ct = getCellCofaces(c,cfs);

    ASSERT(is_in_range(cfs_ct,0,40));

    uint pos = 0;

    for(int i = 0 ; i< cfs_ct;++i)
    {
      cellid_t cf = cfs[i];

      ASSERT(m_ext_rect.contains(cf));

      uint c_dim  = getCellDim(c);
      uint cf_dim = getCellDim(cf);

      for(uint j = c_dim; j < cf_dim ; ++j)
        cf = getCellMaxFacetId(cf);

      if(cf == c)
        est[pos++] = cfs[i];
    }

    return pos;
  }

  bool dataset_t::isPairOrientationCorrect (cellid_t c, cellid_t p) const
  {
    return (getCellDim (c) <getCellDim (p));
  }

//  bool dataset_t::isCellMarked (cellid_t c) const
//  {
//    return ! (m_cell_flags (c) == CELLFLAG_UNKNOWN);
//  }

  bool dataset_t::isCellCritical (cellid_t c) const
  {
    return (m_cell_flags (c) & CELLFLAG_CRITICAL);
  }

  bool dataset_t::isCellPaired (cellid_t c) const
  {
    return (((m_cell_flags(c)>>3) & 0x07) !=0);
  }

  bool dataset_t::isCellVisited (cellid_t c) const
  {
    return (m_cell_flags (c) & CELLFLAG_VISITED);
  }

  void dataset_t::pairCells (cellid_t c,cellid_t p)
  {
    ASSERT(isCellPaired(c) == false);
    ASSERT(isCellPaired(p) == false);

    m_cell_flags (c) |= (pair_to_flag(c,p));
    m_cell_flags (p) |= (pair_to_flag(p,c));

    ASSERT(getCellPairId(c) == p);
    ASSERT(getCellPairId(p) == c);
  }

  void dataset_t::visitCell(cellid_t c)
  {
    ASSERT(isCellVisited(c) == false);
    m_cell_flags (c) |= CELLFLAG_VISITED;
    ASSERT(isCellVisited(c) == true);
  }

//  void  dataset_t::unpairCells ( cellid_t c,cellid_t p )
//  {
//    ASSERT(getCellPairId(c) == p);
//    ASSERT(getCellPairId(p) == c);

//    m_cell_pairs (c) = CELLADJDIR_UNKNOWN;
//    m_cell_pairs (p) = CELLADJDIR_UNKNOWN;

//    m_cell_flags (c) &= (CELLFLAG_MASK^CELLFLAG_PAIRED);
//    m_cell_flags (p) &= (CELLFLAG_MASK^CELLFLAG_PAIRED);

//    ASSERT(isCellPaired(c) == false);
//    ASSERT(isCellPaired(p) == false);
//  }

  void dataset_t::setCellMaxFacet (cellid_t c,cellid_t f)
  {
    ASSERT(getCellDim(c) == getCellDim(f)+1);
    m_cell_flags (c) |= mxfct_to_flag(c,f);
    ASSERT(getCellMaxFacetId(c) == f);
  }

  void dataset_t::markCellCritical (cellid_t c)
  {
    ASSERT(isCellCritical(c) == false);
    m_cell_flags (c) |= CELLFLAG_CRITICAL;
    ASSERT(isCellCritical(c) == true);
  }

  bool dataset_t::isTrueBoundryCell (cellid_t c) const
  {
    return (m_domain_rect.isOnBoundry (c));
  }

  bool dataset_t::isFakeBoundryCell (cellid_t c) const
  {
    return (m_rect.isOnBoundry (c) && (!m_domain_rect.isOnBoundry (c)));
  }

  bool dataset_t::isCellExterior (cellid_t c) const
  {
    return (!m_rect.contains (c));
  }

  void dataset_t::assignMaxFacets_thd(int thid,int dim)
  {
    cellid_t f[20],c,s(0,0,0);

    for(int i = 0 ; i < dim; ++i)
      s[i] = 1;

    while(true)
    {
      rect_t rect  = rect_t(m_ext_rect.lc()+s,m_ext_rect.uc()-s);

      int N = num_cells2(rect);

      for( int i = thid; i < N; i += g_num_threads)
      {
        c = i_to_c2(rect,i);

        ASSERT(getCellDim(c) == dim);

        int f_ct   = getCellFacets(c,f);

        setCellMaxFacet(c,*std::max_element(f,f+f_ct,cmp_ftor));
      }

      if(!next_permutation(s.rbegin(),s.rend()))
        break;
    }
  }

  void  dataset_t::pairCellsWithinEst_thd(int tid, cellid_list_t *ccells)
  {
    int N = num_cells2(m_rect);

    for(int i = tid; i < N; i += g_num_threads)
    {
      cellid_t c = i_to_c2(m_rect,i);

      cellid_t est_arr[40];

      uint est_ct = getCellEst(c,est_arr);

      std::sort(est_arr,est_arr+est_ct,cmp_ftor);

      for(int j = 0; j < est_ct; ++j)
      {
        cellid_t p = est_arr[j];

        if(j+1 < est_ct)
        {
          cellid_t q = est_arr[j+1];

          bool is_adj        = areCellsIncident(p,q);
          bool is_same_bndry = (isTrueBoundryCell(p) == isTrueBoundryCell(q));

          if(is_adj && is_same_bndry)
          {
            pairCells(p,q);
            ++j;
            continue;
          }
        }

        if(m_rect.contains(p))
        {
          markCellCritical(p);
          ccells->push_back(p);
        }
      }
    }
  }

  void  dataset_t::markBoundry_thd(int tid, rect_t bnd, cellid_list_t *ccells)
  {
    int N = num_cells(bnd);

    cellid_t bnd_nrm = bnd.get_normal();

    for(int i = tid ; i < N; i += g_num_threads)
    {
      cellid_t c = i_to_c(bnd,i);

      if(!isCellPaired(c))
        continue;

      cellid_t p = getCellPairId(c);

      if((bnd_nrm != (c-p)) && (bnd_nrm != (p-c)))
        continue;

      markCellCritical(c);
      ccells->push_back(c);

      if(!m_rect.contains(p))
        continue;

      markCellCritical(p);
      ccells->push_back(p);
    }
  }

  void  dataset_t::setupCPs(mscomplex_ptr_t msgraph, cellid_list_t *ccells, int offset)
  {
    for(int i = 0 ; i < ccells->size();++i)
    {
      cellid_t c = ccells->at(i);

      cellid_t v = get_cell_vert(c);

      msgraph->set_critpt(offset++,c,getCellDim(c),get_cell_fn(c),v);

      if(!isCellPaired(c))
        continue;

      cellid_t p = getCellPairId(c);

      if(isPairOrientationCorrect(c,p))
        continue;

      if(!m_rect.contains(p))
        continue;

      ASSERT(ccells->at(i-1) == p);

      msgraph->pair_cps(offset-1,offset-2);
    }
    ccells->clear();
  }

  void  dataset_t::extrema_connect_thd
      (mscomplex_ptr_t msgraph, cp_producer_ptr_t prd)
  {
    for(int i ; prd->next(i);)
    {
      cellid_t c = msgraph->cellid(i);

      eGDIR dir = (msgraph->index(i) == 3)?(GDIR_DES):(GDIR_ASC);

      bfs::connect_cps(shared_from_this(),msgraph,c,dir);
    }
  }

  void  dataset_t::saddle_visit(mscomplex_ptr_t msgraph,eGDIR dir)
  {
    int idx = (dir == GDIR_DES)?(2):(1);

    for(int i = 0 ; i < msgraph->get_num_critpts();++i)
    {
      if(msgraph->index(i) != idx )
        continue;

      bfs::mark_visits(shared_from_this(),msgraph->cellid(i),dir);
    }
  }

  void  dataset_t::saddle_connect_thd
      (mscomplex_ptr_t msgraph, cp_producer_ptr_t prd)
  {
    for(int i ; prd->next(i);)
    {
      cellid_t c = msgraph->cellid(i);

      bfs::connect_thru_visted_pairs(shared_from_this(),msgraph,c,GDIR_DES);
    }
  }

  void  dataset_t::computeMsGraph(mscomplex_ptr_t msgraph)
  {

#ifdef BUILD_EXEC_OPENCL
    opencl::worker w;
    w.assign_gradient(shared_from_this(),msgraph);
#else
    for(int dim = 1 ; dim <= gc_grid_dim; ++dim)
    {
      boost::thread_group group;

      for(int tid = 0 ; tid < g_num_threads; ++tid)
        group.create_thread(bind(&dataset_t::assignMaxFacets_thd,this,tid,dim));

      group.join_all();
    }

    cellid_list_t ccells[g_num_threads];

    {
      boost::thread_group group;

      for(int tid = 0 ; tid < g_num_threads; ++tid)
        group.create_thread(bind(&dataset_t::pairCellsWithinEst_thd,this,tid,ccells+tid));

      group.join_all();
    }
    {
      rect_list_t bnds;

      get_boundry_rects(m_rect,m_ext_rect,bnds);

      for( int i = 0 ; i < bnds.size() ; ++i)
      {
        boost::thread_group group;

        for(int tid = 0 ; tid < g_num_threads; ++tid)
          group.create_thread(bind(&dataset_t::markBoundry_thd,this,tid,bnds[i],ccells+tid));

        group.join_all();
      }
    }

    {
      int offset[g_num_threads+1];

      offset[0] = 0;

      for(int tid = 0 ; tid < g_num_threads; ++tid)
        offset[tid + 1] = offset[tid] + ccells[tid].size();

      msgraph->resize(offset[g_num_threads]);

      boost::thread_group group;

      for(int tid = 0 ; tid < g_num_threads; ++tid)
        group.create_thread(bind(&dataset_t::setupCPs,this,msgraph,ccells+tid,offset[tid]));

      group.join_all();
    }
#endif

    msgraph->build_id_cp_map();

    {
      boost::thread_group saddle_group;

      saddle_group.create_thread(bind(&dataset_t::saddle_visit,this,msgraph,GDIR_DES));
      saddle_group.create_thread(bind(&dataset_t::saddle_visit,this,msgraph,GDIR_ASC));

//      saddle_visit(msgraph,GDIR_DES);
//      saddle_visit(msgraph,GDIR_ASC);

#ifdef BUILD_EXEC_OPENCL
      w.owner_extrema(shared_from_this(),msgraph);
#else
      {
        boost::thread_group group;

        cp_producer_ptr_t prd(new cp_producer_t(msgraph,cp_producer_t::extrema_filter));

        for(int tid = 0 ; tid < g_num_threads-2; ++tid)
          group.create_thread(bind(&dataset_t::extrema_connect_thd,this,msgraph,prd));

        group.join_all();
      }
#endif
      saddle_group.join_all();
    }

    {
      boost::thread_group group;

      cp_producer_ptr_t prd(new cp_producer_t(msgraph,cp_producer_t::twosaddle_filter));

      for(int tid = 0 ; tid < g_num_threads; ++tid)
        group.create_thread(bind(&dataset_t::saddle_connect_thd,this,msgraph,prd));

      group.join_all();
    }
  }

  void  dataset_t::get_mfold
      (mfold_t *mfold, mscomplex_const_ptr_t msc, int i, int dir) const
  {
    try
    {
      ASSERT(msc->is_paired(i) == false);

      if(m_rect.contains(msc->cellid(i)))
        bfs::collect_manifolds
            (shared_from_this(),mfold,msc->cellid(i),(eGDIR)dir);

      for( conn_iter_t j  = msc->m_conn[dir][i].begin();
                       j != msc->m_conn[dir][i].end();++j)
      {
        ASSERT(msc->is_paired(*j) == true);
        ASSERT(msc->index(i) == msc->index(msc->pair_idx(*j)));

        bfs::collect_manifolds
            (shared_from_this(),mfold,msc->cellid(msc->pair_idx(*j)),(eGDIR)dir);
      }
    }
    catch(assertion_error e)
    {
      e.push(_FFL);
      e.push(SVAR(m_rect));
      e.push(SVAR(m_ext_rect));

      throw;
    }
  }

  void  dataset_t::mark_extrema_owner_thd(mscomplex_ptr_t msc,cp_producer_ptr_t p)
  {
    for ( int i ; p->next(i);)
    {
      bfs::mark_owner_extrema(shared_from_this(),msc->cellid(i),msc->surv_extrema(i));
    }
  }

  template<typename T>
  class producer_consumer_t:boost::noncopyable
  {

  protected:
    queue<T>                      m_queue;
    boost::mutex                  m_mutex;
    boost::condition_variable_any m_cond;

  public:
    void put(const T & t)
    {
      boost::mutex::scoped_lock scoped_lock(m_mutex);

      m_queue.push(t);

      m_cond.notify_one();
    }

    T get()
    {
      boost::mutex::scoped_lock scoped_lock(m_mutex);

      while (m_queue.empty())
        m_cond.wait(m_mutex);

      T t = m_queue.front();

      m_queue.pop();

      return t;
    }
  };

  template<typename Titer>
  class mt_producer_t:boost::noncopyable
  {
  protected:
    Titer          m_begin;
    Titer          m_end;
    boost::mutex   m_mutex;

  public:

    typedef Titer iterator;

    class no_more_items:public std::exception
    {
    public:
      no_more_items(){}
      virtual ~no_more_items() throw(){}

      virtual const char*
      what() const throw()
      {return "nothing more to produce";}
    };

    mt_producer_t(Titer begin,Titer end)
    {
      boost::mutex::scoped_lock scoped_lock(m_mutex);

      m_begin = begin;
      m_end   = end;
    }

    Titer next()
    {
      boost::mutex::scoped_lock scoped_lock(m_mutex);

      if(m_begin == m_end)
        throw no_more_items();

      return m_begin++;
    }
  };

}

namespace std
{
template<typename _IIter>
  typename iterator_traits<_IIter>::difference_type
  count(_IIter b, _IIter e)
  {
    typename iterator_traits<_IIter>::difference_type val = 0;

    for( ; b != e; ++b)
      ++val;

    return val;
  }
}

namespace grid
{
  namespace save_mfolds
  {

  typedef dataset_t::mfold_t                                mfold_t;
  typedef boost::shared_ptr<mfold_t>                        mfold_ptr_t;

  struct work_item_t
  {
    int            cp_no; // cp no in the msc
    eGDIR          dir;   // which dir to track
    mfold_ptr_t    mfold;


    work_item_t(int _cp_no,eGDIR _dir)
      :cp_no(_cp_no),dir(_dir),mfold(new mfold_t){}

  };

  typedef std::vector<work_item_t>             work_list_t;
  typedef std::vector<mfold_t>                 mfold_list_t;
  typedef std::map<int,int>                    int_to_int_t;
  typedef mt_producer_t<work_list_t::iterator> wi_producer_t;

  template <typename Titer>  void build_work_list
    (mscomplex_const_ptr_t msc,work_list_t &work_list,Titer begin,Titer end)
  {
    std::set<int> tracked_cps[2];

    for (; begin != end;)
    {
      int i = *begin++;

      for(int d = 0 ; d < 2; ++d)
      {
        if(msc->m_rect.contains(msc->cellid(i)))
        {
          work_list.push_back(work_item_t(i,(eGDIR)d));
          tracked_cps[d].insert(i);
        }

        for( conn_iter_t j  = msc->m_conn[d][i].begin();
                         j != msc->m_conn[d][i].end();++j)
        {
          if(tracked_cps[d].count(msc->pair_idx(*j)) == 0)
          {
            tracked_cps[d].insert(msc->pair_idx(*j));
            work_list.push_back(work_item_t(msc->pair_idx(*j),(eGDIR)d));
          }
        }
      }
    }
  }


  void do_work_on_list
    (dataset_const_ptr_t ds,mscomplex_const_ptr_t msc,wi_producer_t &prd)
  {
    try
    {
      for(;;)
      {
        work_item_t wi = *prd.next();
        bfs::collect_manifolds(ds,wi.mfold.get(),msc->cellid(wi.cp_no),wi.dir);
      }
    }
    catch(wi_producer_t::no_more_items){} // not an error
  }

  mfold_ptr_t consolidate_mfold
    (mscomplex_const_ptr_t msc,int i,eGDIR d,
     work_list_t &wl,int_list_t cpno_to_wino[])
  {
    ASSERT(msc->is_saddle(i)&&(!msc->is_paired(i)));
    ASSERT(cpno_to_wino[d][i] != -1);

    work_item_t iwi = wl[cpno_to_wino[d][i]];

    for( conn_iter_t j  = msc->m_conn[d][i].begin();
                     j != msc->m_conn[d][i].end();++j)
    {
      ASSERT(msc->index(i) == msc->index(msc->pair_idx(*j)));
      ASSERT(cpno_to_wino[d][msc->pair_idx(*j)] != -1);

      work_item_t jwi = wl[cpno_to_wino[d][msc->pair_idx(*j)]];

      iwi.mfold->insert(iwi.mfold->end(),jwi.mfold->begin(),jwi.mfold->end());
    }

    return iwi.mfold;
  }

  int get_header_size(int num_cps)
  {
    return sizeof(rect_t)*3         + // rects
           sizeof(int)              + // num_cps
           sizeof(cellid_t)*num_cps + // cellids
           sizeof(int)*(2*num_cps+1); // offsets
  }


  template<typename T>
  void bin_write(std::ostream & os,const T & d)
  {
    os.write((const char*)(const void*)&d,sizeof(T));
  }

  template <typename Titer>
  void write_header(dataset_const_ptr_t ds,
                    mscomplex_const_ptr_t msc,
                    const int_list_t & offsets,
                    std::ostream & os,
                    Titer begin,Titer end)
  {
    os.seekp(0,ios::beg);

    bin_write(os,ds->m_rect);
    bin_write(os,ds->m_ext_rect);
    bin_write(os,ds->m_domain_rect);

    bin_write(os,std::count(begin,end));

    for(; begin != end;)
      bin_write(os,msc->cellid(*begin++));

    os.write((char*)(void*)offsets.data(),offsets.size()*sizeof(int));
  }

  void save_saddles(std::ostream & os,
            dataset_const_ptr_t ds,
            mscomplex_const_ptr_t msc)
  {
    // figure out which cps we need to track and put them in a work list
    work_list_t  work_list;

    mscomplex_t::filter_iterator_t cp_begin = msc->begin_unpaired_saddle();
    mscomplex_t::filter_iterator_t cp_end   = msc->end_unpaired_saddle();
    mscomplex_t::filter_iterator_t cp_it    = cp_begin;

    int num_cps   = std::count(cp_begin,cp_end);

    build_work_list(msc,work_list,cp_begin,cp_end);

    // advance os by the required space for header

    os.seekp(get_header_size(num_cps),ios::beg);

    // launch a bunch of threads to work on the list
    {
      boost::thread_group group;

      wi_producer_t prd(work_list.begin(),work_list.end());

      for(int tid = 0 ; tid < g_num_threads; ++tid)
        group.create_thread(bind(do_work_on_list,ds,msc,boost::ref(prd)));

      group.join_all();
    }

    // compile a map that maps cpno to workitem (-1 if invalid)

    int_list_t cpno_to_wino[2];

    cpno_to_wino[GDIR_DES].resize(msc->get_num_critpts(),-1);
    cpno_to_wino[GDIR_ASC].resize(msc->get_num_critpts(),-1);

    for(int i = 0 ; i < work_list.size(); ++i)
    {
      work_item_t wi = work_list[i];

      cpno_to_wino[wi.dir][wi.cp_no] = i;
    }

    // consolidate the manifold and write to os

    int_list_t offsets(num_cps*2+1);
    offsets[0] = 0;
    int pos = 0;

    for(cp_it = cp_begin ; cp_it != cp_end; ++cp_it)
    {
      for(int d = 0 ; d < 2; ++d)
      {
        mfold_ptr_t mfold = consolidate_mfold(msc,*cp_it,(eGDIR)d,work_list,cpno_to_wino);

        os.write((char*)(void*)mfold->data(),mfold->size()*sizeof(cellid_t));

        offsets[++pos] = offsets[pos-1]+mfold->size();

        mfold->clear();
      }
    }

    // write the actual header

    write_header(ds,msc,offsets,os,cp_begin,cp_end);
  }

  }

  void  dataset_t::saveManifolds(mscomplex_ptr_t msc,const std::string &bn)
  {
    std::ofstream fs((bn+".mfold.bin").c_str());
    ensure(fs.is_open(),"unable to open file");
    save_mfolds::save_saddles(fs,shared_from_this(),msc);
    fs.close();

#ifdef BUILD_EXEC_OPENCL
    opencl::update_to_surv_extrema(shared_from_this(),msc);
#else
    {
      cp_producer_ptr_t prd(new cp_producer_t(msc,cp_producer_t::extrema_filter));
      boost::thread_group group;
      for(int tid = 0 ; tid < g_num_threads; ++tid)
        group.create_thread(bind(dataset_t::mark_extrema_owner_thd,msc,prd));
      group.join_all();
    }
#endif

    {
      int ex_num = num_cells2(rect_t(m_rect.lc()+1,m_rect.uc()-1));
      std::ofstream fs((bn+".max.raw").c_str());
      ensure(fs.is_open(),"unable to open file");
      fs.write((char*)(void*)m_owner_maxima.data(),sizeof(int)*ex_num);
      fs.close();
    }

    {
      int ex_num = num_cells2(m_rect);
      std::ofstream fs((bn+".min.raw").c_str());
      ensure(fs.is_open(),"unable to open file");
      fs.write((char*)(void*)m_owner_minima.data(),sizeof(int)*ex_num);
      fs.close();
    }

  }

  void dataset_t::log_flags()
  {
    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
    {
      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
      {
        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
        {
          std::cout<<m_cell_flags(c)<<" ";
        }
        std::cout<<std::endl;
      }
      std::cout<<std::endl;
    }
  }

  char get_dir_txt(cellid_t c,cellid_t p)
  {
    int dir = 0;

    if (c[1] != p[1])
      dir = 1;
    else if( c[2] != p[2])
      dir = 2;

    if(dir == 2)
    {
      if( c[dir] > p[dir])
        return 'd';
      else
        return 'u';
    }

    if(c[dir]&1)
    {
      if(c[dir] > p[dir])
      {
        switch(dir)
        {
        case 0: return '>';
        case 1: return 'v';
        }
      }
      else
      {
        switch(dir)
        {
        case 0: return '<';
        case 1: return '^';
        }
      }
    }
    else
    {
      switch(dir)
      {
      case 0: return '-';
      case 1: return '|';
      }
    }

    return '#';
  }

  void dataset_t::log_pairs(std::ostream &os)
  {
    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
    {
      os<<"sheet no:: "<<c[2]<<std::endl;
      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
      {
        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
        {
          if(isCellCritical(c))
          {
            if(isCellPaired(c))
              os<<get_dir_txt(c,getCellPairId(c))<<"c";
            else
              os<<"C ";
          }
          else if(isCellPaired(c))
            os<<get_dir_txt(c,getCellPairId(c))<<" ";
          else
            os<<"? ";
        }
        os<<std::endl;
      }

    }
  }

  void dataset_t::log_visits(std::ostream &os)
  {
    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
    {
      os<<"sheet no:: "<<c[2]<<std::endl;
      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
      {
        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
        {
          if(isCellVisited(c))
            os<<"x ";
          else
            os<<". ";
        }
        os<<std::endl;
      }
    }
  }

  void dataset_t::log_pair_visits(std::ostream &os)
  {
    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
    {
      os<<"sheet no:: "<<c[2]<<std::endl;
      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
      {
        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
        {
          if(isCellVisited(c)&&isCellPaired(c)&&isCellVisited(getCellPairId(c)))
            os<<"x ";
          else
            os<<". ";
        }
        os<<std::endl;
      }
    }
  }

  void dataset_t::log_owner_extrema(eGDIR dir, std::ostream &os)
  {
    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    rect_t ex_rect               = (dir == GDIR_DES)?(rect_t(m_rect.lc()+1,m_rect.uc()-1)):(m_rect);
    owner_array_t &owner_extrema = (dir == GDIR_DES)?(m_owner_maxima):(m_owner_minima);

    int X = ex_rect[0].span()/2+1;
    int Y = ex_rect[1].span()/2+1;
    int N = num_cells2(ex_rect);

    for( int i = 0 ; i < N; ++i)
    {
      cellid_t c = i_to_c2(ex_rect,i);

      os<<owner_extrema(c/2)<<" ";

      if((i+1)%X == 0)
      {
        os<<endl;

        if(((i/X)+1)%Y == 0)
          os<<endl;
      }
    }

//    for(c[2] = ex_rect[2][0] ; c[2] <= ex_rect[2][1]; c[2]+=2)
//    {
//      os<<"sheet no:: "<<c[2]<<std::endl;
//      for(c[1] = ex_rect[1][0] ; c[1] <= ex_rect[1][1]; c[1]+=2)
//      {
//        for(c[0] = ex_rect[0][0] ; c[0] <= ex_rect[0][1]; c[0]+=2)
//        {
//          os<<m_owner_extrema[dir](c/2)<<" ";
//        }
//        os<<std::endl;
//      }
//    }
  }


  void dataset_t::log_pairs(const std::string &s)
  {
    std::ofstream fs(s.c_str());
    ensure(fs.is_open(),"unable to open file");
    log_pairs(fs);
    fs.close();
  }

  void dataset_t::log_pair_visits(const std::string &s)
  {
    std::ofstream fs(s.c_str());
    ensure(fs.is_open(),"unable to open file");
    log_pair_visits(fs);
    fs.close();
  }

  void dataset_t::log_visits(const std::string &s)
  {
    std::ofstream fs(s.c_str());
    ensure(fs.is_open(),"unable to open file");
    log_visits(fs);
    fs.close();
  }


  void dataset_t::log_max_facets()
  {
    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = m_rect[2][0] ; c[2] <= m_rect[2][1]; ++c[2])
    {
      for(c[1] = m_rect[1][0] ; c[1] <= m_rect[1][1]; ++c[1])
      {
        for(c[0] = m_rect[0][0] ; c[0] <= m_rect[0][1]; ++c[0])
        {
          if(getCellDim(c) != 0 )
            std::cout<< getCellMaxFacetId(c);
          else
            std::cout<< "(.,.,.)";
        }
        std::cout<<std::endl;
      }
      std::cout<<std::endl;
    }
  }

  void dataset_t::extract_vdata_subarray(rect_t r,const std::string &filename)
  {
    if(r.lower_corner()%2 != cellid_t::zero ||
       r.upper_corner()%2 != cellid_t::zero )
    {
      throw std::runtime_error("r must specify an aabb with vertex end pts");
    }

    std::ofstream ofs(filename.c_str(),std::ios::out|std::ios::binary);

    if(ofs.is_open() == false)
      throw std::runtime_error("unable to open file");

    static_assert(gc_grid_dim == 3 && "defined for 3-manifolds only");

    cellid_t c;

    for(c[2] = r[2][0] ; c[2] <= r[2][1]; c[2] +=2)
    {
      for(c[1] = r[1][0] ; c[1] <= r[1][1]; c[1] +=2)
      {
        for(c[0] = r[0][0] ; c[0] <= r[0][1]; c[0] +=2)
        {
          cell_fn_t fn = get_cell_fn(c);

          ofs.write((char*)(void*)&fn,sizeof(cell_fn_t));
        }
      }
    }
  }
}
