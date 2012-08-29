#include <stack>
#include <queue>
#include <list>
#include <set>
#include <fstream>
#include <tr1/tuple>

#include <boost/typeof/typeof.hpp>
#include <boost/function.hpp>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/static_assert.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>

#define static_assert BOOST_STATIC_ASSERT

#include <config.h>

#include <grid_dataset.h>
#include <grid_mscomplex.h>
#include <grid_dataset_bfs.h>

#ifdef BUILD_EXEC_OPENCL
#include <grid_dataset_cl.h>
#endif

#ifdef BUILD_EXEC_TBB

#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>

#endif

using namespace std;
namespace br = boost::range;
namespace ba = boost::adaptors;


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
      m_work_rect((r.lc()+e.lc())/2,(r.uc()+e.uc())/2),
      m_vert_fns(cellid_t::zero,boost::fortran_storage_order()),
      m_cell_flags(cellid_t::zero,boost::fortran_storage_order()),
      m_cell_order(cellid_t::zero,boost::fortran_storage_order()),
      m_owner_maxima(cellid_t::zero,boost::fortran_storage_order()),
      m_owner_minima(cellid_t::zero,boost::fortran_storage_order())

  {
    // TODO: assert that the given rect is of even size..
    //       since each vertex is in the even positions
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
    m_cell_order.resize(span);
    m_vert_fns.resize(pt_span);

    uint num_cells = span[0]*span[1]*span[2];
    uint num_pts   = pt_span[0]*pt_span[1]*pt_span[2];

    m_cell_flags.reindex(bl);
    m_cell_order.reindex(bl);
    m_vert_fns.reindex(bl/2);

    std::fill_n(m_cell_flags.data(),num_cells,0);

    ifstream ifs(filename.c_str(),ios::in|ios::binary);
    ensure(ifs.is_open(),"unable to open file");

    ifs.read((char*)(void*)m_vert_fns.data(),sizeof(cell_fn_t)*num_pts);
    ensure(ifs.fail()==false,"failed to read some data");

    ifs.seekg(0,ios::end);
    ensure(uint(ifs.tellg())==num_pts*sizeof(cell_fn_t),"file/piece size mismatch");

    ifs.close();

    m_owner_maxima.resize(m_rect.span()/2);
    m_owner_minima.resize((m_rect.span()/2)+1);

    m_owner_maxima.reindex(m_rect.lc()/2);
    m_owner_minima.reindex(m_rect.lc()/2);

  }

  void  dataset_t::clear()
  {
    m_cell_flags.resize(cellid_t::zero);
    m_cell_order.resize(cellid_t::zero);
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

    for(i[0] = -(l[0]) ; i[0] <= (l[0]) ;i[0]+=1)
    {
      for(i[1] = -(l[1]) ; i[1] <= (l[1]) ;i[1]+=1)
      {
        for(i[2] = -(l[2]) ; i[2] <= (l[2]) ;i[2]+=1)
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
//    ASSERT(isCellVisited(c) == false);
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
//    cellid_t f[20],c,s(0,0,0);

//    for(int i = 0 ; i < dim; ++i)
//      s[i] = 1;

//    while(true)
//    {
//      rect_t rect  = rect_t(m_ext_rect.lc()+s,m_ext_rect.uc()-s);

//      int N = num_cells2(rect);

//      for( int i = thid; i < N; i += g_num_threads)
//      {
//        c = i_to_c2(rect,i);

//        ASSERT(getCellDim(c) == dim);

//        int f_ct   = getCellFacets(c,f);

//        setCellMaxFacet(c,*std::max_element(f,f+f_ct,cmp_ftor));
//      }

//      if(!next_permutation(s.rbegin(),s.rend()))
//        break;
//    }
    ASSERT(false&&"impl broken here");
  }

  void  dataset_t::pairCellsWithinEst_thd(int tid, cellid_list_t *ccells)
  {
    ASSERT(false&&"impl broken here");
//    int N = num_cells2(m_rect);

//    for(int i = tid; i < N; i += g_num_threads)
//    {
//      cellid_t c = i_to_c2(m_rect,i);

//      cellid_t est_arr[40];

//      uint est_ct = getCellEst(c,est_arr);

//      std::sort(est_arr,est_arr+est_ct,cmp_ftor);

//      for(int j = 0; j < est_ct; ++j)
//      {
//        cellid_t p = est_arr[j];

//        if(j+1 < est_ct)
//        {
//          cellid_t q = est_arr[j+1];

//          bool is_adj        = areCellsIncident(p,q);
//          bool is_same_bndry = (isTrueBoundryCell(p) == isTrueBoundryCell(q));

//          if(is_adj && is_same_bndry)
//          {
//            pairCells(p,q);
//            ++j;
//            continue;
//          }
//        }

//        if(m_rect.contains(p))
//        {
//          markCellCritical(p);
//          ccells->push_back(p);
//        }
//      }
//    }
  }

  inline bool is_pairable2
  ( const dataset_t &ds,
    cellid_t p, cellid_t p_mf ,cellid_t q)
  {
    if(!ds.m_ext_rect.contains(q))
      return false;

    if(ds.m_domain_rect.isOnBoundry(p) != ds.m_domain_rect.isOnBoundry(q))
      return false;

    cellid_t q_mf    = ds.getCellMaxFacetId(q);
    cellid_t q_mf_mf = ds.getCellMaxFacetId(q_mf);

    return ((q_mf != p) && (q_mf_mf == p_mf));
  }

  cellid_t invalid_cell(-1,-1,-1);

  inline bool set_pair_edge2(dataset_t& ds, cellid_t e,cellid_t d0,cellid_t d1)
  {
    if(!ds.m_rect.contains(e))
      return false;

    if(ds.isCellPaired(e))
      return false;

    cellid_t e_mf =  ds.getCellMaxFacetId(e);

    cellid_t f0 = e - d0;
    cellid_t f1 = e + d0;
    cellid_t f2 = e - d1;
    cellid_t f3 = e + d1;

    cellid_t f = invalid_cell;

    if (is_pairable2(ds,e,e_mf,f0))
      f = f0;

    if ((is_pairable2(ds,e,e_mf,f1)) &&
        ((f == invalid_cell) || ds.compare_cells<2>(f1,f)))
      f = f1;

    if ((is_pairable2(ds,e,e_mf,f2)) &&
        ((f == invalid_cell) || ds.compare_cells<2>(f2,f)))
      f = f2;

    if ((is_pairable2(ds,e,e_mf,f3)) &&
        ((f == invalid_cell) || ds.compare_cells<2>(f3,f)))
      f = f3;

    if(f != invalid_cell && !ds.isCellPaired(f))
    {
      ds.pairCells(e,f);
      return true;
    }
    return false;
  }

  inline bool set_pair_face2
  ( dataset_t &ds,cellid_t f,cellid_t d)
  {
    if(!ds.m_rect.contains(f))
      return false;

    if(ds.isCellPaired(f))
      return false;

    cellid_t c0 = f - d;
    cellid_t c1 = f + d;

    cellid_t f_mf =  ds.getCellMaxFacetId(f);

    cellid_t c = invalid_cell;

    if (is_pairable2(ds,f,f_mf,c0))
      c = c0;

    if ((is_pairable2(ds,f,f_mf,c1)) &&
        ((c == invalid_cell) || ds.compare_cells<3>(c1,c)))
      c = c1;

    if(!(c == invalid_cell) && !ds.isCellPaired(c))
    {
      ds.pairCells(f,c);
      return true;
    }
    return false;
  }


  inline bool is_pairable3
  ( const dataset_t &ds,
    cellid_t p, cellid_t p_mf ,cellid_t p_mf_mf ,cellid_t q)
  {
    if(!ds.m_ext_rect.contains(q))
      return false;

    if(ds.m_domain_rect.isOnBoundry(p) != ds.m_domain_rect.isOnBoundry(q))
      return false;

    cellid_t q_mf       = ds.getCellMaxFacetId(q);
    cellid_t q_mf_mf    = ds.getCellMaxFacetId(q_mf);
    cellid_t q_mf_mf_mf = ds.getCellMaxFacetId(q_mf_mf);

    return ((q_mf != p) && (q_mf_mf != p_mf) && (q_mf_mf_mf == p_mf_mf));
  }

  inline bool set_pair_face3
  ( dataset_t &ds,cellid_t f,cellid_t d)
  {
    if(!ds.m_rect.contains(f))
      return false;

    if(ds.isCellPaired(f))
      return false;

    cellid_t c0 = f - d;
    cellid_t c1 = f + d;

    cellid_t f_mf    =  ds.getCellMaxFacetId(f);
    cellid_t f_mf_mf =  ds.getCellMaxFacetId(f_mf);

    cellid_t c = invalid_cell;

    if (is_pairable3(ds,f,f_mf,f_mf_mf,c0))
      c = c0;

    if ((is_pairable3(ds,f,f_mf,f_mf_mf,c1)) &&
        ((c == invalid_cell) || ds.compare_cells<3>(c1,c)))
      c = c1;

    if(!(c == invalid_cell) && !ds.isCellPaired(c))
    {
      ds.pairCells(f,c);
      return true;
    }
    return false;
  }

  void  dataset_t::assign_pairs2()
  {
    cellid_t X = cellid_t(1,0,0);
    cellid_t Y = cellid_t(0,1,0);
    cellid_t Z = cellid_t(0,0,1);

    cellid_t XY = cellid_t(1,1,0);
    cellid_t YZ = cellid_t(0,1,1);
    cellid_t ZX = cellid_t(1,0,1);

    int n_p = 0;

    for(iterator_dim b = begin(0),e = end(0); b!= e; ++b)
    {
      cellid_t v = *b;

      if (set_pair_edge2(*this,v+X,Y,Z)) ++n_p;
      if (set_pair_edge2(*this,v+Y,Z,X)) ++n_p;
      if (set_pair_edge2(*this,v+Z,X,Y)) ++n_p;

      if (set_pair_face2(*this,v+XY,Z)) ++n_p;
      if (set_pair_face2(*this,v+YZ,X)) ++n_p;
      if (set_pair_face2(*this,v+ZX,Y)) ++n_p;
    }

    cout<<SVAR(n_p)<<endl;
    n_p = 0;

    for(iterator_dim b = begin(0),e = end(0); b!= e; ++b)
    {
      cellid_t v = *b;

      if (set_pair_face3(*this,v+XY,Z)) ++n_p;
      if (set_pair_face3(*this,v+YZ,X)) ++n_p;
      if (set_pair_face3(*this,v+ZX,Y)) ++n_p;
    }

    cout<<SVAR(n_p)<<endl;
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
    for(int i = 0 ; i < (int)ccells->size();++i)
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

  template<int dim,eGDIR dir>
  inline bool is_required_cp(const mscomplex_t& msc,int i)
  {
    static const int pdim = (dir == GDIR_DES)? (dim - 1):(dim +1);

    return msc.index(i) == dim && msc.m_rect.contains(msc.cellid(i))
        && (!msc.is_paired(i) || msc.index(msc.pair_idx(i)) == pdim);
  }

  template<typename T>
  class producer_consumer_t:boost::noncopyable
  {

  public:
    std::queue<T>        m_queue;
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

  typedef producer_consumer_t<cellid_list_ptr_t> cp_conn_que_t;

  void  compute_saddle_connections( dataset_t &ds,cp_id_producer_t &prd,cp_conn_que_t &que)
  {
    for(BOOST_AUTO(it,prd.next()); prd.is_valid(it); it = prd.next())
    {
      que.put(compute_inc_pairs<2>(*it,ds));
    }
  }

  template<int dim,typename T>
  void extract_extrema_connections(T b,T e,dataset_t &ds,cp_conn_que_t &que)
  {
    const eGDIR sad_dir = (dim == 1)?(GDIR_DES):(GDIR_ASC);
    const eGDIR ex_dir  = (dim == 1)?(GDIR_ASC):(GDIR_DES);

    cellid_list_ptr_t clist(new cellid_list_t);

    for(;b != e; ++b)
    {
      cellid_t c = *b,e1,e2;

      get_adj_extrema<sad_dir>(c,e1,e2);

      if(ds.m_rect.contains(e1))
      {
        e1 = i_to_c2(ds.get_extrema_rect<ex_dir>(),ds.owner_extrema<ex_dir>()(e1/2));
        clist->push_back(c);
        clist->push_back(e1);
      }

      if(ds.m_rect.contains(e2))
      {
        e2 = i_to_c2(ds.get_extrema_rect<ex_dir>(),ds.owner_extrema<ex_dir>()(e2/2));
        clist->push_back(c);
        clist->push_back(e2);
      }
    }

    que.put(clist);
  }

  void  store_connections( mscomplex_t &msc,cp_conn_que_t &que,int n)
  {
    map<cellid_t,int> id_cp_map;

    for( int i = 0 ; i < msc.get_num_critpts(); ++i)
      id_cp_map[msc.cellid(i)] = i;


    while( n-- > 0)
    {
      cellid_list_ptr_t ptr =  que.get();
      cellid_list_t &inc_pairs=*ptr;

      for(int i = 0 ; i < (int)inc_pairs.size(); i +=2)
      {
        cellid_t u = inc_pairs[i],v = inc_pairs[i+1];

        ASSERT(id_cp_map.count(u) ==1);
        ASSERT(id_cp_map.count(v) ==1);

        msc.connect_cps(id_cp_map[u],id_cp_map[v]);
      }
    }
  }

  void check_order(const dataset_t & ds)
  {
    cellid_t f[40];

    for(dataset_t::iterator b = ds.begin(),e = ds.end(); b !=e ; ++b)
    {
      cellid_t c= *b;

      int n_gtr = 0;

      for(cellid_t *f_b = f,*f_e = f+ ds.getCellFacets(c,f); f_b != f_e; ++f_b)
      {
        if(ds.compare_cells<-1>(c,*f_b))
          n_gtr++;
      }

      try{ASSERT(n_gtr <= 1);}catch(assertion_error e)
      {
        cout<<SVAR(c)<<endl;
        cout<<SVAR(ds.get_cell_vert(c))<<endl;
        for(cellid_t *f_b = f,*f_e = f+ ds.getCellEst(ds.get_cell_vert(c),f); f_b != f_e; ++f_b)
        {
          cout<<*f_b<<"\t"<<(int)ds.m_cell_order(*f_b);

          if(ds.compare_cells<-1>(c,*f_b))
            cout<<"*";

          cout<<"\t"<<get_cell_dim(*f_b);

          if(get_cell_dim(*f_b) != 0 )
            cout<<"\t"<<ds.getCellMaxFacetId(*f_b);

          cout<<endl;

        }
        throw;
      }

      int n_lsr = 0;

      for(cellid_t *f_b = f,*f_e = f+ ds.getCellCofacets(c,f); f_b != f_e; ++f_b)
      {
        if(ds.compare_cells<-1>(*f_b,c))
          n_lsr++;
      }

      ASSERT(n_lsr <= 1);

      if(ds.isCellPaired(c))
      {
        cellid_t p = ds.getCellPairId(c);

        if(get_cell_dim(c) < get_cell_dim(p))
        {
          ASSERT(ds.m_cell_order(c) > ds.m_cell_order(p));
        }
      }
    }
  }

#ifdef BUILD_EXEC_TBB

  struct compute_saddle_connections_tbb_ftor
  {
    cellid_list_t  &m_saddles;
    dataset_t      &m_ds;
    mscomplex_t    &m_msc;
    cp_conn_que_t  &m_que;

    void operator()( const tbb::blocked_range<size_t>& r) const
    {
      for( size_t i=r.begin(); i!=r.end(); ++i )
        m_que.put(compute_inc_pairs<2>(m_saddles[i],m_ds));
    }

    compute_saddle_connections_tbb_ftor
      (cellid_list_t & s,dataset_t &ds,mscomplex_t &msc,cp_conn_que_t &que):
      m_saddles(s),m_ds(ds),m_msc(msc),m_que(que){}
  };

  void compute_saddle_connections_tbb(dataset_t &ds, mscomplex_t &msc,cp_conn_que_t &que)
  {
    cellid_list_t saddles;

    br::copy(msc.cpno_range()
             |ba::filtered(bind(is_required_cp<2,DES>,boost::cref(msc),_1))
             |ba::transformed(bind(&mscomplex_t::cellid,boost::cref(msc),_1)),
             back_inserter(saddles));

    tbb::parallel_for(tbb::blocked_range<size_t>(0,saddles.size()),
                      compute_saddle_connections_tbb_ftor(saddles,ds,msc,que),
                      tbb::auto_partitioner());
  }
#endif


  void  dataset_t::computeMsGraph(mscomplex_ptr_t msc)
  {
#ifdef BUILD_EXEC_OPENCL
    opencl::worker w;
    w.assign_gradient(shared_from_this(),msc);
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

    mscomplex_t::filter_t f_2des = bind(is_required_cp<2,GDIR_DES>,boost::cref(*msc),_1);

    mscomplex_t::filter_t f_2asc = bind(is_required_cp<2,GDIR_ASC>,boost::cref(*msc),_1);
    mscomplex_t::filter_t f_1des = bind(is_required_cp<1,GDIR_DES>,boost::cref(*msc),_1);

    mark_reachable<1,GDIR_ASC>
        (msc->cpno_range()
         |ba::filtered(bind(is_required_cp<1,GDIR_ASC>,boost::cref(*msc),_1))
         |ba::transformed(bind(&mscomplex_t::cellid,boost::cref(msc),_1)),*this);

    {
      cp_id_producer_t prd(msc,f_2des);

      cp_conn_que_t que;

      boost::thread_group group;

      group.create_thread(bind(store_connections,boost::ref(*msc),boost::ref(que),prd.count()+2));

#ifdef BUILD_EXEC_TBB
      using boost::ref;
      group.create_thread(bind(compute_saddle_connections_tbb,ref(*this),ref(*msc),ref(que)));

#else
      for(int tid = 0 ; tid < g_num_threads; ++tid)
        group.create_thread(bind(compute_saddle_connections,boost::ref(*this),boost::ref(prd),boost::ref(que)));
      //      int n = prd.count();
      //      compute_saddle_connections(*this,prd,que);
#endif


#ifdef BUILD_EXEC_OPENCL
      w.owner_extrema(shared_from_this());
#else
      #error "please write code"
#endif

      extract_extrema_connections<2>
          (msc->cp_id_fbegin(f_2asc),msc->cp_id_fend(f_2asc),*this,que);

      extract_extrema_connections<1>
          (msc->cp_id_fbegin(f_1des),msc->cp_id_fend(f_1des),*this,que);

//      store_connections(*msc,que,n+2);

      group.join_all();
    }
  }

  void dataset_t::compute_owner_grad()
  {
#ifdef BUILD_EXEC_OPENCL
    opencl::assign_gradient_and_owner_extrema(shared_from_this());
#endif
  }

  typedef tr1::tuple<cellid_t,cellid_list_ptr_t,cellid_list_ptr_t> cp_mfold_qitem_t;

  typedef producer_consumer_t<cp_mfold_qitem_t> cp_mfold_que_t;

  inline cellid_t cellid_of_pair(mscomplex_t &msc,int i)
  {return msc.cellid(msc.pair_idx(i));}

  template<eGDIR dir>
  inline void collect_contrib_cps(mscomplex_t &msc,cellid_list_t &cplist, int i)
  {
    BOOST_AUTO(cop_ftr,bind(cellid_of_pair,boost::ref(msc),_1));
    BOOST_AUTO(in_rect_ftr,bind(&rect_t::contains_point,&msc.m_rect,_1));

    if(msc.m_rect.contains(msc.cellid(i)))
      cplist.push_back(msc.cellid(i));

    br::copy(boost::make_iterator_range(msc.m_conn[dir][i])
             |ba::map_keys|ba::transformed(cop_ftr)|ba::filtered(in_rect_ftr),
             back_inserter(cplist));

  }

  template<int dim,eGDIR dir,typename cmp_t>
  inline cellid_list_ptr_t compute_mfold(dataset_t &ds,mscomplex_t &msc,int i,cmp_t cmp)
  {
    cellid_list_ptr_t mfold(new cellid_list_t);

    cellid_list_t cplist;

    collect_contrib_cps<dir>(msc,cplist,i);

    compute_mfold<dim,dir,false>(cplist.begin(),cplist.end(),ds,*mfold,cmp);

    return mfold;
  }

  void compute_saddle_manifold(dataset_t &ds, mscomplex_t &msc, cp_producer_t &prd,cp_mfold_que_t &que)
  {
    BOOST_AUTO(des2_cmp,boost::bind(&dataset_t::compare_cells<2>,&ds,_1,_2));
    BOOST_AUTO(asc2_cmp,boost::bind(&dataset_t::compare_cells<2>,&ds,_2,_1));
    BOOST_AUTO(des1_cmp,boost::bind(&dataset_t::compare_cells<1>,&ds,_1,_2));
    BOOST_AUTO(asc1_cmp,boost::bind(&dataset_t::compare_cells<1>,&ds,_2,_1));


    for(BOOST_AUTO(it,prd.next()); prd.is_valid(it); it = prd.next())
    {
      int i = *it,dim = msc.index(i);

      cellid_list_ptr_t des,asc;

      if(dim == 1)
      {
        des = compute_mfold<1,GDIR_DES>(ds,msc,i,des1_cmp);
        asc = compute_mfold<1,GDIR_ASC>(ds,msc,i,asc1_cmp);
      }
      else
      {
        des = compute_mfold<2,GDIR_DES>(ds,msc,i,des2_cmp);
        asc = compute_mfold<2,GDIR_ASC>(ds,msc,i,asc2_cmp);
      }

      que.put(tr1::make_tuple(msc.cellid(i),des,asc));
    }
  }

  void store_mfold(std::ostream &os,cellid_list_t &cps,int_list_t &offsets, cp_mfold_que_t &que,int n)
  {
    int offset = 0;

    offsets.push_back(offset);

    while(n-- > 0)
    {
      cp_mfold_qitem_t qitem = que.get();

      BOOST_AUTO(c,tr1::get<0>(qitem));
      BOOST_AUTO(des,tr1::get<1>(qitem));
      BOOST_AUTO(asc,tr1::get<2>(qitem));

      offset += des->size(); offsets.push_back(offset);
      offset += asc->size(); offsets.push_back(offset);

      os.write((char*)(void*)des->data(),des->size()*sizeof(cellid_t));
      os.write((char*)(void*)asc->data(),asc->size()*sizeof(cellid_t));

      cps.push_back(c);
    }
  }

  int get_header_size(int num_cps)
  {
    return sizeof(rect_t)*3         + // rects
           sizeof(int)              + // num_cps
           sizeof(cellid_t)*num_cps + // cellids
           sizeof(int)*(2*num_cps+1); // offsets
  }


  template<typename T>
  inline void bin_write(std::ostream & os,const T & d)
  {os.write((const char*)(const void*)&d,sizeof(T));}

  template<typename T>
  inline void bin_read(std::istream &is, const T &v)
  {is.read((char*)(void*)&v,sizeof(T));}

  template<typename T>
  inline void bin_write_vec(std::ostream &os, std::vector<T> &v)
  {os.write((const char*)(const void*)v.data(),v.size()*sizeof(T));}


  void write_header(std::ostream & os, dataset_t &ds,int_list_t & offsets,cellid_list_t &cps)
  {
    os.seekp(0,ios::beg);

    bin_write(os,ds.m_rect);
    bin_write(os,ds.m_ext_rect);
    bin_write(os,ds.m_domain_rect);

    bin_write(os,(int)cps.size());

    os.write((char*)(void*)cps.data(),cps.size()*sizeof(cellid_t));
    os.write((char*)(void*)offsets.data(),offsets.size()*sizeof(int));
  }

  inline bool need_saddle_mfold(const mscomplex_t& msc,int i)
  {
    return msc.is_unpaired_saddle(i);
  }

  void save_saddle_mfolds(std::ostream &os,dataset_t &ds,mscomplex_t &msc)
  {
    mscomplex_t::filter_t fltr = bind(need_saddle_mfold,boost::cref(msc),_1);

    cp_producer_t prd(msc.shared_from_this(),fltr);

    int num_cps = prd.count();

    cp_mfold_que_t que;

    boost::thread_group group;

    for(int tid = 0 ; tid < g_num_threads; ++tid)
      group.create_thread(bind(compute_saddle_manifold,boost::ref(ds),
                               boost::ref(msc),boost::ref(prd),boost::ref(que)));
//      compute_saddle_manifold(ds,msc,prd,que);

    os.seekp(get_header_size(num_cps),ios::beg);

    cellid_list_t cps;
    int_list_t    offsets;

    store_mfold(os,cps,offsets,que,num_cps);

    write_header(os,ds,offsets,cps);

    group.join_all();//redundant
  }

  void  dataset_t::saveManifolds(mscomplex_ptr_t msc,const std::string &bn)
  {
    std::ofstream fs((bn+".mfold.bin").c_str());
    ensure(fs.is_open(),"unable to open file");
    boost::thread_group group;
    group.create_thread(bind(save_saddle_mfolds,boost::ref(fs),boost::ref(*this),boost::ref(*msc)));

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

    group.join_all();
  }

  template<int dim,eGDIR dir,typename cmp_t,typename Titer>
  inline cellid_list_ptr_t compute_mfold(dataset_t &ds,mscomplex_t &msc,Titer b,Titer e,cmp_t cmp)
  {
    using boost::ref;

    cellid_list_ptr_t mfold(new cellid_list_t);

    cellid_list_t cplist;

    for_each(b,e,bind(collect_contrib_cps<dir>,ref(msc),ref(cplist),_1));

    const bool trp = (dim==1 && dir==GDIR_ASC)||(dim==2 && dir==GDIR_DES);

    compute_mfold<dim,dir,trp>(cplist.begin(),cplist.end(),ds,*mfold,cmp);

    return mfold;
  }


  void dataset_t::saveConnectingOneManifolds(mscomplex_ptr_t msc,const std::string &bn)
  {
    mscomplex_t::memb_filter_t fltr_1 = &mscomplex_t::is_unpaired_one_saddle;
    mscomplex_t::memb_filter_t fltr_2 = &mscomplex_t::is_unpaired_two_saddle;

    mscomplex_t::fiterator_t b_1 = msc->fbegin(fltr_1),e_1 = msc->fend(fltr_1);
    mscomplex_t::fiterator_t b_2 = msc->fbegin(fltr_2),e_2 = msc->fend(fltr_2);

    BOOST_AUTO(des2_cmp,boost::bind(&dataset_t::compare_cells<2>,this,_1,_2));
    BOOST_AUTO(asc2_cmp,boost::bind(&dataset_t::compare_cells<2>,this,_2,_1));
    BOOST_AUTO(des1_cmp,boost::bind(&dataset_t::compare_cells<1>,this,_1,_2));
    BOOST_AUTO(asc1_cmp,boost::bind(&dataset_t::compare_cells<1>,this,_2,_1));

    cellid_list_t cplist;

    for_each(b_2,e_2,bind(collect_contrib_cps<GDIR_DES>,ref(*msc),ref(cplist),_1));

    mark_reachable<2,GDIR_DES>(cplist,*this);

    cellid_list_ptr_t one_des = compute_mfold<1,GDIR_DES>(*this,*msc,b_1,e_1,des1_cmp);
    cellid_list_ptr_t one_asc = compute_mfold<1,GDIR_ASC>(*this,*msc,b_1,e_1,asc1_cmp);
    cellid_list_ptr_t two_asc = compute_mfold<2,GDIR_ASC>(*this,*msc,b_2,e_2,asc2_cmp);

    std::ofstream fs((bn+".onemfold.bin").c_str());
    ensure(fs.is_open(),"unable to open file");

    bin_write(fs,int(one_des->size()));
    bin_write(fs,int(one_asc->size()));
    bin_write(fs,int(0));
    bin_write(fs,int(two_asc->size()));

    bin_write_vec(fs,*one_des);
    bin_write_vec(fs,*one_asc);
    bin_write_vec(fs,*two_asc);
  }

  template<typename T>
  inline void bin_write_marray(std::ostream & os,const T * p, const cellid_t & s)
  {os.write((const char*)(const void*)p,sizeof(T)*s[0]*s[1]*s[2]);}

  template<typename T>
  inline void bin_read_marray(std::istream & is,T * p, const cellid_t & s)
  {is.read((char*)(void*)p,sizeof(T)*s[0]*s[1]*s[2]);}


  void dataset_t::store(std::ostream &os)
  {
    bin_write(os,m_rect);
    bin_write(os,m_ext_rect);
    bin_write(os,m_domain_rect);

    bin_write_marray(os,m_cell_flags.data(),m_ext_rect.span()+1);
    bin_write_marray(os,m_cell_order.data(),m_ext_rect.span()+1);
    bin_write_marray(os,m_vert_fns.data(),(m_ext_rect.span()/2)+1);
    bin_write_marray(os,m_owner_maxima.data(),m_rect.span()/2);
    bin_write_marray(os,m_owner_minima.data(),(m_rect.span()/2)+1);
  }

  void dataset_t::load(std::istream &is)
  {
    bin_read(is,m_rect);
    bin_read(is,m_ext_rect);
    bin_read(is,m_domain_rect);

    m_cell_flags.resize(m_ext_rect.span()+1);
    m_cell_order.resize(m_ext_rect.span()+1);
    m_vert_fns.resize((m_ext_rect.span()/2)+1);
    m_owner_maxima.resize(m_rect.span()/2);
    m_owner_minima.resize((m_rect.span()/2)+1);

    m_cell_flags.reindex(m_ext_rect.lc());
    m_cell_order.reindex(m_ext_rect.lc());
    m_vert_fns.reindex(m_ext_rect.lc()/2);
    m_owner_maxima.reindex(m_rect.lc()/2);
    m_owner_minima.reindex(m_rect.lc()/2);

    bin_read_marray(is,m_cell_flags.data(),m_ext_rect.span()+1);
    bin_read_marray(is,m_cell_order.data(),m_ext_rect.span()+1);
    bin_read_marray(is,m_vert_fns.data(),(m_ext_rect.span()/2)+1);
    bin_read_marray(is,m_owner_maxima.data(),(m_rect.span()/2));
    bin_read_marray(is,m_owner_minima.data(),(m_rect.span()/2)+1);
  }

  template<eGDIR dir>
  void copy_owner_extrema(int_marray_t &dest,const int_marray_t &src,const dataset_t & ds)
  {
    rect_t ex_rct = ds.get_extrema_rect<dir>();

    int num_ex = num_cells2(ex_rct);

    dest.resize(ex_rct.span()/2+1);

    memcpy((void*)dest.data(),(void*)src.data(),num_ex*sizeof(int));
  }

  void  dataset_t::storeOwnerArrays(int_marray_t &omax,int_marray_t &omin) const
  {
    copy_owner_extrema<GDIR_DES>(omax,m_owner_maxima,*this);
    copy_owner_extrema<GDIR_ASC>(omin,m_owner_minima,*this);
  }

  void  dataset_t::loadOwnerArrays(int_marray_t &omax,int_marray_t &omin)
  {
    copy_owner_extrema<GDIR_DES>(m_owner_maxima,omax,*this);
    copy_owner_extrema<GDIR_ASC>(m_owner_minima,omin,*this);
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
    int_marray_t &owner_extrema  = (dir == GDIR_DES)?(m_owner_maxima):(m_owner_minima);

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
