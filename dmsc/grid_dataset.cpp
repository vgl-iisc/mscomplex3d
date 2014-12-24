#include <stack>
#include <queue>
#include <list>
#include <set>
#include <fstream>
#include <tr1/tuple>

#include <boost/typeof/typeof.hpp>
#include <boost/function.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/static_assert.hpp>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/iterator/transform_iterator.hpp>

#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>

#include <omp.h>

#define static_assert BOOST_STATIC_ASSERT

#include <config.h>

#include <grid_dataset.h>
#include <grid_mscomplex.h>
#include <grid_dataset_bfs.h>

#include <grid_dataset_cl.h>

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

  template<int dim,eGDIR dir>
  inline bool is_required_cp(const mscomplex_t& msc,int i)
  {
    static const int pdim = (dir == GDIR_DES)? (dim - 1):(dim +1);

    return msc.index(i) == dim && msc.m_rect.contains(msc.cellid(i))
        && (!msc.is_paired(i) || msc.index(msc.pair_idx(i)) == pdim);
  }

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

      #pragma omp critical(mscomplex_connector_critical_section)
      {
//        if (id_cp_map.size() == 0)
//          init();

        msc->connect_cps(i,j,m);
      }
    }
  };


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

  template <int dim,eGDIR dir>
  inline bool cellid_int_pair_cmp(dataset_ptr_t ds, cellid_int_pair_t c1,cellid_int_pair_t c2)
  {return ds->compare_cells_pp_<dir,dim>(c1.first,c2.first);}

  template <int dim,eGDIR dir>
  inline void compute_inc_pairs_pq
  (cellid_t s, dataset_ptr_t ds, mscomplex_connector_t &msc_connector)
  {
    const int pdim = (dir == DES)?(dim - 1):(dim + 1);
    const bool no_vcheck =
        ((dim==1) && (dir==DES)) || (dim==0) ||
        ((dim==2) && (dir==ASC)) || (dim==3);

    BOOST_AUTO(cmp_dim , bind(cellid_int_pair_cmp< dim, dir>,ds,_1,_2));
    BOOST_AUTO(cmp_pdim, bind(cellid_int_pair_cmp<pdim, dir>,ds,_1,_2));

    priority_queue<cellid_int_pair_t,cellid_int_pair_list_t,typeof(cmp_dim) >     pq(cmp_dim );
    priority_queue<cellid_int_pair_t,cellid_int_pair_list_t,typeof(cmp_pdim)> inc_pq(cmp_pdim);

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

          if (p != c && (no_vcheck||ds->isCellVisited(*b)) && ds->getCellDim(p) == dim )
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

  void computeSaddleConnections(mscomplex_ptr_t msc,dataset_ptr_t ds,
                                mscomplex_connector_t &msc_connector)
  {
    BOOST_AUTO(cps_1asc,msc->cpno_range()|
               ba::filtered(bind(is_required_cp<1,GDIR_ASC>,boost::cref(*msc),_1))|
               ba::transformed(bind(&mscomplex_t::cellid,boost::cref(msc),_1)));

    mark_reachable<1,GDIR_ASC,typeof(cps_1asc)>(cps_1asc,ds);


    cellid_list_t cps_2des;

    br::copy(msc->cpno_range()|
             ba::filtered(bind(is_required_cp<2,GDIR_DES>,boost::cref(*msc),_1))|
             ba::transformed(bind(&mscomplex_t::cellid,boost::cref(msc),_1)),
             std::back_inserter(cps_2des));

    #pragma omp parallel for
    for(int i = 0 ; i < cps_2des.size(); ++i)
    {
      compute_inc_pairs_pq<2,DES>(cps_2des[i],ds,msc_connector);
    }
  }

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



  void  dataset_t::computeMsGraph(mscomplex_ptr_t msc)
  {
    opencl::worker w;
    w.assign_gradient(shared_from_this(),msc);

    mscomplex_connector_t msc_connector(msc);

    #pragma omp sections
    {
      #pragma omp section
      {msc_connector.init();}

      #pragma omp section
      {computeSaddleConnections(msc,shared_from_this(),msc_connector);}

      #pragma omp section
      {
        w.owner_extrema(shared_from_this());
        computeExtremaConnections<DES>(msc,shared_from_this(),msc_connector);
        computeExtremaConnections<ASC>(msc,shared_from_this(),msc_connector);
      }
    }
  }

  void dataset_t::compute_owner_grad()
  {
    opencl::assign_gradient_and_owner_extrema(shared_from_this());
  }

//  typedef tr1::tuple<cellid_t,cellid_list_ptr_t,cellid_list_ptr_t> cp_mfold_qitem_t;

//  typedef producer_consumer_t<cp_mfold_qitem_t> cp_mfold_que_t;

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

//  void compute_saddle_manifold(dataset_t &ds, mscomplex_t &msc, cp_producer_t &prd,cp_mfold_que_t &que)
//  {
//    BOOST_AUTO(des2_cmp,boost::bind(&dataset_t::compare_cells_pp<2>,&ds,_1,_2));
//    BOOST_AUTO(asc2_cmp,boost::bind(&dataset_t::compare_cells_pp<2>,&ds,_2,_1));
//    BOOST_AUTO(des1_cmp,boost::bind(&dataset_t::compare_cells_pp<1>,&ds,_1,_2));
//    BOOST_AUTO(asc1_cmp,boost::bind(&dataset_t::compare_cells_pp<1>,&ds,_2,_1));


//    for(BOOST_AUTO(it,prd.next()); prd.is_valid(it); it = prd.next())
//    {
//      int i = *it,dim = msc.index(i);

//      cellid_list_ptr_t des,asc;

//      if(dim == 1)
//      {
//        des = compute_mfold<1,GDIR_DES>(ds,msc,i,des1_cmp);
//        asc = compute_mfold<1,GDIR_ASC>(ds,msc,i,asc1_cmp);
//      }
//      else
//      {
//        des = compute_mfold<2,GDIR_DES>(ds,msc,i,des2_cmp);
//        asc = compute_mfold<2,GDIR_ASC>(ds,msc,i,asc2_cmp);
//      }

//      que.put(tr1::make_tuple(msc.cellid(i),des,asc));
//    }
//  }

//  void store_mfold(std::ostream &os,cellid_list_t &cps,int_list_t &offsets, cp_mfold_que_t &que,int n)
//  {
//    int offset = 0;

//    offsets.push_back(offset);

//    while(n-- > 0)
//    {
//      cp_mfold_qitem_t qitem = que.get();

//      BOOST_AUTO(c,tr1::get<0>(qitem));
//      BOOST_AUTO(des,tr1::get<1>(qitem));
//      BOOST_AUTO(asc,tr1::get<2>(qitem));

//      offset += des->size(); offsets.push_back(offset);
//      offset += asc->size(); offsets.push_back(offset);

//      os.write((char*)(void*)des->data(),des->size()*sizeof(cellid_t));
//      os.write((char*)(void*)asc->data(),asc->size()*sizeof(cellid_t));

//      cps.push_back(c);
//    }
//  }

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

//  void save_saddle_mfolds(std::ostream &os,dataset_t &ds,mscomplex_t &msc)
//  {
//    mscomplex_t::filter_t fltr = bind(need_saddle_mfold,boost::cref(msc),_1);

//    cp_producer_t prd(msc.shared_from_this(),fltr);

//    int num_cps = prd.count();

//    cp_mfold_que_t que;

//    boost::thread_group group;

//    for(int tid = 0 ; tid < g_num_threads; ++tid)
//      group.create_thread(bind(compute_saddle_manifold,boost::ref(ds),
//                               boost::ref(msc),boost::ref(prd),boost::ref(que)));
////      compute_saddle_manifold(ds,msc,prd,que);

//    os.seekp(get_header_size(num_cps),ios::beg);

//    cellid_list_t cps;
//    int_list_t    offsets;

//    store_mfold(os,cps,offsets,que,num_cps);

//    write_header(os,ds,offsets,cps);

//    group.join_all();//redundant
//  }

  void  dataset_t::saveManifolds(mscomplex_ptr_t msc,const std::string &bn)
  {
//    std::ofstream fs((bn+".mfold.bin").c_str());
//    ensure(fs.is_open(),"unable to open file");
//    boost::thread_group group;
//    group.create_thread(bind(save_saddle_mfolds,boost::ref(fs),boost::ref(*this),boost::ref(*msc)));

    opencl::update_to_surv_extrema(shared_from_this(),msc);

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

//    group.join_all();
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
    m_vert_fns.resize((m_ext_rect.span()/2)+1);
    m_owner_maxima.resize(m_rect.span()/2);
    m_owner_minima.resize((m_rect.span()/2)+1);

    m_cell_flags.reindex(m_ext_rect.lc());
    m_vert_fns.reindex(m_ext_rect.lc()/2);
    m_owner_maxima.reindex(m_rect.lc()/2);
    m_owner_minima.reindex(m_rect.lc()/2);

    bin_read_marray(is,m_cell_flags.data(),m_ext_rect.span()+1);
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
