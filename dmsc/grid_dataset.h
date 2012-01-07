/***************************************************************************
 *   Copyright (C) 2009 by nithin,,,   *
 *   nithin@gauss   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifndef __GRID_DATASET_H_INCLUDED_
#define __GRID_DATASET_H_INCLUDED_

#include <vector>
#include <queue>

#include <boost/enable_shared_from_this.hpp>
#include <boost/function.hpp>

#include <grid.h>

namespace grid
{
  class dataset_t:public boost::enable_shared_from_this<dataset_t>
  {

  public:

    // used as a bit mask.. cells can be critical and paired..in theory they all are
    enum eCellFlags
    {
      CELLFLAG_VISITED   = 0x80,
      CELLFLAG_CRITICAL = 0x40,
      CELLFLAG_MASK     = 0xc0
    };

    // bits [0,3) max facet of a cell
    // bits [3,6) pair of a cell
    // bit 6 ..  mark bit used by bfs to say visted or not
    // bit 7 .. is cell critical or not.


    typedef boost::multi_array<cellid_t,gc_grid_dim>        cellid_array_t;
    typedef boost::multi_array<cell_flag_t,gc_grid_dim>     cellflag_array_t;
    typedef boost::multi_array<cell_fn_t,gc_grid_dim>       varray_t;
    typedef cellid_list_t                                   mfold_t;


  public:

    rect_t             m_rect;
    rect_t             m_ext_rect;
    rect_t             m_domain_rect;
    rect_t             m_work_rect;

    varray_t           m_vert_fns;
    cellflag_array_t   m_cell_flags;
    cellflag_array_t   m_cell_order;

    int_marray_t       m_owner_maxima;
    int_marray_t       m_owner_minima;

  public:

    // initialization of the dataset
    dataset_t ( const rect_t &r,const rect_t &e,const rect_t &d );
    ~dataset_t ();

    void  init(const std::string &filename);
    void  clear();

  // the actual work routines
  public:
    void  computeMsGraph(mscomplex_ptr_t msgraph);

    void  compute_owner_grad();

    void  saveManifolds(mscomplex_ptr_t msgraph,const std::string &);

    void  saveConnectingOneManifolds(mscomplex_ptr_t ,const std::string &);

  // subroutines to main functions
  public:

    void  assignMaxFacets_thd(int tid,int dim);

    void  pairCellsWithinEst_thd(int tid,cellid_list_t * ccells);

    void  assign_pairs2();

    void  markBoundry_thd(int tid,rect_t bnd,cellid_list_t * ccells);

    void  setupCPs(mscomplex_ptr_t msgraph,cellid_list_t * ccells,int offset);

  // dataset interface
  public:

    cellid_t   getCellPairId ( cellid_t ) const;

    cellid_t   getCellMaxFacetId ( cellid_t ) const;

    cellid_t   getCellSecondMaxFacetId ( cellid_t ) const;

    uint   getCellPoints ( cellid_t ,cellid_t  * ) const;

    uint   getCellFacets ( cellid_t ,cellid_t * ) const;

    inline uint   getCellIncCells( cellid_t ,cellid_t * ) const;

    uint   getCellCofacets ( cellid_t ,cellid_t * ) const;

    uint   getCellCofaces ( cellid_t ,cellid_t * ) const;

    uint   getCellEst (cellid_t,cellid_t*) const;

    bool   isPairOrientationCorrect ( cellid_t c, cellid_t p ) const;

    bool   isCellCritical ( cellid_t c ) const;

    bool   isCellPaired ( cellid_t c ) const;

    bool   isCellVisited ( cellid_t c ) const;

    bool   areCellsIncident(cellid_t c1,cellid_t c2) const;

    void   pairCells ( cellid_t c,cellid_t p );

    void   visitCell( cellid_t c);

    void   setCellMaxFacet (cellid_t c,cellid_t f);

    void   markCellCritical ( cellid_t c );

    inline int getCellDim ( cellid_t c ) const;

    bool   isTrueBoundryCell ( cellid_t c ) const;

    bool   isFakeBoundryCell ( cellid_t c ) const;

    bool   isCellExterior ( cellid_t c ) const;

    // misc functions
  public:
    inline cellid_t get_cell_vert(cellid_t c) const
    {
      cellid_t v = c;

      switch(getCellDim(c))
      {
        case 3: v = getCellMaxFacetId(v);
        case 2: v = getCellMaxFacetId(v);
        case 1: v = getCellMaxFacetId(v);
      }
      return v;
    }

    inline cell_fn_t get_cell_fn(cellid_t c)
    {
      return m_vert_fns(get_cell_vert(c)/2);
    }

    inline rect_t get_rect()
    {
      return m_rect;
    }

    inline rect_t get_ext_rect()
    {
      return m_ext_rect;
    }

    inline rect_t get_extrema_rect(eGDIR dir);

    void log_flags();

    void log_pairs(std::ostream &os = std::cout);
    void log_pairs(const std::string &s);

    void log_visits(std::ostream &os = std::cout);
    void log_visits(const std::string &s);

    void log_pair_visits(std::ostream &os = std::cout);
    void log_pair_visits(const std::string &s);

    void log_owner_extrema(eGDIR dir, std::ostream &os = std::cout);
    void log_owner_extrema(eGDIR dir, const std::string &s);

    void log_max_facets();

    void extract_vdata_subarray(rect_t r,const std::string &filename);

    typedef rect_t::pt_iterator iterator;
    inline iterator begin() const;
    inline iterator end() const;

    class iterator_dim;
    inline iterator_dim begin(int d) const;
    inline iterator_dim end(int d) const;

    template <int dim>
    inline bool compare_cells_orig(const cellid_t & c1, const cellid_t &c2) const;

    template <int dim>
    inline bool compare_cells(const cellid_t & c1, const cellid_t &c2) const;

    template<eGDIR dir>
    inline uint get_cets(cellid_t c,cellid_t *cets) const;

    template<eGDIR dir>
    inline uint get_co_cets(cellid_t c,cellid_t *cets) const;

    template <eGDIR dir>
    inline rect_t get_extrema_rect() const;

    template<eGDIR dir>
    inline int_marray_t & owner_extrema() ;

    void stow(std::ostream &os);
    void load(std::istream &is);

    inline void stow(const std::string &f)
    {std::fstream fs(f.c_str(),std::ios::out|std::ios::binary);stow(fs);}
    inline void load(const std::string &f)
    {std::fstream fs(f.c_str(),std::ios::in|std::ios::binary);load(fs);}

  };

  inline int c_to_i2(const rect_t &r,cellid_t c)
  {
    int X = (r[0].span())/2 +1;
    int Y = (r[1].span())/2 +1;

    c = (c-r.lc())/2;

    return (X*Y*c[2] + X*c[1] + c[0]);
  }

  inline cellid_t i_to_c2(const rect_t &r,int i)
  {
    int X = (r[0].span())/2 +1;
    int Y = (r[1].span())/2 +1;

    return r.lc() + cellid_t(2*(i%X),2*((i%(X*Y))/X),2*(i/(X*Y)));
  }

  inline int num_cells2(const rect_t &r)
  {return c_to_i2(r,r.uc()) + 1;}

  inline rect_t shrink(const rect_t &r, const cellid_t& s = cellid_t::one)
  {return rect_t(r.lc()+s,r.uc()-s);}


  template <int dim>
  inline bool dataset_t::compare_cells_orig(const cellid_t & c1, const cellid_t &c2) const
  {
    cellid_t f1 = getCellMaxFacetId(c1);
    cellid_t f2 = getCellMaxFacetId(c2);

    if(f1 != f2)
      return compare_cells<dim-1>(f1,f2);

    f1 = getCellSecondMaxFacetId(c1);
    f2 = getCellSecondMaxFacetId(c2);

    int boundry_ct1 = m_domain_rect.boundryCount(f1);
    int boundry_ct2 = m_domain_rect.boundryCount(f2);

    if(boundry_ct1 != boundry_ct2)
      return (boundry_ct1 <boundry_ct2);

    return compare_cells<dim-1>(f1,f2);
  }

  template <>
  inline bool dataset_t::compare_cells_orig<0>(const cellid_t & c1, const cellid_t &c2) const
  {

    ASSERT(get_cell_dim(c1) == 0);
    ASSERT(get_cell_dim(c2) == 0);

    cell_fn_t f1 = m_vert_fns(c1/2);
    cell_fn_t f2 = m_vert_fns(c2/2);

    if (f1 != f2)
      return f1 < f2;

    return c1 < c2;
  }

  template <int dim>
  inline bool dataset_t::compare_cells(const cellid_t & c1, const cellid_t &c2) const
  {
    ASSERT(m_work_rect.contains(c1));
    ASSERT(m_work_rect.contains(c2));

    cellid_t f1 = getCellMaxFacetId(c1);
    cellid_t f2 = getCellMaxFacetId(c2);


    if(f1 == f2)
      return m_cell_order(c1) < m_cell_order(c2);

    return compare_cells<dim-1>(f1,f2);
  }

  template <>
  inline bool dataset_t::compare_cells<0>(const cellid_t & c1, const cellid_t &c2) const
  {
    return compare_cells_orig<0>(c1,c2);
  }

  template <>
  inline bool dataset_t::compare_cells<-1>(const cellid_t & c1, const cellid_t &c2) const
  {
    cellid_t v1 = get_cell_vert(c1);
    cellid_t v2 = get_cell_vert(c2);

    if(v1 == v2)
      return m_cell_order(c1) < m_cell_order(c2);

    return compare_cells<0>(v1,v2);
  }

  template<>
  inline int_marray_t & dataset_t::owner_extrema<GDIR_DES>()
  {return m_owner_maxima;}

  template<>
  inline int_marray_t & dataset_t::owner_extrema<GDIR_ASC>()
  {return m_owner_minima;}

  template <>
  inline rect_t dataset_t::get_extrema_rect<GDIR_DES>() const
  {return shrink(m_rect);}

  template <>
  inline rect_t dataset_t::get_extrema_rect<GDIR_ASC>() const
  {return (m_rect);}

  template<>
  inline uint dataset_t::get_cets<GDIR_DES>(cellid_t c,cellid_t *cets) const
  {return getCellFacets(c,cets);}

  template<>
  inline uint dataset_t::get_cets<GDIR_ASC>(cellid_t c,cellid_t *cets) const
  {return getCellCofacets(c,cets);}

  template<eGDIR dir>
  inline uint dataset_t::get_co_cets(cellid_t c,cellid_t *cets) const
  {return get_cets<(dir == GDIR_DES)?(GDIR_ASC):(GDIR_DES)>(c,cets);}

  inline rect_t dataset_t::get_extrema_rect(eGDIR dir)
  {return (dir == GDIR_DES)?(rect_t(m_rect.lc()+1,m_rect.uc()-1)):(m_rect);}

  inline int dataset_t::getCellDim ( cellid_t c ) const
  {return get_cell_dim(c);}

  class dataset_t::iterator_dim:public std::iterator
      <std::forward_iterator_tag,cellid_t,int,cellid_t,cellid_t>
  {
  public:
    iterator_dim(rect_t r,cellid_t d,int pn = 0):m_rect(r),m_dir(d),m_i(0),m_pn(pn)
    {
      m_cur_rect  = shrink(m_rect,m_dir);
      m_num_cells = num_cells2(m_cur_rect);
    };
    rect_t   m_rect;
    cellid_t m_dir;

    int     m_num_cells;
    rect_t  m_cur_rect;

    int m_i;
    int m_pn;

    inline iterator_dim& operator++()
    {
      ++m_i;

      if(m_i == m_num_cells)
      {
        std::next_permutation(m_dir.begin(),m_dir.end());

        m_cur_rect  = shrink(m_rect,m_dir);
        m_num_cells = num_cells2(m_cur_rect);

        m_i = 0;
        m_pn++;
      }
      return *this;
    }

    inline bool operator== (const iterator_dim &rhs) const
    {return (m_i == rhs.m_i) && (m_pn == rhs.m_pn);}

    inline bool operator!= (const iterator_dim &rhs) const
    {return !(*this == rhs);}

    inline reference operator*() const {return i_to_c2(m_cur_rect,m_i);}
  };


  inline dataset_t::iterator dataset_t::begin() const
  {return m_rect.pt_begin();}

  inline dataset_t::iterator dataset_t::end() const
  {return m_rect.pt_end();}

  inline dataset_t::iterator_dim dataset_t::begin(int d) const
  {
    switch(d)
    {
      case 0: return iterator_dim(m_rect,cellid_t(0,0,0),0);
      case 1: return iterator_dim(m_rect,cellid_t(0,0,1),0);
      case 2: return iterator_dim(m_rect,cellid_t(0,1,1),0);
      case 3: return iterator_dim(m_rect,cellid_t(1,1,1),0);
    }

    ASSERT(false);

    return iterator_dim(m_rect,cellid_t(-1,-1,-1),0);
  }

  inline dataset_t::iterator_dim dataset_t::end(int d) const
  {
    switch(d)
    {
      case 0: return iterator_dim(m_rect,cellid_t(0,0,0),1);
      case 1: return iterator_dim(m_rect,cellid_t(0,0,1),3);
      case 2: return iterator_dim(m_rect,cellid_t(0,1,1),3);
      case 3: return iterator_dim(m_rect,cellid_t(1,1,1),1);
    }

    ASSERT(false);

    return iterator_dim(m_rect,cellid_t(-1,-1,-1),0);
  }

  inline void get_boundry_rects(const rect_t &r,const rect_t & e,rect_list_t &bnds)
  {
    for( int xyz_dir = 0 ; xyz_dir < 3; ++xyz_dir)
    {
      for( int lr_dir = 0 ; lr_dir < 2; ++lr_dir)
      {
        rect_t bnd = r;

        if(r[xyz_dir][lr_dir] != e[xyz_dir][lr_dir])
        {
          bnd[xyz_dir][0] = r[xyz_dir][lr_dir];
          bnd[xyz_dir][1] = r[xyz_dir][lr_dir];

          bnds.push_back(bnd);
        }
      }
    }
  }

  template <eGDIR dir>
  inline void  get_adj_extrema(cellid_t c, cellid_t & e1,cellid_t & e2)
  {
    const int a = (dir == GDIR_ASC)?(1):(0);

    ASSERT(dir != GDIR_ASC || get_cell_dim(c) == 2 );
    ASSERT(dir != GDIR_DES || get_cell_dim(c) == 1 );

    e1[0] = c[0] + ((c[0]+a)&1);
    e1[1] = c[1] + ((c[1]+a)&1);
    e1[2] = c[2] + ((c[2]+a)&1);

    e2[0] = c[0] - ((c[0]+a)&1);
    e2[1] = c[1] - ((c[1]+a)&1);
    e2[2] = c[2] - ((c[2]+a)&1);

    ASSERT(dir != GDIR_ASC || (get_cell_dim(e1) == 3  && get_cell_dim(e2) == 3));
    ASSERT(dir != GDIR_DES || (get_cell_dim(e1) == 0  && get_cell_dim(e2) == 0));
  }

}

//namespace boost
//{
//  namespace serialization
//  {
//    template<class Archive>
//    void serialize(Archive & ar, grid::dataset_t & d, const unsigned int );
//
//  } // namespace serialization
//}
#endif
