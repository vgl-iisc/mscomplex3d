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

#include <boost/enable_shared_from_this.hpp>
#include <boost/multi_array.hpp>
#include <boost/function.hpp>

#include <grid.h>

namespace grid
{

  class mscomplex_t;

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
    typedef boost::multi_array<int,gc_grid_dim>             owner_array_t;
    typedef cellid_list_t                                   mfold_t;


  public:

    rect_t             m_rect;
    rect_t             m_ext_rect;
    rect_t             m_domain_rect;

    varray_t           m_vert_fns;

    cellflag_array_t   m_cell_flags;

    owner_array_t      m_owner_maxima;
    owner_array_t      m_owner_minima;

    boost::function<bool (cellid_t,cellid_t)> cmp_ftor;
    boost::function<bool (cellid_t,cellid_t)> cmp_ftors[2];


  public:

    // initialization of the dataset
    dataset_t ( const rect_t &r,const rect_t &e,const rect_t &d );
    ~dataset_t ();

    void  init(const std::string &filename);
    void  clear();

  // the actual work routines
  public:
    void  computeMsGraph(mscomplex_ptr_t msgraph);

    void  saveManifolds(mscomplex_ptr_t msgraph,const std::string &);

  // subroutines to main functions
  public:

    void  assignMaxFacets_thd(int tid,int dim);

    void  pairCellsWithinEst_thd(int tid,cellid_list_t * ccells);

    void  markBoundry_thd(int tid,rect_t bnd,cellid_list_t * ccells);

    void  setupCPs(mscomplex_ptr_t msgraph,cellid_list_t * ccells,int offset);

    void  extrema_connect_thd(mscomplex_ptr_t msgraph,cp_producer_ptr_t p);

    void  get_mfold(mfold_t * mfold,mscomplex_const_ptr_t msc,int i,int dir) const;

    void  mark_extrema_owner_thd(mscomplex_ptr_t msgraph,cp_producer_ptr_t p);
  // dataset interface
  public:

    cellid_t   getCellPairId ( cellid_t ) const;

    cellid_t   getCellMaxFacetId ( cellid_t ) const;

    cellid_t   getCellSecondMaxFacetId ( cellid_t ) const;

    inline bool   ptLt ( cellid_t c1,cellid_t c2) const
    {
      cell_fn_t f1 = m_vert_fns(c1/2);
      cell_fn_t f2 = m_vert_fns(c2/2);

      if (f1 != f2)
        return f1 < f2;

      return c1<c2;
    }

    bool   compareCells( cellid_t ,cellid_t ) const;

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
    inline cellid_t get_cell_vert(cellid_t c)
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
  };

  inline rect_t dataset_t::get_extrema_rect(eGDIR dir)
  {
    return (dir == GDIR_DES)?(rect_t(m_rect.lc()+1,m_rect.uc()-1)):(m_rect);
  }

  inline int get_cell_dim ( cellid_t c )
  {
    int dim = 0;

    for (size_t i = 0 ; i < gc_grid_dim;++i)
      dim += ( c[i]&0x01 );

    return (dim);
  }


  inline int dataset_t::getCellDim ( cellid_t c ) const
  {
    return get_cell_dim(c);
  }

  inline int c_to_i(const rect_t &r,cellid_t c)
  {
    cellid_t s = r.span()+1;
    c = (c - r.lc());
    return (s[0]*s[1]*c[2] + s[0]*c[1] + c[0]);
  }

  inline cellid_t i_to_c(const rect_t &r,int i)
  {
    cellid_t s = r.span()+1;
    cellid_t c = r.lc() + (cellid_t(i%s[0],(i%(s[0]*s[1]))/s[0],i/(s[0]*s[1])));
    ASSERT(r.contains(c));
    return c;
  }

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

  inline int num_cells(const rect_t &r)
  {
    return c_to_i(r,r.uc()) + 1;
  }

  inline int num_cells2(const rect_t &r)
  {
    return c_to_i2(r,r.uc()) + 1;
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

  inline void  get_adj_extrema(cellid_t c, cellid_t & e1,cellid_t & e2,eGDIR dir)
  {
    ASSERT(dir != GDIR_DES || get_cell_dim(c) == 2 );
    ASSERT(dir != GDIR_ASC || get_cell_dim(c) == 1 );

    int a = (dir == GDIR_DES)?(1):(0);

    e1[0] = c[0] + ((c[0]+a)&1);
    e1[1] = c[1] + ((c[1]+a)&1);
    e1[2] = c[2] + ((c[2]+a)&1);

    e2[0] = c[0] - ((c[0]+a)&1);
    e2[1] = c[1] - ((c[1]+a)&1);
    e2[2] = c[2] - ((c[2]+a)&1);
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
