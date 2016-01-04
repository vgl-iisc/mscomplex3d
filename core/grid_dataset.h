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
#include <fstream>
#include <stack>

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


public:

  rect_t             m_rect;
  rect_t             m_ext_rect;
  rect_t             m_domain_rect;

  varray_t           m_vert_fns;
  cellflag_array_t   m_cell_flags;

  int_marray_t       m_owner_maxima;
  int_marray_t       m_owner_minima;

public:

  // initialization of the dataset
  dataset_t ( const rect_t &r,const rect_t &e,const rect_t &d );
  ~dataset_t ();

  void  init(const cell_fn_t * dptr,bool is_fortran_order);
  void  init(const std::string &filename);
  void  init_storage(); // assumes rects are defined
  void  clear();

  // save/load data
  inline void save(const std::string &f) const;
  inline void load(const std::string &f);

  void save_bin(std::ostream &os) const;
  void load_bin(std::istream &is);

  // dataset base functions

  cellid_t   getCellPairId ( cellid_t ) const;
  cellid_t   getCellMaxFacetId ( cellid_t ) const;
  cellid_t   getCellSecondMaxFacetId ( cellid_t ) const;
  uint       getCellPoints ( cellid_t ,cellid_t  * ) const;
  uint       getCellCubes ( cellid_t ,cellid_t  * ) const;
  uint       getCellFacets ( cellid_t ,cellid_t * ) const;
  uint       getCellIncCells( cellid_t ,cellid_t * ) const;
  uint       getCellCofacets ( cellid_t ,cellid_t * ) const;
  uint       getCellCofaces ( cellid_t ,cellid_t * ) const;
  uint       getCellEst (cellid_t,cellid_t*) const;
  bool       isPairOrientationCorrect ( cellid_t c, cellid_t p ) const;
  bool       isCellCritical ( cellid_t c ) const;
  bool       isCellPaired ( cellid_t c ) const;
  bool       isCellVisited ( cellid_t c ) const;
  bool       areCellsIncident(cellid_t c1,cellid_t c2) const;
  void       pairCells ( cellid_t c,cellid_t p );
  void       visitCell( cellid_t c);
  void       setCellMaxFacet (cellid_t c,cellid_t f);
  void       markCellCritical ( cellid_t c );
  int        getCellDim ( cellid_t c ) const;
  bool       isTrueBoundryCell ( cellid_t c ) const;
  bool       isFakeBoundryCell ( cellid_t c ) const;
  bool       isCellExterior ( cellid_t c ) const;
  cellid_t   get_cell_vert(cellid_t c) const;


  inline cell_fn_t      get_cell_fn(cellid_t c) const;
  inline rect_t         get_rect() const;
  inline rect_t         get_ext_rect() const;
  inline rect_t         get_extrema_rect(eGDIR dir) const;
  template<eGDIR dir>
  inline uint           get_cets(cellid_t c,cellid_t *cets) const;
  template<eCCTYPE dir>
  inline uint           get_points(cellid_t c,cellid_t *cets) const;
  template<eGDIR dir>
  inline uint           get_co_cets(cellid_t c,cellid_t *cets) const;
  template <eGDIR dir>
  inline rect_t         get_extrema_rect() const;
  template<eGDIR dir>
  inline int_marray_t & owner_extrema() ;


  // Comparator routines that establish total order before/after pairing
  template <int dim>
  inline bool compare_cells_orig(const cellid_t & c1, const cellid_t &c2) const;

  template <int dim>
  inline bool compare_cells_pp(const cellid_t & c1, const cellid_t &c2) const;

  template <eGDIR dir, int dim>
  inline bool compare_cells_pp_(cellid_t c1, cellid_t c2) const;

  // core work algorithms
  void  computeMsGraph(mscomplex_ptr_t msgraph);
  void  getManifold(mfold_t &mfold,const cellid_list_t &rng,int dim,eGDIR dir);


  // additional work algorithms for out of core processing
  void  compute_owner_grad();
  void  saveManifolds(mscomplex_ptr_t msgraph,const std::string &);
  void  storeOwnerArrays(int_marray_t &,int_marray_t &) const;
  void  loadOwnerArrays(int_marray_t &,int_marray_t &);


  // misc functions
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

void get_boundry_rects(const rect_t &r,const rect_t & e,rect_list_t &bnds);
}

#include <grid_dataset_inl.h>
#endif
