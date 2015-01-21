#ifndef GRID_DATASET_INL_H_INCLUDED
#define GRID_DATASET_INL_H_INCLUDED

#include <grid_dataset.h>

namespace grid
{

/*===========================================================================*/

inline cell_fn_t dataset_t::get_cell_fn(cellid_t c) const
{return m_vert_fns(get_cell_vert(c)/2);}

/*---------------------------------------------------------------------------*/

inline rect_t dataset_t::get_rect() const
{return m_rect;}

/*---------------------------------------------------------------------------*/

inline rect_t dataset_t::get_ext_rect() const
{return m_ext_rect;}

/*---------------------------------------------------------------------------*/

inline int c_to_i2(const rect_t &r,cellid_t c)
{ int X = (r[0].span())/2 +1,Y=(r[1].span())/2 +1;c = (c-r.lc())/2;
  return (X*Y*c[2] + X*c[1] + c[0]);}

/*---------------------------------------------------------------------------*/

inline cellid_t i_to_c2(const rect_t &r,int i)
{ int X = (r[0].span())/2 +1,Y = (r[1].span())/2 +1;
  return r.lc() + cellid_t(2*(i%X),2*((i%(X*Y))/X),2*(i/(X*Y)));}

/*---------------------------------------------------------------------------*/

inline int num_cells2(const rect_t &r)
{return c_to_i2(r,r.uc()) + 1;}

/*---------------------------------------------------------------------------*/

inline rect_t shrink(const rect_t &r, const cellid_t& s = cellid_t::one)
{return rect_t(r.lc()+s,r.uc()-s);}

/*---------------------------------------------------------------------------*/

template<> inline int_marray_t & dataset_t::owner_extrema<DES>()
{return m_owner_maxima;}
template<> inline int_marray_t & dataset_t::owner_extrema<ASC>()
{return m_owner_minima;}

/*---------------------------------------------------------------------------*/
template <> inline rect_t dataset_t::get_extrema_rect<DES>() const
{return shrink(m_rect);}
template <> inline rect_t dataset_t::get_extrema_rect<ASC>() const
{return (m_rect);}

/*---------------------------------------------------------------------------*/

template<> inline uint dataset_t::get_cets<GDIR_DES>
(cellid_t c,cellid_t *cets) const {return getCellFacets(c,cets);}
template<> inline uint dataset_t::get_cets<GDIR_ASC>
(cellid_t c,cellid_t *cets) const {return getCellCofacets(c,cets);}

/*---------------------------------------------------------------------------*/

template<eGDIR dir> inline uint dataset_t::get_co_cets
(cellid_t c,cellid_t *cets) const
{return get_cets<(dir== DES)?(ASC):(DES)>(c,cets);}

/*---------------------------------------------------------------------------*/

template<> inline uint dataset_t::get_points<CC_PRIM>
(cellid_t c,cellid_t *pts) const {return getCellPoints(c,pts);}
template<> inline uint dataset_t::get_points<CC_DUAL>
(cellid_t c,cellid_t *pts) const {return getCellCubes(c,pts);}

/*---------------------------------------------------------------------------*/

inline rect_t dataset_t::get_extrema_rect(eGDIR dir) const
{return (dir == GDIR_DES)?(rect_t(m_rect.lc()+1,m_rect.uc()-1)):(m_rect);}

/*---------------------------------------------------------------------------*/

inline int dataset_t::getCellDim ( cellid_t c ) const
{return get_cell_dim(c);}

/*---------------------------------------------------------------------------*/

template <int sad_dim>
inline void  get_adj_extrema(cellid_t c, cellid_t & e1,cellid_t & e2)
{
  const int a = (sad_dim == 2)?(1):(0);

  ASSERT(get_cell_dim(c) == sad_dim );

  e1[0] = c[0] + ((c[0]+a)&1);
  e1[1] = c[1] + ((c[1]+a)&1);
  e1[2] = c[2] + ((c[2]+a)&1);

  e2[0] = c[0] - ((c[0]+a)&1);
  e2[1] = c[1] - ((c[1]+a)&1);
  e2[2] = c[2] - ((c[2]+a)&1);

  //    ASSERT(dir != GDIR_ASC || (get_cell_dim(e1) == 3  && get_cell_dim(e2) == 3));
  //    ASSERT(dir != GDIR_DES || (get_cell_dim(e1) == 0  && get_cell_dim(e2) == 0));
}

/*---------------------------------------------------------------------------*/

inline void dataset_t::save(const std::string &f) const
{
  std::fstream fs(f.c_str(),std::ios::out|std::ios::binary);
  ENSUREV(fs.is_open(),"file not found!!",f);
  save_bin(fs);
}

/*---------------------------------------------------------------------------*/

inline void dataset_t::load(const std::string &f)
{
  std::fstream fs(f.c_str(),std::ios::in|std::ios::binary);
  ENSUREV(fs.is_open(),"file not found!!",f);
  load_bin(fs);
}

/*===========================================================================*/



/*===========================================================================*/

template <int dim>
inline bool dataset_t::compare_cells_orig(const cellid_t & c1, const cellid_t &c2) const
{
  cellid_t f1 = getCellMaxFacetId(c1);
  cellid_t f2 = getCellMaxFacetId(c2);

  if(f1 != f2)
    return compare_cells_orig<dim-1>(f1,f2);

  f1 = getCellSecondMaxFacetId(c1);
  f2 = getCellSecondMaxFacetId(c2);

  int boundry_ct1 = m_domain_rect.boundryCount(f1);
  int boundry_ct2 = m_domain_rect.boundryCount(f2);

  if(boundry_ct1 != boundry_ct2)
    return (boundry_ct1 <boundry_ct2);

  return compare_cells_orig<dim-1>(f1,f2);
}

/*---------------------------------------------------------------------------*/

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

/*---------------------------------------------------------------------------*/

template <int dim>
inline bool dataset_t::compare_cells_pp(const cellid_t & c1, const cellid_t &c2) const
{
  cellid_t oc1 = c1,oc2 = c2;

  if(isCellPaired(c1) && getCellDim(getCellPairId(c1)) == dim+1)
    oc1 = getCellMaxFacetId(getCellPairId(c1));

  if(isCellPaired(c2) && getCellDim(getCellPairId(c2)) == dim+1)
    oc2 = getCellMaxFacetId(getCellPairId(c2));

  return compare_cells_orig<dim>(oc1,oc2);
}

/*---------------------------------------------------------------------------*/

template <> inline bool dataset_t::compare_cells_pp_<GDIR_DES,0>
(cellid_t c1, cellid_t c2) const  {return compare_cells_pp<0>(c1,c2);}
template <> inline bool dataset_t::compare_cells_pp_<GDIR_DES,1>
(cellid_t c1, cellid_t c2) const  {return compare_cells_pp<1>(c1,c2);}
template <> inline bool dataset_t::compare_cells_pp_<GDIR_DES,2>
(cellid_t c1, cellid_t c2) const  {return compare_cells_pp<2>(c1,c2);}
template <> inline bool dataset_t::compare_cells_pp_<GDIR_DES,3>
(cellid_t c1, cellid_t c2) const  {return compare_cells_pp<3>(c1,c2);}
template <> inline bool dataset_t::compare_cells_pp_<GDIR_ASC,0>
(cellid_t c1, cellid_t c2) const  {return compare_cells_pp<0>(c2,c1);}
template <> inline bool dataset_t::compare_cells_pp_<GDIR_ASC,1>
(cellid_t c1, cellid_t c2) const  {return compare_cells_pp<1>(c2,c1);}
template <> inline bool dataset_t::compare_cells_pp_<GDIR_ASC,2>
(cellid_t c1, cellid_t c2) const  {return compare_cells_pp<2>(c2,c1);}
template <> inline bool dataset_t::compare_cells_pp_<GDIR_ASC,3>
(cellid_t c1, cellid_t c2) const  {return compare_cells_pp<3>(c2,c1);}

/*===========================================================================*/
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
